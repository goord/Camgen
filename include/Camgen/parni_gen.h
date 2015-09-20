//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file parni_gen.h
    \brief Monte-Carlo generator class using an adaptive, unequal-spaced grid.
 */

#ifndef CAMGEN_PARNI_GEN_H_
#define CAMGEN_PARNI_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Parni is a Monte-Carlo generator class, which uses an adaptive grid with  *
 * unequal spacings. When calling the update method, the bin weights are     *
 * updated, and when the adapt method is called, the highest-weight bin is   *
 * split, approximating the integrand optimally.                             *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <set>
#include <Camgen/obj_alloc.h>
#include <Camgen/MC_gen.h>
#include <Camgen/MC_int_base.h>
#include <Camgen/parni_it.h>

namespace Camgen
{
    /// Parni adaptive grid for multi-dimensional integration.

    template<class value_t,std::size_t D,class rng_t,class key_t>class parni_generator: public object_allocator<typename parni_bin<value_t,D,rng_t,key_t>::point_type>, 
    									      	     	public MC_generator<value_t>, 
									      	     	public MC_integrator_base<value_t>
    {
	friend class parni_sub_grid<value_t,D,rng_t,key_t>;

	/* Type definitions: */

	typedef object_allocator<typename parni_bin<value_t,D,rng_t>::point_type> base_type;

	public:

	    typedef value_t value_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef parni_bin<value_t,D,rng_t,key_t> bin_type;
	    typedef typename bin_type::point_type point_type;
	    typedef typename bin_type::bin_iterator bin_iterator;
	    typedef typename bin_type::const_bin_iterator const_bin_iterator;
	    typedef typename bin_type::reverse_bin_iterator reverse_bin_iterator;
	    typedef typename bin_type::const_reverse_bin_iterator const_reverse_bin_iterator;
	    typedef typename bin_type::leaf_iterator leaf_iterator;
	    typedef typename bin_type::const_leaf_iterator const_leaf_iterator;
	    typedef typename bin_type::reverse_leaf_iterator reverse_leaf_iterator;
	    typedef typename bin_type::const_reverse_leaf_iterator const_reverse_leaf_iterator;

	    /* Public data members: */
	    /*----------------------*/

	    /// Maximal number of bins.

	    const size_type n_bins;

	    /// Splitting criterium
	    
	    const grid_modes::type mode;

	    /// Constructor with lower-left and upper-right integration domain
	    /// corners and the maximal number of bins as arguments.

	    /* Public constructors: */
	    /*----------------------*/
	    
	    parni_generator(const point_type& xmin_,const point_type& xmax_,size_type n_bins_=500,grid_modes::type mode_=grid_modes::cumulant_weights):n_bins(n_bins_),mode(mode_),xmin(xmin_),xmax(xmax_),root_bin(NULL),reg_bin(root_bin)
	    {
		vector<key_t,D>key;
		key.assign(1);
		root_bin=new bin_type(&xmin,&xmax,key,mode);
	    }

	    /// Constructor with an external point address (to be filled) and
	    /// lower-left and upper-right integration domain corners and the
	    /// maximal number of bins as arguments.
	    
	    parni_generator(point_type* point_,const point_type& xmin_,const point_type& xmax_,size_type n_bins_=500,grid_modes::type mode_=grid_modes::cumulant_weights):base_type(point_),n_bins(n_bins_),mode(mode_),xmin(xmin_),xmax(xmax_),root_bin(NULL),reg_bin(root_bin)
	    {
		vector<key_t,D>key;
		key.assign(1);
		root_bin=new bin_type(&xmin,&xmax,key,mode);
	    }

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    ~parni_generator()
	    {
		delete root_bin;
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Generation method implementation.

	    bool generate()
	    {
		value_type rho=rn_stream::throw_number(0,root_bin->w);
		reg_bin=root_bin->find_weight(rho);
		if(reg_bin->generate(this->object()))
		{
		    this->weight()=reg_bin->weight()*norm();
		    return true;
		}
		this->weight()=(value_type)0;
		return false;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		reg_bin=root_bin->find_point(this->object());
		if(reg_bin==NULL)
                {
                    return false;
                }
		this->weight()=reg_bin->weight()*norm();
		return true;
	    }

	    /// Overridden updating method.

	    void update()
	    {
                if(reg_bin==NULL)
                {
                    return;
                }
		if(accept(this->integrand()))
		{
		    reg_bin->update(this->object(),this->integrand());
		}
	    }

	    /// Overridden adaptation method.

	    void adapt()
	    {
		root_bin->adapt();
		if(!final())
		{
		    std::vector<bin_type*>selection;
		    value_type maxw(0),w;
		    for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		    {
			w=(*it)->w;
			if(maxw<w)
			{
			    selection.resize(1);
			    selection[0]=*it;
			    maxw=w;
			}
			else if(maxw==w)
			{
			    selection.push_back(*it);
			}
		    }
		    if(selection.size()==1)
		    {
			selection[0]->split();
		    }
		    else if(selection.size()>1)
		    {
			selection[rn_stream::throw_dice(selection.size())]->split();
		    }
		}
		else
		{
		    std::vector<bin_type*>min_selection,max_selection;
		    value_type maxw(0),minw(std::numeric_limits<value_type>::infinity()),w;
		    for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		    {
			w=(*it)->w;
			if(maxw<w)
			{
			    max_selection.resize(1);
			    max_selection[0]=*it;
			    maxw=w;
			}
			else if(maxw==w)
			{
			    max_selection.push_back(*it);
			}
			w=(*it)->parent->w;
			if(minw>w)
			{
			    min_selection.resize(1);
			    min_selection[0]=(*it)->parent;
			    minw=w;
			}
			else if(minw==w)
			{
			    min_selection.push_back(*it);
			}
		    }
		    bin_type* minb;
		    if(min_selection.size()==1)
		    {
			minb=min_selection[0];
		    }
		    else if(min_selection.size()>1)
		    {
			minb=min_selection[rn_stream::throw_dice(min_selection.size())];
		    }
		    else
		    {
			return;
		    }
		    bin_type* maxb;
		    if(max_selection.size()==1)
		    {
			maxb=max_selection[0];
		    }
		    else if(max_selection.size()>1)
		    {
			maxb=max_selection[rn_stream::throw_dice(max_selection.size())];
		    }
		    else
		    {
			return;
		    }
		    if(!(minb==maxb or maxb->is_child_of(minb)))
		    {
			if(reg_bin==NULL or reg_bin->is_child_of(minb))
			{
			    reg_bin=minb;
			}
			minb->merge();
			maxb->split();
			if(reg_bin==maxb)
			{
			    reg_bin=rn_stream::throw_coin()?(maxb->child1):(maxb->child2);
			}
		    }
		}
	    }

	    /// Resets the cross section and re-initialises the grid to a single
	    /// bin.

	    void reset()
	    {
		root_bin->reset();
	    }
	    
	    /// Creates a new sub grid generator with the current instance as parent.

	    parni_sub_grid<value_t,D,rng_t,key_t>* create_sub_grid(const point_type& xmin,const point_type& xmax)
	    {
		return new parni_sub_grid<value_t,D,rng_t,key_t>(this,xmin,xmax);
	    }

	    /* Public readout functions: */
	    /*---------------------------*/

	    /// Returns whether the maximal number of bins has been reached.

	    bool final() const
	    {
		return (root_bin->count_bins()>=n_bins);
	    }

	    /// Returns the number of bins in the tree.

	    size_type count_bins() const
	    {
		return root_bin->count_bins();
	    }

	    /// Returns the number of leaves in the tree.

	    size_type count_leaves() const
	    {
		return root_bin->count_leaves();
	    }

	    /// Weight normalisation.

	    const value_type& norm() const
	    {
		return root_bin->w;
	    }

	    /// Printing method.

	    std::ostream& print(std::ostream& os) const
	    {
		for(const_bin_iterator it=begin_bins();it!=end_bins();++it)
		{
		    os<<(**it)<<std::endl;
		}
		return os;
	    }

	    /// Method printing only child bins.

	    std::ostream& print_leaves(std::ostream& os) const
	    {
		for(const_leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		{
		    os<<(**it)<<std::endl;
		}
		return os;
	    }

	    /// Adds the rectangles to the plotscript pointer

	    plot_script* print_rectangles(plot_script* ps) const
	    {
		for(const_leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		{
		    (*it)->print_rectangle(ps);
		}
		return ps;
	    }

	    /* Iterator definitions: */
	    /*-----------------------*/

	    leaf_iterator begin_leaves()
	    {
		return root_bin->begin_leaves();
	    }
	    const_leaf_iterator begin_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->begin_leaves();
	    }
	    reverse_leaf_iterator rbegin_leaves()
	    {
		return root_bin->rbegin_leaves();
	    }
	    const_reverse_leaf_iterator rbegin_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->rbegin_leaves();
	    }
	    leaf_iterator end_leaves()
	    {
		return root_bin->end_leaves();
	    }
	    const_leaf_iterator end_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->end_leaves();
	    }
	    reverse_leaf_iterator rend_leaves()
	    {
		return root_bin->rend_leaves();
	    }
	    const_reverse_leaf_iterator rend_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->rend_leaves();
	    }
	    bin_iterator begin_bins()
	    {
		return root_bin->begin_bins();
	    }
	    const_bin_iterator begin_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->begin_bins();
	    }
	    reverse_bin_iterator rbegin_bins()
	    {
		return root_bin->rbegin_bins();
	    }
	    const_reverse_bin_iterator rbegin_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->rbegin_bins();
	    }
	    bin_iterator end_bins()
	    {
		return root_bin->end_bins();
	    }
	    const_bin_iterator end_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->end_bins();
	    }
	    reverse_bin_iterator rend_bins()
	    {
		return root_bin->rend_bins();
	    }
	    const_reverse_bin_iterator rend_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->rend_bins();
	    }

	    void check()
	    {
		for(bin_iterator it=begin_bins();it!=end_bins();++it)
		{
		    if((*it)->is_leaf())
		    {
			continue;
		    }
		    bin_type* c1=(*it)->child1;
		    bin_type* c2=(*it)->child2;
		    if((*it)->volume()!=(c1->volume()+c2->volume()))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"bin \n"<<**it<<"\n has volume unequal to the sum of its childs \n"<<*c1<<"\n"<<*c2<<endlog;
		    }
		    if((*it)->F0!=(c1->F0+c2->F0))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"bin \n"<<**it<<"\n has F0 unequal to the sum of its childs \n"<<*c1<<"\n"<<*c2<<endlog;
		    }
		    if((*it)->w!=(c1->w+c2->w))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"bin \n"<<**it<<"\n has weight unequal to the sum of its childs \n"<<*c1<<"\n"<<*c2<<endlog;
		    }
		}
		for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		{
		    if(!(*it)->is_leaf())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"leaf iterator \n"<<*it<<"\n is not a leaf bin"<<endlog;
		    }
		}
		value_type w=this->weight();
		evaluate_weight();
		if(w!=this->weight())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"re-evaluated weight "<<this->weight()<<" not equal to original "<<w<<endlog;
		}
	    }

	private:

	    /* Integration domain: */
	    
	    point_type xmin,xmax;

	    /* Root of the binary tree of bins: */

	    bin_type* root_bin;
	    
	    /* Register bin, which has thrown the last point: */
	    
	    bin_type* reg_bin;

	    /* Sub-grid register: */

	    std::set<parni_sub_grid<value_t,D,rng_t,key_t>*> subgrids;

	    /* Private constructors: */
	    /*-----------------------*/

	    /* Default constructor: */

	    parni_generator(size_type n_bins_,grid_modes::type mode_=grid_modes::cumulant_weights):n_bins(n_bins_),mode(mode_),root_bin(NULL),reg_bin(NULL){}

	    /* Constructor with point argument: */

	    parni_generator(point_type* point_,size_type n_bins_,grid_modes::type mode_=grid_modes::cumulant_weights):base_type(point_),n_bins(n_bins_),mode(mode_),root_bin(NULL),reg_bin(NULL){}

	    /* Private methods: */
	    /*------------------*/

	    /* Integrand acceptance function: */

	    static bool accept(const value_type& f)
	    {
		if(f==f and f!=std::numeric_limits<value_type>::infinity() and !(f<(value_type)0))
		{
		    return true;
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"integrand "<<f<<" was not accepted by parni"<<endlog;
		return false;
	    }
	    
	    /* Returns the lowest-level bin containing the argument point: */

	    bin_type* find_leaf(const point_type& x)
	    {
		return root_bin->find_point(x);
	    }

	    /* Returns the lowest-level bin containing the argument point and
	     * sets the second argument to the corresponding left-integral: */

	    bin_type* find_leaf(const point_type& x,value_type& s)
	    {
		s=(value_type)0;
		return root_bin->find_point(x,s);
	    }

	    /* Registers a subgrid: */

	    void register_sub_grid(parni_sub_grid<value_t,1,rng_t,key_t>* subgrid)
	    {
		if(subgrid->grid==this)
		{
		    subgrids.insert(subgrid);
		}
	    }

	    /* Unregisters a subgrid: */

	    void unregister_sub_grid(parni_sub_grid<value_t,1,rng_t,key_t>* subgrid)
	    {
		if(subgrid->grid==this)
		{
		    subgrids.erase(subgrid);
		}
	    }
    };

    /// Parni specialisation for one-dimensional integration.

    template<class value_t,class rng_t,class key_t>class parni_generator<value_t,1,rng_t,key_t>: public object_allocator<typename parni_bin<value_t,1,rng_t,key_t>::point_type>, 
											    	 public MC_generator<value_t>,
											    	 public MC_integrator_base<value_t>
    {
	friend class parni_sub_grid<value_t,1,rng_t,key_t>;

	/* Type definitions: */

	typedef object_allocator<typename parni_bin<value_t,1,rng_t,key_t>::point_type> base_type;

	public:

	    typedef value_t value_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef parni_bin<value_t,1,rng_t,key_t> bin_type;
	    typedef typename bin_type::point_type point_type;
	    typedef typename bin_type::bin_iterator bin_iterator;
	    typedef typename bin_type::const_bin_iterator const_bin_iterator;
	    typedef typename bin_type::reverse_bin_iterator reverse_bin_iterator;
	    typedef typename bin_type::const_reverse_bin_iterator const_reverse_bin_iterator;
	    typedef typename bin_type::leaf_iterator leaf_iterator;
	    typedef typename bin_type::const_leaf_iterator const_leaf_iterator;
	    typedef typename bin_type::reverse_leaf_iterator reverse_leaf_iterator;
	    typedef typename bin_type::const_reverse_leaf_iterator const_reverse_leaf_iterator;

	    /* Public data members: */
	    /*----------------------*/

	    /// Maximal number of bins.

	    const size_type n_bins;

	    /// Parni bin-splitting criterium.
	    
	    const grid_modes::type mode;

	    /// Constructor with lower-left and upper-right integration domain
	    /// corners and the maximal number of bins as arguments.

	    /* Public constructors: */
	    /*----------------------*/
	    
	    parni_generator(const point_type& xmin_,const point_type& xmax_,size_type n_bins_=500,grid_modes::type mode_=grid_modes::cumulant_weights):n_bins(n_bins_),mode(mode_),xmin(xmin_),xmax(xmax_),root_bin(new bin_type(&xmin,&xmax,1,mode)),reg_bin(root_bin){}

	    /// Constructor with an external point address (to be filled) and
	    /// lower-left and upper-right integration domain corners and the
	    /// maximal number of bins as arguments.
	    
	    parni_generator(point_type* point_,const point_type& xmin_,const point_type& xmax_,size_type n_bins_=500,grid_modes::type mode_=grid_modes::cumulant_weights):base_type(point_),n_bins(n_bins_),mode(mode_),xmin(xmin_),xmax(xmax_),root_bin(new bin_type(&xmin,&xmax,1,mode)),reg_bin(root_bin){}

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    ~parni_generator()
	    {
		delete root_bin;
		for(typename std::set<parni_sub_grid<value_t,1,rng_t,key_t>*>::iterator it=subgrids.begin();it!=subgrids.end();++it)
		{
		    (*it)->grid=NULL;
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Generation method implementation.

	    bool generate()
	    {
		value_type rho=rn_stream::throw_number(0,root_bin->w);
		reg_bin=root_bin->find_weight(rho);
		if(reg_bin->generate(this->object()))
		{
		    this->weight()=reg_bin->weight()*norm();
		    return true;
		}
		this->weight()=(value_type)0;
		return false;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		reg_bin=root_bin->find_point(this->object());
		if(reg_bin==NULL)
                {
                    return false;
                }
		this->weight()=reg_bin->weight()*norm();
		return true;
	    }

	    /// Overridden updating method.

	    void update()
	    {
                if(reg_bin==NULL)
                {
                    return;
                }
		if(accept(this->integrand()))
		{
		    reg_bin->update(this->object(),this->integrand());
		}
		for(typename std::set<parni_sub_grid<value_t,1,rng_t,key_t>*>::iterator it=subgrids.begin();it!=subgrids.end();++it)
		{
		    (*it)->refresh_bounds();
		}
	    }

	    /// Overridden adaptation method.
	    
	    void adapt()
	    {
		root_bin->adapt();
		if(!final())
		{
		    std::vector<bin_type*>selection;
		    value_type maxw(0),w;
		    for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		    {
			w=(*it)->w;
			if(maxw<w)
			{
			    selection.resize(1);
			    selection[0]=*it;
			    maxw=w;
			}
			else if(maxw==w)
			{
			    selection.push_back(*it);
			}
		    }
		    if(selection.size()==1)
		    {
			selection[0]->split();
		    }
		    else if(selection.size()>1)
		    {
			selection[rn_stream::throw_dice(selection.size())]->split();
		    }
		}
		else
		{
		    std::vector<bin_type*>min_selection,max_selection;
		    value_type maxw(0),minw(std::numeric_limits<value_type>::infinity()),w;
		    for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		    {
			w=(*it)->w;
			if(maxw<w)
			{
			    max_selection.resize(1);
			    max_selection[0]=*it;
			    maxw=w;
			}
			else if(maxw==w)
			{
			    max_selection.push_back(*it);
			}
			w=(*it)->parent->w;
			if(minw>w)
			{
			    min_selection.resize(1);
			    min_selection[0]=(*it)->parent;
			    minw=w;
			}
			else if(minw==w)
			{
			    min_selection.push_back((*it)->parent);
			}
		    }
		    bin_type* minb;
		    if(min_selection.size()==1)
		    {
			minb=min_selection[0];
		    }
		    else if(min_selection.size()>1)
		    {
			minb=min_selection[rn_stream::throw_dice(min_selection.size())];
		    }
		    else
		    {
			return;
		    }
		    bin_type* maxb;
		    if(max_selection.size()==1)
		    {
			maxb=max_selection[0];
		    }
		    else if(max_selection.size()>1)
		    {
			maxb=max_selection[rn_stream::throw_dice(max_selection.size())];
		    }
		    else
		    {
			return;
		    }
		    if(!(minb==maxb or maxb->is_child_of(minb)))
		    {
			if(reg_bin==NULL or reg_bin->is_child_of(minb))
			{
			    reg_bin=minb;
			}
			minb->merge();
			maxb->split();
			if(reg_bin==maxb)
			{
			    reg_bin=rn_stream::throw_coin()?(maxb->child1):(maxb->child2);
			}
		    }
		}
		for(typename std::set<parni_sub_grid<value_t,1,rng_t,key_t>*>::iterator it=subgrids.begin();it!=subgrids.end();++it)
		{
		    (*it)->refresh_bounds();
		}
	    }

	    /// Adapts the grid between the selected boundaries.

	    void adapt(const point_type& a,const point_type& b)
	    {
		root_bin->adapt();
		if(!final())
		{
		    std::vector<bin_type*>selection;
		    value_type maxw(0),w;
		    for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		    {
			if((*it)->upper_bound()<a or (*it)->lower_bound()>b)
			{
			    continue;
			}
			w=(*it)->w;
			if(maxw<w)
			{
			    selection.resize(1);
			    selection[0]=*it;
			    maxw=w;
			}
			else if(maxw==w)
			{
			    selection.push_back(*it);
			}
		    }
		    if(selection.size()==1)
		    {
			selection[0]->split();
		    }
		    else if(selection.size()>1)
		    {
			selection[rn_stream::throw_dice(selection.size())]->split();
		    }
		}
		else
		{
		    std::vector<bin_type*>min_selection,max_selection;
		    value_type maxw(0),minw(std::numeric_limits<value_type>::infinity()),w;
		    for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		    {
			if((*it)->upper_bound()<a or (*it)->lower_bound()>b)
			{
			    continue;
			}
			w=(*it)->w;
			if(maxw<w)
			{
			    max_selection.resize(1);
			    max_selection[0]=*it;
			    maxw=w;
			}
			else if(maxw==w)
			{
			    max_selection.push_back(*it);
			}
			w=(*it)->parent->w;
			if(minw>w)
			{
			    min_selection.resize(1);
			    min_selection[0]=(*it)->parent;
			    minw=w;
			}
			else if(minw==w)
			{
			    min_selection.push_back((*it)->parent);
			}
		    }
		    bin_type* minb;
		    if(min_selection.size()==1)
		    {
			minb=min_selection[0];
		    }
		    else if(min_selection.size()>1)
		    {
			minb=min_selection[rn_stream::throw_dice(min_selection.size())];
		    }
		    else
		    {
			return;
		    }
		    bin_type* maxb;
		    if(max_selection.size()==1)
		    {
			maxb=max_selection[0];
		    }
		    else if(max_selection.size()>1)
		    {
			maxb=max_selection[rn_stream::throw_dice(max_selection.size())];
		    }
		    else
		    {
			return;
		    }
		    if(!(minb==maxb or maxb->is_child_of(minb)))
		    {
			if(reg_bin==NULL or reg_bin->is_child_of(minb))
			{
			    reg_bin=minb;
			}
			minb->merge();
			maxb->split();
			if(reg_bin==maxb)
			{
			    reg_bin=rn_stream::throw_coin()?(maxb->child1):(maxb->child2);
			}
		    }
		}
		for(typename std::set<parni_sub_grid<value_t,1,rng_t,key_t>*>::iterator it=subgrids.begin();it!=subgrids.end();++it)
		{
		    (*it)->refresh_bounds();
		}
	    }

	    /// Resets the cross section and re-initialises the grid to a single
	    /// bin.

	    void reset()
	    {
		root_bin->reset();
	    }

	    /// Creates a new sub grid generator with the current instance as parent.

	    parni_sub_grid<value_t,1,rng_t,key_t>* create_sub_grid(const point_type& xmin,const point_type& xmax)
	    {
		return new parni_sub_grid<value_t,1,rng_t,key_t>(this,xmin,xmax);
	    }

	    /* Public readout functions: */
	    /*---------------------------*/

	    /// Returns whether the maximal number of bins has been reached.

	    bool final() const
	    {
		return (root_bin->count_bins()>=n_bins);
	    }

	    /// Returns the number of bins in the tree.

	    size_type count_bins() const
	    {
		return root_bin->count_bins();
	    }

	    /// Returns the number of leaves in the tree.

	    size_type count_leaves() const
	    {
		return root_bin->count_leaves();
	    }

	    /// Weight normalisation.

	    const value_type& norm() const
	    {
		return root_bin->w;
	    }

	    /// Printing method.

	    std::ostream& print(std::ostream& os) const
	    {
		for(const_bin_iterator it=begin_bins();it!=end_bins();++it)
		{
		    os<<(**it)<<std::endl;
		}
		return os;
	    }

	    /// Method printing only child bins.

	    std::ostream& print_leaves(std::ostream& os) const
	    {
		for(const_leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		{
		    os<<(**it)<<std::endl;
		}
		return os;
	    }

	    /// Adds the rectangles to the plotscript pointer

	    plot_script* print_rectangles(plot_script* ps) const
	    {
		for(const_leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		{
		    (*it)->print_rectangle(ps);
		}
		return ps;
	    }

	    /* Iterator definitions: */
	    /*-----------------------*/

	    leaf_iterator begin_leaves()
	    {
		return root_bin->begin_leaves();
	    }
	    const_leaf_iterator begin_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->begin_leaves();
	    }
	    reverse_leaf_iterator rbegin_leaves()
	    {
		return root_bin->rbegin_leaves();
	    }
	    const_reverse_leaf_iterator rbegin_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->rbegin_leaves();
	    }
	    leaf_iterator end_leaves()
	    {
		return root_bin->end_leaves();
	    }
	    const_leaf_iterator end_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->end_leaves();
	    }
	    reverse_leaf_iterator rend_leaves()
	    {
		return root_bin->rend_leaves();
	    }
	    const_reverse_leaf_iterator rend_leaves() const
	    {
		return const_cast<const bin_type*>(root_bin)->rend_leaves();
	    }
	    bin_iterator begin_bins()
	    {
		return root_bin->begin_bins();
	    }
	    const_bin_iterator begin_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->begin_bins();
	    }
	    reverse_bin_iterator rbegin_bins()
	    {
		return root_bin->rbegin_bins();
	    }
	    const_reverse_bin_iterator rbegin_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->rbegin_bins();
	    }
	    bin_iterator end_bins()
	    {
		return root_bin->end_bins();
	    }
	    const_bin_iterator end_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->end_bins();
	    }
	    reverse_bin_iterator rend_bins()
	    {
		return root_bin->rend_bins();
	    }
	    const_reverse_bin_iterator rend_bins() const
	    {
		return const_cast<const bin_type*>(root_bin)->rend_bins();
	    }

	    void check()
	    {
		for(bin_iterator it=begin_bins();it!=end_bins();++it)
		{
		    if(!(*it)->is_leaf())
		    {
			if((*it)->volume()!=((*it)->child1->volume()+(*it)->child2->volume()))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"bin \n"<<**it<<"\n has volume unequal to the sum of its childs \n"<<(*(*it)->child1)<<"\n"<<(*(*it)->child2)<<endlog;
			}
			if((*it)->F0!=((*it)->child1->F0+(*it)->child2->F0))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"bin \n"<<**it<<"\n has F0 unequal to the sum of its childs \n"<<(*(*it)->child1)<<"\n"<<(*(*it)->child2)<<endlog;
			}
			if((*it)->w!=((*it)->child1->w+(*it)->child2->w))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"bin \n"<<**it<<"\n has weight unequal to the sum of its childs \n"<<(*(*it)->child1)<<"\n"<<(*(*it)->child2)<<endlog;
			}
		    }
		}
		for(leaf_iterator it=begin_leaves();it!=end_leaves();++it)
		{
		    if(!(*it)->is_leaf())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"leaf iterator \n"<<*it<<"\n is not a leaf bin"<<endlog;
		    }
		}
		value_type w=this->weight();
		evaluate_weight();
		if(w!=this->weight())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"re-evaluated weight "<<this->weight()<<" not equal to original "<<w<<endlog;
		}
	    }

	private:
	    
	    /* Integration domain: */
	    
	    point_type xmin,xmax;

	    /* Root of the binary tree of bins: */

	    bin_type* root_bin;
	    
	    /* Register bin, which has thrown the last point: */
	    
	    bin_type* reg_bin;

	    /* Sub-grid register: */

	    std::set<parni_sub_grid<value_t,1,rng_t,key_t>*>subgrids;

	    /* Private constructors: */
	    /*-----------------------*/

	    /* Default constructor: */

	    parni_generator(size_type n_bins_,grid_modes::type mode_=grid_modes::cumulant_weights):n_bins(n_bins_),mode(mode_),root_bin(NULL),reg_bin(NULL){}

	    /* Constructor with point argument: */

	    parni_generator(point_type* point_,size_type n_bins_,grid_modes::type mode_=grid_modes::cumulant_weights):base_type(point_),n_bins(n_bins_),mode(mode_),root_bin(NULL),reg_bin(NULL){}

	    /* Private methods: */
	    /*------------------*/

	    /* Integrand acceptance function: */

	    static bool accept(const value_type& f)
	    {
		if(f==f and f!=std::numeric_limits<value_type>::infinity() and !(f<(value_type)0))
		{
		    return true;
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"integrand "<<f<<" was not accepted by parni"<<endlog;
		return false;
	    }

	    /* Returns the lowest-level bin containing the argument point: */

	    bin_type* find_leaf(const point_type& x)
	    {
		return root_bin->find_point(x);
	    }

	    /* Returns the lowest-level bin containing the argument point and
	     * sets the second argument to the corresponding left-integral: */

	    bin_type* find_leaf(const point_type& x,value_type& s)
	    {
		s=(value_type)0;
		return root_bin->find_point(x,s);
	    }

	    /* Registers a subgrid: */

	    void register_sub_grid(parni_sub_grid<value_t,1,rng_t,key_t>* subgrid)
	    {
		if(subgrid->grid==this)
		{
		    subgrids.insert(subgrid);
		}
	    }

	    /* Unregisters a subgrid: */

	    void unregister_sub_grid(parni_sub_grid<value_t,1,rng_t,key_t>* subgrid)
	    {
		if(subgrid->grid==this)
		{
		    subgrids.erase(subgrid);
		}
	    }
    };

    /// Parni sub grid generator class, generating points from a parent grid within smaller bounds.

    template<class value_t,std::size_t D,class rng_t,class key_t>class parni_sub_grid: public MC_generator<value_t>
    {
	friend class parni_generator<value_t,D,rng_t,key_t>;

	/* Type definitions: */

	typedef parni_generator<value_t,D,rng_t,key_t> grid_type;
	typedef object_allocator<typename parni_bin<value_t,D,rng_t,key_t>::point_type> base_type;

	public:

	    typedef typename grid_type::value_type value_type;
	    typedef typename grid_type::rn_stream rn_stream;
	    typedef typename grid_type::size_type size_type;
	    typedef typename grid_type::bin_type bin_type;
	    typedef typename grid_type::point_type point_type;

	    /// Destructor.

	    ~parni_sub_grid()
	    {
		if(grid!=NULL)
		{
		    grid->unregister_sub_grid(this);
		}
	    }

	    /// Sets the sub-interval lower bound.
	    
	    bool set_x_min(const point_type& xmin_)
	    {
		if(grid==NULL)
		{
		    return false;
		}
		if(xmin_>=grid->xmin and xmin_<=grid->xmax)
		{
		    xmin=xmin_;
		    return (grid->find_leaf(xmin,smin)!=NULL);
		}
		return false;
	    }

	    /// Returns the sub-interval lower bound.

	    point_type x_min() const
	    {
		return xmin;
	    }

	    /// Sets the sub-interval upper bound.
	    
	    bool set_x_max(const point_type& xmax_)
	    {
		if(grid==NULL)
		{
		    return false;
		}
		if(xmax_>=grid->xmin and xmax_<=grid->xmax)
		{
		    xmax=xmax_;
		    return (grid->find_leaf(xmax,smax)!=NULL);
		}
		return false;
	    }

	    /// Returns the sub-interval upper bound.

	    point_type x_max() const
	    {
		return xmax;
	    }

	    /// Generation method implementation.

	    bool generate()
	    {
		if(grid==NULL)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type rho=rn_stream::throw_number(smin,smax);
		reg_bin=grid->root_bin->find_weight(rho);
		bool q=reg_bin->generate(grid->object(),xmin,xmax);
		if(!q)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->weight()=(smax-smin)*reg_bin->weight();
		grid->reg_bin=reg_bin;
		grid->weight()=grid->norm()*reg_bin->weight();
		return true;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		if(grid==NULL)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!grid->evaluate_weight())
		{
		    this->weight()=0;
		    return false;
		}
		if(grid->object()>=xmin and grid->object()<=xmax)
		{
		    reg_bin=grid->reg_bin;
		    this->weight()=(smax-smin)*reg_bin->weight();
		    return true;
		}
		else
		{
		    this->weight()=0;
		    return false;
		}
	    }

	    /// Refreshes internal parameters depending on the boundary values.

	    void refresh_bounds()
	    {
		if(grid!=NULL)
		{
		    grid->find_leaf(xmin,smin);
		    grid->find_leaf(xmax,smax);
		}
	    }

	private:

	    /// Construct over sub-interval.

	    parni_sub_grid(grid_type* grid_,const point_type& xmin_,const point_type& xmax_):grid(grid_),xmin(xmin_),xmax(xmax_),reg_bin(grid->reg_bin)
	    {
		if(grid!=NULL)
		{
		    grid->register_sub_grid(this);
		}
		refresh_bounds();
	    }

	    /* Parent grid: */

	    grid_type* grid;

	    /* Bounds: */

	    point_type xmin,xmax;

	    /* Register bin: */

	    bin_type* reg_bin;

	    /* Left integral estimates: */

	    value_type smin,smax;
    };
}

#endif /*CAMGEN_PARNI_H_*/

