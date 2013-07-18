//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PARNI_BIN_H_
#define CAMGEN_PARNI_BIN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the parni-bin class. Parni data consists of unequal-sized bins  *
 * with specific weights, stored in the class defined below. The type also       *
 * defines a splitting method and a generation method, throwing uniformly a      *
 * point within the bin.                                                         *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <climits>
#include <iomanip>
#include <sstream>
#include <Camgen/vector.h>
#include <Camgen/utils.h>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/MC_config.h>
#include <Camgen/plt_script.h>
#include <Camgen/rn_strm.h>

namespace Camgen
{
    /* Parni class forward declaration: */

    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class parni;
    
    /* Parni subgrid forward declaration: */
    
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class parni_sub_grid;
    
    /* Parni bin forward declaration: */
    
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class parni_bin;
    
    /* Parni leaf iterators forward declaration: */
    
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class parni_bin_iterator;
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class const_parni_bin_iterator;
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class reverse_parni_bin_iterator;
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class const_reverse_parni_bin_iterator;
    
    /* Parni leaf iterators forward declaration: */
    
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class parni_leaf_iterator;
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class const_parni_leaf_iterator;
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class reverse_parni_leaf_iterator;
    template<class value_t,std::size_t D,class rng_t,class key_t=std::size_t>class const_reverse_parni_leaf_iterator;

    /* Parni bin class: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class parni_bin
    {
	/* Friend declarations: */

	friend class parni<value_t,D,rng_t>;
	friend class parni_sub_grid<value_t,D,rng_t>;
	friend class parni_bin_iterator<value_t,D,rng_t,key_t>;
	friend class const_parni_bin_iterator<value_t,D,rng_t,key_t>;
	friend class reverse_parni_bin_iterator<value_t,D,rng_t,key_t>;
	friend class const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>;
	friend class parni_leaf_iterator<value_t,D,rng_t,key_t>;
	friend class const_parni_leaf_iterator<value_t,D,rng_t,key_t>;
	friend class reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>;
	friend class const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>;

	public:
	    
	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef vector<value_t,D> point_type;
	    typedef key_t key_type;
	    typedef parni_bin_iterator<value_t,D,rng_t,key_t> bin_iterator;
	    typedef const_parni_bin_iterator<value_t,D,rng_t,key_t> const_bin_iterator;
	    typedef reverse_parni_bin_iterator<value_t,D,rng_t,key_t> reverse_bin_iterator;
	    typedef const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t> const_reverse_bin_iterator;
	    typedef parni_leaf_iterator<value_t,D,rng_t,key_t> leaf_iterator;
	    typedef const_parni_leaf_iterator<value_t,D,rng_t,key_t> const_leaf_iterator;
	    typedef reverse_parni_leaf_iterator<value_t,D,rng_t,key_t> reverse_leaf_iterator;
	    typedef const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t> const_reverse_leaf_iterator;

	    /* Static factory methods: */
	    /*-------------------------*/

	    /* Raises the keys by one if zero: */

	    static vector<key_type,D> raise(const vector<key_type,D>& n)
	    {
		vector<key_type,D>b;
		for(size_type i=0;i<D;++i)
		{
		    b[i]=(n[i]<1)?1:n[i];
		}
		return b;
	    }

	    /* Most significant bit readout: */

	    static vector<key_type,D> msb(const vector<key_type,D>& n)
	    {
		vector<key_type,D> a(n),b;
		b.assign(0);
		for(size_type i=0;i<D;++i)
		{
		    while(a[i]>>=1)
		    {
			++b[i];
		    }
		}
		return b;
	    }

	    /* Factory method with input stream argument: */

	    static parni_bin<value_t,D,rng_t,key_t>* create_instance(std::istream& is,const point_type* min_pos_,const point_type* max_pos_,grid_modes::type mode_=grid_modes::cumulant_weights)
	    {
		std::string dummy;
		is>>std::ws;
		if(is.peek()=='N')
		{
		    std::getline(is,dummy);
		    return NULL;
		}
		vector<key_type,D>key_;
		for(size_type i=0;i<D;++i)
		{
		    is>>key_[i];
		}
		parni_bin<value_t,D,rng_t,key_t>* result=new parni_bin<value_t,D,rng_t,key_t>(min_pos_,max_pos_,key_,mode_);
		result->load(is);
		result->child1=create_instance(is,min_pos_,max_pos_,mode_);
		if(result->child1!=NULL)
		{
		    result->child1->parent=result;
		}
		result->child2=create_instance(is,min_pos_,max_pos_,mode_);
		if(result->child2!=NULL)
		{
		    result->child2->parent=result;
		}
		return result;
	    }

	    /* Maximal key value, determines the maximal tree depth: */

	    static const key_type max_depth;

	    /* Public data members: */
	    /*----------------------*/

	    /* Integration interval corners: */

	    const point_type* const min_pos;
	    const point_type* const max_pos;

	    /* Bin key, specifying location and volume: */

	    const vector<key_type,D> key;

	    /* Bin depth, inverse of the volume: */

	    const vector<key_type,D> depth;

	    /* Mode: */

	    const grid_modes::type mode;

	    /* Public constructors/destructors: */
	    /*----------------------------------*/

	    /* Constructor of a rectangular bin with total integration domain
	    lower and upper bound arguments: */

	    parni_bin(const point_type* min_pos_,const point_type* max_pos_,const vector<key_type,D>& key_,grid_modes::type mode_=grid_modes::cumulant_weights):min_pos(min_pos_),max_pos(max_pos_),key(raise(key_)),depth(msb(key)),mode(mode_),F0(0),F1(0),F2(0),fmax(0),fmax1(0),fmax2(0),w(volume()),parent(NULL),child1(NULL),child2(NULL),split_ind(0)
	    {
		value_type x=0,max_edge=edge(0),prev_edge=max_edge;
		bool q=true;
		for(size_type i=1;i<D;++i)
		{
		    x=edge(i);
		    q&=(x==prev_edge);
		    prev_edge=x;
		    if(x>max_edge)
		    {
			split_ind=i;
			max_edge=x;
		    }
		}
		if(q)
		{
		    split_ind=rn_stream::throw_dice(D);
		}
	    }
	    
	    /* Destructor. The bin owns its children: */
	    
	    ~parni_bin()
	    {
		if(child1!=NULL)
		{
		    delete child1;
		}
		if(child2!=NULL)
		{
		    delete child2;
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Splitting method. Returns a boolean denoting whether the split
	     * was successful. */
	    
	    bool split()
	    {
		if(!is_leaf())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"request to split non-leaf bin ignored"<<endlog;
		    return false;
		}
		if(depth[split_ind]==(max_depth-1))
		{
		    return false;
		}
		vector<key_type,D>key1=key;
		(key1[split_ind])<<=1;
		child1=new(std::nothrow)parni_bin<value_t,D,rng_t>(min_pos,max_pos,key1,mode);
		if(child1==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to allocate first child--split aborted"<<endlog;
		    return false;
		}
		child1->parent=this;
		child1->F0=(value_type)0.5*F0;
		child1->F1=(value_type)0.5*F1;

		vector<key_type,D>key2=key;
		(key2[split_ind])<<=1;
		++(key2[split_ind]);
		child2=new(std::nothrow)parni_bin<value_t,D,rng_t>(min_pos,max_pos,key2,mode);
		if(child2==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to allocate second child--split aborted"<<endlog;
		    delete child1;
		    child1=NULL;
		    return false;
		}
		child2->parent=this;
		child2->F0=child1->F0;
		child2->F1=child1->F1;
		if(mode==grid_modes::variance_weights)
		{
		    child1->F2=(value_type)0.5*F2;
		    child2->F2=child1->F2;
		}
		if(mode==grid_modes::maximum_weights)
		{
		    value_type vol=0.5*volume();
		    child1->fmax=fmax1;
		    child1->w=vol*fmax1;
		    child2->fmax=fmax2;
		    child2->w=vol*fmax2;
		    update_weight();
		}
		else
		{
		    child1->w=(value_type)0.5*w;
		    child2->w=child1->w;
		}
		
		return true;
	    }

	    /* Merging method: */

	    void merge()
	    {
		if(child1!=NULL)
		{
		    child1->merge();
		    delete child1;
		    child1=NULL;
		}
		if(child2!=NULL)
		{
		    child2->merge();
		    delete child2;
		    child2=NULL;
		}
	    }

	    /* Resetting method: */

	    void reset()
	    {
		merge();
		F0=0;
		F1=0;
		F2=0;
		fmax=0;
		fmax1=0;
		fmax2=0;
		w=volume();
	    }

	    /* Updates the weight of the bin by the argument, denoting the integrand
	    evaluated at the thrown point: */

	    void update(const point_type& x,const value_type& f)
	    {
		++F0;
		F1+=f;
		if(mode==grid_modes::variance_weights)
		{
		    F2+=(f*f);
		}
		if(mode==grid_modes::maximum_weights)
		{
		    fmax=std::max(f,fmax);
		    if(x[split_ind]<lower_bound(split_ind)+0.5*edge(split_ind))
		    {
			fmax1=std::max(f,fmax1);
		    }
		    else
		    {
			fmax2=std::max(f,fmax2);
		    }
		}
	    }

	    /* Adapts the weight and all subbin weights: */

	    void adapt()
	    {
		if(child1!=NULL)
		{
		    child1->adapt();
		}
		if(child2!=NULL)
		{
		    child2->adapt();
		}
		if(child1==NULL and child2==NULL)
		{
		    value_type vol(volume());
		    switch(mode)
		    {
			case grid_modes::cumulant_weights:
			    w=(F0==(value_type)0)?vol:(vol*F1/F0);
			    break;
			case grid_modes::variance_weights:
			    w=(F0==(value_type)0)?vol:(vol*std::sqrt(F2/F0));
			    break;
			case grid_modes::maximum_weights:
			    w=vol*fmax;
			    break;
			default:
			    w=(F0==(value_type)0)?vol:(vol*F1/F0);
		    }
		}
		if(child1!=NULL and child2!=NULL)
		{
		    F0=child1->F0+child2->F0;
		    F1=child1->F1+child2->F1;
		    if(mode==grid_modes::variance_weights)
		    {
			F2=child1->F2+child2->F2;
		    }
		    if(mode==grid_modes::maximum_weights)
		    {
			fmax1=child1->fmax;
			fmax2=child2->fmax;
			fmax=std::max(fmax1,fmax2);
		    }
		    w=child1->w+child2->w;
		}
	    }

	    /* Generates uniformly a point within the bin: */

	    bool generate(point_type& x)
	    {
		x=lower_bound();
		for(size_type i=0;i<D;++i)
		{
		    x[i]+=rn_stream::throw_number(edge(i));
		}
		return true;
	    }

	    /* Generates uniformly a point within the bin and within argument
	     * bounds: */

	    bool generate(point_type& x,const point_type& xmin,const point_type& xmax)
	    {
		point_type lb=lower_bound();
		point_type ub=upper_bound();
		value_type xminbis,xmaxbis;
		for(size_type i=0;i<D;++i)
		{
		    xminbis=std::max(xmin[i],lb[i]);
		    xmaxbis=std::min(xmax[i],ub[i]);
		    if(xminbis>xmaxbis)
		    {
			return false;
		    }
		    x[i]=rn_stream::throw_number(xminbis,xmaxbis);
		}
		return true;
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /* Total grid volume: */

	    value_type grid_volume() const
	    {
		value_type v((*max_pos)[0]-(*min_pos)[0]);
		for(size_type i=1;i<D;++i)
		{
		    v*=((*max_pos)[i]-(*min_pos)[i]);
		}
		return v;
	    }

	    /* Bin volume: */

	    value_type volume() const
	    {
		value_type v(((*max_pos)[0]-(*min_pos)[0])/(1<<depth[0]));
		for(size_type i=1;i<D;++i)
		{
		    v*=(((*max_pos)[i]-(*min_pos)[i])/(1<<depth[i]));
		}
		return v;
	    }

	    /* i-th edge length: */

	    value_type edge(size_type i) const
	    {
		return ((*max_pos)[i]-(*min_pos)[i])/(1<<depth[i]);
	    }

	    /* Bin lower bound: */

	    point_type lower_bound() const
	    {
		point_type lb(*min_pos);
		value_type denom;
		for(size_type i=0;i<D;++i)
		{
		    denom=(value_type)1/(1<<depth[i]);
		    lb[i]+=((key[i]*denom-1)*((*max_pos)[i]-(*min_pos)[i]));
		}
		return lb;
	    }

	    /* i-th lower bound direction: */

	    value_type lower_bound(size_type i) const
	    {
		value_type denom=(value_type)1/(1<<depth[i]);
		return ((key[i]*denom-1)*((*max_pos)[i]-(*min_pos)[i])+(*min_pos)[i]);
	    }

	    /* Bin upper bound: */

	    point_type upper_bound() const
	    {
		point_type ub(*min_pos);
		value_type denom;
		for(size_type i=0;i<D;++i)
		{
		    denom=(value_type)1/(1<<depth[i]);
		    ub[i]+=(((key[i]+1)*denom-1)*((*max_pos)[i]-(*min_pos)[i]));
		}
		return ub;
	    }

	    /* i-th upper bound direction: */

	    value_type upper_bound(size_type i) const
	    {
		value_type denom=(value_type)1/(1<<depth[i]);
		return (((key[i]+1)*denom-1)*((*max_pos)[i]-(*min_pos)[i])+(*min_pos)[i]);
	    }

	    /* Returns the splitted index: */

	    size_type split_index() const
	    {
		return split_ind;
	    }

	    /* Returns the bin weight: */

	    value_type weight() const
	    {
		return (w==(value_type)0)?1:(volume()/w);
	    }

	    /* Returns whether there is no parent: */

	    bool is_root() const
	    {
		return (parent==NULL);
	    }

	    /* Returns whether there are no children: */

	    bool is_leaf() const
	    {
		return (child1==NULL and child2==NULL);
	    }

	    /* Returns whether the instance is the first child of the parent: */

	    bool is_first() const
	    {
		return (parent==NULL)?false:(this==parent->child1);
	    }

	    /* Returns whether the instance is the second child of the parent: */

	    bool is_second() const
	    {
		return (parent==NULL)?false:(this==parent->child2);
	    }

	    /* Returns the number of descendants: */

	    std::size_t count_bins() const
	    {
		return is_leaf()?1:(1+child1->count_bins()+child2->count_bins());
	    }

	    /* Returns the number of descendent leaves: */

	    std::size_t count_leaves() const
	    {
		return is_leaf()?1:(child1->count_leaves()+child2->count_leaves());
	    }

	    /* Returns whether the argument is an ancester of the instance: */

	    bool is_child_of(const parni_bin<value_t,D,rng_t>* other) const
	    {
		parni_bin<value_type,D,rng_t>* bin=parent;
		while(bin!=NULL and bin!=other and bin!=this)
		{
		    bin=bin->parent;
		}
		return (bin!=NULL and bin==other and bin!=this);
	    }

	    /* Returns whther the argument is a descendant of the instance: */

	    bool is_parent_of(const parni_bin<value_t,D,rng_t>* other) const
	    {
		return other->is_child_of(this);
	    }

	    /* Finds the lowest-level sub-bin containing the point x: */

	    parni_bin<value_t,D,rng_t>* find_point(const point_type& x)
	    {
		value_type lb(lower_bound(split_ind)),vol(edge(split_ind));
                if(x[split_ind]<lb or x[split_ind]>(lb+vol))
                {
                    log(log_level::warning)<<CAMGEN_STREAMLOC<<"point "<<x<<" not found in parni bin\n"<<*this<<"\n--returning NULL."<<endlog;
                    return NULL;
                }
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (x[split_ind]<(lb+0.5*vol))?(child1->find_point(x)):(child2->find_point(x));
	    }

	    /* Finds the lowest-level const sub-bin containing the point x: */

	    const parni_bin<value_t,D,rng_t>* find_point(const point_type& x) const
	    {
		value_type lb(lower_bound(split_ind)),vol(edge(split_ind));
                if(x[split_ind]<lb or x[split_ind]>(lb+vol))
                {
                    log(log_level::warning)<<CAMGEN_STREAMLOC<<"point "<<x<<" not found in parni bin\n"<<*this<<"\n--returning NULL."<<endlog;
                    return NULL;
                }
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (x[split_ind]<(lb+0.5*vol))?(child1->find_point(x)):(child2->find_point(x));
	    }

	    /* Finds the lowest-level sub-bin with weight smaller than the argument: */

	    parni_bin<value_t,D,rng_t>* find_weight(const value_type& w_)
	    {
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (w_<child1->w)?(child1->find_weight(w_)):(child2->find_weight(w_-child1->w));
	    }

	    /* Finds the lowest-level const sub-bin with weight smaller than the argument: */

	    const parni_bin<value_t,D,rng_t>* find_weight(const value_type& w_) const
	    {
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (w_<child1->w)?(child1->find_weight(w_)):(child2->find_weight(w_-child1->w));
	    }

	    /* Recursive printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		std::stringstream ss;
		for(std::size_t i=0;i<D-1;++i)
		{
		    ss<<"["<<lower_bound(i)<<","<<upper_bound(i)<<"] x ";
		}
		ss<<'['<<lower_bound(D-1)<<','<<upper_bound(D-1)<<']';
		os<<key<<std::setw(20*D)<<depth<<std::setw(46*D)<<ss.str()<<std::setw(20)<<w/volume();
		return os;
	    }

	    /* Printing helper function: */

	    plot_script* print_rectangle(plot_script* ps) const
	    {
		if(D==2)
		{
		    if(is_leaf())
		    {
			ps->add_object(new plot_rectangle(lower_bound(0),lower_bound(1),upper_bound(0),upper_bound(1)));
		    }
		}
		return ps;
	    }

	    /* Iterator definitions: */

	    leaf_iterator begin_leaves()
	    {
		return leaf_iterator(first_leaf());
	    }
	    const_leaf_iterator begin_leaves() const
	    {
		return const_leaf_iterator(first_leaf());
	    }
	    reverse_leaf_iterator rbegin_leaves()
	    {
		return reverse_leaf_iterator(last_leaf());
	    }
	    const_reverse_leaf_iterator rbegin_leaves() const
	    {
		return const_reverse_leaf_iterator(last_leaf());
	    }
	    leaf_iterator end_leaves()
	    {
		return leaf_iterator(NULL);
	    }
	    const_leaf_iterator end_leaves() const
	    {
		return const_leaf_iterator(NULL);
	    }
	    reverse_leaf_iterator rend_leaves()
	    {
		return reverse_leaf_iterator(NULL);
	    }
	    const_reverse_leaf_iterator rend_leaves() const
	    {
		return const_reverse_leaf_iterator(NULL);
	    }
	    bin_iterator begin_bins()
	    {
		return bin_iterator(first_bin());
	    }
	    const_bin_iterator begin_bins() const
	    {
		return const_bin_iterator(first_bin());
	    }
	    reverse_bin_iterator rbegin_bins()
	    {
		return reverse_bin_iterator(last_bin());
	    }
	    const_reverse_bin_iterator rbegin_bins() const
	    {
		return const_reverse_bin_iterator(last_bin());
	    }
	    bin_iterator end_bins()
	    {
		return bin_iterator(NULL);
	    }
	    const_bin_iterator end_bins() const
	    {
		return const_bin_iterator(NULL);
	    }
	    reverse_bin_iterator rend_bins()
	    {
		return reverse_bin_iterator(NULL);
	    }
	    const_reverse_bin_iterator rend_bins() const
	    {
		return const_reverse_bin_iterator(NULL);
	    }

	private:

	    /* Bin volume and weight counter variables: */
	    
	    value_type F0,F1,F2,fmax,fmax1,fmax2,w;

	    /* Bin parent: */

	    parni_bin<value_t,D,rng_t>* parent;
	    
	    /* Sub-bins: */
	    
	    parni_bin<value_t,D,rng_t>* child1;
	    parni_bin<value_t,D,rng_t>* child2;

	    /* division index: */

	    size_type split_ind;

	    /* Loads bin data from input stream: */

	    std::istream& load(std::istream& is)
	    {
		safe_read(is,F0);
		safe_read(is,F1);
		if(mode==grid_modes::variance_weights)
		{
		    safe_read(is,F2);
		}
		if(mode==grid_modes::maximum_weights)
		{
		    safe_read(is,fmax);
		    safe_read(is,fmax1);
		    safe_read(is,fmax2);
		}
		safe_read(is,w);
		return is;
	    }

	    /* Recursively saves bin data to output stream: */

	    std::ostream& save(std::ostream& os) const
	    {
		for(size_type i=0;i<D;++i)
		{
		    os<<key[i]<<"\t";
		}
		safe_write(os,F0);
		os<<"\t";
		safe_write(os,F1);
		if(mode==grid_modes::variance_weights)
		{
		    os<<"\t";
		    safe_write(os,F2);
		}
		if(mode==grid_modes::maximum_weights)
		{
		    os<<"\t";
		    safe_write(os,fmax);
		    os<<"\t";
		    safe_write(os,fmax1);
		    os<<"\t";
		    safe_write(os,fmax2);
		}
		os<<"\t";
		safe_write(os,w);
		os<<std::endl;
		if(child1!=NULL)
		{
		    child1->save(os);
		}
		else
		{
		    os<<'N'<<std::endl;
		}
		if(child2!=NULL)
		{
		    child2->save(os);
		}
		else
		{
		    os<<'N'<<std::endl;
		}
		return os;
	    }

	    /* Returns the current instance: */

	    parni_bin<value_t,D,rng_t>* first_bin()
	    {
		return this;
	    }

	    /* Returns the current const instance: */

	    const parni_bin<value_t,D,rng_t>* first_bin() const
	    {
		return this;
	    }
	    
	    /* Returns the left-most leaf of the tree with the current bin as
	     * root bin: */

	    parni_bin<value_t,D,rng_t>* first_leaf()
	    {
		parni_bin<value_t,D,rng_t>* b=this;
		while(b->child1!=NULL)
		{
		    b=b->child1;
		}
		return b;
	    }

	    /* Returns the left-most leaf of the tree with the current bin as
	     * root bin: */

	    const parni_bin<value_t,D,rng_t>* first_leaf() const
	    {
		const parni_bin<value_t,D,rng_t>* b=this;
		while(b->child1!=NULL)
		{
		    b=b->child1;
		}
		return b;
	    }

	    /* Returns the next bin in the left-up-right traversal method: */

	    parni_bin<value_t,D,rng_t>* next_bin()
	    {
		if(!is_leaf())
		{
		    return child1;
		}
		if(is_first())
		{
		    return parent->child2;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2);
	    }

	    /* Returns the const next bin in the left-up-right traversal method: */

	    const parni_bin<value_t,D,rng_t>* next_bin() const
	    {
		if(!is_leaf())
		{
		    return child1;
		}
		if(is_first())
		{
		    return parent->child2;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2);
	    }

	    /* Returns the right-next leaf. If the current instance is not a
	     * leaf, returns the left-most leaf. */

	    parni_bin<value_t,D,rng_t>* next_leaf()
	    {
		if(!is_leaf())
		{
		    return first_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child1)
		{
		    return parent->child2->first_leaf();
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2->first_leaf());
	    }

	    /* Returns the right-next const leaf. If the current instance is not a
	     * leaf, returns the left-most leaf. */

	    const parni_bin<value_t,D,rng_t>* next_leaf() const
	    {
		if(!is_leaf())
		{
		    return first_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child1)
		{
		    return parent->child2->first_leaf();
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2->first_leaf());
	    }

	    /* Returns the current instance: */

	    parni_bin<value_t,D,rng_t>* last_bin()
	    {
		return this;
	    }

	    /* Returns the current const instance: */

	    const parni_bin<value_t,D,rng_t>* last_bin() const
	    {
		return this;
	    }

	    /* Returns the right-most leaf of the binary tree with the current
	     * bin as root: */

	    parni_bin<value_t,D,rng_t>* last_leaf()
	    {
		parni_bin<value_t,D,rng_t>* b=this;
		while(b->child2!=NULL)
		{
		    b=b->child2;
		}
		return b;
	    }
	    const parni_bin<value_t,D,rng_t>* last_leaf() const
	    {
		parni_bin<value_t,D,rng_t>* b=this;
		while(b->child2!=NULL)
		{
		    b=b->child2;
		}
		return b;
	    }

	    /* Returns the previous bin in the left-up-right traversal method: */

	    parni_bin<value_t,D,rng_t>* previous_bin()
	    {
		if(!is_leaf())
		{
		    return child2;
		}
		if(is_second())
		{
		    return parent->child1;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1);
	    }

	    /* Returns the const previous bin in the left-up-right traversal method: */

	    const parni_bin<value_t,D,rng_t>* previous_bin() const
	    {
		if(!is_leaf())
		{
		    return child2;
		}
		if(is_second())
		{
		    return parent->child1;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1);
	    }

	    /* Returns the previous leaf: */

	    parni_bin<value_t,D,rng_t>* previous_leaf()
	    {
		if(!is_leaf())
		{
		    return last_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child2)
		{
		    return parent->child1->last_leaf();
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1->last_leaf());
	    }

	    /* Returns the previous const leaf: */

	    const parni_bin<value_t,D,rng_t>* previous_leaf() const
	    {
		if(!is_leaf())
		{
		    return last_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child2)
		{
		    return parent->child1->last_leaf();
		}
		parni_bin<value_t,D,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1->last_leaf());
	    }

	    /* Updates the weights up the tree with the sum rule: */

	    void update_weight()
	    {
		if(child1!=NULL and child2!=NULL)
		{
		    w=child1->w+child2->w;
		}
		if(parent!=NULL)
		{
		    parent->update_weight();
		}
	    }
    };
    template<class value_t,std::size_t D,class rng_t,class key_type>const key_type parni_bin<value_t,D,rng_t,key_type>::max_depth=sizeof(key_type)*CHAR_BIT-1;

    /* Parni bin class template specialisation for 1-dimensional domains: */

    template<class value_t,class rng_t,class key_t>class parni_bin<value_t,1,rng_t,key_t>
    {
	/* Friend declarations: */

	friend class parni<value_t,1,rng_t>;
	friend class parni_sub_grid<value_t,1,rng_t>;
	friend class parni_bin_iterator<value_t,1,rng_t,key_t>;
	friend class const_parni_bin_iterator<value_t,1,rng_t,key_t>;
	friend class reverse_parni_bin_iterator<value_t,1,rng_t,key_t>;
	friend class const_reverse_parni_bin_iterator<value_t,1,rng_t,key_t>;
	friend class parni_leaf_iterator<value_t,1,rng_t,key_t>;
	friend class const_parni_leaf_iterator<value_t,1,rng_t,key_t>;
	friend class reverse_parni_leaf_iterator<value_t,1,rng_t,key_t>;
	friend class const_reverse_parni_leaf_iterator<value_t,1,rng_t,key_t>;

	public:
	    
	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef value_t point_type;
	    typedef key_t key_type;
	    typedef parni_bin_iterator<value_t,1,rng_t,key_t> bin_iterator;
	    typedef const_parni_bin_iterator<value_t,1,rng_t,key_t> const_bin_iterator;
	    typedef reverse_parni_bin_iterator<value_t,1,rng_t,key_t> reverse_bin_iterator;
	    typedef const_reverse_parni_bin_iterator<value_t,1,rng_t,key_t> const_reverse_bin_iterator;
	    typedef parni_leaf_iterator<value_t,1,rng_t,key_t> leaf_iterator;
	    typedef const_parni_leaf_iterator<value_t,1,rng_t,key_t> const_leaf_iterator;
	    typedef reverse_parni_leaf_iterator<value_t,1,rng_t,key_t> reverse_leaf_iterator;
	    typedef const_reverse_parni_leaf_iterator<value_t,1,rng_t,key_t> const_reverse_leaf_iterator;

	    /* Static factory methods: */
	    /*-------------------------*/

	    /* Most significant bit readout: */

	    static key_type msb(key_type n)
	    {
		key_type a(n),b(0);
		while(a>>=1)
		{
		    ++b;
		}
		return b;
	    }

	    /* Factory method with input stream argument: */

	    static parni_bin<value_t,1,rng_t,key_t>* create_instance(std::istream& is,const point_type* min_pos_,const point_type* max_pos_,grid_modes::type mode_=grid_modes::cumulant_weights)
	    {
		std::string dummy;
		is>>std::ws;
		if(is.peek()=='N')
		{
		    std::getline(is,dummy);
		    return NULL;
		}
		key_type key_(1);
		is>>key_;
		parni_bin<value_t,1,rng_t,key_t>* result=new parni_bin<value_t,1,rng_t,key_t>(min_pos_,max_pos_,key_,mode_);
		result->load(is);
		result->child1=create_instance(is,min_pos_,max_pos_,mode_);
		if(result->child1!=NULL)
		{
		    result->child1->parent=result;
		}
		result->child2=create_instance(is,min_pos_,max_pos_,mode_);
		if(result->child2!=NULL)
		{
		    result->child2->parent=result;
		}
		return result;
	    }

	    /* Maximal key value, determines the maximal tree depth: */

	    static const key_type max_depth;

	    /* Public data members: */
	    /*----------------------*/

	    /* Integration interval corners: */

	    const point_type* const min_pos;
	    const point_type* const max_pos;

	    /* Bin key, specifying location and volume: */

	    const key_type key;

	    /* Bin depth, inverse of the volume: */

	    const key_type depth;

	    /* Grid mode: */

	    const grid_modes::type mode;

	    /* Public constructors/destructors: */
	    /*----------------------------------*/

	    /* Constructor of a rectangular bin with total integration domain
	    lower and upper bound arguments: */

	    parni_bin(const point_type* min_pos_,const point_type* max_pos_,key_type key_,grid_modes::type mode_=grid_modes::cumulant_weights):min_pos(min_pos_),max_pos(max_pos_),key(std::max((key_type)1,key_)),depth(msb(key)),mode(mode_),F0(0),F1(0),F2(0),fmax(0),fmax1(0),fmax2(0),w(volume()),parent(NULL),child1(NULL),child2(NULL){}
	    
	    /* Destructor. The bin owns its children: */
	    
	    ~parni_bin()
	    {
		if(child1!=NULL)
		{
		    delete child1;
		}
		if(child2!=NULL)
		{
		    delete child2;
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Splitting method. Returns a boolean denoting whether the split
	     * was successful. */
	    
	    bool split()
	    {
		if(!is_leaf())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"request to split non-leaf bin ignored"<<endlog;
		    return false;
		}
		if(depth==(max_depth-1))
		{
		    return false;
		}
		child1=new(std::nothrow)parni_bin<value_t,1,rng_t>(min_pos,max_pos,(key<<1),mode);
		if(child1==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to allocate first child--split aborted."<<endlog;
		    return false;
		}
		child1->parent=this;
		child1->F0=(value_type)0.5*F0;
		child1->F1=(value_type)0.5*F1;

		child2=new(std::nothrow)parni_bin<value_t,1,rng_t>(min_pos,max_pos,(key<<1)+1,mode);
		if(child2==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to allocate second child--split aborted."<<endlog;
		    delete child1;
		    child1=NULL;
		    return false;
		}
		child2->parent=this;
		child2->F0=child1->F0;
		child2->F1=child1->F1;

		if(mode==grid_modes::variance_weights)
		{
		    child1->F2=(value_type)0.5*F2;
		    child2->F2=child1->F2;
		}
		if(mode==grid_modes::maximum_weights)
		{
		    value_type vol=0.5*volume();
		    child1->fmax=fmax1;
		    child1->w=vol*fmax1;
		    child2->fmax=fmax2;
		    child2->w=vol*fmax2;
		    update_weight();
		}
		else
		{
		    child1->w=(value_type)0.5*w;
		    child2->w=child1->w;
		}

		return true;
	    }

	    /* Merging method: */

	    void merge()
	    {
		if(child1!=NULL)
		{
		    child1->merge();
		    delete child1;
		    child1=NULL;
		}
		if(child2!=NULL)
		{
		    child2->merge();
		    delete child2;
		    child2=NULL;
		}
	    }

	    /* Resetting method: */

	    void reset()
	    {
		merge();
		F0=0;
		F1=0;
		F2=0;
		fmax=0;
		fmax1=0;
		fmax2=0;
		w=volume();
	    }

	    /* Updates the weight of the bin by the argument, denoting the integrand
	    evaluated at the thrown point: */

	    void update(const point_type& x,const value_type& f)
	    {
		++F0;
		F1+=f;
		if(mode==grid_modes::variance_weights)
		{
		    F2+=(f*f);
		}
		value_type vol(volume());
		if(mode==grid_modes::maximum_weights)
		{
		    fmax=std::max(f,fmax);
		    if(x<lower_bound()+0.5*vol)
		    {
			fmax1=std::max(f,fmax1);
		    }
		    else
		    {
			fmax2=std::max(f,fmax2);
		    }
		}
	    }

	    /* Adapts the weight and all subbin weights: */

	    void adapt()
	    {
		if(child1!=NULL)
		{
		    child1->adapt();
		}
		if(child2!=NULL)
		{
		    child2->adapt();
		}
		if(child1==NULL and child2==NULL)
		{
		    value_type vol(volume());
		    if(F0==(value_type)0)
		    {
			w=vol;
		    }
		    else
		    {
			switch(mode)
			{
			    case grid_modes::cumulant_weights:
				w=vol*F1/F0;
				break;
			    case grid_modes::variance_weights:
				w=vol*std::sqrt(F2/F0);
				break;
			    case grid_modes::maximum_weights:
				w=vol*fmax;
				break;
			    default:
				w=vol;
			}
		    }
		}
		if(child1!=NULL and child2!=NULL)
		{
		    F0=child1->F0+child2->F0;
		    F1=child1->F1+child2->F1;
		    if(mode==grid_modes::variance_weights)
		    {
			F2=child1->F2+child2->F2;
		    }
		    if(mode==grid_modes::maximum_weights)
		    {
			fmax1=child1->fmax;
			fmax2=child2->fmax;
			fmax=std::max(fmax1,fmax2);
		    }
		    w=child1->w+child2->w;
		}
	    }

	    /* Generates uniformly a point within the bin and calls add_point():
	     * */

	    bool generate(point_type& x)
	    {
		x=lower_bound()+rn_stream::throw_number(0,volume());
		return true;
	    }

	    /* Generates uniformly a point within the bin within the boundaries
	     * and calls add_point() if succesful: */

	    bool generate(point_type& x,const point_type& xmin,const point_type& xmax)
	    {
		value_type xminbis=std::max(xmin,lower_bound());
		value_type xmaxbis=std::min(xmax,upper_bound());
		if(xminbis>xmaxbis)
		{
		    return false;
		}
		x=rn_stream::throw_number(xminbis,xmaxbis);
		return true;
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /* Total grid volume: */

	    value_type grid_volume() const
	    {
		return (*max_pos-*min_pos);
	    }

	    /* Bin volume: */

	    value_type volume() const
	    {
		return grid_volume()/(1<<depth);
	    }

	    /* Bin lower bound: */

	    value_type lower_bound() const
	    {
		return (((value_type)key)/(1<<depth)-1)*grid_volume()+*min_pos;
	    }

	    /* Bin upper bound: */

	    value_type upper_bound() const
	    {
		return (((value_type)(key+1))/(1<<depth)-1)*grid_volume()+*min_pos;
	    }

	    /* Returns the bin weight: */

	    value_type weight() const
	    {
		return (w==(value_type)0)?1:(volume()/w);
	    }

	    /* Returns whether there is no parent: */

	    bool is_root() const
	    {
		return (parent==NULL);
	    }

	    /* Returns whether there are no children: */

	    bool is_leaf() const
	    {
		return (child1==NULL and child2==NULL);
	    }

	    /* Returns whether the instance is the first child of the parent: */

	    bool is_first() const
	    {
		return (parent==NULL)?false:(this==parent->child1);
	    }

	    /* Returns whether the instance is the second child of the parent: */

	    bool is_second() const
	    {
		return (parent==NULL)?false:(this==parent->child2);
	    }

	    /* Returns the number of descendants: */

	    std::size_t count_bins() const
	    {
		return is_leaf()?1:(1+child1->count_bins()+child2->count_bins());
	    }

	    /* Returns the number of descendent leaves: */

	    std::size_t count_leaves() const
	    {
		return is_leaf()?1:(child1->count_leaves()+child2->count_leaves());
	    }

	    /* Returns whether the argument is an ancester of the instance: */

	    bool is_child_of(const parni_bin<value_t,1,rng_t>* other) const
	    {
		parni_bin<value_type,1,rng_t>* bin=parent;
		while(bin!=NULL and bin!=other and bin!=this)
		{
		    bin=bin->parent;
		}
		return (bin!=NULL and bin==other and bin!=this);
	    }

	    /* Returns whther the argument is a descendant of the instance: */

	    bool is_parent_of(const parni_bin<value_t,1,rng_t>* other) const
	    {
		return other->is_child_of(this);
	    }

	    /* Finds the lowest-level sub-bin containing the point x: */

	    parni_bin<value_t,1,rng_t>* find_point(const point_type& x)
	    {
		value_type lb(lower_bound()),vol(volume());
                if(x<lb or x>(lb+vol))
                {
                    log(log_level::warning)<<CAMGEN_STREAMLOC<<"point "<<x<<" not found in parni bin\n"<<*this<<"\n--returning NULL."<<endlog;
                    return NULL;
                }
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (x<(lb+0.5*vol))?(child1->find_point(x)):(child2->find_point(x));
	    }

	    /* Finds the lowest-level const sub-bin containing the point x: */

	    const parni_bin<value_t,1,rng_t>* find_point(const point_type& x) const
	    {
		value_type lb(lower_bound()),vol(volume());
                if(x<lb or x>(lb+vol))
                {
                    log(log_level::warning)<<CAMGEN_STREAMLOC<<"point "<<x<<" not found in parni bin\n"<<*this<<"\n--returning NULL."<<endlog;
                    return NULL;
                }
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (x<(lb+0.5*vol))?(child1->find_point(x)):(child2->find_point(x));
	    }

	    /* Finds the lowest-level sub-bin containing the point x and store the
	     * corresponding left integral in s: */

	    parni_bin<value_t,1,rng_t>* find_point(const point_type& x,value_type& s)
	    {
		value_type lb(lower_bound()),vol(volume());
                if(x<lb or x>(lb+vol))
                {
                    log(log_level::warning)<<CAMGEN_STREAMLOC<<"point "<<x<<" not found in parni bin\n"<<*this<<"\n--returning NULL."<<endlog;
		    return NULL;
                }
		if(child1==NULL or child2==NULL)
		{
		    s+=((x-lb)*w/vol);
		    return this;
		}
		return (x<(lb+0.5*vol))?(child1->find_point(x,s)):(child2->find_point(x,s+=(child1->w)));
	    }

	    /* Finds the lowest-level const sub-bin containing the point x and store the
	     * corresponding left integral in s: */

	    const parni_bin<value_t,1,rng_t>* find_point(const point_type& x,value_type& s) const
	    {
		value_type lb(lower_bound()),vol(volume());
                if(x<lb or x>(lb+vol))
                {
                    log(log_level::warning)<<CAMGEN_STREAMLOC<<"point "<<x<<" not found in parni grid--returning NULL."<<endlog;
                    return NULL;
                }
		if(child1==NULL or child2==NULL)
		{
		    s+=((x-lb)*w/vol);
		    return this;
		}
		return (x<(lb+0.5*vol))?(child1->find_point(x,s)):(child2->find_point(x,s+=(child1->w)));
	    }

	    /* Finds the lowest-level sub-bin with weight smaller than the argument: */

	    parni_bin<value_t,1,rng_t>* find_weight(const value_type& w_)
	    {
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (w_<child1->w)?(child1->find_weight(w_)):(child2->find_weight(w_-child1->w));
	    }

	    /* Finds the lowest-level const sub-bin with weight smaller than the argument: */

	    const parni_bin<value_t,1,rng_t>* find_weight(const value_type& w_) const
	    {
		if(child1==NULL or child2==NULL)
		{
		    return this;
		}
		return (w_<child1->w)?(child1->find_weight(w_)):(child2->find_weight(w_-child1->w));
	    }

	    /* Recursive printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		std::stringstream ss;
		ss<<'['<<lower_bound()<<','<<upper_bound()<<']';
		os<<key<<std::setw(20)<<depth<<std::setw(20)<<ss.str()<<std::setw(20)<<w/volume();
		return os;
	    }

	    /* Printing helper function: */

	    plot_script* print_rectangle(plot_script* ps) const
	    {
		if(is_leaf())
		{
		    ps->add_object(new plot_rectangle(lower_bound(),0,upper_bound(),w/volume()));
		}
		return ps;
	    }

	    /* Iterator definitions: */

	    leaf_iterator begin_leaves()
	    {
		return leaf_iterator(first_leaf());
	    }
	    const_leaf_iterator begin_leaves() const
	    {
		return const_leaf_iterator(first_leaf());
	    }
	    reverse_leaf_iterator rbegin_leaves()
	    {
		return reverse_leaf_iterator(last_leaf());
	    }
	    const_reverse_leaf_iterator rbegin_leaves() const
	    {
		return const_reverse_leaf_iterator(last_leaf());
	    }
	    leaf_iterator end_leaves()
	    {
		return leaf_iterator(NULL);
	    }
	    const_leaf_iterator end_leaves() const
	    {
		return const_leaf_iterator(NULL);
	    }
	    reverse_leaf_iterator rend_leaves()
	    {
		return reverse_leaf_iterator(NULL);
	    }
	    const_reverse_leaf_iterator rend_leaves() const
	    {
		return const_reverse_leaf_iterator(NULL);
	    }
	    bin_iterator begin_bins()
	    {
		return bin_iterator(first_bin());
	    }
	    const_bin_iterator begin_bins() const
	    {
		return const_bin_iterator(first_bin());
	    }
	    reverse_bin_iterator rbegin_bins()
	    {
		return reverse_bin_iterator(last_bin());
	    }
	    const_reverse_bin_iterator rbegin_bins() const
	    {
		return const_reverse_bin_iterator(last_bin());
	    }
	    bin_iterator end_bins()
	    {
		return bin_iterator(NULL);
	    }
	    const_bin_iterator end_bins() const
	    {
		return const_bin_iterator(NULL);
	    }
	    reverse_bin_iterator rend_bins()
	    {
		return reverse_bin_iterator(NULL);
	    }
	    const_reverse_bin_iterator rend_bins() const
	    {
		return const_reverse_bin_iterator(NULL);
	    }

	private:

	    /* Bin volume and weight counter variables: */
	    
	    value_type F0,F1,F2,fmax,fmax1,fmax2,w;

	    /* Bin parent: */

	    parni_bin<value_t,1,rng_t>* parent;
	    
	    /* Sub-bins: */
	    
	    parni_bin<value_t,1,rng_t>* child1;
	    parni_bin<value_t,1,rng_t>* child2;

	    /* Loads bin data from input stream: */

	    std::istream& load(std::istream& is)
	    {
		safe_read(is,F0);
		safe_read(is,F1);
		if(mode==grid_modes::variance_weights)
		{
		    safe_read(is,F2);
		}
		if(mode==grid_modes::maximum_weights)
		{
		    safe_read(is,fmax);
		    safe_read(is,fmax1);
		    safe_read(is,fmax2);
		}
		safe_read(is,w);
		return is;
	    }

	    /* Recursively saves bin data to output stream: */

	    std::ostream& save(std::ostream& os) const
	    {
		os<<key<<"\t";
		safe_write(os,F0);
		os<<"\t";
		safe_write(os,F1);
		if(mode==grid_modes::variance_weights)
		{
		    os<<"\t";
		    safe_write(os,F2);
		}
		if(mode==grid_modes::maximum_weights)
		{
		    os<<"\t";
		    safe_write(os,fmax);
		    os<<"\t";
		    safe_write(os,fmax1);
		    os<<"\t";
		    safe_write(os,fmax2);
		}
		os<<"\t";
		safe_write(os,w);
		os<<std::endl;
		if(child1!=NULL)
		{
		    child1->save(os);
		}
		else
		{
		    os<<'N'<<std::endl;
		}
		if(child2!=NULL)
		{
		    child2->save(os);
		}
		else
		{
		    os<<'N'<<std::endl;
		}
		return os;
	    }

	    /* Returns the current instance: */

	    parni_bin<value_t,1,rng_t>* first_bin()
	    {
		return this;
	    }

	    /* Returns the current const instance: */

	    const parni_bin<value_t,1,rng_t>* first_bin() const
	    {
		return this;
	    }
	    
	    /* Returns the left-most leaf of the tree with the current bin as
	     * root bin: */

	    parni_bin<value_t,1,rng_t>* first_leaf()
	    {
		parni_bin<value_t,1,rng_t>* b=this;
		while(b->child1!=NULL)
		{
		    b=b->child1;
		}
		return b;
	    }

	    /* Returns the left-most leaf of the tree with the current bin as
	     * root bin: */

	    const parni_bin<value_t,1,rng_t>* first_leaf() const
	    {
		const parni_bin<value_t,1,rng_t>* b=this;
		while(b->child1!=NULL)
		{
		    b=b->child1;
		}
		return b;
	    }

	    /* Returns the next bin in the left-up-right traversal method: */

	    parni_bin<value_t,1,rng_t>* next_bin()
	    {
		if(!is_leaf())
		{
		    return child1;
		}
		if(is_first())
		{
		    return parent->child2;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2);
	    }

	    /* Returns the const next bin in the left-up-right traversal method: */

	    const parni_bin<value_t,1,rng_t>* next_bin() const
	    {
		if(!is_leaf())
		{
		    return child1;
		}
		if(is_first())
		{
		    return parent->child2;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2);
	    }

	    /* Returns the right-next leaf. If the current instance is not a
	     * leaf, returns the left-most leaf. */

	    parni_bin<value_t,1,rng_t>* next_leaf()
	    {
		if(!is_leaf())
		{
		    return first_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child1)
		{
		    return parent->child2->first_leaf();
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2->first_leaf());
	    }

	    /* Returns the right-next const leaf. If the current instance is not a
	     * leaf, returns the left-most leaf. */

	    const parni_bin<value_t,1,rng_t>* next_leaf() const
	    {
		if(!is_leaf())
		{
		    return first_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child1)
		{
		    return parent->child2->first_leaf();
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_second())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child2->first_leaf());
	    }

	    /* Returns the current instance: */

	    parni_bin<value_t,1,rng_t>* last_bin()
	    {
		return this;
	    }

	    /* Returns the current const instance: */

	    const parni_bin<value_t,1,rng_t>* last_bin() const
	    {
		return this;
	    }

	    /* Returns the right-most leaf of the binary tree with the current
	     * bin as root: */

	    parni_bin<value_t,1,rng_t>* last_leaf()
	    {
		parni_bin<value_t,1,rng_t>* b=this;
		while(b->child2!=NULL)
		{
		    b=b->child2;
		}
		return b;
	    }
	    const parni_bin<value_t,1,rng_t>* last_leaf() const
	    {
		parni_bin<value_t,1,rng_t>* b=this;
		while(b->child2!=NULL)
		{
		    b=b->child2;
		}
		return b;
	    }

	    /* Returns the previous bin in the left-up-right traversal method: */

	    parni_bin<value_t,1,rng_t>* previous_bin()
	    {
		if(!is_leaf())
		{
		    return child2;
		}
		if(is_second())
		{
		    return parent->child1;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1);
	    }

	    /* Returns the const previous bin in the left-up-right traversal method: */

	    const parni_bin<value_t,1,rng_t>* previous_bin() const
	    {
		if(!is_leaf())
		{
		    return child2;
		}
		if(is_second())
		{
		    return parent->child1;
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1);
	    }

	    /* Returns the previous leaf: */

	    parni_bin<value_t,1,rng_t>* previous_leaf()
	    {
		if(!is_leaf())
		{
		    return last_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child2)
		{
		    return parent->child1->last_leaf();
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1->last_leaf());
	    }

	    /* Returns the previous const leaf: */

	    const parni_bin<value_t,1,rng_t>* previous_leaf() const
	    {
		if(!is_leaf())
		{
		    return last_leaf();
		}
		if(parent==NULL)
		{
		    return NULL;
		}
		if(this==parent->child2)
		{
		    return parent->child1->last_leaf();
		}
		parni_bin<value_t,1,rng_t>* b=parent;
		while(b->is_first())
		{
		    b=b->parent;
		}
		return (b->parent==NULL)?NULL:(b->parent->child1->last_leaf());
	    }

	    /* Updates the weights up the tree with the sum rule: */

	    void update_weight()
	    {
		if(child1!=NULL and child2!=NULL)
		{
		    w=child1->w+child2->w;
		}
		if(parent!=NULL)
		{
		    parent->update_weight();
		}
	    }
    };
    template<class value_t,class rng_t,class key_type>const key_type parni_bin<value_t,1,rng_t,key_type>::max_depth=sizeof(key_type)*CHAR_BIT-1;

    /* Output stream operator overload: */

    template<class value_t,std::size_t D,class rng_t>std::ostream& operator << (std::ostream& os,const parni_bin<value_t,D,rng_t>& bin)
    {
	return bin.print(os);
    }
}

#endif /*CAMGEN_PARNI_BIN_H_*/

