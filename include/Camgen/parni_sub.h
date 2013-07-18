//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file parni_sub.h
    \brief Parni sub-grid integrator class.
 */

#ifndef CAMGEN_PARNI_SUB_H_
#define CAMGEN_PARNI_SUB_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The parni sub-grid Monte Carlo generator generates according to a predefined  *
 * parni-grid, within a specified interval.                                      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/parni.h>

namespace Camgen
{
    template<class value_t,class rng_t,class key_t>class parni_sub_grid<value_t,1,rng_t,key_t>: public MC_object_generator<value_t,typename parni_bin<value_t,1,rng_t,key_t>::point_type,1,true>
    {
	friend class parni<value_t,1,rng_t,key_t>;

	/* Type definitions: */

	typedef parni<value_t,1,rng_t,key_t> grid_type;
	typedef typename grid_type::base_type base_type;

	public:

	    typedef typename grid_type::value_type value_type;
	    typedef typename grid_type::rn_stream rn_stream;
	    typedef typename grid_type::size_type size_type;
	    typedef typename grid_type::bin_type bin_type;
	    typedef typename grid_type::point_type point_type;

	    /// Constructor over total interval.
	    
	    parni_sub_grid(grid_type* grid_):base_type(grid_->get_object()),grid(grid_),xmin(grid->xmin),xmax(grid->xmax),reg_bin(grid->reg_bin),smin(0),smax(grid_->norm())
	    {
		grid->register_sub_grid(this);
	    }

	    /// Construct over sub-interval.

	    parni_sub_grid(grid_type* grid_,const point_type& xmin_,const point_type& xmax_):base_type(grid_->get_object()),grid(grid_),xmin(xmin_),xmax(xmax_),reg_bin(grid->reg_bin)
	    {
		refresh_bounds();
		grid->register_sub_grid(this);
	    }

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
		bool q=reg_bin->generate(this->object(),xmin,xmax);
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
		if(this->object()>=xmin and this->object()<=xmax)
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

	    double grid_norm()
	    {
		return grid->norm();
	    }

	    /// Update method. Updates the mother grid.
	    
	    void update()
	    {
		if(grid!=NULL)
		{
		    grid->integrand()=this->integrand();
		    grid->update();
		}
	    }

	    /// Adaptation method. Adapts the mother grid.

	    void adapt()
	    {
		if(grid!=NULL)
		{
		    grid->adapt(xmin,xmax);
		}
	    }

	    /// Refreshes internal parameters depending on the boundary values.

	    void refresh_bounds()
	    {
		grid->find_leaf(xmin,smin);
		grid->find_leaf(xmax,smax);
	    }

	    /* Overridden loading method: */

	    std::istream& load(std::istream& is)
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<parni_sub>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read"<<endlog;
		    return is;
		}
		this->MC_generator<value_t>::load(is);
		safe_read(is,xmin);
		safe_read(is,xmax);
		refresh_bounds();
		std::string endflag;
		do
		{
		    std::getline(is,endflag);
		}
		while(endflag!="</parni_sub>" and !is.eof());
		return is;
	    }

	    /* Overridden saving method: */

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<parni_sub>"<<std::endl;
		this->MC_generator<value_t>::save(os);
		safe_write(os,xmin);
		os<<"\t";
		safe_write(os,xmax);
		os<<std::endl;
		os<<"</parni_sub>"<<std::endl;
		return os;
	    }

	private:

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

#endif /*CAMGEN_PARNI_SUB_H_*/

