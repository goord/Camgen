//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file val_gen_grid.h
    \brief Adaptve propagator sampler class.
 */

#ifndef CAMGEN_VAL_GEN_GRID_H_
#define CAMGEN_VAL_GEN_GRID_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Adaptive value generator class definition, consisting of a smoothing mapping  *
 * with an underlying adaptive grid.                                             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/val_gen.h>
#include <Camgen/parni_gen.h>

namespace Camgen
{
    template<class value_t,class rng_t>class adaptive_value_generator: public value_generator<value_t,rng_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef value_generator<value_t,rng_t> base_type;
	    typedef parni_generator<value_t,1,rng_t> grid_type;
	    typedef parni_sub_grid<value_t,1,rng_t> sub_grid_type;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Constructor: */

	    adaptive_value_generator(base_type* mapping_,size_type N_bins=grid_bins(),grid_modes::type mode=grid_mode()):mapping(mapping_),grid(new grid_type(&r,0,1,N_bins,mode)),rmin(0),rmax(1),subgrid(grid->create_sub_grid(rmin,rmax))
	    {
		mapping->set_value(this->get_value());
	    }

	    /* Destructor: */

	    ~adaptive_value_generator()
	    {
		if(mapping!=NULL)
		{
		    delete mapping;
		}
		if(subgrid!=NULL)
		{
		    delete subgrid;
		}
		if(grid!=NULL)
		{
		    delete grid;
		}
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /* Mapping implementation: */

	    value_type map(const value_type& r) const
	    {
		return mapping->map(r);
	    }

	    /* Inverse mapping implementation: */

	    value_type inverse_map(const value_type& s) const
	    {
		return mapping->inverse_map(s);
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    void set_value(value_type* x_)
	    {
		this->base_type::set_value(x_);
		mapping->set_value(x_);
	    }

	    /* Sets the minimal minimal invariant mass-squared: */

	    bool set_min_lower_bound(const value_type& xmin)
	    {
		return mapping->set_lower_bound(xmin) and this->set_lower_bound(xmin);
	    }

	    /* Sets the maximal maximal invariant mass-squared: */

	    bool set_max_upper_bound(const value_type& xmax)
	    {
		return mapping->set_upper_bound(xmax) and this->set_upper_bound(xmax);
	    }

	    /* Minimal invariant mass refresher: */

	    bool refresh_lower_bound()
	    {
		if(this->lower_bound()<mapping->lower_bound())
		{
		    this->norm=(value_type)0;
		    return false;
		}
		rmin=(this->lower_bound()==mapping->lower_bound())?0:(mapping->inverse_map(this->lower_bound()));
		subgrid->set_x_min(rmin);
		refresh_norm();
		return this->normalisable();
	    }

	    /* Maximal invariant mass refresher: */

	    bool refresh_upper_bound()
	    {
		if(this->upper_bound()>mapping->upper_bound())
		{
		    this->norm=(value_type)0;
		    return false;
		}
		rmax=(this->upper_bound()==mapping->upper_bound())?1:(mapping->inverse_map(this->upper_bound()));
		subgrid->set_x_max(rmax);
		refresh_norm();
		return this->normalisable();
	    }

	    /* Parameter refresher function: */

	    bool refresh_params()
	    {
		return mapping->refresh_params();
	    }

	    /* Invariant mass bounds refresher: */

	    bool refresh_bounds()
	    {
		if(this->lower_bound()<mapping->lower_bound() or this->upper_bound()>mapping->upper_bound())
		{
		    this->norm=(value_type)0;
		    return false;
		}
		rmin=(this->lower_bound()==mapping->lower_bound())?0:(mapping->inverse_map(this->lower_bound()));
		rmax=(this->upper_bound()==mapping->upper_bound())?1:(mapping->inverse_map(this->upper_bound()));
		subgrid->set_x_min(rmin);
		subgrid->set_x_max(rmax);
		refresh_norm();
		return this->normalisable();
	    }

	    /* Generation implementation: */

	    bool generate()
	    {
		if(!this->normalisable())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		subgrid->generate();
		if(r<rmin or r>rmax)
		{
		    std::string t=detailed_type();
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<t<<": "<<r<<" not in interval ["<<rmin<<','<<rmax<<"] corresponding to s-range ["<<this->lower_bound()<<','<<this->upper_bound()<<']'<<endlog;
		    this->weight()=(value_type)0;
		    return false;
		}
		this->value()=mapping->map(r);
		mapping->weight()=mapping->get_fast_weight();
		this->weight()=subgrid->weight()*mapping->weight();
		return true;
	    }

	    /* Weight evaluation method: */

	    bool evaluate_weight()
	    {
		if(this->normalisable() and mapping->evaluate_weight())
		{
		    r=mapping->inverse_map(mapping->value());
		    subgrid->evaluate_weight();
		    this->weight()=subgrid->weight()*mapping->weight();
		    return true;
		}
		this->weight()=(value_type)0;
		return false;
	    }

	    /* Mapping weight factor: */

	    value_type mapping_weight() const
	    {
		return mapping->weight();
	    }

	    /* Adaptive grid weight factor: */

	    value_type grid_weight() const
	    {
		return subgrid->weight();
	    }

	    /* Grid update method: */

	    void update()
	    {
		if(grid->weight()!=(value_type)0)
		{
		    grid->integrand()=this->integrand()/grid->weight();
		    grid->update();
		}
	    }

	    /* Adaptive method: */

	    void adapt()
	    {
		grid->adapt();
	    }

	    /* Overridden resetting method */

	    void reset()
	    {
		grid->reset();
	    }

	    /* Overridden cross section resetting method */

	    void reset_cross_section()
	    {
		this->MC_generator<value_t>::reset_cross_section();
		grid->reset_cross_section();
		subgrid->reset_cross_section();
	    }

	    /* Normalisation refreshing method: */

	    void refresh_norm()
	    {
		if((!mapping->normalisable()) or rmin<(value_type)0 or rmax>(value_type)1 or rmin>rmax)
		{
		    this->norm=(value_type)0;
		}
		else
		{
		    this->norm=(value_type)1;
		}
	    }

	    /* Returns the weight without bound-checking: */

	    value_type get_fast_weight() const
	    {
		subgrid->evaluate_weight();
		return subgrid->weight()*mapping->get_fast_weight();
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Type output: */

	    std::string type() const
	    {
		return mapping->type()+"*";
	    }

	    /* Type output: */

	    std::string detailed_type() const
	    {
		return mapping->detailed_type()+"*";
	    }

	    /* Adds the rectangles to the plotscript pointer: */

	    plot_script* print_rectangles(plot_script* ps) const
	    {
		return grid->print_rectangles(ps);
	    }

	    /* Overridden printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		mapping->print(os);
		os<<"*";
		return os;
	    }

	private:

	    /* Transformation: */

	    base_type* const mapping;

	    /* Adaptive grid: */

	    grid_type* grid;
	    
	    /* Grid veriables: */
	    
	    value_type r,rmin,rmax;

	    /* Subinterval grid instance: */

	    parni_sub_grid<value_t,1,rng_t>* subgrid;
    };
}

#endif /*CAMGEN_VAL_GEN_GRID_H_*/

