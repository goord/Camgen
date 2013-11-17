//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file s_gen_grid.h
    \brief Adaptve propagator sampler class.
 */

#ifndef CAMGEN_S_GEN_GRID_H_
#define CAMGEN_S_GEN_GRID_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The adaptive propagator sampler classes consist of a s-generator map and an   *
 * underlying adaptive grid streaming the random numbers which are mapped to the *
 * s-values.                                                                     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/s_gen.h>
#include <Camgen/parni.h>

namespace Camgen
{
    template<class value_t,class rng_t>class adaptive_s_generator: public s_generator<value_t,rng_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef s_generator<value_t,rng_t> map_type;
	    typedef parni<value_t,1,rng_t> grid_type;
	    typedef parni_sub_grid<value_t,1,rng_t> sub_grid_type;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Constructor: */

	    adaptive_s_generator(map_type* mapping_,size_type N_bins=grid_bins(),grid_modes::type mode=grid_mode()):mapping(mapping_),grid(new grid_type(&r,0,1,N_bins,mode)),rmin(0),rmax(1),subgrid(new sub_grid_type(grid,rmin,rmax))
	    {
		mapping->set_s(this->get_s());
	    }

	    /* Constructor with invariant mass pointer argument: */

	    adaptive_s_generator(value_type* s,map_type* mapping_,size_type N_bins=grid_bins(),grid_modes::type mode=grid_mode()):map_type(s),mapping(mapping_),grid(new grid_type(&r,0,1,N_bins,mode)),rmin(0),rmax(1),subgrid(new sub_grid_type(grid,rmin,rmax))
	    {
		mapping->set_s(s);
	    }

	    /* Destructor: */

	    ~adaptive_s_generator()
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

	    /* Sets the minimal minimal invariant mass-squared: */

	    bool set_s_min_min(const value_type& sminmin)
	    {
		return mapping->set_s_min(sminmin);
	    }

	    /* Sets the maximal maximal invariant mass-squared: */

	    bool set_s_max_max(const value_type& smaxmax)
	    {
		return mapping->set_s_max(smaxmax);
	    }

	    /* Minimal invariant mass refresher: */

	    bool refresh_s_min()
	    {
		if(this->s_min()<mapping->s_min())
		{
		    this->norm=(value_type)0;
		    return false;
		}
		rmin=(this->s_min()==mapping->s_min())?0:(mapping->inverse_map(this->s_min()));
		subgrid->set_x_min(rmin);
		refresh_norm();
		return this->normalisable();
	    }

	    /* Maximal invariant mass refresher: */

	    bool refresh_s_max()
	    {
		if(this->s_max()>mapping->s_max())
		{
		    this->norm=(value_type)0;
		    return false;
		}
		rmax=(this->s_max()==mapping->s_max())?1:(mapping->inverse_map(this->s_max()));
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

	    bool refresh_s_bounds()
	    {
		if(this->s_min()<mapping->s_min() or this->s_max()>mapping->s_max())
		{
		    this->norm=(value_type)0;
		    return false;
		}
		rmin=(this->s_min()==mapping->s_min())?0:(mapping->inverse_map(this->s_min()));
		rmax=(this->s_max()==mapping->s_max())?1:(mapping->inverse_map(this->s_max()));
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
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<t<<": "<<r<<" not in interval ["<<rmin<<','<<rmax<<"] corresponding to s-range ["<<this->s_min()<<','<<this->s_max()<<']'<<endlog;
		    this->weight()=(value_type)0;
		    return false;
		}
		this->s()=mapping->map(r);
		mapping->weight()=mapping->get_fast_weight();
		this->weight()=subgrid->weight()*mapping->weight();
		return true;
	    }

	    /* Weight evaluation method: */

	    bool evaluate_weight()
	    {
		if(this->normalisable() and mapping->evaluate_weight())
		{
		    r=mapping->inverse_map(mapping->s());
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
		grid->integrand()=this->integrand()*mapping->weight();
		grid->update();
	    }

	    /* Total update method: */

	    void update_weight()
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
		subgrid->adapt();
	    }

	    /* Overridden resetting method */

	    void reset()
	    {
		this->MC_generator<value_t>::reset();
		grid->reset();
		subgrid->reset();
	    }

	    /* Overridden cross section resetting method */

	    void reset_cross_section()
	    {
		this->MC_generator<value_t>::reset_cross_section();
		grid->reset_cross_section();
		subgrid->reset_cross_section();
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Double-dispatch integrator functions. Unknown and should not be
	     * used. */
	    
	    value_type integrate_with(const s_generator<value_t,rng_t>* other,const value_type& sqrts) const
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"constrained integral with adaptive grid map requested--returning 0"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const BW_s_generator<value_t,rng_t>* other,const value_type& sqrts) const
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"constrained integral with adaptive grid map requested--returning 0"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const pl_s_generator<value_t,rng_t>* other,const value_type& sqrts) const
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"Constrained integral with adaptive grid map requested--returning 0"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const uni_s_generator<value_t,rng_t>* other,const value_type& sqrts) const
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"Constrained integral with adaptive grid map requested--returning 0"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const Dd_s_generator<value_t,rng_t>* other,const value_type& sqrts) const
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"Constrained integral with adaptive grid map requested--returning 0"<<endlog;
		return 0;
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

	    /* Overridden loading method. */

	    std::istream& load_data(std::istream& is)
	    {
		mapping->load(is);
		if(grid!=NULL)
		{
		    delete grid;
		}
		grid=grid_type::create_instance(&r,is);
		if(subgrid!=NULL)
		{
		    delete subgrid;
		}
		subgrid=new parni_sub_grid<value_t,1,rng_t>(grid);
		subgrid->load(is);
		return is;
	    }

	    /* Overridden saving method. */

	    std::ostream& save_data(std::ostream& os) const
	    {
		mapping->save(os);
		grid->save(os);
		subgrid->save(os);
		return os;
	    }

	protected:

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

	private:

	    /* Transformation: */

	    map_type* const mapping;

	    /* Adaptive grid: */

	    grid_type* grid;
	    
	    /* Grid veriables: */
	    
	    value_type r,rmin,rmax;

	    /* Subinterval grid instance: */

	    parni_sub_grid<value_t,1,rng_t>* subgrid;
    };
}

#endif /*CAMGEN_S_GEN_GRID_H_*/

