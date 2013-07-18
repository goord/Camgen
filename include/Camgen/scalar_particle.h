//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SCALAR_PARTICLE_H_
#define CAMGEN_SCALAR_PARTICLE_H_

#include <Camgen/spin.h>
#include <Camgen/scalar_contr.h>
#include <Camgen/scalar_prop.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * File including all necessary means to do computations involving scalars...*
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class scalar_particle
    {
	public:

	    /* Anti-particle type definition: */

	    typedef scalar_particle<model_t> anti_particle_type;

	    /* Wave function, propagator and propagator refresher function pointer type
	     * definitions: */

	    typedef typename get_basic_types<model_t>::wave_func wave_func;
	    typedef typename get_basic_types<model_t>::prop_func prop_func;
	    typedef typename get_basic_types<model_t>::prop_refresher prop_refresher;

	    /* Contraction type definition: */
	     
	    typedef scalar_contraction<model_t> contraction_type;

	    /* Propagator type definition: */

	    typedef scalar_propagator<model_t> propagator_type; 

	private:

	    /* Private utility type definitions: */

	    typedef typename get_basic_types<model_t>::size_type size_type;
	    typedef typename get_basic_types<model_t>::value_type value_type;
	    typedef typename get_basic_types<model_t>::r_value_type r_value_type;
	    typedef typename get_basic_types<model_t>::tensor_type tensor_type;
	    typedef typename get_basic_types<model_t>::momentum_type momentum_type;
	    typedef typename get_basic_types<model_t>::iterator iterator;
	
	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		propagator_type::initialise();
	    }
	    
	    /* Function returning the index ranges: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& r)
	    {
		r.clear();
		return r;
	    }

	    /* Particle type spin: */

	    static const spin particle_spin;

	    /* Massless helicity spinor wave functions: */

	    static const vector<wave_func,0> pos_hel_incoming_massless_states;
	    static const vector<wave_func,0> neg_hel_incoming_massless_states;
	    static const wave_func zero_hel_incoming_massless_state;
	    
	    static const vector<wave_func,0> pos_hel_outgoing_massless_states;
	    static const vector<wave_func,0> neg_hel_outgoing_massless_states;
	    static const wave_func zero_hel_outgoing_massless_state;
	    
	    /* Massive helicity spinor wave functions: */

	    static const vector<wave_func,0> pos_hel_incoming_massive_states;
	    static const vector<wave_func,0> neg_hel_incoming_massive_states;
	    static const wave_func zero_hel_incoming_massive_state;
	    
	    static const vector<wave_func,0> pos_hel_outgoing_massive_states;
	    static const vector<wave_func,0> neg_hel_outgoing_massive_states;
	    static const wave_func zero_hel_outgoing_massive_state;

	    /* The scalar wave function for discrete helicities: */
	
	    static void apply(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()==0),"tensor iterator out of range");
		(*it)=value_type(1,0);
	    }

	    /* Scalar wave function for continuous helicities: */

	    static void apply(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()==0),"tensor iterator out of range");
		(*it)=h;
	    }
    };
    template<class model_t>const spin scalar_particle<model_t>::particle_spin=spin::integer(0);
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::pos_hel_incoming_massless_states={};
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::neg_hel_incoming_massless_states={};
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::pos_hel_outgoing_massless_states={};
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::neg_hel_outgoing_massless_states={};
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::pos_hel_incoming_massive_states={};
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::neg_hel_incoming_massive_states={};
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::pos_hel_outgoing_massive_states={};
    template<class model_t>const vector<typename scalar_particle<model_t>::wave_func,0> scalar_particle<model_t>::neg_hel_outgoing_massive_states={};
    template<class model_t>const typename scalar_particle<model_t>::wave_func scalar_particle<model_t>::zero_hel_incoming_massless_state=&scalar_particle<model_t>::apply;
    template<class model_t>const typename scalar_particle<model_t>::wave_func scalar_particle<model_t>::zero_hel_outgoing_massless_state=&scalar_particle<model_t>::apply;
    template<class model_t>const typename scalar_particle<model_t>::wave_func scalar_particle<model_t>::zero_hel_incoming_massive_state=&scalar_particle<model_t>::apply;
    template<class model_t>const typename scalar_particle<model_t>::wave_func scalar_particle<model_t>::zero_hel_outgoing_massive_state=&scalar_particle<model_t>::apply;
}

#endif /*CAMGEN_SCALAR_PARTICLE_H_*/

