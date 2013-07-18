//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VECTOR_PARTICLE_H_
#define CAMGEN_VECTOR_PARTICLE_H_

#include <Camgen/spin.h>
#include <Camgen/spinor_fac.h>
#include <Camgen/vec_prop_fg.h>
#include <Camgen/vec_prop_ug.h>
#include <Camgen/R_gauge.h>
#include <Camgen/vector_contr.h>
#include <Camgen/Minkowski.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * File including all necessary means to do computations involving vectors     *
 * particles. The vector_particle class simply includes static storage of      *
 * function pointers to the in- and outgoing positive, negative and zero,      *
 * massless and massive helicity polarisation vectors, and type definitions to *
 * the contraction and propagator classes.                                     *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,class gauge,class spacetime_t=typename model_t::spacetime_type,std::size_t dim=model_t::dimension>class vector_particle;
    
    template<class model_t,class gauge>class vector_particle<model_t,gauge,Minkowski_type,4>
    {
	public:

	    /* Anti-particle type definition: */

	    typedef vector_particle<model_t,gauge,Minkowski_type,4> anti_particle_type;

	    /* Wave function, propagator and propagator refresher function pointer type
	     * definitions: */

	    typedef typename get_basic_types<model_t>::size_type size_type;
	    typedef typename get_basic_types<model_t>::wave_func wave_func;
	    typedef typename get_basic_types<model_t>::prop_func prop_func;
	    typedef typename get_basic_types<model_t>::prop_refresher prop_refresher;

	    /* Contraction type definition: */
	     
	    typedef vector_contraction<model_t> contraction_type;

	    /* Propagator type definition: */

	    typedef gauge propagator_type; 

	private:

	    /* Private utility type definitions: */

	    typedef typename get_Dirac_algebra_type<model_t>::type Dirac_algebra_type;
	    typedef massless_spinor_factory<Dirac_algebra_type,model_t::beam_direction,4> spinor_fac;
	
	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spinor_fac::initialise();
		propagator_type::initialise();
		evaluate<contraction_type>::initialise();
	    }
	    
	    /* Function returning the index ranges: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& r)
	    {
		r.clear();
		r.push_back(model_t::dimension);
		return r;
	    }

	    /* Particle type spin: */

	    static const spin particle_spin;

	    /* Massless helicity spinor wave functions: */

	    static const vector<wave_func,1> pos_hel_incoming_massless_states;
	    static const vector<wave_func,1> neg_hel_incoming_massless_states;
	    static const wave_func zero_hel_incoming_massless_state;
	    
	    static const vector<wave_func,1> pos_hel_outgoing_massless_states;
	    static const vector<wave_func,1> neg_hel_outgoing_massless_states;
	    static const wave_func zero_hel_outgoing_massless_state;
	    
	    /* Massive helicity spinor wave functions: */

	    static const vector<wave_func,1> pos_hel_incoming_massive_states;
	    static const vector<wave_func,1> neg_hel_incoming_massive_states;
	    static const wave_func zero_hel_incoming_massive_state;
	    
	    static const vector<wave_func,1> pos_hel_outgoing_massive_states;
	    static const vector<wave_func,1> neg_hel_outgoing_massive_states;
	    static const wave_func zero_hel_outgoing_massive_state;
    };
    template<class model_t,class gauge>const spin vector_particle<model_t,gauge,Minkowski_type,4>::particle_spin=spin::integer(1);
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::pos_hel_incoming_massless_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_e_p}};
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::neg_hel_incoming_massless_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_e_m}};
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::pos_hel_outgoing_massless_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_e_p}};
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::neg_hel_outgoing_massless_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_e_m}};
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::pos_hel_incoming_massive_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_m_e_p}};
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::neg_hel_incoming_massive_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_m_e_m}};
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::pos_hel_outgoing_massive_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_m_e_p}};
    template<class model_t,class gauge>const vector<typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func,1> vector_particle<model_t,gauge,Minkowski_type,4>::neg_hel_outgoing_massive_states={{&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_m_e_m}};
    template<class model_t,class gauge>const typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func vector_particle<model_t,gauge,Minkowski_type,4>::zero_hel_incoming_massless_state=NULL;
    template<class model_t,class gauge>const typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func vector_particle<model_t,gauge,Minkowski_type,4>::zero_hel_outgoing_massless_state=NULL;
    template<class model_t,class gauge>const typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func vector_particle<model_t,gauge,Minkowski_type,4>::zero_hel_incoming_massive_state=&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_m_e_L;
    template<class model_t,class gauge>const typename vector_particle<model_t,gauge,Minkowski_type,4>::wave_func vector_particle<model_t,gauge,Minkowski_type,4>::zero_hel_outgoing_massive_state=&vector_particle<model_t,gauge,Minkowski_type,4>::spinor_fac::make_m_e_L;
}

#endif /*CAMGEN_VECTOR_PARTICLE_H_*/

