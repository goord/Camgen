//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FERMION_H_
#define CAMGEN_FERMION_H_

#include <Camgen/spin.h>
#include <Camgen/m_spinor_fac.h>
#include <Camgen/spinor_contr.h>
#include <Camgen/fermion_prop.h>
#include <Camgen/Minkowski.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * File including all necessary means to do computations involving fermions. The *
 * fermion classes simply include static storage of function pointers to the in- *
 * and outgoing positive and negative helicity spinor wave functions, and type   *
 * definitions to the contraction and propagator classes.                        *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Class declarations: */

    template<class model_t,class spacetime_t=typename model_t::spacetime_type,std::size_t dim=model_t::dimension>class fermion;
    template<class model_t,class spacetime_t=typename model_t::spacetime_type,std::size_t dim=model_t::dimension>class anti_fermion;

    /* Fermion class: */

    template<class model_t>class fermion<model_t,Minkowski_type,4>
    {
	public:

	    /* Anti-particle type definition: */

	    typedef anti_fermion<model_t,Minkowski_type,4> anti_particle_type;

	    /* Wave function, propagator and propagator refresher function pointer type
	     * definitions: */

	    typedef typename get_basic_types<model_t>::size_type size_type;
	    typedef typename get_basic_types<model_t>::wave_func wave_func;
	    typedef typename get_basic_types<model_t>::prop_func prop_func;
	    typedef typename get_basic_types<model_t>::prop_refresher prop_refresher;

	    /* Contraction type definition: */
	     
	    typedef spinor_contraction<model_t> contraction_type;

	    /* Propagator type definition: */

	    typedef fermion_propagator<model_t> propagator_type; 

	private:

	    /* Private utility type definitions: */

	    typedef typename get_Dirac_algebra_type<model_t>::type Dirac_algebra_type;
	    typedef massless_spinor_factory<Dirac_algebra_type,model_t::beam_direction,4> spinor_fac;
	    typedef massive_spinor_factory<Dirac_algebra_type,typename model_t::spin_vector_type,model_t::beam_direction,4> m_spinor_fac;
	
	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spinor_fac::initialise();
		m_spinor_fac::initialise();
		propagator_type::initialise();
		evaluate<contraction_type>::initialise();
	    }
	    
	    /* Function returning the index ranges: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& r)
	    {
		r.clear();
		r.push_back(Dirac_algebra_type::base_type::index_range);
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
    template<class model_t>const spin fermion<model_t,Minkowski_type,4>::particle_spin=spin::half_integer(1);
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::pos_hel_incoming_massless_states={{&fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_p}};
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::neg_hel_incoming_massless_states={{&fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_m}};
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::pos_hel_outgoing_massless_states={{&fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_p_bar}};
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::neg_hel_outgoing_massless_states={{&fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_m_bar}};
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::pos_hel_incoming_massive_states={{&fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_u_p}};
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::neg_hel_incoming_massive_states={{&fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_u_m}};
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::pos_hel_outgoing_massive_states={{&fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_u_p_bar}};
    template<class model_t>const vector<typename fermion<model_t,Minkowski_type,4>::wave_func,1> fermion<model_t,Minkowski_type,4>::neg_hel_outgoing_massive_states={{&fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_u_m_bar}};
    template<class model_t>const typename fermion<model_t,Minkowski_type,4>::wave_func fermion<model_t,Minkowski_type,4>::zero_hel_incoming_massless_state=NULL;
    template<class model_t>const typename fermion<model_t,Minkowski_type,4>::wave_func fermion<model_t,Minkowski_type,4>::zero_hel_outgoing_massless_state=NULL;
    template<class model_t>const typename fermion<model_t,Minkowski_type,4>::wave_func fermion<model_t,Minkowski_type,4>::zero_hel_incoming_massive_state=NULL;
    template<class model_t>const typename fermion<model_t,Minkowski_type,4>::wave_func fermion<model_t,Minkowski_type,4>::zero_hel_outgoing_massive_state=NULL;
    
    /* Antifermion class: */

    template<class model_t>class anti_fermion<model_t,Minkowski_type,4>
    {
	public:

	    /* Anti-particle type definition: */

	    typedef fermion<model_t,Minkowski_type,4> anti_particle_type;

	    /* Wave function, propagator and propagator refresher function pointer type
	     * definitions: */

	    typedef typename get_basic_types<model_t>::size_type size_type;
	    typedef typename get_basic_types<model_t>::wave_func wave_func;
	    typedef typename get_basic_types<model_t>::prop_func prop_func;
	    typedef typename get_basic_types<model_t>::prop_refresher prop_refresher;

	    /* Contraction type definition: */
	    
	    typedef spinor_contraction<model_t> contraction_type;

	    /* Propagator type definition: */

	    typedef anti_fermion_propagator<model_t> propagator_type;

	private:

	    /* Private utility type definitions: */

	    typedef typename get_Dirac_algebra_type<model_t>::type Dirac_algebra_type;
	    typedef massless_spinor_factory<Dirac_algebra_type,model_t::beam_direction,4> spinor_fac;
	    typedef massive_spinor_factory<Dirac_algebra_type,typename model_t::spin_vector_type,model_t::beam_direction,4> m_spinor_fac;
	
	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spinor_fac::initialise();
		m_spinor_fac::initialise();
		propagator_type::initialise();
		evaluate<contraction_type>::initialise();
	    }
	    
	    /* Function returning the index ranges: */
	    
	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& r)
	    {
		r.clear();
		r.push_back(Dirac_algebra_type::base_type::index_range);
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
    template<class model_t>const spin anti_fermion<model_t,Minkowski_type,4>::particle_spin=spin::half_integer(1);
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::pos_hel_incoming_massless_states={{&anti_fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_p_bar}};
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::neg_hel_incoming_massless_states={{&anti_fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_m_bar}};
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::pos_hel_outgoing_massless_states={{&anti_fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_p}};
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::neg_hel_outgoing_massless_states={{&anti_fermion<model_t,Minkowski_type,4>::spinor_fac::make_u_m}};
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::pos_hel_incoming_massive_states={{&anti_fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_v_p_bar}};
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::neg_hel_incoming_massive_states={{&anti_fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_v_m_bar}};
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::pos_hel_outgoing_massive_states={{&anti_fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_v_p}};
    template<class model_t>const vector<typename anti_fermion<model_t,Minkowski_type,4>::wave_func,1> anti_fermion<model_t,Minkowski_type,4>::neg_hel_outgoing_massive_states={{&anti_fermion<model_t,Minkowski_type,4>::m_spinor_fac::make_v_m}};
    template<class model_t>const typename anti_fermion<model_t,Minkowski_type,4>::wave_func anti_fermion<model_t,Minkowski_type,4>::zero_hel_incoming_massless_state=NULL;
    template<class model_t>const typename anti_fermion<model_t,Minkowski_type,4>::wave_func anti_fermion<model_t,Minkowski_type,4>::zero_hel_outgoing_massless_state=NULL;
    template<class model_t>const typename anti_fermion<model_t,Minkowski_type,4>::wave_func anti_fermion<model_t,Minkowski_type,4>::zero_hel_incoming_massive_state=NULL;
    template<class model_t>const typename anti_fermion<model_t,Minkowski_type,4>::wave_func anti_fermion<model_t,Minkowski_type,4>::zero_hel_outgoing_massive_state=NULL;
}

#endif /*CAMGEN_FERMION_H_*/

