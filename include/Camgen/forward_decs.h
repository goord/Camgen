//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FORWARD_DECS_H_
#define CAMGEN_FORWARD_DECS_H_

#include <cstddef>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Forward declarations of class templates which depend upon each other. *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
	template<class model_t>class anti_fermion_propagator;
	template<class model_t>class anti_fermion_wave_function;
	template<class Feynrule_t>class cfd_evaluate;
	template<class model_t,std::size_t N_in,std::size_t N_out>class CM_algorithm;
	
	/// Colour generator abstract base class template.
	/// The template parameter denote the model type, incoming and outgoing number of particles in the
	/// process, and a boolean denoting whether the colours are continuous. 
	
	template<class model_t,std::size_t N_in,std::size_t N_out>class current_tree;
	template<class model_t,std::size_t N_in,std::size_t N_out>class process_tree;
	template<class model_t,std::size_t N,bool decomp>class current;
	template<class model_t,std::size_t N>class current_base;
	template<class Feynrule_t>class evaluate;
	template<class model_t>class particle_family;
	template<class model_t>class fermion_propagator;
	template<class model_t>class fermion_wave_function;
	template<class model_t>class Feynman_gauge;
	template<class value_t>class flavour_comp;
	template<class value_t>class flavour_eq;
	template<class value_t>class flavour_neq;
	template<class value_t>class flavour_select;
	template<class value_t>class flavourvec_comp;
	template<class model_t>class fusion_class;
	template<class model_t>class has_leg;
	template<class model_t>class has_legs;
	
	/// Helicity generator abstract base class template.
	/// The template parameter denote the model type, incoming and outgoing number of particles in the
	/// process, and a boolean denoting whether the helicities are continuous. 
	
	template<class model_t,std::size_t N,bool decomp>class interaction;
	template<class model_t,std::size_t N>class interaction_base;
	template<class model_t>class m_anti_fermion_wave_function;
	template<class model_t>class m_fermion_wave_function;
	template<class model_t,std::size_t dim=model_t::dimension,class spacetime_t=typename model_t::spacetime_type,class Dirac_t=typename model_t::Dirac_algebra_type>class m_vector_wave_function;
	template<class model_t>class model;
	template<class model_t>class model_wrapper;
	template<class value_t>class name_comp;
	template<class value_t>class name_eq;
	template<class value_t>class name_neq;
	template<class value_t>class name_select;

	/// Single-particle phase space class. The first template parameter
	/// denotes the model type, the second parameter the spacetime
	/// dimensions and the third the coloured flag.
	
	template<class model_t,std::size_t dim=model_t::dimension,bool q=model_t::coloured>class particle_ps;
	template<class model_t>class particle;
	template<class model_t>class R_xi_gauge;
	template<class model_t>class scalar_propagator;
	template<class model_t>class scalar_wave_function;
	template<class model_t,bool cont_hels=model_t::continuous_helicities>class singlet_ps;
	template<class model_t>class unitary_gauge;
	template<class model_t,std::size_t dim=model_t::dimension,class spacetime_t=typename model_t::spacetime_type,class Dirac_t=typename model_t::Dirac_algebra_type>class vector_wave_function;
	template<class model_t>class vertex;
}

#endif /*CAMGEN_FORWARD_DECS_H_*/

