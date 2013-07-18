//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_M_FERMION_WFUNC_H_
#define CAMGEN_M_FERMION_WFUNC_H_

#include <Camgen/m_spinor_fac.h>
#include <Camgen/spinor_contr.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the massive fermion and anti-fermion wave *
 * function classes. The wave function class dispatch function pointers    *
 * depending on the argument integer, which denotes the helicity of the    *
 * requested wave function.                                                *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Massive fermion wave function class: */

    template<class model_t>class m_fermion_wave_function
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	    
	    /* Definition of the wave function pointer type, as defined in the
	     * get_basic_types class template: */
	    
	    typedef typename get_basic_types<model_t>::wave_func wave_func;

	    /* Type definition of the wave function contraction algorithm: */

	    typedef spinor_contraction<model_t> contraction_type;

	private:

	    /* Type definitions of the Dirac algebra and massive spinor factory:
	     * */

	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);
	    typedef massive_spinor_factory<Dirac_algebra_type,typename model_t::spin_vector_type,model_t::beam_direction,model_t::dimension> spinor_fac;

	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spinor_fac::initialise();
	    }

	    /* Function returning the index ranges of a subamplitude that is
	     * filled by the wave functions : */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& r)
	    {
		r.clear();
		r.push_back(Dirac_algebra_type::base_type::index_range);
		return r;
	    }

	    /* Incoming wave function dispatcher: */

	    static wave_func dispatch_incoming_wave_function(int s)
	    {
		switch(s)
		{
		    case -1:
			return &spinor_fac::make_u_m;
			break;
		    case 0:
			return &make_u_L;
			break;
		    case 1:
			return &spinor_fac::make_u_p;
			break;
		}
		CAMGEN_ERROR("requested incoming helicity state out of range");
		return NULL;
	    }

	    /* Outgoing wave function dispatcher: */

	    static wave_func dispatch_outgoing_wave_function(int s)
	    {
		switch(s)
		{
		    case -1:
			return &spinor_fac::make_u_m_bar;
			break;
		    case 0:
			return &make_u_L;
			break;
		    case 1:
			return &spinor_fac::make_u_p_bar;
			break;
		}
		CAMGEN_ERROR("requested incoming helicity state out of range");
		return NULL;
	    }

	private:

	    /* 'Longitudinal fermion' dummy wave function, is being dispatched,
	     * but not supposed to be called: */

	    static void make_u_L(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_MESSAGE("zero helicity for fermion encountered and ignored");
	    }

	    /* 'Longitudinal fermion' dummy wave function for continuous
	     * helicities will not produce any message: */

	    static void make_u_L(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m){}
    };

    /* Massive anti-fermion wave function class: */

    template<class model_t>class m_anti_fermion_wave_function
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	    
	    /* Definition of the wave function pointer type, as defined in the
	     * get_basic_types class template: */
	    
	    typedef typename get_basic_types<model_t>::wave_func wave_func;

	    /* Type definition of the wave function contraction algorithm: */

	    typedef spinor_contraction<model_t> contraction_type;

	private:

	    /* Type definitions of the Dirac algebra and massive spinor factory:
	     * */

	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);
	    typedef massive_spinor_factory<Dirac_algebra_type,typename model_t::spin_vector_type,model_t::beam_direction,model_t::dimension> spinor_fac;

	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spinor_fac::initialise();
	    }

	    /* Function returning the index ranges of a subamplitude that is
	     * filled by the wave functions : */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& r)
	    {
		r.clear();
		r.push_back(Dirac_algebra_type::base_type::index_range);
		return r;
	    }

	    /* Incoming wave function dispatcher: */

	    static wave_func dispatch_incoming_wave_function(int s)
	    {
		switch(s)
		{
		    case -1:
			return &spinor_fac::make_v_m_bar;
			break;
		    case 0:
			return &make_v_L;
			break;
		    case 1:
			return &spinor_fac::make_v_p_bar;
			break;
		}
		CAMGEN_ERROR("requested incoming helicity state out of range");
		return NULL;
	    }

	    /* Outgoing wave function dispatcher: */

	    static wave_func dispatch_outgoing_wave_function(int s)
	    {
		switch(s)
		{
		    case -1:
			return &spinor_fac::make_v_m;
			break;
		    case 0:
			return &make_v_L;
			break;
		    case 1:
			return &spinor_fac::make_v_p;
			break;
		}
		CAMGEN_ERROR("requested incoming helicity state out of range");
		return NULL;
	    }
	private:

	    /* 'Longitudinal anti-fermion' dummy wave function, is being dispatched,
	     * but not supposed to be called: */

	    static void make_v_L(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_MESSAGE("zero helicity for anti-fermion encountered and ignored");
	    }

	    /* 'Longitudinal anti-fermion' dummy wave function for continuous
	     * helicities will not produce any message: */

	    static void make_v_L(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m){}
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_M_FERMION_WFUNC_H_*/

