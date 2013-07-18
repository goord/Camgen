//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_M_VECTOR_WFUNC_H_
#define CAMGEN_M_VECTOR_WFUNC_H_

#include <Camgen/forward_decs.h>
#include <Camgen/Minkowski.h>
#include <Camgen/spinor_fac.h>
#include <Camgen/vector_contr.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Wave function class definition for massive vector particles in Camgen. The   *
 * class is declared in the forward declarations header, and its instantiation   *
 * is only nonempty for models in 4-dimensional Minkowksi spacetime. The wave    *
 * functions are called from the spinor factory class, but an explicit           *
 * instantiation is provided without this invocation if the Dirac algebra type   *
 * in the model is defined as 'unused'. In the latter case the wave functions in *
 * the Pauli basis of gamma matrices have been copied to the class and will      *
 * return their address.                                                         *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Massive vector wave function class for models with a Dirac algebra basis defined
     * (in 4d-Minkowksi spacetime): */

    template<class model_t,class Dirac_alg>class m_vector_wave_function<model_t,4,Minkowski_type,Dirac_alg>
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	    /* Definition of the wave function pointer type, as defined in the
	     * get_basic_types class template: */
	    
	    typedef typename get_basic_types<model_t>::wave_func wave_func;

	    /* Type definition of the wave function contraction algorithm: */

	    typedef vector_contraction<model_t> contraction_type;
	    
	private:

	    /* Type definitions of the spacetime type, Dirac algebra type and
	     * massless spinor factory: */
	    
	    DEFINE_SPACETIME_TYPE(model_t);
	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);
	    typedef massless_spinor_factory<Dirac_algebra_type,model_t::beam_direction,4> spinor_fac;
	
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
		r.push_back(4);
		return r;
	    }

	    /* Incoming wave function dispatcher: */

	    static wave_func dispatch_incoming_wave_function(int s)
	    {
		switch(s)
		{
		    case -1:
			return &spinor_fac::make_m_e_m;
			break;
		    case 0:
			return &spinor_fac::make_m_e_L;
			break;
		    case 1:
			return &spinor_fac::make_m_e_p;
			break;
		}
		CAMGEN_ERROR("requested incoming helicity state out of range");
		return NULL;
	    }

	    /* Outgoing wave function dispatcher: */

	    static wave_func dispatch_outgoing_wave_function(int s)
	    {
		return dispatch_incoming_wave_function(s);
	    }
    };

    /* Specialisation for models with no Dirac algebra basis defined: use the
     * Pauli basis overloaded methods in the spinor factory: */

    template<class model_t>class m_vector_wave_function<model_t,4,Minkowski_type,unused>
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	    /* Definition of the wave function pointer type, as defined in the
	     * get_basic_types class template: */
	    
	    typedef typename get_basic_types<model_t>::wave_func wave_func;

	    /* Type definition of the wave function contraction algorithm: */

	    typedef vector_contraction<model_t> contraction_type;
	    
	    /* Initialisation phase: */

	    static void initialise()
	    {
		Minkowski_type::implementation<r_value_type,4>::initialise();
	    }

	    /* Function returning the index ranges of a subamplitude that is
	     * filled by the wave functions : */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& r)
	    {
		r.clear();
		r.push_back(4);
		return r;
	    }

	    /* Incoming wave function dispatcher: */

	    static wave_func dispatch_incoming_wave_function(int s)
	    {
		switch(s)
		{
		    case -1:
			return &make_m_e_m;
			break;
		    case 0:
			return &make_m_e_L;
			break;
		    case 1:
			return &make_m_e_p;
			break;
		}
		CAMGEN_ERROR("requested incoming helicity state out of range");
		return NULL;
	    }

	    /* Outgoing wave function dispatcher: */

	    static wave_func dispatch_outgoing_wave_function(int s)
	    {
		return dispatch_incoming_wave_function(s);
	    }

	private:

	    static const r_value_type sqrthlf;

	    /* Dirac-algebra independent methods to construct massive helicity
	     * vectors (see the Pauli basis specialisation for the massless
	     * spinor factory): */

	    static void make_m_e_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered");

		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		if(P==(r_value_type)0)
		{
		    it[1]+=value_type(sqrthlf,0);
		    it[2]+=value_type(0,sqrthlf);
		    return;
		}
		it[1]-=(r_value_type)0.5/(pT*P);
		it[2]+=pT*value_type((*p)[1]*(*p)[2]/P,(*p)[3]);
		it[3]+=pT*value_type((*p)[1]*(*p)[3]/P,-(*p)[2]); 

	    }
	    static void make_m_e_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered");

		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		if(P==(r_value_type)0)
		{
		    it[1]+=value_type(sqrthlf,0)*h;
		    it[2]+=value_type(0,sqrthlf)*h;
		    return;
		}
		it[1]-=(r_value_type)0.5*h/(pT*P);
		it[2]+=h*pT*value_type((*p)[1]*(*p)[2]/P,(*p)[3]);
		it[3]+=h*pT*value_type((*p)[1]*(*p)[3]/P,-(*p)[2]); 

	    }
	    static void make_m_e_L(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered");

		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		if(P==(r_value_type)0)
		{
		    it[3]+=value_type(1,0);
		    return;
		}
		r_value_type xi=(*p)[0]/((*m)*P);
		it[0]+=P/(*m);
		it[1]+=xi*(*p)[1];
		it[2]+=xi*(*p)[2];
		it[3]+=xi*(*p)[3];

	    }
	    static void make_m_e_L(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered");

		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		if(P==(r_value_type)0)
		{
		    it[3]+=h;
		    return;
		}
		value_type xi=(*p)[0]*h/((*m)*P);
		it[0]+=P*h/(*m);
		it[1]+=xi*(*p)[1];
		it[2]+=xi*(*p)[2];
		it[3]+=xi*(*p)[3];

	    }
	    static void make_m_e_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered");

		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		if(P==(r_value_type)0)
		{
		    it[1]+=value_type(sqrthlf,0);
		    it[2]+=value_type(0,-sqrthlf);
		    return;
		}
		it[1]-=(r_value_type)0.5/(pT*P);
		it[2]+=pT*value_type((*p)[1]*(*p)[2]/P,-(*p)[3]);
		it[3]+=pT*value_type((*p)[1]*(*p)[3]/P,(*p)[2]); 

	    }
	    static void make_m_e_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered");

		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		if(P==(r_value_type)0)
		{
		    it[1]+=value_type(sqrthlf,0)*h;
		    it[2]-=value_type(0,sqrthlf)*h;
		    return;
		}
		it[1]-=(r_value_type)0.5*h/(pT*P);
		it[2]+=h*pT*value_type((*p)[1]*(*p)[2]/P,-(*p)[3]);
		it[3]+=h*pT*value_type((*p)[1]*(*p)[3]/P,(*p)[2]); 

	    }
    };
    template<class model_t>const typename m_vector_wave_function<model_t,4,Minkowski_type,unused>::r_value_type m_vector_wave_function<model_t,4,Minkowski_type,unused>::sqrthlf=std::sqrt((typename model_t::value_type)0.5); 
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_M_VECTOR_WFUNC_H_*/

