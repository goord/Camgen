//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Specialisation of the massive polarised spinor factory in the Pauli basis of  *
 * gamma matrices. This header cannot be included stand-alone, it can only be    *
 * included through the files Pauli_basis.h or polarised_type.h                  *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t>class massive_spinor_factory<typename Pauli_basis::template implementation<value_t,4,true>,polarised_type,3,4>: public massless_spinor_factory<typename Pauli_basis::template implementation<value_t,4,true>,3,4>
    {
	public:

	    /* Type definition of the parent class: */

	    typedef massless_spinor_factory<typename Pauli_basis::template implementation<value_t,4,true>,3,4> base_type;
	    
	    /* Usual type definitions for Feynman rule classes: */
	    
	    typedef value_t r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::momentum_type momentum_type;

	    /* Type definitions of the spacetime type, Dirac algebra type and
	     * spin vector type: */

	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::Dirac_algebra_type Dirac_algebra_type;
	    typedef typename polarised_type::template implementation<spacetime_type,3> spin_vec_type;

	    /* No private data members, so only the base class initialises: */
	    
	    static void initialise()
	    {
		base_type::initialise();
	    }

	    /* Positive-helicity spinor construction: it is the iterator in the
	     * tensor that will be filled with the spinor wave function, p is a
	     * pointer to the momentum and m a pointer to the mass: */

	    static void make_u_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]+=c*e;
		it[2]+=c*(*p)[3];
		it[3]+=c*value_type((*p)[1],(*p)[2]);
	    }
	    
	    /* Positive-helicity spinor construction with prefactor h: necessary
	     * when dealing with continuous colours:*/

	    static void make_u_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]+=c*e;
		it[2]+=c*(*p)[3];
		it[3]+=c*value_type((*p)[1],(*p)[2]);
	    }
	    
	    /* Positive-helicity row spinor construction: */
	    
	    static void make_u_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]+=c*e;
		it[2]-=c*(*p)[3];
		it[3]+=c*value_type(-(*p)[1],(*p)[2]);
	    }
	    
	    /* Positive-helicity row spinor construction with prefactor: */
	    
	    static void make_u_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]+=c*e;
		it[2]-=c*(*p)[3];
		it[3]+=c*value_type(-(*p)[1],(*p)[2]);
	    }

	    /* Negative-helicity spinor construction: */
	    
	    static void make_u_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[1]+=c*e;
		it[2]+=c*value_type((*p)[1],-(*p)[2]);
		it[3]-=c*(*p)[3];
	    }

	    /* Negative-helicity spinor construction with prefactor: */
	    
	    static void make_u_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[1]+=c*e;
		it[2]+=c*value_type((*p)[1],-(*p)[2]);
		it[3]-=c*(*p)[3];
	    }

	    /* Negative-helicity row-spinor construction: */
	    
	    static void make_u_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[1]+=c*e;
		it[2]-=c*value_type((*p)[1],(*p)[2]);
		it[3]+=c*(*p)[3];
	    }

	    /* Negative-helicity row-spinor construction with prefactor: */
	    
	    static void make_u_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[1]+=c*e;
		it[2]-=c*value_type((*p)[1],(*p)[2]);
		it[3]+=c*(*p)[3];
	    }

	    /* Positive-helicity anti-spinor construction: */

	    static void make_v_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]-=c*value_type(-(*p)[1],(*p)[2]);
		it[1]-=c*(*p)[3];
		it[3]+=c*e;
	    }

	    /* Positive-helicity anti-spinor construction with prefactor: */

	    static void make_v_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]-=c*value_type(-(*p)[1],(*p)[2]);
		it[1]-=c*(*p)[3];
		it[3]+=c*e;
	    }

	    /* Positive-helicity row-anti-spinor construction: */

	    static void make_v_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]+=c*value_type((*p)[1],(*p)[2]);
		it[1]-=c*(*p)[3];
		it[3]-=c*e;
	    }

	    /* Positive-helicity row-anti-spinor construction with prefactor: */

	    static void make_v_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]+=c*value_type((*p)[1],(*p)[2]);
		it[1]-=c*(*p)[3];
		it[3]-=c*e;
	    }

	    /* Negative-helicity anti-spinor construction: */

	    static void make_v_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]-=c*(*p)[3];
		it[1]-=c*value_type((*p)[1],(*p)[2]);
		it[2]-=c*e;
	    }

	    /* Negative-helicity anti-spinor construction with prefactor: */

	    static void make_v_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]-=c*(*p)[3];
		it[1]-=c*value_type((*p)[1],(*p)[2]);
		it[2]-=c*e;
	    }

	    /* Negative-helicity row-anti-spinor construction: */

	    static void make_v_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type e=(*p)[0]+*m;
		value_type c=value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]-=c*(*p)[3];
		it[1]-=c*value_type((*p)[1],-(*p)[2]);
		it[2]+=c*e;
	    }

	    /* Negative-helicity row-anti-spinor construction with prefactor: */

	    static void make_v_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");

		r_value_type e=(*p)[0]+*m;
		value_type c=h*value_type(e-(*p)[1]+(*p)[3],-(*p)[2])/std::sqrt((r_value_type)2*e*(e*((*p)[0]-(*p)[1]+(*p)[3])-(*p)[1]*(*p)[3]));
		it[0]-=c*(*p)[3];
		it[1]-=c*value_type((*p)[1],-(*p)[2]);
		it[2]+=c*e;
	    }
    };
}

