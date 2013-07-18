//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_M_SPINOR_FAC_H_
#define CAMGEN_M_SPINOR_FAC_H_

#include <Camgen/spinor_fac.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Massive spinor factory class template declaration and definition. The       *
 * template parameters are the Dirac algebra type, the spin vector class, the  *
 * beam direction integer and the dimension of the spacetime.                  *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* For general dimensions, the massive spinor factory is empty: */

    template<class Dirac_alg_t,class spin_vec_t,int beam_dir,std::size_t Dim>class massive_spinor_factory: public massless_spinor_factory<Dirac_alg_t,beam_dir,Dim>{};

    /* In four dimensions, the specialisation defines all the necessary spinors by
     * applying the formulas
     *
     *                 u(l,p)=N*(pslash + m).(1 + l*g5.sslash)u0(-l)
     *                 v(l,p)=C.ubar(l,p)
     *
     * where N is a suitable normalisation factor, pslash is the momentum slashed, m is
     * the mass, l is the helicity quantum number, g5 is the chiral gamma matrix, sslash
     * is the spin vector slashed and u0(±) are the basic massless spinors, constructed in
     * the massless spinor factory base class. */

    template<class Dirac_alg_t,class spin_vec_t,int beam_dir>class massive_spinor_factory<Dirac_alg_t,spin_vec_t,beam_dir,4>: public massless_spinor_factory<Dirac_alg_t,beam_dir,4>
    {
	public:

	    /* Type definition of the base class: */

	    typedef massless_spinor_factory<Dirac_alg_t,beam_dir,4> base_type;
	    
	    /* The usual type definitions for Feynman rule classes, derived from the base
	     * type: */
	    
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::r_value_type r_value_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename tensor_type::size_type size_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::momentum_type momentum_type;
	    
	    /* Type definitions of the spacetime, Dirac algebra basis and spin
	     * vector types: */
	    
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef Dirac_alg_t Dirac_algebra_type;
	    typedef typename spin_vec_t::template implementation<spacetime_type,beam_dir> spin_vec_type;
	
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
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 - sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_second(-c,s,temp.begin(),base_type::u_m.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_m[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_second(value_type(1,0),*p,it,temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=(*m)*temp[i];
		}
	    }
	    
	    /* Positive-helicity spinor construction with prefactor h: necessary
	     * when dealing with continuous colours:*/

	    static void make_u_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 - sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_second(-c,s,temp.begin(),base_type::u_m.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_m[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_second(h,*p,it,temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=h*(*m)*temp[i];
		}
	    }
	    
	    /* Positive-helicity row spinor construction: */
	    
	    static void make_u_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 - sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_third(-c,s,base_type::u_m_bar.begin(),temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_m_bar[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_third(value_type(1,0),*p,temp.begin(),it);
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=(*m)*temp[i];
		}
	    }

	    /* Positive-helicity row spinor construction with prefactor: */
	    
	    static void make_u_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 - sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_third(-c,s,base_type::u_m_bar.begin(),temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_m_bar[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_third(h,*p,temp.begin(),it);
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=h*(*m)*temp[i];
		}
	    }

	    /* Negative-helicity spinor construction: */
	    
	    static void make_u_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 + sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_second(c,s,temp.begin(),base_type::u_p.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_p[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_second(value_type(1,0),*p,it,temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=(*m)*temp[i];
		}
	    }

	    /* Negative-helicity spinor construction with prefactor: */
	    
	    static void make_u_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 + sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_second(c,s,temp.begin(),base_type::u_p.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_p[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_second(h,*p,it,temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=h*(*m)*temp[i];
		}
	    }

	    /* Negative-helicity row-spinor construction: */
	    
	    static void make_u_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 + sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_third(c,s,base_type::u_p_bar.begin(),temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_p_bar[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_third(value_type(1,0),*p,temp.begin(),it);
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=(*m)*temp[i];
		}
	    }

	    /* Negative-helicity row-spinor construction with prefactor: */
	    
	    static void make_u_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		/* Create spin vector: */
		
		spin_vec_type::make_spinvec(s,p,m);

		/* Compute normalisation constant: */

		value_type c=spacetime_type::dot(base_type::k_0,*p) + (*m)*(spacetime_type::dot(base_type::k_0,s));
		c=(r_value_type)0.5/std::sqrt(c);

		/* Multiplying the appropriate massless base spinor with (1 + sslash.g5)
		 * and storing the result in temp: */

		temp.reset();
		Dirac_alg_t::gg5_third(c,s,base_type::u_p_bar.begin(),temp.begin());
		for(size_type i=0;i<4;++i)
		{
		    temp[i]+=c*(base_type::u_p_bar[i]);
		}

		/* Multiplying with (pslash + m) to the iterator: */

		Dirac_alg_t::g_third(h,*p,temp.begin(),it);
		for(size_type i=0;i<4;++i)
		{
		    it[i]+=h*(*m)*temp[i];
		}
	    }

	    /* Positive-helicity anti-spinor construction: */

	    static void make_v_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_p_bar(temp2.begin(),p,m);
		Dirac_algebra_type::C_left(it,temp2.begin());
	    }

	    /* Positive-helicity anti-spinor construction with prefactor: */

	    static void make_v_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_p_bar(h,temp2.begin(),p,m);
		Dirac_algebra_type::C_left(it,temp2.begin());
	    }

	    /* Positive-helicity row-anti-spinor construction: */

	    static void make_v_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_p(temp2.begin(),p,m);
		Dirac_algebra_type::Cc_right(it,temp2.begin());
	    }

	    /* Positive-helicity row-anti-spinor construction with prefactor: */

	    static void make_v_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_p(h,temp2.begin(),p,m);
		Dirac_algebra_type::Cc_right(it,temp2.begin());
	    }

	    /* Negative-helicity anti-spinor construction: */
	    
	    static void make_v_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_m_bar(temp2.begin(),p,m);
		Dirac_algebra_type::C_left(it,temp2.begin());
	    }

	    /* Negative-helicity anti-spinor construction with prefactor: */
	    
	    static void make_v_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_m_bar(h,temp2.begin(),p,m);
		Dirac_algebra_type::C_left(it,temp2.begin());
	    }

	    /* Negative-helicity row-anti-spinor construction: */
	    
	    static void make_v_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_m(temp2.begin(),p,m);
		Dirac_algebra_type::Cc_right(it,temp2.begin());
	    }

	    /* Negative-helicity row-anti-spinor construction with prefactor: */
	    
	    static void make_v_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp2.reset();
		make_u_m(h,temp2.begin(),p,m);
		Dirac_algebra_type::Cc_right(it,temp2.begin());
	    }
	private:
	    
	    /* temporary spin vector: */

	    static momentum_type s;
	    
	    /* Temporary spinors: */
	    
	    static tensor_type temp;
	    static tensor_type temp2;
    };
    template<class Dirac_alg_t,class spin_vec_t,int beam_dir>typename massive_spinor_factory<Dirac_alg_t,spin_vec_t,beam_dir,4>::momentum_type massive_spinor_factory<Dirac_alg_t,spin_vec_t,beam_dir,4>::s;
    template<class Dirac_alg_t,class spin_vec_t,int beam_dir>typename massive_spinor_factory<Dirac_alg_t,spin_vec_t,beam_dir,4>::tensor_type massive_spinor_factory<Dirac_alg_t,spin_vec_t,beam_dir,4>::temp(1,4);
    template<class Dirac_alg_t,class spin_vec_t,int beam_dir>typename massive_spinor_factory<Dirac_alg_t,spin_vec_t,beam_dir,4>::tensor_type massive_spinor_factory<Dirac_alg_t,spin_vec_t,beam_dir,4>::temp2(1,4);
}

/* Specialisations in the case of the Pauli and Weyl bases combined with standard spin
 * vector choices: */

#ifdef CAMGEN_PAULI_BASIS_H_

#ifdef CAMGEN_KS_TYPE_H_
#include <Camgen/KSspinPb.h>
#endif /*CAMGEN_KS_TYPE_H_*/

#ifdef CAMGEN_HEL_TYPE_H_
#include <Camgen/helspinPb.h>
#endif /*CAMGEN_HELICITY_TYPE_H_*/

#ifdef CAMGEN_POL_TYPE_H_
#include <Camgen/polspinPb.h>
#endif /*CAMGEN_POL_TYPE_H_*/

#endif /*CAMGEN_PAULI_BASIS_H_*/

#ifdef CAMGEN_WEYL_BASIS_H_

#ifdef CAMGEN_KS_TYPE_H_
#include <Camgen/KSspinWb.h>
#endif /*CAMGEN_KS_TYPE_H_*/

#ifdef CAMGEN_HEL_TYPE_H_
#include <Camgen/helspinWb.h>
#endif /*CAMGEN_HELICITY_TYPE_H_*/

#ifdef CAMGEN_POL_TYPE_H_
#include <Camgen/polspinWb.h>
#endif /*CAMGEN_POL_TYPE_H_*/

#endif /*CAMGEN_WEYL_BASIS_H_*/

#endif /*CAMGEN_M_SPINOR_FAC_H_*/

