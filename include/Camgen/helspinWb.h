//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Specialisation of the massive helicity spinor factory in the Weyl basis of  *
 * gamma matrices. This header cannot be included stand-alone, it can only be  *
 * included through the files Weyl_basis.h or helicity_type.h                  *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t>class massive_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,helicity_type,3,4>: public massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>
    {
	public:
	    
	    /* Type definition of the parent class: */
	    
	    typedef massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4> base_type;
	    
	    /* Usual type definitions in Feynman rule classes: */
	    
	    typedef value_t r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::momentum_type momentum_type;
	    
	    /* Type definitions of the spacetime, Dirac algebra basis and spin
	     * vector types: */
	    
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::Dirac_algebra_type Dirac_algebra_type;
	    typedef typename helicity_type::template implementation<spacetime_type,3> spin_vec_type;

	    /* No private data members, so only the base class initialises: */
	    
	    static void initialise()
	    {
		base_type::initialise();
	    }

	    /* Positive-helicity spinor construction: it is the iterator in the
	     * tensor that will be filled ny the spinor wave function, p is a
	     * pointer to the momentum and m a pointer to the mass: */

	    static void make_u_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		value_type b(c*(-P+(*p)[1]+(*p)[3]),c*(*p)[2]);
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]+=((*p)[0]+P)*a;
		it[3]+=((*p)[0]+P)*b;
	    }

	    /* Positive-helicity spinor construction with prefactor h: necessary
	     * when dealing with continuous colours:*/

	    static void make_u_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		value_type b(c*(-P+(*p)[1]+(*p)[3]),c*(*p)[2]);
		a*=h;
		b*=h;
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]+=((*p)[0]+P)*a;
		it[3]+=((*p)[0]+P)*b;
	    }

	    /* Positive-helicity row spinor construction: */

	    static void make_u_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		value_type b(c*(-P+(*p)[1]+(*p)[3]),-c*(*p)[2]);
		
		it[0]+=((*p)[0]+P)*a;
		it[1]+=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }

	    /* Positive-helicity row spinor construction with prefactor: */

	    static void make_u_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		std::complex<r_value_type>a(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		std::complex<r_value_type>b(c*(-P+(*p)[1]+(*p)[3]),-c*(*p)[2]);
		a*=h;
		b*=h;
		
		it[0]+=((*p)[0]+P)*a;
		it[1]+=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }

	    /* negative-helicity spinor construction: */

	    static void make_u_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		
		it[0]+=((*p)[0]+P)*a;
		it[1]+=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }

	    /* Negative-helicity spinor construction with prefactor: */

	    static void make_u_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		a*=h;
		b*=h;
		
		it[0]+=((*p)[0]+P)*a;
		it[1]+=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }

	    /* Negative-helicity row-spinor construction: */

	    static void make_u_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),-c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]+=((*p)[0]+P)*a;
		it[3]+=((*p)[0]+P)*b;
	    }

	    /* Negative-helicity row-spinor construction with prefactor: */

	    static void make_u_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),-c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		a*=h;
		b*=h;
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]+=((*p)[0]+P)*a;
		it[3]+=((*p)[0]+P)*b;
	    }

	    /* Positive-helicity anti-spinor construction: */

	    static void make_v_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		
		it[0]-=((*p)[0]+P)*a;
		it[1]-=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }

	    /* Positive-helicity anti-spinor construction with prefactor: */

	    static void make_v_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		a*=h;
		b*=h;
		
		it[0]-=((*p)[0]+P)*a;
		it[1]-=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }

	    /* Positive-helicity row-anti-spinor construction: */

	    static void make_v_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),-c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]-=((*p)[0]+P)*a;
		it[3]-=((*p)[0]+P)*b;
	    }

	    /* Positive-helicity row-anti-spinor construction with prefactor: */

	    static void make_v_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]-(*p)[3]),-c*(*p)[2]);
		value_type b(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		a*=h;
		b*=h;
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]-=((*p)[0]+P)*a;
		it[3]-=((*p)[0]+P)*b;
	    }

	    /* Negative-helicity anti-spinor construction: */

	    static void make_v_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		value_type b(c*(-P+(*p)[1]+(*p)[3]),c*(*p)[2]);
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]-=((*p)[0]+P)*a;
		it[3]-=((*p)[0]+P)*b;
	    }

	    /* Negative-helicity anti-spinor construction with prefactor: */

	    static void make_v_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]+(*p)[3]),c*(*p)[2]);
		value_type b(c*(-P+(*p)[1]+(*p)[3]),c*(*p)[2]);
		a*=h;
		b*=h;
		
		it[0]+=(*m)*a;
		it[1]+=(*m)*b;
		it[2]-=((*p)[0]+P)*a;
		it[3]-=((*p)[0]+P)*b;
	    }

	    /* Negative-helicity row-anti-spinor construction: */

	    static void make_v_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		value_type b(c*(-P+(*p)[1]+(*p)[3]),-c*(*p)[2]);

		it[0]-=((*p)[0]+P)*a;
		it[1]-=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }

	    /* Negative-helicity row-anti-spinor construction with prefactor: */

	    static void make_v_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p == NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m == NULL),"NULL mass pointer encountered in");

		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		CAMGEN_MESSAGE_IF((P==0 or P==-(*p)[0] or P==(*p)[1]),"inappropriate spinor base for given momentum");

		r_value_type c=std::sqrt((r_value_type)0.25/(P*((*p)[0]+P)*(P-(*p)[1])));
		value_type a(c*(P-(*p)[1]+(*p)[3]),-c*(*p)[2]);
		value_type b(c*(-P+(*p)[1]+(*p)[3]),-c*(*p)[2]);
		a*=h;
		b*=h;
		it[0]-=((*p)[0]+P)*a;
		it[1]-=((*p)[0]+P)*b;
		it[2]+=(*m)*a;
		it[3]+=(*m)*b;
	    }
    };
}
