//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Specialisation of the massless spinor factory class in the Weyl basis of      *
 * Dirac gamma matrices. This header cannot be compiled stand-alone, instead it  *
 * should be included through the files Weyl_basis.h or spinor_fac.h.            *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t>class massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>
    {
	public:
	    
	    /* The usual type definitions for Feynman rule classes, derived from the base
	     * type: */
	    
	    typedef value_t r_value_type;
	    typedef std::complex<value_t> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::size_type size_type;
	    typedef typename tensor_type::iterator iterator;
	    
	    /* Type definitions of the spacetime, Dirac algebra basis and spin
	     * vector types: */
	    
	    typedef typename Weyl_basis::template implementation<r_value_type,4,true> Dirac_algebra_type;
	    typedef typename Minkowski_type::template implementation<r_value_type,4> spacetime_type;
	    typedef vector<r_value_type,4> momentum_type;
	    
	    /* Initialisation function: */

	    static void initialise()
	    {
		if(!initialised)
		{

		    Dirac_algebra_type::initialise();
		    spacetime_type::template make_base_vectors<3>(k_0,k_1,k_2);
		    value_type c(1,0);
		    
		    u_p[2]=c;
		    u_p[3]=c;
		    u_m[0]=c;
		    u_m[1]=-c;

		    u_p_bar[0]=c;
		    u_p_bar[1]=c;
		    u_m_bar[2]=c;
		    u_m_bar[3]=-c;
		    
		    initialised=true;

		}
	    }

	    /* Testing whether the class is initialised: */

	    static bool is_initialised()
	    {
		return initialised;
	    }

	    /* Positive-helicity massless spinor construction: */

	    static void make_u_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type c=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[0]-(*p)[1]));

		it[2]+=c*value_type((*p)[0]-(*p)[1]+(*p)[3],(*p)[2]);
		it[3]+=c*value_type(-(*p)[0]+(*p)[1]+(*p)[3],(*p)[2]);	
		
	    }

	    /* Positive-helicity massless spinor construction with prefactor: */

	    static void make_u_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		value_type c=h/std::sqrt((r_value_type)2*((*p)[0]-(*p)[1]));
		
		it[2]+=c*value_type((*p)[0]-(*p)[1]+(*p)[3],(*p)[2]);
		it[3]+=c*value_type(-(*p)[0]+(*p)[1]+(*p)[3],(*p)[2]);	
		
	    }

	    /* Positive-helicity massless row-spinor construction: */

	    static void make_u_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type c=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[0]-(*p)[1]));
		
		it[0]+=c*value_type((*p)[0]-(*p)[1]+(*p)[3],-(*p)[2]);
		it[1]+=c*value_type(-(*p)[0]+(*p)[1]+(*p)[3],-(*p)[2]);
		
	    }

	    /* Positive-helicity massless row-spinor construction with prefactor: */

	    static void make_u_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		value_type c=h/std::sqrt((r_value_type)2*((*p)[0]-(*p)[1]));
		
		it[0]+=c*value_type((*p)[0]-(*p)[1]+(*p)[3],-(*p)[2]);
		it[1]+=c*value_type(-(*p)[0]+(*p)[1]+(*p)[3],-(*p)[2]);
		
	    }

	    /* Negative-helicity massless spinor construction: */

	    static void make_u_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type c=(r_value_type)1/std::sqrt(2*((*p)[0]-(*p)[1]));
		
		it[0]-=c*value_type((*p)[0]-(*p)[1]-(*p)[3],(*p)[2]);
		it[1]-=c*value_type((*p)[0]-(*p)[1]+(*p)[3],-(*p)[2]);
		
	    }

	    /* Negative-helicity massless spinor construction with prefactor: */

	    static void make_u_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		value_type c=h/std::sqrt((r_value_type)2*((*p)[0]-(*p)[1]));

		it[0]-=c*value_type((*p)[0]-(*p)[1]-(*p)[3],(*p)[2]);
		it[1]-=c*value_type((*p)[0]-(*p)[1]+(*p)[3],-(*p)[2]);
		
	    }

	    /* Negative-helicity massless row-spinor construction: */

	    static void make_u_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type c=(r_value_type)1/std::sqrt(2*((*p)[0]-(*p)[1]));
		
		it[2]-=c*value_type((*p)[0]-(*p)[1]-(*p)[3],-(*p)[2]);
		it[3]-=c*value_type((*p)[0]-(*p)[1]+(*p)[3],(*p)[2]);
		
	    }

	    /* Negative-helicity massless row-spinor construction with prefactor: */

	    static void make_u_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		value_type c=h/std::sqrt((r_value_type)2*((*p)[0]-(*p)[1]));

		it[2]-=c*value_type((*p)[0]-(*p)[1]-(*p)[3],-(*p)[2]);
		it[3]-=c*value_type((*p)[0]-(*p)[1]+(*p)[3],(*p)[2]);
		
	    }

	    /* Longitudinal massless polarisation vector construction: */

	    static void make_e_L(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_MESSAGE("longitudinal massless polarisation vector encountered");
		return;
	    }

	    /* Positive-helicity massless polarisation vector construction: */

	    static void make_e_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1] or (*p)[0]==(*p)[2]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type a=(r_value_type)0.5/std::sqrt(((*p)[0]-(*p)[1])*((*p)[0]-(*p)[2]));

		it[0]+=value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),a*(*p)[3]);
		it[1]+=value_type(a*(-(*p)[0]+(*p)[1]+(*p)[2]),a*(*p)[3]);
		it[2]+=value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),a*(*p)[3]);
		it[3]+=value_type(a*(*p)[3],a*((*p)[0]-(*p)[1]-(*p)[2]));
		
	    }

	    /* Positive-helicity massless polarisation vector construction with prefactor: */

	    static void make_e_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1] or (*p)[0]==(*p)[2]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type a=(r_value_type)0.5/std::sqrt(((*p)[0]-(*p)[1])*((*p)[0]-(*p)[2]));
		
		it[0]+=h*value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),a*(*p)[3]);
		it[1]+=h*value_type(a*(-(*p)[0]+(*p)[1]+(*p)[2]),a*(*p)[3]);
		it[2]+=h*value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),a*(*p)[3]);
		it[3]+=h*value_type(a*(*p)[3],a*((*p)[0]-(*p)[1]-(*p)[2]));
		
	    }

	    /* Longitudinal massless polarisation vector construction with prefactor: */

	    static void make_e_L(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m){}

	    /* Negative-helicity massless polarisation vector construction: */

	    static void make_e_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1] or (*p)[0]==(*p)[2]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type a=(r_value_type)0.5/std::sqrt(((*p)[0]-(*p)[1])*((*p)[0]-(*p)[2]));

		it[0]+=value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),-a*(*p)[3]);
		it[1]+=value_type(a*(-(*p)[0]+(*p)[1]+(*p)[2]),-a*(*p)[3]);
		it[2]+=value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),-a*(*p)[3]);
		it[3]+=value_type(a*(*p)[3],a*(-(*p)[0]+(*p)[1]+(*p)[2]));
		
	    }

	    /* Negative-helicity massless polarisation vector construction with prefactor: */

	    static void make_e_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF(((*p)[0]==(*p)[1] or (*p)[0]==(*p)[2]),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");
		
		r_value_type a=(r_value_type)0.5/std::sqrt(((*p)[0]-(*p)[1])*((*p)[0]-(*p)[2]));

		it[0]+=h*value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),-a*(*p)[3]);
		it[1]+=h*value_type(a*(-(*p)[0]+(*p)[1]+(*p)[2]),-a*(*p)[3]);
		it[2]+=h*value_type(a*((*p)[0]-(*p)[1]+(*p)[2]),-a*(*p)[3]);
		it[3]+=h*value_type(a*(*p)[3],a*(-(*p)[0]+(*p)[1]+(*p)[2]));
		
	    }

	    /* Positive-helicity massive polarisation vector construction: */

	    static void make_m_e_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		
		if(P==(r_value_type)0)
		{
		    it[1]=value_type(sqrthlf,0);
		    it[2]=value_type(0,sqrthlf);
		    return;
		}

		CAMGEN_MESSAGE_IF(((*p)[2]==0 and (*p)[3]==0),"inappropriate spinor base for given momentum");
		
		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));

		it[1]=-(r_value_type)0.5/(pT*P);
		it[2]=pT*value_type((*p)[1]*(*p)[2]/P,(*p)[3]);
		it[3]=pT*value_type((*p)[1]*(*p)[3]/P,-(*p)[2]); 
		
	    }

	    /* Positive-helicity massive polarisation vector construction with prefactor: */

	    static void make_m_e_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");

		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		
		if(P==(r_value_type)0)
		{
		    r_value_type sqrthlf=std::sqrt((r_value_type)0.5);
		    it[1]+=value_type(sqrthlf,0)*h;
		    it[2]+=value_type(0,sqrthlf)*h;
		    return;
		}

		CAMGEN_MESSAGE_IF(((*p)[2]==0 and (*p)[3]==0),"inappropriate spinor base for given momentum");
		
		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));

		it[1]-=(r_value_type)0.5*h/(pT*P);
		it[2]+=h*pT*value_type((*p)[1]*(*p)[2]/P,(*p)[3]);
		it[3]+=h*pT*value_type((*p)[1]*(*p)[3]/P,-(*p)[2]); 
		
	    }

	    /* Longitudinal polarisation vector construction: */

	    static void make_m_e_L(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		if(P==(r_value_type)0)
		{
		    it[3]=value_type(1,0);
		    return;
		}

		CAMGEN_MESSAGE_IF((*m==0),"zero-valued mass encountered");

		r_value_type xi=(*p)[0]/((*m)*P);

		it[0]=P/(*m);
		it[1]=xi*(*p)[1];
		it[2]=xi*(*p)[2];
		it[3]=xi*(*p)[3];
		
	    }

	    /* Longitudinal polarisation vector construction with prefactor: */

	    static void make_m_e_L(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		
		if(P==(r_value_type)0)
		{
		    it[3]+=h;
		    return;
		}
		
		CAMGEN_MESSAGE_IF((*m==0),"zero-valued mass encountered");

		value_type xi=(*p)[0]*h/((*m)*P);
		
		it[0]+=P*h/(*m);
		it[1]+=xi*(*p)[1];
		it[2]+=xi*(*p)[2];
		it[3]+=xi*(*p)[3];
		
	    }

	    /* Negative-helicity massive polarisation vector construction: */

	    static void make_m_e_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));

		if(P==(r_value_type)0)
		{
		    r_value_type sqrthlf=std::sqrt((r_value_type)0.5);
		    it[1]=value_type(sqrthlf,0);
		    it[2]=value_type(0,-sqrthlf);
		    return;
		}

		CAMGEN_MESSAGE_IF(((*p)[2]==0 and (*p)[3]==0),"inappropriate spinor base for given momentum");
		
		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));
		
		it[1]=-(r_value_type)0.5/(pT*P);
		it[2]=pT*value_type((*p)[1]*(*p)[2]/P,-(*p)[3]);
		it[3]=pT*value_type((*p)[1]*(*p)[3]/P,(*p)[2]); 
		
	    }

	    /* Negative-helicity massive polarisation vector construction with prefactor: */

	    static void make_m_e_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((it.range()<4),"tensor iterator out of range");
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		
		r_value_type P=std::sqrt((*p)[0]*(*p)[0]-(*m)*(*m));
		
		if(P==(r_value_type)0)
		{
		    r_value_type sqrthlf=std::sqrt((r_value_type)0.5);
		    it[1]+=value_type(sqrthlf,0)*h;
		    it[2]-=value_type(0,sqrthlf)*h;
		    return;
		}

		CAMGEN_MESSAGE_IF(((*p)[2]==0 and (*p)[3]==0),"inappropriate spinor base for given momentum");
		
		r_value_type pT=(r_value_type)1/std::sqrt((r_value_type)2*((*p)[2]*(*p)[2]+(*p)[3]*(*p)[3]));
		
		it[1]-=(r_value_type)0.5*h/(pT*P);
		it[2]+=h*pT*value_type((*p)[1]*(*p)[2]/P,-(*p)[3]);
		it[3]+=h*pT*value_type((*p)[1]*(*p)[3]/P,(*p)[2]); 
		
	    }

	    /* Output method: */

	    static std::ostream& print(std::ostream& os)
	    {
		os<<"base vectors:"<<k_0<<","<<k_1<<","<<k_2<<std::endl;
		os<<"base spinors:"<<std::endl<<std::endl;
		for(size_type i=0;i<4;++i)
		{
		    os<<"\t\t\t";
		    os<<std::left<<"[";
		    os<<std::setw(2);
		    shortprint(os,u_p[i]);
		    os<<"]"<<"\t\t\t";
		    os<<std::left<<"[";
		    os<<std::setw(2);
		    shortprint(os,u_m[i]);
		    os<<"]"<<std::endl;
		}
		os<<"base row spinors:"<<std::endl;
		os<<u_p_bar<<"        "<<u_m_bar<<std::endl;
		return os;
	    }

	    /* Checking function: */

	    static void check_base_spinors()
	    {

		if(!initialised)
		{
		    initialise();
		}
		value_type ap;
		value_type am;
		value_type bp;
		value_type bm;
		for(size_type i=0;i<4;++i)
		{
		    for(size_type j=0;j<4;++j)
		    {
			ap=u_p[i]*u_p_bar[j];
			am=u_m[i]*u_m_bar[j];
			bp=0;
			bm=0;
			for(size_type mu=0;mu<4;++mu)
			{
			    for(size_type nu=0;nu<4;++nu)
			    {
				bp+=Dirac_algebra_type::ggamma_L(mu,i,j)*spacetime_type::metric(mu,nu)*k_0[nu];
				bm+=Dirac_algebra_type::ggamma_R(mu,i,j)*spacetime_type::metric(mu,nu)*k_0[nu];
			    }
			}
			bool q=false;
			if(!equals(ap,bp))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"positive-helicity base spinors do not multiply to gL.slash(k_0)"<<endlog;
			    q=true;
			}
			if(!equals(am,bm))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"negative-helicity base spinors do not multiply to gR.slash(k_0)"<<endlog;
			    q=true;
			}
			if(q)
			{
			    return;
			}
		    }
		}

	    }

	    /* Spinor product s(p,q) = bar{u}(+,p).u(-,q): */

	    static value_type s(const momentum_type& p,const momentum_type& q)
	    {
		
		tensor_type s1(1,4);
		make_u_p_bar(s1.begin(),&p,NULL);
		tensor_type s2(1,4);
		make_u_m(s2.begin(),&q,NULL);
		
		
		return s1[0]*s2[0]+s1[1]*s2[1]+s1[2]*s2[2]+s1[3]*s2[3];
	    }

	    /* Spinor product t(p,q) = bar{u}(-,p).u(+,q): */

	    static value_type t(const momentum_type& p,const momentum_type& q)
	    {
		
		tensor_type s1(1,4);
		make_u_m_bar(s1.begin(),&p,NULL);
		tensor_type s2(1,4);
		make_u_p(s2.begin(),&q,NULL);
		
		
		return s1[0]*s2[0]+s1[1]*s2[1]+s1[2]*s2[2]+s1[3]*s2[3];
	    }
	
	protected:
	    
	    /* Protected static data members: */
	    
	    /* Base spinors: */
	    
	    static tensor_type u_p;
	    static tensor_type u_p_bar;
	    static tensor_type u_m;
	    static tensor_type u_m_bar;
	    
	    /* Secondary base spinors for massless polarisation vectors: */
	    
	    static tensor_type u_2_p;
	    static tensor_type u_2_m;
	    
	    /* Base momenta: */
	    
	    static momentum_type k_0;
	    static momentum_type k_1;
	    static momentum_type k_2;
	
	private:

	    /* Private static data members: */

	    /* Initialisation tag: */

	    static bool initialised;

	    /* temporary spinors: */

	    static tensor_type temp1;
	    static tensor_type temp2;

	    /* Useful constants: */

	    static const value_t sqrthlf;
    };

    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::tensor_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::u_p(1,4);
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::tensor_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::u_m(1,4);
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::tensor_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::u_p_bar(1,4);
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::tensor_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::u_m_bar(1,4);
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::tensor_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::u_2_p(1,4);
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::tensor_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::u_2_m(1,4);
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::momentum_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::k_0;
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::momentum_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::k_1;
    template<class value_t>typename massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::momentum_type massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::k_2;
    template<class value_t>bool massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::initialised=false;
    template<class value_t>const value_t massless_spinor_factory<typename Weyl_basis::template implementation<value_t,4,true>,3,4>::sqrthlf=std::sqrt((value_t)0.5);
}

