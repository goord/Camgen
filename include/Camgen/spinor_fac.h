//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SPINOR_FAC_H_
#define CAMGEN_SPINOR_FAC_H_

#include <Camgen/forward_decs.h>
#include <Camgen/debug.h>
#include <Camgen/utils.h>
#include <Camgen/tensor.h>
#include <Camgen/vector.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Massless spinor factory class definition. It constructs the massless basic  *
 * spinors, used to construct massless spinor wave functions, (massive)        *
 * polarisation vectors and massive spinors (by the massive spinor factory     *
 * subclasses). The template parameters are the Dirac algebra type, the beam   *
 * direction integer and the spacetime dimensions.                             *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* For general dimensions, the massless spinor factory is empty: */

    template<class Dirac_alg_t,int beam_dir,std::size_t D>class massless_spinor_factory;

    /* In four dimensions, the base spinors are constructed fulfiling the
     * equations
     *
     * 			u0(-)bar{u0}(-) = omega(-).slash{k_0}
     * 			          u0(+) = slash(k_1)u0(-)
     *
     * where omega(±) = 1 ± g5 and k_0 and k_1 are the base momentum constructed
     * by the spacetime class. Then massless spinors are constructed by
     *
     * 		      u(l,p) = (1/sqrt(2p.k_0))*slash{p}.u0(-l)
     *
     * The massive polarisation vectors use the momentum-splitting function
     * defined in the spacetime class. Let p_1 and p_2 the momenta p splits
     * into. Then the factory uses the massive helicity vectors
     *
     * 		 e(±,p)(mu) = (1/(sqrt{2}*m)*bar{u(±,p_1)g(mu)u(±,p_2)
     * 		 		 e(0,p) = (p_1-p_2)/m
     *
     * For the massless helicity vectors, a third massless base vector k_2 is
     * constructed, not parallel to the momentum, and yields
     *
     * 		 e(±,p)(mu) = (±1/s(±,p,k_2))*bar{u}(±,p).g(mu).u(±,k_2)
     * 
     * where s(l,p,q)=bar{u}(l,p).u(-l,q) is the spinor product. */ 		 		 

    template<class Dirac_alg_t,int beam_dir>class massless_spinor_factory<Dirac_alg_t,beam_dir,4>
    {
	public:
	    
	    /* The usual type definitions for Feynman rule classes, derived from the base
	     * type: */
	    
	    typedef typename Dirac_alg_t::base_type::r_value_type r_value_type;
	    typedef typename Dirac_alg_t::base_type::value_type value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::size_type size_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef Dirac_alg_t Dirac_algebra_type;
	    
	    /* Type definitions of the spacetime, Dirac algebra basis and spin
	     * vector types: */
	    
	    typedef typename Dirac_alg_t::spacetime_type spacetime_type;
	    typedef vector<r_value_type,4> momentum_type;

	    /* Initialisation function: */

	    static void initialise()
	    {
		if(!initialised)
		{
		    
		    /* If not done yet, construct Dirac matrices: */

		    Dirac_algebra_type::initialise();

		    /* Construct the three massless base vectors for the
		     * massless spinors and helicity vectors: */

		    spacetime_type::template make_base_vectors<beam_dir>(k_0,k_1,k_2);
		    
		    /* Compute left-hand side of the equation defining a basic
		     * massless spinor (M denotes the Dirac conjugation matrix): */

		    value_type kslash[4][4];
		    if(spacetime_type::diagonal)
		    {
			for(size_type i=0;i<4;++i)
			{
			    for(size_type j=0;j<4;++j)
			    {
				kslash[i][j]=0;
				for(size_type mu=0;mu<4;++mu)
				{
				    for(size_type k=0;k<4;++k)
				    {
					kslash[i][j]+=(spacetime_type::metric(mu,mu)*k_0[mu]*Dirac_algebra_type::gamma(mu,i,k)*Dirac_algebra_type::M(k,j));
				    }
				}
			    }
			}
		    }
		    else
		    {
			for(size_type i=0;i<4;++i)
			{
			    for(size_type j=0;j<4;++j)
			    {
				kslash[i][j]=0;
				for(size_type mu=0;mu<4;++mu)
				{
				    for(size_type nu=0;nu<4;++nu)
				    {
					for(size_type k=0;k<4;++k)
					{
					    kslash[i][j]+=(spacetime_type::metric(mu,nu)*k_0[mu]*Dirac_algebra_type::gamma(nu,i,k)*Dirac_algebra_type::M(k,j));
					}
				    }
				}
			    }
			}
		    }

		    /* Compute the positive-helicity basic spinor: */

		    for(size_type i=0;i<4;++i)
		    {
			u_p[i]=0;
			for(size_type j=0;j<4;++j)
			{
			    u_p[i]+=(Dirac_algebra_type::gamma_L(i,j)*kslash[j][i]);
			}
			u_p[i]=std::sqrt(u_p[i]);
			if(u_p[i] != (value_type)0)
			{
			    for(size_type k=i+1;k<4;++k)
			    {
				for(size_type j=0;j<4;++j)
				{
				    u_p[k]+=((Dirac_algebra_type::gamma_L(k,j)*kslash[j][i])/u_p[i]);
				}
			    }
			    break;
			}
		    }

		    /* Compute the negative-helicity basic spinor: */

		    for(size_type i=0;i<4;++i)
		    {
			u_m[i]=0;
			for(size_type j=0;j<4;++j)
			{
			    u_m[i]+=(Dirac_algebra_type::gamma_R(i,j)*kslash[j][i]);
			}
			u_m[i]=sqrt(u_m[i]);
			if(u_m[i] != (value_type)0)
			{
			    for(size_type k=i+1;k<4;++k)
			    {
				for(size_type j=0;j<4;++j)
				{
				    u_m[k]+=((Dirac_algebra_type::gamma_R(k,j)*kslash[j][i])/u_m[i]);
				}
			    }
			    break;
			}
		    }

		    /* Compute the corresponding row spinors: */

		    Dirac_algebra_type::make_bar(u_p_bar.begin(),u_p.begin());
		    Dirac_algebra_type::make_bar(u_m_bar.begin(),u_m.begin());

		    /* Compute the spinors u(±,k_2), needed for the massless
		     * polarization vectors: */

		    Dirac_algebra_type::g_second(value_type((r_value_type)1/std::sqrt((r_value_type)2),0),k_2,u_2_p.begin(),u_m.begin());
		    temp3.reset();
		    Dirac_algebra_type::make_bar(temp3.begin(),u_2_p.begin());
		    Dirac_algebra_type::C_left(u_2_m.begin(),temp3.begin());
		    
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
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_0)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		value_type c(std::sqrt((r_value_type)0.5/spacetime_type::dot(*p,k_0)),0);
		Dirac_algebra_type::g_second(c,*p,it,u_m.begin());

	    }

	    /* Positive-helicity massless spinor construction with prefactor: */

	    static void make_u_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_0)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		value_type c(std::sqrt((r_value_type)0.5/spacetime_type::dot(*p,k_0)),0);
		Dirac_algebra_type::g_second(h*c,*p,it,u_m.begin());

	    }

	    /* Positive-helicity massless row-spinor construction: */

	    static void make_u_p_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_0)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		value_type c(std::sqrt((r_value_type)0.5/spacetime_type::dot(*p,k_0)),0);
		Dirac_algebra_type::g_third(c,*p,u_m_bar.begin(),it);

	    }

	    /* Positive-helicity massless row-spinor construction with prefactor: */

	    static void make_u_p_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_0)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		value_type c(std::sqrt((r_value_type)0.5/spacetime_type::dot(*p,k_0)),0);
		Dirac_algebra_type::g_third(h*c,*p,u_m_bar.begin(),it);

	    }

	    /* Negative-helicity massless spinor construction: */

	    static void make_u_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp3.reset();
		make_u_p_bar(temp3.begin(),p,m);
		Dirac_algebra_type::C_left(it,temp3.begin());
	    }

	    /* Negative-helicity massless spinor construction with prefactor: */

	    static void make_u_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp3.reset();
		make_u_p_bar(h,temp3.begin(),p,m);
		Dirac_algebra_type::C_left(it,temp3.begin());
	    }

	    /* Negative-helicity massless row-spinor construction: */

	    static void make_u_m_bar(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp3.reset();
		make_u_p(temp3.begin(),p,m);
		Dirac_algebra_type::Cc_right(it,temp3.begin());
	    }

	    /* Negative-helicity massless row-spinor construction with prefactor: */

	    static void make_u_m_bar(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		temp3.reset();
		make_u_p(h,temp3.begin(),p,m);
		Dirac_algebra_type::Cc_right(it,temp3.begin());
	    }

	    /* Positive-helicity massless polarisation vector construction: */

	    static void make_e_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_2)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		temp1.reset();
		make_u_p_bar(temp1.begin(),p,m);
		value_type c((r_value_type)0.5/std::sqrt(spacetime_type::dot(*p,k_2)));
		Dirac_algebra_type::g_first(c,it,temp1.begin(),u_2_p.begin());

	    }

	    /* Positive-helicity massless polarisation vector construction with
	     * prefactor: */

	    static void make_e_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_2)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		temp1.reset();
		make_u_p_bar(h,temp1.begin(),p,m);
		value_type c((r_value_type)0.5/std::sqrt(spacetime_type::dot(*p,k_2)));
		Dirac_algebra_type::g_first(c,it,temp1.begin(),u_2_p.begin());

	    }

	    /* Longitudinal massless polarisation vector construction: */

	    static void make_e_L(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_MESSAGE("longitudinal massless vector polarisation construction called");
		return;
	    }

	    /* Longitudinal massless polarisation vector construction with
	     * prefactor: */

	    static void make_e_L(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_MESSAGE("longitudinal massless vector polarisation construction called");
		return;
	    }

	    /* Negative-helicity massless polarisation vector construction: */

	    static void make_e_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_2)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		temp1.reset();
		make_u_m_bar(temp1.begin(),p,m);
		value_type c((r_value_type)0.5/std::sqrt(spacetime_type::dot(*p,k_2)));
		Dirac_algebra_type::g_first(c,it,temp1.begin(),u_2_m.begin());

	    }

	    /* Negative-helicity massless polarisation vector construction with
	     * prefactor: */

	    static void make_e_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((spacetime_type::dot(*p,k_2)==(r_value_type)0),"inappropriate spinor base for given momentum");
		CAMGEN_MESSAGE_IF((m!=NULL),"nonzero mass pointer encountered");

		temp1.reset();
		make_u_m_bar(h,temp1.begin(),p,m);
		value_type c((r_value_type)0.5/std::sqrt(spacetime_type::dot(*p,k_2)));
		Dirac_algebra_type::g_first(c,it,temp1.begin(),u_2_m.begin());

	    }

	    /* Positive-helicity massive polarisation vector construction: */

	    static void make_m_e_p(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((*m==0),"zero mass value encountered");

		temp1.reset();
		temp2.reset();
		momentum_type p1,p2;
		spacetime_type::split_momentum(p1,p2,*p,*m);
		make_u_p(temp2.begin(),&p2,NULL);
		make_u_p_bar(temp1.begin(),&p1,NULL);
		Dirac_algebra_type::g_first((r_value_type)1/(std::sqrt((r_value_type)2)*(*m)),it,temp1.begin(),temp2.begin());

	    }

	    /* Positive-helicity massive polarisation vector construction with
	     * prefactor: */

	    static void make_m_e_p(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((*m==0),"zero mass value encountered");

		temp1.reset();
		temp2.reset();
		momentum_type p1,p2;
		spacetime_type::split_momentum(p1,p2,*p,*m);
		make_u_p(temp2.begin(),&p2,NULL);
		make_u_p_bar(temp1.begin(),&p1,NULL);
		Dirac_algebra_type::g_first(h/(std::sqrt((r_value_type)2)*(*m)),it,temp1.begin(),temp2.begin());

	    }

	    /* Longitudinal polarisation vector construction: */

	    static void make_m_e_L(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((*m==0),"zero mass value encountered");
		
		momentum_type p1,p2;
		spacetime_type::split_momentum(p1,p2,*p,*m);
		it[0]=(p1[0]-p2[0])/(*m);
		it[1]=(p1[1]-p2[1])/(*m);
		it[2]=(p1[2]-p2[2])/(*m);
		it[3]=(p1[3]-p2[3])/(*m);

	    }

	    /* Longitudinal polarisation vector construction with prefactor: */

	    static void make_m_e_L(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((*m==0),"zero mass value encountered");
		
		momentum_type p1,p2;
		spacetime_type::split_momentum(p1,p2,*p,*m);
		it[0]+=h*(p1[0]-p2[0])/(*m);
		it[1]+=h*(p1[1]-p2[1])/(*m);
		it[2]+=h*(p1[2]-p2[2])/(*m);
		it[3]+=h*(p1[3]-p2[3])/(*m);

	    }
	    
	    /* Negative-helicity massive polarisation vector construction: */

	    static void make_m_e_m(iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((*m==0),"zero mass value encountered");

		temp1.reset();
		temp2.reset();
		momentum_type p1,p2;
		spacetime_type::split_momentum(p1,p2,*p,*m);
		make_u_m(temp2.begin(),&p2,NULL);
		make_u_m_bar(temp1.begin(),&p1,NULL);
		Dirac_algebra_type::g_first((r_value_type)1/(std::sqrt((r_value_type)2)*(*m)),it,temp1.begin(),temp2.begin());

	    }

	    /* Negative-helicity massive polarisation vector construction with
	     * prefactor: */

	    static void make_m_e_m(const value_type& h,iterator it,const momentum_type* p,const r_value_type* m)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		CAMGEN_ERROR_IF((m==NULL),"NULL mass pointer encountered");
		CAMGEN_MESSAGE_IF((!initialised),"call without initialised class");
		CAMGEN_MESSAGE_IF((*m==0),"zero mass value encountered");
		
		temp1.reset();
		temp2.reset();
		momentum_type p1,p2;
		spacetime_type::split_momentum(p1,p2,*p,*m);
		make_u_m(temp2.begin(),&p2,NULL);
		make_u_m_bar(temp1.begin(),&p1,NULL);
		Dirac_algebra_type::g_first(h/(std::sqrt((r_value_type)2)*(*m)),it,temp1.begin(),temp2.begin());

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
		os<<"additional spinors:"<<std::endl<<std::endl;;
		for(size_type i=0;i<4;++i)
		{
		    os<<"\t\t\t";
		    os<<std::left<<"[";
		    os<<std::setw(2);
		    shortprint(os,u_2_p[i]);
		    os<<"]"<<"\t\t\t";
		    os<<std::left<<"[";
		    os<<std::setw(2);
		    shortprint(os,u_2_m[i]);
		    os<<"]"<<std::endl;
		}
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
	    static tensor_type temp3;
    };

    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::u_p(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::u_m(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::u_p_bar(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::u_m_bar(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::u_2_p(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::u_2_m(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::temp1(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::temp2(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::tensor_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::temp3(1,4);
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::momentum_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::k_0;
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::momentum_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::k_1;
    template<class Dirac_alg_t,int beam_dir>typename massless_spinor_factory<Dirac_alg_t,beam_dir,4>::momentum_type massless_spinor_factory<Dirac_alg_t,beam_dir,4>::k_2;
    template<class Dirac_alg_t,int beam_dir>bool massless_spinor_factory<Dirac_alg_t,beam_dir,4>::initialised=false;
}

#ifdef CAMGEN_PAULI_BASIS_H_
#include <Camgen/spinfacPb.h>
#endif /*CAMGEN_PAULI_BASIS_H_*/

#ifdef CAMGEN_WEYL_BASIS_H_
#include <Camgen/spinfacWb.h>
#endif /*CAMGEN_WEYL_BASIS_H_*/

#endif /*CAMGEN_SPINOR_FAC*/

