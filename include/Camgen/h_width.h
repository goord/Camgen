//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file h_width.h
  \brief Evaluates the Higgs width with off-shell vector-boson effects.
  */

#ifndef CAMGEN_H_WIDTH_H_
#define CAMGEN_H_WIDTH_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility functions to evaluate the standard model Higgs decay width due to   *
 * off-shell vector bosons. A small Monte Carlo is evaluated up to some aimed  *
 * precision.                                                                  *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/Breit_Wigner.h>
#include <Camgen/ps_vol.h>
#include <Camgen/stdrand.h>

namespace Camgen
{
    /// Computes the decay with of a scalar with mass mh into vector bosons with
    /// masses m1,m2 and widths w1,w2 assuming a standard SVV-coupling
    /// proportional to c. The last argument denotes the aimed precision.
    
    template<class value_type>value_type offshell_VV_width(std::complex<value_type>& ymc,const value_type& mh,const value_type& m1,const value_type& w1,const value_type& m2,const value_type& w2,const value_type& eps=0.01)
    {
	value_type mw1=m1*w1;
	value_type mw12=mw1*mw1;
	value_type m12=m1*m1;
	value_type mw2=m2*w2;
	value_type mw22=mw2*mw2;
	value_type m22=m2*m2;
	value_type mh2=mh*mh;
	BW_s_generator<value_type,std::random>gen1(&m1,&w1),gen2(&m2,&w2);
	gen1.set_m_range(0,mh);
	gen2.set_m_range(0,mh);
	const value_type& s1=gen1.s();
	const value_type& s2=gen2.s();
	value_type wght(0),wsum(0),w2sum(0),err(1);
	std::size_t N=0,n=0;
	do
	{
	    bool q,q1,q2;
	    n=0;
	    do
	    {
		q1=gen1.generate();
		q2=gen2.generate();
		q=(q1 and q2 and (std::sqrt(s1)+std::sqrt(s2)<=mh));
		++n;
	    }
	    while(!q and (n<max_s_pairs()));
	    if(!q)
	    {
		continue;
	    }
	    else
	    {
		++N;
		value_type num=((value_type)2+std::pow(mh2-s1-s2,(int)2)/(4*m12*m22))*std::sqrt(Kallen(mh2,s1,s2));
		value_type denom=(std::pow(s1-m12,(int)2)+mw12)*(std::pow(s2-m22,(int)2)+mw22);
		wght=integrate(&gen1,&gen2,mh)*gen1.weight()*gen2.weight()*num/denom;
		wsum+=wght;
		w2sum+=(wght*wght);
	    }
	    err=(wsum==0 or N<2)?1:(std::sqrt(w2sum-wsum*wsum/N)/wsum);
	}
	while(err>eps);
	value_type pi=std::acos(-(value_type)1);
	return mw1*mw2*std::norm(ymc)*(wsum/N)/std::pow(2*pi*mh,(int)3);
    }

    /// Returns the Higgs to fermions width

    template<class model_type>typename model_type::value_type Gamma_H_ff(const typename model_type::value_type& eps=0.01)
    {
	typedef typename model_type::value_type value_type;
	value_type result(0);
	value_type MH=model_type::M_h0;
	value_type ml[3]={model_type::M_e,model_type::M_mu,model_type::M_tau};
	value_type pi=std::acos(-(value_type)1);
	value_type prefactor=model_type::G_F*MH/(std::sqrt((value_type)2)*(value_type)4*pi);
	for(int i=0;i<3;++i)
	{
	    if(ml[i]>(value_type)0 and ml[i]<=(value_type)0.5*MH)
	    {
		result+=prefactor*ml[i]*ml[i]*std::pow((value_type)1-std::pow(2*ml[i]/MH,(int)2),(value_type)1.5);
	    }
	}
	value_type mq[6]={model_type::M_d,model_type::M_u,model_type::M_s,model_type::M_c,model_type::M_b,model_type::M_t};
	prefactor*=(model_type::N_c);
	for(int i=0;i<6;++i)
	{
	    if(mq[i]>(value_type)0 and mq[i]<=(value_type)0.5*MH)
	    {
		result+=prefactor*mq[i]*mq[i]*std::pow((value_type)1-std::pow(2*mq[i]/MH,(int)2),(value_type)1.5);
	    }
	}
    }

    /// Returns the Higgs to vector bosons width

    template<class model_type>typename model_type::value_type Gamma_H_VV(const typename model_type::value_type& eps=0.01)
    {
	typedef typename model_type::value_type value_type;
	value_type result(0);
	result+=offshell_VV_width(model_type::hWW,model_type::M_h0,model_type::M_W,model_type::W_W,model_type::M_W,model_type::W_W,eps);
	result+=(0.5*offshell_VV_width(model_type::hZZ,model_type::M_h0,model_type::M_Z,model_type::W_Z,model_type::M_Z,model_type::W_Z,eps));
	return result;
    }

    template<class value_type>class H_branchings
    {
	public:

	    value_type M_h0;
	    value_type W_h0;

	    value_type Br_l[3];
	    value_type Br_q[6];
	    value_type Br_W,Br_Z;

	    H_branchings():M_h0(0),W_h0(0),Br_W(0),Br_Z(0)
	    {
		for(int i=0;i<3;++i)
		{
		    Br_l[i]=(value_type)0;
		}
		for(int i=0;i<6;++i)
		{
		    Br_q[i]=(value_type)0;
		}
	    }
    };

    /// Evaluates and substitutes the higgs decay width with fermionic and
    /// (off-shell) vector boson contributions.

    template<class model_type>H_branchings<typename model_type::value_type>Gamma_H(const typename model_type::value_type& eps=0.01)
    {
	typedef typename model_type::value_type value_type;
	H_branchings<value_type>result;
	value_type MH=model_type::M_h0;
	result.M_h0=MH;

	value_type ml[3]={model_type::M_e,model_type::M_mu,model_type::M_tau};
	value_type pi=std::acos(-(value_type)1);
	value_type prefactor=model_type::G_F*MH/(std::sqrt((value_type)2)*(value_type)4*pi);
	for(int i=0;i<3;++i)
	{
	    if(ml[i]>(value_type)0 and ml[i]<=(value_type)0.5*MH)
	    {
		result.Br_l[i]=prefactor*ml[i]*ml[i]*std::pow((value_type)1-std::pow(2*ml[i]/MH,(int)2),(value_type)1.5);
		result.W_h0+=(result.Br_l[i]);
	    }
	}
	value_type mq[6]={model_type::M_d,model_type::M_u,model_type::M_s,model_type::M_c,model_type::M_b,model_type::M_t};
	prefactor*=(model_type::N_c);
	for(int i=0;i<6;++i)
	{
	    if(mq[i]>(value_type)0 and mq[i]<=(value_type)0.5*MH)
	    {
		result.Br_q[i]+=prefactor*mq[i]*mq[i]*std::pow((value_type)1-std::pow(2*mq[i]/MH,(int)2),(value_type)1.5);
		result.W_h0+=(result.Br_q[i]);
	    }
	}
	result.Br_W=offshell_VV_width(model_type::hWW,MH,model_type::M_W,model_type::W_W,model_type::M_W,model_type::W_W,eps);
	result.W_h0+=result.Br_W;
	result.Br_Z=(0.5*offshell_VV_width(model_type::hZZ,MH,model_type::M_Z,model_type::W_Z,model_type::M_Z,model_type::W_Z,eps));
	result.W_h0+=result.Br_Z;
	for(int i=0;i<3;++i)
	{
	    result.Br_l[i]/=result.W_h0;
	}
	for(int i=0;i<6;++i)
	{
	    result.Br_q[i]/=result.W_h0;
	}
	result.Br_W/=result.W_h0;
	result.Br_Z/=result.W_h0;
	return result;
    }

    /// Evaluates and substitutes the higgs decay width with fermionic and
    /// (off-shell) vector boson contributions.

    template<class model_type>typename model_type::value_type evaluate_Higgs_width(const typename model_type::value_type& eps=0.01)
    {
	typedef typename model_type::value_type value_type;
	model_type::W_h0=(value_type)0;
	value_type MH=model_type::M_h0;
	value_type ml[3]={model_type::M_e,model_type::M_mu,model_type::M_tau};
	value_type pi=std::acos(-(value_type)1);
	value_type prefactor=model_type::G_F*MH/(std::sqrt((value_type)2)*(value_type)4*pi);
	for(int i=0;i<3;++i)
	{
	    if(ml[i]>(value_type)0 and ml[i]<=(value_type)0.5*MH)
	    {
		model_type::W_h0+=prefactor*ml[i]*ml[i]*std::pow((value_type)1-std::pow(2*ml[i]/MH,(int)2),(value_type)1.5);
	    }
	}
	value_type mq[6]={model_type::M_d,model_type::M_u,model_type::M_s,model_type::M_c,model_type::M_b,model_type::M_t};
	prefactor*=(model_type::N_c);
	for(int i=0;i<6;++i)
	{
	    if(mq[i]>(value_type)0 and mq[i]<=(value_type)0.5*MH)
	    {
		model_type::W_h0+=prefactor*mq[i]*mq[i]*std::pow((value_type)1-std::pow(2*mq[i]/MH,(int)2),(value_type)1.5);
	    }
	}
	model_type::W_h0+=offshell_VV_width(model_type::hWW,MH,model_type::M_W,model_type::W_W,model_type::M_W,model_type::W_W,eps);
	model_type::W_h0+=(0.5*offshell_VV_width(model_type::hZZ,MH,model_type::M_Z,model_type::W_Z,model_type::M_Z,model_type::W_Z,eps));
	return model_type::W_h0;
    }
}

#endif /*CAMGEN_H_WIDTH_H_*/

