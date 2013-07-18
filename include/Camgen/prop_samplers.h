//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PROP_SAMPLERS_H_
#define CAMGEN_PROP_SAMPLERS_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Classes sampling invariant mass-squared values according to some propagator *
 * density: BW_sampler generates a Breit-Wigner distribution, pl_sampler a     *
 * power-law and d_sampler a delta function (on-shell particles). Furthermore  *
 * the file includes all implementations of the overloaded phase-space         *
 * integration function for s-type phase space branchings into the defined     *
 * sampled propagator distributions.                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <complex>
#include <Camgen/MC_tool_base.h>
#include <Camgen/dilog.h>

/* Gauss-Legendre quadrature order definition: */

#define GAUSS_INT_ORDER 20

namespace Camgen
{
    /* Utility helper function: */

    template<class value_t>std::complex<value_t> BWhelper(const std::complex<value_t>& z1,const std::complex<value_t>& z2)
    {
	return std::log(z1)*std::log(z2)+Li2(z1)+Li2(z2);
    }
    
    /* Phase-space integral of a constrained s-type decay into 2 Breit-Wigner
     * propagators: */

    template<class value_t,class rng_t>value_t integrate(const BW_sampler<value_t,rng_t>& bw1,const BW_sampler<value_t,rng_t>& bw2,const value_t& s)
    {
	value_t w=std::sqrt(s)-bw1.sqrtsmin-bw2.sqrtsmin;
	if(w<(value_t)0)
	{
	    return (value_t)0;
	}
	std::complex<value_t>result;
	result = BWhelper(bw1.nu/(bw1.nu+bw2.nu-w),bw2.nu/(bw1.nu+bw2.nu-w))+BWhelper(bw1.nu/(bw1.nu-bw2.nu+w),bw2.nu/(-bw1.nu+bw2.nu-w))+BWhelper(bw1.nu/(bw1.nu-bw2.nu-w),bw2.nu/(-bw1.nu+bw2.nu+w))+BWhelper(bw1.nu/(bw1.nu+bw2.nu+w),bw2.nu/(bw1.nu+bw2.nu+w)) ;
	result+=(BWhelper(bw1.nubar/(bw1.nubar+bw2.nubar-w),bw2.nubar/(bw1.nubar+bw2.nubar-w))+BWhelper(bw1.nubar/(bw1.nubar-bw2.nubar+w),bw2.nubar/(-bw1.nubar+bw2.nubar-w))+BWhelper(bw1.nubar/(bw1.nubar-bw2.nubar-w),bw2.nubar/(-bw1.nubar+bw2.nubar+w))+BWhelper(bw1.nubar/(bw1.nubar+bw2.nubar+w),bw2.nubar/(bw1.nubar+bw2.nubar+w)));
	result-=(BWhelper(bw1.nu/(bw1.nu+bw2.nubar-w),bw2.nubar/(bw1.nu+bw2.nubar-w))+BWhelper(bw1.nu/(bw1.nu-bw2.nubar+w),bw2.nubar/(-bw1.nu+bw2.nubar-w))+BWhelper(bw1.nu/(bw1.nu-bw2.nubar-w),bw2.nubar/(-bw1.nu+bw2.nubar+w))+BWhelper(bw1.nu/(bw1.nu+bw2.nubar+w),bw2.nubar/(bw1.nu+bw2.nubar+w)));
	result-=(BWhelper(bw1.nubar/(bw1.nubar+bw2.nu-w),bw2.nu/(bw1.nubar+bw2.nu-w))+BWhelper(bw1.nubar/(bw1.nubar-bw2.nu+w),bw2.nu/(-bw1.nubar+bw2.nu-w))+BWhelper(bw1.nubar/(bw1.nubar-bw2.nu-w),bw2.nu/(-bw1.nubar+bw2.nu+w))+BWhelper(bw1.nubar/(bw1.nubar+bw2.nu+w),bw2.nu/(bw1.nubar+bw2.nu+w)));
	return -(value_t)0.25*result.real()/(bw1.mw*bw2.mw);
    }
    
    /* Phase-space integral of a constrained s-type decay into 2 power-law
     * propagators: */
    
    template<class value_t,class rng_t>value_t integrate(const pl_sampler<value_t,rng_t>& pl1,const pl_sampler<value_t,rng_t>& pl2,const value_t& s)
    {
	value_t zmin1=std::sqrt(pl1.smin()/s);
	value_t zmin2=std::sqrt(pl2.smin()/s);
	if(zmin1+zmin2>(value_t)1)
	{
	    return (value_t)0;
	}
	value_t zplus2=(value_t)1-zmin2;
	value_t a=(value_t)2*((value_t)1-*(pl1.nu));
	value_t b=(value_t)2*((value_t)1-*(pl2.nu));
	value_t integral=(value_t)0;
	value_t up =(value_t)0.5*(zplus2+zmin1);
	value_t low=(value_t)0.5*(zplus2-zmin1);
	value_t x;
	for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
	{
	    x=Legendre<value_t,GAUSS_INT_ORDER>::root(n)*low+up;
	    integral+=Legendre<value_t,GAUSS_INT_ORDER>::weight(n)*std::pow(x,a-(value_t)1)*std::pow(((value_t)1-x),b);
	}
	integral*=low;
	value_t prefactor=(value_t)4*std::pow(s,(value_t)0.5*(a+b))/b;
	return prefactor*(integral-std::pow(zmin2,b)*(std::pow(zplus2,a)-std::pow(zmin1,a))/a);
    }

    /* Phase-space integral of a constrained s-type decay into 2 on-shell
     * particles: */

    template<class value_t,class rng_t>value_t integrate(const d_sampler<value_t,rng_t>& d1,const d_sampler<value_t,rng_t>& d2,const value_t& s)
    {
	return (*(d1.m)+*(d2.m)<=std::sqrt(s))?(value_t)1:(value_t)0;
    }

    /* Phase-space integrals of constrained s-type decays into a Breit-Wigner
     * and a power-law propagator: */

    template<class value_t,class rng_t>value_t integrate(const BW_sampler<value_t,rng_t>& bw,const pl_sampler<value_t,rng_t>& pl,const value_t& s)
    {
	value_t sqrts=std::sqrt(s);
	value_t sqrtsminbw=std::sqrt(bw.smin());
	value_t sqrtsminpl=std::sqrt(pl.smin());
	if(sqrts<=sqrtsminbw+sqrtsminpl)
	{
	    return (value_t)0;
	}
	value_t nuplus=(value_t)1-*(pl.nu);
	value_t low=std::pow(sqrts-sqrtsminpl,(value_t)2*nuplus);
	value_t mid=std::pow(sqrts-*(bw.m),(value_t)2*nuplus);
	value_t integral=(value_t)0;
	value_t umin2=((value_t)1-*(pl.nu))*pl.r_min;
	value_t x,y,z,I1,I2;

	for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
	{
	    x=(value_t)0.5*(Legendre<value_t,GAUSS_INT_ORDER>::root(n)+(value_t)1);
	    y=std::pow(sqrts-std::pow(low+(mid-low)*x,(value_t)0.5*pl.inv_exp),2);
	    I1=std::atan((y-bw.mm)/bw.mw)-bw.r_min;
	    z=std::pow(sqrts-std::pow(umin2-(umin2-mid)*x,(value_t)0.5*pl.inv_exp),2);
	    I2=std::atan((z-bw.mm)/bw.mw)-bw.r_min;
	    integral+=Legendre<value_t,GAUSS_INT_ORDER>::weight(n)*((mid-low)*I1+(umin2-mid)*I2);
	}
	return -(value_t)0.5*integral/(bw.mw*nuplus);
    }

    template<class value_t,class rng_t>value_t integrate(const pl_sampler<value_t,rng_t>& pl,const BW_sampler<value_t,rng_t>& bw,const value_t& s)
    {
	return integrate(bw,pl,s);
    }

    /* Phase-space integrals of constrained s-type decays into a Breit-Wigner
     * propagator and an on-shell particle: */
    
    template<class value_t,class rng_t>value_t integrate(const BW_sampler<value_t,rng_t>& bw,const d_sampler<value_t,rng_t>& d,const value_t& s)
    {
	if(s<d.mm){return (value_t)0;}
	value_t seff=(d.m==NULL)?s:s+d.mm-(value_t)2*(*(d.m))*std::sqrt(s);
	return (bw.smin()<seff)?(std::atan((seff-bw.mm)/bw.mw)-bw.r_min)/bw.mw:(value_t)0;
    }

    template<class value_t,class rng_t>value_t integrate(const d_sampler<value_t,rng_t>& d,const BW_sampler<value_t,rng_t>& bw,const value_t& s)
    {
	return integrate(bw,d,s);
    }

    /* Phase-space integrals of constrained s-type decays into a power-law
     * propagator and an on-shell particle: */
    
    template<class value_t,class rng_t>value_t integrate(const d_sampler<value_t,rng_t>& d,const pl_sampler<value_t,rng_t>& pl,const value_t& s)
    {
	if(s<d.mm){return (value_t)0;}
	value_t seff=(d.m==NULL)?s:s+d.mm-(value_t)2*(*(d.m))*std::sqrt(s);
	return (pl.smin()<seff)?(std::pow(seff,(value_t)1-*(pl.nu))/((value_t)1-*(pl.nu))-pl.r_min):(value_t)0;
    }

    template<class value_t,class rng_t>value_t integrate(const pl_sampler<value_t,rng_t>& pl,const d_sampler<value_t,rng_t>& d,const value_t& s)
    {
	return integrate(d,pl,s);
    }

    /* Phase-space integral of constrained s-type decays into a pair of
     * uniformly generated momenta: */

    template<class value_t,class rng_t>value_t integrate(const u_sampler<value_t,rng_t>& u1,const u_sampler<value_t,rng_t>& u2,const value_t& s)
    {
	value_t x=u1.smin_val+u2.smin_val;
	return s*s/(value_t)6-(value_t)0.5*x*x-s*x+(value_t)4*(u1.smin_val*std::sqrt(s*u1.smin_val)+u2.smin_val*std::sqrt(s*u2.smin_val))/(value_t)3;
    }

    /* Phase-space integral of constrained s-type decays into a uniformly
     * generated momentum and an on-shell particle: */

    template<class value_t,class rng_t>value_t integrate(const u_sampler<value_t,rng_t>& u,const d_sampler<value_t,rng_t>& d,const value_t& s)
    {
	value_t seff=(d.m==NULL)?s:s+d.mm-(value_t)2*(*(d.m))*std::sqrt(s);
	return std::max(seff-u.smin_val,(value_t)0);
    }
    
    template<class value_t,class rng_t>value_t integrate(const d_sampler<value_t,rng_t>& d,const u_sampler<value_t,rng_t>& u,const value_t& s)
    {
	return integrate(u,d,s);
    }

    /* Phase-space intergral of constrained s-type decays into a uniformly
     * generated momentum and a power-law propagator: */

    template<class value_t,class rng_t>value_t integrate(const u_sampler<value_t,rng_t>& u,const pl_sampler<value_t,rng_t>& pl,const value_t& s)
    {
	value_t sqrts=std::sqrt(s);
	value_t a=sqrts-std::sqrt(u.smin_val);
	value_t b=std::sqrt(pl.smin());
	value_t fac=(value_t)2-(value_t)2*(*pl.nu);
	value_t A=std::pow(a,fac);
	value_t B=std::pow(b,fac);
	value_t result=(value_t)2*(s-u.smin_val)*(A-B)/fac;
	++fac;
	A*=a;
	B*=b;
	result-=((value_t)4*sqrts*(A-B)/fac);
	++fac;
	A*=a;
	B*=b;
	result+=((value_t)2*(A-B)/fac);
	return result;
    }

    template<class value_t,class rng_t>value_t integrate(const pl_sampler<value_t,rng_t>& pl,const u_sampler<value_t,rng_t>& u,const value_t& s)
    {
	return integrate(u,pl,s);
    }

    /* Phase-space integral of constrained s-type decay into s uniformly
     * generated momentum and a Breit-Wigner resonance: */

    template<class value_t,class rng_t>value_t integrate(const u_sampler<value_t,rng_t>& u,const BW_sampler<value_t,rng_t>& bw,const value_t& s)
    {
	value_t sqrts=std::sqrt(s);
	value_t sqrts1=sqrts-std::sqrt(u.smin_val);
	value_t sqrts2=std::sqrt(bw.smin());
	value_t seff=sqrts1*sqrts1;
	std::complex<value_t>lterm=bw.nu*std::log(((bw.nu+sqrts1)*(bw.nu-sqrts2))/((bw.nu-sqrts1)*(bw.nu+sqrts2)));
	std::complex<value_t>lbarterm=bw.nubar*std::log(((bw.nubar+sqrts1)*(bw.nubar-sqrts2))/((bw.nubar-sqrts1)*(bw.nubar+sqrts2)));

	return ((s-u.smin_val+bw.mm)/bw.mw)*(std::atan((seff-bw.mm)/bw.mw)-bw.r_min)+(value_t)0.5*std::log(((seff-bw.mm)*(seff+bw.mm)+bw.mw*bw.mw)/(std::pow(bw.smin()-bw.mm,2)+bw.mw*bw.mw))+sqrts*(lterm.real()+lbarterm.real())/bw.mw;
    }
}

#undef GAUSS_INT_ORDER

#endif /*CAMGEN_PROP_SAMPLERS_H_*/

