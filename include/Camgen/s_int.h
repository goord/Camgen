//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file s_int.h
\brief Class template integrating a pair of invariant mass pdf's with energy conservation constraint.
*/

#ifndef CAMGEN_S_INT_H_
#define CAMGEN_S_INT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Class template definition with overloaded static members evaluating the   *
* integrals of the product of 2 invariant mass squared densities within the *
* domain                                                                    *
*                         sqrt(s1)+sqrt(s2)<sqrt(s)                         *
*                                                                           *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <complex>
#include <Camgen/dilog.h>

/* Gauss-Legendre quadrature order definition: */

#define GAUSS_INT_ORDER 20

/* Use numeric double Breit-Wigner evaluation: */

#define NUMERIC_BW_INTEGRAL 1

namespace Camgen
{
    template<class value_t,class rng_t>class s_integrator
    {
	public:

	    typedef value_t value_type;
	    typedef std::complex<value_type> c_value_type;
	    typedef Dd_s_generator<value_t,rng_t> Dirac_delta_type;
	    typedef BW_s_generator<value_t,rng_t> Breit_Wigner_type;
	    typedef pl_s_generator<value_t,rng_t> power_law_type;
	    typedef uni_s_generator<value_t,rng_t> uniform_type;

	    /* Phase-space integral of a constrained s-type decay into 2 (massless)
	     * generic inversion-method channels by Gauss quadrature: */

	    template<class gen_t1,class gen_t2>static value_type integrate(const inversion_s_generator<value_t,rng_t,gen_t1>* gen1,const inversion_s_generator<value_t,rng_t,gen_t2>* gen2,const value_type& sqrts)
	    {
		value_type A=(value_type)0.5*std::abs(gen1->max_cumulant-gen1->min_cumulant);
		value_type B=(value_type)0.5*(gen1->max_cumulant+gen1->min_cumulant);
		value_type integral(0);
		value_type x,y;
		for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
		{
		    x=static_cast<const gen_t1*>(gen1)->inverse_cdf(Legendre<value_type,GAUSS_INT_ORDER>::root(n)*A+B);
		    y=sqrts-std::sqrt(x);
		    integral+=Legendre<value_type,GAUSS_INT_ORDER>::weight(n)*(static_cast<const gen_t2*>(gen2)->cdf(y*y));
		}
		return (A*integral-static_cast<const gen_t1*>(gen1)->norm*static_cast<const gen_t2*>(gen2)->min_cumulant);
	    }

	    /* Phase-space integral involving a Breit-Wigner: swap arguments such that the
	     * peak is integrated analytically, and split the outer interval at
	     * the mass threshold: */

	    template<class gen_t>static value_type integrate(const Breit_Wigner_type* gen1,const inversion_s_generator<value_t,rng_t,gen_t>* gen2,const value_type& sqrts)
	    {
		value_type min=std::min(gen2->min_cumulant,gen2->max_cumulant);
		value_type midpt=sqrts-gen1->mass();
		value_type mid=static_cast<const gen_t*>(gen2)->cdf(midpt*midpt);
		value_type max=std::max(gen2->min_cumulant,gen2->max_cumulant);

		if((min<mid) and (mid<max))
		{
		    value_type A=(value_type)0.5*(max-mid);
		    value_type B=(value_type)0.5*(max+mid);
		    value_type C=(value_type)0.5*(mid-min);
		    value_type D=(value_type)0.5*(mid+min);
		    value_type integral(0);
		    value_type x,y,z;
		    for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
		    {
			x=Legendre<value_type,GAUSS_INT_ORDER>::root(n);
			y=std::pow(sqrts-std::sqrt(static_cast<const gen_t*>(gen2)->inverse_cdf(A*x+B)),(int)2);
			z=std::pow(sqrts-std::sqrt(static_cast<const gen_t*>(gen2)->inverse_cdf(C*x+D)),(int)2);
			integral+=Legendre<value_type,GAUSS_INT_ORDER>::weight(n)*(A*gen1->cdf(y)+C*gen1->cdf(z));
		    }
		    return integral-gen2->norm*gen1->min_cumulant;
		}
		else
		{
		    value_type A=(value_type)0.5*std::abs(gen1->max_cumulant-gen1->min_cumulant);
		    value_type B=(value_type)0.5*(gen1->max_cumulant+gen1->min_cumulant);
		    value_type integral(0);
		    value_type x,y;
		    for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
		    {
			x=gen1->inverse_cdf(Legendre<value_type,GAUSS_INT_ORDER>::root(n)*A+B);
			y=sqrts-std::sqrt(x);
			integral+=Legendre<value_type,GAUSS_INT_ORDER>::weight(n)*(static_cast<const gen_t*>(gen2)->cdf(y*y));
		    }
		    return (A*integral-gen1->norm*gen2->min_cumulant);
		}
	    }

	    /* Utility helper functions: */

	    static value_type h(const value_type& sqrts,const c_value_type& nu1,const value_type& sqrtsmin1,const c_value_type& nu2,const value_type& sqrtsmin2)
	    {
		c_value_type z=sqrts+nu1+nu2;
		c_value_type z1=sqrtsmin1+nu1;
		c_value_type z2=sqrtsmin2+nu2;
		return (Li2((z-z1)/z)-Li2(z2/z)+std::log(z-z1)*std::log(z1/z)-std::log(z2)*(std::log((z-z2)/z)-std::log(((z-z2)/z1)))).real();
	    }

	    static value_type g(const value_type& sqrts,const c_value_type& nu1,const value_type& sqrtsmin1,const c_value_type& nu2,const value_type& sqrtsmin2)
	    {
		return ( h(sqrts,nu1,sqrtsmin1,nu2,sqrtsmin2)
			+h(sqrts,nu1,sqrtsmin1,-nu2,sqrtsmin2)
			+h(sqrts,-nu1,sqrtsmin1,nu2,sqrtsmin2)
			+h(sqrts,-nu1,sqrtsmin1,-nu2,sqrtsmin2));
	    }

	    /* Phase-space integral of a constrained s-type decay into 2 Breit-Wigner
	     * propagators: */

	    static value_type integrate(const Breit_Wigner_type* gen1,const Breit_Wigner_type* gen2,const value_type& sqrts)
	    {
#if NUMERIC_BW_INTEGRAL

		value_type min=std::min(gen2->min_cumulant,gen2->max_cumulant);
		value_type midpt=sqrts-gen1->mass();
		value_type mid=gen2->cdf(midpt*midpt);
		value_type max=std::max(gen2->min_cumulant,gen2->max_cumulant);

		if((min<mid) and (mid<max))
		{
		    value_type A=(value_type)0.5*(max-mid);
		    value_type B=(value_type)0.5*(max+mid);
		    value_type C=(value_type)0.5*(mid-min);
		    value_type D=(value_type)0.5*(mid+min);
		    value_type integral(0);
		    value_type x,y,z;
		    for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
		    {
			x=Legendre<value_type,GAUSS_INT_ORDER>::root(n);
			y=std::pow(sqrts-std::sqrt(gen2->inverse_cdf(A*x+B)),(int)2);
			z=std::pow(sqrts-std::sqrt(gen2->inverse_cdf(C*x+D)),(int)2);
			integral+=Legendre<value_type,GAUSS_INT_ORDER>::weight(n)*(A*gen1->cdf(y)+C*gen1->cdf(z));
		    }
		    return integral-gen2->norm*gen1->min_cumulant;
		}
		else
		{
		    value_type A=(value_type)0.5*std::abs(gen1->max_cumulant-gen1->min_cumulant);
		    value_type B=(value_type)0.5*(gen1->max_cumulant+gen1->min_cumulant);
		    value_type integral(0);
		    value_type x,y;
		    for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
		    {
			x=gen1->inverse_cdf(Legendre<value_type,GAUSS_INT_ORDER>::root(n)*A+B);
			y=sqrts-std::sqrt(x);
			integral+=Legendre<value_type,GAUSS_INT_ORDER>::weight(n)*(gen2->cdf(y*y));
		    }
		    return (A*integral-gen1->norm*gen2->min_cumulant);
		}
#else
		return (value_type)0.25*(g(sqrts,gen1->nu,gen1->m_min(),gen2->nu,gen2->m_min())
			-g(sqrts,gen1->nu,gen1->m_min(),gen2->nubar,gen2->m_min())
			-g(sqrts,gen1->nubar,gen1->m_min(),gen2->nu,gen2->m_min())
			+g(sqrts,gen1->nubar,gen1->m_min(),gen2->nubar,gen2->m_min()))/(gen1->mw*gen2->mw);
#endif
	    }

	    /* Phase-space integrals involving a Dirac delta distribution and a inversion
	     * method channel: */

	    template<class gen_t>static value_type integrate(const inversion_s_generator<value_t,rng_t,gen_t>* gen1,const Dirac_delta_type* gen2,const value_type& sqrts)
	    {
		if(sqrts<gen2->mass()){return (value_type)0;}
		value_type seff=pow(sqrts-gen2->mass(),(int)2);
		value_type result=(gen1->s_min()<=seff)?std::abs(static_cast<const gen_t*>(gen1)->cdf(std::min(seff,gen1->s_max()))-gen1->min_cumulant):(value_type)0;
		return result;
	    }

	    /* Phase-space integral of a constrained state of two on-shell particles: */

	    static value_type integrate(const Dirac_delta_type* gen1,const Dirac_delta_type* gen2,const value_type& sqrts)
	    {
		return (gen1->normalisable() and gen2->normalisable() and (sqrts>=(gen1->mass()+gen2->mass())))?(value_type)1:(value_type)0;
	    }

	    /* Phase-space integral of a constrained state of a uniform
	     * invariant mass distribution and an inversion-method generated
	     * invariant mass spectrum: */

	    template<class gen_t>static value_type integrate(const uniform_type* gen1,const inversion_s_generator<value_t,rng_t,gen_t>* gen2,const value_type& sqrts)
	    {
		value_type A=(value_type)0.5*(gen1->s_max()-gen1->s_min());
		value_type B=(value_type)0.5*(gen1->s_max()+gen1->s_min());
		value_type x,y;
		value_type integral(0);
		for(unsigned n=0;n<GAUSS_INT_ORDER;++n)
		{
		    x=Legendre<value_type,GAUSS_INT_ORDER>::root(n)*A+B;
		    y=sqrts-std::sqrt(x);
		    integral+=Legendre<value_type,GAUSS_INT_ORDER>::weight(n)*(static_cast<const gen_t*>(gen2)->cdf(y*y)-static_cast<const gen_t*>(gen2)->min_cumulant);
		}
		return A*integral;
	    }

	    /* Phase-space integral of constrained s-type decays into a pair of
	     * uniformly generated momenta: */

	    static value_type integrate(const uniform_type* gen1,const uniform_type* gen2,const value_type& sqrts)
	    {
		value_type n=gen1->s_max()-gen1->s_min();
		value_type mmax1=sqrts-gen2->m_min();
		return ((sqrts*sqrts-gen2->s_min())*n + (value_type)0.5*(gen1->s_min()+gen1->s_max())*n - (value_type)4*sqrts*(mmax1*gen1->s_max()-gen1->m_min()*gen1->s_min())/(value_type)3);
	    }

	    /* Phase-space integral of constrained s-type decays into a uniformly
	     * generated momentum and an on-shell particle: */

	    static value_type integrate(const uniform_type* gen1,const Dirac_delta_type* gen2,const value_type& sqrts)
	    {
		if(gen2->mass()<gen2->m_min() or sqrts<(gen1->m_min()+gen2->mass()))
		{
		    return 0;
		}
		return (std::min(std::pow(sqrts-gen2->mass(),(int)2),gen1->s_max())-gen1->s_min());
	    }

	    /* Phase-space intergral of constrained s-type decays into a uniformly
	     * generated momentum and a power-law propagator: */

	    static value_type integrate(const power_law_type* gen1,const uniform_type* gen2,const value_type& sqrts)
	    {
		if(gen1->m!=NULL)
		{
		    return integrate(gen2,static_cast<const inversion_s_generator<value_t,rng_t,power_law_type>*>(gen1),sqrts);
		}
		value_type term1=(sqrts*sqrts-gen2->s_min())*gen1->norm;
		value_type term2;
		if(gen1->exp()==(value_type)1.5)
		{
		    term2=-(value_type)2*sqrts*std::log(gen1->s_max()/gen1->s_min());
		}
		else
		{
		    value_type mu=(value_type)1.5-gen1->exp();
		    term2=-(value_type)2*sqrts*(std::pow(gen1->s_max(),mu)-std::pow(gen1->s_min(),mu))/mu;
		}
		value_type term3;
		if(gen1->exp()==(value_type)2)
		{
		    term3=std::log(gen1->s_max()/gen1->s_min());
		}
		else
		{
		    value_type rho=(value_type)2-gen1->exp();
		    term3=(std::pow(gen1->s_max(),rho)-std::pow(gen1->s_min(),rho))/rho;
		}
		return (term1+term2+term3);
	    }

	    /* Phase-space integral of constrained s-type decay into s uniformly
	     * generated momentum and a Breit-Wigner resonance: */

	    static value_type integrate(const Breit_Wigner_type* gen1,const uniform_type* gen2,const value_type& sqrts)
	    {
		value_type mmax1=sqrts-gen2->m_min();
		return (sqrts*sqrts-gen2->s_min()+gen1->mm)*gen1->norm
		    +0.5*std::log((std::pow(gen1->s_max()-gen1->mm,2)+gen1->mw2)/(std::pow(gen1->s_min()-gen1->mm,2)+gen1->mw2))
		    -(sqrts/gen1->xi_plus)*std::log((gen1->s_max()-(gen1->xi_plus)*mmax1+gen1->alpha)*(gen1->s_min()+gen1->xi_plus*gen1->m_min()+gen1->alpha)/((gen1->s_max()+gen1->xi_plus*mmax1+gen1->alpha)*(gen1->s_min()-gen1->xi_plus*gen1->m_min()+gen1->alpha)))
		    -(value_type)2*sqrts*(std::atan((2*mmax1+gen1->xi_plus)/gen1->xi_min)
			    +std::atan((2*mmax1-gen1->xi_plus)/gen1->xi_min)
			    -std::atan((2*gen1->m_min()+gen1->xi_plus)/gen1->xi_min)
			    -std::atan((2*gen1->m_min()-gen1->xi_plus)/gen1->xi_min))/gen1->xi_min;
	    }
    };
}

#undef GAUSS_INT_ORDER

#endif /*CAMGEN_S_INT_H_*/

