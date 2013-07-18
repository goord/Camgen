//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file dilog.h
    \brief Complex dilogarithm function definition.
 */

#ifndef CAMGEN_DILOG_H_
#define CAMGEN_DILOG_H_

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Dilogarithm function template definitions.  *
 *                                             *
 * * * * * * * * * * * * * * * * * * * * * * * */                                            

#include <complex>
#include <Camgen/Legendre.h>

/* Order of the dilogarithm series evaluation: */

#define DILOG_SERIES_ORDER 50

/* Order of the Gauss quadrature integration: */

#define INTEGRATION_ORDER 15

/* Numerical constant \pi^2/6: */

#define PI_SQ_OVER_SIX 1.644934066848226

namespace Camgen
{

    /* Function template evaluating the series x^k/k^2 up to the DILOG_SERIES_ORDER
     * integer: */

    template<class value_t>value_t dilog_series(const value_t& x)
    {
	value_t temp=x;
	value_t result=temp;
	for(unsigned n=2;n<DILOG_SERIES_ORDER;++n)
	{
	    temp*=x;
	    result+=temp/(value_t)(n*n);
	}
	return result;
    }

    /* Helper function for real-valued dilogarithm evaluation; if 1/2 < x < 1, the
     * dilogarithm series converges slowly, and the fomula
     *
     * Li2(1-z)=-Li2(z)-log(z)log(1-z)+\pi^2/6
     *
     * is applied before series expansion: */

    template<class value_t>value_t dilog_helper(const value_t& x)
    {
	return (x<0.5)?dilog_series(x):(-dilog_series((value_t)1-x)-std::log(x)*std::log((value_t)1-x)+(value_t)PI_SQ_OVER_SIX);
    }

    /* Helper function for complex-valued dilogarithm evaluation; if |z| < 0.5, the
     * series expansion converges rapidly enough, if 0.5 < |z| < 1, Gauss quadrature
     * is used to evaluate the integral: */

    template<class value_t>std::complex<value_t> dilog_helper(const std::complex<value_t>& z)
    {
	if(std::abs(z)<(value_t)0.5)
	{
	    return dilog_series(z);
	}
	else
	{
	    std::complex<value_t>result(0,0);
	    std::complex<value_t>w(1,0);
	    value_t t;
	    for(unsigned i=0;i<INTEGRATION_ORDER;++i)
	    {
		t=(value_t)0.5*(Legendre<value_t,INTEGRATION_ORDER>::root(i)+(value_t)1);
		result-=(value_t)0.5*Legendre<value_t,INTEGRATION_ORDER>::weight(i)*std::log(w-t*z)/t;
	    }
	    return result;
	}
    }

    /* Dilogarithm function template for real arguments. If -1 < x < 1, the helper
     * function is directly applied, if x < -1, the formula
     *
     * Li2(1/z)=-Li2(z)-log(-z)^2/2-\pi^2/6
     *
     * is applied. */

    /// Dilogarithm definition.
    
    template<class value_t>value_t Li2(const value_t& x)
    {
	if(x==(value_t)0)
	{
	    return x;
	}
	if(x<-(value_t)1)
	{
	    value_t l=std::log(-x);
	    return -dilog_helper((value_t)1/x)-(value_t)0.5*l*l-(value_t)PI_SQ_OVER_SIX;
	}
	else if(x<-(value_t)0.5)
	{
	    return (value_t)0.5*dilog_helper(x*x)-dilog_helper(-x);
	}
	else if(x<1)
	{
	    return dilog_helper(x);
	}
	return HUGE_VAL;
    }

    /* Dilogarithm function template for complex arguments. If |z| < 1, the helper
     * function is directly applied else the formula
     *
     * Li2(1/z)=-Li2(z)-log(-z)^2/2-\pi^2/6
     *
     * is applied. */

    
    /// Dilogarithm function overload from complex arguments.
    
    template<class value_t>std::complex<value_t> Li2(const std::complex<value_t>& z)
    {
	if(std::abs(z)<(value_t)1)
	{
	    return dilog_helper(z);
	}
	else if(!(z.imag()==(value_t)0 and z.real()>(value_t)1))
	{
	    std::complex<value_t>l=std::log(-z);
	    return -dilog_helper(std::complex<value_t>(1,0)/z)-(value_t)0.5*l*l-std::complex<value_t>(PI_SQ_OVER_SIX,0);
	}
	return HUGE_VAL;
    }
}

#undef DILOG_SERIES_ORDER
#undef PI_SQ_OVER_SIX
#undef INTEGRATION_ORDER

#endif /*CAMGEN_DILOG_H_*/

