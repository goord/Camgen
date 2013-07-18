//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file Legendre.h
    \brief Legendre polynomial class templates, equipped with Gauss quadature functionality.
 */

#ifndef CAMGEN_LEGENDRE_H_
#define CAMGEN_LEGENDRE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Legendre polynomial class templates definitions and implementations. The      *
 * classes define functions to evaluate the polynomial and its derivative at any *
 * point, and to compute the roots and the Gauss-quadrature weights, and a       *
 * member template integrating any given function.                               *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdlib>
#include <cmath>

namespace Camgen
{
    /// Rank-N Legendre polynomials with argument type value_t.

    template<class value_t,std::size_t N>class Legendre
    {
	public:

	    /// Value type definition.

	    typedef value_t value_type;

	    /// Polynomial rank definition.

	    static const std::size_t rank=N;

	    /// Symmetry tag definition.

	    static const bool symmetric=(N%2==0);

	    /// Initialisation method. Computes the roots with Newton-Raphson and the
	    /// corresponding Gauss quadrature weights.

	    static void initialise(std::size_t NR_iterations=10)
	    {
		if(!initialised)
		{
		    value_type pi=std::acos(-(value_type)1);
		    for(std::size_t n=0;n<N;++n)
		    {
			roots[n]=std::cos(value_type(4*n+3)*pi/value_type(4*N+2));
			for(std::size_t i=0;i<NR_iterations;++i)
			{
			    roots[n]-=value(roots[n])/derivative(roots[n]);
			}
			value_type p=derivative(roots[n]);
			weights[n]=(value_type)2/(((value_type)1-roots[n]*roots[n])*p*p);
		    }
		    initialised=true;
		}
	    }

	    /// Legendre polynomial evaluation in point x.

	    static value_type value(const value_type& x)
	    {
		return ((2*N-1)*x*Legendre<value_t,N-1>::value(x)-(N-1)*Legendre<value_t,N-2>::value(x))/(value_type)N;
	    }

	    /// Legendre polynomial derivative evaluation in point x.

	    static value_type derivative(const value_type& x)
	    {
		return ((2*N-1)*(x*Legendre<value_t,N-1>::derivative(x)+Legendre<value_t,N-1>::value(x))-(N-1)*Legendre<value_t,N-2>::derivative(x))/(value_type)N;			
	    }

	    /// Returns the n-th root.

	    static value_type root(std::size_t n)
	    {
		initialise();
		return (n<N)?roots[n]:(value_type)0;
	    }

	    /// Returns the n-th Gauss quadrature integration weight.

	    static value_type weight(std::size_t n)
	    {
		initialise();
		return (n<N)?weights[n]:(value_type)0;
	    }

	    /// Integrates the template argument function on the interval [a,b]
	    /// by Gauss quadrature.

	    template<value_type f(const value_type&)>static value_type Gauss_integrate(const value_type& a,const value_type& b)
	    {
		initialise();
		value_type s=(value_type)0.5*(a+b);
		value_type d=(value_type)0.5*(b-a);
		value_type x;
		value_type result=0;
		for(unsigned i=0;i<N;++i)
		{
		    x=roots[i]*d+s;
		    result+=weights[i]*f(x);
		}
		
		return d*result;
	    }

	private:

	    /* Initialisation flag: */

	    static bool initialised;

	    /* Roots of the Legendre polynomial: */

	    static value_type roots[N];
	    
	    /* Weights of Gauss quadrature: */
	    
	    static value_type weights[N];
    };
    template<class value_t,std::size_t N>bool Legendre<value_t,N>::initialised=false;
    template<class value_t,std::size_t N>typename Legendre<value_t,N>::value_type Legendre<value_t,N>::roots[]={0};
    template<class value_t,std::size_t N>typename Legendre<value_t,N>::value_type Legendre<value_t,N>::weights[]={0};

    /* Zeroth-order Legendre polynomial specialisation: */

    template<class value_t>class Legendre<value_t,0>
    {
	public:

	    /* Value type definition. */

	    typedef value_t value_type;
	    
	    /* Polynomial rank definition. */
	    
	    static const std::size_t rank=0;

	    /* Symmetry tag definition. */

	    static const bool symmetric=true;

	    /* Initialisation method. */

	    static void initialise(){}	
	    
	    /* Legendre polynomial evaluation in point x. */
	    
	    static value_type value(const value_type& x)
	    {
		return (value_type)1;
	    }

	    /* Evaluation of the derivative in point x. */

	    static value_type derivative(const value_type& x)
	    {
		return (value_type)0;

	    }

	    /* Returns the n-th root. */

	    static value_type root(std::size_t n)
	    {
		return (value_type)0;
	    }

	    /* Returns the n-th Gauss quadrature weight. */

	    static value_type weight(std::size_t n)
	    {
		return (value_type)0;
	    }

	    /* Integrates the template argument function on the interval [a,b]
	     * by Gauss quadrature. */

	    template<value_type f(const value_type&)>static value_type Gauss_integrate(const value_type& a,const value_type& b)
	    {
		return (value_type)0;
	    }
    };

    /* First-order Legendre polynomial specialisation: */

    template<class value_t>class Legendre<value_t,1>
    {
	public:

	    /* Value type definition. */

	    typedef value_t value_type;
	    
	    /* Polynomial rank definition. */
	    
	    static const std::size_t rank=1;

	    /* Symmetry tag definition. */

	    static const bool symmetric=false;

	    /* Initialisation method. */

	    static void initialise(){}	
	    
	    /* Legendre polynomial evaluation in point x. */
	    
	    static value_type value(const value_type& x)
	    {
		return x;
	    }

	    /* Evaluation of the derivative in point x. */

	    static value_type derivative(const value_type& x)
	    {
		return (value_type)1;

	    }

	    /* Returns the n-th root. */

	    static value_type root(std::size_t n)
	    {
		return (value_type)0;
	    }

	    /* Returns the n-th Gauss quadrature weight. */

	    static value_type weight(std::size_t n)
	    {
		return (n==0)?(value_type)1:(value_type)0;
	    }

	    /* Integrates the template argument function on the interval [a,b]
	     * by Gauss quadrature. */

	    template<value_type f(const value_type&)>static value_type Gauss_integrate(const value_type& a,const value_type& b)
	    {
		return f((value_type)0);
	    }
    };
}

#endif /* CAMGEN_LEGENDRE_H_*/

