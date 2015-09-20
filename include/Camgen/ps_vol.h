//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_vol.h
    \brief Massless phase space volume calculating classes.
 */

#ifndef CAMGEN_PS_VOL_H_
#define CAMGEN_PS_VOL_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Routines for computing the volume of phase space (with no cuts imposed).  *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <Camgen/combs.h>

namespace Camgen
{
    /// Kallen phase space weight factor definition.

    template<class value_type>value_type Kallen(value_type x,value_type y,value_type z)
    {
	return x*x+y*y+z*z-(value_type)2*(x*(y+z)+y*z);
    }

    /// Static class computing the massless phase space volume with numerical
    /// type value_t, final-state multiplicity N and in dimension D.
    
    template<class value_t,std::size_t N,std::size_t D=4>class massless_ps
    {
	public:

	    /* TYpe definitions. */

	    typedef value_t value_type;
	    static const std::size_t multiplicity=N;
	    static const std::size_t dimension=D;

	    /// Evaluates the massless phase space volume at CM energy Ecm.
	    
	    static value_type volume(const value_type& Ecm)
	    {
		return prefactor*intpow(Ecm,(D-2)*N-D);
	    }

	    /// Evaluates the massless phase space volume at CM energy Ecm for
	    /// final-state multiplicity n.
	    
	    static value_type volume(const value_type& Ecm,std::size_t n)
	    {
		return (value_type)2*intpow(sqrtpi,(n-1)*(D-2))*intpow(Ecm,(D-2)*n-D)/(gamma2(n*(D-2))*gamma2((n-1)*(D-2))*intpow(value_type(D-2),n));
	    }

	private:

	    /* Digamma function: */

	    static value_type gamma2(std::size_t n)
	    {
		if(n%2==0)
		{
		    return factorial(n/2-1);
		}
		else
		{
		    std::size_t k=n/2;
		    return double_factorial(2*k-1)*sqrtpi/(1<<k);
		}
	    }

	    /* Integer power function: */

	    static value_type intpow(const value_type& x,std::size_t n)
	    {
		value_type result=(value_type)1;
		for(std::size_t i=0;i<n;++i)
		{
		    result*=x;
		}
		return result;
	    }
	    static const value_type pi;
	    static const value_type pi2;
	    static const value_type sqrtpi;
	    static const value_type gammas;
	    static const value_type prefactor;
    };
    template<class value_t,std::size_t N,std::size_t D>const std::size_t massless_ps<value_t,N,D>::dimension;
    template<class value_t,std::size_t N,std::size_t D>const std::size_t massless_ps<value_t,N,D>::multiplicity;
    template<class value_t,std::size_t N,std::size_t D>const typename massless_ps<value_t,N,D>::value_type massless_ps<value_t,N,D>::pi=std::acos(-(value_t)1);
    template<class value_t,std::size_t N,std::size_t D>const typename massless_ps<value_t,N,D>::value_type massless_ps<value_t,N,D>::pi2=massless_ps<value_t,N,D>::pi*massless_ps<value_t,N,D>::pi;
    template<class value_t,std::size_t N,std::size_t D>const typename massless_ps<value_t,N,D>::value_type massless_ps<value_t,N,D>::sqrtpi=std::sqrt(massless_ps<value_t,N,D>::pi);
    template<class value_t,std::size_t N,std::size_t D>const typename massless_ps<value_t,N,D>::value_type massless_ps<value_t,N,D>::gammas=massless_ps<value_t,N,D>::intpow(value_t(D-2),N)*massless_ps<value_t,N,D>::gamma2(N*(D-2))*massless_ps<value_t,N,D>::gamma2((N-1)*(D-2));
    template<class value_t,std::size_t N,std::size_t D>const typename massless_ps<value_t,N,D>::value_type massless_ps<value_t,N,D>::prefactor=(value_t)2*massless_ps<value_t,N,D>::intpow(massless_ps<value_t,N,D>::sqrtpi,(D-2)*(N-1))/massless_ps<value_t,N,D>::gammas;
    

    /* Phase space volume specialisation for 2-dimensional theories: */

    template<class value_t,std::size_t N>class massless_ps<value_t,N,2>
    {
	public:
	    
	    typedef value_t value_type;
	    static const std::size_t multiplicity=N;
	    static const std::size_t dimension=2;

	    static value_type volume(const value_type& Ecm)
	    {
		return (value_type)0;
	    }
	    static value_type volume(const value_type& Ecm,std::size_t n)
	    {
		return (value_type)0;
	    }
    };

    /* Phase space volume specialisation for 1-dimensional theories: */

    template<class value_t,std::size_t N>class massless_ps<value_t,N,1>
    {
	public:
	    
	    typedef value_t value_type;
	    static const std::size_t multiplicity=N;
	    static const std::size_t dimension=1;

	    static value_type volume(const value_type& Ecm)
	    {
		return (value_type)0;
	    }
	    static value_type volume(const value_type& Ecm,std::size_t n)
	    {
		return (value_type)0;
	    }
    };

    /* Phase space volume specialisation for 0-dimensional theories: */

    template<class value_t,std::size_t N>class massless_ps<value_t,N,0>
    {
	public:
	    
	    typedef value_t value_type;
	    static const std::size_t multiplicity=N;
	    static const std::size_t dimension=0;

	    static value_type volume(const value_type& Ecm)
	    {
		return (value_type)0;
	    }
	    static value_type volume(const value_type& Ecm,std::size_t n)
	    {
		return (value_type)0;
	    }
    };
}

#endif /*CAMGEN_PS_VOL_H_*/

