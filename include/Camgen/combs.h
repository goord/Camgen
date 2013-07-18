//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_COMBS_H_
#define CAMGEN_COMBS_H_

/* * * * * * * * * * * * * * * * * * * * * *
 * Some useful combinatorical functions... *
 *                                         *
 * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Integer power: */

    inline unsigned power(unsigned n,unsigned m)
    {
	unsigned k=1;
	for(unsigned i=0;i<m;++i)
	{
	    k*=n;
	}
	return k;
    }

    /* Factorial: */

    inline unsigned factorial(unsigned N)
    {
	unsigned k=1;
	for(unsigned i=1;i<N+1;i++)
	{
	    k*=i;
	}
	return k;
    }

    /* Double factorial: */

    inline unsigned double_factorial(unsigned N)
    {
	unsigned k=1;
	for(unsigned i=((N%2==0)?2:1);i<N+1;i+=2)
	{
	    k*=i;
	}
	return k;
    }

    /* Factorial running from a lower bound N-k+1: */

    inline unsigned lower_factorial(unsigned N,unsigned k)
    {
	if(k<=N)
	{
	    unsigned n=1;
	    for(unsigned i=N-k+1;i<N+1;i++)
	    {
		n*=i;
	    }
	    return n;
	}
	else
	{
	    return 0;
	}
    }

    /* Binomial coefficients: */

    inline unsigned binomial(unsigned N,unsigned k)
    {
	if(k>N)
	{
	    return 0;
	}
	if(k==0)
	{
	    return 1;
	}
	if(k>N/2)
	{
	    return binomial(N,N-k);
	}
	unsigned result=N-k+1;
	for(unsigned i=2;i<k+1;++i)
	{
	    result*=(N-k+i);
	    result/=i;
	}
	return result;
    }

    /* Sum of binomial coefficients from 0 to k: */

    inline unsigned binomial_sum(unsigned N,unsigned k)
    {
	int n=0;
	for(unsigned i=0;i<k+1;++i)
	{
	    n+=binomial(N,i);
	}
	return n;
    }

    /* Stirling numbers: */

    inline unsigned Stirling_nr(unsigned N,unsigned k)
    {
	if(k>N)
	{
	    return 0;
	}
	if(N==0)
	{
	    return 1;
	}
	if(k==0)
	{
	    return 0;
	}
	return Stirling_nr(N-1,k-1)+k*Stirling_nr(N-1,k);
    }
}

#endif /*CAMGEN_COMBS_H_*/

