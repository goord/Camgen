//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/lower_binom.h>
#include <Camgen/combs.h>

namespace Camgen
{
    /* Constructor: */

    lower_binomial::lower_binomial(unsigned N,unsigned n)
    {
	/* Subtract the total binomial sum until n is smaller than it: */

	n%=(1<<N);
	
	/* Start summing binomials until their sum exceeds n: */

	unsigned m=0;
	unsigned k=0;
	do
	{
	    k+=binomial(N,m);
	    ++m;
	}
	while(k<=n and m != N);

	/* Subtract the last term again: */

	k-=binomial(N,m-1);
	
	/* Store the obtained integers in their respective data members: */
	
	order=m-2;
	binom=binomial(N,order);
	binom_sum=k;
	difference=n-binom_sum;
    }
}


