//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_LOWER_BINOM_H_
#define CAMGEN_LOWER_BINOM_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Combinatorical helper class, used for bitstring-to-integer conversion in  *
 * Camgen.                                                                  *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* The lower binomial class' contructor takes two integer arguments: the
     * first denotes the number of objects to choose from (as in the usual
     * binomal function), the second is an arbitrary number. The constructor
     * will sum the binomials N over k for k as long as the result does not
     * exceed the second argument. The last binomial in the sum is stored in
     * 'binom', the total sum in 'binom_sum', the number k in 'order' and the
     * difference of the second argument and the binomial sum in 'difference':*/

    class lower_binomial
    {
	public:
	    /* Data: */

	    unsigned binom;
	    unsigned binom_sum;
	    unsigned order;
	    unsigned difference;

	    /* Constructor: */

	    lower_binomial(unsigned,unsigned);
    };
}

#endif /*CAMGEN_LOWER_BINOM_H_*/
