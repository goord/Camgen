//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file stdrand.h
    \brief std random number generator wrapper class.
 */

#ifndef CAMGEN_STDRAND_H_
#define CAMGEN_STDRAND_H_

#include <cstdlib>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Wrapper class for the standard library random number generator. *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace std
{
    /// Wrapper class for the standard library random number stream.

    class random
    {
	public:

	    /* Numerical output type: */

	    typedef int result_type;

	    /* Minimal and maximal thrown values: */

	    static const result_type min_value=0;
	    static const result_type max_value=RAND_MAX;

	    /// Random number seed.

	    static int seed;

	    /// Constructor.

	    random();

	    /// Throwing operator.

	    result_type operator ()(void);
    };
}

#endif /*CAMGEN_STDRAND_H_*/

