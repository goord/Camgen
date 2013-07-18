//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/stdrand.h>

namespace std
{
    /* Minimal and maximal thrown values: */

    const random::result_type random::min_value;
    const random::result_type random::max_value;
    
    /* Random number generator seed: */
    
    random::result_type random::seed=290881;

    /* Constructor: */

    random::random()
    {
	srand(seed);
    }

    /* Throwing operator: */

    random::result_type random::operator()(void)
    {
	return rand();
    }
}

