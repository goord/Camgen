//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file rcarry.h
    \brief RCARRY pseudo-random number generator.
 */

#ifndef CAMGEN_RCARRY_H_
#define CAMGEN_RCARRY_H_

namespace Camgen
{
    /// RCARRY random number stream, based on a 24-integer registry.

    class rcarry
    {
	public:

	    /* Numerical output type: */

	    typedef int result_type;

	    /* Minimal and maximal thrown values: */

	    static const result_type min_value=0;
	    static const result_type max_value=(2<<24);

	    /// Register of seeds.

	    static result_type seeds[24];

	    /// Constructor.

	    rcarry();

	    /// Throwing operator.

	    result_type operator ()(void);

	private:

	    /* Register of values: */

	    result_type reg[24];
	    
	    /* Carry-bit: */
	    
	    bool carry;
    };
}

#endif /*CAMGEN_RCARRY_H_*/

