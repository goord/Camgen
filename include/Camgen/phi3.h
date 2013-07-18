//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file phi3.h
    \brief phi^3 theory model definition
 */

#ifndef CAMGEN_PHI3_H_
#define CAMGEN_PHI3_H_

#include <Camgen/unused.h>
#include <Camgen/Minkowski.h>
#include <Camgen/model.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Phi-to-the-third model declaration header. The particle content is a neutral  *
 * scalar phi, the vertex a single sss-vertex.					 *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class phi3: public model<phi3>
    {
	public:

	    /* Numerical type definition: */

	    typedef double value_type;
	    
	    /* Spacetime dimensions definition: */
	    
	    static const std::size_t dimension=4;

	    /* Number of distinct colours a particle can have: */

	    static const bool coloured=false;
	    
	    /* Spacetime signature definition: */
	    
	    typedef Minkowski_type spacetime_type;

	    /* Dirac algebra basis definition: */

	    typedef unused Dirac_algebra_type;
	    
	    /* Boolean denoting whether to mix helicities: */
	    
	    static const bool continuous_helicities=false;

	    /* Coupling constant: */

	    static std::complex<value_type>mu;

	    /* Scalar mass: */

	    static value_type m;

	    /* Constructor: */

	    phi3();
    };
}

#endif /*CAMGEN_PHI3_H_*/

