//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file phi3.h
    \brief phi^3 + phi^4 theory model definition
 */

#ifndef CAMGEN_PHI34_H_
#define CAMGEN_PHI34_H_

#include <Camgen/unused.h>
#include <Camgen/Minkowski.h>
#include <Camgen/model.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Phi-to-the-third+fourth model declaration header. The particle content is a *
 * neutral scalar phi, the vertices a triple-scalar and four-scalar vertex.    *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class phi34: public model<phi34>
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

	    /* Coupling constants: */

	    static std::complex<value_type>mu;
	    static std::complex<value_type>lambda;

	    /* Scalar mass: */

	    static value_type m;

	    /* Constructor: */

	    phi34();
    };
}

#endif /*CAMGEN_PHI34_H_*/

