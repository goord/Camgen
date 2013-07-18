//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_WZM_H_
#define CAMGEN_WZM_H_

#include <Camgen/unused.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/model.h>
#include <Camgen/hel_type.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the superymmetric Wess-Zumino model in Camgen. The model  *
 * consists of two conjugate scalars and a single Majorana fermion.          *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class WZM:public model<WZM>
    {
	public:

	    /* Compile-time data: */

	    typedef double value_type;
	    static const std::size_t dimension=4;
	    static const bool coloured=false;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    static const bool continuous_helicities=true;
	    typedef helicity_type spin_vector_type;
	    static const int beam_direction=3;
	    
	    /* Scalar-fermion (unbroken susy) mass: */
	    
	    static value_type m;

	    /* Coupling constant: */

	    static value_type g;
	    
	    /* Three-scalar coupling: */
	    
	    static std::complex<value_type>phi3;

	    /* Four-scalar coupling: */

	    static std::complex<value_type>phi4;

	    /* Scalar-fermion coupling: */

	    static std::complex<value_type>phipsipsi;

	    /* Constructor: */

	    WZM();

	    /* Mass input: */

	    static void set_mass(const value_type&);
	    
	    /* Coupling input: */
	    
	    static void set_coupling(const value_type&);

	    /* Function computing the couplings from g and m: */

	    static void refresh_couplings();
    };
}

#endif /*CAMGEN_WZM_H_*/


