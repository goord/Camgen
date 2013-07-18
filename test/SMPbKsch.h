//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef SMPBKSCH_H_
#define SMPBKSCH_H_

#include <Camgen/SM_base.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/KS_type.h>
#include <Camgen/col_flow.h>

namespace Camgen
{
    class SMPbKsch: public SM_base<SMPbKsch,double>
    {
	public:

	    /* Compile-time data: */

	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    typedef KS_type spin_vector_type;
	    typedef colour_flow colour_treatment;

	    static const std::size_t dimension=4;
	    static const bool coloured=true;
	    static const std::size_t N_c=SM_params::N_c;
	    static const bool continuous_helicities=true;
	    static const bool continuous_colours=false;
	    static const int beam_direction=3;

	    SMPbKsch();
    };
}

#endif /*SMPBKSCH_H_*/

