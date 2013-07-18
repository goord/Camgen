//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef SUSY_QEDWB_H_
#define SUSY_QEDWB_H_

#include <Camgen/susy_QED_base.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Weyl_basis.h>
#include <Camgen/hel_type.h>

namespace Camgen
{
    class susy_QEDWb: public susy_QED_base<susy_QEDWb,double>
    {
	public:
	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    typedef Weyl_basis Dirac_algebra_type;
	    typedef helicity_type spin_vector_type;
	    
	    static const std::size_t dimension=4;
	    static const bool coloured=false;
	    static const bool continuous_helicities=true;
	    static const int beam_direction=3;

	    susy_QEDWb();
    };
}

#endif /*SUSY_QEDWB_H_*/

