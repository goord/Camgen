//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef QEDWBKSDH_H_
#define QEDWBKSDH_H_

#include <Camgen/QED_base.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Weyl_basis.h>
#include <Camgen/KS_type.h>

namespace Camgen
{
    class QEDWbKsdh: public QED_base<QEDWbKsdh,double>
    {
	public:
	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    typedef Weyl_basis Dirac_algebra_type;
	    typedef KS_type spin_vector_type;
	    
	    static const std::size_t dimension=4;
	    static const bool coloured=false;
	    static const bool continuous_helicities=false;
	    static const int beam_direction=3;

	    QEDWbKsdh();
    };
}

#endif /*QEDWBKSDH_H_*/

