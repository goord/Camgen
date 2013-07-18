//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef QEDPBKSDH_H_
#define QEDPBKSDH_H_

#include <Camgen/QED_base.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/KS_type.h>

namespace Camgen
{
    class QEDPbKsdh: public QED_base<QEDPbKsdh,double>
    {
	public:
	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    typedef KS_type spin_vector_type;
	    
	    static const std::size_t dimension=4;
	    static const bool coloured=false;
	    static const bool continuous_helicities=false;
	    static const int beam_direction=3;

	    QEDPbKsdh();
    };
}

#endif /*QEDPBKSDH_H_*/

