//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef QEDPBPSDH_H_
#define QEDPBPSDH_H_

#include <Camgen/QED_base.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/pol_type.h>

namespace Camgen
{
    class QEDPbpsdh: public QED_base<QEDPbpsdh,double>
    {
	public:
	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    typedef polarised_type spin_vector_type;
	    
	    static const std::size_t dimension=4;
	    static const bool coloured=false;
	    static const bool continuous_helicities=false;
	    static const int beam_direction=3;

	    QEDPbpsdh();
    };
}

#endif /*QEDPBPSDH_H_*/

