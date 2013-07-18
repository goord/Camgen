//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/phi3.h>
#include <Camgen/scalar_particle.h>
#include <Camgen/sss.h>

namespace Camgen
{
    /* Static constant integer member definitions: */

    const std::size_t phi3::dimension;
    const bool phi3::coloured;
    const bool phi3::continuous_helicities;

    /* Initialisation of static floating-point members: */

    std::complex<phi3::value_type> phi3::mu=(phi3::value_type)0.5;
    phi3::value_type phi3::m=(phi3::value_type)20;

    /* Constructor: */

    phi3::phi3()
    {
	add_scalar("phi",&m);
	add_vertex<sss>("phi","phi","phi",&mu);
    }
}


