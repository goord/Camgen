//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/phi34.h>
#include <Camgen/scalar_particle.h>
#include <Camgen/sss.h>
#include <Camgen/ssss.h>

namespace Camgen
{
    /* Static constant integer member definitions: */

    const std::size_t phi34::dimension;
    const bool phi34::coloured;
    const bool phi34::continuous_helicities;

    /* Initialisation of static floating-point members: */

    std::complex<phi34::value_type> phi34::mu=(phi34::value_type)0.5;
    std::complex<phi34::value_type> phi34::lambda=(phi34::value_type)0.25;
    phi34::value_type phi34::m=(phi34::value_type)20;

    /* Constructor: */

    phi34::phi34()
    {
	add_scalar("phi",&m);
	add_vertex<sss>("phi","phi","phi",&mu);
	add_vertex<ssss>("phi","phi","phi","phi",&mu);
    }
}


