//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Tests for ROOT interface to event/process generators  *
 *                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef std::random rn_engine;
    typedef typename model_type::value_type value_type;
    typedef std::size_t size_type;

    license_print::disable();
}
