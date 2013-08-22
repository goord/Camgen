//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/QED.h>
#include <Camgen/ps_gen_base.h>

using namespace Camgen;

int main()
{
    typedef QED model_type;

    ps_generator_base<model_type>* ps=new momentum_collection<model_type>(3);
}
