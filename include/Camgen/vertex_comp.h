//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VERTEX_COMP_H_
#define CAMGEN_VERTEX_COMP_H_

#include <Camgen/vertex.h>

/* * * * * * * * * * * * * * * * * * * * * * * * *
 * Vertex pointer comparison function objects... *
 *                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Binary comparison function object: */

    template<class model_t>class vertex_comp
    {
	public:
	    vertex_comp(){}
	    bool operator () (const vertex<model_t>* v1,const vertex<model_t>* v2) const
	    {
		return *v1 < *v2;
	    }
    };
}

#endif /*CAMGEN_VERTEX_COMP_H_*/


