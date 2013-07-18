//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FASTYM4V_H_
#define CAMGEN_FASTYM4V_H_

#include <Camgen/asymtvv.h>
#include <Camgen/YM3v.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the fast Yang-Mills 4-vertex. This is a dummy Feynman rule     *
 * type, and its assignment to e.g. a four-gluon vertex results in the creation  *
 * of an auxiliary tensor field in the adjoint representation, coupling to the   *
 * gluons via an antisymmetric tensor vertex (asymtvv) with a structure constant *
 * as colour structure.                                                          *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
	template<class model_t>class fast_YM4v;
}

#endif /*CAMGEN_FASTYM4V_H_*/
