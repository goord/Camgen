//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DIRAC_DIM_H_
#define CAMGEN_DIRAC_DIM_H_

#include <cstddef>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Compile-time determination of the size of spinors: via recursive templates    *
 * the number 2^([D/2]) is computed, where [D/2] is the lower bound of the half  *
 * of the spacetime dimension.                                                   *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<std::size_t dim>class Dirac_dim
    {
	public:
	    static const std::size_t value=2*Dirac_dim<dim-2>::value;
    };

    template<std::size_t dim>const std::size_t Dirac_dim<dim>::value;

    template<>class Dirac_dim<1>
    {
	public:
	    static const std::size_t value=1;
    };

    template<>class Dirac_dim<0>
    {
	public:
	    static const std::size_t value=1;
    };
}

#endif /*CAMGEN_DIRAC_DIM_H_*/

