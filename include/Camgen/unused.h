//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_UNUSED_H_
#define CAMGEN_UNUSED_H_

#include <cstddef>

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Multi-purpose class denoting an empty type. *
 *                                             *
 * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class unused
    {
	public:

	    /* For convenience, we define some constant static properties... */

	    static const std::size_t index_range=0;
	    static const bool decomposes=false;
    };

    /* Dummy zero-dimensional spacetime type definition: */

    template<class value_t>class unused_spacetime
    {
	public:

	    typedef value_t value_type;
	    static const std::size_t dimension=0;
    };
    template<class value_t>const std::size_t unused_spacetime<value_t>::dimension;
}

#endif /*CAMGEN_UNUSED_H_*/

