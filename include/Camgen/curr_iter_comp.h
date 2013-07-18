//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_CURR_ITER_COMP_H_
#define CAMGEN_CURR_ITER_COMP_H_

#include <Camgen/current.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Comparison object declaration and definition for current iterators. It  *
 * simply orders the ietartors according to the currents they point to.    *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,std::size_t N_bits>class current_iter_comp
    {
	public:
	    typedef typename std::vector< current< model_t,N_bits,get_colour_treatment<model_t,model_t::coloured>::decomposes> >::iterator iterator;

	    bool operator () (iterator it1,iterator it2) const
	    {
		return (*it1 < *it2);
	    }
    };
}

#endif /*CAMGEN_CURR_ITER_COMP_H_*/

