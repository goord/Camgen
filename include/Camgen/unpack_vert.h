//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_UNPACK_VERT_H_
#define CAMGEN_UNPACK_VERT_H_

#include <Camgen/comp_vert.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Helper class template to 'unpack' chained compositions of vertex colour *
 * structures.                                                             *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class Feynrule_t1,class Feynrule_t2,std::size_t N=Feynrule_t1::multiplicity>class unpack_vertices
    {
	public:

	    typedef typename unpack_vertices<typename Feynrule_t1::min_back,compose_vertices<typename Feynrule_t1::back,Feynrule_t2>,N-1>::type type;
    };

    template<class Feynrule_t1,class Feynrule_t2>class unpack_vertices<Feynrule_t1,Feynrule_t2,1>
    {
	public:

	    typedef compose_vertices<Feynrule_t1,Feynrule_t2> type;
    };
}

#endif /*CAMGEN_UNPACK_VERT_H_*/


