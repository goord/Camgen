//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_COMPOSE_H_
#define CAMGEN_COMPOSE_H_

#include <cstdlib>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility class template for composition of contractions  *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class T1,class T2,class T3>class compose
    {
	public:
	    
	    static const std::size_t multiplicity=3;

	    typedef T1 front;
	    typedef T3 back;
	    typedef compose<T1,T2,void> min_back;
    };
    
    template<class T1,class T2>class compose<T1,T2,void>
    {
	public:
	    
	    static const std::size_t multiplicity=2;

	    typedef T1 front;
	    typedef T2 back;
	    typedef compose<T1> min_back;
    };
    
    template<class T1>class compose<T1,void,void>
    {
	public:
	    
	    static const std::size_t multiplicity=1;

	    typedef T1 front;
	    typedef T1 back;
    };
}

#endif /*CAMGEN_COMPOSE_H_*/

