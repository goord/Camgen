//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_STATIC_COMP_H_
#define CAMGEN_STATIC_COMP_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * *
 * Compile-time equal/unequal comparison classes...*
 *                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class A,class B>class static_eq
    {
	public:
	    static const bool value=false;
    };
    template<class A,class B>const bool static_eq<A,B>::value;

    template<class A>class static_eq<A,A>
    {
	public:
	    static const bool value=true;
    };
    template<class A>const bool static_eq<A,A>::value;

    template<class A,class B>class static_neq
    {
	public:
	    static const bool value=true;
    };
    template<class A,class B>const bool static_neq<A,B>::value;

    template<class A>class static_neq<A,A>
    {
	public:
	    static const bool value=false;
    };
    template<class A>const bool static_neq<A,A>::value;
}

#endif /*CAMGEN_STATIC_COMP_H_*/

