//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_NAME_COMP_H_
#define CAMGEN_NAME_COMP_H_

#include <string>
#include <Camgen/forward_decs.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Comparison objects where the criterium is the value returned by a member  *
 * function called 'name()'. This ordering can be used to sort particle and  *
 * vertex addresses by name and perform stl algorithms in containers of such *
 * pointers.                                                                 *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Name ordering class: */

    template<class T>class name_comp
    {
	public:
	    name_comp(){}
	    bool operator () (const T* phi1,const T* phi2) const
	    {
		return (phi1->get_name() < phi2->get_name());
	    }
    };

    /* Name equality class: */

    template<class T>class name_eq
    {
	public:
	    name_eq(){}
	    bool operator () (const T* phi1,const T* phi2) const
	    {
		return (phi1->get_name()==phi2->get_name());
	    }
    };

    /* Name unequality class: */

    template<class T>class name_neq
    {
	public:
	    name_neq(){}
	    bool operator () (const T* phi1,const T* phi2) const
	    {
		return (phi1->get_name() != phi2->get_name());
	    }
    };

    /* Name selection class: */

    template<class T>class name_select
    {
	public:
	    const std::string& name;
	    name_select(const std::string& x):name(x){}
	    name_select(const T* f):name(f->get_name()){}
	    bool operator () (const T* phi) const
	    {
		return (phi->get_name()==name);
	    }
    };
}

#endif /*CAMGEN_NAME_COMP_H_*/

