//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FUNDAM_REP_H_
#define CAMGEN_FUNDAM_REP_H_

#include <Camgen/rep.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the fundamental representation class template. By default, it  *
 * has an implementation class with an empty fill_generators() function. For a   *
 * nontrivial specialisation, see the su(n)-file.                                *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Fundamental representation wrapper class definition for default groups: */
    
    template<class group_t>class fundamental_rep
    {
	public:

	    /* Reference typedef to the group wrapper class: */

	    typedef group_t group_type;
	    
	    /* Dimension of the representation, default zero: */
	    
	    static const std::size_t dimension=0;

	    /* Defining the multiplicity w.r.t. representation composition: */

	    static const std::size_t multiplicity=1;

	    /* Implementation class template, derived from the representation
	     * class template. The template argument type is the value type of
	     * the generator metrices:*/

	    template<class value_t>class implementation: public representation<value_t,implementation<value_t>,0,group_t::dimension>
	    {
		public:

		    /* Reference typedef to the base class: */

		    typedef representation<value_t,implementation<value_t>,0,group_t::dimension> base_type;

		    /* Reference typedef to the group implementation class: */

		    typedef typename group_t::template implementation<value_t> group_implementation;

		    /* Static constant boolean denoting whether the
		     * representation is unitary: */
		    
		    static const bool hermitian=false;
	    };
    };
    template<class group_t>const std::size_t fundamental_rep<group_t>::dimension;
    template<class group_t>template<class value_t>const bool fundamental_rep<group_t>::implementation<value_t>::hermitian;
}

#endif /*CAMGEN_FUNDAM_REP_H_*/

