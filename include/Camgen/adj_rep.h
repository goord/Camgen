//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_ADJ_REP_H_
#define CAMGEN_ADJ_REP_H_

#include <Camgen/rep.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Adjoint representation class definition. Consists of a wrapper, containing    *
 * the symmetry group type and the dimension of the representation, i.e. the     *
 * group rank. Wrapped up is the implementation, derived from the representation *
 * base class. Fills generators with group structure constants.                  *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class group_t>class adjoint_rep
    {
	public:

	    /* Defining group type, as usual given by the template parameter, */
	    
	    typedef group_t group_type;

	    /* Defining dimension, necessary to determine the corresponding
	     * index range, */
	    
	    static const std::size_t dimension=group_type::rank;

	    /* Defining the multiplicity w.r.t. representation composition: */

	    static const std::size_t multiplicity=1;

	    /* Defining the implementation of the generator as a metafunction,
	     * i.e. nested class template*/
	    
	    template<class value_t>class implementation: public representation<value_t,implementation<value_t>,dimension,dimension>
	    {
		public:

		    /* Defining the base representation type */

		    typedef representation<value_t,implementation<value_t>,dimension,dimension> base_type;

		    /* Defining the group implementation metafunction as the
		     * group type, */

		    typedef typename group_type::template implementation<value_t> group_implementation;

		    /* For optimisation purposes, hermiticity of the adjoint
		     * representation is passed as type trait information, */

		    static const bool hermitian=true;

		    /* The generator arrays in the representation base class are
		     * filled with the following static function, according to
		     * the formula T^a_{bc} = i*f_{abc} */

		    static void fill_generators()
		    {
			group_implementation::initialise();
			std::complex<value_t>I(0,1);
			for(std::size_t a=0;a<dimension;++a)
			{
			    for(std::size_t b=0;b<dimension;++b)
			    {
				for(std::size_t c=0;c<dimension;++c)
				{
				    base_type::T[a][b][c]=I*(group_implementation::structure_constant(a,b,c));
				}
			    }
			}
		    }
	    };
    };
    template<class group_t>const std::size_t adjoint_rep<group_t>::dimension;
    template<class group_t>template<class value_t>const bool adjoint_rep<group_t>::implementation<value_t>::hermitian;
}

#endif /*CAMGEN_ADJ_REP_H_*/

