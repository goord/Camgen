//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SCALAR_PROP_FG_H_
#define CAMGEN_SCALAR_PROP_FG_H_

#include <vector>
#include <Camgen/width_scheme.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Scalar propagator declaration and definition. There is also a specialisation  *
 * for zero-dimensional models, omitting the square of the momentum.             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class scalar_propagator: public width_scheme<model_t>
    {

	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	
	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);

	public:	

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
	    }

	    /* Propagating subamplitude tensor size: */

	    static const size_type size=1;

	    /* Function returning the scalar index ranges: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& v)
	    {
		v.clear();
		return v;
	    }

	    /* Function computing the propagator denominator and storing the result in
	     * the 'denom' variable: */

	    static void refresh(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		width_scheme<model_t>::evaluate(p,m,w);
	    }

	    /* Function propagating the scalar that 'first' points to: */

	    static void evaluate(iterator first)
	    {
		CAMGEN_ERROR_IF((first.range()==0),"tensor iterator out of range");
		(*first)*=(width_scheme<model_t>::denominator);
	    }

	    /* Function propagating all scalars between 'first' and 'last': */

	    static void evaluate(iterator first,iterator last)
	    {
		CAMGEN_ERROR_IF((first.range()==0),"tensor iterator out of range");
		CAMGEN_ERROR_IF(((last-first)<0),"bounding iterator smaller than running iterator");
		do
		{
		    (*first)*=(width_scheme<model_t>::denominator);
		    ++first;
		}
		while(first!=last);
	    }
    };
    template<class model_t>const typename scalar_propagator<model_t>::size_type scalar_propagator<model_t>::size;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SCALAR_PROP_FG_H_*/

