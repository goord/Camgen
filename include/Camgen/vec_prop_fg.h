//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VEC_PROP_FG_H_
#define CAMGEN_VEC_PROP_FG_H_

#include <vector>
#include <Camgen/width_scheme.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Feynman-gauge vector propagator class declaration and definition. *
 *                                                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class Feynman_gauge: public width_scheme<model_t>
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

	    /* Propagating subamplitude size: */

	    static const size_type size=model_t::dimension;

	    /* Function returning the vector index ranges: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& v)
	    {
		v.clear();
		v.push_back(model_t::dimension);
		return v;
	    }

	    /* Function computing the propagator denominator and storing the result in
	     * the 'denom' variable: */

	    static void refresh(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		
		width_scheme<model_t>::evaluate(p,m,w);
	    }

	    /* Function propagating the vector that 'first' points to: */
	    
	    static void evaluate(iterator first)
	    {
		CAMGEN_ERROR_IF((first.range()<model_t::dimension),"tensor iterator out of range");
		
		for(int mu=0;mu<model_t::dimension;++mu)
		{
		    (*first)*=(-width_scheme<model_t>::denominator);
		    ++first;
		}
	    }

	    /* Function propagating all vectors between 'first' and 'last': */

	    static void evaluate_range(iterator first,iterator last)
	    {
		CAMGEN_ERROR_IF((first.range()<model_t::dimension),"tensor iterator out of range");
		CAMGEN_ERROR_IF(((last-first)<0),"bounding iterator smaller than running iterator");
		CAMGEN_ERROR_IF(((last-first)%model_t::dimension != 0),"first and last iterator differ by non-integer multiple of vector size");
		
		do
		{
		    (*first)*=(-width_scheme<model_t>::denominator);
		    ++first;
		}
		while(first!=last);
	    }
    };
    template<class model_t>const typename Feynman_gauge<model_t>::size_type Feynman_gauge<model_t>::size;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VEC_PROP_FG_H_*/

