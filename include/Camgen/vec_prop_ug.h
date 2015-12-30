//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VEC_PROP_UG_H_
#define CAMGEN_VEC_PROP_UG_H_

#include <vector>
#include <Camgen/width_scheme.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Unitary-gauge vector propagator class declaration and definition. *
 *                                                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class unitary_gauge: public width_scheme<model_t>
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Specetime type definition: */

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

	    /* Function computing the propagator denominator and gauge term
	     * denominator and storing them in the 'bigdenom' and 'smalldenom'
	     * variables and storing the momentum address: */

	    static void refresh(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		
		width_scheme<model_t>::evaluate(p,m,w);
	    }

	    /* Function propagating the vector 'first' points to: */

	    static void evaluate(iterator first)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;
		
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		CAMGEN_ERROR_IF((first.range()<model_t::dimension),"tensor iterator out of range");
		
		value_type a=spacetime_type::dot(*p,first)/width_scheme<model_t>::gauge_mass2;
		for(int mu=0;mu<model_t::dimension;++mu)
		{
		    (*first)-=a*(*p)[mu];
		    (*first)*=(-width_scheme<model_t>::denominator);
		    ++first;
		}
	    }

	    /* Function propagating the vectors between 'first' an 'last': */

	    static void evaluate_range(iterator first,iterator last)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;

		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		CAMGEN_ERROR_IF((first.range()<model_t::dimension),"tensor iterator out of range");
		CAMGEN_ERROR_IF(((last-first)<0),"bounding iterator smaller than running iterator");
		CAMGEN_ERROR_IF(((last-first)%model_t::dimension != 0),"first and last iterator differ by non-integer multiple of vector size");
		value_type a;
		do
		{
		    a=spacetime_type::dot(*p,first)/width_scheme<model_t>::gauge_mass2;
		    for(size_type mu=0;mu<model_t::dimension;++mu)
		    {
			(*first)-=a*(*p)[mu];
			(*first)*=(-width_scheme<model_t>::denominator);
			++first;
		    }
		}
		while(first!=last);
	    }
    };
    template<class model_t>const typename unitary_gauge<model_t>::size_type unitary_gauge<model_t>::size;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VEC_PROP_UG_H_*/

