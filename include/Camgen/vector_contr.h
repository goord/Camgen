//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VECTOR_CONTR_H_
#define CAMGEN_VECTOR_CONTR_H_

#include <vector>
#include <Camgen/eval.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Vector particle contraction class declaration and definition. *
 *                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Vector contraction class definition: */

    template<class model_t>class vector_contraction
    {
	public:

	    /* Fermionic property: */

	    static const bool fermionic=false;
	    
	    /* Metric tensor size: */
	    
	    static const std::size_t tensor_size=(model_t::dimension)*(model_t::dimension);

	    /* Contracted subamplitude size: */

	    static const std::size_t size=model_t::dimension;
    };
    template<class model_t>const bool vector_contraction<model_t>::fermionic;
    template<class model_t>const std::size_t vector_contraction<model_t>::tensor_size;
    template<class model_t>const std::size_t vector_contraction<model_t>::size;

    /* Evaluate class template specialisation: */

    template<class model_t>class evaluate< vector_contraction<model_t> >
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

	    /* Function returning the vector index ranges: */

	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		r.resize(1);
		r[0]=model_t::dimension;
		return r;
	    }

	    /* Contraction functions, applying the dot product defined in the
	     * spacetime class: */

	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		return spacetime_type::dot(it1,it2);
	    }
	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		result+=spacetime_type::dot(it1,it2);
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VECTOR_CONTR_H_*/

