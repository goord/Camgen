//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SCALAR_CONTR_H_
#define CAMGEN_SCALAR_CONTR_H_

#include <vector>
#include <Camgen/eval.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Scalar particle contraction class declaration and definition. *
 *                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Scalar contraction class: */

    template<class model_t>class scalar_contraction
    {
	public:

	    /* Fermionic property boolean: */

	    static const bool fermionic=false;
	    
	    /* Contraction tensor size: */
	    
	    static const std::size_t tensor_size=1;

	    /* Function returning the contracted tensor sizes: */

	    static const std::size_t size=1;
    };
    template<class model_t>const bool scalar_contraction<model_t>::fermionic;
    template<class model_t>const std::size_t scalar_contraction<model_t>::tensor_size;
    template<class model_t>const std::size_t scalar_contraction<model_t>::size;

    /* Evaluate class template specialisation: */

    template<class model_t>class evaluate< scalar_contraction<model_t> >
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	    /* Trivial initialisation phase: */

	    static void initialise(){}
	    
	    /* Scalar rank vector filling function: */
	    
	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		r.clear();
		return r;
	    }

	    /* Scalar contraction=multiplication functions: */

	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		return (*it1)*(*it2);
	    }
	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		result+=(*it1)*(*it2);
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SCALAR_CONTR_H_*/

