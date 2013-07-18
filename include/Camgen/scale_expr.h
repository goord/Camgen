//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file scale_expr.h
    \brief Abstract interface base class for factorisation and renormalisation scale expressions.
 */

#ifndef CAMGEN_SCALE_EXPR_H_
#define CAMGEN_SCALE_EXPR_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for factorisation/renormalisation scale expressions.*
 * Subclasses should implement the get_scale() method.                     *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

namespace Camgen
{
    /// Base type providing values for the factorisation and renormalisation
    /// scales.
    
    template<class value_t>class scale_expression
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;

	    /* Destructor: */

	    virtual ~scale_expression(){}

	    /// Factorisation scale. By default calls QCD_scale().

	    virtual value_type F_scale()
	    {
		return QCD_scale();
	    }

	    /// Renormalisation scale. By default calls QCD_scale().

	    virtual value_type R_scale()
	    {
		return QCD_scale();
	    }

	    /// Overall QCD scale. By default return the Z pole mass 91.1876

	    virtual value_type QCD_scale()
	    {
		return (value_type)91.1876;
	    }
    };

    template<class value_t>class scale_wrapper: public scale_expression<value_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;

	    /* Pointer member: */

	    const value_type* const Q;

	    /// Constructor.

	    scale_wrapper(const value_type* Q_):Q(Q_){}

	    /* QCD scale implementation: */

	    value_type QCD_scale()
	    {
		return *Q;
	    }
    };
}

#endif /*CAMGEN_SCALE_EXPR_H_*/

