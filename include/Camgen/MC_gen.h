//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file MC_gen.h
    \brief Abstract base class template for Monte Carlo generator classes.
 */

#ifndef CAMGEN_MC_GEN_H_
#define CAMGEN_MC_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for the Monte Carlo generators. *
 *                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/num_utils.h>

namespace Camgen
{
    /// Abstract base class for event generators, containing a method to generate an event and a method to evaluate the
    /// event weight.

    template<class value_t>class MC_generator
    {
	public:

	    /* Numerical value type: */

	    typedef value_t value_type;

	    /// Constructor.

	    MC_generator():w(0){}

	    /// Destructor.

	    virtual ~MC_generator(){}

	    /* Public modifiers: */
	    /*-------------------*/
	    
	    /// 'Bare generation' method (purely virtual). Should update the weight
	    /// data member.
	    
	    virtual bool generate()=0;

            /// Unweighted generation method. Returns failure by default.

            virtual bool generate_unweighted()
            {
                return false;
            }

	    /// Event weight calculation method.

	    virtual bool evaluate_weight()=0;

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Returns the Monte Carlo weight reference.

	    value_type& weight()
	    {
		return w;
	    }

	    /// Returns the Monte Carlo const weight reference.

	    const value_type& weight() const
	    {
		return w;
	    }

	    /// Returns whether the weight is valid.

	    bool validate_weight() const
	    {
		return is_finite_number(w);
	    }

	private:

	    /* Weight member: */

	    value_type w;
    };
}

#endif /*CAMGEN_MC_GEN_H_*/

