//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file MC_int_base.h
    \brief Abstract base class template for Monte Carlo integrator classes.
 */

#ifndef CAMGEN_MC_INT_BASE_H_
#define CAMGEN_MC_INT_BASE_H_

#include <Camgen/num_utils.h>

namespace Camgen
{
    /// Abstract base class for integrator classes.

    template<class value_t>class MC_integrator_base
    {
	public:

	    /* Numerical value type: */

	    typedef value_t value_type;

	    /// Constructor.

	    MC_integrator_base():f(1){}

	    /// Destructor.

	    virtual ~MC_integrator_base(){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Updates internal parameters with the current integrand value.

	    virtual void update()=0;

	    /// Adapts internal parameters with the current integrand value.

	    virtual void adapt(){}

	    /// Resets internal parameters to default values.

	    virtual void reset(){}

	    /// Substitutes the integrand with the argument and performs update().

	    void update(const value_type& f_)
	    {
		f=f_;
		update();
	    }

	    /* Public readout functions: */
	    /*---------------------------*/

	    /// Intergrand value reference.

	    value_type& integrand()
	    {
		return f;
	    }

	    /// Integrand value const reference.

	    const value_type& integrand() const
	    {
		return f;
	    }

	    /// Returns whether the integrand is valid.

	    bool validate_integrand() const
	    {
		return is_finite_number(f);
	    }

	private:

	    /* Integrand: */

	    value_type f;
    };
}

#endif /*CAMGEN_MC_INT_BASE_H_*/


