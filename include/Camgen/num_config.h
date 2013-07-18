//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file num_config.h
    \brief Floating-point numeric configuration definitions.
 */

#ifndef CAMGEN_NUM_CONFIG_H_
#define CAMGEN_NUM_CONFIG_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The numeric configuration class template defines various constants that     *
 * determine the behaviour of floating point arithmetic and equality testing in*
 * Camgen.                                                                     *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Forward class template declaration: */

    template<class T>class numeric_configuration;

    /* Specialisation for doubles: */

    template<>class numeric_configuration<double>
    {
	public:

	    const static double epsilon_abs;

	    const static double epsilon_rel;
    };

    /// Numeric configuration class template.

    template<class value_type>class numeric_configuration
    {
	public:

	    /// Absolute maximal difference in equality testing.

	    const static value_type epsilon_abs;

	    /// Relative maximal difference in equality testing.

	    const static value_type epsilon_rel;
    };
    template<class value_type>const value_type numeric_configuration<value_type>::epsilon_abs=static_cast<value_type>(numeric_configuration<double>::epsilon_abs);
    template<class value_type>const value_type numeric_configuration<value_type>::epsilon_rel=static_cast<value_type>(numeric_configuration<double>::epsilon_rel);
}

#endif /*CAMGEN_NUM_CONFIG_H_*/

