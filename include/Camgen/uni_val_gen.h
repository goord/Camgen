//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file uni_val_gen.h
    \brief Uniform Monte Carlo scalar value generator.
 */

#ifndef CAMGEN_UNI_VAL_GEN_H_
#define CAMGEN_UNI_VAL_GEN_H_

#include <Camgen/val_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Uniform scalar Monte Carlo value generator class. *
 *                                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Uniform scalar generator class.

    template<class value_t,class rng_t>class uniform_value_generator: public value_generator<value_t,rng_t>
    {
	typedef value_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Default constructor: */

	    uniform_value_generator(){}

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Mapping implementation: */

	    value_type map(const value_type& r) const
	    {
		return ((this->upper_bound()-this->lower_bound())*r+this->lower_bound());
	    }

	    /* Inverse mapping implementation: */

	    value_type inverse_map(const value_type& s) const
	    {
		return ((s-this->lower_bound())/(this->upper_bound()-this->lower_bound()));
	    }

	    /* Clone method implementation: */

	    uniform_value_generator<value_t,rng_t>* clone() const
	    {
		return new uniform_value_generator<value_t,rng_t>(*this);
	    }

	    /* Normalisation refresher: */

	    void refresh_norm()
	    {
		this->norm=(this->upper_bound()-this->lower_bound());
	    }

	    /* Fast weight evaluation method implementation: */
	    
	    value_type get_fast_weight() const
	    {
		return this->norm;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Returns the channel type: */

	    std::string type() const
	    {
		return "uni";
	    }
    };
}

#endif /*CAMGEN_UNI_S_H_*/

