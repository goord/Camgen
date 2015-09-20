//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file inv_gen.h
    \brief Base class for invariant mass generators by inversion method.
 */

#ifndef CAMGEN_INV_VAL_GEN_H_
#define CAMGEN_INV_VAL_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Base class template for one-dimensional Monte Carlo generators by inversion.*
 * The base class carries the derived type as the first template parameter     *
 * (CRTP). The derived type should provide implementations of the generated    *
 * density, its cumulant and the inverse of the cumulant.                      *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/val_gen.h>

namespace Camgen
{
    /// Base class template for invariant mass samplers of type sub_type.

    template<class value_t,class rng_t,class sub_type>class inversion_value_generator: public value_generator<value_t,rng_t>
    {
	typedef value_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Public constructors: */
	    /*----------------------*/

	    /// Default constructor.

	    inversion_value_generator(){}

	    /// Copy constructor.

	    inversion_value_generator(const inversion_value_generator<value_t,rng_t,sub_type>& other):value_generator<value_t,rng_t>(other),min_cumulant(other.min_cumulant),max_cumulant(other.max_cumulant),min_cumulant_lim(other.min_cumulant_lim),max_cumulant_lim(other.max_cumulant_lim){}

	    /// Virtual destructor.

	    ~inversion_value_generator(){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Re-evaluates all internal parameters depending on mininvmass.

	    bool refresh_lower_bound()
	    {
		min_cumulant=(this->lower_bound()==-std::numeric_limits<value_t>::infinity())?min_cumulant_lim:static_cast<const sub_type*>(this)->cdf(this->lower_bound());
		refresh_norm();
		return this->normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on maxinvmass.

	    bool refresh_upper_bound()
	    {
		max_cumulant=(this->upper_bound()==std::numeric_limits<value_type>::infinity())?max_cumulant_lim:static_cast<const sub_type*>(this)->cdf(this->upper_bound());
		refresh_norm();
		return this->normalisable();
	    }

	    /// Normalisation refresher.

	    virtual void refresh_norm()
	    {
		this->norm=std::abs(max_cumulant-min_cumulant);
	    }

	    /// Weight evaluation method.

	    value_type get_fast_weight() const
	    {
		return this->norm/static_cast<const sub_type*>(this)->pdf(this->value());
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /// Virtual clone method implementation.

	    virtual inversion_value_generator<value_t,rng_t,sub_type>* clone() const
	    {
		return static_cast<inversion_value_generator<value_t,rng_t,sub_type>*>(new sub_type(*static_cast<const sub_type*>(this)));
	    }

	    /// Mapping implementation.

	    value_type map(const value_type& r) const
	    {
		return static_cast<const sub_type*>(this)->inverse_cdf((max_cumulant-min_cumulant)*r+min_cumulant);
	    }

	    /// Inverse mapping implementation.

	    value_type inverse_map(const value_type& x) const
	    {
		return (static_cast<const sub_type*>(this)->cdf(x)-min_cumulant)/(max_cumulant-min_cumulant);
	    }

	protected:

	    /// Minimal and maximal cumulants.

	    value_type min_cumulant,max_cumulant;

	    /// Limiting values of the cumulants.

	    value_type min_cumulant_lim,max_cumulant_lim;
    };
}

#endif /*CAMGEN_INV_VAL_GEN_H_*/

