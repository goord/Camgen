//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file inv_gen.h
    \brief Base class for invariant mass generators by inversion method.
 */

#ifndef CAMGEN_INV_GEN_H_
#define CAMGEN_INV_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Base class template for invariant mass generators by the inversion method.  *
 * The base class carries the derived type as the first template parameter     *
 * (CRTP). The derived type should provide implementations of the generated    *
 * density, its cumulant and the inverse of the cumulant.                      *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/s_gen.h>
#include <Camgen/s_int.h>

namespace Camgen
{
    /// Base class template for invariant mass samplers of type sub_type.

    template<class value_t,class rng_t,class sub_type>class inversion_s_generator: public s_generator<value_t,rng_t>
    {
	friend class s_integrator<value_t,rng_t>;
	typedef s_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Public constructors: */
	    /*----------------------*/

	    /// Default constructor.

	    inversion_s_generator(){}

	    /// Constructor with invariant mass address.

	    inversion_s_generator(value_type* inv_mass_):base_type(inv_mass_){}

	    /// Copy constructor.

	    inversion_s_generator(const inversion_s_generator<value_t,rng_t,sub_type>& other):s_generator<value_t,rng_t>(other),min_cumulant(other.min_cumulant),max_cumulant(other.max_cumulant),min_cumulant_lim(other.min_cumulant_lim),max_cumulant_lim(other.max_cumulant_lim){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Re-evaluates all internal parameters depending on mininvmass.

	    bool refresh_s_min()
	    {
		min_cumulant=(this->s_min()==-std::numeric_limits<value_t>::infinity())?min_cumulant_lim:static_cast<const sub_type*>(this)->cdf(this->s_min());
		refresh_norm();
		return this->normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on maxinvmass.

	    bool refresh_s_max()
	    {
		max_cumulant=(this->s_max()==std::numeric_limits<value_type>::infinity())?max_cumulant_lim:static_cast<const sub_type*>(this)->cdf(this->s_max());
		refresh_norm();
		return this->normalisable();
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /// Virtual clone method implementation.

	    virtual inversion_s_generator<value_t,rng_t,sub_type>* clone() const
	    {
		return static_cast<inversion_s_generator<value_t,rng_t,sub_type>*>(new sub_type(*static_cast<const sub_type*>(this)));
	    }

	    /// Mapping implementation.

	    value_type map(const value_type& r) const
	    {
		return static_cast<const sub_type*>(this)->inverse_cdf((max_cumulant-min_cumulant)*r+min_cumulant);
	    }

	    /// Inverse mapping implementation.

	    value_type inverse_map(const value_type& s) const
	    {
		return (static_cast<const sub_type*>(this)->cdf(s)-min_cumulant)/(max_cumulant-min_cumulant);
	    }

	    /// Double-dispatch integrator functions.
	    
	    virtual value_type integrate_with(const s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return gen->integrate_with(this,sqrts);
	    }
	    virtual value_type integrate_with(const BW_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::template integrate<sub_type>(gen,this,sqrts);
	    }
	    virtual value_type integrate_with(const pl_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		typedef pl_s_generator<value_t,rng_t> gen_t2;
		return s_integrator<value_t,rng_t>::template integrate<sub_type,gen_t2>(this,static_cast<const inversion_s_generator<value_t,rng_t,gen_t2>*>(gen),sqrts);
	    }
	    virtual value_type integrate_with(const uni_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::template integrate<sub_type>(gen,this,sqrts);
	    }
	    virtual value_type integrate_with(const Dd_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::template integrate<sub_type>(this,gen,sqrts);
	    }

	protected:

	    /// Minimal and maximal cumulants.

	    value_type min_cumulant,max_cumulant;

	    /// Limiting values of the cumulants.

	    value_type min_cumulant_lim,max_cumulant_lim;

	    /// Normalisation refresher.

	    virtual void refresh_norm()
	    {
		this->norm=std::abs(max_cumulant-min_cumulant);
	    }

	    /// Weight evaluation method.

	    value_type get_fast_weight() const
	    {
		return this->norm/static_cast<const sub_type*>(this)->pdf(this->s());
	    }
    };
}

#endif /*CAMGEN_INV_GEN_H_*/

