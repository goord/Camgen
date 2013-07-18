//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file uniform.h
    \brief Momentum channel generator with uniform invariant mass generation.
 */

#ifndef CAMGEN_UNI_S_H_
#define CAMGEN_UNI_S_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and implementation of the uniform type invariant mass sampling  *
 * momentum channel.                                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Uniform invariant mass generating momentum channel class.

    template<class value_t,class rng_t>class uni_s_generator: public s_generator<value_t,rng_t>
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

	    /* Default constructor: */

	    uni_s_generator(){}

	    /* Constructor with momentum address argument: */

	    uni_s_generator(value_type* s_):base_type(s_){}

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Mapping implementation: */

	    value_type map(const value_type& r) const
	    {
		return ((this->s_max()-this->s_min())*r+this->s_min());
	    }

	    /* Inverse mapping implementation: */

	    value_type inverse_map(const value_type& s) const
	    {
		return ((s-this->s_min())/(this->s_max()-this->s_min()));
	    }

	    /// Double-dispatch integrator functions.
	    
	    value_type integrate_with(const s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return gen->integrate_with(this,sqrts);
	    }
	    value_type integrate_with(const BW_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(gen,this,sqrts);
	    }
	    value_type integrate_with(const pl_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(gen,this,sqrts);
	    }
	    value_type integrate_with(const uni_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(this,gen,sqrts);
	    }
	    value_type integrate_with(const Dd_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(this,gen,sqrts);
	    }

	    /* Clone method implementation: */

	    uni_s_generator<value_t,rng_t>* clone() const
	    {
		return new uni_s_generator<value_t,rng_t>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Returns the channel type: */

	    std::string type() const
	    {
		return "uni";
	    }

	protected:

	    /* Normalisation refresher: */

	    void refresh_norm()
	    {
		this->norm=(this->s_max()-this->s_min());
	    }

	    /* Fast weight evaluation method implementation: */
	    
	    value_type get_fast_weight() const
	    {
		return this->norm;
	    }
    };
}

#endif /*CAMGEN_UNI_S_H_*/

