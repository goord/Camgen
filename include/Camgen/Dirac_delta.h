//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DIRAC_DELTA_H_
#define CAMGEN_DIRAC_DELTA_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and implementation of the Dirac delta type invariant mass *
 * sampler, for on-shell particle propagator sampling.                   *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <sstream>
#include <Camgen/utils.h>
#include <Camgen/s_gen.h>
#include <Camgen/s_int.h>

namespace Camgen
{
    /* On-shell momentum generator channel class. */

    template<class value_t,class rng_t>class Dd_s_generator: public s_generator<value_t,rng_t>
    {
	friend class s_integrator<value_t,rng_t>;
	typedef s_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Public data members: */
	    /*----------------------*/

	    /* Mass address: */

	    const value_type* const m;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Default constructor with (optional) mass address argument. */

	    Dd_s_generator(const value_type* m_=NULL):m(m_),mm(squareval(m)){}

	    /* Constructor with momentum and (optional) mass address argument.
	     * */

	    Dd_s_generator(value_type* s_,const value_type* m_=NULL):base_type(s_),m(m_),mm(squareval(m)){}

	    /* Copy constructor: */

	    Dd_s_generator(const Dd_s_generator<value_t,rng_t>& other):base_type(other),m(other.m),mm(other.mm){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Parameter refresher function. */

	    bool refresh_params()
	    {
		mm=squareval(m);
		refresh_norm();
		return this->normalisable();
	    }

	    /* Generation method. */

	    bool generate()
	    {
		if(!this->normalisable())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->s()=mm;
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Virtual cloning method. */

	    virtual Dd_s_generator<value_t,rng_t>* clone() const
	    {
		return new Dd_s_generator<value_t,rng_t>(*this);
	    }

	    /* Mass output. */

	    value_type mass() const
	    {
		return (m==NULL)?(value_type)0:(*m);
	    }

	    /* Mass-squared output. */

	    value_type mass2() const
	    {
		return mm;
	    }

	    /* Mapping implementation: */

	    value_type map(const value_type& r) const
	    {
		return mm;
	    }

	    /* Inverse mapping implementation: */

	    value_type inverse_map(const value_type& s) const
	    {
		return s;
	    }

	    /* Double-dispatch integrator functions. */
	    
	    value_type integrate_with(const s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return gen->integrate_with(this,sqrts);
	    }
	    value_type integrate_with(const BW_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		typedef BW_s_generator<value_t,rng_t> type;
		return s_integrator<value_t,rng_t>::template integrate<type>(static_cast<const inversion_s_generator<value_t,rng_t,type>*>(gen),this,sqrts);
	    }
	    value_type integrate_with(const pl_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		typedef pl_s_generator<value_t,rng_t> type;
		return s_integrator<value_t,rng_t>::template integrate<type>(static_cast<const inversion_s_generator<value_t,rng_t,type>*>(gen),this,sqrts);
	    }
	    value_type integrate_with(const uni_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(gen,this,sqrts);
	    }
	    value_type integrate_with(const Dd_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(this,gen,sqrts);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Overridden printing function: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<"Dd("<<mass()<<')';
		return os;
	    }

	    /* Outputs the type of the generator. */

	    std::string type() const
	    {
		return "Dd";
	    }

	    /* Outputs the type of the generator. */

	    std::string detailed_type() const
	    {
		std::stringstream ss;
		ss<<"Dd("<<mass()<<')';
		return ss.str();
	    }

	protected:

	    /* Normalisation refresher: */

	    void refresh_norm()
	    {
		bool q=(mm>=this->s_min() and mm<=this->s_max());
		this->norm=q?(value_type)1:(value_type)0;
	    }

	    /* Fast weight evaluation method: */

	    value_type get_fast_weight() const
	    {
		return equals(this->s(),mm)?(value_type)1:(value_type)0;
	    }

	private:

	    /* Square utility: */

	    static value_type squareval(const value_type* x)
	    {
		return (x==NULL)?(value_type)0:(*x)*(*x);
	    }

	    /* Squared mass data member: */

	    value_type mm;
    };
}

#endif /*CAMGEN_DIRAC_DELTA_H_*/

