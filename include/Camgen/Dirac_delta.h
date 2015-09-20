//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DIRAC_DELTA_H_
#define CAMGEN_DIRAC_DELTA_H_

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Dirac-delta type value generator definition *
 *                                             *
 * * * * * * * * * * * * * * * * * * * * * * * */

#include <sstream>
#include <Camgen/utils.h>
#include <Camgen/val_gen.h>

namespace Camgen
{
    /* On-shell momentum generator channel class. */

    template<class value_t,class rng_t>class Dirac_delta: public value_generator<value_t,rng_t>
    {
	typedef value_generator<value_t,rng_t> base_type;

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

	    Dirac_delta(const value_type* m_=NULL):m(m_),mm(squareval(m)){}

	    /* Copy constructor: */

	    Dirac_delta(const Dirac_delta<value_t,rng_t>& other):base_type(other),m(other.m),mm(other.mm){}

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
		this->value()=mm;
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Virtual cloning method. */

	    virtual Dirac_delta<value_t,rng_t>* clone() const
	    {
		return new Dirac_delta<value_t,rng_t>(*this);
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

	    /* Normalisation refresher: */

	    void refresh_norm()
	    {
		bool q=(mm>=this->lower_bound() and mm<=this->upper_bound());
		this->norm=q?(value_type)1:(value_type)0;
	    }

	    /* Weight evaluation method: */

	    bool evaluate_weight()
	    {
		if(!(this->normalisable()))
		{
		    return false;
		}
		this->weight()=get_fast_weight();
		return true;
	    }

	    /* Fast weight evaluation method: */

	    value_type get_fast_weight() const
	    {
		return equals(this->value(),mm)?(value_type)1:(value_type)0;
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

