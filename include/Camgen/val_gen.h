//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file val_gen.h
    \brief Base class for 1D continuous Monte Carlo generators.
 */

#ifndef CAMGEN_VAL_GEN_H_
#define CAMGEN_VAL_GEN_H_

#include <limits>
#include <Camgen/utils.h>
#include <Camgen/logstream.h>
#include <Camgen/debug.h>
#include <Camgen/MC_gen.h>
#include <Camgen/MC_int_base.h>
#include <Camgen/rn_strm.h>
#include <Camgen/histogram.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and implementation of the value Monte Carlo generator base class. *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Momentum channel base class for phase space generators.

    template<class value_t,class rng_t>class value_generator: public MC_generator<value_t>, public MC_integrator_base<value_t>
    {
	typedef MC_generator<value_t> base_type;

	public:

	    /* Type definitions: */
	    
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;

	    /* Public constructors: */
	    /*----------------------*/

	    /// Default constructor.
	    
	    value_generator():min_val(-std::numeric_limits<value_t>::infinity()),max_val(std::numeric_limits<value_t>::infinity()),norm(0),val(new value_type),alloc_val(true){}

	    /// Copy constructor.

	    value_generator(const value_generator<value_t,rng_t>& other):MC_generator<value_t>(other),min_val(other.min_val),max_val(other.max_val),norm(0),alloc_val(other.alloc_val)
	    {
		val=alloc_val?(new value_type(*(other.val))):(other.val);
	    }

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    virtual ~value_generator()
	    {
		if(alloc_val)
		{
		    delete val;
		}
	    }

	    /* Public modifiers: */
	    /*--------------------*/

	    /// Sets the invariant mass address to the argument pointer.
	    
	    virtual void set_value(value_type* val_)
	    {
		if(alloc_val)
		{
		    delete val;
		    alloc_val=false;
		}
		val=val_;
	    }

	    /// Sets the minimal generated invariant mass-squared.

	    bool set_lower_bound(const value_type& min_val_)
	    {
		min_val=min_val_;
		return refresh_lower_bound();
	    }

	    /// Sets the maximal generated invariant mass-squared.

	    bool set_upper_bound(const value_type& max_val_)
	    {
		max_val=max_val_;
		return refresh_upper_bound();
	    }

	    /// Sets the invariant mass-squared range.

	    bool set_bounds(const value_type& min_val_,const value_type& max_val_)
	    {
		min_val=min_val_;
		max_val=max_val_;
		return refresh_bounds();
	    }

	    /// Generation method.

	    virtual bool generate()
	    {
		if(!normalisable())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		*val=map(rn_stream::throw_number());
		this->weight()=get_fast_weight();
		return true;
	    }

	    /// Weight evaluation method.

	    virtual bool evaluate_weight()
	    {
		if(!normalisable() or !within_bounds())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->weight()=get_fast_weight();
		return true;
	    }

	    /// Trivial update method.

	    virtual void update(){}

	    /// Re-evaluates internal parameters from external pointer values. Returns
	    /// false if the input parameters are invalid.

	    virtual bool refresh_params()
	    {
		refresh_norm();
		return normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on the lower bound. Returns false if
	    /// the input parameters are invalid.

	    virtual bool refresh_lower_bound()
	    {
		refresh_norm();
		return normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on the upper bound. Returns false if
	    /// the input parameters are invalid.

	    virtual bool refresh_upper_bound()
	    {
		refresh_norm();
		return normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on uper and lower
	    /// bounds. Returns false if the input parameters are invalid.

	    virtual bool refresh_bounds()
	    {
		refresh_lower_bound();
		refresh_upper_bound();
		return normalisable();
	    }

	    /* Public constant methods: */
	    /*--------------------------*/

	    /// Instance cloning method. Returns the NULL pointer if not
	    /// overridden.

	    virtual value_generator<value_t,rng_t>* clone() const
	    {
		return NULL;
	    }

	    /// Maps the first argument, a random number in [0,1] to a value.

	    virtual value_type map(const value_type&) const=0;

	    /// Maps the value argument to a random number in [0,1].

	    virtual value_type inverse_map(const value_type&) const=0;

	    /// Returns whether the requested density is nomalisable within
	    /// current interval.

	    bool normalisable() const
	    {
		return (norm>0 and norm!=std::numeric_limits<value_t>::infinity());
	    }

	    /// Returns whether the value range is bounded.

	    bool is_bounded() const
	    {
		return !((min_val==-std::numeric_limits<value_t>::infinity()) or (max_val==std::numeric_limits<value_t>::infinity()));
	    }

	    /// Returns the probability density at the given invariant mass.

	    value_type density() const
	    {
		if(this->within_bounds())
		{
		    return this->norm/this->weight();
		}
		return 0;
	    }

	    /// Returns a reference to the invariant mass-squared.
	    
	    value_type& value()
	    {
		return *val;
	    }

	    /// Returns a constant reference to the invariant mass-squared.
	    
	    const value_type& value() const
	    {
		return *val;
	    }

	    /// Returns the invariant mass-squared address.

	    value_type* get_value()
	    {
		return val;
	    }

	    /// Returns the constant invariant mass address.

	    const value_type* get_value() const
	    {
		return val;
	    }

	    /// Returns the minimal invariant mass-squared value.

	    const value_type& lower_bound() const
	    {
		return min_val;
	    }

	    /// Returns the maximal invariant mass-squared value.

	    const value_type& upper_bound() const
	    {
		return max_val;
	    }

	    /// Returns whether the invariant mass is within limits.

	    bool within_bounds() const
	    {
		return (greater_equal(*val,min_val) and smaller_equal(*val,max_val));
	    }

	    /// Evaluates the (private) normalisation factor.

	    virtual void refresh_norm()=0;

	    /// Returns the weight of the current configuration without bound-
	    /// or norm-checking.

	    virtual value_type get_fast_weight() const=0;

	    /* Testing methods: */
	    /*------------------*/

	    /// Checking function.

	    bool check()
	    {
		if(!within_bounds())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed check: invariant mass "<<value()<<" not within range ["<<lower_bound()<<','<<upper_bound()<<"]"<<endlog;
		    return false;
		}
		value_type w=this->weight();
		this->evaluate_weight();
		if(w!=this->weight())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed check: recomputed weight "<<this->weight()<<" does not equal generation weight "<<w<<endlog;
		    return false;
		}
		return true;
	    }

	    /// Testrun of N_events generations returning a histogram with
	    /// N_bins bins.

	    //TODO: get this out.
	    histogram<value_type> testrun(size_type N_events,size_type N_bins)
	    {
		histogram<value_type>hist(val,&(this->weight()),N_events);
		for(size_type n=0;n<N_events;++n)
		{
		    if(this->generate())
		    {
			check();
		    }
		    else
		    {
			this->weight()=(value_type)0;
		    }
		    hist.store();
		}
		hist.make(N_bins);
		return hist;
	    }

	    /* Serialization: */
	    /*----------------*/
	    
	    /// Returns a string containing the type.

	    virtual std::string type() const=0;

	    /// Returns the detailed type of the generator.

	    virtual std::string detailed_type() const
	    {
		return type();
	    }

	    /// Printing method:

	    virtual std::ostream& print(std::ostream& os) const
	    {
		os<<type();
		return os;
	    }

	protected:
	    
	    /// Minimal and maximal invariant mass-squared.

	    value_type min_val,max_val;

	    /// Normalisation constant.

	    value_type norm;

	private:

	    /* Pointer to the generated invariant mass: */

	    value_type* val;

	    /* Boolean denoting whether the invariant mass address was allocated
	     * dynamically: */

	    bool alloc_val;
    };
}

#endif /*CAMGEN_VAL_GEN_H_*/

