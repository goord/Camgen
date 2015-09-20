//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ss_gen.h
    \brief Constrained invariant mass pair sampler class.
 */

#ifndef CAMGEN_SS_GEN_H_
#define CAMGEN_SS_GEN_H_

#include <Camgen/val_gen.h>
#include <Camgen/multi_channel.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Monte Carlo class generating a pair of invariant masses for which the sum is  *
 * smaller than a given value sqrt(s).                                           *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Abstract base class for the contrained composition of two invariant mass generators.

    template<class value_t,class rng_t>class s_pair_generator: public virtual MC_generator<value_t>
    {
	typedef MC_generator<value_t> base_type;

	public:

	    /* Type definitions: */
	    
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef value_generator<value_t,rng_t> s_generator_type;

	    /// Total invariant mass-squared, upper bound to the sum of both
	    /// generated invariants.
	    
	    const value_type& s;
	    
	    /// First invariant mass generator.
	    
	    s_generator_type* const s1_generator;

	    /// Second invariant mass generator.

	    s_generator_type* const s2_generator;

	    /// Constructor.

	    s_pair_generator(s_generator_type* s1_generator_,s_generator_type* s2_generator_,const value_type& s_):s(s_),s1_generator(s1_generator_),s2_generator(s2_generator_){}

	    /// Copy constructor.
	    
	    s_pair_generator(const s_pair_generator<value_t,rng_t>& other):MC_generator<value_t>(other),s(other.s),s1_generator(other.s1_generator),s2_generator(other.s2_generator){}

	    /// Destructor.

	    virtual ~s_pair_generator(){}

	    /// First invariant mass.

	    const value_type& s1() const
	    {
		return s1_generator->value();
	    }

	    /// First minimal invariant mass.

	    const value_type& s1_min() const
	    {
		return s1_generator->lower_bound();
	    }

	    /// First maximal invariant mass.

	    const value_type& s1_max() const
	    {
		return s1_generator->upper_bound();
	    }

	    /// Second invariant mass.

	    const value_type& s2() const
	    {
		return s2_generator->value();
	    }

	    /// Second minimal invariant mass.

	    const value_type& s2_min() const
	    {
		return s2_generator->lower_bound();
	    }

	    /// Second maximal invariant mass.

	    const value_type& s2_max() const
	    {
		return s2_generator->upper_bound();
	    }
    };

    /// Generates a constained pair of invariant masses asymmetrically.

    template<class value_t,class rng_t>class s_generator_composition: public s_pair_generator<value_t,rng_t>
    {
	typedef s_pair_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */
	    
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef value_generator<value_t,rng_t> s_generator_type;

	    /* Constructor: */

	    s_generator_composition(s_generator_type* s1_generator_,s_generator_type* s2_generator_,const value_type& s_):base_type(s1_generator_,s2_generator_,s_){}

	    /* Copy constructor: */
	    
	    s_generator_composition(const s_generator_composition<value_t,rng_t>& other):base_type(other){}

	    /* Generation method implementation: */

	    bool generate()
	    {
		value_type s1max=this->s2_min()+this->s-2*std::sqrt(this->s2_min()*this->s);
		this->s1_generator->set_upper_bound(s1max);
		if(!this->s1_generator->generate())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type s2max=this->s1()+this->s-2*std::sqrt(this->s1()*this->s);
		this->s2_generator->set_upper_bound(s2max);
		if(!this->s2_generator->generate())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->weight()=this->s1_generator->weight()*this->s2_generator->weight();
		return true;
	    }

	    /* Weight evaluation method implementation: */

	    bool evaluate_weight()
	    {
		value_type s1max=this->s2_min()+this->s-2*std::sqrt(this->s2_min()*this->s);
		this->s1_generator->set_upper_bound(s1max);
		if(!this->s1_generator->evaluate_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type s2max=this->s1()+this->s-2*std::sqrt(this->s1()*this->s);
		this->s2_generator->set_upper_bound(s2max);
		if(!this->s2_generator->evaluate_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		
		this->weight()=this->s1_generator->weight()*this->s2_generator->weight();
		return true;
	    }
    };

    /// Generates a pair of invariant masses, with zero weights if they do not
    /// respect the momentum conservation constraint:
    
    template<class value_t,class rng_t>class symmetric_s_pair_generator: public s_pair_generator<value_t,rng_t>
    {
	typedef s_pair_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */
	    
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef value_generator<value_t,rng_t> s_generator_type;

	    /* Constructor: */

	    symmetric_s_pair_generator(s_generator_type* s1_generator_,s_generator_type* s2_generator_,const value_type& s_):base_type(s1_generator_,s2_generator_,s_){}

	    /* Copy constructor: */
	    
	    symmetric_s_pair_generator(const symmetric_s_pair_generator<value_t,rng_t>& other):base_type(other){}

	    /* Generation method implementation: */

	    bool generate()
	    {
		value_type s1max=this->s2_min()+this->s-2*std::sqrt(this->s2_min()*this->s);
		this->s1_generator->set_upper_bound(s1max);
		value_type s2max=this->s1_min()+this->s-2*std::sqrt(this->s1_min()*this->s);
		this->s2_generator->set_upper_bound(s2max);

		bool q=true;
		q&=this->s1_generator->generate();
		q&=this->s2_generator->generate();
		
		if(!q)
		{
		    this->weight()=(value_type)0;
		    return false;
		}

		value_type ssum=this->s1()+this->s2()+2*std::sqrt(this->s1()*this->s2());
		if(this->s<ssum)
		{
		    this->weight()=(value_type)0;
		    return false;
		}

		this->weight()=this->s1_generator->weight()*this->s2_generator->weight();
		return true;
	    }

	    /* Weight evaluation method implementation: */

	    bool evaluate_weight()
	    {
		value_type s1max=this->s2_min()+this->s-2*std::sqrt(this->s2_min()*this->s);
		this->s1_generator->set_upper_bound(s1max);
		value_type s2max=this->s1_min()+this->s-2*std::sqrt(this->s1_min()*this->s);
		this->s2_generator->set_upper_bound(s2max);

		bool q=true;
		q&=this->s1_generator->evaluate_weight();
		q&=this->s2_generator->evaluate_weight();
		
		if(!q)
		{
		    this->weight()=(value_type)0;
		    return false;
		}

		value_type ssum=this->s1()+this->s2()+2*std::sqrt(this->s1()*this->s2());
		if(this->s<ssum)
		{
		    this->weight()=(value_type)0;
		    return false;
		}

		this->weight()=this->s1_generator->weight()*this->s2_generator->weight();
		return true;
	    }
    };

    /// Generates a constrained pair of invariant masses with multichannel
    /// adaptive symmetry.

    template<class value_t,class rng_t>class symmetric_s_generator_composition: public s_pair_generator<value_t,rng_t>, multi_channel<value_t,rng_t>
    {
	public:

	    /* Type definitions: */
	    
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef value_generator<value_t,rng_t> s_generator_type;

	    /// Constructor.

	    symmetric_s_generator_composition(s_generator_type* s1_generator_,s_generator_type* s2_generator_,const value_type& s_):s_pair_generator<value_t,rng_t>(s1_generator_,s2_generator_,s_),s1s2_generator(s1_generator_,s2_generator_,s_),s2s1_generator(s2_generator_,s1_generator_,s_)
	    {
		this->add_generator(&s1s2_generator);
		this->add_generator(&s2s1_generator);
	    }

	    /// Copy constructor.
	    
	    symmetric_s_generator_composition(const symmetric_s_generator_composition<value_t,rng_t>& other):s_pair_generator<value_t,rng_t>(other),multi_channel<value_t,rng_t>(other),s1s2_generator(other.s1s2_generator),s2s1_generator(other.s2s1_generator){}

	private:

	    s_generator_composition<value_t,rng_t> s1s2_generator;
	    s_generator_composition<value_t,rng_t> s2s1_generator;
    };
}

#endif /*CAMGEN_SS_GEN_H_*/

