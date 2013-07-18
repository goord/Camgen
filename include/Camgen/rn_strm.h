//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file rn_strm.h
    \brief random-number generator wrapper class.
 */

#ifndef CAMGEN_RN_STRM_H_
#define CAMGEN_RN_STRM_H_

#include <iostream>
#include <cmath>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Wrapper class for random number generators. The first template parameter      *
 * denotes the numerical type the integers will be casted to, the second         *
 * template parameter denotes the random number generator class used. This class *
 * should contain the const static data members min_value and max_value,         *
 * denoting the minimum and maximum integers thrown and the operator (void) for  *
 * a throw.                                                                      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Random number generator wrapper. Converts an integer random number
    /// generator rng_t to a stream of numerical type value_t.

    template<class value_t,class rng_t>class random_number_stream
    {
	public:

	    /// Numerical type to convert the random numbers to.

	    typedef value_t value_type;
	    
	    /// Reference type definition to the random number generator.
	    
	    typedef rng_t rn_engine;

	    /// Result (integer type) of the random number generator.

	    typedef typename rng_t::result_type int_type;

	    /// Minimal value thrown by rng_t.

	    static const int_type min_value=rng_t::min_value;
	    
	    /// Maximal value thrown by rng_t.
	    
	    static const int_type max_value=rng_t::max_value;
    
	    /// Range of values thrown by rng_t (assuming max_value is positive).

	    static const int_type range_value=max_value-min_value;
	    
	    /// Range of values thrown by rng_t, as a value_type.
	    
	    static const value_type range;
	    
	    /// Returns the number of calls to rng.

	    static std::size_t number_of_throws()
	    {
		return counter;
	    }

	    /// Returns a random integer between min- and max.

	    static int_type throw_integer()
	    {
		init();
		++counter;
		return (*rng)();
	    }

	    /// Returns randomly 0 or 1.

	    static int_type throw_coin()
	    {
		init();
		++counter;
		return ((*rng)()%2);
	    }

	    /// Returns a random integer from 0 to max.

	    static int_type throw_dice(int_type max)
	    {
		init();
		++counter;
		if(max>0)
		{
		    return (*rng)()%max;
		}
		else
		{
		    return -((*rng)()%(-max));
		}
	    }

	    /// Returns a random integer between min and max.

	    static int_type throw_dice(int_type min,int_type max)
	    {
		init();
		++counter;
                return (*rng)()%std::abs(max-min)+std::min(min,max);
	    }

	    /// Returns a random floating-point number between 0 and 1.

	    static value_type throw_number()
	    {
		init();
		++counter;
		return (value_type)((*rng)()-min_value)/range;
	    }

	    /// Returns a random floating-point number between 0 and the argument.

	    static value_type throw_number(const value_type& max)
	    {
		init();
		++counter;
		return (value_type)((*rng)()-min_value)*max/range;
	    }
	    
	    /// Returns a random floating-point number between min and max.
	    
	    static value_type throw_number(const value_type& min,const value_type& max)
	    {
		init();
		++counter;
		return (value_type)((*rng)()-min_value)*std::abs(max-min)/range+std::min(min,max);
	    }

	    /// Resets the engine.

	    static void reset_engine()
	    {
		delete rng;
		rng=new rn_engine;
		counter=0;
	    }

	    value_type operator()(const value_type& min,const value_type& max) const
	    {
		return throw_number(min,max);
	    }
	    value_type operator()(const value_type& max) const
	    {
		return throw_number(max);
	    }
	    value_type operator()() const
	    {
		return throw_number();
	    }
	
	private:

	    /* Lazy instantiator: */

	    static void init()
	    {
		if(rng==NULL)
		{
		    rng=new rn_engine;
		}
	    }

	    /* Static random number generator instance: */

	    static rn_engine* rng;
	    
	    /* Call counter: */
	    
	    static std::size_t counter;
    };

    template<class value_t,class rng_t>const typename random_number_stream<value_t,rng_t>::int_type random_number_stream<value_t,rng_t>::min_value;
    template<class value_t,class rng_t>const typename random_number_stream<value_t,rng_t>::int_type random_number_stream<value_t,rng_t>::max_value;
    template<class value_t,class rng_t>const typename random_number_stream<value_t,rng_t>::int_type random_number_stream<value_t,rng_t>::range_value;
    template<class value_t,class rng_t>typename random_number_stream<value_t,rng_t>::rn_engine* random_number_stream<value_t,rng_t>::rng=NULL;
    template<class value_t,class rng_t>const typename random_number_stream<value_t,rng_t>::value_type random_number_stream<value_t,rng_t>::range=random_number_stream<value_t,rng_t>::range_value;
    template<class value_t,class rng_t>std::size_t random_number_stream<value_t,rng_t>::counter=0;
}

#endif /*CAMGEN_RN_STRM_H_*/

