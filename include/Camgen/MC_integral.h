//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file MC_integral.h
    \brief Monte Carlo cross section class template.
 */

#ifndef MC_INTEGRAL_H_
#define MC_INTEGRAL_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Class template holding a cross section, its error and the error on the error. *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <iostream>

namespace Camgen
{
    /// Cross section data holder class template.

    template<class value_t>class MC_integral
    {
	public:

	    typedef value_t value_type;

	    /// Value of the variable.

	    value_type value;
	    
	    /// Error on the value.
	    
	    value_type error;

	    /// Error on the error.

	    value_type error_error;

	    /// Constructor.

	    MC_integral(const value_type& value_=(value_type)0,const value_type& error_=(value_type)0,const value_type& error_error_=(value_type)0):value(value_),error(error_),error_error(error_error_){}

	    /// Copy constructor.

	    MC_integral(const MC_integral<value_t>& other):value(other.value),error(other.error),error_error(other.error_error){}
	    /// Returns the error percentage

	    value_type error_percent() const
	    {
		return (value_type)100*error/value;
	    }

	    /// Returns the error on the error percentage

	    value_type error_error_percent() const
	    {
		return (value_type)100*error_error/error;
	    }

	    /// Adds a cross section to the data

	    MC_integral<value_t>& operator += (const MC_integral<value_t>& other)
	    {
		value+=other.value;
		error=std::sqrt(error*error+other.error*other.error);
		error_error=std::sqrt(error_error*error_error+other.error_error*other.error_error);
	    }

	    /// Subtracts a cross section from the data

	    MC_integral<value_t>& operator -= (const MC_integral<value_t>& other)
	    {
		value-=other.value;
		error=std::sqrt(error*error+other.error*other.error);
		error_error=std::sqrt(error_error*error_error+other.error_error*other.error_error);
	    }

	    /// Multiplies the cross section by a constant.

	    MC_integral<value_t>& operator *= (const value_t& c)
	    {
		value*=c;
		error*=c;
		error_error*=c;
	    }

	    /// Divides the cross section by a constant.

	    MC_integral<value_t>& operator /= (const value_t& c)
	    {
		value/=c;
		error/=c;
		error_error/=c;
	    }
    };

    /// Addition operator.

    template<class value_t> MC_integral<value_t> operator + (const MC_integral<value_t>& xs1,const MC_integral<value_t>& xs2)
    {
	MC_integral<value_t> result(xs1);
	return result+=xs2;
    }

    /// Subtraction operator.

    template<class value_t> MC_integral<value_t> operator - (const MC_integral<value_t>& xs1,const MC_integral<value_t>& xs2)
    {
	MC_integral<value_t> result(xs1);
	return result-=xs2;
    }

    /// Multiplication operators.

    template<class value_t> MC_integral<value_t> operator * (const value_t& c,const MC_integral<value_t>& xs)
    {
	MC_integral<value_t> result(xs);
	return result*=c;
    }
    template<class value_t> MC_integral<value_t> operator * (const MC_integral<value_t>& xs,const value_t& c)
    {
	MC_integral<value_t> result(xs);
	return result*=c;
    }

    /// Division operator.

    template<class value_t> MC_integral<value_t> operator / (const MC_integral<value_t>& xs,const value_t& c)
    {
	MC_integral<value_t> result(xs);
	return result/=c;
    }

    /// Output stream operator.

    template<class value_t>std::ostream& operator << (std::ostream& os,const MC_integral<value_t>& xs)
    {
	os<<xs.value<<" ± "<<xs.error<<" ± "<<xs.error_error;
	return os;
    }

    /// Comparison for cross sctions.

    template<class value_t>bool equals(const MC_integral<value_t>& first,const MC_integral<value_t>& second)
    {
	value_t tolerance(2.575829303549);
	value_t first_error=std::abs(first.error)+tolerance*std::abs(first.error_error);
	value_t second_error=std::abs(second.error)+tolerance*std::abs(second.error_error);
	if(std::abs(first.value-second.value)<=tolerance*(first_error+second_error))
	{
	    return true;
	}
	return false;
    }
}

#endif /*CAMGEN_MC_INTEGRAL_H_*/
