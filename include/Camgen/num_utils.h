//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file num_utils.h
    \brief Numerical general utility functions.
 */

#ifndef CAMGEN_NUM_UTILS_H_
#define CAMGEN_NUM_UTILS_H_

/* * * * * * * * * * * * * * * * *
 * Numeric utilities for Camgen. *
 *                               *
 * * * * * * * * * * * * * * * * */

#include <limits>
#include <Camgen/num_config.h>

namespace Camgen
{
    /// Checks equality up to precision defined in num_config.h

    template<class T>bool equals(const T& a,const T& b)
    {
	if(std::abs(a-b)<numeric_configuration<T>::epsilon_abs)
	{
	    return true;
	}
	return std::abs(a-b)<=numeric_configuration<T>::epsilon_rel*std::max(std::abs(a),std::abs(b));
    }

    /// Compares less up to precision defined in num_config.h

    template<class T>bool smaller(const T& a,const T& b)
    {
	return (!equals(a,b) and a<b);
    }

    /// Compares greater up to precision defined in num_config.h

    template<class T>bool greater(const T& a,const T& b)
    {
	return (!equals(a,b) and a>b);
    }

    /// Compares less-or-equal up to precision defined in num_config.h

    template<class T>bool smaller_equal(const T& a,const T& b)
    {
	return (equals(a,b) or a<b);
    }

    /// Compares greater-or-equal up to precision defined in num_config.h

    template<class T>bool greater_equal(const T& a,const T& b)
    {
	return (equals(a,b) or a>b);
    }

    /// Checks equality of two sequences up to precision defined in num_config.h: */

    template<class T>bool equal_sequences(const T& first,const T& second)
    {
	if(first.size()==second.size())
	{
	    typename T::const_iterator it1=first.begin();
	    for(typename T::const_iterator it2=second.begin();it2 != second.end();++it2)
	    {
		if(!equals(*it1,*it2))
		{
		    return false;
		}	
		++it1;
	    }
	    return true;
	}
	else
	{
	    return false;
	}
    }

    /// Returns whether the argument is a number.
    
    template<class T>bool is_number(const T& x)
    {
	return x==x;
    }

    /// Returns whether the argument is a finite number.

    template<class T>bool is_finite_number(const T& x)
    {
	return (is_number(x) and x!=std::numeric_limits<T>::infinity() and x!=-std::numeric_limits<T>::infinity());
    }
}

#endif /*CAMGEN_NUM_UTILS_H_*/
