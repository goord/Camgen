//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file spin.h
    \brief spin class interface header.
 */

#ifndef CAMGEN_SPIN_H_
#define CAMGEN_SPIN_H_

#include <iostream>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the spin type class and methods. The spin type is implemented  *
 * as an integer, corresponding to a half-integer spin for odd numbers and       *
 * integer spin for even integers.                                               *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Declaration: */

    class spin;

    /* Declaration of arithmetic operators: */

    spin operator + (const spin s1,const spin s2);
    spin operator + (const spin s1,const int s2);
    spin operator + (const int s1,const spin s2);
    spin operator - (const spin s1,const spin s2);
    spin operator - (const spin s1,const int s2);
    spin operator - (const int s1,const spin s2);
    spin operator - (const spin s);
    
    /* Cast function template: */
    
    template<class T>T spin_cast(const spin s);

    /// Particle spin class in Camgen.
    /// The class stores integers and half-integers by twice the value, which is
    /// regular integer.

    class spin
    {
	/* Declaring the cast function template friend: */

	template<class T>friend T spin_cast(const spin);
	
	/* Declaring the arithmetic operators friends: */
	
	friend spin operator + (const spin s1,const spin s2);
	friend spin operator + (const spin s1,const int s2);
	friend spin operator + (const int s1,const spin s2);
	friend spin operator - (const spin s1,const spin s2);
	friend spin operator - (const spin s1,const int s2);
	friend spin operator - (const int s1,const spin s2);
	friend spin operator - (const spin s);

	public:

	    /// Trivial contructor, creating the zero spin.

	    spin();
	    
	    /// Copy constructor.
	    
	    spin(const spin&);

	    /// Constructor of an integer spin from an integer.

	    spin(const int&);
	    
	    /// Function returning an integer-spin instance.
	    
	    static spin integer(int);
	    
	    /// Function returning an half-integer-spin instance.
	    /// The value of the resulting spin is half the argument integer.
	    /// Hence i the argument is even, the result spin is bosonic.
	    
	    static spin half_integer(int);
	    
	    /// Returns whether the spin object blongs to a fermionic particle.
	    
	    bool fermionic() const;
	    
	    /// Returns whether the spin object blongs to a bosonic particle.
	    
	    bool bosonic() const;
	    
	    /// Returns twice the spin value.
	    
	    int twice() const;

	    /// Returns the Lorentz rank of a particle with this spin.

	    unsigned Lorentz_rank() const;

	    /// Returns the Dirac rank of a particle with this spin.

	    unsigned Dirac_rank() const;

	    /// Returns the smallest integer bigger than or equal to the spin value.

	    int ceil() const;

	    /// Returns the largest integer smaller than or equal to the spin value.

	    int floor() const;

	    /// Assignment operator.

	    spin& operator = (const spin);
	    
	    /// Assignment operator from integer.
	    
	    spin& operator = (const int);

	    /// Pre-increment operator.

	    spin& operator ++ ();

	    /// Post-increment operator.

	    spin operator ++ (int) const;

	    /// Pre-decrement operator.

	    spin& operator -- ();

	    /// Post-decrement operator.

	    spin operator -- (int) const;
	    
	    /// Inrements value by another spin.
	    
	    spin& operator += (const spin);

	    /// Increments the value by an integer.

	    spin& operator += (const int);
	    
	    /// Decrements value by another spin.
	    
	    spin& operator -= (const spin);

	    /// Decrements the value by an integer.

	    spin& operator -= (const int);

	    /// Comparison operator.

	    bool operator < (const spin) const;

	    /// Comparison operator.

	    bool operator < (const int) const;

	    /// Comparison operator.

	    bool operator <= (const spin) const;

	    /// Comparison operator.

	    bool operator <= (const int) const;

	    /// Comparison operator.

	    bool operator == (const spin) const;

	    /// Comparison operator.

	    bool operator == (const int) const;

	    /// Comparison operator.

	    bool operator != (const spin) const;

	    /// Comparison operator.

	    bool operator != (const int) const;

	    /// Comparison operator.

	    bool operator > (const spin) const;

	    /// Comparison operator.

	    bool operator > (const int) const;

	    /// Comparison operator.

	    bool operator >= (const spin) const;

	    /// Comparison operator.

	    bool operator >= (const int) const;

	    /// Output function.

	    std::ostream& print(std::ostream&) const;

	private:

	    /* Data integer: */

	    int val;
    };

    /// Overloaded streaming operator.

    std::ostream& operator << (std::ostream& os,const spin s);

    /// Comparison operator between integer and spin.

    bool operator < (const int,const spin);

    /// Comparison operator between integer and spin.

    bool operator <= (const int,const spin);

    /// Comparison operator between integer and spin.

    bool operator == (const int,const spin);

    /// Comparison operator between integer and spin.

    bool operator != (const int,const spin);

    /// Comparison operator between integer and spin.

    bool operator > (const int,const spin);

    /// Comparison operator between integer and spin.

    bool operator >= (const int,const spin);

    /// Spin casting function.

    template<class T>T spin_cast(const spin s)
    {
	return ((T)s.val)/((T)2);
    }
}

#endif /*CAMGEN_SPIN_H_*/

