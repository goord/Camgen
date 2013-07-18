//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <cstdio>
#include <Camgen/spin.h>

namespace Camgen
{
    /* Trivial contructor, creating the zero spin: */

    spin::spin():val(0){}

    /* Creating an integer spin from an integer: */

    spin::spin(const int& n):val(n<<1){}
	    
    /* Copy constructor: */
	    
    spin::spin(const spin& other):val(other.val){}

    /* Function returning an integer spin: */

    spin spin::integer(int n)
    {
	return spin(n);
    }

    /* Function returning a half-integer spin (however, if the argument is even, the
     * result is an integer spin with value n/2: */

    spin spin::half_integer(int n)
    {
	spin s;
	s.val=n;
	return s;
    }

    /* Fermionic spin property: */
	    
    bool spin::fermionic() const
    {
	return ((val & 1)==1);
    }

    /* Bosonic spin property: */

    bool spin::bosonic() const
    {
	return ((val & 1)==0);
    }

    /* Function returning twice the spin value: */

    int spin::twice() const
    {
	return val;
    }

    /* Function computing the Lorentz rank of a tensor describing a particle
     * with the current spin: */

    unsigned spin::Lorentz_rank() const
    {
	if(val < 0)
	{
	    return -val/2;
	}
	return val/2;
    }

    /* Function computing the Dirac rank of a tensor describing a particle with
     * the current spin: */

    unsigned spin::Dirac_rank() const
    {
	if(val < 0)
	{
	    return (-val)%2;
	}
	return val%2;
    }

    /* Smallest integer bigger than or equal to the spin: */

    int spin::ceil() const
    {
	if(val>0)
	{
	    return (val+1)/2;
	}
	if(val<0)
	{
	    return (val-1)/2;
	}
	return 0;
    }

    /* Biggest integer smaller than or equal to the spin: */

    int spin::floor() const
    {
	return val/2;
    }

    /* Assignment operator: */

    spin& spin::operator = (const spin s)
    {
	val=s.val;
	return *this;
    }
	    
    /* Assignment operator from integer: */
	    
    spin& spin::operator = (const int n)
    {
	val=(n<<1);
	return *this;
    }

    /* Pre-increment operator: */

    spin& spin::operator ++ ()
    {
	val+=2;
	return *this;
    }

    /* Post-decrement operator: */

    spin spin::operator ++ (int) const
    {
	spin clone(*this);
	return ++clone;
    }

    /* Pre-decrement operator: */

    spin& spin::operator -- ()
    {
	val-=2;
	return *this;
    }

    /* Post-decrement operator: */

    spin spin::operator -- (int) const
    {
	spin clone(*this);
	return --clone;
    }
	    
    /* Self-addition by another spin: */
	    
    spin& spin::operator += (const spin s)
    {
	val+=(s.val);
	return *this;
    }

    /* Self-addition by an integer: */

    spin& spin::operator += (const int n)
    {
	val+=(n<<1);
	return *this;
    }

    /* Self-subtraction by a spin: */

    spin& spin::operator -= (const spin s)
    {
	val-=(s.val);
	return *this;
    }

    /* Self-subtraction by an integer: */

    spin& spin::operator -= (const int n)
    {
	val-=(n<<1);
	return *this;
    }

    /* Comparison operators: */

    bool spin::operator < (const spin other) const
    {
	return val<other.val;
    }
    bool spin::operator <= (const spin other) const
    {
	return val<=other.val;
    }
    bool spin::operator == (const spin other) const
    {
	return val==other.val;
    }
    bool spin::operator != (const spin other) const
    {
	return val!=other.val;
    }
    bool spin::operator >= (const spin other) const
    {
	return val>=other.val;
    }
    bool spin::operator > (const spin other) const
    {
	return val>other.val;
    }
    bool spin::operator < (const int other) const
    {
	return val<(other<<1);
    }
    bool spin::operator <= (const int other) const
    {
	return val<=(other<<1);
    }
    bool spin::operator == (const int other) const
    {
	return val==(other<<1);
    }
    bool spin::operator != (const int other) const
    {
	return val!=(other<<1);
    }
    bool spin::operator >= (const int other) const
    {
	return val>=(other<<1);
    }
    bool spin::operator > (const int other) const
    {
	return val>(other<<1);
    }

    /* Printing function: */

    std::ostream& spin::print(std::ostream& os) const
    {
	if(bosonic())
	{
	    os<<(val>>1);
	}
	else
	{
	    char buf[33];
	    std::sprintf(buf,"%d/2",val);
	    os<<buf;
	}
	return os;
    }

    /* Outstream operator: */

    std::ostream& operator << (std::ostream& os,const spin s)
    {
	return s.print(os);
    }

    /* Comparison operator between integer and spin: */

    bool operator < (const int n,const spin s)
    {
	return s>n;
    }
    bool operator <= (const int n,const spin s)
    {
	return s>=n;
    }
    bool operator == (const int n,const spin s)
    {
	return s==n;
    }
    bool operator != (const int n,const spin s)
    {
	return s!=n;
    }
    bool operator >= (const int n,const spin s)
    {
	return s<=n;
    }
    bool operator > (const int n,const spin s)
    {
	return s<n;
    }

    /* Arithmetic operators: */

    spin operator + (const spin s1,const spin s2)
    {
	return spin(s1.val+s2.val);
    }
    spin operator + (const spin s1,const int s2)
    {
	return spin(s1.val+(s2<<1));
    }
    spin operator + (const int s1,const spin s2)
    {
	return spin((s1<<1)+s2.val);
    }
    spin operator - (const spin s1,const spin s2)
    {
	return spin(s1.val-s2.val);
    }
    spin operator - (const spin s1,const int s2)
    {
	return spin(s1.val-(s2<<1));
    }
    spin operator - (const int s1,const spin s2)
    {
	return spin((s1<<1)-s2.val);
    }
    spin operator - (const spin other)
    {
	return spin(-other.val);
    }
}

