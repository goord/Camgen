//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_BIT_STRING_H_
#define CAMGEN_BIT_STRING_H_

#include <iostream>
#include <bitset>
#include <Camgen/combs.h>
#include <Camgen/lower_binom.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Static-size bitstring class definition and declaration. The bit-string type   *
 * is derived from the std::bitset class template, and some algorithmic          *
 * functionality is added. The bitstrings are used to label the momentum channel *
 * of an off-shell current in Camgen. The ordering of the strings is            *
 * lexicographical w.r.t. the bits set, which is not optimal from a              *
 * computational point of view, but is more insightful to the user.              *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition of the static-size bitstring class: */
    
    template<std::size_t N>class bit_string: public std::bitset<N>
    {
	public:

	    /* Inherited type definitions: */

	    typedef std::size_t size_type;

	    /* Default constructor, all bits will be initialised to 0: */

	    bit_string(){}

	    /* Copy constructor: */

	    bit_string(const std::bitset<N>& b):std::bitset<N>(b){}
	    
	    /* Ordered constructor, with higher complexity: */
	    
	    bit_string(size_type m)
	    {
		if(m>0)
		{
		    /* Counting the number of bitstrings of length N: */
		    
		    size_type M=(1<<N);

		    /* Subtract the number of possible bitstrings until m is
		     * smaller then M: */

		    m%=M;
		    
		    /* Create a lower-binomial class object, holding the order
		     * (corresponding to the number of bits set) and the
		     * difference (corresponding to the offset of the bitstring
		     * to-be-constructed w.r.t. the first C-level bitstring: */ 
		    
		    lower_binomial B(N,m);
		    size_type C=B.order;
		    size_type D=B.difference;
		    
		    if(m!=(M-1))
		    {
			size_type M=N;
			
			/* Bit-flipping procedure. C counts the number of bits
			 * that still need to be set. M denotes the length of
			 * the part that needs to be changed. D denotes the
			 * offset w.r.t. the bitstring with C bits starting at
			 * the (N-M)-th place: */

			do
			{
			    if(C==0)
			    {
				this->flip(D+N-M);
				M=0;
			    }
			    else if(D<binomial(M-1,C))
			    {
				this->flip(N-M);
				--M;
				--C;
			    }
			    else
			    {
				this->reset(N-M);
				D-=binomial(M-1,C);
				--M;
			    }

			}
			while(M>0);
		    }

		    /* If m equals the binomial sum minus one (m=0 yields a
		     * string with the first bit set), set all bits: */

		    else
		    {
			this->flip();
		    }
		}
	    }

	    /* Assignment operator: */

	    bit_string<N>& operator = (const std::bitset<N>& b)
	    {
		if(this != &b)
		{
		    this->std::bitset<N>::operator=(b);
		}
		return *this;
	    }

	    /* Bit shift operators, reversed w.r.t. operators in the std::bitset
	     * base class: */

	    bit_string<N>& operator <<=(size_type n)
	    {
		std::bitset<N>::operator>>=(n);
		return *this;
	    }
	    bit_string<N> operator << (size_type n)
	    {
		bit_string result(*this);
		return result<<=n;
	    }
	    bit_string<N>& operator >>=(size_type n)
	    {
		std::bitset<N>::operator<<=(n);
		return *this;
	    }
	    bit_string<N> operator >> (size_type n)
	    {
		bit_string result(*this);
		return result>>=n;
	    }

	    /* Count the unset bits: */

	    size_type count_zeros() const
	    {
		return N-this->count();
	    }
	    
	    /* Computing the next bit string in the lexicographical ordering: */
	    
	    bit_string<N>& next()
	    {
		size_type bits=0;
		size_type zeros=0;

		/* Compute the number of bits at the end of the string, and the
		 * number of unset bits preceding those, and store the results
		 * in 'bits' and 'zeros' respectively: */

		inspect_back(bits,zeros);

		size_type C=this->count();
		
		/* If we have the trivial string, set the first bit: */
		
		if(C==0)
		{
		    this->set(0);
		    return *this;
		}

		/* If the last bit is zero, reset the last bit set and set the
		 * proceeding one: */

		if(bits==0)
		{
		    this->flip(N-zeros-1);
		    this->flip(N-zeros);
		    return *this;
		}
		
		/* If all bits are sitting at the end of the string, construct
		 * the first next-level bitstring: */
		
		if(bits==C and bits!=N)
		{
		    for(size_type i=0;i<bits+1;++i)
		    {
			this->set(i);
		    }
		    for(size_type i=bits+1;i<N;++i)
		    {
			this->reset(i);
		    }
		    return *this;
		}

		/* If all bits are set, reset the bit string: */

		if(bits==N)
		{
		    this->flip();
		    return *this;
		}
		
		/* Else, reset the first bit preceding the zeros, set all the
		 * subsequent zero bits and reset all the remaining set bits: */
		
		this->flip(N-bits-zeros-1);
		for(size_type i=N-bits-zeros;i<N-zeros+1;++i)
		{
		    this->set(i);
		}
		for(size_type i=N-zeros+1;i<N;++i)
		{
		    this->reset(i);
		}
		return *this;
	    }
	    
	    /* Computing the previous bitstring in the lexicographical ordering
	     * by computing the next conjugate bitstring: */
	    
	    bit_string<N>& previous()
	    {
		this->flip();
		next();
		this->flip();
		return *this;
	    }

	    /* Compute the rank of the bitstring w.r.t. to the lexicographical
	     * ordering: */

	    size_type to_integer() const
	    {
		std::bitset<N>B(*this);

		/* Check if any bit is set: */

		if(B.any())
		{

		    /* Start with the binomial sum of the lower level
		     * bitstrings: */

		    size_type n=binomial_sum(N,B.count()-1);
		    
		    /* Repeatedly shift the copy to the left, adding smaller based binomials
		     * for each zero first bit: */
		    
		    for(int M=N-1;M!=0;--M)
		    {
			if(B.any())
			{
			    if(!B[0])
			    {
				n+=binomial((unsigned)M,B.count()-1);
			    }
			    B>>=1;
			}
			else
			{
			    break;
			}
		    }
		    return n;
		}
		return 0;
	    }

	    /* Pre-increment operator, a substitute for next(): */

	    bit_string<N> operator ++ ()
	    {
		return next();
	    }

	    /* Post-increment operator: returns a copy and increments the
	     * original: */

	    bit_string<N>& operator ++ (int)
	    {
		bit_string<N>clone(*this);
		next();
		return clone;	
	    }

	    /* Multi-increment operator: */

	    bit_string<N>& operator += (int n)
	    {
		if(n>=0)
		{
		    for(int i=0;i<n;++i)
		    {
			next();
		    }
		    return *this;
		}
		for(int i=0;i<-n;++i)
		{
		    previous();
		}
		return *this;
	    }

	    /* Pre-decrement operator, a substitute for previous(): */

	    bit_string<N> operator -- ()
	    {
		return previous();	
	    }

	    /* Post-decrement operator: returns a copy and decrements the
	     * original: */

	    bit_string<N>& operator -- (int)
	    {
		bit_string<N>clone(*this);
		previous();
		return clone;	
	    }

	    /* Multi-decrement operator: */

	    bit_string<N>& operator -= (int n)
	    {
		if(n>=0)
		{
		    for(int i=0;i<n;++i)
		    {
			previous();
		    }
		    return *this;
		}
		for(int i=0;i<-n;++i)
		{
		    next();
		}
		return *this;
	    }

	    /* Convolution operator: whenever a bit is set, a bit of the
	     * argument is copied to that place: */

	    bit_string<N>& operator *=(const std::bitset<N>& b)
	    {
		size_type j=0;
		for(size_type i=0;i<N;++i)
		{
		    if((*this)[i])
		    {
			(*this)[i]=b[j];
			++j;
		    }
		}
		return *this;
	    }

	    /* Copy-convolution operator: */

	    bit_string<N> operator *(const std::bitset<N>& b) const
	    {
		bit_string<N>clone(*this);
		return clone*=(b);
	    }

	    /* Comprison operators: the less than operation is essentially a
	     * lexicographical comparison: */

	    bool operator < (const std::bitset<N>& b) const
	    {
		if(this->count()==b.count())
		{
		    size_type i=0;
		    while(((*this)[i]==b[i]) and i<N)
		    {
			++i;
		    }
		    if(i<N)
		    {
			return ((*this)[i] and !b[i]);
		    }
		    else
		    {
			return false;
		    }
		}
		else
		{
		    return (this->count() < b.count());
		}
	    }
	    bool operator <= (const std::bitset<N>& b) const
	    {
		return !(b<(*this));
	    }
	    bool operator > (const std::bitset<N>& b) const
	    {
		return b<(*this);
	    }
	    bool operator >= (const std::bitset<N>& b) const
	    {
		return !(*this < b);
	    }

	    /* Output method: */

	    std::ostream& print(std::ostream& os) const
	    {
		for(size_type i=0;i<N;++i)
		{
		    os<<(*this)[i];
		}
		return os;
	    }

	    /* Input method: */

	    std::istream& read(std::istream& is)
	    {
		std::bitset<N> b;
		is>>b;
		for(size_type i=0;i<N;++i)
		{
		    (*this)[i]=b[N-i-1];
		}
		return is;
	    }

	private:

	    /* Private utility method counting subsequently the number of
	     * adjacent set bits and the number zero adjacent bits at the end of
	     * the string: */

	    void inspect_back(size_type& bits,size_type& zeros)
	    {
		bool q1=true;
		bool q2=true;
		size_type M=N;
		do
		{
		    --M;
		    if(q1)
		    {
			if((*this)[M])
			{
			    ++bits;
			}
			else
			{
			    q1=false;
			}
		    }
		    else
		    {
			if(!(*this)[M])
			{
			    ++zeros;
			}
			else
			{
			    q2=false;
			}
		    }
		}
		while(q2 and M != 0);
		++zeros;
	    }
    };

    /* Output operator (note that the order of the bits in the output stream i
     * reversed w.r.t. std::bitset parent output): */

    template<std::size_t N>std::ostream& operator << (std::ostream& os,const bit_string<N>& b)
    {
	return b.print(os);
    }

    /* Input operator (note that the order of the bits in the input stream i
     * reversed w.r.t. std::bitset parent output): */

    template<std::size_t N>std::istream& operator >> (std::istream& is,bit_string<N>& b)
    {
	return b.read(is);
    }
}

#endif /*CAMGEN_BIT_STRING_H_*/

