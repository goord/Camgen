//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file vector.h
    \brief vector class template interface and implementation header.
 */

#ifndef CAMGEN_VECTOR_H_
#define CAMGEN_VECTOR_H_

#include <cstddef>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <Camgen/rn_strm.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the vector type in Camgen. It is used for the  *
 * momenta and spin vectors in the program. The size of the vector is a          *
 * compile-time constant. The implementation is similar to the boost::array      *
 * class...                                                                      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Fixed-size vector class.
    /** Because Camgen uses many vector objects which have their sizes
     * determined at compile-time, this class increases the performance as it
     * allows compilers to perform optimisations. */

    template<class T,std::size_t N>class vector
    {
	public:

	    /* Type definitions. The iterators are simply pointers, and the reverse
	     * iterators are defined with the help of the stl reverse_iterator class
	     * templates: */

	    typedef T 						value_type;
	    typedef T* 						iterator;
	    typedef const T*					const_iterator;
	    typedef std::reverse_iterator<iterator>		reverse_iterator;
	    typedef std::reverse_iterator<const_iterator> 	const_reverse_iterator;
	    typedef T&						reference;
	    typedef const T&					const_reference;
	    typedef std::size_t					size_type;
	    typedef std::ptrdiff_t				difference_type;

	    /// Compile-time vector size.

	    static const size_type static_size=N;

	    /// Returns the vector size.

	    static size_type size()
	    {
		return static_size;
	    }

	    /// Returns true if the size is zero, false otherwise.

	    static bool empty()
	    {
		return false;
	    }

	    /// Iterator to the first entry.

	    iterator begin()
	    {
		return data;
	    }

	    /// Constant iterator to the first entry.

	    const_iterator begin() const
	    {
		return data;
	    }

	    /// Iterator past the last entry.
	    
	    iterator end()
	    {
		return data+N;
	    }

	    /// Constant iterator past the last entry.

	    const_iterator end() const
	    {
		return data+N;
	    }

	    /// Reverse iterator to the last entry.
	    
	    reverse_iterator rbegin()
	    {
		return reverse_iterator(end());
	    }

	    /// Constant reverse iterator to the last entry.

	    const_reverse_iterator rbegin() const
	    {
		return const_reverse_iterator(end());
	    }

	    /// Reverse iterator before the first entry.

	    reverse_iterator rend()
	    {
		return reverse_iterator(begin());
	    }

	    /// Constant reverse iterator before the first entry.

	    const_reverse_iterator rend() const
	    {
		return const_reverse_iterator(begin());
	    }

	    /// Returns a reference to the n-th entry.

	    reference operator [] (size_type n)
	    {
		return data[n];
	    }

	    /// Returns a constant reference to the n-th entry.

	    const_reference operator [] (size_type n) const
	    {
		return data[n];
	    }

	    /// Returns a reference to the n-th entry.
	    /// Throws an std::out_of_range exception when the argument is
	    /// greater than or equel to the vector size. 

	    reference at(size_type n)
	    {
		if(n<static_size)
		{
		    return data[n];
		}
		else
		{
		    throw std::out_of_range("Attempt to access vector component out of range");
		}
	    }

	    /// Returns a constant reference to the n-th entry.
	    /// Throws an std::out_of_range exception when the argument is
	    /// greater than or equel to the vector size. 

	    const_reference at(size_type n) const
	    {
		if(n<static_size)
		{
		    return data[n];
		}
		else
		{
		    throw std::out_of_range("Attempt to access vector component out of range");
		}
	    }

	    /// Returns a reference to the first entry.

	    reference front()
	    {
		return data[0];
	    }

	    /// Returns a constant reference to the first entry.

	    const_reference front() const
	    {
		return data[0];
	    }

	    /// Returns a reference to the lst entry.

	    reference back()
	    {
		return data[static_size-1];
	    }

	    /// Returns a constant reference to the last entry.

	    const_reference back() const
	    {
		return data[static_size-1];
	    }

	    /// Assigns the value c to all entries and returns the instance.

	    vector<T,N>& assign(const value_type& c)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]=c;
		}
		return *this;
	    }

	    /// Swaps the content with vector other.

	    void swap(vector<T,N>& other)
	    {
		std::swap_ranges(data,data+N,other.data);
	    }

	    /// Adds c to all entries and returns the instance.

	    vector<T,N>& operator += (const value_type& c)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]+=c;
		}
		return *this;
	    }

	    /// Subtracts c from all entries and returns the instance.

	    vector<T,N>& operator -= (const value_type& c)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]-=c;
		}
		return *this;
	    }

	    /// Multiplies all entries with c and returns the instance.

	    vector<T,N>& operator *= (const value_type& c)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]*=c;
		}
		return *this;
	    }

	    /// Divides all entries by c and returns the instance.

	    vector<T,N>& operator /= (const value_type& c)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]/=c;
		}
		return *this;
	    }

	    /// Multiplies the instance element-wise with v and returns the instance.

	    vector<T,N>& operator *= (const vector<T,N>& v)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]*=v[i];
		}
		return *this;
	    }

	    /// Divides the instance element-wise by v and returns the instance.

	    vector<T,N>& operator /= (const vector<T,N>& v)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]/=v[i];
		}
		return *this;
	    }

	    /// Element-wise adds v to the instance and returns the instance.

	    vector<T,N>& operator += (const vector<T,N>& v)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]+=v[i];
		}
		return *this;
	    }

	    /// Element-wise subtracts v to from instance and returns the instance.

	    vector<T,N>& operator -= (const vector<T,N>& v)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]-=v[i];
		}
		return *this;
	    }

	    /// Fills all entries with random numbers provided by gen and returns the instance.

	    template<class rng_t>vector<T,N>& fill(random_number_stream<T,rng_t>& gen)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]=gen();
		}
		return *this;
	    }

	    /// Fills all entries with random numbers provided by gen within [0,max] and returns the instance.

	    template<class rng_t>vector<T,N>& fill(random_number_stream<T,rng_t>& gen,const T& max)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]=gen(max);
		}
		return *this;
	    }

	    /// Fills all entries with random numbers provided by gen within [min,max] and returns the instance.

	    template<class rng_t>vector<T,N>& fill(random_number_stream<T,rng_t>& gen,const T& min,const T& max)
	    {
		for(size_type i=0;i<static_size;++i)
		{
		    data[i]=gen(min,max);
		}
		return *this;
	    }
	    value_type data[N];
    };
    template<class T,std::size_t N>const std::size_t vector<T,N>::static_size;

    /* Explicit template specialisation for N==0, avoiding a compile time error reporting
     * the zero-sized array data[0]: */

    template<class T>class vector<T,0>
    {
	public:

	    /* Type definitions. The iterators are simply pointers, and the reverse
	     * iterators are defined with the help of the stl reverse_iterator class
	     * templates: */

	    typedef T 						value_type;
	    typedef T* 						iterator;
	    typedef const T*					const_iterator;
	    typedef std::reverse_iterator<iterator>		reverse_iterator;
	    typedef std::reverse_iterator<const_iterator> 	const_reverse_iterator;
	    typedef T&						reference;
	    typedef const T&					const_reference;
	    typedef std::size_t					size_type;
	    typedef std::ptrdiff_t				difference_type;

	    /* Compile-time constant zero size: */

	    static const size_type static_size=0;

	    /* Run-time constant zero size: */

	    static size_type size()
	    {
		return 0;
	    }

	    /* Empty query return true: */

	    static bool empty()
	    {
		return true;
	    }

	    /* Since no data exists, the iterators return a reinterpret cast of the vector
	     * object itself: */

	    iterator begin()
	    {
		return iterator(reinterpret_cast<T*>(this));
	    }
	    const_iterator begin() const
	    {
		return iterator(reinterpret_cast<const T*>(this));
	    }
	    iterator end()
	    {
		return begin();
	    }
	    const_iterator end() const
	    {
		return begin();
	    }
	    reverse_iterator rbegin()
	    {
		return reverse_iterator(end());
	    }
	    const_reverse_iterator rbegin() const
	    {
		return const_reverse_iterator(end());
	    }
	    reverse_iterator rend()
	    {
		return reverse_iterator(begin());
	    }
	    const_reverse_iterator rend() const
	    {
		return const_reverse_iterator(begin());
	    }

	    /* Random access operators throw exceptions: */

	    reference operator [] (size_type n)
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }
	    const_reference operator [] (size_type n) const
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }
	    reference at(size_type n)
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }
	    const_reference at(size_type n) const
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }
	    /* First and last element access results in exception throwing: */

	    reference front()
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }
	    const_reference front() const
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }
	    reference back()
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }
	    const_reference back() const
	    {
		throw std::out_of_range("Attempt to access element of empty vector");
	    }

	    /* Assignment, swapping and arithmetic operators are all trivial: */

	    vector<T,0>& assign(const value_type& c)
	    {
		return *this;
	    }
	    void swap(vector<T,0>& other){}
	    vector<T,0>& operator += (const value_type& c)
	    {
		return *this;
	    }
	    vector<T,0>& operator -= (const value_type& c)
	    {
		return *this;
	    }
	    vector<T,0>& operator *= (const value_type& c)
	    {
		return *this;
	    }
	    vector<T,0>& operator /= (const value_type& c)
	    {
		return *this;
	    }
	    vector<T,0>& operator += (const vector<T,0>& other)
	    {
		return *this;
	    }
	    vector<T,0>& operator -= (const vector<T,0>& other)
	    {
		return *this;
	    }
	    vector<T,0>& operator *= (const vector<T,0>& other)
	    {
		return *this;
	    }
	    vector<T,0>& operator /= (const vector<T,0>& other)
	    {
		return *this;
	    }
	    template<class rng_t>vector<T,0>& fill(random_number_stream<T,rng_t>& gen)
	    {
		return *this;
	    }
	    template<class rng_t>vector<T,0>& fill(random_number_stream<T,rng_t>& gen,const T& max)
	    {
		return *this;
	    }
	    template<class rng_t>vector<T,0>& fill(random_number_stream<T,rng_t>& gen,const T& min,const T& max)
	    {
		return *this;
	    }
    };
    template<class T>const std::size_t vector<T,0>::static_size;

    /// Element-wise equality operator.

    template<class T,std::size_t N>bool operator == (const vector<T,N>& p,const vector<T,N>& q)
    {
	return std::equal(p.begin(),p.end(),q.begin());
    }

    /// Element-wise unequality operator.
    
    template<class T,std::size_t N>bool operator != (const vector<T,N>& p,const vector<T,N>& q)
    {
	return !(p==q);
    }

    /// Lexicographical comparison operator.

    template<class T,std::size_t N>bool operator < (const vector<T,N>& p,const vector<T,N>& q)
    {
	return std::lexicographical_compare(p.begin(),p.end(),q.begin(),q.end());
    }

    /// Lexicographical comparison operator.

    template<class T,std::size_t N>bool operator > (const vector<T,N>& p,const vector<T,N>& q)
    {
	return q<p;
    }

    /// Lexicographical comparison operator.

    template<class T,std::size_t N>bool operator <= (const vector<T,N>& p,const vector<T,N>& q)
    {
	return !(q<p);
    }

    /// Lexicographical comparison operator.

    template<class T,std::size_t N>bool operator >= (const vector<T,N>& p,const vector<T,N>& q)
    {
	return !(q>p);
    }

    /// Returns p element-wise multiplied by c.

    template<class T,std::size_t N>vector<T,N> operator * (const T& c,const vector<T,N>& q)
    {
	vector<T,N>result(q);
	return result*=c;
    }

    /// Returns p element-wise multiplied by c.

    template<class T,std::size_t N>vector<T,N> operator * (const vector<T,N>& q,const T& c)
    {
	vector<T,N>result(q);
	return result*=c;
    }

    /// Returns p element-wise divided by c.

    template<class T,std::size_t N>vector<T,N> operator / (const vector<T,N>& p,const T& c)
    {
	vector<T,N>result(p);
	return result/=c;
    }

    /// Returns p element-wise multiplied by q.

    template<class T,std::size_t N>vector<T,N> operator * (const vector<T,N>& p,const vector<T,N>& q)
    {
	vector<T,N>result(q);
	return result*=q;
    }

    /// Returns the element-wise addition of p and q.

    template<class T,std::size_t N>vector<T,N> operator + (const vector<T,N>& p,const vector<T,N>& q)
    {
	vector<T,N>result(p);
	return result+=q;
    }

    /// Returns the element-wise subtraction of p and q.

    template<class T,std::size_t N>vector<T,N> operator - (const vector<T,N>& p,const vector<T,N>& q)
    {
	vector<T,N>result(p);
	return result-=q;
    }

    /// Returns the opposite of p.
    
    template<class T,std::size_t N>vector<T,N> operator - (const vector<T,N>& p)
    {
	vector<T,N>result(p);
	for(typename vector<T,N>::size_type i=0;i<N;++i)
	{
	    result[i]=-result[i];
	}
	return result;
    }

    /// Returns the (Euclidean) inner product of p and q

    template<class T,std::size_t N>T dot(const vector<T,N>& p,const vector<T,N>& q)
    {
	T s(0);
	for(typename vector<T,N>::size_type i=0;i<N;++i)
	{
	    s+=(p[i]*q[i]);
	}
	return s;
    }

    /// Prints the content of v to the output stream os.

    template<class T,std::size_t N>std::ostream& operator << (std::ostream& os,const vector<T,N>& v)
    {
	os<<"[";
	if(N>0)
	{
	    for(typename vector<T,N>::size_type i=0;i<v.size()-1;++i)
	    {
		os<<v[i]<<",";
	    }
	    os<<v[v.size()-1];
	}
	os<<"]";
	return os;
    }
}

#endif /*CAMGEN_VECTOR*/


