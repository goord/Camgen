//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_TENS_H_
#define CAMGEN_TENS_H_

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstdarg>
#include <cstddef>
#include <stdexcept>
#include <Camgen/rn_strm.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Tensor class declaration and definition. The is the main data container type  *
 * in Camgen, since all internal polarisation vectors and spinors are in the    *
 * tensor format. This implementation has only one template parameter: the       *
 * numerical type. The implementation is essentially that of <valarray>: the     *
 * data is stored in a linear sequence (stl vector), but the access members      *
 * emulate a tensorial shape. Traversing the tensor in all directions is         *
 * optimised because the block sizes are stored in the class. However, this      *
 * makes the creation of tensors rather slow.                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Iterator declarations: */

    template<class T>class tensor_iter;
    template<class T>class const_tensor_iter;
    template<class T>class reverse_tensor_iter;
    template<class T>class const_reverse_tensor_iter;

    /* Tensor class declaration and definition: */

    template<class T>class tensor
    {
	
	public:

	    /* The usual container type definitions: */

	    typedef T 							value_type;
	    typedef T*							pointer;
	    typedef T&							reference;
	    typedef const T&						const_reference;
	    typedef typename std::vector<T>::size_type 			size_type;
	    typedef typename std::vector<T>::difference_type		difference_type;
	    typedef tensor_iter<T> 					iterator;
	    typedef const_tensor_iter<T> 				const_iterator;
	    typedef reverse_tensor_iter<T> 				reverse_iterator;
	    typedef const_reverse_tensor_iter<T> 			const_reverse_iterator;

	    /* The corresponding complex tensor class and the iterators are
	     * declared friends: */

	    friend class tensor< std::complex<T> >;
	    friend class tensor_iter<T>;
	    friend class const_tensor_iter<T>;
	    friend class reverse_tensor_iter<T>;
	    friend class const_reverse_tensor_iter<T>;

	    /* Trivial constructor, resulting in a size 1 (scalar) tensor with
	     * value 0: */
	    
	    tensor():data(1),tensor_rank(0)
	    {
	        data[0]=(value_type)0;
	    }

	    /* Constructor of a rank n tensor with index ranges r0,r1,...,rn: */

	    tensor(size_type n,const size_type r0,...)
	    {
		if(n != 0)
		{
		    va_list R;
		    size_type r=r0;
		    va_start(R,r0);
		    for(size_type i=0;i<n;++i)
		    {
			index_ranges.push_back(r);
			r=va_arg(R,size_type);
		    }
		    va_end(R);
		    clean_index_ranges();
		    make_blocks();
		    resize_data();
		}
	    }

	    /* Constructor of a multi-rank tensor with index ranges given by the
	     * integer vector R: */

	    tensor(const std::vector<size_type>& R):index_ranges(R)
	    {
		clean_index_ranges();
		make_blocks();
		resize_data();
	    }

	    /* Inserting n indices with range r at position pos: */

	    tensor<T>& insert_index(size_type r,size_type pos=0,size_type n=1)
	    {
		if(n!=0 and r!=0)
		{
		    if(tensor_rank==0)
		    {
			index_ranges.resize(n,r);
		    }
		    else
		    {
			index_ranges.insert(index_ranges.begin()+pos,n,r);
		    }
		    clean_index_ranges();
		    make_blocks();
		    resize_data();
		}
		return *this;
	    }

	    /* Adding n indices of 'type' G at position pos: */

	    template<class G>tensor<T>& add_index(size_type pos=0,size_type n=1)
	    {
		return add_index(G::index_range,pos,n);
	    }

	    /* Copy constructor: */

	    tensor(const tensor<T>& other):data(other.data),tensor_rank(other.tensor_rank),index_ranges(other.index_ranges),block_sizes(other.block_sizes){}
	    
	    /* Assignment operator: */
	    
	    tensor<T>& operator = (const tensor<T>& other)
	    {
		data=other.data;
		tensor_rank=other.tensor_rank;
		index_ranges=other.index_ranges;
		block_sizes=other.block_sizes;
		return *this;
	    }

	    /* Resize tensor shape (this may delete entries): */

	    tensor<T>& resize(size_type n,size_type r0,...)
	    {
		index_ranges.clear();
		if(n != 0)
		{
		    va_list R;
		    size_type r=r0;
		    va_start(R,r0);
		    for(size_type i=0;i<n;++i)
		    {
			index_ranges.push_back(r);
			r=va_arg(R,size_type);
		    }
		    va_end(R);
		}
		clean_index_ranges();
		make_blocks();
		resize_data();
		return *this;
	    }

	    /* Resize tensor shape with integer vector (this may delete
	     * entries): */

	    tensor<T>& resize(const std::vector<size_type>& R)
	    {
		index_ranges=R;
		clean_index_ranges();
		make_blocks();
		resize_data();
		return *this;
	    }

	    /* Resize tensor shape by the shape of another tensor: */

	    tensor<T>& resize(const tensor<T>& other)
	    {
		data.resize(other.data.size());
		block_sizes=other.block_sizes;
		index_ranges=other.index_ranges;
		tensor_rank=other.tensor_rank;
		return *this;
	    }

	    /* Resize tensor shape by the shape of a complex tensor: */

	    tensor<T>& resize(const tensor< std::complex<T> >& other)
	    {
		data.resize(other.data.size());
		block_sizes=other.block_sizes;
		index_ranges=other.index_ranges;
		tensor_rank=other.tensor_rank;
		return *this;
	    }

	    /* Clear all data, block sizes and index ranges included: */

	    tensor<T>& clear()
	    {
		data.resize(0);
		tensor_rank=0;
		block_sizes.clear();
		index_ranges.clear();
		return *this;
	    }

	    /* Reset all entries to zero: */

	    tensor<T>& reset()
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]=(T)0;
		}
		return *this;
	    }

	    /* Get the total tensor rank: */

	    size_type rank() const
	    {
		return tensor_rank;
	    }

	    /* Get the range of an index: */

	    size_type index_range(size_type i) const
	    {
		if(i<tensor_rank)
		{
			return index_ranges[i];
		}
		return 0;
	    }

	    /* Get (read-only) the vector of index ranges: */

	    const std::vector<size_type>& get_index_ranges() const
	    {
		return index_ranges;
	    }

	    /* Get the block size (i.e. the number of entries in a subtensor) of
	     * an index: */

	    size_type block_size(size_type i) const
	    {
		if(i<tensor_rank)
		{
			return block_sizes[i];
		}
		return 0;
	    }

	    /* Get (read-only) the vector of block sizes: */

	    const std::vector<size_type>& get_block_sizes() const
	    {
		return block_sizes;
	    }

	    /* Get the total number of components: */

	    size_type size() const
	    {
		return data.size();
	    }

	    /* Random element access by tensor notation: */

	    reference operator () (const size_type i0, ...)
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<tensor_rank;++j)
		{
		    if(i < index_ranges[j])
		    {
		    	n+=(i*block_sizes[j]);
		    	i=va_arg(I,size_type);
		    }
		    else
		    {
			throw std::out_of_range("tensor index out of range");
		    }
		}
		va_end(I);
		return data[n];
	    }

	    const_reference operator () (const size_type i0, ...) const
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<tensor_rank;++j)
		{
		    if(i < index_ranges[j])
		    {
		    	n+=(i*block_sizes[j]);
		    	i=va_arg(I,size_type);
		    }
		    else
		    {
			throw std::out_of_range("tensor index out of range");
		    }
		}
		va_end(I);
		return data[n];
	    }

	    /* Random element access in linear notation: */

	    reference operator [] (const size_type& n)
	    {
		return data[n];
	    }
	    const_reference operator [] (const size_type& n) const
	    {
		return data[n];
	    }

	    /* Random element access with exception throwing: */

	    reference at(const size_type& n)
	    {
		if(n<data.size())
		{
			return data[n];
		}
		else
		{
		    throw std::out_of_range("tensor element access out of range");
		}
	    }
	    const_reference at(const size_type& n) const
	    {
		if(n<data.size())
		{
			return data[n];
		}
		else
		{
		    throw std::out_of_range("tensor element access out of range");
		}
	    }

	    /* Multiplication by scalar: */

	    tensor<T>& operator *= (const T& c)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]*=c;
		}
		return *this;
	    }

	    tensor<T> operator * (const T& c) const
	    {
		tensor<T>result(*this);
		return result*=c;
	    }

	    /* Multiplication by complex scalar: */

	    tensor< std::complex<T> > operator * (const std::complex<T>& c) const
	    {
		tensor< std::complex<T> >result(*this);
		return result*c;
	    }

	    tensor<T> operator * (const std::pair<std::vector<T>,std::size_t>& v) const
	    {
		tensor<T>result(*this);
		return result*=v;
	    }

	    /* Division by a scalar: */

	    tensor<T>& operator /= (const T& c)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]/=c;
		}
		return *this;
	    }

	    tensor<T> operator / (const T& c) const
	    {
		tensor<T>result;
		return result/=c;
	    }

	    /* Division by a complex scalar: */

	    tensor< std::complex<T> > operator / (const std::complex<T>& c) const
	    {
		tensor<std::complex<T> >result(*this);
		return result/c;
	    }

	    /* Addition with another tensor: */

	    tensor<T>& operator += (const tensor<T>& other)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]+=other[i];
		}
		return *this;
	    }
	    tensor<T> operator + (const tensor<T>& other) const
	    {
		tensor<T>result(*this);
		return result+=other;
	    }
	    tensor< std::complex<T> > operator + (const tensor< std::complex<T> >& other) const
	    {
		tensor< std::complex<T> >result(*this);
		return result+=other;
	    }

	    /* Subtraction with another tensor: */

	    tensor<T>& operator -= (const tensor<T>& other)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]-=other[i];
		}
		return *this;
	    }
	    tensor<T> operator - (const tensor<T>& other) const
	    {
		tensor<T>result(*this);
		return result-=other;
	    }
	    tensor< std::complex<T> > operator - (const tensor< std::complex<T> >& other) const
	    {
		tensor< std::complex<T> >result(*this);
		return result-=other;
	    }

	    /* Comparison operators: */

	    bool operator == (const tensor<T>& other) const
	    {
		return (index_ranges==other.index_ranges and data==other.data);
	    }
	    bool operator != (const tensor<T>& other) const
	    {
		return !(*this==other);
	    }

	    /* Filling function with number x: */
	    
	    tensor<T>& fill(const T& x)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=x;
		}
		return *this;
	    }
	    
	    /* Filling function with random numbers provided by rng: */

	    template<class rng_t>tensor<T>& fill(random_number_stream<T,rng_t>& gen)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=gen();
		}
		return *this;
	    }
	    
	    /* Filling function with random numbers provided by rng within
	     * [0,max]: */

	    template<class rng_t>tensor<T>& fill(random_number_stream<T,rng_t>& gen,const value_type& max)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=gen(max);
		}
		return *this;
	    }
	    
	    /* Filling function with random numbers provided by rng within
	     * [min,max]: */

	    template<class rng_t>tensor<T>& fill(random_number_stream<T,rng_t>& gen,const value_type& min,const value_type& max)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=gen(min,max);
		}
		return *this;
	    }

	    /* Output method: */

	    std::ostream& print(std::ostream& os) const
	    {
		if(tensor_rank>0)
		{
		    os<<"rank:"<<tensor_rank<<std::endl;
		    for(size_type i=0;i<tensor_rank-1;++i)
		    {
			os<<index_ranges[i]<<" x ";
		    }
		    os<<index_ranges.back()<<std::endl;
		    if(tensor_rank==1)
		    {
			os<<std::setw(2);
			os<<"[";
			for(size_type i=0;i<index_ranges[0]-1;++i)
			{
			    os<<data[i]<<"  ";
			}
			os<<data.back()<<" ]"<<std::endl;
		    }
		    else if(tensor_rank==2)
		    {
			for(size_type i=0;i<index_ranges[0];++i)
			{
			    os<<std::setw(2);
			    os<<std::left<<"[";
			    for(size_type j=0;j<index_ranges[1]-1;++j)
			    {
				os<<data[i*block_sizes[0]+j]<<"  ";
			    }
			    os<<data[i*block_sizes[0]+index_ranges[1]-1]<<" ]"<<std::endl;
			}
		    }
		    else if(tensor_rank==3)
		    {
			for(size_type i=0;i<index_ranges[0];++i)
			{
			    for(size_type j=0;j<index_ranges[1];++j)
			    {
				os<<std::setw(2);
				os<<std::left<<"[";
				for(size_type k=0;k<index_ranges[2]-1;++k)
				{
				    os<<data[i*block_sizes[0]+j*block_sizes[1]+k]<<"  ";
				}
				os<<data[i*block_sizes[0]+j*block_sizes[1]+index_ranges[2]-1]<<" ]";
			    }
			    os<<std::endl;
			}
		    }
		    else
		    {
			os<<"no output implementation defined for tensors of rank>3"<<std::endl;
		    }
		}
		else
		{
		    os<<"["<<data[0]<<"]"<<std::endl;
		}
		return os;
	    }

	    /* Tensor iterators: */

	    iterator begin()
	    {
		return iterator(this,0);
	    }
	    iterator end()
	    {
		return iterator(this,data.size());
	    }
	    const_iterator begin() const
	    {
		return const_iterator(this,0);
	    }
	    const_iterator end() const
	    {
		return const_iterator(this,data.size());
	    }

	    /* Reverse tensor iterators: */

	    reverse_iterator rbegin()
	    {
		return reverse_iterator(this,0);
	    }
	    reverse_iterator rend()
	    {
		return iterator(this,data.size());
	    }
	    const_reverse_iterator rbegin() const
	    {
		return const_reverse_iterator(this,0);
	    }
	    const_reverse_iterator rend() const
	    {
		return const_reverse_iterator(this,data.size());
	    }
	private:
	    
	    /* Linearly stored tensor components: */

	    std::vector<value_type>data;
	    
	    /* Number of indices: */
	    
	    size_type tensor_rank;

	    /* Ranges of the indices: */

	    std::vector<size_type>index_ranges;
	    
	    /* Products of index ranges, i.e. sizes of subtensors: */
	    
	    std::vector<size_type>block_sizes;

	    /* Zero index range removal function: */

	    void clean_index_ranges()
	    {
		size_type zero(0);
		std::remove(index_ranges.begin(),index_ranges.end(),zero);
		tensor_rank=index_ranges.size();
	    }

	    /* Compute block sizes from index ranges: */

	    void make_blocks()
	    {
		if(tensor_rank != 0)
		{
		    block_sizes.resize(tensor_rank);
		    block_sizes.back()=1;
		    for(size_type i=0;i<tensor_rank-1;++i)
		    {
			block_sizes[tensor_rank-2-i]=block_sizes[tensor_rank-1-i]*index_ranges[tensor_rank-1-i];
		    }
		}
		else
		{
		    block_sizes.clear();
		}
	    }

	    /* Resize the data vector in correspondence with the given index
	     * ranges: */

	    void resize_data()
	    {
		if(tensor_rank != 0)
		{
		    data.resize(block_sizes[0]*index_ranges[0]);
		}
		else
		{
		    data.resize(1);
		}
	    }
    };

    /* Complex tensor specialisation: */

    template<class T>class tensor< std::complex<T> >
    {
	public:

	    /* The usual container type definitions: */

	    typedef std::complex<T> 						value_type;
	    typedef std::complex<T>*						pointer;
	    typedef std::complex<T>&						reference;
	    typedef const std::complex<T>&					const_reference;
	    typedef typename std::vector<std::complex<T> >::size_type 		size_type;
	    typedef typename std::vector<std::complex<T> >::difference_type	difference_type;
	    typedef tensor_iter< std::complex<T> > 				iterator;
	    typedef const_tensor_iter< std::complex<T> > 			const_iterator;
	    typedef reverse_tensor_iter< std::complex<T> > 			reverse_iterator;
	    typedef const_reverse_tensor_iter< std::complex<T> > 		const_reverse_iterator;

	    /* Corresponding real tensor class and iterators are declared
	     * friends: */

	    friend class tensor<T>;
	    friend class tensor_iter< std::complex<T> >;
	    friend class const_tensor_iter< std::complex<T> >;
	    friend class reverse_tensor_iter< std::complex<T> >;
	    friend class const_reverse_tensor_iter< std::complex<T> >;

	    /* Trivial constructor, resulting in a size 1 (scalar) tensor with
	     * value 0: */
	    
	    tensor():data(1),tensor_rank(0)
	    {
		data[0]=(value_type)0;
	    }

	    /* Constructor of a complex rank n tensor with index ranges r0,r1,...,rn: */

	    tensor(size_type n,const size_type r0,...)
	    {
		if(n != 0)
		{
		    va_list R;
		    size_type r=r0;
		    va_start(R,r0);
		    for(size_type i=0;i<n;++i)
		    {
			index_ranges.push_back(r);
			r=va_arg(R,size_type);
		    }
		    va_end(R);
		    clean_index_ranges();
		    make_blocks();
		    resize_data();
		}
	    }

	    /* Constructor of a multi-rank tensor with index ranges given by the
	     * integer vector R: */

	    tensor(const std::vector<size_type>& R):index_ranges(R)
	    {
		clean_index_ranges();
		make_blocks();
		resize_data();
	    }

	    /* Copy constructors: */

	    tensor(const tensor<std::complex<T> >& other):data(other.data),tensor_rank(other.tensor_rank),index_ranges(other.index_ranges),block_sizes(other.block_sizes){}
	    tensor(const tensor<T>& tens):tensor_rank(tens.tensor_rank),index_ranges(tens.index_ranges),block_sizes(tens.block_sizes)
	    {
		data.resize(tens.data.size());
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]=value_type(tens.data[i],0);
		}
	    }

	    /* Assignment operators: */

	    tensor<std::complex<T> >& operator = (const tensor<std::complex<T> >& other)
	    {
		data=other.data;
		tensor_rank=other.tensor_rank;
		index_ranges=other.index_ranges;
		block_sizes=other.block_sizes;
		return *this;
	    }
	    tensor<std::complex<T> >& operator = (const tensor<T>& other)
	    {
		tensor_rank=other.tensor_rank;
		index_ranges=other.index_ranges;
		block_sizes=other.block_sizes;
		data.resize(other.data.size());
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]=std::complex<T>(other.data[i],0);
		}
		return *this;
	    }

	    /* Inserting n indices with range r at position pos: */

	    tensor< std::complex<T> >& insert_index(size_type r,size_type pos=0,size_type n=1)
	    {
		if(n!=0 and r!=0)
		{
		    if(tensor_rank==0)
		    {
			index_ranges.resize(n,r);
		    }
		    else
		    {
			index_ranges.insert(index_ranges.begin()+pos,n,r);
		    }
		    clean_index_ranges();
		    make_blocks();
		    resize_data();
		}
		return *this;
	    }

	    /* Adding n indices of 'type' G at position pos: */

	    template<class G>tensor< std::complex<T> >& add_index(size_type pos=0,size_type n=1)
	    {
		return add_index(G::index_range,pos,n);
	    }

	    /* Resize tensor shape (this may delete entries): */

	    tensor<std::complex<T> >& resize(size_type n,size_type r0,...)
	    {
		index_ranges.clear();
		if(n != 0)
		{
		    va_list R;
		    size_type r=r0;
		    va_start(R,r0);
		    for(size_type i=0;i<n;++i)
		    {
			index_ranges.push_back(r);
			r=va_arg(R,size_type);
		    }
		    va_end(R);
		}
		clean_index_ranges();
		make_blocks();
		resize_data();
		return *this;
	    }

	    /* Resize tensor shape with integer vector (this may delete
	     * entries): */

	    tensor<std::complex<T> >& resize(const std::vector<size_type>& R)
	    {
		index_ranges=R;
		clean_index_ranges();
		make_blocks();
		resize_data();
		return *this;
	    }

	    /* Resize tensor shape by the shape of another complex tensor: */

	    tensor<std::complex<T> >& resize(const tensor<std::complex<T> >& other)
	    {
		data.resize(other.data.size());
		block_sizes=other.block_sizes;
		index_ranges=other.index_ranges;
		tensor_rank=other.tensor_rank;
		return *this;
	    }

	    /* Resize tensor shape by the shape of another real tensor: */

	    tensor<std::complex<T> >& resize(const tensor<T>& other)
	    {
		data.resize(other.data.size());
		block_sizes=other.block_sizes;
		index_ranges=other.index_ranges;
		tensor_rank=other.tensor_rank;
		return *this;
	    }

	    /* Clear all data, block sizes and index ranges included: */

	    tensor<std::complex<T> >& clear()
	    {
		data.resize(0);
		tensor_rank=0;
		block_sizes.clear();
		index_ranges.clear();
		return *this;
	    }

	    /* Reset all entries to zero: */

	    tensor<std::complex<T> >& reset()
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]=std::complex<T>(0,0);
		}
		return *this;
	    }

	    /* Get the total tensor rank: */

	    size_type rank() const
	    {
		return tensor_rank;
	    }

	    /* Get the range of an index: */

	    size_type index_range(size_type i) const
	    {
		return index_ranges[i];
	    }

	    /* Get (read-only) the vector of index ranges: */

	    const std::vector<size_type>& get_index_ranges() const
	    {
		return index_ranges;
	    }

	    /* Get the block size (i.e. the number of entries in a subtensor) of
	     * an index: */

	    size_type block_size(size_type i) const
	    {
		return block_sizes[i];
	    }

	    /* Get (read-only) the vector of block sizes: */

	    const std::vector<size_type>& get_block_sizes() const
	    {
		return block_sizes;
	    }

	    /* Get the total number of components: */

	    size_type size() const
	    {
		return data.size();
	    }
	    
	    /* Random element access by tensor notation: */

	    reference operator () (const size_type i0, ...)
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<tensor_rank;++j)
		{
		    if(i < index_ranges[j])
		    {
		    	n+=(i*block_sizes[j]);
		    	i=va_arg(I,size_type);
		    }
		    else
		    {
			throw std::out_of_range("tensor index out of range");
		    }
		}
		va_end(I);
		return data[n];
	    }

	    const_reference operator () (const size_type i0, ...) const
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<tensor_rank;++j)
		{
		    if(i < index_ranges[j])
		    {
		    	n+=(i*block_sizes[j]);
		    	i=va_arg(I,size_type);
		    }
		    else
		    {
			throw std::out_of_range("tensor index out of range");
		    }
		}
		va_end(I);
		return data[n];
	    }

	    /* Random element access in linear notation: */

	    reference operator [] (const size_type& n)
	    {
		return data[n];
	    }
	    const_reference operator [] (const size_type& n) const
	    {
		return data[n];
	    }

	    /* Random element access with exception throwing: */

	    reference at(const size_type& n)
	    {
		if(n<data.size())
		{
			return data[n];
		}
		else
		{
		    throw std::out_of_range("tensor element access out of range");
		}
	    }
	    const_reference at(const size_type& n) const
	    {
		if(n<data.size())
		{
			return data[n];
		}
		else
		{
		    throw std::out_of_range("tensor element access out of range");
		}
	    }

	    /* Multiplication by complex scalar: */

	    tensor<std::complex<T> >& operator *= (const std::complex<T> & c)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]*=c;
		}
		return *this;
	    }
	    tensor<std::complex<T> > operator * (const std::complex<T>& c) const
	    {
		tensor<std::complex<T> >result(*this);
		return result*=c;
	    }

	    /* Multiplication by real scalar: */

	    tensor<std::complex<T> >& operator *= (const T& c)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]*=c;
		}
		return *this;
	    }			
	    tensor<std::complex<T> > operator * (const T& c) const
	    {
		tensor<std::complex<T> >result(*this);
		return result*=c;
	    }

	    /* Division by scalar: */

	    tensor<std::complex<T> >& operator /= (const std::complex<T>& c)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]/=c;
		}
		return *this;
	    }
	    tensor<std::complex<T> > operator / (const std::complex<T>& c) const
	    {
		tensor<std::complex<T> >result(*this);
		return result/=c;
	    }

	    /* Division by real scalar: */

	    tensor<std::complex<T> >& operator /= (const T& c)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]/=c;
		}
		return *this;
	    }
	    tensor<std::complex<T> > operator / (const T& c) const
	    {
		tensor<std::complex<T> >result(*this);
		return result/=c;
	    }

	    /* Addition with another tensor: */

	    tensor<std::complex<T> >& operator += (const tensor<std::complex<T> >& other)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]+=other[i];
		}
		return *this;
	    }
	    tensor<std::complex<T> >& operator += (const tensor<T>& other)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]+=std::complex<T>(other[i],0);
		}
		return *this;
	    }
	    tensor<std::complex<T> > operator + (const tensor<std::complex<T> >& other) const
	    {
		tensor<std::complex<T> >result(*this);
		return result+=other;
	    }
	    tensor<std::complex<T> > operator + (const tensor<T>& other) const
	    {
		tensor<std::complex<T> >result(*this);
		return result+=other;
	    }

	    /* Subtraction with another tensor: */

	    tensor<std::complex<T> >& operator -= (const tensor<std::complex<T> >& other)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]-=other[i];
		}
		return *this;
	    }
	    tensor<std::complex<T> >& operator -= (const tensor<T>& other)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i]-=std::complex<T>(other[i],0);
		}
		return *this;
	    }
	    tensor<std::complex<T> > operator - (const tensor<std::complex<T> >& other) const
	    {
		tensor<std::complex<T> >result(*this);
		return result-=other;
	    }
	    tensor<std::complex<T> > operator - (const tensor<T>& other) const
	    {
		tensor<std::complex<T> >result(*this);
		return result+=other;
	    }

	    /* Comparison operators: */

	    bool operator == (const tensor<std::complex<T> >& other) const
	    {
		return (index_ranges==other.index_ranges and data==other.data);
	    }
	    bool operator != (const tensor<std::complex<T> >& other) const
	    {
		return !(*this==other);
	    }

	    /* Filling function with number x: */
	    
	    tensor< std::complex<T> >& fill(const std::complex<T>& x)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=x;
		}
		return *this;
	    }
	    
	    /* Filling function with random numbers provided by rng: */

	    template<class rng_t>tensor< std::complex<T> >& fill(random_number_stream<T,rng_t>& gen)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=std::complex<T>(gen(),gen());
		}
		return *this;
	    }
	    
	    /* Filling function with random numbers provided by rng within
	     * [0,max]: */

	    template<class rng_t>tensor< std::complex<T> >& fill(random_number_stream<T,rng_t>& gen,const value_type& max)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=std::complex<T>(gen(max.real()),gen(max.imag()));
		}
		return *this;
	    }
	    
	    /* Filling function with random numbers provided by rng within
	     * [min,max]: */

	    template<class rng_t>tensor< std::complex<T> >& fill(random_number_stream<T,rng_t>& gen,const value_type& min,const value_type& max)
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    this->data[i]=std::complex<T>(gen(min.real(),max.real()),gen(min.imag(),max.imag()));
		}
		return *this;
	    }

	    /* Output method: */

	    std::ostream& print(std::ostream& os) const
	    {
		if(tensor_rank>0)
		{
		    os<<"rank:"<<tensor_rank<<std::endl;
		    for(size_type i=0;i<tensor_rank-1;++i)
		    {
			os<<index_ranges[i]<<" x ";
		    }
		    os<<index_ranges.back()<<std::endl;
		    if(tensor_rank==1)
		    {
			os<<std::setw(2);
			os<<"[";
			for(size_type i=0;i<index_ranges[0]-1;++i)
			{
			    os<<data[i]<<"  ";
			}
			os<<data.back()<<" ]"<<std::endl;
		    }
		    else if(tensor_rank==2)
		    {
			for(size_type i=0;i<index_ranges[0];++i)
			{
			    os<<std::setw(2);
			    os<<std::left<<"[";
			    for(size_type j=0;j<index_ranges[1]-1;++j)
			    {
				os<<data[i*block_sizes[0]+j]<<"  ";
			    }
			    os<<data[i*block_sizes[0]+index_ranges[1]-1]<<" ]"<<std::endl;
			}
		    }
		    else if(tensor_rank==3)
		    {
			for(size_type i=0;i<index_ranges[0];++i)
			{
			    for(size_type j=0;j<index_ranges[1];++j)
			    {
				os<<std::setw(2);
				os<<std::left<<"[";
				for(size_type k=0;k<index_ranges[2]-1;++k)
				{
				    os<<data[i*block_sizes[0]+j*block_sizes[1]+k]<<"  ";
				}
				os<<data[i*block_sizes[0]+j*block_sizes[1]+index_ranges[2]-1]<<" ]";
			    }
			    os<<std::endl;
			}
		    }
		    else
		    {
			os<<"no output implementation defined for tensors of rank>3"<<std::endl;
		    }
		}
		else
		{
		    os<<"["<<data[0]<<"]"<<std::endl;
		}
		return os;
	    }

	    /* Iterators: */

	    iterator begin()
	    {
		return iterator(this,0);
	    }
	    iterator end()
	    {
		return iterator(this,data.size());
	    }
	    const_iterator begin() const
	    {
		return const_iterator(this,0);
	    }
	    const_iterator end() const
	    {
		return const_iterator(this,data.size());
	    }

	    /* Reverse iterators: */

	    reverse_iterator rbegin()
	    {
		return reverse_iterator(this,0);
	    }
	    iterator rend()
	    {
		return reverse_iterator(this,data.size());
	    }
	    const_reverse_iterator rbegin() const
	    {
		return const_reverse_iterator(this,0);
	    }
	    const_reverse_iterator rend() const
	    {
		return const_reverse_iterator(this,data.size());
	    }
	private:
	    
	    /* Linearly stored tensor components: */

	    std::vector<value_type>data;
	    
	    /* Number of indices: */
	    
	    size_type tensor_rank;

	    /* Ranges of the indices: */

	    std::vector<size_type>index_ranges;
	    
	    /* Products of index ranges, i.e. sizes of subtensors: */
	    
	    std::vector<size_type>block_sizes;

	    /* Zero index range removal function: */

	    void clean_index_ranges()
	    {
		size_type zero(0);
		std::remove(index_ranges.begin(),index_ranges.end(),zero);
		tensor_rank=index_ranges.size();
	    }

	    /* Compute block sizes from index ranges: */

	    void make_blocks()
	    {
		if(tensor_rank != 0)
		{
		    block_sizes.resize(tensor_rank);
		    block_sizes.back()=1;
		    for(size_type i=0;i<tensor_rank-1;++i)
		    {
			block_sizes[tensor_rank-2-i]=block_sizes[tensor_rank-1-i]*index_ranges[tensor_rank-1-i];
		    }
		}
		else
		{
		    block_sizes.clear();
		}
	    }

	    /* Resize the data vector in correspondence with the given index
	     * ranges: */

	    void resize_data()
	    {
		if(tensor_rank != 0)
		{
		    data.resize(block_sizes[0]*index_ranges[0]);
		}
		else
		{
		    data.resize(1);
		}
	    }
    };

    /* Arithmetic operators with tensors: */

    template<class T>Camgen::tensor<T> operator - (const Camgen::tensor<T>& tens)
    {
	return tens*=(-(T)1);
    }

    template<class T>Camgen::tensor<T> operator * (const T& c,const Camgen::tensor<T>& tens)
    {
	return tens*c;
    }

    template<class T>Camgen::tensor<std::complex<T> >operator * (const std::complex<T>& c,const Camgen::tensor<T>& tens)
    {
	return tens*c;
    }

    template<class T>Camgen::tensor<std::complex<T> >operator * (const T& c,const Camgen::tensor<std::complex<T> >& tens)
    {
	return tens*c;
    }

    template<class T>Camgen::tensor<std::complex<T> >operator * (const std::complex<T>& c,const Camgen::tensor<std::complex<T> >& tens)
    {
	return tens*c;
    }

    /* Output operator: */

    template<class T>std::ostream& operator << (std::ostream& os, const Camgen::tensor<T>& tens)
    {
	return tens.print(os);
    }
}

namespace std
{
    template<class T>Camgen::tensor< complex<T> >conj(const Camgen::tensor< complex<T> >& t)
    {
	Camgen::tensor< complex<T> >result;
	result.resize(t);
	for(typename Camgen::tensor< complex<T> >::size_type i=0;i<result.size();++i)
	{
	    result[i]=std::conj(t[i]);
	}
	return result;
    }
}

#endif /*CAMGEN_TENS_H_*/

