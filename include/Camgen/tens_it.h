//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_TENS_IT_H_
#define CAMGEN_TENS_IT_H_

#include <iterator>
#include <Camgen/tens.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Tensor iterator definitions. The tensor iterators are STL-compliant           *
 * random-access iterators with the extra functionality to move in arbitrary     *
 * directions and to compute the tensor indices corresponding to their position. *
 * There are four types: the normal, constant, reverse and constant reverse      *
 * tensor iterators.                                                             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Normal random-access tensor iterator definition: */

    template<class T>class tensor_iter
    {
	public:

	    /* The usual iterator category-compliant type definitions: */

	    typedef std::random_access_iterator_tag 		iterator_category;
	    typedef T 						value_type;
	    typedef T* 						pointer;
	    typedef T& 						reference;
	    typedef typename tensor<T>::difference_type		difference_type;
	    typedef typename tensor<T>::size_type 		size_type;

	    /* Declare the constant iterator a friend to allow conversion to
	     * const_iterator: */

	    friend class const_tensor_iter<T>;

	    /* Trivial constructor: */

	    tensor_iter():container(NULL){}
	    
	    /* Regular constructor: */
	    
	    tensor_iter(tensor<T>* t,const size_type n=0):container(t),offset(n){}
	    
	    /* Copy constructor: */
	    
	    tensor_iter(const tensor_iter<T>& other):container(other.container),offset(other.offset){}
	    
	    /* Dereferencing operator: */

	    reference operator * () const
	    {
		return container->data[offset];
	    }
	    
	    /* Dereferencing-to-member operator: */
	    
	    pointer operator -> () const
	    {
		return &(container->data[offset]);
	    }

	    /* Dereferencing with offset operator: */

	    reference operator [] (difference_type n) const
	    {
		return container->data[this->offset+n];
	    }

	    /* Assignment operator: */

	    tensor_iter<T>& operator = (const tensor_iter<T>& other)
	    {
		container=other.container;
		offset=other.offset;
		return *this;
	    }

	    /* Iterator position readout: */

	    size_type get_offset() const
	    {
		return offset;
	    }

	    /* Range checking utility: */

	    size_type range() const
	    {
		return container->size()-offset;
	    }

	    /* Pre-increment operator: */

	    tensor_iter<T>& operator ++ ()
	    {
		++offset;
		return *this;
	    }

	    /* Post-increment operator: */

	    tensor_iter<T> operator ++ (int)
	    {
		tensor_iter<T>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Multiple pre-increment operator: */

	    tensor_iter<T>& operator += (difference_type n)
	    {
		offset+=n;
		return *this;
	    }

	    /* Multiple copy-increment operator: */

	    tensor_iter<T> operator + (difference_type n) const
	    {
		return tensor_iter<T>(container,offset+n);
	    }

	    /* Pre-decrement operator: */

	    tensor_iter<T>& operator -- ()
	    {
		--offset;
		return *this;
	    }

	    /* Post-decrement operator: */

	    tensor_iter<T> operator -- (int)
	    {
		tensor_iter<T>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Multiple pre-decrement operator: */

	    tensor_iter<T>& operator -= (difference_type n)
	    {
		offset-=n;
		return *this;
	    }

	    /* Signed linear distance between iterators: */

	    difference_type operator - (tensor_iter<T> other) const
	    {
		return offset-other.offset;
	    }
	    difference_type operator - (const_tensor_iter<T> other) const
	    {
		return -(other-(*this));
	    }

	    /* Multiple copy-decrement operator: */

	    tensor_iter<T> operator - (difference_type n) const
	    {
		return tensor_iter<T>(container,offset-n);
	    }

	    /* Pre-increment in index direction dir: */

	    tensor_iter<T>& forward(size_type dir)
	    {
		offset+=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-increment in index direction dir: */

	    tensor_iter<T>& forward(size_type dir,difference_type n)
	    {
		offset+=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Pre-decrement in index direction dir: */

	    tensor_iter<T>& backward(size_type dir)
	    {
		offset-=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-decrement in index direction dir: */

	    tensor_iter<T>& backward(size_type dir,difference_type n)
	    {
		offset-=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Comparison operators: */

	    bool operator == (tensor_iter<T> it) const
	    {
		return (container==it.container and offset==it.offset);
	    }
	    bool operator != (tensor_iter<T> it) const
	    {
		return (container != it.container or offset != it.offset);
	    }
	    bool operator < (tensor_iter<T> it) const
	    {
		return (container==it.container and offset<it.offset);
	    }
	    bool operator > (tensor_iter<T> it) const
	    {
		return (container==it.container and offset>it.offset);
	    }
	    bool operator <= (tensor_iter<T> it) const
	    {
		return (container==it.container and offset<=it.offset);
	    }
	    bool operator >= (tensor_iter<T> it) const
	    {
		return (container==it.container and offset>=it.offset);
	    }

	    /* Send the iterator back to the beginning: */

	    tensor_iter<T>& reset()
	    {
		offset=0;
		return *this;
	    }

	    /* Compute the m-th index corresponding to the iterator position: */

	    size_type index(size_type m) const
	    {
		size_type n=offset;
		for(size_type i=0;i<m;++i)
		{
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		return n/(container->block_sizes[m]);
	    }

	    /* Index output facility: */

	    void print_indices() const
	    {
		std::cout<<"[  ";
		size_type n=offset;
		for(int i=0;i<container->tensor_rank;++i)
		{
		    std::cout<<n/(container->block_sizes[i])<<"  ";
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		std::cout<<"]";
	    }

	protected:

	    /* Tensor object holding the data: */

	    tensor<T>* container;
	    
	    /* Position of the iterator within the tensor: */
	    
	    size_type offset;
    };

    /* Reverse random-access tensor iterator definition: */

    template<class T>class reverse_tensor_iter
    {
	public:

	    /* The usual iterator category-compliant type definitions: */

	    typedef std::random_access_iterator_tag 		iterator_category;
	    typedef T 						value_type;
	    typedef T* 						pointer;
	    typedef T& 						reference;
	    typedef typename tensor<T>::difference_type		difference_type;
	    typedef typename tensor<T>::size_type 		size_type;

	    /* Declaring the corresponding const iterator friend to allow
	     * conversion: */

	    friend class const_reverse_tensor_iter<T>;

	    /* Trivial constructor: */

	    reverse_tensor_iter():container(NULL){}
	    
	    /* Regular constructor: */
	    
	    reverse_tensor_iter(tensor<T>* t,const size_type n):container(t),offset(t->size()-1-n){}
	    
	    /* Copy-constructor: */
	    
	    reverse_tensor_iter(const reverse_tensor_iter<T>& other):container(other.container),offset(other.offset){}
	    
	    /* Dereferencing operator: */
	    
	    reference operator * () const
	    {
		return container->data[offset];
	    }

	    /* Dereferencing-to-member operator: */

	    pointer operator -> () const
	    {
		return &(container->data[offset]);
	    }

	    /* Dereferencing with offset operator: */

	    reference operator [] (difference_type n) const
	    {
		return container->data[this->offset-n];
	    }

	    /* Assignment operator: */

	    reverse_tensor_iter<T>& operator = (const reverse_tensor_iter<T>& other)
	    {
		container=other.container;
		offset=other.offset;
		return *this;
	    }

	    /* Iterator position readout: */

	    size_type get_offset() const
	    {
		return container->size()-offset;
	    }

	    /* Range checking utility: */

	    size_type range() const
	    {
		return offset;
	    }

	    /* Pre-increment operator: */

	    reverse_tensor_iter<T>& operator ++ ()
	    {
		--offset;
		return *this;
	    }

	    /* Post-increment operator: */

	    reverse_tensor_iter<T> operator ++ (int)
	    {
		reverse_tensor_iter<T>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Multiple pre-increment operator: */

	    reverse_tensor_iter<T>& operator += (difference_type n)
	    {
		offset-=n;
		return *this;
	    }

	    /* Multiple copy-increment operator: */

	    reverse_tensor_iter<T> operator + (difference_type n) const
	    {
		return reverse_tensor_iter<T>(container,offset+n);
	    }

	    /* Pre-decrement operator: */

	    reverse_tensor_iter<T>& operator -- ()
	    {
		++offset;
		return *this;
	    }

	    /* Post-decrement operator: */

	    reverse_tensor_iter<T> operator -- (int)
	    {
		reverse_tensor_iter<T>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Multiple pre-decrement operator: */

	    reverse_tensor_iter<T>& operator -= (difference_type n)
	    {
		offset+=n;
		return *this;
	    }

	    /* Signed linear distance between iterators: */

	    difference_type operator - (reverse_tensor_iter<T> other) const
	    {
		return offset-other.offset;
	    }
	    difference_type operator - (const_reverse_tensor_iter<T> other) const
	    {
		return -(other-(*this));
	    }

	    /* Multiple copy-decrement operator: */

	    reverse_tensor_iter<T> operator - (difference_type n) const
	    {
		return reverse_tensor_iter<T>(container,offset-n);
	    }

	    /* Pre-increment in index direction dir: */

	    reverse_tensor_iter<T>& forward(size_type dir)
	    {
		offset-=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-increment in index direction dir: */

	    reverse_tensor_iter<T>& forward(size_type dir,difference_type n)
	    {
		offset-=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Pre-decrement in index direction dir: */

	    reverse_tensor_iter<T>& backward(size_type dir)
	    {
		offset+=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-decrement in index direction dir: */

	    reverse_tensor_iter<T>& backward(size_type dir,difference_type n)
	    {
		offset+=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Comparison operators: */

	    bool operator == (reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset==it.offset);
	    }
	    bool operator != (reverse_tensor_iter<T> it) const
	    {
		return (container != it.container or offset != it.offset);
	    }
	    bool operator < (reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset>it.offset);
	    }
	    bool operator > (reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset<it.offset);
	    }
	    bool operator <= (reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset>=it.offset);
	    }
	    bool operator >= (reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset<=it.offset);
	    }

	    /* Sending the iterator to the back of the tensor: */

	    reverse_tensor_iter<T>& reset()
	    {
		offset=container->size()-1;
		return *this;
	    }

	    /* m-th tensor index corresponding to the iterator position: */

	    size_type index(size_type m) const
	    {
		size_type n=offset;
		for(int i=0;i<m;++i)
		{
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		return n/(container->block_sizes[m]);
	    }

	    /* Index output facility: */

	    void print_indices() const
	    {
		std::cout<<"[  ";
		size_type n=offset;
		for(int i=0;i<container->tensor_rank;++i)
		{
		    std::cout<<n/(container->block_sizes[i])<<"  ";
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		std::cout<<"]";
	    }

	protected:

	    /* Tensor object holding the data: */

	    tensor<T>* container;
	    
	    /* Position of the iterator within the tensor: */
	    
	    size_type offset;
    };

    /* Constant iterator definition: */

    template<class T>class const_tensor_iter
    {
	public:

	    /* The usual iterator category-compliant type definitions: */

	    typedef std::random_access_iterator_tag 	iterator_category;
	    typedef T 					value_type;
	    typedef const T* 				pointer;
	    typedef const T& 				reference;
	    typedef typename tensor<T>::difference_type	difference_type;
	    typedef typename tensor<T>::size_type 		size_type;

	    /* Trivial constructor: */

	    const_tensor_iter():container(NULL){}
	    
	    /* Regular constructor: */
	    
	    const_tensor_iter(const tensor<T>* t,const size_type n):container(t),offset(n){}
	    
	    /* Copy constructor: */
	    
	    const_tensor_iter(const const_tensor_iter<T>& other):container(other.container),offset(other.offset){}
	    
	    /* Copy constructor from a non-const tensor tensor iterator: */
	    
	    const_tensor_iter(const tensor_iter<T>& other):container(other.container),offset(other.offset){}
	    
	    /* Dereferencing operator: */
	    
	    reference operator * () const
	    {
		return container->data[offset];
	    }

	    /* Dereference-to-member operator: */

	    pointer operator -> () const
	    {
		return &(container->data[offset]);
	    }

	    /* Dereferencing with offset operator: */

	    reference operator [] (difference_type n) const
	    {
		return container->data[this->offset+n];
	    }

	    /* Assignment operator: */
	    
	    const_tensor_iter<T>& operator = (const const_tensor_iter<T>& other)
	    {
		container=other.container;
		offset=other.offset;
		return *this;
	    }

	    /* Assignment from non-const iterator operator: */

	    const_tensor_iter<T>& operator = (const tensor_iter<T>& other)
	    {
		container=other.container;
		offset=other.offset;
		return *this;
	    }

	    /* Iterator position readout: */

	    size_type get_offset() const
	    {
		return offset;
	    }

	    /* Range checking utility: */

	    size_type range() const
	    {
		return container->size()-offset;
	    }

	    /* Pre-increment operator: */

	    const_tensor_iter<T>& operator ++ ()
	    {
		++offset;
		return *this;
	    }

	    /* Post-increment operator: */

	    const_tensor_iter<T> operator ++ (int)
	    {
		const_tensor_iter<T>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Multiple pre-increment operator: */

	    const_tensor_iter<T>& operator += (difference_type n)
	    {
		offset+=n;
		return *this;
	    }

	    /* Multiple copy-increment operator: */

	    const_tensor_iter<T> operator + (difference_type n) const
	    {
		return const_tensor_iter<T>(container,offset+n);
	    }

	    /* Pre-decrement operator: */

	    const_tensor_iter<T>& operator -- ()
	    {
		--offset;
		return *this;
	    }

	    /* Post-decrement operator: */

	    const_tensor_iter<T> operator -- (int)
	    {
		const_tensor_iter<T>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Multiple pre-decrement operator: */

	    const_tensor_iter<T>& operator -= (difference_type n)
	    {
		offset-=n;
		return *this;
	    }

	    /* Signed linear distance between iterators: */

	    difference_type operator - (const_tensor_iter<T> other) const
	    {
		return offset-other.offset;
	    }
	    difference_type operator - (tensor_iter<T> other) const
	    {
		return offset-other.offset;
	    }

	    /* Multiple copy-decrement operator: */

	    const_tensor_iter<T> operator - (difference_type n) const
	    {
		return const_tensor_iter<T>(container,offset-n);
	    }

	    /* Pre-increment in index-direction dir: */

	    const_tensor_iter<T>& forward(size_type dir)
	    {
		offset+=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-increment in index-direction dir: */

	    const_tensor_iter<T>& forward(size_type dir,difference_type n)
	    {
		offset+=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Pre-decrement in index direction dir: */

	    const_tensor_iter<T>& backward(size_type dir)
	    {
		offset+=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-increment in index direction dir: */

	    const_tensor_iter<T>& backward(size_type dir,difference_type n)
	    {
		offset+=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Comparison operators: */

	    bool operator == (const_tensor_iter<T> it) const
	    {
		return (container==it.container and offset==it.offset);
	    }
	    bool operator != (const_tensor_iter<T> it) const
	    {
		return (container != it.container or offset != it.offset);
	    }
	    bool operator < (const_tensor_iter<T> it) const
	    {
		return (container==it.container and offset<it.offset);
	    }
	    bool operator > (const_tensor_iter<T> it) const
	    {
		return (container==it.container and offset>it.offset);
	    }
	    bool operator <= (const_tensor_iter<T> it) const
	    {
		return (container==it.container and offset<=it.offset);
	    }
	    bool operator >= (const_tensor_iter<T> it) const
	    {
		return (container==it.container and offset>=it.offset);
	    }

	    /* Sending the iterator back to the first tensor entry: */

	    const_tensor_iter<T>& reset()
	    {
		offset=0;
		return *this;
	    }

	    /* m-th Tensor index corresponding to the iterator position: */

	    size_type index(size_type m) const
	    {
		size_type n=offset;
		for(int i=0;i<m;++i)
		{
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		return n/(container->block_sizes[m]);
	    }

	    /* Index output facility: */

	    void print_indices() const
	    {
		std::cout<<"[  ";
		size_type n=offset;
		for(int i=0;i<container->tensor_rank;++i)
		{
		    std::cout<<n/(container->block_sizes[i])<<"  ";
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		std::cout<<"]";
	    }

	protected:

	    /* Tensor object holding the data: */

	    const tensor<T>* container;
	    
	    /* Position of the iterator within the tensor: */
	    
	    size_type offset;
    };

    /* Constant reverse tensor iterator definition: */

    template<class T>class const_reverse_tensor_iter
    {
	public:

	    /* The usual iterator category-compliant type definitions: */

	    typedef std::random_access_iterator_tag 	iterator_category;
	    typedef T 					value_type;
	    typedef const T* 				pointer;
	    typedef const T& 				reference;
	    typedef typename tensor<T>::difference_type	difference_type;
	    typedef typename tensor<T>::size_type 		size_type;

	    /* Trivial constructor: */

	    const_reverse_tensor_iter():container(NULL){}
	    
	    /* Regular constructor: */

	    const_reverse_tensor_iter(const tensor<T>* t,const size_type n):container(t),offset(t->size()-1-n){}
	    
	    /* Copy constructor: */

	    const_reverse_tensor_iter(const const_reverse_tensor_iter<T>& other):container(other.container),offset(other.offset){}
	    
	    /* Copy constructor from non-const reverse iterator: */
	    
	    const_reverse_tensor_iter(const reverse_tensor_iter<T>& other):container(other.container),offset(other.offset){}
	    
	    /* Dereferencing operator: */
	    
	    reference operator * () const
	    {
		return container->data[offset];
	    }

	    /* Dereference-to-member operator: */

	    pointer operator -> () const
	    {
		return &(container->data[offset]);
	    }
	    
	    /* Dereference with offset operator: */
	    
	    reference operator [] (difference_type n) const
	    {
		return container->data[this->offset-n];
	    }

	    /* Assignment operator: */

	    const_reverse_tensor_iter<T>& operator = (const const_reverse_tensor_iter<T>& other)
	    {
		container=other.container;
		offset=other.offset;
		return *this;
	    }

	    /* Assignment from non-const reverse iterator: */

	    const_reverse_tensor_iter<T>& operator = (const reverse_tensor_iter<T>& other)
	    {
		container=other.container;
		offset=other.offset;
		return *this;
	    }

	    /* Iterator position readout: */

	    size_type get_offset() const
	    {
		return container->size()-offset;
	    }

	    /* Range checking utility: */

	    size_type range() const
	    {
		return container->size()-offset;
	    }

	    /* Pre-increment operator: */
	    
	    const_reverse_tensor_iter<T>& operator ++ ()
	    {
		--offset;
		return *this;
	    }

	    /* Post-increment operator: */

	    const_reverse_tensor_iter<T> operator ++ (int)
	    {
		const_reverse_tensor_iter<T>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Multiple pre-increment operator: */

	    const_reverse_tensor_iter<T>& operator += (difference_type n)
	    {
		offset-=n;
		return *this;
	    }

	    /* Multiple copy-increment operator: */

	    const_reverse_tensor_iter<T> operator + (difference_type n) const
	    {
		return const_reverse_tensor_iter<T>(container,offset+n);
	    }
	    
	    /* Pre-decrement operator: */
	    
	    const_reverse_tensor_iter<T>& operator -- ()
	    {
		++offset;
		return *this;
	    }

	    /* Post-decrement operator: */

	    const_reverse_tensor_iter<T> operator -- (int)
	    {
		const_reverse_tensor_iter<T>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Multiple post-decrement operator: */

	    const_reverse_tensor_iter<T>& operator -= (difference_type n)
	    {
		offset+=n;
		return *this;
	    }

	    /* Signed linear distance between iterators: */

	    difference_type operator - (const_reverse_tensor_iter<T> other) const
	    {
		return offset-other.offset;
	    }
	    difference_type operator - (reverse_tensor_iter<T> other) const
	    {
		return offset-other.offset;
	    }

	    /* Multiple copy-decrement operator: */

	    const_reverse_tensor_iter<T> operator - (difference_type n) const
	    {
		return const_reverse_tensor_iter<T>(container,offset-n);
	    }

	    /* Pre-increment in index-direction dir: */

	    const_reverse_tensor_iter<T>& forward(size_type dir)
	    {
		offset-=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-increment in index-direction dir: */

	    const_reverse_tensor_iter<T>& forward(size_type dir,difference_type n)
	    {
		offset-=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Pre-decrement in index-direction dir: */

	    const_reverse_tensor_iter<T>& backward(size_type dir)
	    {
		offset+=(container->block_sizes[dir]);
		return *this;
	    }

	    /* n-Fold pre-decrement in index-direction dir: */

	    const_reverse_tensor_iter<T>& backward(size_type dir,difference_type n)
	    {
		offset+=(n*(container->block_sizes[dir]));
		return *this;
	    }

	    /* Comparison operators: */

	    bool operator == (const_reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset==it.offset);
	    }
	    bool operator != (const_reverse_tensor_iter<T> it) const
	    {
		return (container != it.container or offset != it.offset);
	    }
	    bool operator < (const_reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset>it.offset);
	    }
	    bool operator > (const_reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset<it.offset);
	    }
	    bool operator <= (const_reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset>=it.offset);
	    }
	    bool operator >= (const_reverse_tensor_iter<T> it) const
	    {
		return (container==it.container and offset<=it.offset);
	    }

	    /* Sending the iterator to the back of the tensor: */

	    const_reverse_tensor_iter<T>& reset()
	    {
		offset=container->size()-1;
		return *this;
	    }

	    /* m-th tensor index corresponding to the iterator position: */

	    size_type index(size_type m) const
	    {
		size_type n=offset;
		for(int i=0;i<m;++i)
		{
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		return n/(container->block_sizes[m]);
	    }

	    /* Iterator output facility: */

	    void print_indices() const
	    {
		std::cout<<"[  ";
		size_type n=offset;
		for(int i=0;i<container->tensor_rank;++i)
		{
		    std::cout<<n/(container->block_sizes[i])<<"  ";
		    n-=((n/(container->block_sizes[i]))*(container->block_sizes[i]));
		}
		std::cout<<"]";
	    }

	protected:

	    /* Tensor object holding the data: */

	    tensor<T>* container;
	    
	    /* Position of the iterator within the tensor: */
	    
	    size_type offset;
    };
}

#endif /*CAMGEN_TENS_IT_H_*/

