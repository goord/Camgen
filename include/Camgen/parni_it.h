//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PARNI_IT_H_
#define CAMGEN_PARNI_IT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Iterator definitions to traverse the parni binary tree. *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iterator>
#include <Camgen/parni_bin.h>

namespace Camgen
{
    /* Bin iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class parni_bin_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    parni_bin_iterator<value_t,D,rng_t,key_t>& operator = (const parni_bin_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    parni_bin_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_bin();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    parni_bin_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    parni_bin_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_bin();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    parni_bin_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    parni_bin_iterator(parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    parni_bin_iterator(const parni_bin_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    parni_bin<value_t,D,rng_t,key_t>* bin;
    };
    
    /* Const bin iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class const_parni_bin_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef const parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef const parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef const parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    const_parni_bin_iterator<value_t,D,rng_t,key_t>& operator = (const const_parni_bin_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    const_parni_bin_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_bin();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    const_parni_bin_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		const_parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    const_parni_bin_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_bin();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    const_parni_bin_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		const_parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const const_parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const const_parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    const_parni_bin_iterator(const parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    const_parni_bin_iterator(const const_parni_bin_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    const parni_bin<value_t,D,rng_t,key_t>* bin;
    };

    /* Leaf iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class parni_leaf_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    parni_leaf_iterator<value_t,D,rng_t,key_t>& operator = (const parni_leaf_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    parni_leaf_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_leaf();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    parni_leaf_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    parni_leaf_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_leaf();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    parni_leaf_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    parni_leaf_iterator(parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    parni_leaf_iterator(const parni_leaf_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    parni_bin<value_t,D,rng_t,key_t>* bin;
    };

    /* const Leaf iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class const_parni_leaf_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef const parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef const parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef const parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef const typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    const_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator = (const const_parni_leaf_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    const_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_leaf();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    const_parni_leaf_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		const_parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    const_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_leaf();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    const_parni_leaf_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		const_parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const const_parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const const_parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    const_parni_leaf_iterator(const parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    const_parni_leaf_iterator(const const_parni_leaf_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    const parni_bin<value_t,D,rng_t,key_t>* bin;
    };

    /* Reverse bin iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class reverse_parni_bin_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& operator = (const reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_bin();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    reverse_parni_bin_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		reverse_parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_bin();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    reverse_parni_bin_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		reverse_parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    reverse_parni_bin_iterator(parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    reverse_parni_bin_iterator(const reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    parni_bin<value_t,D,rng_t,key_t>* bin;
    };
    
    /* Const bin iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class const_reverse_parni_bin_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef const parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef const parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef const parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& operator = (const const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_bin();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_bin();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    const_reverse_parni_bin_iterator(const parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    const_reverse_parni_bin_iterator(const const_reverse_parni_bin_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    const parni_bin<value_t,D,rng_t,key_t>* bin;
    };

    /* Leaf iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class reverse_parni_leaf_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator = (const reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_leaf();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    reverse_parni_leaf_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_leaf();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    reverse_parni_leaf_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    reverse_parni_leaf_iterator(parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    reverse_parni_leaf_iterator(const reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    parni_bin<value_t,D,rng_t,key_t>* bin;
    };

    /* const Leaf iterator class definition: */

    template<class value_t,std::size_t D,class rng_t,class key_t>class const_reverse_parni_leaf_iterator
    {
	friend class parni<value_t,D,rng_t,key_t>;
	friend class parni_bin<value_t,D,rng_t,key_t>;

	public:

	    /* Type definitions: */

	    typedef std::forward_iterator_tag iterator_category;
	    typedef const parni_bin<value_t,D,rng_t,key_t>* value_type;
	    typedef const parni_bin<value_t,D,rng_t,key_t>** pointer;
	    typedef const parni_bin<value_t,D,rng_t,key_t>*& reference;
	    typedef const typename parni_bin<value_t,D,rng_t,key_t>::size_type size_type;

	    /* Assignment operator: */

	    const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator = (const const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& other)
	    {
		bin=other.bin;
		return *this;
	    }

	    /* Dereferencing operator: */

	    value_type operator * () const
	    {
		return bin;
	    }

	    /* Pre-increment operator: */

	    const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator ++ ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->previous_leaf();
		}
		return *this;
	    }

	    /* Post-increment operator: */

	    const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t> operator ++ (int)
	    {
		const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		++(*this);
		return clone;
	    }

	    /* Pre-decrement operator: */

	    const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& operator -- ()
	    {
		if(bin!=NULL)
		{
		    bin=bin->next_leaf();
		}
		return *this;
	    }

	    /* Post-decrement operator: */

	    const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t> operator -- (int)
	    {
		const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>clone(*this);
		--(*this);
		return clone;
	    }

	    /* Comparison operator: */

	    bool operator == (const const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin==other.bin);
	    }
	    bool operator != (const const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& other) const
	    {
		return (bin!=other.bin);
	    }

	private:

	    /* Constructor: */

	    const_reverse_parni_leaf_iterator(const parni_bin<value_t,D,rng_t,key_t>* bin_=NULL):bin(bin_){}

	    /* Copy constructor: */

	    const_reverse_parni_leaf_iterator(const const_reverse_parni_leaf_iterator<value_t,D,rng_t,key_t>& it):bin(it.bin){}

	    /* Container element pointer: */

	    const parni_bin<value_t,D,rng_t,key_t>* bin;
    };
}

#endif /*CAMGEN_PARNI_IT_H_*/

