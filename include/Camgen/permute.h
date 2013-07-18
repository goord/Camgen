//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PERMUTE_H_
#define CAMGEN_PERMUTE_H_

#include <utility>
#include <functional>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Permutation class declaration and definition. Computes and stores the       *
 * permutation taking a given vector into a given other one in terms of swaps. *
 * Then this permutation can be used to reshuffle any vector.                  *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 

namespace Camgen
{
    /* Definition with generic comparison type: */

    template<class container_type,class compare=std::less<typename container_type::value_type> >class permutation: public std::vector<std::pair<typename container_type::size_type,typename container_type::size_type> >
    {
	public:

	    /* Some useful typ definitions: */

	    typedef typename container_type::size_type size_type;
	    typedef typename container_type::value_type value_type;
	    typedef typename container_type::iterator iterator;
	    typedef typename container_type::const_iterator const_iterator;
	    typedef std::vector<std::pair<size_type,size_type> > base_type;

	private:

	    typedef std::pair<size_type,size_type> swap_type;

	public:

	    /* The trivial constructor: */

	    permutation():container_size(0){}

	    /* A constructor taking one vector as argument; the permutation shall be the
	     * one taking v into the sorted v: */

	    permutation(const container_type& v):container_size(v.size())
	    {

		/* Create 2 copies of v and a comparison object: */

		container_type test=v;
		container_type sorted_test=v;
		compare comp;
		
		/* Sort one copy: */
		
		std::sort(sorted_test.begin(),sorted_test.end(),comp);
		
		/* Start constructing the swap pairs: */
		
		const_iterator it=sorted_test.begin();
		size_type i=0;
		do
		{
		    /* If the entry in the sorted vector equals the one in the unsorted
		     * range, move on: */

		    if(*it==test[i])
		    {
			++it;
			++i;
		    }
		    else
		    {
			/* Look for the element of the sorted range in the unsorted range
			 * (stop if the sorted range contains duplicates):
			 * */

			size_type j=i;
			do
			{
			    ++j;
			}
			while(test[j] != *it or sorted_test[i] == sorted_test[j]);

			/* Add the found swap to the data of the class: */

			std::pair<size_type,size_type>p(i,j);
			this->push_back(p);

			/* Perform the swap to the unsorted range: */

			std::swap(test[i],test[j]);
			++it;
			++i;
		    }
		}
		while(it != sorted_test.end());
		CAMGEN_ERROR_IF((test != sorted_test),"permutation construction failed");
	    }

	    /* A constructor taking two vectors as arguments; the permutation shall be the
	     * one taking v into w: */

	    permutation(const container_type& v,const container_type& w):container_size(v.size())
	    {

		/* Creating a copy of v: */

		container_type test(v);

		/* Start constructing the swap pairs: */

		const_iterator it=w.begin();

		size_type i=0;
		
		do
		{
		    /* If the entry in w equals the one in v, move on: */

		    if(*it==test[i])
		    {
			++it;
			++i;
		    }
		    else
		    {
			/* Look for the element of w in v (stop if w contains duplicates):*/

			size_type j=i;
			do
			{
			    ++j;
			}
			while(test[j] != *it or w[i]==w[j]);

			/* Add the found swap to the class data: */

			this->push_back(std::pair<size_type,size_type>(i,j));

			/* Perform the swap on v's copy: */

			std::swap(test[i],test[j]);
			++it;
			++i;
		    }
		}
		while(it != w.end());
		CAMGEN_ERROR_IF((test!=w),"permutation construction failed");
	    }

	    /* Inverting the permutation = reversing the order of the swaps: */

	    permutation<container_type,compare>& invert()
	    {
		std::reverse(this->begin(),this->end());
		return *this;
	    }

	    /* Parity of the permutation = number of swaps modulo 2*/
	    
	    int parity() const
	    {
		if(this->size()%2==0)
		{
		    return 1;
		}
		else
		{
		    return -1;
		}
	    }

	    /* Application operator: apply swaps to any random-access container: */

	    template<class random_access_container>random_access_container& operator () (random_access_container& w) const
	    {
		CAMGEN_ERROR_IF((w.size() < container_size),"permutation swaps may be out of range");
		for(size_type i=0;i<this->size();++i)
		{
		    std::swap(w[(*this)[i].first],w[(*this)[i].second]);
		}
		return w;
	    }
	    template<class random_access_container>random_access_container& reverse_permute(random_access_container& w) const
	    {
		CAMGEN_ERROR_IF((w.size() < container_size),"permutation swaps may be out of range");
		for(typename base_type::const_reverse_iterator it=this->rbegin();it!=this->rend();++it)
		{
		    std::swap(w[it->first],w[it->second]);
		}
		return w;
	    }

	private:

	    /* (Minimum) size of vectors that may be permuted: */

	    size_type container_size;
    };
}

/* Output operator: */

template<class container_type,class compare>std::ostream& operator << (std::ostream& os,const Camgen::permutation<container_type,compare>& s)
{
    for(int i=0;i<s.size();++i)
    {
	os<<"("<<s[i].first<<" "<<s[i].second<<")";
    }
    return os;
}


#endif /*PERMUTE_H_*/

