//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PART_SEQ_H_
#define CAMGEN_PART_SEQ_H_

#include <algorithm>
#include <Camgen/particle.h>
#include <Camgen/model.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Particle/jet sequence definition and declaration. When including a jet (or    *
 * some defined set of particles) the existing sequence is copied for all the    *
 * particles in the set, resulting in a 2D-array of particles, each row          *
 * corresponding to a subprocess. There is also a cleaning function defined that *
 * identifies identical subprocesses and erases the copies.                      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Particle vector comparison class, based on a lexicographical comparison
     * of the sorted sequences with element-wise comparison type comp_t: */

    template<class model_t,std::size_t N,class comp_t=flavour_comp< particle<model_t> > >class part_arr_comp
    {
	public:
	    bool operator () (const vector<const particle<model_t>*,N>& a1,const vector<const particle<model_t>*,N>& a2)
	    {
		/* Copy particle vectors in ascending order to the first and
		 * second temporaries: */

		std::partial_sort_copy(a1.begin(),a1.end(),first.begin(),first.end(),comp);
		std::partial_sort_copy(a2.begin(),a2.end(),second.begin(),second.end(),comp);
		
		/* Lexicographical comparison of the ordered sequences: */
		
		return std::lexicographical_compare(first.begin(),first.end(),second.begin(),second.end(),comp);
	    }

	private:

	    /* Temporary particle sequences to hold ordered vectors: */

	    vector<const particle<model_t>*,N>first;
	    vector<const particle<model_t>*,N>second;

	    /* Elementwise comparison object: */

	    comp_t comp;
    };

    /* Particle vector equality class, based on pairwise equality of the sorted
     * sequences with an element-wise comparison type comp_t: */

    template<class model_t,std::size_t N,class comp_t=flavour_comp< particle<model_t> > >class part_arr_eq
    {
	public:
	    bool operator () (const vector<const particle<model_t>*,N>& a1,const vector<const particle<model_t>*,N>& a2)
	    {
		/* Copy particle vectors in ascending order to the first and
		 * second temporaries: */

		std::partial_sort_copy(a1.begin(),a1.end(),first.begin(),first.end(),comp);
		std::partial_sort_copy(a2.begin(),a2.end(),second.begin(),second.end(),comp);
	
		return (first==second);	
	    }

	private:

	    /* Temporary particle sequences to hold ordered vectors: */

	    vector<const particle<model_t>*,N>first;
	    vector<const particle<model_t>*,N>second;

	    /* Elementwise comparison object: */

	    comp_t comp;
    };

    /* Particle sequence container, implemented as a stl (dynamic-sized) vector
     * of Camgen (static-sized) vectors: */

    template<class model_t,std::size_t N>class particle_sequence: public std::vector<vector<const particle<model_t>*,N> >
    {
	private:

	    /* Base type definition: */

	    typedef std::vector< vector<const particle<model_t>*,N> > base_type;
	public:

	    /* The usual container-compliant type definitions: */

	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::reference reference;
	    typedef typename base_type::const_reference const_reference;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;
	    typedef typename base_type::reverse_iterator reverse_iterator;
	    typedef typename base_type::const_reverse_iterator const_reverse_iterator;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::difference_type difference_type;
	    typedef typename base_type::pointer pointer;

	    /* Trivial constructor: */

	    particle_sequence():counter(0){}
	    
	    /* Assignment operator form particle name: */
	    
	    particle_sequence<model_t,N>& operator = (const std::string& name)
	    {
		/* Clear container: */

		this->clear();

		/* Obtain the model instance: */

		model<model_t>::get_instance();

		/* If the argument is a particle name, create a size 1 stl
		 * vector with element the vector with first entry the particle
		 * pointer and the rest NULL pointers: */

		if(model_wrapper<model_t>::contains_particle(name))
		{
		    value_type first;
		    first.assign(NULL);
		    first[0]=model_wrapper<model_t>::get_particle(name);
		    this->push_back(first);
		    ++counter;
		    return *this;
		}

		/* If the argument is jet name, create n vectors with first
		 * elements all addresses of particles in the jet, and the rest
		 * NULL pointers: */

		else if(model_wrapper<model_t>::contains_jet(name))
		{
		    value_type first;
		    first.assign(NULL);
		    const jet_definition<model_t>* jet=model_wrapper<model_t>::get_jet(name);
		    for(typename jet_definition<model_t>::const_iterator it=jet->begin();it != jet->end();++it)
		    {
			first[0]=*it;
			this->push_back(first);
		    }
		    ++counter;
		    return *this;
		}

		/* If the name is not recognised, complain...*/

		log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle "<<name<<" not found in model"<<endlog;
		return *this;
	    }

	    /* Overloaded comma operator, making the user interface easy: */

	    particle_sequence<model_t,N>& operator , (const std::string& name)
	    {

		/* If the string is recognised as a particle name, insert a
		 * particle address: */

		if(model_wrapper<model_t>::contains_particle(name))
		{
		    this->particle_sequence<model_t,N>::operator , (model_wrapper<model_t>::get_particle(name));
		    return *this;
		}

		/* If the string is recognised as a family name, insert a family
		 * address: */

		if(model_wrapper<model_t>::contains_family(name))
		{
		    this->particle_sequence<model_t,N>::operator , (model_wrapper<model_t>::get_family(name));
		    return *this;
		}

		/* If the name is not recognised, complain...*/

		log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle "<<name<<" not found in model"<<endlog;
		return *this;
	    }

	    /* Function determining whether the container is full: */

	    bool full() const
	    {
		return (counter==N);
	    }

	    /* Cleaning algorithm, deleting copies of the same final state: */

	    void clean()
	    {
		if(full())
		{
		    /* Comparison object creation: */

		    part_arr_comp<model_t,N>comp;

		    /* Sorting the container: */

		    std::sort(this->begin(),this->end(),comp);

		    /* Removal of copies: */

		    part_arr_eq<model_t,N>pred;
		    typename std::vector<vector<const particle<model_t>*,N> >::iterator it=std::unique(this->begin(),this->end(),pred);
		    this->resize(it-this->begin());
		}
	    }

	private:

	    /* Counter utility data type: */

	    size_type counter;

	    /* Insertion of particle pointer: */

	    particle_sequence<model_t,N>& operator , (const particle<model_t>* f)
	    {
		if(counter < N)
		{
		    for(size_type i=0;i<this->size();++i)
		    {
			(*this)[i][counter]=f;
		    }
		    ++counter;
		}
		return *this;
	    }

	    /* Insertion of jet definition pointer: */

	    particle_sequence<model_t,N>& operator , (const jet_definition<model_t>* f)
	    {
		if(counter<N)
		{
		    size_type s=this->size();
		    for(size_type i=0;i<s;++i)
		    {
			(*this)[i][counter]=f->front();
		    }
		    for(size_type l=1;l<f->size();++l)
		    {
			for(size_type i=0;i<s;++i)
			{
			    this->push_back((*this)[i]);
			    (this->back())[counter]=(*f)[l];
			}
		    }
		    ++counter;
		}
		return *this;
	    }
    };

    /* Outstream operator overload: */

    template<class model_t,std::size_t N>std::ostream& operator << (std::ostream& os,const particle_sequence<model_t,N>& p)
    {
	for(std::size_t i=0;i<p.size();++i)
	{
	    os<<"[";
	    for(std::size_t j=0;j<N-1;++j)
	    {
		os<<p[i][j]->get_name()<<",";
	    }
	    os<<p[i][N-1]->get_name()<<"]"<<std::endl;
	}
	return os;
    }
}

#endif /*CAMGEN_PART_SEQ_H_*/

