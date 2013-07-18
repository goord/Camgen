//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FUSION_CLASS_H_
#define CAMGEN_FUSION_CLASS_H_

#include <map>
#include <Camgen/type_holders.h>
#include <Camgen/permute.h>
#include <Camgen/flav_comp.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the fusion class template. Given an n-vertex and a list of      *
 * particle addresses of length n-1, and a produced particle pointer, it         *
 * computes the permutation to bring the legs of the vertex minus the produced   *
 * particles to the given sequence. The fusion class objects are stored in the   *
 * model wrapper as vertices are included, and represent all the ways a particle *
 * of a given type may be produced. This listings are the searched through as    *
 * the fusion tree is constructed.                                               *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class fusion_class
    {
	public:

	    /* Particle and vertex type definitions: */

	    typedef typename get_basic_types<model_t>::particle_type particle_type;
	    typedef typename get_basic_types<model_t>::vertex_type vertex_type;
	    typedef std::multimap<std::vector<const particle_type*>,fusion_class<model_t>,flavourvec_comp<particle_type> > fusion_map_type;
	    typedef typename std::vector<const particle_type*>::iterator particle_iterator;
	    typedef typename std::vector<const particle_type*>::size_type size_type;

	    /* Static function inserting fusion rules originating from vertex
	     * 'vert' into the fusion_map, using the private regular fusion
	     * class constructor: */

	    static void insert_fusion_rules(const vertex_type* vert,fusion_map_type& fusion_map)
	    {
		/* Copy the sorted legs of the vertex and omit all identical
		 * legs: */

		std::vector<const particle_type*>myvec(vert->get_rank());
		particle_iterator it=std::unique_copy(vert->begin_sorted_legs(),vert->end_sorted_legs(),myvec.begin());
		myvec.resize(it-myvec.begin());

		/* For every unique leg, store all permutations of the remaining
		 * legs in fusion classes and insert them in the fusion map: */

		std::vector<const particle_type*>v;
		flavour_comp<particle_type>comp;
		for(particle_iterator it=myvec.begin();it != myvec.end();++it)
		{
		    v=vert->get_sorted_legs();
		    v.erase(std::find(v.begin(),v.end(),*it));
		    do
		    {
			fusion_class<model_t>fusion_object(*it,vert,v);
			fusion_map.insert(std::pair<std::vector<const particle_type*>,fusion_class<model_t> >(v,fusion_object));
		    }
		    while(std::next_permutation(v.begin(),v.end(),comp));
		}
	    }
	    
	    /* Trivial constructor: */

	    fusion_class():prod_part(NULL),vert(NULL){}

	    /* Produced particle readout: */

	    const particle_type* get_produced_particle() const
	    {
		return prod_part;
	    }

	    /* Vertex readout: */

	    const vertex_type* get_vertex() const
	    {
		return vert;
	    }

	    /* Permutation object readout: */

	    const permutation<std::vector<const particle_type*> >& get_ordering() const
	    {
		return ordering;
	    }

	    /* Function applying the permutation to another vector: */

	    template<class T>std::vector<T>& apply_ordering(std::vector<T>& v) const
	    {
		ordering(v);
		return v;
	    }

	    /* Function returning the leg integer: */

	    size_type get_leg() const
	    {
		return leg;
	    }

	    /* Comparison operator: */

	    bool operator < (const fusion_class<model_t>& other) const
	    {
		if(prod_part->get_flavour()==other.prod_part->get_flavour())
		{
		    return (*vert<*(other.vert));
		}
		else
		{
		    return (prod_part->get_flavour() < other.prod_part->get_flavour());
		}
	    }

	    /* Output method: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<vert->get_name()<<"--->"<<prod_part->get_name()<<"\t"<<leg<<std::endl;
		return os;
	    }

	private:

	    /* Produced particle address: */

	    const particle_type* prod_part;

	    /* Interaction vertex: */

	    const vertex_type* vert;
	    
	    /* Permutation converting the sorted set of incoming legs to the
	     * defined set of incoming legs: */
	    
	    permutation< std::vector<const particle_type*> > ordering;
	    
	    /* Leg vector index: */
	    
	    size_type leg;
	    
	    /* Regular constructor: stores the anti-particle of the first argument as the
	     * produced particle, copies the vertex pointer, and computes the permutation
	     * of the particle vector l w.r.t. the legs of the vertex (omitting the
	     * produced particle): */
	    
	    fusion_class(const particle_type* p,const vertex_type* v,std::vector<const particle_type*> partvec):prod_part(p->get_anti_particle()),vert(v)
	    {
		std::vector<const particle_type*>w=v->get_legs();
		typename std::vector<const particle_type*>::iterator it=std::find(w.begin(),w.end(),p);
		leg=it-w.begin();
		w.erase(it);
		ordering=permutation< std::vector<const particle_type*> >(partvec,w);
	    }
    };

    /* Streaming operator overload: */

    template<class model_t>std::ostream& operator << (std::ostream& os,const fusion_class<model_t>& f)
    {
	return f.print(os);
    }
}

#endif /*CAMGEN_FUSION_CLASS_H_*/

