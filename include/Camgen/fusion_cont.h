//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FUSION_CONT_H_
#define CAMGEN_FUSION_CONT_H_

#include <utility>
#include <Camgen/fusion_class.h>
#include <Camgen/flav_comp.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Unary function objects for fusion classes, requesting whether a fusion  *
 * contains a particle as incoming or a vertex.                            *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Function object requesting whether a particle (constructor argument) is
     * contained in the vector of incoming particles: */
    
    template<class model_t>class fusion_contains_particle
    {
	public:
	    typedef const particle<model_t>* particle_type;

	    fusion_contains_particle(const particle_type& phi_):phi(phi_){}
	    bool operator() (const std::pair< std::vector<const particle<model_t>*>,fusion_class<model_t> >& fusion) const
	    {
		if(std::binary_search(fusion.first.begin(),fusion.first.end(),phi,comp))
		{
		    return true;
		}
		return (fusion.second.get_produced_particle()==phi);
	    }
	private:
	    const particle_type& phi;
	    flavour_comp< particle<model_t> >comp;
    };

    /* Function object requested whether a vertex (constructor) pointer is
     * member of the fusion class (operator argument): */

    template<class model_t>class fusion_contains_vertex
    {
	public:
	    typedef const vertex<model_t>* vertex_type;

	    fusion_contains_vertex(const vertex_type& vert_):vert(vert_){}
	    bool operator() (const std::pair< std::vector<const particle<model_t>*>,fusion_class<model_t> >& fusion) const
	    {
		return (fusion.second.get_vertex()==vert);
	    }
	private:
	    const vertex_type& vert;
    };
}

#endif /*CAMGEN_FUSION_CONT_H_*/


