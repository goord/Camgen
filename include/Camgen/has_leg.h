//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_HAS_LEG_H_
#define CAMGEN_HAS_LEG_H_

#include <algorithm>
#include <vector>
#include <Camgen/forward_decs.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Predicate classes determining whether a vertex gives an interaction with a    *
 * given particle, or vector of particles. These can be used in queries through  *
 * containers of vertex pointers, looking or the interactions of a given         *
 * particle.                                                                     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Binary function object determining whether a vertex contains a particle
     * in its legs: */

    template<class model_t>class has_leg
    {
	public:
	    has_leg(const particle<model_t>*phi_):phi(phi_){}
	    bool operator () (const vertex<model_t>* v)
	    {
		return std::binary_search(v->sorted_legs.begin(),v->sorted_legs.end(),phi,comp);
	    }
	private:
	    const particle<model_t>* phi;
	    flavour_comp< particle<model_t> > comp;
    };

    /* Binary function object determining whether a vertex legs are a
     * permutation of a given sorted particle vector: */

    template<class model_t>class has_legs
    {
	public:
	    has_legs(const std::vector<const particle<model_t>*>& phi_):phi(phi_){}
	    bool operator () (const vertex<model_t>* v)
	    {
		return v->legs==phi;
	    }
	private:
	    const std::vector<const particle<model_t>*>& phi;
    };
}

#endif /*CAMGEN_HAS_LEG_H_*/

