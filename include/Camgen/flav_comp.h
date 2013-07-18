//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FLAV_COMP_H_
#define CAMGEN_FLAV_COMP_H_

#include <vector>
#include <algorithm>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Flavour comparison object classes. These make all stl algorithms available in *
 * containers of particle pointers, sorted by flavour.                           *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Less-than comparison class: */

    template<class particle_t>class flavour_comp
    {
	public:
	    flavour_comp(){}
	    bool operator () (const particle_t* phi1,const particle_t* phi2) const
	    {
		return (phi1->flavour < phi2->flavour);
	    }
    };

    /* Equality comparison class: */

    template<class particle_t>class flavour_eq
    {
	public:
	    flavour_eq(){}
	    bool operator () (const particle_t* phi1,const particle_t* phi2) const
	    {
		return (phi1->flavour==phi2->flavour);
	    }
    };

    /* Unequality comparison class: */

    template<class particle_t>class flavour_neq
    {
	public:
	    flavour_neq(){}
	    bool operator () (const particle_t* phi1,const particle_t* phi2) const
	    {
		return (phi1->flavour!=phi2->flavour);
	    }
    };

    /* Selection class for finding algorithms: */

    template<class particle_t>class flavour_select
    {
	public:
	    int flavour;
	    flavour_select(int x):flavour(x){}
	    flavour_select(const particle_t* f):flavour(f->flavour){}
	    bool operator () (const particle_t* phi) const
	    {
		return (phi->flavour==flavour);
	    }
    };

    /* Less-than extension to vectors: */

    template<class particle_t>class flavourvec_comp
    {
	public:
	    bool operator () (const std::vector<const particle_t*>& v1,const std::vector<const particle_t*>& v2)
	    {
		return std::lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end(),elem_comp);				
	    }
	private:
	    flavour_comp<particle_t>elem_comp;
    };
}

#endif /*CAMGEN_FLAV_COMP_H_*/
