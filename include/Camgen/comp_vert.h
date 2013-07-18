//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_COMP_VERT_H_
#define CAMGEN_COMP_VERT_H_

#include <string>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Vertex composition class template. The first template argument denotes an *
 * addition internal symmetry tensor factor, the second argument the inner   *
 * vertex, which should always contain a spacetime vertex factor.            *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class Feynrule_t>std::size_t rank_finder(std::size_t n)
    {
	const std::size_t* p=std::find(Feynrule_t::contractions,Feynrule_t::contractions+Feynrule_t::rank,n);
	return (p==Feynrule_t::contractions+Feynrule_t::rank)?1:Feynrule_t::ranges[p-Feynrule_t::contractions];
    }
    
    template<class Feynrule_t1,class Feynrule_t2>class compose_vertices
    {
	public:

	    /* Necessary compile-time data. Model type, vertex rank, number of couplings
	     * (params),momentum dependence (p_dependent), and fermionic property are all
	     * taken from the inner vertex contraction class: */ 

	    typedef typename Feynrule_t2::model_type model_type;
	    static const std::size_t rank=Feynrule_t2::rank;
	    static const std::size_t params=Feynrule_t2::params;
	    static const bool p_dependent=Feynrule_t2::p_dependent;
	    static const bool fermionic=Feynrule_t2::fermionic;

	    /* The vertex tensor size is simply the product of the sizes: */

	    static const std::size_t tensor_size=(Feynrule_t1::tensor_size)*(Feynrule_t2::tensor_size);

	    /* Sizes of the interacting current tensors: */

	    static const std::size_t sizes[4];

	    /* Output to model logfile: */

	    static const std::string formula;
    };
    template<class Feynrule_t1,class Feynrule_t2>const std::size_t compose_vertices<Feynrule_t1,Feynrule_t2>::rank;
    template<class Feynrule_t1,class Feynrule_t2>const std::size_t compose_vertices<Feynrule_t1,Feynrule_t2>::params;
    template<class Feynrule_t1,class Feynrule_t2>const bool compose_vertices<Feynrule_t1,Feynrule_t2>::p_dependent;
    template<class Feynrule_t1,class Feynrule_t2>const bool compose_vertices<Feynrule_t1,Feynrule_t2>::fermionic;
    template<class Feynrule_t1,class Feynrule_t2>const std::size_t compose_vertices<Feynrule_t1,Feynrule_t2>::tensor_size;
    template<class Feynrule_t1,class Feynrule_t2>const std::size_t compose_vertices<Feynrule_t1,Feynrule_t2>::sizes[4]={rank_finder<Feynrule_t1>(0)*Feynrule_t2::sizes[0],rank_finder<Feynrule_t1>(1)*Feynrule_t2::sizes[1],rank_finder<Feynrule_t1>(2)*Feynrule_t2::sizes[2],rank_finder<Feynrule_t1>(3)*Feynrule_t2::sizes[3]};
    template<class Feynrule_t1,class Feynrule_t2>const std::string compose_vertices<Feynrule_t1,Feynrule_t2>::formula=Feynrule_t1::formula + Feynrule_t2::formula;
}

#endif /*CAMGEN_COMP_VERT_H_*/

