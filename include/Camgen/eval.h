//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_EVAL_H_
#define CAMGEN_EVAL_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declarations and definitions of the evaluate and cfd_evaluate class           *
 * templates. The evaluate class takes a vertex type (possibly, a composed type) *
 * as its template argument, and incorporates the corresponding recursive        *
 * equations. These are however absent in the default definition, and therefore  *
 * each usable vertex class should also contain a specialisation of the evaluate *
 * class template with the static functions first, second, third and fourth,     *
 * denoting the recursions and a static vector-valued member function            *
 * get_index_ranges, filling the argument vector with the correct index ranges   *
 * of the amplitudes involved in the combination.                                *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */                                                                               

namespace Camgen
{

    /* Definition and declaration of the by default trivial evaluate class
     * template. If this class is instantiated in Camgen, it shall result in a
     * compile error, because the static members get_index_ranges, first,
     * second, this and fourth are not present. Every specialisation should
     * hence contain these additional static members: */

    template<class Feynrule_t>class evaluate
    {
	public:
	    DEFINE_BASIC_TYPES(typename Feynrule_t::model_type);
	    static void initialise(){}
    };

    /* Declaration and definition of the cfd_evaluate class, which is invoked
     * instead of the evaluate class if the colour_flow treatment is used with
     * discrete colours. The recursive static member function contain an extra
     * argument, a tensor iterator std::set called produced_iters in
     * CFD_ARG_LIST, denoting the propagating colours in the produced amplitude.
     * By default, the class template invokes the evaluate classes' recursive
     * member functions and inserts the produced iterator in this set: */

    template<class Feynrule_t>class cfd_evaluate
    {
	public:
	    DEFINE_BASIC_TYPES(typename Feynrule_t::model_type);

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(const size_type n,std::vector<size_type>& r)
	    {
		return evaluate<Feynrule_t>::get_index_ranges(n,r);
	    }
	    static void first(CFD_ARG_LIST)
	    {
		produced_iters.insert(iters[0]);
		evaluate<Feynrule_t>::first(factor,couplings,iters,momenta);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		produced_iters.insert(iters[1]);
		evaluate<Feynrule_t>::second(factor,couplings,iters,momenta);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		produced_iters.insert(iters[2]);
		evaluate<Feynrule_t>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		produced_iters.insert(iters[3]);
		evaluate<Feynrule_t>::fourth(factor,couplings,iters,momenta);
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_EVAL_H_*/

