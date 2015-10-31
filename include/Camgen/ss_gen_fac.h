//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ss_gen_fac.h
    \brief Invariant mass pair sampler factory
 */

#ifndef CAMGEN_SS_GEN_FAC_H_
#define CAMGEN_SS_GEN_FAC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Factory class creating s-pair sampler instances depending on* 
 * the s_pair_generator_mode flag in MC_config.                *
 *                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/ss_gen.h>
#include <Camgen/MC_config.h>

namespace Camgen
{

    /// Factory for s-pair MC samplers, depending on the current configuration.

    template<class value_t,class rng_t>class s_pair_generator_factory
    {
	public:

	    typedef value_t value_type;
	    typedef value_generator<value_t,rng_t> s_generator_type;
	    typedef s_pair_generator<value_t,rng_t> s_pair_generator_type;
	    typedef Dirac_delta<value_t,rng_t> Dirac_delta_type;

	    /* Factory method: */

	    static s_pair_generator_type* create(s_generator_type* s1_generator,s_generator_type* s2_generator,const value_type& s)
	    {
		s_generator_type* s1;
		s_generator_type* s2;
		if(dynamic_cast<Dirac_delta_type*>(s1_generator)==NULL and dynamic_cast<Dirac_delta_type*>(s2_generator)!=NULL)
		{
		    s1=s2_generator;
		    s2=s1_generator;
		}
		else
		{
		    s1=s1_generator;
		    s2=s2_generator;
		}
		switch(s_pair_generation_mode())
		{
		    case s_pair_generation_modes::symmetric:
			return new symmetric_s_generator_composition<value_t,rng_t>(s1,s2,s);
		    case s_pair_generation_modes::asymmetric:
			return new s_generator_composition<value_t,rng_t>(s1,s2,s);
		    case s_pair_generation_modes::hit_and_miss:
			return new symmetric_s_pair_generator<value_t,rng_t>(s1,s2,s);
		    default:
			return NULL;
		}
	    }
    };
}

#endif /*CAMGEN_SS_GEN_FAC_H_*/

