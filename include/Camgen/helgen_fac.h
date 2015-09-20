//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_HELGEN_FAC_H_
#define CAMGEN_HELGEN_FAC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Helicity generator factory. Creates uniform or importance sampling  *
 * helicity generators or sequential spin summation objects.           *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/MC_config.h>
#include <Camgen/uni_hels.h>
#include <Camgen/long_hels.h>

namespace Camgen
{
    /// Helicity generator factory.

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class helicity_generator_factory
    {
	public:

	    typedef typename model_t::value_type value_type;
	    typedef helicity_generator<typename model_t::value_type,N_in,N_out,model_t::continuous_helicities> generator_type;
	    typedef helicity_generators::type generator_tag;

	    static const bool ch=(model_t::dimension==0)?false:(model_t::continuous_helicities);

	    /// Factory method returning a helicity generator instance from a
	    /// tree iterator.

	    static generator_type* create_generator(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		return create_generator(it,helicity_generator_type());
	    }

	    /// Factory method returning a helicity generator instance from a
	    /// tree iterator and generator tag.

	    static generator_type* create_generator(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it,generator_tag tag)
	    {
		switch(tag)
		{
		    case helicity_generators::uniform:
			return uniform_helicities<value_type,N_in,N_out,rng_t,ch>::template create_instance<model_t>(it);
		    case helicity_generators::longitudinal:
			return longitudinal_helicities<value_type,N_in,N_out,rng_t,ch>::template create_instance<model_t>(it);
		    case helicity_generators::summation:
			return helicity_summer<value_type,N_in,N_out,ch>::template create_instance<model_t>(it);
		    default:
			return NULL;
		}
	    }
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const bool helicity_generator_factory<model_t,N_in,N_out,rng_t>::ch;
}

#endif /*CAMGEN_HELGEN_FAC_H_*/

