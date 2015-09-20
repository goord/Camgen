//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file uni_ps_fac.h
    \brief Factory for uniform phase space generators.
 */

#ifndef CAMGEN_UNI_PSGEN_FAC_H_
#define CAMGEN_UNI_PSGEN_FAC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Factory for uniform phase space generators, creating instances from *
 * configuration data or file streams.                                 *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/psgen_fac_base.h>
#include <Camgen/isgen_fac.h>
#include <Camgen/rambo.h>

namespace Camgen
{
    /// Phase space generator factory.

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class uniform_ps_generator_factory: public ps_generator_factory_base<model_t,N_in,N_out,initial_state_factory<model_t,N_in,rng_t>,uniform_ps_generator_factory<model_t,N_in,N_out,rng_t> >
    {
	public:

	    /* Type definitions. */

	    typedef ps_generator<model_t,N_in,N_out> generator_type;
	    typedef phase_space_generators::type fs_type;
	    typedef typename initial_state_factory<model_t,N_in,rng_t>::generator_type is_generator;

	    /// Creates a uniform phase space generator for the given initial
	    /// state generator and final state type (latter must be uniform).

	    static generator_type* create_instance(is_generator* isgen,fs_type fs_type_)
	    {
		switch(fs_type_)
		{
		    case phase_space_generators::uniform:
			return new rambo<model_t,N_in,N_out,rng_t>(isgen);
			break;
		    default:
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"momentum generator type "<<fs_type_<<" cannot be instantiated by this factory..."<<endlog;
			return NULL;
		}
	    }
    };
}

#endif /*CAMGEN_UNI_PSGEN_FAC_H_*/

