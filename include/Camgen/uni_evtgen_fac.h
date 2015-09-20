//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file uni_procgen_fac.h
    \brief Factory for uniform single-process generators.
 */

#ifndef CAMGEN_UNI_PROCGEN_FAC_H_
#define CAMGEN_UNI_PROCGEN_FAC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Child class of the process generator base which uses the phase space, *
 * helicity and colour generator factories to create sub-generators.     *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/procgen_fac_base.h>
#include <Camgen/evtgen_fac_base.h>
#include <Camgen/uni_psgen_fac.h>
#include <Camgen/helgen_fac.h>
#include <Camgen/colgen_fac.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class uniform_process_generator_factory: public process_generator_factory_base<model_t,N_in,N_out,rng_t>
    {
	typedef process_generator_factory_base<model_t,N_in,N_out,rng_t> base_type;

	public:

	    /* Public type definitions: */
	    /*--------------------------*/

	    typedef typename base_type::process_generator_type process_generator_type;
	    typedef typename base_type::momentum_generator_type momentum_generator_type;
	    typedef typename base_type::helicity_generator_type helicity_generator_type;
	    typedef typename base_type::colour_generator_type colour_generator_type;
	    typedef typename base_type::CM_tree_iterator CM_tree_iterator;

	    /* Public factory methods: */
	    /*--------------------------*/

	    /// Creates a phase space generator for the given single-process amplitude.

	    momentum_generator_type* create_momentum_generator(CM_tree_iterator amplitude)
	    {
		uniform_ps_generator_factory<model_t,N_in,N_out,rng_t>factory;
		return factory.create_generator(amplitude);
	    }

	    /// Creates a helicity generator for the given single-process
	    /// amplitude.
	    
	    helicity_generator_type* create_helicity_generator(CM_tree_iterator amplitude)
	    {
		return helicity_generator_factory<model_t,N_in,N_out,rng_t>::create_generator(amplitude);
	    }
	    
	    /// Creates a colour generator for the given single-process
	    /// amplitude.
	    
	    colour_generator_type* create_colour_generator(CM_tree_iterator amplitude)
	    {
		return colour_generator_factory<model_t,N_in,N_out,rng_t>::create_generator(amplitude);
	    }
    };

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class uniform_event_generator_factory: public event_generator_factory_base<model_t,N_in,N_out,rng_t>
    {
	typedef event_generator_factory_base<model_t,N_in,N_out,rng_t> base_type;

	public:

	    /* Public type definitions: */
	    /*--------------------------*/

	    typedef typename base_type::process_generator_type process_generator_type;
	    typedef typename base_type::momentum_generator_type momentum_generator_type;
	    typedef typename base_type::helicity_generator_type helicity_generator_type;
	    typedef typename base_type::colour_generator_type colour_generator_type;
	    typedef typename base_type::CM_tree_iterator CM_tree_iterator;
	    typedef typename base_type::event_generator_type event_generator_type;
	    typedef typename base_type::amplitude_type amplitude_type;

	    /// Constructor
	    
	    uniform_event_generator_factory():base_type(new uniform_process_generator_factory<model_t,N_in,N_out,rng_t>()){}

	    /// Destructor

	    ~uniform_event_generator_factory()
	    {
		delete this->process_generator_factory;
	    }
    };
}

#endif /*CAMGEN_UNI_PROCGEN_FAC_H_*/
