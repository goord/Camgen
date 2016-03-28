//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evtgen_fac.h
    \brief Multi-process event generator factory.
 */

#ifndef CAMGEN_EVTGEN_FAC_H_
#define CAMGEN_EVTGEN_FAC_H_

/* * * * * * * * * * * * * * * * *
 * Factory for event generators. *
 *                               *
 * * * * * * * * * * * * * * * * */

#include <Camgen/evtgen_fac_base.h>
#include <Camgen/procgen_fac.h>

namespace Camgen
{
    /// Multi-process event generator factory base.
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class event_generator_factory: public event_generator_factory_base<model_t,N_in,N_out,rng_t>
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
	    
	    event_generator_factory():base_type(new process_generator_factory<model_t,N_in,N_out,rng_t>()){}

	    /// Destructor

	    ~event_generator_factory()
	    {
		delete this->process_generator_factory;
	    }
    };
}

#endif /*CAMGEN_PROCGEN_FAC_BASE_H_*/

