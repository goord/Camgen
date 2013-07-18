//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file gen_conf.h
    \brief Event generator configuration base class.
 */

#ifndef CAMGEN_GEN_CONF_H_
#define CAMGEN_GEN_CONF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Base class for automatic event generator configuration/initialisation.  *
 * Derived objects may be used to create and initialise event generator    *
 * instances.                                                              *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/MC_config.h>
#include <Camgen/ps_cut.h>
#include <Camgen/ps_gen_base.h>
#include <Camgen/scale_expr.h>

namespace Camgen
{
    template<class model_t>class generator_configuration: public phase_space_cut,
    							  public scale_expression<typename model_t::value_type>,
    							  public ps_generator_viewer<model_t>
    {
	public:

	    typedef model_t model_type;
	    typedef typename ps_generator_viewer<model_t>::momentum_type momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef typename ps_generator_viewer<model_t>::spacetime_type spacetime_type;
	    typedef helicity_generators::type helicity_generator_tag;
	    typedef colour_generators::type colour_generator_tag;
	    typedef initial_states::type initial_state_tag;
	    typedef phase_space_generators::type ps_generator_tag;

	    generator_configuration(){}

	    virtual void configure()=0;
    };

}

#endif /*CAMGEN_GEN_CONF_H_*/

