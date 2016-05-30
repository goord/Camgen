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

#include <Camgen/ps_cut.h>
#include <Camgen/ps_gen_base.h>

namespace Camgen
{
    template<class model_t>class generator_configuration: public ps_generator_viewer<model_t>
    {
	public:

	    typedef model_t model_type;
	    typedef typename ps_generator_viewer<model_t>::momentum_type momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef typename ps_generator_viewer<model_t>::spacetime_type spacetime_type;

	    generator_configuration(){}

	    virtual void configure()=0;
    };
}

#endif /*CAMGEN_GEN_CONF_H_*/

