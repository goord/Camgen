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

#include <Camgen/sub_proc.h>

namespace Camgen
{
    /// Abstract base class for event generator configuration.
    
    template<class model_t,std::size_t N_in,std::size_t N_out>class generator_configuration
    {
	public:

	    typedef sub_process<model_t,N_in,N_out> process_type;

            /// Virtual destructor.

            virtual ~generator_configuration(){}

            /// Virtual subprocess-dependent configuration function.

	    virtual void configure(const process_type&)
            {
                configure();
            }

            /// Abstract subprocess-independent configuration method.

	    virtual void configure()=0;
    };
}

#endif /*CAMGEN_GEN_CONF_H_*/

