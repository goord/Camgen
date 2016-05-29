//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_gen_base.h
    \brief Utility base class for event-generating objects in Camgen.
 */

#ifndef CAMGEN_EVT_GEN_BASE_H_
#define CAMGEN_EVT_GEN_BASE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility class combining event-, scale- and cut-container objects in Camgen. *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_container.h>
#include <Camgen/cut_container.h>
#include <Camgen/scale_container.h>

namespace Camgen
{
    /// Utility base class combining event container, scale container and cut containers.

    template<class model_t,std::size_t N_in,std::size_t N_out>class event_generator_base: public event_container<model_t,N_in,N_out>,
                                                                                          public cut_container<model_t,N_in,N_out>,
                                                                                          public scale_container<model_t,N_in,N_out>
    {
        /// Returns true if the current event passes the cuts.

        bool pass() const
        {
            return this->pass(*(this->get_event()));
        }

        /// Returns the factorization scale for the current event.
        
        bool F_scale() const
        {
            return this->F_scale(*(this->get_event()));
        }

        /// Returns the renormalization scale for the current event.

        bool R_scale() const
        {
            return this->R_scale(*(this->get_event()));
        }
    };
}

#endif /*CAMGEN_EVT_GEN_BASE_H_*/
