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
#include <Camgen/evt_data.h>

namespace Camgen
{
    /// Default factorisation/renormalisation scale expression class.

    template<class model_t,std::size_t N_in,std::size_t N_out>class pT_scale: public ps_function<model_t,N_in,N_out>
    {
        public:

            typedef ps_function<model_t,N_in,N_out> base_type;
            typedef typename base_type::event_type event_type;
            typedef typename base_type::value_type value_type;
            typedef typename event_type::size_type size_type;

            /// Operator returning the qcd scale for the argument event, the square root of the average pT-squared of
            /// the outgoing particles.

            value_type operator()(const event_type& evt) const
            {
		value_type s(0);
		for(size_type i=0;i<N_out;++i)
		{
		    s+=(evt.pT2_out(i));
		}
		return std::sqrt(s/N_out);
            }
    };

    /// Utility base class combining event container, scale container and cut containers.

    template<class model_t,std::size_t N_in,std::size_t N_out>class event_generator_base: public event_container<model_t,N_in,N_out>,
                                                                                          public cut_container<model_t,N_in,N_out>,
                                                                                          public scale_container<model_t,N_in,N_out>
    {
        public:

            typedef typename event_container<model_t,N_in,N_out>::event_type event_type;
            typedef typename cut_container<model_t,N_in,N_out>::ps_cut_type ps_cut_type;
            typedef typename scale_container<model_t,N_in,N_out>::functor_type functor_type;
            typedef typename scale_container<model_t,N_in,N_out>::function_type function_type;

            /// Returns true if the current event passes the cuts.

            virtual bool pass_cuts() const
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

        protected:

            /// Implementation of the event creation factory method.
            //
            fillable_event<model_t,N_in,N_out>* create_event() const
            {
                return new event_data<model_t,N_in,N_out>();
            }

            /// Implementation of the factorization scale creation method, creating a pT-scale function object.

            ps_function<model_t,N_in,N_out>* create_F_scale() const
            {
                return new pT_scale<model_t,N_in,N_out>();
            }

            /// Implementation of the renormalisation scale creation method, creating a pT-scale function object.

            ps_function<model_t,N_in,N_out>* create_R_scale() const
            {
                return new pT_scale<model_t,N_in,N_out>();
            }
    };
}

#endif /*CAMGEN_EVT_GEN_BASE_H_*/
