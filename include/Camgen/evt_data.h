//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_data.h
    \brief Event class implementation.
 */

#ifndef CAMGEN_EVT_DATA_H_
#define CAMGEN_EVT_DATA_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Implementation of the event class, containing a set of incoming and outgoing momenta. *
 *                                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_fill.h>

namespace Camgen
{
    /// Event interface implementation containing arrays of incoming and outgoing momenta, owned by the instance.

    template<class model_t,std::size_t N_in,std::size_t N_out>class event_data: public fillable_event<model_t,N_in,N_out>
    {
        public:

            /* Type definitions: */

            typedef event<model_t,N_in,N_out> base_type;
	    typedef typename base_type::model_type model_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::spacetime_type spacetime_type;
            typedef typename base_type::particle_type particle_type;

            /// Constructor. Allocates all momenta to zero.

            event_data()
            {
                for(size_type i=0;i<N_in;++i)
                {
                    pin[i].assign((value_type)0);
                }
                for(size_type i=0;i<N_out;++i)
                {
                    pout[i].assign((value_type)0);
                }
            }

            /// Implementation of the weight accessor.

            value_type w() const
            {
                return weight;
            }

            /// Implementation of the cross section accessor.

            MC_integral<value_type> xsec() const
            {
                return cross_section;
            }

            /// Returns the partonic invariant mass.

            value_type Ecm_hat() const
            {
                return ecm;
            }

            /// Implementation of the incoming momentum accessor.

            momentum_type p_in(size_type i) const
            {
                return pin[i];
            }

            /// Implementation of the outgoing momentum accessor.

            momentum_type p_out(size_type i) const
            {
                return pout[i];
            }

            /// Sets the partonic invariant mass.

            void set_Ecm_hat(const value_type& e)
            {
               ecm=e;
            }

            /// Sets the event weight.

            void set_w(const value_type& weight_)
            {
                weight=weight_;
            }

            /// Updates the process cross_section.

            void set_xsec(const MC_integral<value_type> cross_section_)
            {
                cross_section=cross_section_;
            }

            /// Copies the argument vector to the i-th incoming momentum.

            void set_p_in(const momentum_type& p,size_type i)
            {
                pin[i]=p;
            }

            /// Copies the argument vector to the i-th outgoing momentum.

            void set_p_out(const momentum_type& p,size_type i)
            {
                pout[i]=p;
            }

        private:

            /* Momentum collections: */

            vector<momentum_type,N_in> pin;
            vector<momentum_type,N_out> pout;

            /* Invariant mass: */

            value_type ecm;

            /* Event weight: */

            value_type weight;

            /* Process cross_section: */

            MC_integral<value_type> cross_section;
    };
}

#endif /*CAMGEN_EVT_DATA_H_*/

