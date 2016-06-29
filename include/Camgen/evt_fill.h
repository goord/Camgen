//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_fill.h
    \brief Abstract base class for fillable events in Camgen
 */

#ifndef CAMGEN_EVT_FILL_H_
#define CAMGEN_EVT_FILL_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract event subclass containing the event modifiers.     *
 * Should not be exposed in the Camgen API, only for internal  *
 * usage.                                                      *
 *                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/event.h>

namespace Camgen
{
    /// Editable event class. Contains methods for copying momenta to the event.

    template<class model_t,std::size_t N_in,std::size_t N_out>class fillable_event: public event<model_t,N_in,N_out>
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

            /// Sets the current sub-process.

            virtual void set_process(const sub_process<model_t,N_in,N_out>* sub_proc_,int i=1)
            {
                this->sub_proc=sub_proc_;
                this->procid=i;
            }

            /// Copies the argument vector to the i-th incoming momentum.

            virtual void set_p_in(const momentum_type&,size_type)=0;

            /// Copies the argument vector to the i-th outgoing momentum.

            virtual void set_p_out(const momentum_type&,size_type)=0;

            /// Sets the event weight.

            virtual void set_w(const value_type&)=0;

            /// Sets the cross section.

            virtual void set_xsec(const MC_integral<value_type>)=0;

            /// Sets the partonic invariant mass.

            virtual void set_Ecm_hat(const value_type&)=0;

            /// Sets the partonic invariant mass.

            virtual void set_Ecm_beams(const value_type&)=0;

            /// Resets to default (zero)momenta.

            void reset()
            {
                momentum_type p;
                p.assign((value_type)0);
                for(size_type i=0;i<N_in;++i)
                {
                    set_p_in(p,i);
                }
                for(size_type i=0;i<N_out;++i)
                {
                    set_p_out(p,i);
                }
                set_w((value_type)1);
                set_xsec(MC_integral<value_type>(0));
            }
    };
}

#endif /*CAMGEN_EVT_FILL_H_*/

