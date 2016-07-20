//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PS_INTERFACE_H_
#define CAMGEN_PS_INTERFACE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility base class for parton-shower interfaces.  *
 *                                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_queue.h>
#include <Camgen/evt_gen_base.h>
#include <Camgen/MC_gen.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out>class parton_shower_interface
    {
        public:

            typedef event_generator_base<model_t,N_in,N_out> generator_type; 
            typedef typename generator_type::event_type event_type;
            typedef typename event_type::value_type value_type;

	    /// Les-Houches accord weight switch: 
	    /// 1: weighted events to PS, unweighted out, no xsec provided.
	    /// 2: weighted events to PS, unweighted out, xsecs provided.
	    /// 3: unweighted events to PS, unweighted out.
	    /// 4: weighted events to PS, weighted out.

	    const int weight_switch;

	    /// Generation switch. If set true, the matrix element generation is
	    /// triggered by the parton shower, otherwise, it has to be performed
	    /// by the user.

	    bool generate_events;

            /// Constructor.

            parton_shower_interface(generator_type* gen_,int weight_switch_,bool generate_events_=true):weight_switch(weight_switch_),generate_events(generate_events_),gen(gen_){}

            /// Generates the next event: if generate_events is set true, lets the generator instance generate a new
            /// event, otherwise returns the oldest event in the internal queue.
            
            const event_type* next_event()
            {
                if(generate_events)
                {
                    generate_event();
                    return const_cast<const generator_type*>(gen)->get_event_ptr();
                }
                else
                {
                    if(evt_queue.empty())
                    {
                        return NULL;
                    }
                    return evt_queue.front();
                }
            }

            /// Copies the existing generator's event to the queue.

            const event_type* fill()
            {
                evt_queue.fill_event(gen->get_event());
                return evt_queue.back();
            }

            /// Pops the oldest event from the queue.

            void pop_event()
            {
                if(!evt_queue.empty())
                {
                    evt_queue.pop();
                }
            }

            /// Returns the most actual event.

            const event_type* get_current_event() const
            {
                return const_cast<const generator_type*>(gen)->get_event_ptr();
            }

        protected:

            /* Event generator instance. */

            generator_type* const gen;

        private:

            /* Internal event queue */

            event_queue<model_t,N_in,N_out> evt_queue;

            /// Generates new event according to the given strategy argument.

            bool generate_event()
            {
                MC_generator<value_type>* mcgen=dynamic_cast<MC_generator<value_type>*>(gen);
                if(mcgen==NULL)
                {
                    return false;
                }
                if(std::abs(weight_switch)!=3)
                {
                    return mcgen->generate();
                }
                return mcgen->generate_unweighted();
            }
    };
}

#endif /*CAMGEN_PS_INTERFACE_H_*/

