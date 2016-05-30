//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_container.h
    \brief Abstract base class for event-containing objects in Camgen
 */

#ifndef CAMGEN_EVT_CONTAINER_H_
#define CAMGEN_EVT_CONTAINER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility abstract base class for event-owning (phase space) generators *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_fill.h>

namespace Camgen
{
    /// Utility abstract base class for event containers. Implementors must provide a default event creation method.

    template<class model_t,std::size_t N_in,std::size_t N_out> class event_container
    {
        public:

            /* Type definitions: */

            typedef event<model_t,N_in,N_out> event_type;

            /// Constructor. Creates the owned event.

            event_container():evt(NULL),alloc_evt(false){}

            /// Destructor. Deletes the event if necessary.

            virtual ~event_container()
            {
                if(alloc_evt)
                {
                    delete evt;
                    evt=NULL;
                    alloc_evt=false;
                }
            }

            /// Returns the event (read-only).

            const event_type* get_event() const
            {
                return evt;
            }

            /// Sets the event instance to the argument pointer, deleting the current instance whenever necessary. 
            /// If NULL, re-allocates the event. 

            void set_event(fillable_event<model_t,N_in,N_out>* evt_)
            {
                if(evt_==NULL)
                {
                    if(!alloc_evt)
                    {
                        evt=create_event();
                        alloc_evt=true;
                    }
                }
                else
                {
                    if(alloc_evt)
                    {
                        delete evt;
                        alloc_evt=false;
                    }
                    evt=evt_;
                }
                after_event_set(evt);
            }

            /// Allocates the event by calling the create_event implementation.

            void allocate_event()
            {
                set_event(NULL);
            }

            /// Utility method for post-event setter operations.

            virtual void after_event_set(fillable_event<model_t,N_in,N_out>* evt_){}

        protected:

            /// Abstract event factory method.

            virtual fillable_event<model_t,N_in,N_out>* create_event() const=0;

            /// Returns the event (read/write)

            fillable_event<model_t,N_in,N_out>* get_event()
            {
                return evt;
            }

        private:

            //TODO: Use smart pointers!

            /* Event pointer. */

            fillable_event<model_t,N_in,N_out>* evt;

            /* Flag denoting whether the evt pointer is owned by this class instance. */

            bool alloc_evt;
    };
}

#endif /*CAMGEN_EVT_CONTAINER_H_*/

