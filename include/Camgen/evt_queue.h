//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_queue.h
  \brief Event interface copying input events into a queue on the heap.
*/

#ifndef CAMGEN_EVT_QUEUE_H_
#define CAMGEN_EVT_QUEUE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Partial implementation of the event stream interface that *
 * pushes the events to a queue in memory.                   *
 *                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_stream.h>
#include <queue>

namespace Camgen
{
    /// Partial event interface implementation copying events to a queue in memory.

    template<class model_t,std::size_t N_in,std::size_t N_out>class event_queue: public event_stream<model_t,N_in,N_out>
    {
        typedef event_stream<model_t,N_in,N_out> base_type;

        public:

            /* Public type definitions: */

            typedef typename base_type::event_type event_type;
            typedef typename std::queue<const event_type*>::size_type size_type;

            /// Virtual destructor.

            virtual ~event_queue()
            {
                while(!empty())
                {
                    delete evt_queue.front();
                    evt_queue.pop();
                }
            }

            /// Returns whether the queue contains events.

            bool empty() const
            {
                return evt_queue.empty();
            }

            /// Returns the number of events in the queue.

            size_type size() const
            {
                return evt_queue.size();
            }

            /// Returns the front of the queue (oldest event).

            const event_type* front() const
            {
                return evt_queue.front();
            }

            /// Returns the back of the queue (newest event).

            const event_type* back() const
            {
                return evt_queue.back();
            }

            /// Pops the front of the queue.

            bool pop()
            {
                if(empty())
                {
                    return false;
                }
                delete evt_queue.front();
                evt_queue.pop();
                return true;
            }

            /// Event size implementation.

            size_type event_size() const
            {
                return (N_in+N_out)*sizeof(typename event_type::momentum_type);
            }

            /// Write implementation, does nothing.

            virtual bool write()
            {
                while(!empty())
                {
                    std::cout<<*front()<<std::endl;
                    pop();
                }
                return true;
            }

        protected:

            /// Pushes to the back of the queue.

            bool fill_event(const event_type& evt)
            {
                evt_queue.push(evt.clone());
                return true;
            }

        private:

            std::queue<const event_type*> evt_queue;
    };
}

#endif /*CAMGEN_EVT_QUEUE_H_*/

