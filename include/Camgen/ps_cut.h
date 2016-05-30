//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_cut.h
    \brief Abstract interface base class for phase space cut objects.
 */

#ifndef CAMGEN_PS_CUT_H_
#define CAMGEN_PS_CUT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for phase space cut objects. Derived classes should *
 * implement the pass() method.                                            *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <Camgen/event.h>

namespace Camgen
{
    /// Abstract base class for phase space cuts.

    template<class model_t,std::size_t N_in,std::size_t N_out> class ps_cut
    {
        public:

            /* Type definitions: */
            
            typedef event<model_t,N_in,N_out> event_type;

            /// Virtual destructor.

            virtual ~ps_cut(){}

            /// Abstract function for passing the imposed cuts.

            virtual bool operator()(const event_type&)=0;
    };

    /// Phase space cut implementation wrapping a function pointer.

    template<class model_t,std::size_t N_in,std::size_t N_out> class ps_cut_wrapper
    {
        public:

            /* Type definitions: */
            
            typedef event<model_t,N_in,N_out> event_type;
            typedef bool (*function_type)(const event_type&);

            /// Constructor, taking a function pointer as argument.

            ps_cut_wrapper(function_type function_):function(function_){}

            /// Abstract function for passing the imposed cuts.

            virtual bool operator()(const event_type& evt)
            {
                return function(evt);
            }

        protected:

            /* Wrapped function pointer: */

            function_type function;
    };
}

#endif /*CAMGEN_PS_CUT_H_*/
