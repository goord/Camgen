//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_func.h
    \brief Abstract base class for phase space functions in Camgen
 */

#ifndef CAMGEN_PS_FUNC_H_
#define CAMGEN_PS_FUNC_H_

#include <Camgen/event.h>

namespace Camgen
{
    /// Abstract base class for real-valued functions of events in Camgen.
    
    template<class model_t,std::size_t N_in,std::size_t N_out>class ps_function
    {
        public:

            /* Type definitions: */
            
            typedef event<model_t,N_in,N_out> event_type;
            typedef typename event_type::value_type value_type;

            /// Virtual destructor.

            virtual ~ps_function(){}

            /// Function evaluation operator.

            virtual value_type operator()(const event_type&) const=0;
    };

    /// Phase space function implementation wrapping a function pointer.

    template<class model_t,std::size_t N_in,std::size_t N_out>class ps_function_wrapper: public ps_function<model_t,N_in,N_out>
    {
        public:

            /* Type definitions: */
            
            typedef event<model_t,N_in,N_out> event_type;
            typedef typename event_type::value_type value_type;
            typedef value_type (*function_type)(const event_type&);

            /// Constructor, taking a function pointer as argument.

            ps_function_wrapper(function_type function_):function(function_){}

            /// Evaluates the function pointer.

            value_type operator()(const event_type& evt) const
            {
                return function(evt);
            }

        protected:

            /* Wrapped function pointer: */

            function_type function;
    };
}

#endif /*CAMGEN_PS_FUNC_H_*/

