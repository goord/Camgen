//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file cut_container.h
    \brief Base class for cut-containing objects in Camgen.
 */

#ifndef CAMGEN_CUT_CONTAINER_H_
#define CAMGEN_CUT_CONTAINER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Base class for objects that need to hold cuts to veto events. *
 *                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <vector>
#include <Camgen/ps_cut.h>

namespace Camgen
{
    /// Base class for objects that need to own phase space cuts for vetoing events,

    template<class model_t,std::size_t N_in,std::size_t N_out>class cut_container
    {
        public:

            /* Type definitions: */

            typedef ps_cut<model_t,N_in,N_out> ps_cut_type;
            typedef typename ps_cut_type::event_type event_type;
            typedef typename event_type::value_type value_type;
            typedef typename event_type::size_type size_type;
            typedef typename ps_cut_wrapper<model_t,N_in,N_out>::function_type function_type;

            /// Constructor, by default no cuts are present.
            
            cut_container(){}

            /// Destructor.

            virtual ~cut_container()
            {
                for(size_type i=0;i<cuts.size();++i)
                {
                    if(alloc_cuts[i])
                    {
                        delete cuts[i];
                        cuts[i]=NULL;
                        alloc_cuts[i]=false;
                    }
                }
            }

            /// Adds a non-owned phase space cut.

            void add_cut(const ps_cut_type* psc)
            {
                cuts.push_back(psc);
                alloc_cuts.push_back(false);
                after_cut_inserted(psc);
            }

            /// Adds a phase space cut function.

            void add_cut(function_type f)
            {
                const ps_cut_type* psc=new ps_cut_wrapper<model_t,N_in,N_out>(f);
                cuts.push_back(psc);
                alloc_cuts.push_back(true);
                after_cut_inserted(psc);
            }

            /// Returns true if the argument event passes all cuts.

            bool pass(const event_type& evt) const
            {
                for(size_type i=0;i<cuts.size();++i)
                {
                    if(!(*(cuts[i]))(evt))
                    {
                        return false;
                    }
                }
                return true;
            }

        protected:

            /// Virtual function for implementing side-effects of cut-insertions.

            virtual void after_cut_inserted(const ps_cut_type* psc){}

        private:

            /* Phase space cuts: */

            std::vector<const ps_cut_type*> cuts;

            /* Flags denoting which cuts are owned by this instance: */

            std::vector<bool> alloc_cuts;
    }; 
}

#endif /*CAMGEN_CUT_CONTAINER_H_*/

