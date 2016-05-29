//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file scale_container.h
    \brief Base class containing expressions for renormalization and factorization scales.
 */

#ifndef CAMGEN_SCALE_CONTAINER_H_
#define CAMGEN_SCALE_CONTAINER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Base class containing phase space functions for factorization and renormalization scales. *
 *                                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <Camgen/ps_func.h>

namespace Camgen
{
    /// Base class responsible for ownership management of factorization/renormalization scale expressions.
    
    template<class model_t,std::size_t N_in,std::size_t N_out>class scale_container
    {
	public:

	    /* Type definitions: */

            typedef ps_function<model_t,N_in,N_out> functor_type;
            typedef typename functor_type::event_type event_type;
            typedef typename functor_type::value_type value_type;
            typedef value_type (*function_type)(const event_type&);

            /// Constructor. By default, both scales are the Z-mass.
            
            scale_container():mu_F(NULL),alloc_mu_F(false),mu_R(NULL),alloc_mu_R(false){}

	    /// Destructor.

	    virtual ~scale_container()
            {
                delete_F_scale();
                delete_R_scale();
            }

            /// Sets the factorization expression to the argument phase space expression.

            void set_F_scale(const functor_type* f)
            {
                set_F_scale(f,false);
            }

            /// Sets the factorization expression to the given event function.

            void set_F_scale(function_type f)
            {
                set_F_scale(new ps_function_wrapper<model_t,N_in,N_out>(f),true);
            }

            /// Sets the factorization expression to the newly created one by the implementer.

            void allocate_F_scale()
            {
                set_F_scale(create_F_scale(),true);
            }

            /// Sets the renormalization expression to the argument phase space expression.

            void set_R_scale(const functor_type* f)
            {
                set_R_scale(f,false);
            }

            /// Sets the renormalization expression to the given event function.

            void set_R_scale(function_type f)
            {
                set_R_scale(new ps_function_wrapper<model_t,N_in,N_out>(f),true);
            }

            /// Sets the renormalization expression to the newly created one by the implementer.

            void allocate_R_scale()
            {
                set_R_scale(create_R_scale(),true);
            }

            /// Sets both factorization and renormalization scales to the argument expression.

            void set_QCD_scale(const functor_type* f)
            {
                set_F_scale(f);
                set_R_scale(f);
            }

	    /// Returns the factorisation scale for the given event.

	    value_type F_scale(const event_type& evt) const
	    {
		return mu_F==NULL?default_scale():(*mu_F)(evt);
	    }

	    /// Returns the renormalization scale for the given event.

	    value_type R_scale(const event_type& evt) const
	    {
		return mu_R==NULL?default_scale():(*mu_R)(evt);
	    }

	    /// Default scale value: Z-boson pole mass.

	    virtual value_type default_scale() const
	    {
		return (value_type)91.1876;
	    }

        protected:

            /// Helper function for setting the factorization scale.

            void set_F_scale(const functor_type* f,bool owner)
            {
                delete_F_scale();
                mu_F=f;
                alloc_mu_F=owner and (mu_F!=NULL);
                after_F_scale_set(mu_F);
            }

            /// Implements side effects after setting the factorization scale.

            virtual void after_F_scale_set(const functor_type* f){}

            /// Creates a new factorization scale expression. NULL if not overridden.

            virtual functor_type* create_F_scale() const
            {
                return NULL;
            }

            /// Helper function for setting the renormalization scale.

            void set_R_scale(const functor_type* f,bool owner)
            {
                delete_R_scale();
                mu_R=f;
                alloc_mu_R=owner and (mu_R!=NULL);
                after_R_scale_set(mu_R);
            }

            /// Implements side effects after setting the renormalization scale.

            virtual void after_R_scale_set(const functor_type* f){}

            /// Creates a new renormalization scale expression. NULL if not overridden.

            virtual functor_type* create_R_scale() const
            {
                return NULL;
            }

        private:

            /* Factorization scale phase space function: */
            // TODO: use smart pointer!

            const functor_type* mu_F;
            bool alloc_mu_F;

            /* Renormalization scale phase space function: */
            // TODO: use smart pointer!

            const functor_type* mu_R;
            bool alloc_mu_R;

            /* Helper function for cleaning up the factorization scale functor: */

            void delete_F_scale()
            {
                if(alloc_mu_F)
                {
                    delete mu_F;
                    alloc_mu_F=false;
                }
                mu_F=NULL;
            }

            /* Helper function for cleaning up the renormalization scale functor: */

            void delete_R_scale()
            {
                if(alloc_mu_R)
                {
                    delete mu_R;
                    alloc_mu_R=false;
                }
                mu_R=NULL;
            }
    };
}

#endif /*CAMGEN_SCALE_CONTAINER_H_*/

