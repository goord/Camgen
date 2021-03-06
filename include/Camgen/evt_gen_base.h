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
                return (value_type)91.1876;
            }
    };

    template<class model_t,std::size_t N_out>class pT_scale<model_t,2,N_out>: public ps_function<model_t,2,N_out>
    {
        public:

            typedef ps_function<model_t,2,N_out> base_type;
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
            typedef typename event_type::size_type size_type;
            typedef typename event_type::value_type value_type;

            virtual ~event_generator_base(){}

            /// Returns true if the current event passes the cuts.

            virtual bool pass_cuts() const
            {
                return this->pass(this->get_event());
            }

            /// Returns the factorization scale for the current event.

            bool F_scale() const
            {
                return this->scale_container<model_t,N_in,N_out>::F_scale(this->get_event());
            }

            /// Returns the renormalization scale for the current event.

            bool R_scale() const
            {
                return this->scale_container<model_t,N_in,N_out>::R_scale(this->get_event());
            }

            /// Returns the number of subprocesses.

            virtual size_type processes() const
            {
                return 1;
            }

            /// Returns the id oft he i-th subprocess.

            virtual int process_id(size_type i) const
            {
                return i==0?1:-1;
            }

            /// Returns the cross section of the subprocess with the given id.

            virtual MC_integral<value_type> process_xsec(int proc_id) const
            {
                const event_generator_base<model_t,N_in,N_out>* sub_gen=get_sub_generator(proc_id);
                return sub_gen==NULL?MC_integral<value_type>():(sub_gen->process_xsec(1));
            }

            /// Returns the maximal weight of the subprocess with the given id.

            virtual value_type process_maxw(int proc_id) const
            {
                const event_generator_base<model_t,N_in,N_out>* sub_gen=get_sub_generator(proc_id);
                return sub_gen==NULL?(value_type)0:(sub_gen->process_maxw(1));
            }

        protected:

            /// Returns the generator with the given id.

            virtual const event_generator_base<model_t,N_in,N_out>* get_sub_generator(int proc_id) const
            {
                return NULL;
            }

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
