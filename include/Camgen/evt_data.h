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
                bids.assign(0);
                pdfgs.assign(-1);
                pdfss.assign(-1);
                muF=(value_type)0;
                cols.assign(0);
                anticols.assign(0);
            }

            /// Clone implementation.

            event<model_t,N_in,N_out>* clone() const
            {
                return new event_data<model_t,N_in,N_out>(*this);
            }

            /// Implementation of the weight accessor.

            value_type w() const
            {
                return weight;
            }

            /// Implementation of the maximal weight accessor.

            value_type max_w() const
            {
                return max_weight;
            }

            /// Implementation of the cross section accessor.

            MC_integral<value_type> xsec() const
            {
                return cross_section;
            }

            /// Implementation of the cross section accessor.

            MC_integral<value_type> process_xsec() const
            {
                return process_cross_section;
            }

            /// Returns the partonic invariant mass-squared.

            value_type s_hat() const
            {
                return ecm*ecm;
            }

            /// Returns the partonic invariant mass.

            value_type Ecm_hat() const
            {
                return ecm;
            }

            /// Returns the hadronic invariant mass.

            value_type E_beam(size_type i) const
            {
                return ebeam[i];
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

            /// Implementation of the beam id accessor.

            int beam_id(int i) const
            {
                return (i<0 and i>=-(int)N_in)?bids[-i-1]:0;
            }

            /// Implementation of the pdf group id accessor.

            int pdfg(int i) const
            {
                return (i<0 and i>=-(int)N_in)?pdfgs[-i-1]:-1;
            }

            /// Implementation of the pdf set id accessor.

            int pdfs(int i) const
            {
                return (i<0 and i>=-(int)N_in)?pdfss[-i-1]:-1;
            }

            /// Implementation of the factorization scale accessor.

            value_type mu_F() const
            {
                return muF;
            }

            /// Implementation of the incoming colour connection accessor.

            int c_in(size_type i) const
            {
                return (i<N_in)?cols[i]:0;
            }

            /// Implementation of the outgoing colour connection accessor.

            int c_out(size_type i) const
            {
                return (i<N_out)?cols[i+N_in]:0;
            }

            /// Implementation of the incoming anti-colour connection accessor.

            int cbar_in(size_type i) const
            {
                return (i<N_in)?anticols[i]:0;
            }

            /// Implementation of the outgoing anti-colour connection accessor.

            int cbar_out(size_type i) const
            {
                return (i<N_out)?anticols[i+N_in]:0;
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

            /// Sets the event weight.

            void set_max_w(const value_type& max_weight_)
            {
                max_weight=max_weight_;
            }

            /// Updates the total cross_section.

            void set_xsec(const MC_integral<value_type> cross_section_)
            {
                cross_section=cross_section_;
            }

            /// Updates the process cross_section.

            void set_process_xsec(const MC_integral<value_type> cross_section_)
            {
                process_cross_section=cross_section_;
            }

            /// Copies the argument vector to the i-th incoming momentum.

            void set_p_in(size_type i,const momentum_type& p)
            {
                pin[i]=p;
            }

            /// Copies the argument vector to the i-th outgoing momentum.

            void set_p_out(size_type i,const momentum_type& p)
            {
                pout[i]=p;
            }

            /// Sets the incoming beam energy.

            void set_beam_energy(int i,const value_type& e)
            {
                if(i<0 and i>=-(int)N_in)
                {
                    ebeam[-i-1]=e;
                }
            }

            /// Implementation of the beam id setter.

            void set_beam_id(int i,int id)
            {
                if(i<0 and i>=-(int)N_in)
                {
                    bids[-i-1]=id;
                }
            }

            /// Implementation of the pdf group id setter.

            void set_pdfg(int i,int id)
            {
                if(i<0 and i>=-(int)N_in)
                {
                    pdfgs[-i-1]=id;
                }
            }

            /// Implementation of the pdf set id setter.

            void set_pdfs(int i,int id)
            {
                if(i<0 and i>=-(int)N_in)
                {
                    pdfss[-i-1]=id;
                }
            }

            /// Implementation of the factorization scale setter.

            void set_mu_F(const value_type& mu)
            {
                muF=mu;
            }

            /// Implementation of the colour connection setter.

            void set_colour_connection(const vector<int,N_in+N_out>& c,const vector<int,N_in+N_out>& cbar)
            {
                cols=c;
                anticols=cbar;
            }

        private:

            /* Momentum collections: */

            vector<momentum_type,N_in> pin;
            vector<momentum_type,N_out> pout;

            /* Beam ids: */
            
            vector<int,N_in> bids;

            /* Cernlib pd group and set numbers: */

            vector<int,N_in> pdfgs,pdfss;

            /* Invariant mass: */

            value_type ecm;

            /* Collider invariant mass: */

            vector<value_type,N_in> ebeam;

            /* Event weight: */

            value_type weight;

            /* Maximal event weight: */

            value_type max_weight;

            /* Total cross_section: */

            MC_integral<value_type> cross_section;

            /* Process cross_section: */

            MC_integral<value_type> process_cross_section;

            /* Factorization scale: */

            value_type muF;

            /* Color connection: */

            vector<int,N_in+N_out> cols,anticols;
    };
}

#endif /*CAMGEN_EVT_DATA_H_*/

