//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file Pythia_if.h
    \brief Les-Houches interface output for Pythia 8.
 */

#ifndef CAMGEN_PYTHIA_IF_H_
#define CAMGEN_PYTHIA_IF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Les-Houches-type output for event generators in camgen to Pythia 8  *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/ps_interface.h>

#if 1

#include <Pythia8/LesHouches.h>

namespace Camgen
{
    /// Les-Houches event file output interface class.

    template<class model_t,std::size_t N_out>class Pythia_interface: public Pythia8::LHAup,
                                                                     public parton_shower_interface<model_t,2,N_out>
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
            typedef parton_shower_interface<model_t,2,N_out> base_type;
            typedef typename base_type::generator_type generator_type;
            typedef typename base_type::event_type event_type;
	    typedef typename event_type::value_type value_type;
	    typedef typename event_type::size_type size_type;
	    typedef typename event_type::momentum_type momentum_type;

	    /// Process id.

	    const unsigned proc_id;

	    /// Constructor.

	    Pythia_interface(generator_type* gen_,int weight_switch_,bool generate_events_=true,unsigned proc_id_=1):base_type(gen_,weight_switch_,generate_events_),proc_id(proc_id_){}

	    /// Destructor.

	    ~Pythia_interface(){}

	    /// Initialisation method.

	    bool setInit()
	    {
		this->setStrategy(this->weight_switch);

                const event_type* e=this->get_current_event();
		this->setBeamA(e->beam_id(-1),e->beam_energy(-1),e->pdfg(-1),e->pdfs(-1));
		this->setBeamB(e->beam_id(-2),e->beam_energy(-2),e->pdfg(-2),e->pdfs(-2));
		MC_integral<value_type>sigma=e->xsec();
		this->addProcess(proc_id,sigma.value,sigma.error,e->max_w());
		return true;
	    }

	    /// Fills the event common block.

	    bool setEvent(int idProcess=0)
	    {
                const event_type* e=this->next_event();
                if(e==NULL)
                {
                    return false;
                }

		value_type w=(this->weight_switch==3)?1:e->w();
		this->setProcess(proc_id,w,(double)(e->mu_F()),(double)(model_t::alpha),(double)(model_t::alpha_s));

		for(unsigned i=0;i<2;++i)
		{
		    this->addParticle(e->id_in(i),-1,0,0,
                                      e->c(i),
                                      e->cbar(i),
                                      (double)(e->p_in(i,1)),
                                      (double)(e->p_in(i,2)),
                                      (double)(e->p_in(i,3)),
                                      (double)(e->p_in(i,0)),
                                      (double)(e->M_in(i)),
                                      (double)0,
                                      (double)9);
		}
		for(unsigned i=0;i<N_out;++i)
		{
		    this->addParticle(e->id_out(i),1,1,2,
                                      e->c(i+2),
                                      e->cbar(i+2),
                                      (double)(e->p_out(i,1)),
                                      (double)(e->p_out(i,2)),
                                      (double)(e->p_out(i,3)),
                                      (double)(e->p_out(i,0)),
                                      (double)(e->M_out(i)),
                                      (double)0,
                                      (double)9);
		}
                this->pop_event();
		return true;
	    }
    };
}

#else

namespace Camgen
{
    /// Les-Houches event file output interface class.

    template<class model_t,std::size_t N_out>class Pythia_interface: public parton_shower_interface<model_t,2,N_out>
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
            typedef parton_shower_interface<model_t,2,N_out> base_type;
            typedef typename base_type::generator_type generator_type;
            typedef typename base_type::event_type event_type;
	    typedef typename event_type::value_type value_type;
	    typedef typename event_type::size_type size_type;
	    typedef typename event_type::momentum_type momentum_type;

	    /// Process id.

	    const unsigned proc_id;

	    /// Constructor.

	    Pythia_interface(generator_type* gen_,int weight_switch_,bool generate_events_=true,unsigned proc_id_=1):base_type(gen_,weight_switch_,generate_events_),proc_id(proc_id_){}

	    /// Destructor.

	    ~Pythia_interface(){}

	    /// Initialisation method.

	    bool setInit()
	    {
                return false;
	    }

	    /// Fills the event common block.

	    bool setEvent(int idProcess=0)
	    {
                return false;
	    }
    };
}
#endif /*HAVE_PYTHIA8_H*/

#endif /*CAMGEN_PYTHIA_IF_H_*/

