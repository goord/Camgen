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

#include <fstream>
#include <Camgen/if_output.h>

/* External component, should be in the c-flags when compiling this code: */

#include <LesHouches.h>

namespace Camgen
{
    /// Les-Houches event file output interface class.

    template<class model_t,std::size_t N_out>class Pythia_interface: public Pythia8::LHAup, 
                                                                     public interface_base<model_t,2,N_out>
    {
	typedef interface_base<model_t,2,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::event_type event_type;
	    typedef typename event_type::value_type value_type;
	    typedef typename event_type::size_type size_type;
	    typedef typename event_type::momentum_type momentum_type;

	    /// Process id.

	    const unsigned proc_id;

	    /// Les-Houches accord weight switch: 
	    /// 1: weighted events to PS, unweighted out, no xsec provided.
	    /// 2: weighted events to PS, unweighted out, xsecs provided.
	    /// 3: unweighted events to PS, unweighted out.
	    /// 4: weighted events to PS, weighted out.

	    const int weight_switch;

	    /// Generation switch. If set true, the matrix element generation is
	    /// triggered by the parton shower, otherwise, it has to be performed
	    /// by the user.

	    bool generation_switch;

	    /// Constructor.

	    Pythia_interface(const event_type& e,int weight_switch_,unsigned proc_id_=1):proc_id(proc_id_),weight_switch(weight_switch_),generation_switch(true)
            {
                this->set_event(e);
            }

	    /// Constructor with description.

	    Pythia_interface(const event_type& e,int weight_switch_,const std::string& descr_,unsigned proc_id_=1):proc_id(proc_id_),weight_switch(weight_switch_),generation_switch(true),gen(gen_){}

	    /// Destructor.

	    ~Pythia_interface(){}

	    /// Initialisation method.

	    bool setInit()
	    {
		if(gen->n_in()!=2)
		{
		    return false;
		}
		this->setStrategy(weight_switch);
		this->setBeamA(gen->beam_id(-1),gen->beam_energy(-1),gen->pdfg(-1),gen->pdfs(-1));
		this->setBeamB(gen->beam_id(-2),gen->beam_energy(-2),gen->pdfg(-2),gen->pdfs(-2));
		MC_integral<value_type>sigma=gen->xsec();
		this->addProcess(proc_id,sigma.value,sigma.error,gen->max_w());
		return true;
	    }

	    /// Fills the event common block.

	    bool setEvent(int idProcess=0)
	    {
		if(generation_switch)
		{
		    gen->next_event(weight_switch);
		}
		value_type w=(weight_switch==3)?1:this->gen->w();
		this->setProcess(proc_id,w,(double)(gen->mu_F()),(double)(model_t::alpha),(double)(model_t::alpha_s));
		std::vector<int>c,cbar;
		gen->fill_colours(c,cbar);
		size_type nin=gen->n_in();
		size_type npart=nin+gen->n_out();
		if(c.size()!=npart)
		{
		    log(log_level::warning)<<"Colour vectors were incorrectly filled--automatic resizing performed."<<endlog;
		    c.resize(npart,0);
		    cbar.resize(npart,0);
		}

		for(unsigned i=0;i<nin;++i)
		{
		    this->addParticle(gen->id_in(i),-1,0,0,c[i],cbar[i],(double)(gen->p_in(i,1)),(double)(gen->p_in(i,2)),(double)(gen->p_in(i,3)),(double)(gen->p_in(i,0)),(double)(gen->m_in(i)),(double)0,(double)9);
		}
		for(unsigned i=0;i<this->gen->n_out();++i)
		{
		    this->addParticle(gen->id_out(i),1,1,2,c[i+nin],cbar[i+nin],(double)(gen->p_out(i,1)),(double)(gen->p_out(i,2)),(double)(gen->p_out(i,3)),(double)(gen->p_out(i,0)),(double)(gen->m_out(i)),(double)0,(double)9);
		}
		return true;
	    }

	    /// Required disk space per event.

	    size_type event_size() const
	    {
		return ((gen->n_in()+gen->n_out())*(5*sizeof(value_type)+6*sizeof(int))+4*sizeof(value_type));
	    }

	private:

	    /* Matrix-element event generator: */

	    generator_type* gen;
    };
}

#endif /*CAMGEN_PYTHIA_IF_H_*/

