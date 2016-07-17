//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file Herwig_if.h
    \brief Les-Houches interface output for Herwig++.
 */

#ifndef CAMGEN_HERWIG_IF_H_
#define CAMGEN_HERWIG_IF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Les-Houches-type output for event generators in camgen to Herwig++  *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/ps_interface.h>

/* External component, should be in the c-flags when compiling this code: */

#include <LesHouches/LesHouchesReader.h>

namespace Camgen
{
    /// Les-Houches event file output interface class.

    template<class model_t,std::size_t N_out>class Herwig_interface: public ThePEG::LesHouchesReader,
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

	    Herwig_interface(generator_type* gen_,int weight_switch_,bool generate_events_=true,unsigned proc_id_=1):base_type(gen_,weight_switch_,generate_events_),proc_id(proc_id_){}

	    /// Destructor.

	    ~Herwig_interface(){}

	    /// Initialisation method.

	    void open()
	    {
		(this->heprup).IDWTUP=this->weight_switch;

                const event_type* e=this->get_current_event();
		(this->heprup).IDBMUP.first=e->beam_id(-1);
		(this->heprup).IDBMUP.second=e->beam_id(-2);
		(this->heprup).PDFGUP.first=e->pdfg(-1);
		(this->heprup).PDFGUP.second=e->pdfg(-2);
		(this->heprup).PDFSUP.first=e->pdfs(-1);
		(this->heprup).PDFSUP.second=e->pdfs(-2);
		MC_integral<value_type>sigma=e->xsec();
		(this->heprup).resize(1);
		(this->heprup).LPRUP[0]=proc_id;
		(this->heprup).XSECUP[0]=(double)sigma.value;
		(this->heprup).XERRUP[0]=(double)sigma.error;
		(this->heprup).XMAXUP[0]=e->max_w();
	    }

	    /// Fills the event common block.

	    bool setEvent(int idProcess=0)
	    {
                const event_type* e=this->next_event();
                if(e==NULL)
                {
                    return false;
                }

		this->hepeup.resize(2+N_out);
		this->hepeup.IDPRUP=proc_id;
		value_type w=(this->weight_switch==3)?1:e->w();
		this->hepeup.XWGTUP=(double)w;
//
//		this->hepeup.XPDUP=...; ...not necessary (yet).
//				
		this->hepeup.SCALUP=(double)(e->mu_F());
		this->hepeup.AQEDUP=(double)(model_type::alpha);
		this->hepeup.AQCDUP=(double)(model_type::alpha_s);
		
		for(unsigned i=0;i<2;++i)
		{
		    this->hepeup.IDUP[i]=e->id_in(i);
		    this->hepeup.ISTUP[i]=-1;
		    this->hepeup.MOTHUP[i]=std::pair<int,int>(0,0);
		    this->hepeup.ICOLUP[i]=std::pair<int,int>(e->c(i),e->cbar(i));
		    this->hepeup.PUP[i][0]=(double)e->p_in(i,1);
		    this->hepeup.PUP[i][1]=(double)e->p_in(i,2);
		    this->hepeup.PUP[i][2]=(double)e->p_in(i,3);
		    this->hepeup.PUP[i][3]=(double)e->p_in(i,0);
		    this->hepeup.PUP[i][4]=(double)e->M_in(i);
		    this->hepeup.VTIMUP[0]=(double)0;
		    this->hepeup.SPINUP[0]=(double)9;
		}
		for(unsigned i=0;i<N_out;++i)
		{
		    unsigned j=i+2;
		    this->hepeup.IDUP[j]=e->id_out(i);
		    this->hepeup.ISTUP[j]=1;
		    this->hepeup.MOTHUP[j]=std::pair<int,int>(1,2);
		    this->hepeup.ICOLUP[j]=std::pair<int,int>(e->c(j),e->cbar(j));
		    this->hepeup.PUP[j][0]=(double)e->p_out(i,1);
		    this->hepeup.PUP[j][1]=(double)e->p_out(i,2);
		    this->hepeup.PUP[j][2]=(double)e->p_out(i,3);
		    this->hepeup.PUP[j][3]=(double)e->p_out(i,0);
		    this->hepeup.PUP[j][4]=(double)e->M_out(i);
		    this->hepeup.VTIMUP[0]=(double)0;
		    this->hepeup.SPINUP[0]=(double)9;
		}
		return true;
	    }
    };
}

#endif /*CAMGEN_HERWIG_IF_H_*/

