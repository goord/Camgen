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

#include <fstream>

/* External component, should be in the c-flags when compiling this code: */

#include <LesHouches/LesHouchesReader.h>

namespace Camgen
{
    /// Les-Houches event file output interface class.

    template<class model_t>class Herwig_interface: public ThePEG::LesHouchesReader
    {
	typedef interface_base<model_t> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef ps_generator_base<model_t> generator_type;
	    typedef typename model_type::value_type value_type;
	    typedef typename generator_type::size_type size_type;
	    typedef vector<value_type,model_t::dimension> momentum_type;

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

	    Herwig_interface(generator_type* gen_,int weight_switch_,unsigned proc_id_=1):proc_id(proc_id_),weight_switch(weight_switch_),generation_switch(true),gen(gen_){}

	    /// Constructor with description.

	    Herwig_interface(generator_type* gen_,int weight_switch_,const std::string& descr_,unsigned proc_id_=1):proc_id(proc_id_),weight_switch(weight_switch_),generation_switch(true),gen(gen_){}

	    /// Destructor.

	    ~Herwig_interface(){}

	    /// Initialisation method.

	    void open()
	    {
		if(gen->n_in()!=2)
		{
		    return;
		}
		(this->heprup).IDWTUP=weight_switch;
		(this->heprup).IDBMUP.first=gen->beam_id(-1);
		(this->heprup).IDBMUP.second=gen->beam_id(-2);
		(this->heprup).PDFGUP.first=gen->pdfg(-1);
		(this->heprup).PDFGUP.second=gen->pdfg(-2);
		(this->heprup).PDFSUP.first=gen->pdfs(-1);
		(this->heprup).PDFSUP.second=gen->pdfs(-2);
		MC_integral<value_type>sigma=gen->xsec();
		(this->heprup).resize(1);
		(this->heprup).LPRUP[0]=proc_id;
		(this->heprup).XSECUP[0]=(double)sigma.value;
		(this->heprup).XERRUP[0]=(double)sigma.error;
		(this->heprup).XMAXUP[0]=gen->max_w();
	    }

	    /// Fills the event common block.

	    bool setEvent(int idProcess=0)
	    {
		if(generation_switch)
		{
		    gen->next_event(weight_switch);
		}
		this->hepeup.resize(gen->n_in()+gen->n_out());
		this->hepeup.IDPRUP=proc_id;
		value_type w=(weight_switch==3)?1:this->gen->w();
		this->hepeup.XWGTUP=(double)w;
//
//		this->hepeup.XPDUP=...; ...not necessary (yet).
//				
		this->hepeup.SCALUP=(double)(gen->mu_F());
		this->hepeup.AQEDUP=(double)(model_type::alpha);
		this->hepeup.AQCDUP=(double)(model_type::alpha_s);
		
		std::vector<int>c,cbar;
		gen->fill_colours(c,cbar);
		size_type nin=gen->n_in();
		size_type npart=nin+gen->n_out();
		if(c.size()!=npart)
		{
		    LOG(warning)<<"Colour vectors were incorrectly filled...automatic resizing performed."<<ENDLOG(warning);
		    c.resize(npart,0);
		    cbar.resize(npart,0);
		}
		for(unsigned i=0;i<nin;++i)
		{
		    this->hepeup.IDUP[i]=gen->id_in(i);
		    this->hepeup.ISTUP[i]=-1;
		    this->hepeup.MOTHUP[i]=std::pair<int,int>(0,0);
		    this->hepeup.ICOLUP[i]=std::pair<int,int>(c[i],cbar[i]);
		    this->hepeup.PUP[i][0]=(double)gen->p_in(i,1);
		    this->hepeup.PUP[i][1]=(double)gen->p_in(i,2);
		    this->hepeup.PUP[i][2]=(double)gen->p_in(i,3);
		    this->hepeup.PUP[i][3]=(double)gen->p_in(i,0);
		    this->hepeup.PUP[i][4]=(double)gen->m_in(i);
		    this->hepeup.VTIMUP[0]=(double)0;
		    this->hepeup.SPINUP[0]=(double)9;
		}
		for(unsigned i=0;i<this->gen->n_out();++i)
		{
		    unsigned j=i+nin;
		    this->hepeup.IDUP[j]=gen->id_out(i);
		    this->hepeup.ISTUP[j]=1;
		    this->hepeup.MOTHUP[j]=std::pair<int,int>(1,2);
		    this->hepeup.ICOLUP[j]=std::pair<int,int>(c[j],cbar[j]);
		    this->hepeup.PUP[j][0]=(double)gen->p_out(i,1);
		    this->hepeup.PUP[j][1]=(double)gen->p_out(i,2);
		    this->hepeup.PUP[j][2]=(double)gen->p_out(i,3);
		    this->hepeup.PUP[j][3]=(double)gen->p_out(i,0);
		    this->hepeup.PUP[j][4]=(double)gen->m_out(i);
		    this->hepeup.VTIMUP[0]=(double)0;
		    this->hepeup.SPINUP[0]=(double)9;
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

#endif /*CAMGEN_HERWIG_IF_H_*/

