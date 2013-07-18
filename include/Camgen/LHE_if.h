//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file LHE_if.h
    \brief Les-Houches interface output for phase space generators.
 */

#ifndef CAMGEN_LHE_IF_H_
#define CAMGEN_LHE_IF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Les-Houches event format output interface for process/event generators in *
 * Camgen.                                                                   *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <fstream>
#include <Camgen/if_output.h>

namespace Camgen
{
    /// Les-Houches event file output interface class.

    template<class model_t>class LHE_interface: public interface_base<model_t>
    {
	typedef interface_base<model_t> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef ps_generator_base<model_t> generator_type;
	    typedef typename model_type::value_type value_type;
	    typedef typename generator_type::size_type size_type;
	    typedef vector<value_type,model_t::dimension> momentum_type;

	    /// Output file name.

	    const std::string file_name;

	    /// Process id.

	    const unsigned proc_id;

	    /// Les-Houches accord weight switch: 
	    /// 1: weighted events to PS, unweighted out, no xsec provided.
	    /// 2: weighted events to PS, unweighted out, xsecs provided.
	    /// 3: unweighted events to PS, unweighted out.
	    /// 4: weighted events to PS, weighted out.

	    const int weight_switch;

	    /// Description string.

	    std::string description;

	    /// Constructor.

	    LHE_interface(generator_type* gen,const std::string& file_name_,int weight_switch_,unsigned proc_id_=1):base_type(gen),file_name(file_name_),proc_id(proc_id_),weight_switch(weight_switch_)
	    {
		open_file();
		write_init();
	    }

	    /// Constructor with description.

	    LHE_interface(generator_type* gen,const std::string& file_name_,int weight_switch_,const std::string& descr_,unsigned proc_id_=1):base_type(gen),file_name(file_name_),proc_id(proc_id_),weight_switch(weight_switch_),description(descr_)
	    {
		open_file();
		write_init();
	    }

	    /// Destructor.

	    ~LHE_interface(){}

	    /// Writes and closes the datafile.

	    bool write()
	    {
		if(ofs.is_open())
		{
		    ofs.close();
		}
		return !(ofs.is_open());
	    }

	    /// Fills the event common block.

	    bool fill_event()
	    {
		if(!ofs.is_open())
		{
		    return false;
		}
		std::vector<int>c,cbar;
		this->gen->fill_colours(c,cbar);
		size_type nin=this->gen->n_in();
		size_type npart=nin+this->gen->n_out();
		if(c.size()!=npart)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"colour vectors were incorrectly filled--automatic resizing performed"<<endlog;
		    c.resize(npart,0);
		    cbar.resize(npart,0);
		}
		ofs<<"<event>"<<std::endl;
		value_type w=(weight_switch==3)?1:this->gen->w();
		ofs<<npart<<"\t1\t"<<w<<"\t"<<this->gen->mu_F()<<"\t"<<model_t::alpha<<"\t"<<model_t::alpha_s<<std::endl;
		for(unsigned i=0;i<nin;++i)
		{
		    ofs<<std::setw(4)<<this->gen->id_in(i);
		    ofs<<std::setw(6)<<-1;
		    ofs<<std::setw(6)<<0;
		    ofs<<std::setw(6)<<0;
		    ofs<<std::setw(6)<<c[i];
		    ofs<<std::setw(6)<<cbar[i];
		    ofs<<std::setw(20)<<this->gen->p_in(i,1);
		    ofs<<std::setw(20)<<this->gen->p_in(i,2);
		    ofs<<std::setw(20)<<this->gen->p_in(i,3);
		    ofs<<std::setw(20)<<this->gen->p_in(i,0);
		    ofs<<std::setw(20)<<this->gen->m_in(i);
		    ofs<<std::setw(6)<<"0.";
		    ofs<<std::setw(6)<<"9.";
		    ofs<<std::endl;
		}
		for(unsigned i=0;i<this->gen->n_out();++i)
		{
		    ofs<<std::setw(4)<<this->gen->id_out(i);
		    ofs<<std::setw(6)<<1;
		    ofs<<std::setw(6)<<1;
		    ofs<<std::setw(6)<<2;
		    ofs<<std::setw(6)<<c[i+nin];
		    ofs<<std::setw(6)<<cbar[i+nin];
		    ofs<<std::setw(20)<<this->gen->p_out(i,1);
		    ofs<<std::setw(20)<<this->gen->p_out(i,2);
		    ofs<<std::setw(20)<<this->gen->p_out(i,3);
		    ofs<<std::setw(20)<<this->gen->p_out(i,0);
		    ofs<<std::setw(20)<<this->gen->m_out(i);
		    ofs<<std::setw(6)<<"0.";
		    ofs<<std::setw(6)<<"9.";
		    ofs<<std::endl;
		}
		ofs<<"</event>"<<std::endl;
		return true;
	    }

	    /// Required disk space per event.

	    size_type event_size() const
	    {
		return ((this->gen->n_in()+this->gen->n_out())*(5*sizeof(value_type)+6*sizeof(int))+4*sizeof(value_type));
	    }

	protected:

	    /// Opens the (temporary) datafile.

	    bool open_file()
	    {
		std::string fname=this->file_name+".LHE";
		ofs.open(fname.c_str());
		if(ofs.is_open())
		{
		    ofs.precision(10);
		    ofs.setf(std::ios::scientific,std::ios::floatfield);
		    return true;
		}
		return false;
	    }

	    /// Initialisation block writing.

	    void write_init()
	    {
		if(!ofs.is_open())
		{
		    return;
		}
		ofs<<"<LesHouchesEvents version=\"1.0\">"<<std::endl;
		ofs<<"<!--"<<std::endl;
		ofs<<"File written by Camgen "<<VERSION<<" on "<<__DATE__<<" at "<<__TIME__<<std::endl;
		ofs<<"-->"<<std::endl;
		if(description.size()>0)
		{
		    ofs<<"# "<<description<<std::endl;
		}
		ofs<<"<init>"<<std::endl;
		ofs<<"\t"<<this->gen->beam_id(-1)<<"\t"<<this->gen->beam_id(-2);
		ofs<<"\t"<<this->gen->beam_energy(-1)<<"\t"<<this->gen->beam_energy(-2);
		ofs<<"\t"<<this->gen->pdfg(-1)<<"\t"<<this->gen->pdfg(-2);
		ofs<<"\t"<<this->gen->pdfs(-1)<<"\t"<<this->gen->pdfs(-2)<<std::endl;
		int ws(weight_switch);
		if(std::abs(ws)>4)
		{
		    ws=-4;
		}
		MC_integral<value_type>sigma=this->gen->xsec();
		ofs<<"\t"<<ws<<"\t 1 \t"<<sigma.value<<"\t"<<sigma.error<<"\t"<<this->gen->max_w()<<"\t"<<proc_id<<std::endl;
		ofs<<"</init>"<<std::endl;
	    }

	private:

	    /* Output file stream: */

	    std::ofstream ofs;
    };
}

#endif /*CAMGEN_LHE_IF_H_*/

