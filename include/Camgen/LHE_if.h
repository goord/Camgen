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
#include <Camgen/if_base.h>
#include <Camgen/if_output.h>

namespace Camgen
{
    /// Les-Houches event file output interface class.

    template<class model_t,std::size_t N_in,std::size_t N_out>class LHE_interface: public interface_base<model_t,N_in,N_out>
    {
	typedef interface_base<model_t,N_in,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::event_type event_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;

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

	    LHE_interface(const std::string& file_name_,int weight_switch_,unsigned proc_id_=1):file_name(file_name_),proc_id(proc_id_),weight_switch(weight_switch_),initialised(false)
	    {
		open_file();
	    }

	    /// Constructor with description.

	    LHE_interface(const std::string& file_name_,int weight_switch_,const std::string& descr_,unsigned proc_id_=1):file_name(file_name_),proc_id(proc_id_),weight_switch(weight_switch_),description(descr_),initialised(false)
	    {
		open_file();
	    }

	    /// Destructor.

	    ~LHE_interface(){}

	    /// Writes and closes the datafile.

	    bool write()
	    {
		if(ofs.is_open())
		{
		    ofs.close();
                    initialised=false;
		}
		return !(ofs.is_open());
	    }

	    /// Required disk space per event.

	    size_type event_size() const
	    {
		return ((N_in+N_out)*(5*sizeof(value_type)+6*sizeof(int))+4*sizeof(value_type));
	    }

	protected:

	    /// Fills the event common block.

	    bool fill_event(const event_type& evt)
	    {
		if(!ofs.is_open())
		{
		    return false;
		}
                if(!initialised)
                {
                    write_init(evt);
                    initialised=true;
                }
		ofs<<"<event>"<<std::endl;
		value_type w=(weight_switch==3)?1:evt.w();
		ofs<<(N_in+N_out)<<"\t1\t"<<w<<"\t"<<evt.mu_F()<<"\t"<<model_t::alpha<<"\t"<<model_t::alpha_s<<std::endl;
		for(unsigned i=0;i<N_in;++i)
		{
		    ofs<<std::setw(4)<<evt.id_in(i);
		    ofs<<std::setw(6)<<-1;
		    ofs<<std::setw(6)<<0;
		    ofs<<std::setw(6)<<0;
		    ofs<<std::setw(6)<<evt.c_in(i);
		    ofs<<std::setw(6)<<evt.cbar_in(i);
		    ofs<<std::setw(20)<<evt.p_in(i,1);
		    ofs<<std::setw(20)<<evt.p_in(i,2);
		    ofs<<std::setw(20)<<evt.p_in(i,3);
		    ofs<<std::setw(20)<<evt.p_in(i,0);
		    ofs<<std::setw(20)<<evt.M_in(i);
		    ofs<<std::setw(6)<<"0.";
		    ofs<<std::setw(6)<<"9.";
		    ofs<<std::endl;
		}
		for(unsigned i=0;i<N_out;++i)
		{
		    ofs<<std::setw(4)<<evt.id_out(i);
		    ofs<<std::setw(6)<<1;
		    ofs<<std::setw(6)<<1;
		    ofs<<std::setw(6)<<2;
		    ofs<<std::setw(6)<<evt.c_out(i);
		    ofs<<std::setw(6)<<evt.cbar_out(i);
		    ofs<<std::setw(20)<<evt.p_out(i,1);
		    ofs<<std::setw(20)<<evt.p_out(i,2);
		    ofs<<std::setw(20)<<evt.p_out(i,3);
		    ofs<<std::setw(20)<<evt.p_out(i,0);
		    ofs<<std::setw(20)<<evt.M_out(i);
		    ofs<<std::setw(6)<<"0.";
		    ofs<<std::setw(6)<<"9.";
		    ofs<<std::endl;
		}
		ofs<<"</event>"<<std::endl;
		return true;
	    }

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

            //TODO: Read version from config.h 
	    void write_init(const event_type& evt)
	    {
		ofs<<"<LesHouchesEvents version=\"1.0\">"<<std::endl;
		ofs<<"<!--"<<std::endl;
		ofs<<"File written by Camgen 1.0 "<<" on "<<__DATE__<<" at "<<__TIME__<<std::endl;
		ofs<<"-->"<<std::endl;
		if(description.size()>0)
		{
		    ofs<<"# "<<description<<std::endl;
		}
		ofs<<"<init>"<<std::endl;
		ofs<<"\t"<<evt.beam_id(-1)<<"\t"<<evt.beam_id(-2);
		ofs<<"\t"<<evt.beam_energy(-1)<<"\t"<<evt.beam_energy(-2);
		ofs<<"\t"<<evt.pdfg(-1)<<"\t"<<evt.pdfg(-2);
		ofs<<"\t"<<evt.pdfs(-1)<<"\t"<<evt.pdfs(-2)<<std::endl;
		int ws(weight_switch);
		if(std::abs(ws)>4)
		{
		    ws=-4;
		}
		MC_integral<value_type>sigma=evt.xsec();
		ofs<<"\t"<<ws<<"\t 1 \t"<<sigma.value<<"\t"<<sigma.error<<"\t"<<evt.max_w()<<"\t"<<proc_id<<std::endl;
		ofs<<"</init>"<<std::endl;
	    }

	private:

	    /* Output file stream: */

	    std::ofstream ofs;

            /* Initialization flag: */

            bool initialised;
    };
}

#endif /*CAMGEN_LHE_IF_H_*/

