//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file process_split_if.h
  \brief Per-process file event generator output interface.
  */

#ifndef CAMGEN_PROCESS_SPLIT_IF_H_
#define CAMGEN_PROCESS_SPLIT_IF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * File process-split event generator output interface. Creates a process  *
 * generator with separate file for each subprocess.                       *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <map>
#include <Camgen/if_base.h>
#include <Camgen/if_engine.h>
#include <Camgen/if_output.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out>class process_split_interface: public interface_base<model_t,N_in,N_out>
    {
	typedef interface_base<model_t,N_in,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef event_generator_base<model_t,N_in,N_out> generator_type;
	    typedef typename base_type::event_type event_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;

	    struct sub_interface
	    {
		interface_output<model_t,N_in,N_out>* output;
		interface_engine<model_t,N_in,N_out>* engine;
	    };

	    /* Public constructors/destructors: */
	    /*----------------------------------*/

	    /// Constructor with process generator instance.

	    process_split_interface(generator_type* gen,interface_output<model_t,N_in,N_out>* output_,interface_engine<model_t,N_in,N_out>* engine_=NULL):output(output_),engine(engine_)
	    {
		std::string namebase(output->file_name);
		if(namebase.size()==0)
		{
		    namebase="proc";
		}
		namebase.append("_");
		int digits=0;
		size_type step=1;
		while(step<gen->processes())
		{
		    ++digits;
		    step*=10;
		}
		if(digits==0)
		{
		    digits=1;
		}
		for(size_type i=0;i<gen->processes();++i)
		{
		    size_type n=gen->process_id(i);
		    std::stringstream ss;
		    ss<<namebase<<std::setw(digits)<<std::setfill('0')<<n;
		    std::string fname(ss.str());
		    interface_output<model_t,N_in,N_out>* sub_output=output->create(fname);
		    sub_output->open_file();
		    interface_engine<model_t,N_in,N_out>* sub_engine=NULL;
		    if(engine!=NULL)
		    {
			sub_engine=engine->clone();
			sub_engine->add_branches(sub_output);
		    }
		    sub_output->branch(&(this->w),"weight");
		    sub_output->branch(&(this->proc_id),"proc_id");
		    sub_interface sub_int={sub_output,sub_engine};
		    sub_interfaces[n]=sub_int;
		}
	    }

	    /// Destructor.

	    ~process_split_interface()
	    {
		for(typename std::map<int,sub_interface>::iterator it=sub_interfaces.begin();it!=sub_interfaces.end();++it)
		{
		    if((it->second).engine!=NULL)
		    {
			delete (it->second).engine;
		    }
		    delete (it->second).output;
		}
		if(engine!=NULL)
		{
		    delete engine;
		}
		delete output;
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Writes the output file and closes this file. No more events can be
	    /// streamed anymore.

	    bool write()
	    {
		bool q=true;
		for(typename std::map<int,sub_interface>::iterator it=sub_interfaces.begin();it!=sub_interfaces.end();++it)
		{
		    q&=(it->second.output->close_file());
		}
		return q;
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Returns the event size.

	    size_type event_size() const
	    {
		typename std::map<int,sub_interface>::const_iterator it=sub_interfaces.begin();
		if(it==sub_interfaces.end())
		{
		    return 0;
		}
		interface_engine<model_t,N_in,N_out>* e=it->second.engine;
		if(e==NULL)
		{
		    return 0;
		}
		return (e->event_size()+sizeof(value_type));
	    }

	    /// Returns the size of the i-th output file.

	    size_type file_size(size_type i) const
	    {
		typename std::map<int,sub_interface>::const_iterator it=sub_interfaces.find(i);
		if(it==sub_interfaces.end())
		{
		    return 0;
		}
		interface_engine<model_t,N_in,N_out>* e=it->second.engine;
		if(e==NULL)
		{
		    return 0;
		}
		return (e->event_size()+sizeof(value_type))*(it->output_evts);
	    }

	    /// Prints statistics to a file.

	    bool write_statistics() const
	    {
		if(output==NULL)
		{
		    return false;
		}
		std::ofstream ofs;
		std::string fname(output->file_name+"_stats.dat");
		ofs.open(fname.c_str());
		if(ofs.is_open())
		{
		    this->print_statistics(ofs);
		    ofs.close();
		    return true;
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Failed to open file "<<fname<<endlog;
		}
		return false;
	    }

	protected:

	    /// Reads the current event.

	    bool fill_event(const event_type& evt)
	    {
		typename std::map<int,sub_interface>::iterator it=sub_interfaces.find(evt.process_id());
		if(it==sub_interfaces.end())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"process id "<<evt.process_id()<<" not found in map"<<endlog;
		    return false;
		}
		if(it->second.engine!=NULL)
		{
		    it->second.engine->fill(evt);
		}
		return it->second.output->write_event(evt);
	    }

	private:

	    /* Interface output object (unused): */

	    interface_output<model_t,N_in,N_out>* output;

	    /* Interface engine instance: */

	    interface_engine<model_t,N_in,N_out>* engine;

	    /* Collection of sub-process interfaces: */

	    std::map<int,sub_interface> sub_interfaces;
    };
}

#endif /*EVTS_IF_H_*/

