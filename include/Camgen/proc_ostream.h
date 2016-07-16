//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file proc_ostream.h
  \brief Per-process event output stream.
  */

#ifndef CAMGEN_PROC_OSTREAM_H_
#define CAMGEN_PROC_OSTREAM_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * File process-split event output stream. Creates a process *
 * generator with separate file for each subprocess.         *
 *                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <map>
#include <Camgen/evt_stream.h>
#include <Camgen/evt_output_conf.h>

namespace Camgen
{
    /// Per-process event output streaming class. Creates separate files for each sub-process encountered in the event
    /// stream.
    
    template<class model_t,std::size_t N_in,std::size_t N_out>class process_output_stream: public event_stream<model_t,N_in,N_out>
    {
	typedef event_stream<model_t,N_in,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::event_type event_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;

            /* Utility struct: */

	    struct proc_stream
	    {
		event_output<model_t,N_in,N_out>* output;
		event_output_configuration<model_t,N_in,N_out>* config;
	    };

	    /* Public constructors/destructors: */
	    /*----------------------------------*/

	    /// Constructor with output definintion instance and configuration object:

	    process_output_stream(event_output<model_t,N_in,N_out>* output_,event_output_configuration<model_t,N_in,N_out>* config_):output(output_),config(config_){}

	    /// Constructor with output definition. Uses the default configuration.

	    process_output_stream(event_output<model_t,N_in,N_out>* output_):output(output_),config(new standard_output_configuration<model_t,N_in,N_out>()){}

	    /// Destructor.

	    ~process_output_stream()
	    {
		for(typename std::map<int,proc_stream>::iterator it=sub_streams.begin();it!=sub_streams.end();++it)
		{
		    if((it->second).config!=NULL)
		    {
			delete (it->second).config;
		    }
		    delete (it->second).output;
		}
		if(config!=NULL)
		{
		    delete config;
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
		for(typename std::map<int,proc_stream>::iterator it=sub_streams.begin();it!=sub_streams.end();++it)
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
		typename std::map<int,proc_stream>::const_iterator it=sub_streams.begin();
		if(it==sub_streams.end())
		{
		    return 0;
		}
		event_output_configuration<model_t,N_in,N_out>* e=it->second.config;
		if(e==NULL)
		{
		    return 0;
		}
		return (e->event_size()+sizeof(value_type));
	    }

	    /// Returns the size of the i-th output file.

	    size_type file_size(size_type i) const
	    {
		typename std::map<int,proc_stream>::const_iterator it=sub_streams.find(i);
		if(it==sub_streams.end())
		{
		    return 0;
		}
		event_output_configuration<model_t,N_in,N_out>* e=it->second.config;
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
		typename std::map<int,proc_stream>::iterator it=sub_streams.find(evt.process_id());
		if(it==sub_streams.end())
		{
                    proc_stream ps=add_process(evt.process_id());
                    if(ps.config!=NULL)
                    {
                        ps.config->fill(evt);
                    }
                    return ps.output->write_event(evt);
		}
		if(it->second.config!=NULL)
		{
		    it->second.config->fill(evt);
		}
		return it->second.output->write_event(evt);
	    }

	private:

	    /* Interface output object (unused): */

	    event_output<model_t,N_in,N_out>* output;

	    /* Interface config instance: */

	    event_output_configuration<model_t,N_in,N_out>* config;

	    /* Collection of sub-process interfaces: */

	    std::map<int,proc_stream> sub_streams;

            /* Adds a sub-process stream: */

            proc_stream& add_process(size_type id)
            {
		std::string namebase(output->file_name);
		if(namebase.size()==0)
		{
		    namebase="proc";
		}
                std::stringstream ss;
                ss<<namebase<<'_'<<id;
                std::string fname(ss.str());
                
                event_output<model_t,N_in,N_out>* sub_output=output->create(fname);
                sub_output->open_file();
                event_output_configuration<model_t,N_in,N_out>* sub_config=NULL;
                if(config!=NULL)
                {
                    sub_config=config->clone();
                    sub_config->add_branches(sub_output);
                }
                sub_output->branch(&(this->w),"weight");
                sub_output->branch(&(this->proc_id),"proc_id");
                proc_stream sub_int={sub_output,sub_config};
                sub_streams[id]=sub_int;
                return sub_streams[id];
            }
    };
}

#endif /*EVTS_IF_H_*/

