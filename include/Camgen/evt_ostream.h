//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_ostream.h
  \brief Event output streaming class.
  */

#ifndef CAMGEN_EVT_OSTREAM_H_
#define CAMGEN_EVT_OSTREAM_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Event output stream implementation. Takes an output config instance and   * 
 * event output instance as constructor arguments and streams according to   *
 * joint behaviour.                                                          *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_stream.h>
#include <Camgen/evt_output_conf.h>

namespace Camgen
{
    /// Event output stream class. Based upon the event output class and output configuration class, writes event data
    /// to disk.
    
    template<class model_t,std::size_t N_in,std::size_t N_out>class event_output_stream: public event_stream<model_t,N_in,N_out>
    {
	typedef event_stream<model_t,N_in,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::event_type event_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;

	    /* Public constructors/destructors: */
	    /*----------------------------------*/

	    /// Constructor with event output and configuration instances.

	    event_output_stream(event_output<model_t,N_in,N_out>* output_,event_output_configuration<model_t,N_in,N_out>* config_):output(output_),config(config_)
            {
		output->open_file();
		if(config!=NULL)
		{
		    config->add_branches(output);
		}
		output->branch(&(this->w),"weight");
		output->branch(&(this->proc_id),"proc_id");
            }

	    /// Constructor with event output class. Uses the standard configuration.

	    event_output_stream(event_output<model_t,N_in,N_out>* output_):output(output_)
            {
                config=new standard_output_configuration<model_t,N_in,N_out>();
		output->branch(&(this->w),"weight");
		output->branch(&(this->proc_id),"proc_id");
            }

	    /// Destructor.

	    ~event_output_stream()
	    {
		if(config!=NULL)
		{
		    delete config;
		}
		delete output;
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Writes the output file and closes this file. No events can be
	    /// streamed anymore.

	    bool write()
	    {
		return output->close_file();
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Returns the event size.

	    size_type event_size() const
	    {
		if(config!=NULL)
		{
		    return (config->event_size()+sizeof(value_type));
		}
		return 0;
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
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to open file "<<fname<<endlog;
		}
		return false;
	    }

	protected:


	    /// Reads the event.

	    bool fill_event(const event_type& evt)
	    {
		if(config!=NULL)
		{
		    config->fill(evt);
		}
		return output->write_event(evt);
	    }

	private:

	    /* Interface output instance: */

	    event_output<model_t,N_in,N_out>* output;

	    /* Interface config instance: */

	    event_output_configuration<model_t,N_in,N_out>* config;
    };
}

#endif /*CAMGEN_EVT_OSTREAM_H_*/

