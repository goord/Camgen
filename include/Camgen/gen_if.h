//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file gen_if.h
  \brief Phase space/event generator output interface.
  */

#ifndef CAMGEN_GEN_IF_H_
#define CAMGEN_GEN_IF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Generic interface for event generators. Takes a generator base instance,  *
 * interface config instance and interface output instance as constructor    *
 * arguments.                                                                *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/if_base.h>
#include <Camgen/evt_output_conf.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out>class generator_interface: public interface_base<model_t,N_in,N_out>
    {
	typedef interface_base<model_t,N_in,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::event_type event_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;

	    /* Public constructors/destructors: */
	    /*----------------------------------*/

	    /// Constructor with process generator instance.

	    generator_interface(event_output<model_t,N_in,N_out>* output_,event_output_configuration<model_t,N_in,N_out>* config_=NULL):output(output_),config(config_)
            {
		output->open_file();
		if(config!=NULL)
		{
		    config->add_branches(output);
		}
		output->branch(&(this->w),"weight");
		output->branch(&(this->proc_id),"proc_id");
            }

	    /// Destructor.

	    ~generator_interface()
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

#endif /*CAMGEN_GEN_IF_H_*/

