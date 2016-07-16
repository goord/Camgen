//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_output.h
  \brief Abstract base class for event output classes.
  */

#ifndef CAMGEN_EVT_OUTPUT_H_
#define CAMGEN_EVT_OUTPUT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for event output classes. Derived classes should  *
 * implement methods adding branches holding variables, opening closing  *
 * files and writing an event to disk.                                   *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <string>
#include <Camgen/event.h>

namespace Camgen
{
    /// Base class for output interface engines.

    template<class model_t,std::size_t N_in,std::size_t N_out>class event_output
    {
	public:

	    /* Type definitions: */

	    typedef event<model_t,N_in,N_out> event_type;
	    typedef typename event_type::size_type size_type;
	    typedef typename event_type::value_type value_type;
	    typedef typename event_type::momentum_type momentum_type;
	    typedef typename std::vector<value_type>::iterator weight_iterator;
	    typedef typename std::vector<value_type>::const_iterator const_weight_iterator;

	    /// Output file name.

	    const std::string file_name;

	    /// Output file description.

	    std::string description;

	    /// Constructor with file name argument.

	    event_output(const std::string& file_name_):file_name(file_name_){}

	    /// Constructor with file name and description arguments.

	    event_output(const std::string& file_name_,const std::string description_):file_name(file_name_),description(description_){}

	    /// Destructor.

	    virtual ~event_output(){}

	    /// Creation method, no data copying expected.

	    virtual event_output<model_t,N_in,N_out>* create(const std::string& file_name_) const=0;

	    /// Abstract method opening the file.

	    virtual bool open_file()=0;

	    /// Abstract method closing the file.

	    virtual bool close_file()=0;

	    /// Virtual method adding a branch holding a Lorentz vector.

	    virtual bool branch(const momentum_type*,const std::string&)
            {
                return true;
            }

	    /// Virtual method adding a branch holding a floating-point number.

	    virtual bool branch(const value_type*,const std::string&)
            {
                return true;
            }

	    /// Virtual method adding a branch holding an integer.

	    virtual bool branch(const int*,const std::string&)
            {
                return true;
            }

	    /// Virtual method adding a branch holding a boolean.

	    virtual bool branch(const bool*,const std::string&)
            {
                return true;
            }

	    /// Virtual method writing the event to the output file.

	    virtual bool write_event()
            {
                return true;
            }

	    /// Virtual method writing the event to the output file, with process
	    /// generator instance argument.

	    virtual bool write_event(const event_type& evt)
	    {
	    	return write_event();
	    }
    };
}

#endif /*CAMGEN_EVT_OUTPUT_H_*/

