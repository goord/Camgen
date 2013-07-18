//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file if_output.h
  \brief Abstract base class for interface output types.
  */

#ifndef CAMGEN_IF_OUTPUT_H_
#define CAMGEN_IF_OUTPUT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for interface output types. Derived classes should  *
 * implement methods adding branches holding variables, opening closing    *
 * files and writing an event to disk.                                     *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/ps_gen_base.h>

namespace Camgen
{
    /// Base class for output interface engines.

    template<class model_t>class interface_output
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_type::value_type value_type;
	    typedef vector<value_type,model_t::dimension> momentum_type;
	    typedef typename std::vector<value_type>::iterator weight_iterator;
	    typedef typename std::vector<value_type>::const_iterator const_weight_iterator;

	    /// Output file name.

	    const std::string file_name;

	    /// Output file description.

	    std::string description;

	    /// Constructor with file name argument.

	    interface_output(const std::string& file_name_):file_name(file_name_){}

	    /// Constructor with file name and description arguments.

	    interface_output(const std::string& file_name_,const std::string description_):file_name(file_name_),description(description_){}

	    /// Destructor.

	    virtual ~interface_output(){}

	    /// Creation method, no data copying expected.

	    virtual interface_output<model_t>* create(const std::string& file_name_) const=0;

	    /// Abstract method opening the file.

	    virtual bool open_file()=0;

	    /// Abstract method closing the file.

	    virtual bool close_file()=0;

	    /// Virtual method adding a branch holding a Lorentz vector.

	    virtual bool branch(const momentum_type*,const std::string&){}

	    /// Virtual method adding a branch holding a floating-point number.

	    virtual bool branch(const value_type*,const std::string&){}

	    /// Virtual method adding a branch holding an integer.

	    virtual bool branch(const int*,const std::string&){}

	    /// Virtual method adding a branch holding a boolean.

	    virtual bool branch(const bool*,const std::string&){}

	    /// Virtual method writing the event to the output file.

	    virtual bool write_event(){}

	    /// Virtual method writing the event to the output file, with process
	    /// generator instance argument.

	    virtual bool write_event(const ps_generator_base<model_t>* gen)
	    {
		return write_event();
	    }
    };
}

#endif /*CAMGEN_IF_OUTPUT_H_*/

