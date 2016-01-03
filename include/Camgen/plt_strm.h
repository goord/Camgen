//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file plt_strm.h
    \brief data stream class definitions for plot-scripting class.
 */

#ifndef CAMGEN_PLT_STRM_H_
#define CAMGEN_PLT_STRM_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <Camgen/file_utils.h>

namespace Camgen
{
    /// Abstract base type of the plot streaming classes.

    class plot_stream
    {
	public:

	    /// Trivial constructor.

	    plot_stream();
	    
	    /// Trivial destructor.
	    
	    virtual ~plot_stream();

	    /// Plot command writing method.

	    virtual void write(std::ostream& os,bool strip_paths=false) const=0;

	    /// Virtual clone method.

	    virtual plot_stream* clone() const=0;
    };

    /// Plot command streaming class for gnuplot functions.

    class function_stream: public plot_stream
    {
	public:

	    /// Function name.

	    const std::string function;

	    /// Title of the plotted variable.

	    std::string title;
	    
	    /// Style of the plotted variable.
	    
	    std::string style;

	    /// Constructor with function name argument.

	    function_stream(const std::string&);

	    /// Destructor.

	    ~function_stream();

	    /// Streams the gnuplot command.

	    void write(std::ostream& os,bool strip_paths=false) const;

	    /// Clone method implementation.

	    function_stream* clone() const;
    };

    /// Data output stream wrapper class.

    class data_wrapper
    {
	friend class data_stream;

	public:

	    /// Number of temporary data files currently open.

	    static int tempfiles;

	    /// Temporary datafile-mode constructor with 2 leafs.

	    data_wrapper(const void* var1,const void* var2);

	    /// Temporary datafile-mode constructor with 3 leafs.

	    data_wrapper(const void* var1,const void* var2,const void* var3);

	    /// Permanent datafile-mode constructor with 2 leafs.

	    data_wrapper(const std::string&,const void* var1,const void* var2);

	    /// Permanent datafile-mode constructor with 3 leafs.

	    data_wrapper(const std::string&,const void* var1,const void* var2,const void* var3);

	    /// Adds a leaf to the list.

	    void add_leaf(const void* var);

	    /// Destructor. Closes the data file, and removes the file if it was
	    /// temporary.

	    ~data_wrapper();

	    /// Adds a new line to the output file with the values pointed to by
	    /// the internal adresses.

	    bool fill();

	    /// Writes the data and closes the file.

	    bool write();

	    /// Returns the name of the datafile.

	    const std::string& filename() const;

	private:

	    /* Initialization utility function: */

	    void init();

	    /* Valid input flag: */

	    bool valid;

	    /* Temporary/permanent output file flag: */

	    bool tempfile;
	    
	    /* Data file name: */
	    
	    std::string fname;

	    /* Pointers to written to the data file: */

	    std::vector<const double*> data;
	    
	    /* Output stream instance: */
	    
	    std::ofstream datastream;

	    /* Users counter: */

	    unsigned users;

	    /* Open data file flag: */

	    bool open;
    };

    /// Plot command streaming class for data-plotting.

    class data_stream: public plot_stream
    {
	public:

	    /// Plot command.
	    
	    const std::string command;
	    
	    /// Title of the plotted variable.
	    
	    std::string title;

	    /// Gnuplot style of the plotted variable.

	    std::string style;

	    /// Constructor with data buffer.

	    data_stream(data_wrapper*);

	    /// Constructor with data buffer and column string.

	    data_stream(data_wrapper*,const std::string&);

	    /// Constructor with data buffer and 2 column strings.

	    data_stream(data_wrapper*,const std::string&,const std::string&);

	    /// Constructor with data buffer and 3 column strings.

	    data_stream(data_wrapper*,const std::string&,const std::string&,const std::string&);

	    /// Copy constructor.

	    data_stream(const data_stream&);
	    
	    /// Destructor.
	    
	    ~data_stream();

	    /// Streams the gnuplot command.

	    void write(std::ostream& os,bool strip_paths=false) const;

	    /// Clone method implementation.

	    data_stream* clone() const;

	private:

	    data_wrapper* data;
    };
}

#endif /*CAMGEN_PLT_STRM_H_*/

