//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file log_strm.h
    \brief Log stream object for Camgen.
 */

#ifndef CAMGEN_LOGSTREAM_H_
#define CAMGEN_LOGSTREAM_H_

#include <string>
#include <iostream>
#include <fstream>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Basic log stream definition for Camgen. Redirects the stream to the standard  *
 * output stream and/or a log file stream. Two extern global log streams are     *
 * defined: one for the matrix element calculation and one for the Monte Carlo   *
 * component.                                                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Log level type enumerable.

    struct log_level
    {
	enum type
	{
	    message=0,
	    warning=1,
	    error=2,
	    abort=3
	};
    };

    /// Placeholder class for ending log messages.

    class end_msg
    {
	public:

	    end_msg();
    };

    /// Ending message placeholder object.
    
    extern end_msg endlog;

    /// Log stream class declaration with standard output and file streams.

    class logstream: public std::ostream
    {
	public:

	    template<class T>friend logstream& operator << (logstream&,const T&);
	    friend logstream& operator << (logstream&,const end_msg&);
	    friend logstream& operator << (logstream&,std::ostream& (*fp)(std::ostream&));

	    /// Level above which logging is enabled.
	    
	    log_level::type enable_level;

	    /// Level of the current log message.

	    log_level::type level;

	    /// Flag controlling whether the log messages are streamed to
	    /// standard output.

	    bool streaming;

	    /// Flag controlling whether the log messages are streamed to file
	    /// (if such a file is created).
	    
	    bool writing;
	    
	    /// Flag controlling whether the logging prompts at message streams.

	    bool prompting;

	    /// Maximum of warnings allowed before exiting. A negative value means
	    /// that there is no limit.

	    int max_warnings;

	    /// Variable to keep track of all the warnings issued.

	    int warnings;

	    /// Default constructor.

	    logstream();

	    /// Constructor with log file name argument.

	    logstream(const std::string&);

	    /// Constructor with log file name and input/output stream
	    /// arguments.

	    logstream(const std::string&,std::istream&,std::ostream&);

	    /// Destructor.

	    ~logstream();

	    /// Opens the log file.

	    logstream& open(const std::string&);

	    /// Closes the log file.

	    logstream& close();

	    /// Sets the log level.

	    logstream& operator()(log_level::type);

	    /// Assertion: if failing, logging is disabled until end_message is
	    /// called.

	    logstream& assert(bool);

	    /// Returns whether the current configuration enables logging.

	    bool logging() const;

	private:

	    // Ends the log message and prompts (if prompting is set).

	    logstream& end_message();

	    // Helper method providing the log file header.

	    void init_logfile();

	    // Input stream (for prompting).
	    
	    std::istream& is;

	    // Output stream.
	    
	    std::ostream& os;

	    // Output file stream.

	    std::ofstream fs;

	    // Custom boolean for assertions.

	    bool enabled;
    };

    /// Logstream instance used in Camgen.

    extern logstream log;

    // Overloaded output stream operator.

    template<class T>logstream& operator << (logstream& ls,const T& message)
    {
	if(!ls.logging())
	{
	    return ls;
	}
	if(ls.streaming)
	{
	    ls.os<<message;
	}
	if(ls.writing and ls.fs.is_open())
	{
	    ls.fs<<message;
	}
	return ls;
    }

    logstream& operator << (logstream& ls,std::ostream& (*fp)(std::ostream&));

    logstream& operator << (logstream& ls,const end_msg&);
}

#endif /*CAMGEN_LOGSTREAM_H_*/

