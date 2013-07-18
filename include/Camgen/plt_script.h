//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file plt_script.h
    \brief gnuplot-script producing class.
 */

#ifndef CAMGEN_PLT_SCRIPT_H_
#define CAMGEN_PLT_SCRIPT_H_

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <Camgen/plt_obj.h>
#include <Camgen/plt_strm.h>

namespace Camgen
{
    /// Gnuplot-scripting class.

    class plot_script
    {
	public:

	    /// Filename.

	    std::string file;

	    /// Terminal.

	    const char* terminal;

	    /// Plot title.

	    std::string title;

	    /// X-axis label.
	    
	    std::string xlabel;

	    /// X-range.

	    double xmin,xmax;

	    /// Logarithmic x-axis flag.

	    bool xlog;
	    
	    /// Y-axis label.
	    
	    std::string ylabel;

	    /// Y-range.

	    double ymin,ymax;
	    
	    /// Logarithmic y-axis flag.
	    
	    bool ylog;

	    /// Histogram style flag.

	    bool histo;
	    
	    /// Autoscale y-axis flag.
	    
	    bool autoscale;

	    /// Tics on the x-axis.

	    std::vector< std::pair<double,std::string> >xtics;

	    /// Tics on the y-axis.

	    std::vector< std::pair<double,std::string> >ytics;

	    /// Gridline flag.

	    bool grid;

	    /// Surface flag.
	    
	    bool splot;

	    /// Constructor with title argument.

	    plot_script(const std::string&,const char* term=NULL);

	    /// Constructor with title and x-ranges arguments.

	    plot_script(const std::string&,const double&,const double&,const char* term=NULL);
	    
	    /// Constructor with title and x- and y-range arguments.
	    
	    plot_script(const std::string&,const double&,const double&,const double&,const double&,const char* term=NULL);

	    /// Destructor.

	    ~plot_script();

	    /// Adds a plot object to the script.

	    void add_object(const plot_object*);
	    
	    /// Adds a plot command to the script.
	    
	    void add_plot(const plot_stream*);

	    /// Adds an x-tic with label.

	    void add_x_tic(const double&,const std::string&);

	    /// Adds an x-tic without label.

	    void add_x_tic(const double&);

	    /// Adds a y-tic with label.

	    void add_y_tic(const double&,const std::string&);

	    /// Adds a y-tic without label.

	    void add_y_tic(const double&);

	    /// Writes the script to the output stream.

	    std::ostream& write(std::ostream& os) const;

	    /// Resets settings for multiplots.

	    std::ostream& reset(std::ostream&) const;

	    /// Writes the script to an output file with extension ".gp".

	    bool write() const;

	    /// Invokes gnuplot to execute the plot command. If gnuplot cannot
	    /// be invoked, writes the plot script to disc.
	    
	    bool plot() const;

	    /// Returns the file extension for the terminal type argument.

	    static std::string file_extension(const std::string&);

	private:

	    /// Plot objects.
	    
	    std::vector<const plot_object*>objects;

	    /// Plot command streams.

	    std::vector<const plot_stream*>streams;

	    /// Graphics file extensions map.

	    static std::map<std::string,std::string> filetypes;

	    /// Initializer for file extensions map.

	    static void init_filetypes();
    };
}

#endif /*CAMGEN_PLT_SCRIPT_H_*/

