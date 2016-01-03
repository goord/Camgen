//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_MULTI_PLOT_H_
#define CAMGEN_MULTI_PLOT_H_

/*! \file multi_plot.h
    \brief gnuplot multiplot wrapper class.
 */

#include <vector>
#include <iostream>
#include <Camgen/plt_script.h>

namespace Camgen
{
    /// Multiplot-scripting class.

    class multi_plot
    {
	public:

	    /// Plot filename.

	    std::string file;

	    /// Terminal.

	    const char* terminal;

	    /// Number of rows.
	    
	    const std::size_t rows;

	    /// Number of columns.

	    const std::size_t cols;

	    /// Multiplot title.
	    
	    std::string title;

	    /// Constructor.

	    multi_plot(std::size_t rows_,std::size_t cols_,const std::string& file_,const char* term=NULL);

	    /// Adds a plot script.

	    void add_plot(plot_script*,std::size_t,std::size_t);

	    /// Returns a plot script.

	    const plot_script* get_plot(std::size_t,std::size_t) const;

	    /// Writes the multiplot script

	    std::ostream& write(std::ostream& os,bool strip_paths=false) const;

	    /// Writes the script to a file with extension ".gp".

	    bool write() const;

	    /// Invokes gnuplot to execute the plot script.

	    bool plot() const;

	private:

	    std::vector<plot_script*> plots;

	    plot_script*& get_plot_script(std::size_t,std::size_t); 
    };
}

#endif /*CAMGEN_MULTI_PLOT_H_*/

