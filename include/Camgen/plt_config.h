//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file plt_conf.h
    \brief Static configuration data for plotting in Camgen.
 */

#ifndef CAMGEN_PLT_CONFIG_H_
#define CAMGEN_PLT_CONFIG_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Plotting configuration class, containing information about temporary  *
 * data files and a path to gnuplot for testing output.                  *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdio>

namespace Camgen
{
    /// Plotting configuration data holder.

    class plot_config
    {
	public:

	    /// Path to gnuplot program (NULL of not existent).

	    static const char* gnuplot_path;

	    /// Maximum number of temporary data files.

	    static int max_tmp_files();

	    /// Creates a unique temporary file name (dynamically
	    //allocated--needs to be deleted after usage).

	    static char* tmp_file();
	    
	    /// Opens pipe to gnuplot, returns NULL if gnuplot is not installed.
	    
	    static FILE* open_pipe();

	    /// Closes the pipe to gnuplot.

	    static int close_pipe(FILE*);
    };
}

#endif /*CAMGEN_PLT_CONFIG_H_*/

