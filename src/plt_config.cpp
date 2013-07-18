//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <cstdlib>
#include <limits>
#include <cstring>
#include <iostream>
#include <config.h>
#include <Camgen/plt_config.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
 #include <io.h>
 #define MAX_TMP_FILES  27
 #define MAKE_TMP_FILE(fname) _mktemp(fname)!=NULL
 #define FNAME_TEMPLATE "camgendataiXXXXXX"
 #define PIPE_OPEN _popen
 #define PIPE_CLOSE _pclose
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__) 
 #include <unistd.h>
 #define MAX_TMP_FILES  64
 #define FNAME_TEMPLATE "/tmp/camgendataiXXXXXX"
 #define MAKE_TMP_FILE(fname) mkstemp(fname)!=-1
 #define PIPE_OPEN popen
 #define PIPE_CLOSE pclose
#endif

namespace Camgen
{

#ifdef GNUPLOTPATH
    const char* plot_config::gnuplot_path=GNUPLOTPATH;
#else
    const char* plot_config::gnuplot_path=NULL;
#endif

    int plot_config::max_tmp_files()
    {
#ifdef MAX_TMP_FILES
	return MAX_TMP_FILES;
#else
	return std::numeric_limits<int>::max();
#endif
    }

    char* plot_config::tmp_file()
    {
#ifdef FNAME_TEMPLATE
#ifdef MAKE_TMP_FILE
	char* fname=new char[sizeof(FNAME_TEMPLATE)];
	std::strcpy(fname,FNAME_TEMPLATE);
	if(MAKE_TMP_FILE(fname))
	{
	    return fname;
	}
	else
	{
	    delete[] fname;
	    return NULL;
	}
#else
	return NULL;
#endif

#else
	return NULL;
#endif
    }

    FILE* plot_config::open_pipe()
    {
#ifdef GNUPLOTPATH
#ifdef PIPE_OPEN
	return PIPE_OPEN(GNUPLOTPATH,"w");
#endif
#endif
	return NULL;
    }

    int plot_config::close_pipe(FILE* file)
    {
#ifdef GNUPLOTPATH
#ifdef PIPE_CLOSE
	return PIPE_CLOSE(file);
#endif
#endif
	return -1;
    }
}

#ifdef MAX_TMP_FILE
#undef MAX_TMP_FILE
#endif
#ifdef FNAME_TEMPLATE
#undef FNAME_TEMPLATE
#endif
#ifdef MAKE_TMP_FILE
#undef MAKE_TMP_FILE
#endif
#ifdef PIPE_OPEN
#undef PIPE_OPEN
#endif
#ifdef PIPE_CLOSE
#undef PIPE_CLOSE
#endif

