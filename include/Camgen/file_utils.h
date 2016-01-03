//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file file_utils.h
    \brief filesystem utilities
 */

#ifndef CAMGEN_FILE_UTILS_H_
#define CAMGEN_FILE_UTILS_H_ 

#include <string>

namespace Camgen
{
    /// Class holding filesystem utility functions
    
    class file_utils
    {
	public:

	    /// Returns the file, including extension, that the path refers to.

	    static std::string get_file(const std::string& path);

	    /// Returns the file name, without extension, that the path refers to.

	    static std::string get_file_name(const std::string& path);

	    /// Returns the full path of the directory.

	    static std::string get_directory(const std::string& path);

	    /// Creates a directory with the given name and returns its full path.

	    static int create_directory(const std::string& name,int mode=0777);

	private:

	    static int mkdir(const std::string& s,int m);

	    static bool match_directory_separator(char c);

	    static bool match_extension_separator(char c);

    };
}

#endif /*CAMGEN_FILE_UTILS_H_*/
