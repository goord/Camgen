//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/file_utils.h>
#include <algorithm>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
 #include <direct.h>
 #define MKDIR(dirname,mode) _mkdir(dirname)
 #define STATSTRUCT _stat
 #define STATFUNC(path,buf) _stat(path,buf)
 #define S_ISDIR(code) (((code) & S_IFMT) == S_IFDIR)
#elif defined(__APPLE__) || defined(__linux__) || defined(__unix__) 
 #include <unistd.h>
 #define MKDIR(dirname,mode) mkdir(dirname,mode)
 #define STATSTRUCT stat
 #define STATFUNC(path,buf) stat(path,buf) 
#else
 #error Unknown target platform
#endif

namespace Camgen
{
    std::string file_utils::get_file(const std::string& path)
    {
	typename std::string::const_reverse_iterator it=std::find_if(path.rbegin(),path.rend(),file_utils::match_directory_separator);
	if(it==path.rend())
	{
	    return path;
	}
	return std::string(it.base(),path.end());
    }

    std::string file_utils::get_file_name(const std::string& path)
    {
	std::string file=get_file(path);
	typename std::string::reverse_iterator it=std::find_if(file.rbegin(),file.rend(),file_utils::match_extension_separator);
	if(it==file.rend())
	{
	    return file;
	}
	return std::string(file.begin(),it.base());
    }

    std::string file_utils::get_directory(const std::string& path)
    {
	typename std::string::const_reverse_iterator it=std::find_if(path.rbegin(),path.rend(),file_utils::match_directory_separator);
	if(it==path.rend())
	{
	    return "";
	}
	return std::string(path.begin(),it.base());
    }

    int file_utils::create_directory(const std::string& path, int mode)
    {
	int status=0;
	typename std::string::const_iterator it1=path.begin();
	typename std::string::const_iterator it2;

	while(status==0 and (it2=std::find_if(it1,path.end(),file_utils::match_directory_separator))!=path.end())
	{
	    if(it2!=it1)
	    {
		status=file_utils::mkdir(std::string(path.begin(),it2),mode);
	    }
	    it1=++it2;
	}
	if(status==0)
	{
	    status=file_utils::mkdir(path,mode);
	}
	return status;
    }

    int file_utils::mkdir(const std::string& path,int mode)
    {
	struct STATSTRUCT status;
	int retval=0;
	const char* cp=path.c_str();
	if(STATFUNC(cp,&status)!=0)
	{
	    if(MKDIR(cp,mode)!=0 and errno!=EEXIST)
	    {
		retval=-1;
	    }
	}
	else if(!S_ISDIR(status.st_mode))
	{
	    errno=ENOTDIR;
	    retval=-1;
	}
	return retval;
    }

    bool file_utils::match_directory_separator(char c)
    {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	return c=='\\' or c=='/';
#elif defined(__APPLE__) || defined(__linux__) || defined(__unix__)
	return c=='/';
#else
 #error Unknown target platform
#endif
    }

    bool file_utils::match_extension_separator(char c)
    {
	return c=='.';
    }
}

