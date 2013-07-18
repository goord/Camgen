//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DECOMP_PROC_H_
#define CAMGEN_DECOMP_PROC_H_

#include <cstring>
#include <iostream>
#include <Camgen/debug.h>
#include <Camgen/vector.h>
#include <Camgen/logstream.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility functions for the decomposition of process strings with particle  *
 * families into a list of subprocesses.                                     *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<std::size_t N_in,std::size_t N_out>bool decompose_process(const std::string& proc,vector<std::string,N_in+N_out>& vec)
    {
	typename std::string::size_type counter=0;
	bool in=true;
	typename std::string::const_iterator it=proc.begin();
	while(it != proc.end())
	{
	    if(counter==(N_in+N_out))
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"detected more than "<<N_out<<" final state particles for "<<proc<<endlog;
		return true;
	    }
	    if(*it==' ')
	    {
		++it;
		continue;
	    }
	    if(*it==',')
	    {
		if(counter==(N_in-1) and in)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"detected more than "<<N_in<<" initial state particles for "<<proc<<endlog;
		    it=std::find(it,proc.end(),'>');
		}
		else
		{
		    ++counter;
		    ++it;
		}
		continue;
	    }
	    if(*it=='>')
	    {
		if(counter!=(N_in-1))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"detected less than "<<N_in<<" initial state particles for "<<proc<<endlog;
		    return false;
		}
		++counter;
		++it;
		in=false;
		continue;
	    }
	    vec[counter].push_back(*it);
	    ++it;
	}
	if(counter<(N_in+N_out-1))
	{
	    log(log_level::warning)<<"detected less than "<<N_out<<" final state particles for "<<proc<<endlog;
	    return false;
	}
	return true;
    }
}

#endif /*CAMGEN_DECOMP_PROC_H_*/

