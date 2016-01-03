//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/plt_strm.h>
#include <Camgen/plt_config.h>

namespace Camgen
{
    /* Plot stream abstract base class constructor: */

    plot_stream::plot_stream(){}
    
    /* Plot stream abstract base class destructor: */
    
    plot_stream::~plot_stream(){}

    /* Function plot stream constructor with gnuplot function character array
     * argument: */

    function_stream::function_stream(const std::string& function_):function(function_){}

    /* Function stream destructor: */

    function_stream::~function_stream(){}

    /* Function stream plot command writing method: */

    void function_stream::write(std::ostream& os,bool strip_paths) const
    {
	os<<function;
	if(style.size()!=0)
	{
	    os<<" with "<<style;
	}
	if(title.size()==0)
	{
	    os<<" notitle";
	}
	else
	{
	    os<<" title \""<<title<<"\"";
	}
    }

    /* Clone method implementation: */

    function_stream* function_stream::clone() const
    {
	return new function_stream(*this);
    }

    /* Static number of temporary data files: */
    
    int data_wrapper::tempfiles=0;

    /* Temporary-file mode data wrapper constructor with 2 leafs: */

    data_wrapper::data_wrapper(const void* var1,const void* var2):valid(true),tempfile(true),users(0),open(false)
    {
	data.resize(2);
	data[0]=static_cast<const double*>(var1);
	data[1]=static_cast<const double*>(var2);
        init();	
    }

    /* Temporary-file mode data wrapper constructor with 3 leafs: */

    data_wrapper::data_wrapper(const void* var1,const void* var2,const void* var3):valid(true),tempfile(true),users(0),open(false)
    {
	data.resize(3);
	data[0]=static_cast<const double*>(var1);
	data[1]=static_cast<const double*>(var2);
	data[2]=static_cast<const double*>(var3);
        init();	
    }

    /* Permanent-file mode data wrapper constructor with 2 leafs: */

    data_wrapper::data_wrapper(const std::string& fname_,const void* var1,const void* var2):valid(true),tempfile(false),fname(fname_),users(0),open(false)
    {
	data.resize(2);
	data[0]=static_cast<const double*>(var1);
	data[1]=static_cast<const double*>(var2);
        init();	
    }

    /* Permanent-file mode data wrapper constructor with 3 leafs: */

    data_wrapper::data_wrapper(const std::string& fname_,const void* var1,const void* var2,const void* var3):valid(true),tempfile(false),fname(fname_),users(0),open(false)
    {
	data.resize(3);
	data[0]=static_cast<const double*>(var1);
	data[1]=static_cast<const double*>(var2);
	data[2]=static_cast<const double*>(var3);
        init();	
    }

    void data_wrapper::add_leaf(const void* var)
    {
	data.push_back(static_cast<const double*>(var));
    }

    /* Data buffer destructor: */

    data_wrapper::~data_wrapper()
    {
	if(open)
	{
	    write();
	}
	if(valid)
	{
	    if(tempfile)
	    {
		remove(fname.c_str());
		--tempfiles;
	    }
	}
    }

    /* Filling method: */

    bool data_wrapper::fill()
    {
	if(valid and open)
	{
	    for(std::vector<const double*>::size_type i=0;i<data.size();++i)
	    {
		datastream<<*(data[i])<<"\t";
	    }
	    datastream<<std::endl;
	}
	return (valid and open);
    }

    /* Writing method: */

    bool data_wrapper::write()
    {
	if(valid and open)
	{
	    datastream.close();
	    open=false;
	    return true;
	}
	return false;
    }

    /* Data file output: */

    const std::string& data_wrapper::filename() const
    {
	return fname;
    }

    /* Initialisation utility method: */

    void data_wrapper::init()
    {
	if(tempfile)
	{
	    char* tmpfile=plot_config::tmp_file();
	    if(tmpfile==NULL)
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"error generating temporary data file name--no data will be written"<<endlog;
		valid=false;
		return;
	    }
	    if(tempfiles==plot_config::max_tmp_files()-1)
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"maximum number of temporary data files reached--no data will be written."<<endlog;
		valid=false;
		return;
	    }
	    fname=std::string(tmpfile);
	    delete[] tmpfile;
	}
	datastream.open(fname.c_str());
	if(!datastream.is_open())
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"error opening temporary data stream--no data will be written."<<endlog;
	    valid=false;
	    return;
	}
	open=true;
	++tempfiles;
    }

    /* Data stream constructor with data file: */

    data_stream::data_stream(data_wrapper* data_):command("using 1:2"),data(data_)
    {
	++(data->users);
    }

    /* Constructor with data file name and column function string: */

    data_stream::data_stream(data_wrapper* data_,const std::string& col2):command("using 1:"+col2),data(data_)
    {
	++(data->users);
    }

    /* Constructor with data file name and column function string: */

    data_stream::data_stream(data_wrapper* data_,const std::string& col1,
	    					 const std::string& col2):command("using "+col1+":"+col2),data(data_)
    {
	++(data->users);
    }

    /* Constructor with data file name and column function string: */

    data_stream::data_stream(data_wrapper* data_,const std::string& col1,
	    					 const std::string& col2,
						 const std::string& col3):command("using "+col1+":"+col2+":"+col3),data(data_)
    {
	++(data->users);
    }

    /* Copy constructor: */

    data_stream::data_stream(const data_stream& other):command(other.command),title(other.title),style(other.style),data(other.data)
    {
	++(data->users);
    }

    /* Destructor implementation: */

    data_stream::~data_stream()
    {
	--(data->users);
	if(data->users==0)
	{
	    delete data;
	}
    }
    
    /* Plot command writing method implementation: */

    void data_stream::write(std::ostream& os,bool strip_paths) const
    {
	std::string fpath=strip_paths?file_utils::get_file(data->filename()):data->filename();
	os<<"\""<<fpath<<"\" "<<command;
	if(style.size()!=0)
	{
	    os<<" with "<<style;
	}
	if(title.size()==0)
	{
	    os<<" notitle";
	}
	else
	{
	    os<<" title \""<<title<<"\"";
	}
    }

    /* Clone method implementation: */

    data_stream* data_stream::clone() const
    {
	return new data_stream(*this);
    }
}

