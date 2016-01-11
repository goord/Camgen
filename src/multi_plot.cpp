//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <config.h>
#include <cstdio>
#include <sstream>
#include <Camgen/multi_plot.h>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/file_utils.h>
#include <Camgen/plt_config.h>

namespace Camgen
{
    multi_plot::multi_plot(std::size_t rows_,std::size_t cols_,const std::string& file_,const char* term):file(file_),terminal(term),rows(rows_),cols(cols_)
    {
	plots.resize(rows*cols,NULL);
    }

    void multi_plot::add_plot(plot_script* plt,std::size_t i,std::size_t j)
    {
	if(i<rows and j<cols)
	{
	    plots[i*cols+j]=plt;
	}
    }

    const plot_script* multi_plot::get_plot(std::size_t i,std::size_t j) const
    {
	if(i<rows and j<cols)
	{
	    return static_cast<const plot_script*>(plots[i*cols+j]);
	}
	return NULL;
    }

    std::ostream& multi_plot::write(std::ostream& os,bool strip_paths) const
    {
	if(terminal!=NULL)
	{
	    std::string term(terminal);
	    std::string::iterator it=term.begin();
	    while(std::isalpha(*it) and it!=term.end())
	    {
		++it;
	    }
	    std::string file_ext=plot_script::file_extension(std::string(term.begin(),it));
	    std::string outputfile=strip_paths?file_utils::get_file(file):file;
	    if(file_ext.size()>0)
	    {
		os<<"set terminal "<<term<<std::endl;
		os<<"set output \""<<outputfile<<'.'<<file_ext<<"\""<<std::endl;
	    }
	    else
	    {
		os<<"set terminal table"<<std::endl;
		os<<"set output \""<<outputfile<<".dat"<<"\""<<std::endl;
	    }
	}
	os<<"set multiplot layout "<<rows<<","<<cols;
	if(title.size()!=0)
	{
	    os<<" title \""<<title<<"\"";
	}
	os<<std::endl;
	double origin[2]={0,0};
	double rowsize=(double)1/(double)rows;
	double colsize=(double)1/(double)cols;
	for(std::size_t m=0;m<rows;++m)
	{
	    origin[0]=0;
	    for(std::size_t n=0;n<cols;++n)
	    {
		const plot_script* p=get_plot(m,n);
		if(p!=NULL)
		{
		    os<<"set size "<<colsize<<","<<rowsize<<std::endl;
		    os<<"set origin "<<origin[0]<<","<<origin[1]<<std::endl;
		    p->write(os,strip_paths);
		    os<<std::endl;
		    p->reset(os);
		    os<<std::endl;
		}
		origin[0]+=colsize;
	    }
	    origin[1]+=rowsize;
	}
	return os;
    }

    bool multi_plot::write() const
    {
	std::ofstream fs;
	fs.open((file+".gp").c_str());
	if(!fs.is_open())
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to open new file "<<file<<".gp...ignoring call."<<endlog;
	    return false;
	}
	write(fs,true);
	fs.close();
	return true;
    }

    bool multi_plot::plot() const
    {
#ifdef GNUPLOTPATH
	FILE* gnustream=plot_config::open_pipe();
	if(gnustream==NULL)
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to open connection to gnuplot...proceeding by writing the script "<<file<<".gp"<<endlog;
	    return write();	    
	}
	std::ostringstream ss;
	write(ss,false);
	std::fputs(ss.str().c_str(),gnustream);
	if(plot_config::close_pipe(gnustream)==-1)
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"closing connection to gnuplot failed"<<endlog;
	    return false;
	}
	return true;
#else
	return write();
#endif
    }
}

