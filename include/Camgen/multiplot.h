//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_MULTIPLOT_H_
#define CAMGEN_MULTIPLOT_H_

/*! \file multiplot.h
    \brief gnuplot multiplot wrapper class.
 */

#include <cstdio>
#include <sstream>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/plt_config.h>
#include <Camgen/plt_script.h>

namespace Camgen
{
    /// Multiplot-scripting class.

    template<std::size_t M,std::size_t N>class multi_plot
    {
	public:

	    /// Plot filename.

	    std::string file;

	    /// Terminal.

	    const char* terminal;

	    /// Multiplot title.
	    
	    std::string title;

	    /// Constructor.

	    multi_plot(const std::string& file_,const char* term=NULL):file(file_),terminal(term)
	    {
		for(std::size_t i=0;i<M;++i)
		{
		    for(std::size_t j=0;j<N;++j)
		    {
			plots[i][j]=NULL;
		    }
		}
	    }

	    /// Adds a plot script.

	    void add_plot(plot_script* plt,std::size_t i,std::size_t j)
	    {
		if(i<M and j<N)
		{
		    plots[i][j]=plt;
		}
	    }

	    /// Writes the multiplot script

	    std::ostream& write(std::ostream& os) const
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
		    if(file_ext.size()>0)
		    {
			os<<"set terminal "<<term<<std::endl;
			os<<"set output \""<<file<<'.'<<file_ext<<"\""<<std::endl;
		    }
		    else
		    {
			os<<"set terminal table"<<std::endl;
			os<<"set output \""<<file<<".dat"<<"\""<<std::endl;
		    }
		}
		os<<"set multiplot layout "<<M<<","<<N;
		if(title.size()!=0)
		{
		    os<<" title \""<<title<<"\"";
		}
		os<<std::endl;
		double origin[2]={0,0};
		double rowsize=(double)1/(double)M;
		double colsize=(double)1/(double)N;
		for(std::size_t m=0;m<M;++m)
		{
		    origin[0]=0;
		    for(std::size_t n=0;n<N;++n)
		    {
			if(plots[m][n]!=NULL)
			{
			    os<<"set size "<<colsize<<","<<rowsize<<std::endl;
			    os<<"set origin "<<origin[0]<<","<<origin[1]<<std::endl;
			    plots[m][n]->write(os);
			    os<<std::endl;
			    plots[m][n]->reset(os);
			    os<<std::endl;
			}
			origin[0]+=colsize;
		    }
		    origin[1]+=rowsize;
		}
		return os;
	    }

	    /// Writes the script to a file with extension ".gp".

	    bool write() const
	    {
		std::ofstream fs;
		fs.open((file+".gp").c_str());
		if(!fs.is_open())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to open new file "<<file<<".gp...ignoring call."<<endlog;
		    return false;
		}
		write(fs);
		fs.close();
		return true;
	    }

	    /// Invokes gnuplot to execute the plot script.

	    bool plot() const
	    {
		FILE* gnustream=plot_config::open_pipe();
		if(gnustream==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to open connection to gnuplot...proceeding by writing the script "<<file<<".gp"<<endlog;
		    return write();	    
		}
		std::ostringstream ss;
		write(ss);
		std::fputs(ss.str().c_str(),gnustream);
		if(plot_config::close_pipe(gnustream)==-1)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"closing connection to gnuplot failed"<<endlog;
		    return false;
		}
		return true;
	    }

	private:

	    plot_script* plots[M][N];
    };
}

#endif /*CAMGEN_MULTIPLOT_H_*/
