//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <cstdio>
#include <sstream>
#include <locale>
#include <algorithm>
#include <Camgen/plt_config.h>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/plt_script.h>

namespace Camgen
{
    /* Static data initialization: */

    std::map<std::string,std::string> plot_script::filetypes;

    /* Constructor with filename argument: */

    plot_script::plot_script(const std::string& file_,const char* term):file(file_),terminal(term),xmin(-10),xmax(10),xlog(false),ymin(0),ymax(0),ylog(false),histo(false),autoscale(true),grid(false),splot(false){init_filetypes();}
    
    /* Constructor with filename and x-range arguments: */

    plot_script::plot_script(const std::string& file_,const double& xmin_,const double& xmax_,const char* term):file(file_),terminal(term),xmin(xmin_),xmax(xmax_),xlog(false),ymin(0),ymax(0),ylog(false),histo(false),autoscale(true),grid(false),splot(false){init_filetypes();}
    
    /* Constructor with filename, x-range and y-range arguments: */

    plot_script::plot_script(const std::string& file_,const double& xmin_,const double& xmax_,const double& ymin_,const double& ymax_,const char* term):file(file_),terminal(term),xmin(xmin_),xmax(xmax_),xlog(false),ymin(ymin_),ymax(ymax_),ylog(false),histo(false),autoscale(false),grid(false),splot(false){init_filetypes();}
    
    /* Destructor: */
    
    plot_script::~plot_script()
    {
	for(std::vector<const plot_stream*>::size_type i=0;i<streams.size();++i)
	{
	    delete streams[i];
	}
	for(std::vector<const plot_object*>::size_type i=0;i<objects.size();++i)
	{
	    delete objects[i];
	}
    }

    /* Adds a plot object: */

    void plot_script::add_object(const plot_object* pobj)
    {
	if(std::find(objects.begin(),objects.end(),pobj)==objects.end())
	{
	    objects.push_back(pobj);
	}
	else
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to add identical plotobject instance to plotscript ignored"<<endlog;
	}
    }

    /* Adds a plot command: */

    void plot_script::add_plot(const plot_stream* pstr)
    {
	if(std::find(streams.begin(),streams.end(),pstr)==streams.end())
	{
	    streams.push_back(pstr);
	}
	else
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to add identical plotstream instance to plotscript ignored"<<endlog;
	}
    }

    /* Adds an x-tic at position x with label label: */

    void plot_script::add_x_tic(const double& x,const std::string& label)
    {
	xtics.push_back(std::pair<double,std::string>(x,label));
    }

    /* Adds an x-tic at position x without label: */

    void plot_script::add_x_tic(const double& x)
    {
	std::ostringstream sstrm;
	sstrm<<x;
	xtics.push_back(std::pair<double,std::string>(x,sstrm.str()));
    }

    /* Adds a y-tic at position y with label label: */

    void plot_script::add_y_tic(const double& y,const std::string& label)
    {
	ytics.push_back(std::pair<double,std::string>(y,label));
    }

    /* Adds an y-tic at position y without label: */

    void plot_script::add_y_tic(const double& y)
    {
	std::ostringstream sstrm;
	sstrm<<y;
	xtics.push_back(std::pair<double,std::string>(y,sstrm.str()));
    }

    /* Writes the plot script to the output stream argument: */

    std::ostream& plot_script::write(std::ostream& os) const
    {
	if(terminal!=NULL)
	{
	    std::string term(terminal);
	    std::string::iterator it=term.begin();
	    while(std::isalpha(*it) and it!=term.end())
	    {
		++it;
	    }
	    std::string ftype(term.begin(),it);
	    std::map<std::string,std::string>::const_iterator it2=filetypes.find(ftype);
	    if(it2!=filetypes.end())
	    {
		os<<"set terminal "<<term<<std::endl;
		os<<"set output \""<<file<<'.'<<it2->second<<"\""<<std::endl;
	    }
	    else
	    {
		os<<"set terminal table"<<std::endl;
		os<<"set output \""<<file<<".dat"<<"\""<<std::endl;
	    }
	}
	if(title.size()!=0)
	{
	    os<<"set title \""<<title<<"\""<<std::endl;
	}
	if(xlabel.size()!=0)
	{
	    os<<"set xlabel \""<<xlabel<<"\""<<std::endl;
	}
	if(ylabel.size()!=0)
	{
	    os<<"set ylabel \""<<ylabel<<"\""<<std::endl;
	}
	if(histo)
	{
	    os<<"set style data histogram"<<std::endl;
	}
	if(xlog and !ylog)
	{
	    os<<"set logscale x"<<std::endl;
	}
	if(!xlog and ylog)
	{
	    os<<"set logscale y"<<std::endl;
	}
	if(xlog and ylog)
	{
	    os<<"set logscale xy"<<std::endl;
	}
	os<<"set xrange ["<<xmin<<":"<<xmax<<"]"<<std::endl;
	if(autoscale)
	{
	    os<<"set autoscale"<<std::endl;
	}
	else
	{
	    os<<"set yrange ["<<ymin<<":"<<ymax<<"]"<<std::endl;
	}
	if(xtics.size()!=0)
	{
	    os<<"set xtics (";
	    for(std::vector< std::pair<double,std::string> >::size_type i=0;i<xtics.size()-1;++i)
	    {
		os<<"\""<<xtics[i].second<<"\" "<<xtics[i].first<<",";
	    }
	    os<<"\""<<xtics.back().second<<"\" "<<xtics.back().first<<")"<<std::endl;
	}
	if(ytics.size()!=0)
	{
	    os<<"set ytics (";
	    for(std::vector< std::pair<double,std::string> >::size_type i=0;i<ytics.size()-1;++i)
	    {
		os<<"\""<<ytics[i].second<<"\" "<<ytics[i].first<<",";
	    }
	    os<<"\""<<ytics.back().second<<"\" "<<ytics.back().first<<")"<<std::endl;
	}
	if(grid)
	{
	    os<<"set grid"<<std::endl;
	}
	for(std::vector<const plot_object*>::size_type i=0;i<objects.size();++i)
	{
	    objects[i]->write(os,int(i+1));
	    os<<std::endl;
	}
	if(streams.size()!=0)
	{
	    if(splot)
	    {
		os<<"splot      ";
	    }
	    else
	    {
		os<<"plot      ";
	    }
	    streams[0]->write(os);
	    for(std::vector<const plot_stream*>::size_type i=1;i<streams.size();++i)
	    {
		os<<",	\\"<<std::endl<<"          ";
		streams[i]->write(os);
	    }
	}
	return os;
    }

    /* Resets log settings etc. for multiplots. */

    std::ostream& plot_script::reset(std::ostream& os) const
    {
	if(xlog and !ylog)
	{
	    os<<"unset logscale x"<<std::endl;
	}
	if(ylog and !xlog)
	{
	    os<<"unset logscale y"<<std::endl;
	}
	if(xlog and ylog)
	{
	    os<<"unset logscale xy"<<std::endl;
	}
	if(grid)
	{
	    os<<"unset grid"<<std::endl;
	}
	os<<"unset xtics"<<std::endl;
	os<<"set xtics"<<std::endl;
	os<<"unset ytics"<<std::endl;
	os<<"set ytics"<<std::endl;
	for(unsigned i=1;i<=objects.size();++i)
	{
	    os<<"unset object "<<i<<std::endl;
	}
	return os;
    }

    bool plot_script::write() const
    {
	std::ofstream fs;
	fs.open((file+".gp").c_str());
	if(!fs.is_open())
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to open new file "<<file<<".gp--ignoring call"<<endlog;
	    return false;
	}
	write(fs);
	fs.close();
	return true;
    }

    bool plot_script::plot() const
    {
#ifdef GNUPLOTPATH
	FILE* gnustream=plot_config::open_pipe();
	if(gnustream==NULL)
	{
	    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to open connection to gnuplot--proceeding by writing the script "<<file<<".gp..."<<endlog;
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
#else
	return write();
#endif
    }

    std::string plot_script::file_extension(const std::string& term_type)
    {
	if(filetypes.size()==0)
	{
	    init_filetypes();
	}
	if(filetypes.find(term_type)!=filetypes.end())
	{
	    return filetypes[term_type];
	}
	return "";
    }

    void plot_script::init_filetypes()
    {
	if(filetypes.size()==0)
	{
	    filetypes["postscript"]="eps";
	    filetypes["pdf"]="pdf";
	    filetypes["png"]="png";
	    filetypes["gif"]="gif";
	    filetypes["jpeg"]="jpeg";
	    filetypes["table"]="dat";
	}
    }
}

