//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <iostream>
#include <Camgen/plt_config.h>
#include <Camgen/file_utils.h>
#include <Camgen/Minkowski.h>
#include <Camgen/stdrand.h>
#include <Camgen/histogram.h>
#include <Camgen/norm_gen.h>
#include <Camgen/uni_sphere.h>
#include <Camgen/val_gen_fac.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Testing facility for momentum channel sampling algorithms.                     *
*                                                                                *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing channel generation algorithms...................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    
    typedef double value_type;
    typedef std::size_t size_type;
    typedef random_number_stream<value_type,std::random> rn_stream;

    ///////////////////////////////////////////////////////////////

    /* Initial power value for power-law generation: */

    value_type nu(0.5);
    
    /* Initial mass value for power-law/Breit-Wigner generation: */
    
    value_type m(0);

    /* Initial width value for Breit-Wigner generation: */

    value_type w(10);
    
    /* Initial minimal generated invariant mass: */
    
    value_type smin(-5);
    
    /* Initial maximal generated invariant mass: */
    
    value_type smax(10);
    
    /* Initial number of Monte Carlo throws: */
    
    unsigned N_events=50000;
    
    /* Initial number of bins written to tape: */
    
    unsigned N_bins=500;
    
    ///////////////////////////////////////////////////////////////
    
    bool have_gp=plot_config::gnuplot_path!=NULL;
    file_utils::create_directory("test_output/MC_gen_test");
    
    {
	std::cerr<<"Checking normally distributed numbers.....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/Gauss";
	normal_generator<value_type,std::random>gen;
	value_type rho;
	value_type p;
	histogram<value_type>hist(&rho,&p,N_events);
	for(size_type i=0;i<N_events;++i)
	{
	    rho=gen.generate();
	    p=gen.weight();
	    hist.store();
	}
	hist.make(N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	value_type norm=N_events*hist.binsize();
	plt->add_plot(gen.pdf_plot(norm));
	plt->plot();
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    {
	std::cerr<<"Checking uniform circle generation.....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/circle";
	vector<value_type,2>pt;
	uniform_sphere<value_type,1,std::random>gen(&pt);
	data_wrapper* data=new data_wrapper(&pt[0],&pt[1]);
	for(size_type n=0;n<N_events/100;++n)
	{
	    gen.generate();
	    data->fill();
	}
	data->write();
	data_stream* dataplt=new data_stream(data);
	dataplt->style="points";
	plot_script* plt=new plot_script(filename,"postscript color");
	plt->add_plot(dataplt);
	plt->plot();
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    {
	std::cerr<<"Checking uniform sphere generation..........";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/sphere";
	vector<value_type,3>pt;
	uniform_sphere<value_type,2,std::random>gen(&pt);
	data_wrapper* data=new data_wrapper(&pt[0],&pt[1],&pt[2]);
	for(size_type n=0;n<N_events/10;++n)
	{
	    gen.generate();
	    data->fill();
	}
	data->write();
	data_stream* dataplt=new data_stream(data,"1","2","3");
	dataplt->style="points";
	plot_script* plt=new plot_script(filename,"postscript color");
	plt->splot=true;
	plt->add_plot(dataplt);
	plt->plot();
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    {
	std::cerr<<"Checking uniform 3-sphere generation.....";
	std::cerr.flush();
	vector<value_type,4>pt;
	uniform_sphere<value_type,3,std::random>gen(&pt);
	for(size_type n=0;n<N_events;++n)
	{
	    gen.generate();
	    value_type radius(0);
	    for(size_type i=0;i<pt.size();++i)
	    {
		radius+=(pt[i]*pt[i]);
	    }
	    if(!equals(radius,(value_type)1))
	    {
		Camgen::log(log_level::warning)<<"uniform sphere point has incorrect length "<<radius<<endlog;
	    }
	}
	std::cerr<<".....done."<<std::endl;
    }

    {
	std::cerr<<"Checking uniform 4-sphere generation.....";
	std::cerr.flush();
	vector<value_type,5>pt;
	uniform_sphere<value_type,4,std::random>gen(&pt);
	for(size_type n=0;n<N_events;++n)
	{
	    gen.generate();
	    value_type radius(0);
	    for(size_type i=0;i<pt.size();++i)
	    {
		radius+=(pt[i]*pt[i]);
	    }
	    if(!equals(radius,(value_type)1))
	    {
		Camgen::log(log_level::warning)<<"uniform sphere point has incorrect length "<<radius<<endlog;
	    }
	}
	std::cerr<<".....done."<<std::endl;
    }

    {
	std::cerr<<"Checking uniform 5-sphere generation.....";
	std::cerr.flush();
	vector<value_type,6>pt;
	uniform_sphere<value_type,5,std::random>gen(&pt);
	for(size_type n=0;n<N_events;++n)
	{
	    gen.generate();
	    value_type radius(0);
	    for(size_type i=0;i<pt.size();++i)
	    {
		radius+=(pt[i]*pt[i]);
	    }
	    if(!equals(radius,(value_type)1))
	    {
		Camgen::log(log_level::warning)<<"uniform sphere point has incorrect length "<<radius<<endlog;
	    }
	}
	std::cerr<<".....done."<<std::endl;
    }

    {
	std::cerr<<"Checking uniform 6-sphere generation.....";
	std::cerr.flush();
	vector<value_type,7>pt;
	uniform_sphere<value_type,6,std::random>gen(&pt);
	for(size_type n=0;n<N_events;++n)
	{
	    gen.generate();
	    value_type radius(0);
	    for(size_type i=0;i<pt.size();++i)
	    {
		radius+=(pt[i]*pt[i]);
	    }
	    if(!equals(radius,(value_type)1))
	    {
		Camgen::log(log_level::warning)<<"uniform sphere point has incorrect length "<<radius<<endlog;
	    }
	}
	std::cerr<<".....done."<<std::endl;
    }

    {
	std::cerr<<"Checking uniform 7-sphere generation.....";
	std::cerr.flush();
	vector<value_type,8>pt;
	uniform_sphere<value_type,7,std::random>gen(&pt);
	for(size_type n=0;n<N_events;++n)
	{
	    gen.generate();
	    value_type radius(0);
	    for(size_type i=0;i<pt.size();++i)
	    {
		radius+=(pt[i]*pt[i]);
	    }
	    if(!equals(radius,(value_type)1))
	    {
		Camgen::log(log_level::warning)<<"uniform sphere point has incorrect length "<<radius<<endlog;
	    }
	}
	std::cerr<<".....done."<<std::endl;
    }

    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law1";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    ///////////////////////////////

    m=3;
    nu=1.7;
    smin=20.;
    smax=100.;
    N_bins=500;

    ///////////////////////////////

    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law2";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    ///////////////////////////////

    smin=-(value_type)4;
    smax=-(value_type)1;
    N_bins=100;
    
    ///////////////////////////////
    
    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law3";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////

    nu=-2.5;
    smin=-(value_type)4;
    smax=(value_type)2;
    N_bins=500;
    
    ///////////////////////////////

    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law4";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////

    nu=(value_type)1;
    smin=-(value_type)4;
    smax=-(value_type)0.5;
    N_bins=250;
    
    ///////////////////////////////

    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law5";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////

    smin=(value_type)10;
    smax=(value_type)20;
    N_bins=500;
    
    ///////////////////////////////

    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law6";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////

    nu=(value_type)2.5;
    smin=(value_type)1;
    smax=(value_type)5;
    N_bins=250;
    
    ///////////////////////////////

    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law7";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////
    
    smin=-(value_type)2;
    smax=-(value_type)0.5;
    N_bins=100;
    
    ///////////////////////////////

    {
	std::cerr<<"Checking power-law generator with (m,nu) = ("<<m<<","<<nu<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/pow_law8";
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(0,"0");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////

    value_type mmin(0);
    value_type mmax(200);
    value_type numin(-4);
    value_type numax(4);
    N_bins=100;
    
    ///////////////////////////////
    
    {
	std::cerr<<"Checking 1M random power-law generators within random bounds.....";
	std::cerr.flush();
	power_law<value_type,std::random>* channel=new power_law<value_type,std::random>(&m,&nu);
	for(size_type n=0;n<10000;++n)
	{
	    nu=rn_stream::throw_number(numin,numax);
	    m =rn_stream::throw_number(mmin,mmax);
	    channel->refresh_params();
	    value_type mmax2=mmax*mmax;
	    if(nu>=(value_type)1)
	    {
		smin=rn_stream::throw_number(-mmax2,mmax2);
		if(smin<m*m)
		{
		    smax=rn_stream::throw_number(smin,m*m);
		}
		else
		{
		    smax=rn_stream::throw_number(smin,mmax2);
		}
	    }
	    else
	    {
		smin=rn_stream::throw_number(-mmax2,mmax2);
		value_type x=rn_stream::throw_number(-mmax2,mmax2);
		smax=std::max(smin,x);
		smin=std::min(smin,x);
	    }
	    if(!channel->set_bounds(smin,smax))
	    {
		return 1;
	    }
	    if(channel->normalisable())
	    {
		for(size_type i=0;i<100;++i)
		{
		    channel->generate();
		    if(!channel->check())
		    {
			Camgen::log(log_level::warning)<<"Power-law generation with mass "<<m<<" and exponent "<<nu<<" failed within ["<<smin<<","<<smax<<"]....."<<endlog;
			return 1;
		    }
		}
	    }
	    else
	    {
		Camgen::log(log_level::warning)<<"Power-law generator with mass "<<m<<" and exponent "<<nu<<" is not normalizable within ["<<smin<<","<<smax<<"]....."<<endlog;
		return 1;
	    }
	}
	delete channel;
	std::cerr<<".....done."<<std::endl;
    }
    
    ///////////////////////////////

    m=(value_type)7;
    smin=(value_type)0;
    smax=(value_type)200;
    N_bins=250;
    
    ///////////////////////////////

    {
	std::cerr<<"Checking Breit-Wigner generator with (m,w) = ("<<m<<","<<w<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/Breit_Wigner1";
	Breit_Wigner<value_type,std::random>* channel=new Breit_Wigner<value_type,std::random>(&m,&w);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(m*m,"m2");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////
    
    smin=-(value_type)200;
    smax=(value_type)60;
    N_bins=500;
    
    ///////////////////////////////
    
    {
	std::cerr<<"Checking Breit-Wigner generator with (m,w) = ("<<m<<","<<w<<") within range ["<<smin<<","<<smax<<"].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/Breit_Wigner2";
	Breit_Wigner<value_type,std::random>* channel=new Breit_Wigner<value_type,std::random>(&m,&w);
	if(!channel->set_bounds(smin,smax))
	{
	    return 1;
	}
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_x_tic(smin,"smin");
	plt->add_x_tic(smax,"smax");
	plt->add_x_tic(m*m,"m2");
	plt->grid=true;
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
    
    ///////////////////////////////
    
    N_events=100000;
    N_bins=5000;
    
    ///////////////////////////////
    
    {
	std::cerr<<"Checking Breit-Wigner generator with (m,w) = ("<<m<<","<<w<<") within range [-inf,+inf].....";
	std::cerr.flush();
	std::string filename="test_output/MC_gen_test/Breit_Wigner3";
	Breit_Wigner<value_type,std::random>* channel=new Breit_Wigner<value_type,std::random>(&m,&w);
	histogram<value_type>hist=channel->testrun(N_events,N_bins);
	plot_script* plt=hist.plot(filename,"postscript color");
	plt->add_plot(channel->pdf_plot(N_events*hist.binsize()));
	plt->xmin=-1000;
	plt->xmax=1000;
	plt->plot();
	delete channel;
	delete plt;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
}

