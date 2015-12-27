//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/plt_config.h>
#include <Camgen/stdrand.h>
#include <Camgen/histogram.h>
#include <Camgen/multiplot.h>
#include <Camgen/power_law.h>
#include <Camgen/Breit_Wigner.h>
#include <Camgen/val_gen_grid.h>
#include <Camgen/parni_int.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Facility testing adaptive grids in parni for various intgrands in various *
 * dimensions.                                                               *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing adaptive integration algorithm..................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    typedef double value_type;
    typedef std::size_t size_type;
    
    //////////////////////////////////////////////////////////////////

    size_type N_events = 100000;
    size_type N_batch  = 100;
    size_type N_bins   = 100;
    grid_modes::type mode = grid_modes::maximum_weights;
    std::string fext = plot_config::gnuplot_path==NULL?".dat/.gp":".eps";

    //////////////////////////////////////////////////////////////////

    size_type offset=N_events/N_bins;
    value_type pi=std::acos(-(value_type)1);
    
    {
	value_type x=1,xmin=0,xmax=10,m=5,w=0.75;
	value_type wght;
	std::cerr<<"Checking 1D parni on Cauchy distribution centered at "<<m<<" with width "<<w<<" within ["<<xmin<<","<<xmax<<"]..........";
        std::cerr.flush();
	std::string filename("plots/parni1D_Cauchy");
	parni_integrator<value_type,1,std::random>* gen=new parni_integrator<value_type,1,std::random>(&x,xmin,xmax,N_bins,mode);
	histogram<value_type>hist(&x,&wght);
	value_type result=(std::atan((xmax-m)/w)-std::atan((xmin-m)/w))/w;
	value_type evt,c1,c2;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    wght=gen->weight();
	    hist.store();
	    gen->integrand()=(value_type)1/(std::pow(x-m,(int)2)+w*w);
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		evt=(value_type)n;
		MC_integral<value_type>integral=gen->cross_section();
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist.make(N_bins);

	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin;
	plot1->xmax=xmax;
	gen->print_rectangles(plot1);
	std::ostringstream sstrm;
	sstrm<<"1/((x-"<<m<<")**2 + "<<w*w<<")";
	function_stream* fstr=new function_stream(sstrm.str());
	fstr->title="density";
	fstr->style="lines";
	plot1->add_plot(fstr);
	
	plot_script* plot2=hist.plot(filename+"_x");
	
	plot_script* plot3=new plot_script(filename+"_err");
	plot3->grid=true;
	plot3->xmin=0;
	plot3->xmax=N_events;
	plot3->ymin=gen->cross_section().error/result;
	plot3->ymax=1;
	plot3->autoscale=false;
	plot3->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot3->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"3");
	interr->title="MC error";
	interr->style="lines";
	plot3->add_plot(interr);
	multi_plot<3,1>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot3,0,0);
	plotall.add_plot(plot2,1,0);
	plotall.add_plot(plot1,2,0);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete gen;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }

    {
	value_type x=1,xmin=0.1,xmax=0.9;
	value_type wght;
	std::cerr<<"Checking 1D parni 1/(x(1-x)) within ["<<xmin<<","<<xmax<<"]..........";
        std::cerr.flush();
	std::string filename("plots/parni1D_beta11");
	parni_integrator<value_type,1,std::random>* gen=new parni_integrator<value_type,1,std::random>(&x,xmin,xmax,N_bins,mode);
	histogram<value_type>hist(&x,&wght);
	value_type result=std::log((1-xmin)*xmax/(1-xmax)/xmin);
	value_type evt,c1,c2;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    wght=gen->weight();
	    hist.store();
	    gen->integrand()=((value_type)1)/(x*(1-x));
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		evt=(value_type)n;
		MC_integral<value_type>integral=gen->cross_section();
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist.make(N_bins);

	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin;
	plot1->xmax=xmax;
	gen->print_rectangles(plot1);
	function_stream* fstr=new function_stream("1/(x*(1-x))");
	fstr->title="density";
	fstr->style="lines";
	plot1->add_plot(fstr);
	
	plot_script* plot2=hist.plot(filename+"_x");
	
	plot_script* plot3=new plot_script(filename+"_err");
	plot3->grid=true;
	plot3->xmin=0;
	plot3->xmax=N_events;
	plot3->ymin=gen->cross_section().error/result;
	plot3->ymax=1;
	plot3->autoscale=false;
	plot3->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot3->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"3");
	interr->title="MC error";
	interr->style="lines";
	plot3->add_plot(interr);
	multi_plot<3,1>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot3,0,0);
	plotall.add_plot(plot2,1,0);
	plotall.add_plot(plot1,2,0);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete gen;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    {
	value_type x=1,xmin=-2*pi,xmax=-xmin;
	value_type wght;
	std::cerr<<"Checking 1D parni cos(x)^2*exp(-x^2/8) within ["<<xmin<<","<<xmax<<"]..........";
        std::cerr.flush();
	std::string filename("plots/parni1D_cosgauss");
	parni_integrator<value_type,1,std::random>* gen=new parni_integrator<value_type,1,std::random>(&x,xmin,xmax,N_bins,mode);
	histogram<value_type>hist(&x,&wght);
	value_type result=2.501415718;
	value_type evt,c1,c2;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    wght=gen->weight();
	    hist.store();
	    gen->integrand()=std::pow(std::cos(x),(int)2)*std::exp(-x*x/(value_type)8);
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		evt=(value_type)n;
		MC_integral<value_type>integral=gen->cross_section();
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist.make(N_bins);

	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin;
	plot1->xmax=xmax;
	gen->print_rectangles(plot1);
	function_stream* fstr=new function_stream("(cos(x)**2)*exp(-x**2/8)");
	fstr->title="density";
	fstr->style="lines";
	plot1->add_plot(fstr);
	
	plot_script* plot2=hist.plot(filename+"_x");
	
	plot_script* plot3=new plot_script(filename+"_err");
	plot3->grid=true;
	plot3->xmin=0;
	plot3->xmax=N_events;
	plot3->ymin=gen->cross_section().error/result;
	plot3->ymax=1;
	plot3->autoscale=false;
	plot3->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot3->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"3");
	interr->title="MC error";
	interr->style="lines";
	plot3->add_plot(interr);
	multi_plot<3,1>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot3,0,0);
	plotall.add_plot(plot2,1,0);
	plotall.add_plot(plot1,2,0);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete gen;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    {
	value_type x=0,xmin=0,xmax=1;
	value_type smin=1,smax=5;
	value_type wght;
	std::cerr<<"Checking power-law(1D parni) on (x-1)e^{-x} within ["<<smin<<","<<smax<<"]..........";
        std::cerr.flush();
	std::string filename("plots/parni1Dpl");
	value_type nu(1);
	parni_integrator<value_type,1,std::random>* gen=new parni_integrator<value_type,1,std::random>(&x,xmin,xmax,N_bins,mode);
	power_law<value_type,std::random>* genmap=new power_law<value_type,std::random>(NULL,&nu);
	genmap->refresh_params();
	genmap->set_lower_bound(smin);
	genmap->set_upper_bound(smax);
	value_type& s=genmap->value();
	histogram<value_type>hist(&s,&wght);
	value_type result=0.3342;
	value_type evt,c1,c2,f;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	MC_integrator<value_type>* genmap_int=new MC_generator_wrapper<value_type>(genmap);
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    s=genmap->map(x);
	    genmap_int->evaluate_weight();
	    wght=gen->weight()*genmap_int->weight();
	    hist.store();
	    f=(s-(value_type)1)/std::exp(s);
	    gen->integrand()=genmap_int->weight()*f;
	    gen->update();
	    genmap_int->integrand()=gen->weight()*f;
	    genmap_int->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		evt=(value_type)n;
		MC_integral<value_type>integral=genmap_int->cross_section();
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist.make(N_bins);
	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin;
	plot1->xmax=xmax;
	gen->print_rectangles(plot1);
	function_stream* fstr=new function_stream("0");
	plot1->add_plot(fstr);
	function_stream* fstr2=new function_stream("0.6");
	plot1->add_plot(fstr2);
	
	plot_script* plot2=hist.plot(filename+"_x");
	
	plot_script* plot3=new plot_script(filename+"_err");
	plot3->grid=true;
	plot3->xmin=0;
	plot3->xmax=N_events;
	plot3->ymin=genmap_int->cross_section().error/result;
	plot3->ymax=1;
	plot3->autoscale=false;
	plot3->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot3->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"3");
	interr->title="MC error";
	interr->style="lines";
	plot3->add_plot(interr);
	multi_plot<3,1>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot3,0,0);
	plotall.add_plot(plot2,1,0);
	plotall.add_plot(plot1,2,0);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete gen;
	delete genmap_int;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    {
	value_type x=0,xmin=0,xmax=1;
	value_type smin=0,smax=4;
	value_type wght;
	std::cerr<<"Checking Breit-Wigner(1D parni) on x/((x-2)^2+0.1*x) within ["<<smin<<","<<smax<<"]..........";
        std::cerr.flush();
	std::string filename("plots/parni1DBW");
	value_type m=std::sqrt(2),w=0.1/m;
	parni_integrator<value_type,1,std::random>* gen=new parni_integrator<value_type,1,std::random>(&x,xmin,xmax,N_bins,mode);
	Breit_Wigner<value_type,std::random>* genmap=new Breit_Wigner<value_type,std::random>(&m,&w);
	MC_integrator<value_type>* genmap_int=new MC_generator_wrapper<value_type>(genmap);
	genmap->refresh_params();
	genmap->set_lower_bound(smin);
	genmap->set_upper_bound(smax);
	value_type& s=genmap->value();
	histogram<value_type>hist(&s,&wght);
	value_type result=11.91253775;
	value_type evt,c1,c2,f;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    s=genmap->map(x);
	    genmap_int->evaluate_weight();
	    wght=gen->weight()*genmap_int->weight();
	    hist.store();
	    f=s/(std::pow(s-m*m,(int)2)+m*w*s);
	    gen->integrand()=genmap_int->weight()*f;
	    gen->update();
	    genmap_int->integrand()=gen->weight()*f;
	    genmap_int->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		evt=(value_type)n;
		MC_integral<value_type>integral=genmap_int->cross_section();
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist.make(N_bins);
	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin;
	plot1->xmax=xmax;
	gen->print_rectangles(plot1);
	function_stream* fstr=new function_stream("0");
	plot1->add_plot(fstr);
	function_stream* fstr2=new function_stream("100");
	plot1->add_plot(fstr2);
	
	plot_script* plot2=hist.plot(filename+"_x");
	
	plot_script* plot3=new plot_script(filename+"_err");
	plot3->grid=true;
	plot3->xmin=0;
	plot3->xmax=N_events;
	plot3->ymin=genmap_int->cross_section().error/result;
	plot3->ymax=1;
	plot3->autoscale=false;
	plot3->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot3->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"3");
	interr->title="MC error";
	interr->style="lines";
	plot3->add_plot(interr);
	multi_plot<3,1>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot3,0,0);
	plotall.add_plot(plot2,1,0);
	plotall.add_plot(plot1,2,0);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete gen;
	delete genmap_int;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    {
	value_type x=1,x0min=0,xmin=1,xmax=5,x0max=10,m=5,w=0.75;
	value_type wght;
	std::cerr<<"Checking 1D parni subgrid on Cauchy centered at "<<m<<" with width "<<w<<" within ["<<xmin<<","<<xmax<<"]..........";
        std::cerr.flush();
	std::string filename("plots/subparni1D_Cauchy");
	parni_integrator<value_type,1,std::random>* gen=new parni_integrator<value_type,1,std::random>(&x,x0min,x0max,N_bins,mode);
	N_batch=20;
	for(size_type n=0;n<200;++n)
	{
	    gen->generate();
	    gen->integrand()=(value_type)1/(std::pow(x-m,(int)2)+w*w);
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0)
	    {
		gen->adapt();
	    }
	}

	MC_integrator<value_type>* subgen=new MC_generator_wrapper<value_type>(gen->create_sub_grid(xmin,xmax));
	histogram<value_type>hist(&x,&wght);
	value_type subresult=(std::atan((xmax-m)/w)-std::atan((xmin-m)/w))/w;
	value_type result=(std::atan((x0max-m)/w)-std::atan((x0min-m)/w))/w;
	value_type evt,c1,c2,c3,c4;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	intdata->add_leaf(&c3);
	intdata->add_leaf(&c4);
	N_batch=50;
	for(size_type n=0;n<2*N_events;++n)
	{
	    subgen->generate();
	    wght=subgen->weight();
	    hist.store();
	    subgen->integrand()=(value_type)1/(std::pow(x-m,(int)2)+w*w);
	    subgen->update();
	    subgen->refresh_cross_section();
	    gen->refresh_cross_section();
	    if(n%N_batch==0)
	    {
		subgen->adapt();
	    }
	    if(n%offset==0)
	    {
		evt=(value_type)n;
		MC_integral<value_type>integral=gen->cross_section();
		MC_integral<value_type>subintegral=subgen->cross_section();
		c1=std::abs(subintegral.value-subresult)/subresult;
		c2=subintegral.error/subresult;
		c3=std::abs(integral.value-result)/result;
		c4=integral.error/result;
		intdata->fill();
	    }
	}
	hist.make(N_bins);
	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=x0min;
	plot1->xmax=x0max;
	gen->print_rectangles(plot1);
	std::ostringstream sstrm;
	sstrm<<"1/((x-"<<m<<")**2 + "<<w*w<<")";
	function_stream* fstr=new function_stream(sstrm.str());
	fstr->title="density";
	fstr->style="lines";
	plot1->add_plot(fstr);
	
	plot_script* plot2=hist.plot(filename+"_x");
	
	plot_script* plot3=new plot_script(filename+"_err");
	plot3->grid=true;
	plot3->xmin=0;
	plot3->xmax=N_events;
	plot3->ymin=subgen->cross_section().error/result;
	plot3->ymax=1;
	plot3->autoscale=false;
	plot3->ylog=true;
	
	data_stream* intsubval=new data_stream(intdata,"2");
	intsubval->title="delta(subMC,exact)";
	intsubval->style="lines";
	plot3->add_plot(intsubval);
	data_stream* intsuberr=new data_stream(intdata,"3");
	intsuberr->title="subMC error";
	intsuberr->style="lines";
	plot3->add_plot(intsuberr);
	data_stream* intval=new data_stream(intdata,"4");
	intval->title="delta(totMC,exact)";
	intval->style="lines";
	plot3->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"5");
	interr->title="totMC error";
	interr->style="lines";
	plot3->add_plot(interr);
	multi_plot<3,1>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot3,0,0);
	plotall.add_plot(plot2,1,0);
	plotall.add_plot(plot1,2,0);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete subgen;
	delete gen;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    {
	value_type smin=-2,smax=-0.25,s0min=-5,s0max=-0.1;
	value_type wght;
	std::cerr<<"Checking power-law(1D parni) on -x*exp(-x^2) within ["<<smin<<","<<smax<<"]..........";
        std::cerr.flush();
	std::string filename("plots/parni1Dplsub");
	value_type nu(1);
	power_law<value_type,std::random>* genmap=new power_law<value_type,std::random>(NULL,&nu);
	MC_integrator<value_type>* genmap_int=new MC_generator_wrapper<value_type>(genmap);
	adaptive_value_generator<value_type,std::random>* gen=new adaptive_value_generator<value_type,std::random>(genmap,N_bins,mode);
	MC_integrator<value_type>* gen_int=new MC_generator_wrapper<value_type>(gen);
	gen->set_min_lower_bound(s0min);
	gen->set_max_upper_bound(s0max);
	gen->set_bounds(smin,smax);
	value_type& s=gen->value();
	histogram<value_type>hist(&s,&wght);
	value_type result=0.4605;
	value_type evt,c1,c2,f;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	for(size_type n=0;n<N_events;++n)
	{
	    gen_int->generate();
	    wght=gen_int->weight();
	    hist.store();
	    f=-s*std::exp(-s*s);
	    gen_int->integrand()=f;
	    gen_int->update();
	    gen_int->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen_int->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		evt=(value_type)n;
		MC_integral<value_type>integral=gen_int->cross_section();
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist.make(N_bins);
	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=0;
	plot1->xmax=1;
	gen->print_rectangles(plot1);
	function_stream* fstr=new function_stream("0");
	plot1->add_plot(fstr);
	function_stream* fstr2=new function_stream("2");
	plot1->add_plot(fstr2);
	
	plot_script* plot2=hist.plot(filename+"_x");
	
	plot_script* plot3=new plot_script(filename+"_err");
	plot3->grid=true;
	plot3->xmin=0;
	plot3->xmax=N_events;
	plot3->ymin=gen_int->cross_section().error/result;
	plot3->ymax=1;
	plot3->autoscale=false;
	plot3->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot3->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"3");
	interr->title="MC error";
	interr->style="lines";
	plot3->add_plot(interr);
	multi_plot<3,1>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot3,0,0);
	plotall.add_plot(plot2,1,0);
	plotall.add_plot(plot1,2,0);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete gen_int;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    {
	N_bins=500;
	offset=N_events/N_bins;
	vector<value_type,2>x;
	vector<value_type,2>xmin;
	xmin.assign(-1);
	vector<value_type,2>xmax;
	xmax.assign(1);
	value_type wght;
	value_type sigma=0.5;
	std::cerr<<"Checking 2D parni on 2D Gauss with sigma "<<sigma<<" within ["<<xmin[0]<<","<<xmax[0]<<"] x ["<<xmin[1]<<","<<xmax[1]<<"]..........";
        std::cerr.flush();
	std::string filename("plots/parni2D_Gauss");
	parni_integrator<value_type,2,std::random>* gen=new parni_integrator<value_type,2,std::random>(&x,xmin,xmax,N_bins,mode);
	histogram<value_type>hist1(&x[0],&wght);
	histogram<value_type>hist2(&x[1],&wght);
	value_type evt,c1,c2;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	data_wrapper* xdata=new data_wrapper(&x[0],&x[1]);
	value_type result=0.7780675800;
	value_type sigma2=sigma*sigma;
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    xdata->fill();
	    wght=gen->weight();
	    hist1.store();
	    hist2.store();
	    gen->integrand()=std::exp((-x[0]*x[0]-x[1]*x[1])/sigma2);
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		MC_integral<value_type>integral=gen->cross_section();
		evt=(value_type)n;
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist1.make(N_bins/5,2,2);
	hist2.make(N_bins/5,2,2);
	
	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin[0];
	plot1->xmax=xmax[0];
	plot1->ymin=xmin[1];
	plot1->ymax=xmax[1];
	gen->print_rectangles(plot1);
	data_stream* points=new data_stream(xdata,"1","2");
	points->style="dots";
	plot1->add_plot(points);

	plot_script* plot2=hist1.plot(filename+"_x1");
	plot_script* plot3=hist2.plot(filename+"_x2");

	plot_script* plot4=new plot_script(filename+"_err");
	plot4->grid=true;
	plot4->xmin=0;
	plot4->xmax=N_events;
	plot4->ymin=gen->cross_section().error/result;
	plot4->ymax=1;
	plot4->autoscale=false;
	plot4->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"1","2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot4->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"1","3");
	interr->title="MC error";
	interr->style="lines";
	plot4->add_plot(interr);
	
	multi_plot<2,2>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot2,0,0);
	plotall.add_plot(plot3,0,1);
	plotall.add_plot(plot1,1,0);
	plotall.add_plot(plot4,1,1);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete plot4;
	delete gen;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    
    {
	N_events=50000;
	N_bins=1000;
	offset=N_events/N_bins;
	vector<value_type,2>x;
	vector<value_type,2>xmin;
	xmin.assign(-pi);
	vector<value_type,2>xmax;
	xmax.assign(pi);
	value_type wght;
	std::cerr<<"Checking 2D parni on double cos^2 within ["<<xmin[0]<<","<<xmax[0]<<"] x ["<<xmin[1]<<","<<xmax[1]<<"].........";
        std::cerr.flush();
	std::string filename("plots/parni2D_cos2");
	parni_integrator<value_type,2,std::random>* gen=new parni_integrator<value_type,2,std::random>(&x,xmin,xmax,N_bins,mode);
	histogram<value_type>hist1(&x[0],&wght);
	histogram<value_type>hist2(&x[1],&wght);
	value_type evt,c1,c2;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	data_wrapper* xdata=new data_wrapper(&x[0],&x[1]);
	value_type result=pi*pi;
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    xdata->fill();
	    wght=gen->weight();
	    hist1.store();
	    hist2.store();
	    gen->integrand()=std::pow(std::cos(x[0])*std::cos(x[1]),(int)2);
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		MC_integral<value_type>integral=gen->cross_section();
		evt=(value_type)n;
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist1.make(N_bins/5,2,2);
	hist2.make(N_bins/5,2,2);
	
	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin[0];
	plot1->xmax=xmax[0];
	plot1->ymin=xmin[1];
	plot1->ymax=xmax[1];
	gen->print_rectangles(plot1);
	data_stream* points=new data_stream(xdata,"1","2");
	points->style="dots";
	plot1->add_plot(points);

	plot_script* plot2=hist1.plot(filename+"_x1");
	plot_script* plot3=hist2.plot(filename+"_x2");

	plot_script* plot4=new plot_script(filename+"_err");
	plot4->grid=true;
	plot4->xmin=0;
	plot4->xmax=N_events;
	plot4->ymin=gen->cross_section().error/result;
	plot4->ymax=1;
	plot4->autoscale=false;
	plot4->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"1","2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot4->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"1","3");
	interr->title="MC error";
	interr->style="lines";
	plot4->add_plot(interr);
	
	multi_plot<2,2>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot2,0,0);
	plotall.add_plot(plot3,0,1);
	plotall.add_plot(plot1,1,0);
	plotall.add_plot(plot4,1,1);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete plot4;
	delete gen;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
    
    {
	N_events=50000;
	N_bins=1000;
	offset=N_events/N_bins;
	vector<value_type,2>x;
	vector<value_type,2>xmin;
	xmin.assign(-1);
	vector<value_type,2>xmax;
	xmax.assign(1);
	value_type wght;
	value_type rmax=1,w=0.1,mu=0.5;
	std::cerr<<"Checking 2D parni on 2D Cauchy with radius "<<mu<<" and width "<<w<<" within D(0,"<<rmax<<")..........";
        std::cerr.flush();
	std::string filename("plots/parni2D_Cauchy");
	parni_integrator<value_type,2,std::random>* gen=new parni_integrator<value_type,2,std::random>(&x,xmin,xmax,N_bins,mode);
	histogram<value_type>hist1(&x[0],&wght);
	histogram<value_type>hist2(&x[1],&wght);
	value_type evt,c1,c2;
	data_wrapper* intdata=new data_wrapper(&evt,&c1,&c2);
	data_wrapper* xdata=new data_wrapper(&x[0],&x[1]);
	value_type mu2=mu*mu;
	value_type w2=w*w;
	value_type rmax2=rmax*rmax;
	value_type result=pi*(std::atan((rmax2-mu2)/w)+std::atan(mu2/w))/w;
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    xdata->fill();
	    wght=gen->weight();
	    hist1.store();
	    hist2.store();
	    value_type r2=x[0]*x[0]+x[1]*x[1];
	    gen->integrand()=(r2<rmax2)?((value_type)1/(std::pow(r2-mu2,(int)2)+w2)):(value_type)0;
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	    if(n%offset==0 and n!=0)
	    {
		MC_integral<value_type>integral=gen->cross_section();
		evt=(value_type)n;
		c1=std::abs(integral.value-result)/result;
		c2=integral.error/result;
		intdata->fill();
	    }
	}
	hist1.make(N_bins/5,2,2);
	hist2.make(N_bins/5,2,2);
	
	plot_script* plot1=new plot_script(filename+"_grid");
	plot1->xmin=xmin[0];
	plot1->xmax=xmax[0];
	plot1->ymin=xmin[1];
	plot1->ymax=xmax[1];
	gen->print_rectangles(plot1);
	data_stream* points=new data_stream(xdata,"1","2");
	points->style="dots";
	plot1->add_plot(points);

	plot_script* plot2=hist1.plot(filename+"_x1");
	plot_script* plot3=hist2.plot(filename+"_x2");

	plot_script* plot4=new plot_script(filename+"_err");
	plot4->grid=true;
	plot4->xmin=0;
	plot4->xmax=N_events;
	plot4->ymin=gen->cross_section().error/result;
	plot4->ymax=1;
	plot4->autoscale=false;
	plot4->ylog=true;
	
	data_stream* intval=new data_stream(intdata,"1","2");
	intval->title="delta(MC,exact)";
	intval->style="lines";
	plot4->add_plot(intval);
	data_stream* interr=new data_stream(intdata,"1","3");
	interr->title="MC error";
	interr->style="lines";
	plot4->add_plot(interr);
	
	multi_plot<2,2>plotall(filename,"postscript enhanced color");
	plotall.add_plot(plot2,0,0);
	plotall.add_plot(plot3,0,1);
	plotall.add_plot(plot1,1,0);
	plotall.add_plot(plot4,1,1);
	plotall.plot();
	delete plot1;
	delete plot2;
	delete plot3;
	delete plot4;
	delete gen;
	std::cerr<<"...........done, file "<<filename+fext<<" written."<<std::endl;
    }
}

