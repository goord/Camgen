//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/plt_config.h>
#include <Camgen/plt_script.h>
#include <Camgen/file_utils.h>
#include <Camgen/SM.h>

/* * * * * * * * * * * * * * * * * * * * *
 * Tests for standard-model parameters.  *
 *                                       *
 * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    license_print::disable();

    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing standard-model parameters........................."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    bool have_gp=plot_config::gnuplot_path!=NULL;
    std::string fext=have_gp?".eps":".dat/.gp";
    file_utils::create_directory("test_output/SM_params_tests");

    {
	std::cerr<<"Checking standard model LO Higgs width..........";
	std::string filename("test_output/SM_params_tests/Higgs_width");
	double mass,width;
	data_wrapper* data=have_gp?new data_wrapper(&mass,&width):new data_wrapper(filename+".dat",&mass,&width);
	std::cerr.flush();
	for(mass=114;mass<500;mass+=10)
	{
	    SM::set_Higgs_mass(mass);
	    width=SM::W_h0;
	    data->fill();
	}
	data->write();
	data_stream* datastr=new data_stream(data,"1","2");
	datastr->title="SM Higgs width (LO)";
	datastr->style="lines";
	plot_script* plot=new plot_script(filename,"postscript color");
	plot->add_plot(datastr);
	plot->ylog=true;
	plot->xlabel="Mh(GeV)";
	plot->ylabel="Gamma(GeV)";
	plot->plot();
	delete plot;
	std::cerr<<".........done, file "<<filename+fext<<" written."<<std::endl;
    }

    {
	std::cerr<<"Checking standard model strong coupling..........";
	std::string filename("test_output/SM_params_tests/alpha_s");
	double mu,alpha;
	data_wrapper* data=have_gp?new data_wrapper(&mu,&alpha):new data_wrapper(filename+".dat",&mu,&alpha);
	for(mu=5;mu<500;mu+=10)
	{
	    SM::set_QCD_scale(mu);
	    alpha=SM::alpha_s;
	    data->fill();
	}
	data->write();
	data_stream* datastr=new data_stream(data,"1","2");
	datastr->title="strong coupling (LO)";
	datastr->style="lines";
	plot_script* plot=new plot_script(filename,"postscript color");
	plot->add_plot(datastr);
	plot->xlabel="mu(GeV)";
	plot->ylabel="alpha_s";
	plot->plot();
	delete plot;
	std::cerr<<".........done, file "<<filename+fext<<" written."<<std::endl;
    }
}

