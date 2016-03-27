//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>

/* * * * * * * * * * * * * * * *
 * Tests for event generators. *
 *                             *
 * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef std::random rn_engine;
    
    license_print::disable();

    initial_states::type isgen_type=initial_states::partonic;
    phase_space_generators::type psgen_type=phase_space_generators::recursive;
    
    set_initial_state_type(isgen_type);
    set_phase_space_generator_type(psgen_type);
    set_s_pair_generation_mode(s_pair_generation_modes::symmetric);
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    file_utils::create_directory("test_output/evt_gen_test");
    
    std::size_t n_init_evts=1000;
    std::size_t n_evts=10000;
    std::size_t n_bins=100;

    set_subprocess_events(n_init_evts);
    
    std::cout<<"-----------------------------------------------------------"<<std::endl;
    std::cout<<"testing event generator for decays........................."<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;

    {
        
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > gamma,gamma");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
        event_generator_factory<model_type,1,2,rn_engine> factory;
        event_generator<model_type,1,2,rn_engine>* evt_gen=factory.create_generator(algo);
        for(std::size_t n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
            }
        }
        evt_gen->refresh_cross_section();
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
        
	Camgen::log.enable_level=log_level::error;
	std::string process("gamma > l+,l-");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
        event_generator_factory<model_type,1,2,rn_engine> factory;
        event_generator<model_type,1,2,rn_engine>* evt_gen=factory.create_generator(algo);
        for(std::size_t n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
            }
        }
        evt_gen->refresh_cross_section();
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > l+,l-,l+,l-");
	std::cerr<<"Checking event correlation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w1(0);
	std::size_t n1(0);
	do
	{
	    evt_gen->generate();
	    w1=evt_gen->weight();
	    n1++;
	}
	while(w1==0);
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w2(0);
	std::size_t n2(0);
	do
	{
	    evt_gen->generate();
	    w2=evt_gen->weight();
	    n2++;
	}
	while(w2==0);
	if(n1!=n2 or w1!=w2)
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    std::cout<<"-----------------------------------------------------------"<<std::endl;
    std::cout<<"testing event generator for scattering....................."<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,dbar > l+,l-");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
        event_generator_factory<model_type,2,2,rn_engine> factory;
        event_generator<model_type,2,2,rn_engine>* evt_gen=factory.create_generator(algo);
	evt_gen->set_beam_energy(-1,E1);
	evt_gen->set_beam_energy(-2,E2);
        for(std::size_t n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
            }
        }
        evt_gen->refresh_cross_section();
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("q,qbar > t,tbar");
	double E1=100;
	double E2=100;
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
        event_generator_factory<model_type,2,2,rn_engine> factory;
        event_generator<model_type,2,2,rn_engine>* evt_gen=factory.create_generator(algo);
	evt_gen->set_beam_energy(-1,E1);
	evt_gen->set_beam_energy(-2,E2);
        for(std::size_t n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
            }
        }
        evt_gen->refresh_cross_section();
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > h0,nu,nu_bar");
	std::cerr<<"Checking event correlation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
        event_generator_factory<model_type,2,3,rn_engine> factory;
        event_generator<model_type,2,3,rn_engine>* evt_gen=factory.create_generator(algo);
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w1(0);
	std::size_t n1(0);
	do
	{
	    evt_gen->generate();
	    w1=evt_gen->weight();
	    n1++;
	}
	while(w1==0);
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w2(0);
	std::size_t n2(0);
	do
	{
	    evt_gen->generate();
	    w2=evt_gen->weight();
	    n2++;
	}
	while(w2==0);
	if(n1!=n2 or w1!=w2)
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    return 0;
}
