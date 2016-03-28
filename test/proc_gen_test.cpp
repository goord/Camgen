//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <proc_gen_tester.h>

/* * * * * * * * * * * * * * * * *
 * Tests for process generators. *
 *                               *
 * * * * * * * * * * * * * * * * */

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
    file_utils::create_directory("test_output/proc_gen_test");
    
    std::size_t n_evts=10000;
    std::size_t n_bins=100;
    
    std::cout<<"-----------------------------------------------------------"<<std::endl;
    std::cout<<"testing process generator for decays......................."<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > gamma,gamma");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,2,rn_engine> tester(algo.get_tree_iterator(),"");
	if(!tester.run(n_evts,n_bins,false,true,true))
	{
	    return 1;
	}
	if(!(tester.uni_gen->cross_section().value==0 and tester.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(tester.proc_gen->cross_section().value==0 and tester.proc_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("gamma > e+,e-");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,2,rn_engine> tester(algo.get_tree_iterator(),"");
	if(!tester.run(n_evts,n_bins,false,true,true))
	{
	    return 1;
	}
	if(!(tester.uni_gen->cross_section().value==0 and tester.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(tester.proc_gen->cross_section().value==0 and tester.proc_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::cerr<<"Checking event correlation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_factory<model_type,1,4,rn_engine> proc_gen_fac;
	process_generator<model_type,1,4,rn_engine>* proc_gen=proc_gen_fac.create_generator(algo.get_tree_iterator());
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w1(0);
	std::size_t n1(0);
	do
	{
	    proc_gen->generate();
	    w1=proc_gen->weight();
	    n1++;
	}
	while(w1==0);
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w2(0);
	std::size_t n2(0);
	do
	{
	    proc_gen->generate();
	    w2=proc_gen->weight();
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

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("W- > mu-,nu_mubar");
	std::string fname("test_output/proc_gen_test/W_ln");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,2,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("t > b,mu+,nu_mu");
	std::string fname("test_output/proc_gen_test/t_bln");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,3>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,3,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_output/proc_gen_test/h_WW_2l2n");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=140;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_output/proc_gen_test/h_WW*_2l2n");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    tester.proc_gen->print_ps(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > mu-,nu_mubar,mu+,nu_mu");
	std::string fname("test_output/proc_gen_test/h_2l2n2");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    tester.proc_gen->print_ps(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > mu-,mu+,e-,e+");
	std::string fname("test_output/proc_gen_test/h_ZZ_4l");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=140;
	model_type::refresh_widths();
	std::string process("h0 > nu_mu,nu_mubar,nu_e,nu_ebar");
	std::string fname("test_output/proc_gen_test/h_ZZ*_4n");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    std::cout<<"-----------------------------------------------------------"<<std::endl;
    std::cout<<"testing process generator for scattering..................."<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > u,dbar");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,2,rn_engine> tester(algo.get_tree_iterator(),"");
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins,false,true,true))
	{
	    return 1;
	}
	if(!(tester.uni_gen->cross_section().value==0 and tester.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(tester.proc_gen->cross_section().value==0 and tester.proc_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > t,tbar");
	double E1=100;
	double E2=100;
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,2,rn_engine> tester(algo.get_tree_iterator(),"");
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins,false,true,true))
	{
	    return 1;
	}
	if(!(tester.uni_gen->cross_section().value==0 and tester.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(tester.proc_gen->cross_section().value==0 and tester.proc_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_e,h0,nu_ebar");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking event correlation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator_factory<model_type,2,3,rn_engine> proc_gen_fac;
	process_generator<model_type,2,3,rn_engine>* proc_gen=proc_gen_fac.create_generator(algo.get_tree_iterator());
	proc_gen->set_beam_energy(-1,E1);
	proc_gen->set_beam_energy(-2,E2);
	proc_gen->refresh_Ecm();
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w1(0);
	std::size_t n1(0);
	do
	{
	    proc_gen->generate();
	    w1=proc_gen->weight();
	    n1++;
	}
	while(w1==0);
	random_number_stream<model_type::value_type,rn_engine>::reset_engine();
	model_type::value_type w2(0);
	std::size_t n2(0);
	do
	{
	    proc_gen->generate();
	    w2=proc_gen->weight();
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

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > mu+,mu-");
	std::string fname("test_output/proc_gen_test/ee_ll");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,2,rn_engine> tester(algo.get_tree_iterator(),fname);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > e+,e-");
	double E1=10;
	double E2=10;
	double pTmin=1;
	std::string fname("test_output/proc_gen_test/ee_ee");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,2,rn_engine> tester(algo.get_tree_iterator(),fname);
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	tester.set_pT_min(pTmin,1);
	tester.set_pT_min(pTmin,2);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > Z,Z");
	double E1=100;
	double E2=100;
	std::string fname("test_output/proc_gen_test/ee_ZZ");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,2,rn_engine> tester(algo.get_tree_iterator(),fname);
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_e,nu_ebar,h0");
	double E1=250;
	double E2=250;
	std::string fname("test_output/proc_gen_test/ee_hll");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,3,rn_engine> tester(algo.get_tree_iterator(),fname);
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > mu-,nu_mubar,W+");
	double E1=250;
	double E2=250;
	std::string fname("test_output/proc_gen_test/ee_lnuW");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,3,rn_engine> tester(algo.get_tree_iterator(),fname);
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > mu-,nu_mubar,u,dbar");
	double E1=100;
	double E2=100;
	std::string fname("test_output/proc_gen_test/ee_lnuqq");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > mu-,nu_mubar,W+,Z");
	double E1=500;
	double E2=500;
	std::string fname("test_output/proc_gen_test/ee_munWZ");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_e,nu_ebar,Z,Z");
	double E1=500;
	double E2=500;
	std::string fname("test_output/proc_gen_test/ee_nnZZ");
	std::cerr<<"Checking process generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,2,4,rn_engine> tester(algo.get_tree_iterator(),fname);
	tester.set_beam_energy(-1,E1);
	tester.set_beam_energy(-2,E2);
	if(!tester.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }
}
