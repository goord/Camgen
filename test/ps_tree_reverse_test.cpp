//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/plt_config.h>
#include <Camgen/file_utils.h>
#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/uni_hels.h>
#include <Camgen/QCD_cols.h>
#include <ps_gen_tester.h>
#include <SMPbKsch.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Tests for reverse invariant mass sampling within recursive *
* Monte Carlo trees.                                         *
*                                                            *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef SM::value_type value_type;
    typedef std::random rn_engine;
    
    license_print::disable();

    initial_states::type isgen_type=initial_states::partonic;
    phase_space_generators::type psgen_type=phase_space_generators::recursive_backward_s;
    set_phase_space_generator_type(psgen_type);
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    file_utils::create_directory("test_output/ps_tree_reverse_test");
    
    std::size_t n_evts=10000;
    std::size_t n_bins=100;
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing recursive phase space MC with backward s-sampling in decays......"<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > gamma,gamma");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,2,true>* helgen=uniform_helicities<value_type,1,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,2,rn_engine>test(algo.get_tree_iterator(),"",isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins,false,true,true))
	{
	    std::cerr<<test.ps_gen->cross_section()<<','<<test.uni_gen->cross_section()<<std::endl;
	    return 1;
	}
	if(!(test.uni_gen->cross_section().value==0 and test.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(test.ps_gen->cross_section().value==0 and test.ps_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("gamma > e+,e-");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,2,true>* helgen=uniform_helicities<value_type,1,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,2,rn_engine>test(algo.get_tree_iterator(),"",isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins,false,true,true))
	{
	    std::cerr<<test.ps_gen->cross_section()<<','<<test.uni_gen->cross_section()<<std::endl;
	    return 1;
	}
	if(!(test.uni_gen->cross_section().value==0 and test.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(test.ps_gen->cross_section().value==0 and test.ps_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("W- > mu-,nu_mubar");
	std::string fname("test_results/ps_tree_reverse_test/W_ln~");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,2,true>* helgen=uniform_helicities<value_type,1,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,2,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("t > b,mu+,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/t_bln~");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,3,true>* helgen=uniform_helicities<value_type,1,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/h_WW_2l2n~");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=140;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/h_WW*_2l2n~");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > mu-,nu_mubar,mu+,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/h_2l2n2~");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > mu-,mu+,e-,e+");
	std::string fname("test_results/ps_tree_reverse_test/h_ZZ_4l~");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=140;
	model_type::refresh_widths();
	std::string process("h0 > mu-,mu+,e-,e+");
	std::string fname("test_results/ps_tree_reverse_test/h_ZZ*_4l~");
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_m_min(1,2,10);
	test.set_m_min(3,4,10);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing recursive phase space MC with backward s-sampling in scattering.."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > u,dbar");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(algo.get_tree_iterator(),"",isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins,false,true,true))
	{
	    std::cerr<<test.ps_gen->cross_section()<<','<<test.uni_gen->cross_section()<<std::endl;
	    return 1;
	}
	if(!(test.uni_gen->cross_section().value==0 and test.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(test.ps_gen->cross_section().value==0 and test.ps_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > t,tbar");
	double E1=100;
	double E2=100;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(algo.get_tree_iterator(),"",isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins,false,true,true))
	{
	    std::cerr<<test.ps_gen->cross_section()<<','<<test.uni_gen->cross_section()<<std::endl;
	    return 1;
	}
	if(!(test.uni_gen->cross_section().value==0 and test.uni_gen->cross_section().error==0))
	{
	    return 1;
	}
	if(!(test.ps_gen->cross_section().value==0 and test.ps_gen->cross_section().error==0))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
    
    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > mu+,mu-");
	std::string fname("test_results/ps_tree_reverse_test/ee_ll~");
	double E1=250;
	double E2=250;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
	
    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,mu- > nu_ebar,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/emu_2n~");
	double E1=250;
	double E2=250;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
    
    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_e,nu_ebar");
	std::string fname("test_results/ps_tree_reverse_test/ee_2n~");
	double E1=250;
	double E2=250;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
    
    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,mu- > nu_ebar,h0,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/emu_h2n~");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_ebar,h0,nu_e");
	std::string fname("test_results/ps_tree_reverse_test/ee_h2n~");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
    
    {
	Camgen::log.enable_level=log_level::error;
	std::string process("nu_e,nu_mu > nu_e,Z,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/nn_Z2n~");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_ebar,mu+,mu-,nu_e");
	std::string fname("test_results/ps_tree_reverse_test/ee_2mu2ne~");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_m_min(2,3,mmin);
	test.set_m_min(1,4,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_mubar,mu+,mu-,nu_mu");
	std::string fname("test_results/ps_tree_reverse_test/ee_2mu2nmu~");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_m_min(2,3,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_e,W+,W-,nu_ebar");
	std::string fname("test_results/ps_tree_reverse_test/ee_2W2n~");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_m_min(1,4,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_e,Z,Z,nu_ebar");
	std::string fname("test_results/ps_tree_reverse_test/ee_2Z2n~");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_m_min(1,4,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > nu_e,W+,W-,nu_ebar");
	std::string fname("test_results/ps_tree_reverse_test/ee_2W2n~");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_m_min(1,4,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
}


