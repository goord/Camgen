//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/plt_config.h>
#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/uni_hels.h>
#include <Camgen/qcd_cols.h>
#include <ps_gen_tester.h>
#include <SMPbKsch.h>

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef SM::value_type value_type;
    typedef std::random rn_engine;
    
    license_print::disable();

    initial_states::type isgen_type=initial_states::partonic;
    phase_space_generators::type psgen_type=phase_space_generators::recursive;
    set_phase_space_generator_type(psgen_type);
    set_s_pair_generation_mode(s_pair_generation_modes::hit_and_miss);
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    
    std::size_t n_evts=10000;
    std::size_t n_bins=100;
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing ps_tree MC modules in decays....................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > gamma,gamma");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,2,true>* helgen=uniform_helicities<value_type,1,2,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,2,rn_engine>test(hdecay.get_tree_iterator(),"",isgen_type,psgen_type);
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
	Camgen::log.enable_level=log_level::message;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("gamma > e+,e-");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,2,true>* helgen=uniform_helicities<value_type,1,2,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,2,rn_engine>test(hdecay.get_tree_iterator(),"",isgen_type,psgen_type);
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
	Camgen::log.enable_level=log_level::message;
    }

    {
	std::string process("W- > mu-,nu_mubar");
	std::string fname("plots/W_ln");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>Wdecay(process);
	Wdecay.load();
	Wdecay.construct();
	helicity_generator<value_type,1,2,true>* helgen=uniform_helicities<value_type,1,2,rn_engine,true>::create_instance<model_type>(Wdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,2,rn_engine>test(Wdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

    {
	std::string process("t > b,mu+,nu_mu");
	std::string fname("plots/t_bln");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,3>tdecay(process);
	tdecay.load();
	tdecay.construct();
	helicity_generator<value_type,1,3,true>* helgen=uniform_helicities<value_type,1,3,rn_engine,true>::create_instance<model_type>(tdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,3,rn_engine>test(tdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

    {
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("plots/h_WW_2l2n");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

    {
	model_type::M_h0=140;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("plots/h_WW*_2l2n");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

    {
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > mu-,nu_mubar,mu+,nu_mu");
	std::string fname("plots/h_2l2n2");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

    {
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > mu-,mu+,e-,e+");
	std::string fname("plots/h_ZZ_4l");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

    {
	model_type::M_h0=140;
	model_type::refresh_widths();
	std::string process("h0 > mu-,mu+,e-,e+");
	std::string fname("plots/h_ZZ*_4l");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_m_min(1,2,10);
	test.set_m_min(3,4,10);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing ps_tree MC modules in scattering................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > u,dbar");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
	std::string process("e+,e- > mu+,mu-");
	std::string fname("plots/ee_ll");
	double E1=250;
	double E2=250;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>bahba(process);
	bahba.load();
	bahba.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(bahba.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(bahba.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	if(!test.run(n_evts,n_bins))
	{
	    test.print_generator(std::cerr);
	    std::cerr<<std::endl;
	    test.print_status(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }
	
    {
	std::string process("e+,mu- > nu_ebar,nu_mu");
	std::string fname("plots/emu_2n");
	double E1=250;
	double E2=250;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }
    
    {
	std::string process("e+,e- > nu_e,nu_ebar");
	std::string fname("plots/ee_2n");
	double E1=250;
	double E2=250;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }
    
    {
	std::string process("e+,mu- > nu_ebar,h0,nu_mu");
	std::string fname("plots/emu_h2n");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }

    {
	std::string process("e+,e- > nu_ebar,h0,nu_e");
	std::string fname("plots/ee_h2n");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }
    
    {
	std::string process("nu_e,nu_mu > nu_e,Z,nu_mu");
	std::string fname("plots/nn_Z2n");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }

    {
	std::string process("e+,e- > nu_ebar,mu+,mu-,nu_e");
	std::string fname("plots/ee_2mu2ne");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }

    {
	std::string process("e+,e- > nu_mubar,mu+,mu-,nu_mu");
	std::string fname("plots/ee_2mu2nmu");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }

    {
	std::string process("e+,e- > nu_e,W+,W-,nu_ebar");
	std::string fname("plots/ee_2W2n");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }

    {
	std::string process("e+,e- > nu_e,Z,Z,nu_ebar");
	std::string fname("plots/ee_2Z2n");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }

    {
	std::string process("e+,e- > nu_e,W+,W-,nu_ebar");
	std::string fname("plots/ee_2W2n");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
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
    }

    {
	std::string process("e+,e- > nu_e,u,dbar,mu-,nu_mubar,nu_ebar");
	std::string fname("plots/ee_qqbarl3n");
	double E1=500;
	double E2=500;
	double mmin=10;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,6,true>* helgen=uniform_helicities<value_type,2,6,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	colour_generator<value_type,2,6,false>* colgen=colour_flow_QCD<value_type,2,6,3,rn_engine,false>::create_instance<model_type>(algo.get_tree_iterator());
	colgen->generate();
	ps_generator_tester<model_type,2,6,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_m_min(1,6,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }
}

