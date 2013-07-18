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
#include <ps_gen_tester.h>
#include <SMPbKsch.h>

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef SM::value_type value_type;
    typedef std::random rn_engine;
    
    license_print::disable();
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing ps_tree MC modules in decays....................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    initial_states::type isgen_type=initial_states::partonic;
    phase_space_generators::type psgen_type=phase_space_generators::recursive;
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    
    {
	std::string process("W- > mu-,nu_mubar");
	std::string fname("plots/Wdecay");
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,1,2>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,2,true>* helgen=uniform_helicities<value_type,1,2,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,2,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.run(10000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }

    {
	std::string process("t > b,mu+,nu_mu");
	std::string fname("plots/tdecay");
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<SMPbKsch,1,3>tdecay(process);
	tdecay.load();
	tdecay.construct();
	helicity_generator<value_type,1,3,true>* helgen=uniform_helicities<value_type,1,3,rn_engine,true>::create_instance<SMPbKsch>(tdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<SMPbKsch,1,3,rn_engine>test(tdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.run(10000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }

    {
	SM::M_h0=180;
	SM::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("plots/hdecay");
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,1,4>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.run(10000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }

    {
	model_type::M_h0=200;
	std::string process("h0 > mu-,nu_mubar,mu+,nu_mu");
	std::string fname("plots/hdecay2");
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,1,4>hdecay(process);
	hdecay.load();
	hdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(hdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(hdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.run(50000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }

    {
	double mmin=1;
	std::string process("W+ > nu_e,nu_ebar,mu+,nu_mu");
	std::string fname("plots/Wdecay2");
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,1,4>Wdecay(process);
	Wdecay.load();
	Wdecay.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(Wdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(Wdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_m_min(1,2,mmin);
	test.set_m_min(3,4,mmin);
	test.run(300000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }

    {
	model_type::W_Z=100.;
	model_type::W_W=100.;
	double mmin=25;
	std::string process("Z > e+,e-,mu+,nu_mu,tau-,nu_taubar");
	std::string fname("plots/Zdecay");
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,1,6>Zdecay(process);
	Zdecay.load();
	Zdecay.construct();
	value_type E_beam=1000;
	helicity_generator<value_type,1,6,true>* helgen=uniform_helicities<value_type,1,6,rn_engine,true>::create_instance<model_type>(Zdecay.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,6,rn_engine>test(Zdecay.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E_beam);
	test.set_m_min(1,2,mmin);
	test.set_m_min(3,4,mmin);
	test.set_m_min(3,6,mmin);
	test.run(100000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing ps_tree MC modules in scattering................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    
    {
	std::string process("e+,e- > mu+,mu-");
	std::string fname("plots/bahbaZ");
	double E1=250;
	double E2=250;
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,2,2>bahba(process);
	bahba.load();
	bahba.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(bahba.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(bahba.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.run(10000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }
	
    {
	std::string process("e+,mu- > nu_ebar,nu_mu");
	std::string fname("plots/tWchannel");
	double E1=250;
	double E2=250;
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.run(10000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }
    
    {
	std::string process("e+,e- > nu_e,nu_ebar");
	std::string fname("plots/sZuWchannel");
	double E1=250;
	double E2=250;
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,2,true>* helgen=uniform_helicities<value_type,2,2,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,2,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.run(10000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }
    
    {
	std::string process("e+,mu- > nu_ebar,h0,nu_mu");
	std::string fname("plots/Hprod1");
	double E1=500;
	double E2=500;
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.run(20000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }

    {
	std::string process("e+,e- > nu_ebar,h0,nu_e");
	std::string fname("plots/Hprod2");
	double E1=500;
	double E2=500;
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.run(20000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }

    {
	std::string process("e+,e- > W+,Z,W-");
	std::string fname("plots/WWZprod");
	double E1=500;
	double E2=500;
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.run(20000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }
    
    {
	std::string process("e+,e- > nu_ebar,W+,W-,nu_e");
	std::string fname("plots/WWnunuprod");
	double E1=500;
	double E2=500;
	double mmin=1;
	std::cout<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_m_min(1,4,mmin);
	test.run(20000,100);
	std::cout<<"done, files "<<fname+fext<<" written."<<std::endl;
    }
}

