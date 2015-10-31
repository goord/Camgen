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

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef SM::value_type value_type;
    typedef std::random rn_engine;
    
    license_print::disable();

    initial_states::type isgen_type=initial_states::proton_proton;
    phase_space_generators::type psgen_type;
    set_s_pair_generation_mode(s_pair_generation_modes::asymmetric);
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    
    std::size_t n_evts=20000;
    std::size_t n_bins=100;

    std::cout<<"--------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing hadronic tree MC with forward shat/y-sampling....................."<<std::endl;
    std::cout<<"--------------------------------------------------------------------------"<<std::endl;

    psgen_type=phase_space_generators::recursive;
    set_phase_space_generator_type(psgen_type);

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,dbar > e+,e-");
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
	std::string process("u,ubar > t,tbar");
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
	Camgen::log.enable_level=log_level::error;
	std::string process("u,dbar > c,sbar");
	std::string fname("plots/qq_qq");
	double E1=500;
	double E2=500;
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
	Camgen::log.enable_level=log_level::warning;
    }


    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,ubar > g,g");
	std::string fname("plots/qq_gg");
	double E1=500;
	double E2=500;
	double mmin=20;
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
	test.set_m_min(1,2,mmin);
	test.set_pT_min(mmin);
	if(!test.run(n_evts,n_bins))
	{
	    test.print_generator(std::cerr);
	    std::cerr<<std::endl;
	    test.print_status(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,ubar > d,dbar");
	std::string fname("plots/uubar_ddbar");
	double E1=500;
	double E2=500;
	double mmin=20;
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
	test.set_m_min(1,2,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    test.print_generator(std::cerr);
	    std::cerr<<std::endl;
	    test.print_status(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,dbar > c,sbar,g");
	std::string fname("plots/qq_qqg");
	double E1=500;
	double E2=500;
	double mmin=20;
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
	test.set_pT_min(mmin,3);
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
	std::string process("u,dbar > c,sbar,mu+,mu-");
	std::string fname("plots/qq_qqll");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	double mmin=20;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_pT_min(mmin,3);
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
	std::string process("g,g > t,tbar,h0");
	std::string fname("plots/gg_tth");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    test.print_generator(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > t,tbar,b,bbar");
	std::string fname("plots/gg_ttbb");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
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
	std::string process("g,g > t,tbar,e-,nu_ebar,mu+,nu_mu");
	std::string fname("plots/gg_ttllnn");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,6,true>* helgen=uniform_helicities<value_type,2,6,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,6,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
    std::cout<<"--------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing hadronic tree MC with forward xx-sampling........................."<<std::endl;
    std::cout<<"--------------------------------------------------------------------------"<<std::endl;

    isgen_type=initial_states::proton_proton_xx;


    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,dbar > c,sbar");
	std::string fname("plots/qqxx_qq");
	double E1=500;
	double E2=500;
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
	Camgen::log.enable_level=log_level::warning;
    }


    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,ubar > g,g");
	std::string fname("plots/qqxx_gg");
	double E1=500;
	double E2=500;
	double mmin=20;
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
	test.set_m_min(1,2,mmin);
	test.set_pT_min(mmin);
	if(!test.run(n_evts,n_bins))
	{
	    test.print_generator(std::cerr);
	    std::cerr<<std::endl;
	    test.print_status(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,ubar > d,dbar");
	std::string fname("plots/uubarxx_ddbar");
	double E1=500;
	double E2=500;
	double mmin=20;
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
	test.set_m_min(1,2,mmin);
	if(!test.run(n_evts,n_bins))
	{
	    test.print_generator(std::cerr);
	    std::cerr<<std::endl;
	    test.print_status(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("u,dbar > c,sbar,g");
	std::string fname("plots/qqxx_qqg");
	double E1=500;
	double E2=500;
	double mmin=20;
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
	test.set_pT_min(mmin,3);
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
	std::string process("u,dbar > c,sbar,mu+,mu-");
	std::string fname("plots/qqxx_qqll");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	double mmin=20;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_pT_min(mmin,3);
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
	std::string process("g,g > t,tbar,h0");
	std::string fname("plots/ggxx_tth");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    test.print_generator(std::cerr);
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > t,tbar,b,bbar");
	std::string fname("plots/ggxx_ttbb");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,4,true>* helgen=uniform_helicities<value_type,2,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
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
	std::string process("g,g > t,tbar,e-,nu_ebar,mu+,nu_mu");
	std::string fname("plots/ggxx_ttllnn");
	set_beam_energy(-1,500);
	set_beam_energy(-2,500);
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,2,6,true>* helgen=uniform_helicities<value_type,2,6,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,6,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
}

