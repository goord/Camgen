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
    phase_space_generators::type psgen_type=phase_space_generators::recursive_backward_shat;
    set_phase_space_generator_type(psgen_type);
    set_s_pair_generation_mode(s_pair_generation_modes::hit_and_miss);
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    
    std::size_t n_evts=20000;
    std::size_t n_bins=100;

    std::cout<<"--------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing hadronic tree MC with forward s-sampling.........................."<<std::endl;
    std::cout<<"--------------------------------------------------------------------------"<<std::endl;

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
	std::string process("u,dbar > c,sbar");
	std::string fname("plots/qq_qq");
	double E1=500;
	double E2=500;
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
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }


    {
	std::string process("u,ubar > g,g");
	std::string fname("plots/qq_gg");
	double E1=500;
	double E2=500;
	double mmin=20;
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
    }

    {
	std::string process("u,ubar > d,dbar");
	std::string fname("plots/uubar_ddbar");
	double E1=500;
	double E2=500;
	double mmin=20;
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
    }

    {
	std::string process("u,dbar > c,sbar,g");
	std::string fname("plots/qq_qqg");
	double E1=500;
	double E2=500;
	double mmin=20;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>bahba(process);
	bahba.load();
	bahba.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(bahba.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(bahba.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_pT_min(mmin,3);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }

/*    {
	std::string process("g,g > t,tbar,h0");
	std::string fname("plots/gg_tth");
	double E1=500;
	double E2=500;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>bahba(process);
	bahba.load();
	bahba.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(bahba.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(bahba.get_tree_iterator(),fname,isgen_type,psgen_type);
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
	std::string process("u,dbar > c,sbar,mu+,mu-");
	std::string fname("plots/qq_qqg");
	double E1=500;
	double E2=500;
	double mmin=20;
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>bahba(process);
	bahba.load();
	bahba.construct();
	helicity_generator<value_type,2,3,true>* helgen=uniform_helicities<value_type,2,3,rn_engine,true>::create_instance<model_type>(bahba.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,2,3,rn_engine>test(bahba.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.set_beam_energy(-1,E1);
	test.set_beam_energy(-2,E2);
	test.set_pT_min(mmin,3);
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
    }*/
}

