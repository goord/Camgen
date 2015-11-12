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

    initial_states::type isgen_type=initial_states::partonic;
    phase_space_generators::type psgen_type=phase_space_generators::recursive;
    set_phase_space_generator_type(psgen_type);
    set_s_pair_generation_mode(s_pair_generation_modes::hit_and_miss);
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    
    std::size_t n_evts=10000;
    std::size_t n_adapt=100;
    std::size_t n_bins=100;
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing ps_tree multichannel adaptation in decays........................"<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;


    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=165;
	model_type::refresh_widths();
	std::string process("h0 > mu-,nu_mubar,mu+,nu_mu");
	std::string fname("plots/h_WW_2l2n");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.adapt_batch_size=n_adapt;
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	test.print_generator(std::cerr);
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	model_type::M_h0=183;
	model_type::refresh_widths();
	std::string process("h0 > mu-,nu_mubar,mu+,nu_mu");
	std::string fname("plots/h_ZZ_2l2n");
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	helicity_generator<value_type,1,4,true>* helgen=uniform_helicities<value_type,1,4,rn_engine,true>::create_instance<model_type>(algo.get_tree_iterator());
	helgen->generate();
	ps_generator_tester<model_type,1,4,rn_engine>test(algo.get_tree_iterator(),fname,isgen_type,psgen_type);
	test.adapt_batch_size=n_adapt;
	if(!test.run(n_evts,n_bins))
	{
	    return 1;
	}
	test.print_generator(std::cerr);
	std::cerr<<"done, files "<<fname+fext<<" written."<<std::endl;
	delete helgen;
	Camgen::log.enable_level=log_level::warning;
    }
}
