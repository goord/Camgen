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
    set_s_pair_generation_mode(s_pair_generation_modes::hit_and_miss);
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
	std::cerr<<"Checking phase space tree decomposition for "<<process<<"............";
	std::cerr.flush();

	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	process_generator_tester<model_type,1,2,rn_engine> tester(algo.get_tree_iterator(),"");
	if(!tester.run(n_evts,n_bins,false,true,true))
	{
	    std::cerr<<tester.proc_gen->cross_section()<<','<<tester.uni_gen->cross_section()<<std::endl;
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
    

}
