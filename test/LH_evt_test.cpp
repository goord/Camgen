//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>
#include <Camgen/LH_evt_stream.h>

/* * * * * * * * * * * * * * * * * * * *
 * Tests for Les-Houches event record. *
 *                                     *
 * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef std::random rn_engine;
    typedef std::size_t size_type;
    typedef double value_type;

    license_print::disable();

    file_utils::create_directory("test_output/LH_evt_test");
    size_type n_evts=100;

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > u,ubar,d,dbar");
	std::string fname("test_output/LH_evt_test/h_uudd");
	std::cerr<<"Checking Les Houches event record for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct_trees();
        process_generator_factory<model_type,1,4,rn_engine> factory;
        process_generator<model_type,1,4,rn_engine>* proc_gen=factory.create_generator(algo.get_tree_iterator());
        event_stream<model_type,1,4>* lh_if=new LH_event_stream<model_type,1,4>(fname,1);
        for(size_type i=0;i<n_evts;++i)
        {
            proc_gen->generate();
            lh_if->fill(proc_gen->get_event());
        }
        lh_if->write();
        std::cerr<<"done, file "<<fname<<".LHE written."<<std::endl;
        delete lh_if;
        delete proc_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	model_type::M_h0=170;
	model_type::refresh_widths();
	std::string process("h0 > q,q,qbar,qbar");
	std::string fname("test_output/LH_evt_test/h_4q");
	std::cerr<<"Checking Les Houches event record for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct_trees();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
        event_stream<model_type,1,4>* lh_if=new LH_event_stream<model_type,1,4>(fname,2);
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            lh_if->fill(evt_gen->get_event());
        }
        lh_if->write();
        std::cerr<<"done, file "<<fname<<".LHE written."<<std::endl;
        delete lh_if;
        delete evt_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	std::string process("e+,e- > q,qbar,Z");
	std::string fname("test_output/LH_evt_test/ee_qqZ");
        value_type E1=100;
        value_type E2=100;
	std::cerr<<"Checking Les Houches event record for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct_trees();
        set_beam_energy(-1,E1);
        set_beam_energy(-2,E2);
        event_generator_factory<model_type,2,3,rn_engine> factory;
        event_generator<model_type,2,3,rn_engine>* evt_gen=factory.create_generator(algo);
        event_stream<model_type,2,3>* lh_if=new LH_event_stream<model_type,2,3>(fname,2);
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            lh_if->fill(evt_gen->get_event());
        }
        lh_if->write();
        std::cerr<<"done, file "<<fname<<".LHE written."<<std::endl;
        delete lh_if;
        delete evt_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::proton_proton);
        set_phase_space_generator_type(phase_space_generators::recursive);
	std::string process("p,p > q,qbar,t,tbar");
	std::string fname("test_output/LH_evt_test/pp_qqtt");
        value_type E1=1000;
        value_type E2=1000;
	std::cerr<<"Checking Les Houches event record for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct_trees();
        set_beam_energy(-1,E1);
        set_beam_energy(-2,E2);
        event_generator_factory<model_type,2,4,rn_engine> factory;
        event_generator<model_type,2,4,rn_engine>* evt_gen=factory.create_generator(algo);
        event_stream<model_type,2,4>* lh_if=new LH_event_stream<model_type,2,4>(fname,2);
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            lh_if->fill(evt_gen->get_event());
        }
        lh_if->write();
        std::cerr<<"done, file "<<fname<<".LHE written."<<std::endl;
        delete lh_if;
        delete evt_gen;
        Camgen::log.enable_level=log_level::warning;
    }
}


