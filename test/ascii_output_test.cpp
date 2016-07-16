//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <config.h>
#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>
#include <Camgen/ascii_file.h>
#include <Camgen/proc_split_if.h>
#include <Camgen/evt_ostream.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Tests for ROOT interface to event/process generators  *
 *                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

template<class model_t,std::size_t N_in,std::size_t N_out>class test_output: public event_output_configuration<model_t,N_in,N_out>
{
    public:

        typedef event_output_configuration<model_t,N_in,N_out> base_type;
        typedef typename base_type::event_type event_type;
        typedef typename base_type::size_type size_type;
        typedef typename base_type::value_type value_type;
        typedef typename base_type::momentum_type momentum_type;

        momentum_type p2;
        value_type alpha12;

        event_output_configuration<model_t,N_in,N_out>* clone() const
        {
            return new test_output(*this);
        }

        void fill(const event_type& evt)
        {
            p2=evt.p(2);
            alpha12=evt.alpha(1,2);
        }

    protected:

        void add_variables()
        {
            this->add_variable(p2,"p2");
            this->add_variable(alpha12,"alpha12");
        }
};

template<class model_t,std::size_t N_out>class test_output<model_t,2,N_out>: public event_output_configuration<model_t,2,N_out>
{
    public:

        typedef event_output_configuration<model_t,2,N_out> base_type;
        typedef typename base_type::event_type event_type;
        typedef typename base_type::size_type size_type;
        typedef typename base_type::value_type value_type;
        typedef typename base_type::momentum_type momentum_type;

        momentum_type p2;
        value_type pT1;
        value_type eta1;
        value_type alpha12;

        event_output_configuration<model_t,2,N_out>* clone() const
        {
            return new test_output(*this);
        }

        void fill(const event_type& evt)
        {
            p2=evt.p(2);
            pT1=evt.pT(1);
            eta1=evt.eta(1);
            alpha12=evt.alpha(1,2);
        }

    protected:

        void add_variables()
        {
            this->add_variable(p2,"p2");
            this->add_variable(pT1,"pT1");
            this->add_variable(eta1,"eta1");
            this->add_variable(alpha12,"alpha12");
        }
};

int main()
{
    typedef SM model_type;
    typedef std::random rn_engine;
    typedef std::size_t size_type;
    typedef double value_type;

    license_print::disable();

    file_utils::create_directory("test_output/ascii_test");
    size_type n_evts=1000;

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_output/ascii_test/h_WW_2l2n");
	std::cerr<<"Checking ascii interface for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
        process_generator_factory<model_type,1,4,rn_engine> factory;
        process_generator<model_type,1,4,rn_engine>* proc_gen=factory.create_generator(algo.get_tree_iterator());
        event_output_stream<model_type,1,4>* evt_os=new event_output_stream<model_type,1,4>(new ascii_file<model_type,1,4>(fname,"test tree"),new test_output<model_type,1,4>());
        for(size_type i=0;i<n_evts;++i)
        {
            proc_gen->generate();
            evt_os->fill(proc_gen->get_event());
        }
        evt_os->write_statistics();
        evt_os->write();
        std::cerr<<"done, file "<<fname<<".dat written."<<std::endl;
        delete evt_os;
        delete proc_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	model_type::M_h0=170;
	model_type::refresh_widths();
	std::string process("h0 > l+,l-,l+,l-");
	std::string fname("test_output/ascii_test/h_WW_4l");
	std::cerr<<"Checking ascii interface for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
        event_output_stream<model_type,1,4>* evt_os=new event_output_stream<model_type,1,4>(new ascii_file<model_type,1,4>(fname,"test tree"),new test_output<model_type,1,4>());
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_os->fill(evt_gen->get_event());
        }
        evt_os->write_statistics();
        evt_os->write();
        std::cerr<<"done, file "<<fname<<".dat written."<<std::endl;
        delete evt_os;
        delete evt_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	std::string process("u,ubar > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_output/ascii_test/uubar_2l2n");
        value_type E1=100;
        value_type E2=100;
	std::cerr<<"Checking ascii interface for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct_trees();
        set_beam_energy(-1,E1);
        set_beam_energy(-2,E2);
        process_generator_factory<model_type,2,4,rn_engine> factory;
        process_generator<model_type,2,4,rn_engine>* proc_gen=factory.create_generator(algo.get_tree_iterator());
        event_output_stream<model_type,2,4>* evt_os=new event_output_stream<model_type,2,4>(new ascii_file<model_type,2,4>(fname,"test tree"),new test_output<model_type,2,4>());
        for(size_type i=0;i<n_evts;++i)
        {
            proc_gen->generate();
            evt_os->fill(proc_gen->get_event());
        }
        evt_os->write_statistics();
        evt_os->write();
        std::cerr<<"done, file "<<fname<<".dat written."<<std::endl;
        delete evt_os;
        delete proc_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	std::string process("p,p > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_output/ascii_test/pp_2l2n");
        value_type E1=300;
        value_type E2=300;
	std::cerr<<"Checking ascii interface for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
        set_beam_energy(-1,E1);
        set_beam_energy(-2,E2);
        event_generator_factory<model_type,2,4,rn_engine> factory;
        event_generator<model_type,2,4,rn_engine>* evt_gen=factory.create_generator(algo);
        event_output_stream<model_type,2,4>* evt_os=new event_output_stream<model_type,2,4>(new ascii_file<model_type,2,4>(fname,"test tree"),new test_output<model_type,2,4>());
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_os->fill(evt_gen->get_event());
        }
        evt_os->write_statistics();
        evt_os->write();
        std::cerr<<"done, file "<<fname<<".dat written."<<std::endl;
        delete evt_os;
        delete evt_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	model_type::M_h0=170;
	model_type::refresh_widths();
	std::string process("h0 > l+,l-,nu,nubar");
	std::string fname("test_output/ascii_test/h_WW_4l");
	std::cerr<<"Checking process-split ascii interface for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct_trees();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
        event_stream<model_type,1,4>* evt_os=new process_split_interface<model_type,1,4>(evt_gen,new ascii_file<model_type,1,4>(fname,"test tree"),new test_output<model_type,1,4>());
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_os->fill(evt_gen->get_event());
        }
        evt_os->write();
        std::cerr<<"done, files "<<fname<<"_x.dat written."<<std::endl;
        delete evt_os;
        delete evt_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
        set_initial_state_type(initial_states::partonic);
        set_phase_space_generator_type(phase_space_generators::recursive);
	std::string process("e+,e- > l+,l-,nu,nubar");
	std::string fname("test_output/ascii_test/ee_WW_4l");
	std::cerr<<"Checking process-split ascii interface for "<<process<<"............";
        value_type E1=100;
        value_type E2=100;
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct_trees();
        set_beam_energy(-1,E1);
        set_beam_energy(-2,E2);
        event_generator_factory<model_type,2,4,rn_engine> factory;
        event_generator<model_type,2,4,rn_engine>* evt_gen=factory.create_generator(algo);
        event_stream<model_type,2,4>* evt_os=new process_split_interface<model_type,2,4>(evt_gen,new ascii_file<model_type,2,4>(fname,"test tree"),new test_output<model_type,2,4>());
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_os->fill(evt_gen->get_event());
        }
        evt_os->write();
        std::cerr<<"done, files "<<fname<<"_x.dat written."<<std::endl;
        delete evt_os;
        delete evt_gen;
        Camgen::log.enable_level=log_level::warning;
    }

    return 0;
}
