//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>

/* * * * * * * * * * * * * * * *
 * Tests for event generators. *
 *                             *
 * * * * * * * * * * * * * * * */

// TODO: add test for auto_update...

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef std::random rn_engine;
    typedef typename model_type::value_type value_type;
    typedef std::size_t size_type;

    license_print::disable();

    initial_states::type isgen_type=initial_states::partonic;
    phase_space_generators::type psgen_type=phase_space_generators::recursive;

    set_initial_state_type(isgen_type);
    set_phase_space_generator_type(psgen_type);
    set_s_pair_generation_mode(s_pair_generation_modes::symmetric);
    std::string fext=plot_config::gnuplot_path==NULL?".dat/.gp":".eps";
    file_utils::create_directory("test_output/evt_gen_test");

    size_type n_init_evts=1000;
    size_type n_evts=10000;
    //    std::size_t n_bins=100;

    set_subprocess_events(n_init_evts);

    std::cout<<"-----------------------------------------------------------"<<std::endl;
    std::cout<<"testing event generator for decays........................."<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("h0 > gamma,gamma");
        std::cerr<<"Checking process generation for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,1,2>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,1,2,rn_engine> factory;
        event_generator<model_type,1,2,rn_engine>* evt_gen=factory.create_generator(algo);
        for(size_type n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
                evt_gen->refresh_cross_section();
            }
        }
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("gamma > l+,l-");
        std::cerr<<"Checking process generation for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,1,2>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,1,2,rn_engine> factory;
        event_generator<model_type,1,2,rn_engine>* evt_gen=factory.create_generator(algo);
        for(size_type n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
                evt_gen->refresh_cross_section();
            }
        }
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("h0 > l+,l-,l+,l-");
        std::cerr<<"Checking event correlation for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,1,4>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
        random_number_stream<model_type::value_type,rn_engine>::reset_engine();
        model_type::value_type w1(0);
        size_type n1(0);
        do
        {
            evt_gen->generate();
            w1=evt_gen->weight();
            n1++;
        }
        while(w1==0);
        random_number_stream<model_type::value_type,rn_engine>::reset_engine();
        model_type::value_type w2(0);
        size_type n2(0);
        do
        {
            evt_gen->generate();
            w2=evt_gen->weight();
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
        std::string process("h0 > l+,l-,l+,l-");
        std::cerr<<"Checking call counting for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,1,4>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_gen->refresh_cross_section();
        }
        size_type calls=0;
        typename event_generator<model_type,1,4,rn_engine>::const_process_iterator it;
        for(it=evt_gen->begin_processes();it!=evt_gen->end_processes();++it)
        {
            calls+=it->generator->calls();
        }
        if(calls!=n_evts)
        {
            return 1;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("h0 > l+,l-,nu,nubar");
        std::cerr<<"Checking process ordering for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,1,4>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_gen->update();
        }
        evt_gen->refresh_cross_section();
        evt_gen->adapt_processes();
        value_type alpha(1);
        typename event_generator<model_type,1,4,rn_engine>::const_process_iterator it;
        for(it=evt_gen->begin_processes();it!=evt_gen->end_processes();++it)
        {
            if(it->alpha>alpha)
            {
                return 1;
            }
            if(it->generator->cross_section().value<evt_gen->process_threshold*evt_gen->cross_section().value)
            {
                return 1;
            }
            alpha=it->alpha;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("h0 > l+,l-,nu,nubar");
        std::cerr<<"Checking process sampling frequency for "<<process<<"..........";
        std::cerr.flush();
        CM_algorithm<model_type,1,4>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,1,4,rn_engine> factory;
        event_generator<model_type,1,4,rn_engine>* evt_gen=factory.create_generator(algo);
        evt_gen->pre_initialise(n_evts);
        for(size_type i=0;i<n_evts;++i)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
                evt_gen->refresh_cross_section();
            }
        }
        for(event_generator<model_type,1,4,rn_engine>::const_process_iterator it=evt_gen->begin_processes();it!=evt_gen->end_processes();++it)
        {
            value_type a_min=it->alpha-MC_integral<value_type>::tolerance/std::sqrt(n_evts);
            value_type a_max=it->alpha+MC_integral<value_type>::tolerance/std::sqrt(n_evts);
            value_type a=((value_type)(it->generator->calls()))/evt_gen->calls();
            if(a<a_min or a>a_max)
            {
                return 1;
            }
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    std::cout<<"-----------------------------------------------------------"<<std::endl;
    std::cout<<"testing event generator for scattering....................."<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("u,dbar > l+,l-");
        value_type E1=500;
        value_type E2=500;
        std::cerr<<"Checking process generation for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,2,2>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,2,2,rn_engine> factory;
        event_generator<model_type,2,2,rn_engine>* evt_gen=factory.create_generator(algo);
        evt_gen->set_beam_energy(-1,E1);
        evt_gen->set_beam_energy(-2,E2);
        for(size_type n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
                evt_gen->refresh_cross_section();
            }
        }
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("q,qbar > t,tbar");
        value_type E1=100;
        value_type E2=100;
        std::cerr<<"Checking process generation for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,2,2>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,2,2,rn_engine> factory;
        event_generator<model_type,2,2,rn_engine>* evt_gen=factory.create_generator(algo);
        evt_gen->set_beam_energy(-1,E1);
        evt_gen->set_beam_energy(-2,E2);
        for(size_type n=0;n<n_evts;++n)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
                evt_gen->refresh_cross_section();
            }
        }
        if(!(evt_gen->cross_section().value==0 and evt_gen->cross_section().error==0))
        {
            return 1;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("e+,e- > h0,nu,nubar");
        value_type E1=250;
        value_type E2=250;
        std::cerr<<"Checking event correlation for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,2,3>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,2,3,rn_engine> factory;
        event_generator<model_type,2,3,rn_engine>* evt_gen=factory.create_generator(algo);
        evt_gen->set_beam_energy(-1,E1);
        evt_gen->set_beam_energy(-2,E2);
        random_number_stream<model_type::value_type,rn_engine>::reset_engine();
        model_type::value_type w1(0);
        size_type n1(0);
        do
        {
            evt_gen->generate();
            w1=evt_gen->weight();
            n1++;
        }
        while(w1==0);
        random_number_stream<model_type::value_type,rn_engine>::reset_engine();
        model_type::value_type w2(0);
        size_type n2(0);
        do
        {
            evt_gen->generate();
            w2=evt_gen->weight();
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
        std::string process("u,ubar > l+,l-,nu,nubar");
        value_type E1=250;
        value_type E2=250;
        std::cerr<<"Checking process ordering for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,2,4>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,2,4,rn_engine> factory;
        event_generator<model_type,2,4,rn_engine>* evt_gen=factory.create_generator(algo);
        evt_gen->set_beam_energy(-1,E1);
        evt_gen->set_beam_energy(-2,E2);
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_gen->update();
            evt_gen->refresh_cross_section();
        }
        evt_gen->adapt_processes();
        value_type alpha(1);
        typename event_generator<model_type,2,4,rn_engine>::const_process_iterator it;
        for(it=evt_gen->begin_processes();it!=evt_gen->end_processes();++it)
        {
            if(it->alpha>alpha)
            {
                return 1;
            }
            if(it->generator->cross_section().value<evt_gen->process_threshold*evt_gen->cross_section().value)
            {
                return 1;
            }
            alpha=it->alpha;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("e+,e- > l+,l-,nu,nubar");
        value_type E1=250;
        value_type E2=250;
        std::cerr<<"Checking process sampling frequency for "<<process<<"..........";
        std::cerr.flush();
        CM_algorithm<model_type,2,4>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,2,4,rn_engine> factory;
        event_generator<model_type,2,4,rn_engine>* evt_gen=factory.create_generator(algo);
        evt_gen->set_beam_energy(-1,E1);
        evt_gen->set_beam_energy(-2,E2);
        evt_gen->pre_initialise(n_evts);
        for(size_type i=0;i<n_evts;++i)
        {
            if(evt_gen->generate())
            {
                evt_gen->update();
                evt_gen->refresh_cross_section();
            }
        }
        for(event_generator<model_type,2,4,rn_engine>::const_process_iterator it=evt_gen->begin_processes();it!=evt_gen->end_processes();++it)
        {
            value_type a_min=it->alpha-MC_integral<value_type>::tolerance/std::sqrt(n_evts);
            value_type a_max=it->alpha+MC_integral<value_type>::tolerance/std::sqrt(n_evts);
            value_type a=((value_type)(it->generator->calls()))/evt_gen->calls();
            if(a<a_min or a>a_max)
            {
                return 1;
            }
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    {
        Camgen::log.enable_level=log_level::error;
        std::string process("e+,e- > l+,l-,nu,nubar");
        value_type E1=250;
        value_type E2=250;
        std::cerr<<"Checking call counting for "<<process<<"............";
        std::cerr.flush();
        CM_algorithm<model_type,2,4>algo(process);
        algo.load();
        algo.construct_trees();
        event_generator_factory<model_type,2,4,rn_engine> factory;
        event_generator<model_type,2,4,rn_engine>* evt_gen=factory.create_generator(algo);
        evt_gen->set_beam_energy(-1,E1);
        evt_gen->set_beam_energy(-2,E2);
        for(size_type i=0;i<n_evts;++i)
        {
            evt_gen->generate();
            evt_gen->refresh_cross_section();
        }
        size_type calls=0;
        typename event_generator<model_type,2,4,rn_engine>::const_process_iterator it;
        for(it=evt_gen->begin_processes();it!=evt_gen->end_processes();++it)
        {
            calls+=it->generator->calls();
        }
        if(calls!=n_evts)
        {
            std::cerr<<std::endl<<calls<<"!="<<n_evts<<std::endl;
            return 1;
        }
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

    return 0;
}
