//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/rambo.h>
#include <Camgen/utils.h>
#include <Camgen/psgen_fac.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Tests checking uniform phase space generation for partonic initial states.*
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 

using namespace Camgen;

/* Function checking event generation of the argument generator: */

template<class model_t,std::size_t N_in,std::size_t N_out>bool check_generator(ps_generator<model_t,N_in,N_out>* ps_gen,std::size_t N_evts,bool skip_check=false)
{
    typedef typename ps_generator<model_t,N_in,N_out>::size_type size_type;
    typedef typename ps_generator<model_t,N_in,N_out>::value_type value_type;

    ps_gen->refresh_Ecm();

    for(size_type i=0;i<N_evts;++i)
    {
	bool ps_generated=ps_gen->generate();

	if(ps_generated)
	{
	    if(!skip_check)
	    {
		if(!ps_gen->check())
		{
		    return false;
		}
		value_type w=ps_gen->weight();
		ps_gen->evaluate_weight();
		value_type w2=ps_gen->weight();
		if(!equals(w,w2))
		{
		    return false;
		}
	    }
	}
	else
	{
	    return false;
	}
    }
    return true;
}

int main()
{
    typedef SM model_type;
    typedef std::random rn_engine;
    
    license_print::disable();

    initial_states::type isgen_type=initial_states::partonic;
    phase_space_generators::type psgen_type=phase_space_generators::uniform;
    set_phase_space_generator_type(psgen_type);
    
    std::size_t n_evts=10000;
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing uniform ps generation in decays.................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > gamma,gamma");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,2>* ps_gen=ps_generator_factory<model_type,1,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	if(!check_generator(ps_gen,n_evts,true))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("gamma > e+,e-");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,2>* ps_gen=ps_generator_factory<model_type,1,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	if(check_generator(ps_gen,n_evts,true))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("W- > mu-,nu_mubar");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,2>* ps_gen=ps_generator_factory<model_type,1,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("W- > mu-,nu_mubar");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,2>* ps_gen=ps_generator_factory<model_type,1,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,250.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("t > b,mu+,nu_mu");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,3>* ps_gen=ps_generator_factory<model_type,1,3,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("t > b,mu+,nu_mu");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,3>* ps_gen=ps_generator_factory<model_type,1,3,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,750.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,4>* ps_gen=ps_generator_factory<model_type,1,4,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,4>* ps_gen=ps_generator_factory<model_type,1,4,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,1750.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu,Z");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,5>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,5>* ps_gen=ps_generator_factory<model_type,1,5,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu,Z");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,5>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,5>* ps_gen=ps_generator_factory<model_type,1,5,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,1750.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu,b,bbar");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,6>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,6>* ps_gen=ps_generator_factory<model_type,1,6,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu,b,bbar");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,6>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,1,6>* ps_gen=ps_generator_factory<model_type,1,6,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,1750.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing uniform ps generation in scattering.............................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("e+,e- > u,dbar");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* ps_gen=ps_generator_factory<model_type,2,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,100.0);
	ps_gen->set_beam_energy(-2,100.0);
	if(!check_generator(ps_gen,n_evts,true))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > t,tbar");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* ps_gen=ps_generator_factory<model_type,2,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,100.0);
	ps_gen->set_beam_energy(-2,100.0);
	if(check_generator(ps_gen,n_evts,true))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* ps_gen=ps_generator_factory<model_type,2,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,100.0);
	ps_gen->set_beam_energy(-2,100.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* ps_gen=ps_generator_factory<model_type,2,2,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,100.0);
	ps_gen->set_beam_energy(-2,500.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,g");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* ps_gen=ps_generator_factory<model_type,2,3,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,100.0);
	ps_gen->set_beam_energy(-2,100.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,g");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* ps_gen=ps_generator_factory<model_type,2,3,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,500.0);
	ps_gen->set_beam_energy(-2,100.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,g,g");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* ps_gen=ps_generator_factory<model_type,2,4,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,100.0);
	ps_gen->set_beam_energy(-2,100.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,g,g");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* ps_gen=ps_generator_factory<model_type,2,4,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,100.0);
	ps_gen->set_beam_energy(-2,500.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,t,tbar");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* ps_gen=ps_generator_factory<model_type,2,4,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,500.0);
	ps_gen->set_beam_energy(-2,500.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,t,tbar");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* ps_gen=ps_generator_factory<model_type,2,4,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,1000.0);
	ps_gen->set_beam_energy(-2,500.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,t,tbar,g");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,5>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,5>* ps_gen=ps_generator_factory<model_type,2,5,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,500.0);
	ps_gen->set_beam_energy(-2,500.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > g,g,t,tbar,g");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,5>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,5>* ps_gen=ps_generator_factory<model_type,2,5,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,500.0);
	ps_gen->set_beam_energy(-2,1000.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > t,tbar,t,tbar,g");
	std::cerr<<"Checking uniform phase space generation for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,5>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,5>* ps_gen=ps_generator_factory<model_type,2,5,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,500.0);
	ps_gen->set_beam_energy(-2,500.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }

    {
	Camgen::log.enable_level=log_level::error;
	std::string process("g,g > t,tbar,t,tbar,g");
	std::cerr<<"Checking uniform phase space generation for boosted "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,2,5>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,5>* ps_gen=ps_generator_factory<model_type,2,5,rn_engine>::create_generator(algo.get_tree_iterator(),isgen_type,psgen_type);
	ps_gen->set_beam_energy(-1,1500.0);
	ps_gen->set_beam_energy(-2,500.0);
	if(!check_generator(ps_gen,n_evts))
	{
	    return 1;
	}
	std::cerr<<"done."<<std::endl;
	Camgen::log.enable_level=log_level::warning;
    }
}

