//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evtgen_fac_base.h
    \brief Multi-process event generator factory base class.
 */

#ifndef CAMGEN_EVTGEN_FAC_BASE_H_
#define CAMGEN_EVTGEN_FAC_BASE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Base class for event generator factories. Takes a single-process factory  *
 * as engine for creating many subprocess generators.                        *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_gen.h>
#include <Camgen/procgen_fac_base.h>

namespace Camgen
{
    /// Multi-process event generator factory base.
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class event_generator_factory_base
    {
	typedef process_generator_factory_base<model_t,N_in,N_out,rng_t> procgen_fac_type;

	public:

	    /* Public type definitions: */
	    /*--------------------------*/

	    typedef typename procgen_fac_type::process_generator_type process_generator_type;
	    typedef typename procgen_fac_type::momentum_generator_type momentum_generator_type;
	    typedef typename procgen_fac_type::helicity_generator_type helicity_generator_type;
	    typedef typename procgen_fac_type::colour_generator_type colour_generator_type;
	    typedef typename procgen_fac_type::CM_tree_iterator CM_tree_iterator;
	    typedef event_generator<model_t,N_in,N_out,rng_t> event_generator_type;
	    typedef CM_algorithm<model_t,N_in,N_out> amplitude_type;

	    procgen_fac_type* const process_generator_factory;
	    
	    /* Public abstract methods: */
	    /*--------------------------*/

	    /// Constructor
	    
	    event_generator_factory_base(procgen_fac_type* proc_gen_fac):process_generator_factory(proc_gen_fac){}

	    /// Creates a multi-process event generator.

	    event_generator_type* create_generator(amplitude_type& amplitude)
	    {
		event_generator_type* result=new event_generator_type(amplitude);
		configure(result,NULL);
		return result;
	    }

	    /// Creates a multi-process event generator with custom
	    /// configuration.

	    event_generator_type* create_generator(CM_tree_iterator amplitude, generator_configuration<model_t,N_in,N_out>& conf)
	    {
		event_generator_type* result=new event_generator_type(amplitude);
		configure(result,&conf);
		return result;
	    }

	    /// Destructor.

	    virtual ~event_generator_factory_base(){}

	    /// Pre-initialises the argument generator to give a first estimate
	    /// of the cross section.

	    static void pre_initialise(event_generator_type* evtgen,bool verbose=false)
	    {
		if(evtgen!=NULL)
		{
		    evtgen->pre_initialise(pre_init_events(),verbose);
		}
	    }

	    /// Initialises the argument generator.

	    static void initialise(event_generator_type* evtgen,bool verbose=false)
	    {
		if(evtgen!=NULL)
		{
		    evtgen->initialise(init_channel_iterations(),init_channel_batch(),init_grid_iterations(),init_grid_batch(),subprocess_events(),verbose);
		}
	    }

	    /// Initialiser method using generator configuration data.

	    static void initialise(event_generator_type* evtgen,generator_configuration<model_t,N_in,N_out>& conf,bool verbose=false)
	    {
		if(evtgen==NULL)
		{
		    return;
		}

		typename event_generator_type::process_container& procs=evtgen->procs;
		typename event_generator_type::process_iterator& sub_proc=evtgen->sub_proc;

		for(sub_proc=evtgen->procs.begin();sub_proc!=evtgen->procs.end();++sub_proc)
		{
		    if(verbose)
		    {
			std::stringstream ss;
			(*sub_proc)->print_process(ss);
			std::cout<<std::endl<<"init subprocess "<<ss.str()<<std::endl;
		    }
		    conf.configure();
		    procgen_fac_type::initialise(&(*sub_proc),verbose);
		    sub_proc.loop_process(subprocess_events(),verbose);
		}
		evtgen->refresh_cross_section();
		evtgen->adapt_processes();
	    }

	protected:

	    /* Configures the argument generator: */

	    void configure(event_generator_type* evtgen,generator_configuration<model_t,N_in,N_out>* conf)
	    {
		typedef typename model_t::value_type value_type;

		amplitude_type& amplitude=evtgen->algorithm;
		typename event_generator_type::process_container& procs=evtgen->procs;


		procs.clear();
		procs.reserve(amplitude.n_trees());
		if(amplitude.reset_process())
		{
		    int id=1;
		    do
		    {
			if(!amplitude.get_tree_iterator()->is_empty())
			{
                            process_generator_type* proc_gen=process_generator_factory->create_generator(amplitude.get_tree_iterator(),*conf,id++);
			    typename event_generator_type::subprocess_type subproc={proc_gen,1.0};
                            procs.push_back(subproc);
			}
		    }
		    while(amplitude.next_process());


                    //TODO: Get this code to evt-gen
                    evtgen->allocate_event();
		    if(procs.size()!=0)
		    {
			value_type alpha=(value_type)1/(value_type)procs.size();
			for(std::size_t i=0;i<procs.size();++i)
			{
			    procs[i].alpha=alpha;
			}
			evtgen->set_sub_process(procs.begin());
		    }
		    else
		    {
			evtgen->set_sub_process(procs.end());
			return;
		    }
		}
		else
		{
                    evtgen->allocate_event();
                    evtgen->set_sub_process(procs.end());
		    return;
		}
		if(conf!=NULL)
		{
		    conf->configure();
		}
	    }
    };

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>event_generator<model_t,N_in,N_out,rng_t>* pre_initialise(event_generator<model_t,N_in,N_out,rng_t>* evtgen,bool verbose=false)
    {
	event_generator_factory_base<model_t,N_in,N_out,rng_t>::pre_initialise(evtgen,verbose);
	return evtgen;
    }

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>event_generator<model_t,N_in,N_out,rng_t>* initialise(event_generator<model_t,N_in,N_out,rng_t>* evtgen,bool verbose=false)
    {
	event_generator_factory_base<model_t,N_in,N_out,rng_t>::initialise(evtgen,verbose);
	return evtgen;
    }

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>event_generator<model_t,N_in,N_out,rng_t>* initialise(event_generator<model_t,N_in,N_out,rng_t>* evtgen,generator_configuration<model_t,N_in,N_out>& conf,bool verbose=false)
    {
	event_generator_factory_base<model_t,N_in,N_out,rng_t>::initialise(evtgen,conf,verbose);
	return evtgen;
    }
}

#endif /*CAMGEN_PROCGEN_FAC_BASE_H_*/

