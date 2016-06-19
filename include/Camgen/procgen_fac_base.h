//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file procgen_fac_base.h
    \brief Single-process event generator factory abstract base class.
 */

#ifndef CAMGEN_PROCGEN_FAC_BASE_H_
#define CAMGEN_PROCGEN_FAC_BASE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for event generator factories. Child classes should   *
 * implement creation methods for phase space, helicity and colour generator *
 * instances.                                                                *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/MC_config.h>
#include <Camgen/proc_gen.h>

namespace Camgen
{
    /// Single-process event generator factory base.
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class process_generator_factory_base
    {
	public:

	    /* Public type definitions: */
	    /*--------------------------*/

	    typedef process_generator<model_t,N_in,N_out,rng_t> process_generator_type;
	    typedef typename process_generator_type::momentum_generator_type momentum_generator_type;
	    typedef typename process_generator_type::helicity_generator_type helicity_generator_type;
	    typedef typename process_generator_type::colour_generator_type colour_generator_type;
	    typedef typename process_generator_type::CM_tree_iterator CM_tree_iterator;

	    /* Public abstract methods: */
	    /*--------------------------*/

	    /// Creates a phase space generator for the given single-process amplitude.

	    virtual momentum_generator_type* create_momentum_generator(CM_tree_iterator amplitude)=0;

	    /// Creates a helicity generator for the given single-process
	    /// amplitude.
	    
	    virtual helicity_generator_type* create_helicity_generator(CM_tree_iterator amplitude)=0;
	    
	    /// Creates a colour generator for the given single-process
	    /// amplitude.
	    
	    virtual colour_generator_type* create_colour_generator(CM_tree_iterator amplitude)=0;

	    /// Virtual destructor

	    virtual ~process_generator_factory_base(){}

	    /// Creates a single process generator.

	    process_generator_type* create_generator(CM_tree_iterator amplitude, int id=0)
	    {
		return create_generator(amplitude,NULL,id);
	    }

	    /// Creates a single process generator with custom configuration.

	    process_generator_type* create_generator(CM_tree_iterator amplitude, generator_configuration<model_t,N_in,N_out>& conf, int id=0)
	    {
		return create_generator(amplitude,&conf,id);
	    }

	    /// Pre-initialises the argument generator to give a first estimate
	    /// of the cross section.

	    static void pre_initialise(process_generator_type* proc_gen,bool verbose=false)
	    {
		proc_gen->pre_initialise(pre_init_events(),verbose);
	    }

	    /// Initialises with number of channel and grid iterations/batch
	    /// sizes defined in the static configuration.
	    
	    static void initialise(process_generator_type* proc_gen,bool verbose=false)
	    {
		proc_gen->initialise(init_channel_iterations(),init_channel_batch(),init_grid_iterations(),init_grid_batch(),verbose);
	    }

	    /// Initialises with the number of channel and grid iterations/batch
	    /// sizes defined in the configuration object.

	    static void initialise(process_generator_type* procgen,generator_configuration<model_t,N_in,N_out>& conf,bool verbose=false)
	    {
		conf.configure(procgen->get_process());
		initialise(procgen,verbose);
	    }

	protected:

	    /* Utility method creating a process generator with configuration
	     * pointer argument: */

	    process_generator_type* create_generator(CM_tree_iterator amplitude, generator_configuration<model_t,N_in,N_out>* conf, int id=0)
	    {
		process_generator_type* result=new process_generator_type(amplitude,id);
                if(result!=NULL)
                {
                    result->allocate_event();
                    configure(result,conf);
                }
		return result;
	    }

	    /* Adopts parameters and generator instances from the configuration object. 
	     * Sets the phase space, helicity and colour generators to the given arguments: */
	    
	    void configure(process_generator_type* procgen,generator_configuration<model_t,N_in,N_out>* conf)
	    {
		if(conf!=NULL)
		{
		    basic_cuts::clear();
		    conf->configure(procgen->get_process());
		}

		procgen->set_ps_generator(create_momentum_generator(procgen->amplitude));
		procgen->set_helicity_generator(create_helicity_generator(procgen->amplitude));
		procgen->set_colour_generator(create_colour_generator(procgen->amplitude));

		if(Camgen::auto_channel_adapt())
		{
		    procgen->set_auto_channel_adapt(auto_channel_batch());
		}
		if(Camgen::auto_grid_adapt())
		{
		    procgen->set_auto_grid_adapt(auto_grid_batch());
		}
		procgen->max_rejects=max_init_rejects();
		if(discarded_weight_fraction()>0 and discarded_weight_fraction()<1)
		{
		    procgen->bin_weights(weight_bins());
		    procgen->reduce_max_weight(discarded_weight_fraction());
		}
		for(std::size_t i=1;i<=N_out;++i)
		{
		    for(std::size_t j=i+1;j<=N_out;++j)
		    {
			procgen->set_m_min(i,j,basic_cuts::m_min(i,j));
		    }
		}
		procgen->set_pdf_alpha_s(use_pdf_alpha_s());
	    }
    };

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>process_generator<model_t,N_in,N_out,rng_t>* pre_initialise(process_generator<model_t,N_in,N_out,rng_t>* procgen,bool verbose=false)
    {
	if(procgen==NULL)
	{
	    return NULL;
	}
	process_generator_factory_base<model_t,N_in,N_out,rng_t>::pre_initialise(procgen,verbose);
	return procgen;
    }

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>process_generator<model_t,N_in,N_out,rng_t>* initialise(process_generator<model_t,N_in,N_out,rng_t>* procgen,bool verbose=false)
    {
	process_generator_factory_base<model_t,N_in,N_out,rng_t>::initialise(procgen,verbose);
	return procgen;
    }

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>process_generator<model_t,N_in,N_out,rng_t>* initialise(process_generator<model_t,N_in,N_out,rng_t>* procgen,generator_configuration<model_t,N_in,N_out>& conf,bool verbose=false)
    {
	process_generator_factory_base<model_t,N_in,N_out,rng_t>::initialise(procgen,conf,verbose);
	return procgen;
    }
}

#endif /*CAMGEN_PROCGEN_FAC_BASE_H_*/

