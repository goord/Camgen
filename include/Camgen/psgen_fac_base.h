//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file psgen_fac_base.h
    \brief Base class for phase space generator factories.
 */

#ifndef CAMGEN_PSGEN_FAC_BASE_H_
#define CAMGEN_PSGEN_FAC_BASE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for phase space generator factories. Implementors must  *
 * provide a phase space generator creation method from the runtime            *
 * configuration or from a fiel stream.                                        *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/MC_config.h>
#include <Camgen/ps_gen.h>

namespace Camgen
{
    /// Base class for phase space factories. Last template argument denotes the
    /// derived type (static polymorphism)
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class is_factory_type,class fs_factory_type>class ps_generator_factory_base
    {
	public:

	    /* Type definitions: */

	    typedef ps_generator<model_t,N_in,N_out> generator_type;
	    typedef initial_state<model_t,N_in> is_generator_type;
	    typedef initial_states::type is_type;
	    typedef phase_space_generators::type fs_type;
	    typedef typename CM_algorithm<model_t,N_in,N_out>::tree_iterator tree_iterator;

	    /// Creates a phase space generator instance for given amplitude
	    /// iterator using runtime configuration data.

	    static generator_type* create_generator(tree_iterator it)
	    {
		return create_generator(it,initial_state_type(),phase_space_generator_type());
	    }

	    /// Creates a phase space generator for the given amplitude
	    /// iterator, initial state type and final state type.

	    static generator_type* create_generator(tree_iterator it,is_type is_type_,fs_type fs_type_)
	    {
		bool backward_shat=(fs_type_==phase_space_generators::recursive_backward_shat);
		is_generator_type* isgen=is_factory_type::create_instance(is_type_,backward_shat);
		if(isgen==NULL)
		{
                    log(log_level::error)<<CAMGEN_STREAMLOC<<"initial state generator could not be constructed -- returning NULL"<<endlog;
		    return NULL;
		}
		generator_type* result=fs_factory_type::create_instance(isgen,fs_type_);
		if(result!=NULL)
		{
                    result->allocate_event();
		    if(!result->set_amplitude(it))
		    {
			return NULL;
		    }
		    for(int i=1;i<=(int)N_in;++i)
		    {
			result->set_beam_energy(-i,beam_energy(i));
		    }
		    result->refresh_Ecm();
		    for(int i=1;i<=(int)N_out;++i)
		    {
			for(int j=i+1;j<=(int)N_out;++j)
			{
			    result->set_m_min(i,j,basic_cuts::m_min(i,j));
			}
		    }
		    result->refresh_m_min();
		}
		return result;
	    }
    };
}

#endif /*CAMGEN_PSGEN_FAC_BASE_H_*/

