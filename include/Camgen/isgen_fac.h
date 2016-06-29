//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_ISGEN_FAC_H_
#define CAMGEN_ISGEN_FAC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Initial state factory specialisation for the istream class, creating an *
 * initial state instance from a filestream.                               *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/MC_config.h>
#include <Camgen/part_is.h>
#include <Camgen/had_is.h>

namespace Camgen
{
    /// Initial state factory (not specialised).

    template<class model_t,std::size_t N,class rng_t,class spacetime_t=typename model_t::spacetime_type>class initial_state_factory
    {
	public:

	    /* Type definitions: */

	    typedef initial_state<model_t,N> generator_type;
	    typedef initial_states::type initial_state_type;

	    /// Factory method creating an initial state generator instance from
	    /// the static configuration data.

	    static generator_type* create_instance()
	    {
		bool q=true;
		return create_instance(initial_state_type(),q);
	    }

	    /// Factory method from initial state type. Returns NULL.

	    static generator_type* create_instance(initial_state_type is_type,bool& backward_shat)
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"initial state type not defined for current spacetime type--returning NULL"<<endlog;
		return NULL;
	    }
    };

    /// Initial state factory specialisation for decays in Minkowski spacetimes.

    template<class model_t,class rng_t>class initial_state_factory<model_t,1,rng_t,Minkowski_type>
    {
	public:

	    /* Type definitions: */

	    typedef initial_state<model_t,1> generator_type;
	    typedef initial_states::type initial_state_type;

	    /// Factory method creating an initial state generator instance from
	    /// the static configuration data.

	    static generator_type* create_instance()
	    {
		return create_instance(initial_state_type(),true);
	    }

	    /// Factory method from initial state type. Returns NULL unless the
	    /// type is partonic.

	    static generator_type* create_instance(initial_state_type is_type,bool backward_shat=false)
	    {
		if(is_type==initial_states::partonic and !backward_shat)
		{
		    return new partonic_is<model_t,1,Minkowski_type>;
		}
		log(log_level::error)<<CAMGEN_STREAMLOC<<"initial state type "<<is_type<<" is not a valid single-particle initial state"<<endlog;
		return NULL;
	    }
    };

    /// Initial state factory specialisation for 2-particle initial states in
    /// Minkowski spacetimes.

    template<class model_t,class rng_t>class initial_state_factory<model_t,2,rng_t,Minkowski_type>
    {
	public:

	    typedef initial_state<model_t,2> generator_type;
	    typedef initial_states::type initial_state_type;

	    /// Factory method creating an initial state generator instance from
	    /// the static configuration data.

	    static generator_type* create_instance()
	    {
		return create_instance(initial_state_type(),backward_shat_sampling());
	    }

	    /// Factory method with initial state type argument.

	    static generator_type* create_instance(initial_state_type is_type,bool backward_shat)
	    {
		switch(is_type)
		{
		    case initial_states::partonic:
			return new partonic_is<model_t,2,Minkowski_type>;
		    case initial_states::eplus_eminus:
			return new partonic_is<model_t,2,Minkowski_type>;
		    case initial_states::proton_proton:
			if(backward_shat)
			{
			    return new hadronic_is_y<model_t,rng_t,Minkowski_type>(false,false);
			}
			else
			{
			    return new hadronic_is_sy<model_t,rng_t,Minkowski_type>(false,false);
			}
		    case initial_states::proton_antiproton:
			if(backward_shat)
			{
			    return new hadronic_is_y<model_t,rng_t,Minkowski_type>(false,true);
			}
			else
			{
			    return new hadronic_is_sy<model_t,rng_t,Minkowski_type>(false,true);
			}
		    case initial_states::antiproton_proton:
			if(backward_shat)
			{
			    return new hadronic_is_y<model_t,rng_t,Minkowski_type>(true,false);
			}
			else
			{
			    return new hadronic_is_sy<model_t,rng_t,Minkowski_type>(true,false);
			}
		    case initial_states::antiproton_antiproton:
			if(backward_shat)
			{
			    return new hadronic_is_y<model_t,rng_t,Minkowski_type>(true,true);
			}
			else
			{
			    return new hadronic_is_sy<model_t,rng_t,Minkowski_type>(true,true);
			}
		    case initial_states::proton_proton_xx:
			return new hadronic_is_xx<model_t,rng_t,Minkowski_type>(false,false);
		    case initial_states::proton_antiproton_xx:
			return new hadronic_is_xx<model_t,rng_t,Minkowski_type>(false,true);
		    case initial_states::antiproton_proton_xx:
			return new hadronic_is_xx<model_t,rng_t,Minkowski_type>(true,false);
		    case initial_states::antiproton_antiproton_xx:
			return new hadronic_is_xx<model_t,rng_t,Minkowski_type>(true,true);
		    default:
                        log(log_level::error)<<CAMGEN_STREAMLOC<<"initial state type "<<is_type<<" is not a valid two-particle initial state"<<endlog;
			return NULL;
		}
	    }
    };
}

#endif /*CAMGEN_ISGEN_FAC_H_*/
