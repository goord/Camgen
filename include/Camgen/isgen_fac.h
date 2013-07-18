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
	    typedef initial_states::type generator_tag;

	    /// Factory method creating an initial state generator instance from
	    /// the static configuration data.

	    static generator_type* create_instance()
	    {
		bool q=true;
		return create_instance(initial_state_type(),q);
	    }

	    /// Factory method from initial state tag. Returns NULL.

	    static generator_type* create_instance(generator_tag tag,bool& backward_shat)
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"initial state type not defined for current spacetime type--returning NULL"<<endlog;
		return NULL;
	    }

	    /// Factory method from input file stream. Returns NULL.

	    static generator_type* create_instance(std::istream& is)
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<isgen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read"<<endlog;
		    return NULL;
		}
		generator_type* result;
		std::string is_type;
		is>>is_type;
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"initial state type "<<is_type<<" not recognised for current spacetime type--returning NULL"<<endlog;
		std::string endflag;
		do
		{
		    std::getline(is,endflag);
		}
		while(endflag!="</isgen>" and !is.eof());
		return NULL;
	    }
    };

    /// Initial state factory specialisation for decays in Minkowski spacetimes.

    template<class model_t,class rng_t>class initial_state_factory<model_t,1,rng_t,Minkowski_type>
    {
	public:

	    /* Type definitions: */

	    typedef initial_state<model_t,1> generator_type;
	    typedef initial_states::type generator_tag;

	    /// Factory method creating an initial state generator instance from
	    /// the static configuration data.

	    static generator_type* create_instance()
	    {
		return create_instance(initial_state_type(),true);
	    }

	    /// Factory method from initial state tag. Returns NULL unless the
	    /// tag is partonic.

	    static generator_type* create_instance(generator_tag tag,bool backward_shat=false)
	    {
		if(tag==initial_states::partonic and !backward_shat)
		{
		    return new partonic_is<model_t,1,Minkowski_type>;
		}
		return NULL;
	    }

	    /// Factory method from input stream.
	    
	    static generator_type* create_instance(std::istream& is)
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<isgen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read"<<endlog;
		    return NULL;
		}
		generator_type* result;
		std::string is_type;
		is>>is_type;
		if(is_type=="partonic")
		{
		    result=new partonic_is<model_t,1>;
		    result->load(is);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"initial state type "<<is_type<<" not recognised for current spacetime type--returning NULL"<<endlog;
		    result=NULL;
		}
		std::string endflag;
		do
		{
		    std::getline(is,endflag);
		}
		while(endflag!="</isgen>" and !is.eof());
		return result;
	    }
    };

    /// Initial state factory specialisation for 2-particle initial states in
    /// Minkowski spacetimes.

    template<class model_t,class rng_t>class initial_state_factory<model_t,2,rng_t,Minkowski_type>
    {
	public:

	    typedef initial_state<model_t,2> generator_type;
	    typedef initial_states::type generator_tag;

	    /// Factory method creating an initial state generator instance from
	    /// the static configuration data.

	    static generator_type* create_instance()
	    {
		return create_instance(initial_state_type(),backward_shat_sampling());
	    }

	    /// Factory method with initial state tag argument.

	    static generator_type* create_instance(generator_tag tag,bool backward_shat)
	    {
		switch(tag)
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
			return NULL;
		}
	    }

	    /// Factory method with input stream argument.
	    
	    static generator_type* create_instance(std::istream& is)
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<isgen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read"<<endlog;
		    return NULL;
		}
		generator_type* result;
		std::string is_type;
		is>>is_type;
		if(is_type=="partonic")
		{
		    result=new partonic_is<model_t,2>;
		    result->load(is);
		}
		else if(is_type=="pp_xx")
		{
		    result=new hadronic_is_xx<model_t,rng_t>(false,false);
		    result->load(is);
		}
		else if(is_type=="ppbar_xx")
		{
		    result=new hadronic_is_xx<model_t,rng_t>(false,true);
		    result->load(is);
		}
		else if(is_type=="pbarp_xx")
		{
		    result=new hadronic_is_xx<model_t,rng_t>(true,false);
		    result->load(is);
		}
		else if(is_type=="pbarpbar_xx")
		{
		    result=new hadronic_is_xx<model_t,rng_t>(true,true);
		    result->load(is);
		}
		else if(is_type=="pp_sy")
		{
		    result=new hadronic_is_sy<model_t,rng_t>(false,false);
		    result->load(is);
		}
		else if(is_type=="ppbar_sy")
		{
		    result=new hadronic_is_sy<model_t,rng_t>(false,true);
		    result->load(is);
		}
		else if(is_type=="pbarp_sy")
		{
		    result=new hadronic_is_sy<model_t,rng_t>(true,false);
		    result->load(is);
		}
		else if(is_type=="pbarpbar_sy")
		{
		    result=new hadronic_is_sy<model_t,rng_t>(true,true);
		    result->load(is);
		}
		else if(is_type=="pp_y")
		{
		    result=new hadronic_is_y<model_t,rng_t>(false,false);
		    result->load(is);
		}
		else if(is_type=="ppbar_y")
		{
		    result=new hadronic_is_y<model_t,rng_t>(false,true);
		    result->load(is);
		}
		else if(is_type=="pbarp_y")
		{
		    result=new hadronic_is_y<model_t,rng_t>(true,false);
		    result->load(is);
		}
		else if(is_type=="pbarpbar_y")
		{
		    result=new hadronic_is_y<model_t,rng_t>(true,true);
		    result->load(is);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"initial state type "<<is_type<<" not recognised for current spacetime type--returning NULL"<<endlog;
		    result=NULL;
		}
		std::string endflag;
		do
		{
		    std::getline(is,endflag);
		}
		while(endflag!="</isgen>" and !is.eof());
		return result;
	    }
    };
}

#endif /*CAMGEN_ISGEN_FAC_H_*/
