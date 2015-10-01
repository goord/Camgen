//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PS_BRANCHING_FAC_H_
#define CAMGEN_PS_BRANCHING_FAC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Factory class for creating momentum channels, particle channels, propagator *
 * samplers and phase space tree branching objects that populate the recursive *
 * phase space tree.                                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/ps_channel.h>
#include <Camgen/p_channel.h>
#include <Camgen/s_branch.h>
#include <Camgen/t_branch.h>
#include <Camgen/sum_branch.h>

namespace Camgen
{

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t>class ps_channel_factory
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_type::value_type value_type;
	    typedef momentum_channel<model_t,N_in,N_out,rng_t> momentum_channel_type;
	    typedef particle_channel<model_t,N_in,N_out,rng_t> particle_channel_type;
	    typedef typename particle_channel_type::s_generator_type s_generator_type;
	    typedef ps_branching<model_t,N_in,N_out,rng_t> branching_type;
	    typedef particle<model_t> particle_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;
	    typedef typename momentum_channel_type::momentum_type momentum_type;
	    typedef typename momentum_channel_type::size_type size_type;

	    static s_generator_type* create_s_generator(const particle_channel_type* channel,const value_type* Ecmhat=NULL)
	    {
		const particle_type* phi=channel->particle_type;
		const value_type* mass=(phi==NULL)?NULL:(phi->get_mass_address());
		const value_type* width=(phi==NULL or channel->spacelike())?NULL:(phi->get_width_address());
		const value_type* exponent=&(channel->s_sampling_exponent);

		bool adaptive=(channel->timelike())?adaptive_s_sampling():adaptive_t_sampling();
		bool dirac_delta=false;	
		if(channel->on_shell()) // final state particles: Dirac delta
		{
		    width=NULL;
		    exponent=NULL;
		    adaptive=false;
		}
		bool total_s=channel->timelike() and channel->bitstring().count()==N_out;
		if(total_s and Ecmhat!=NULL) // total s-channel with fixed invariant mass: Dirac delta
		{
		    width=NULL;
		    mass=Ecmhat;
		    exponent=NULL;
		    adaptive=false;
		}
		return value_generator_factory<value_type,rng_t>::create_instance(mass,width,exponent,adaptive);
	    }

	    static particle_channel_type* find_or_add_particle_channel(momentum_channel_type* ps_channel,const particle_type* phi,const value_type* Ecmhat=NULL)
	    {
		if(phi!=NULL and phi->get_pdg_id()==0 and !phi->is_auxiliary()) //phi is pure gauge->do nothing.
		{
		    return NULL;
		}

		particle_channel_type* found_channel=ps_channel->find_particle_channel(phi);
		if(found_channel!=NULL)
		{
		    return found_channel;
		}

		particle_channel_type* new_channel=new particle_channel_type(ps_channel,phi);
		new_channel->s_sampling_exponent=sampling_exponent(new_channel);
		new_channel->set_s_generator(create_s_generator(new_channel,Ecmhat));
		ps_channel->add_particle_channel(new_channel);
		return new_channel;
	    }

	    static value_type sampling_exponent(particle_channel_type* channel)
	    {
		if(channel->particle_type==NULL)
		{
		    if(channel->timelike() and channel->bitstring().count()==N_out)
		    {
			return shat_exponent();
		    }
		    else
		    {
			return auxiliary_exponent();
		    }
		}
		else
		{
		    if(channel->particle_type->is_auxiliary())
		    {
			return auxiliary_exponent();
		    }
		    if(channel->timelike())
		    {
			return timelike_exponent(channel->particle_type->get_name());
		    }
		    else
		    {
			return spacelike_exponent(channel->particle_type->get_name());
		    }
		}
	    }

	protected:

	    /* Static utility function for read-in from file: */

	    static particle_channel_type* find_channel(const std::string& name, const std::list<momentum_channel_type*>& ps_channels)
	    {
		for(typename std::list<momentum_channel_type*>::const_iterator it=ps_channels.begin();it!=ps_channels.end();++it)
		{
		    for(typename momentum_channel_type::particle_channel_list::const_iterator it2=(*it)->begin_particle_channels();it2!=(*it)->end_particle_channels();++it2)
		    {
			if(name==(*it2)->name)
			{
			    return *it2;
			}
		    }
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle channel "<<name<<" not found in list--returning NULL"<<endlog;
		return NULL;
	    }
    };

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t>class ps_branching_factory: public ps_channel_factory<model_t,N_in,N_out,rng_t,spacetime_t>{};

    template<class model_t,std::size_t N_out,class rng_t>class ps_branching_factory<model_t,1,N_out,rng_t,Minkowski_type>: public ps_channel_factory<model_t,1,N_out,rng_t,Minkowski_type>
    {
	typedef ps_channel_factory<model_t,1,N_out,rng_t,Minkowski_type> base_type;
	
	public:

	    /* Type definitions: */

	    typedef typename base_type::model_type model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_channel_type momentum_channel_type;
	    typedef typename base_type::particle_channel_type particle_channel_type;
	    typedef typename base_type::branching_type branching_type;
	    typedef typename base_type::particle_type particle_type;
	    typedef typename base_type::bit_string_type bit_string_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::size_type size_type;

	    typedef typename Minkowski_type::template implementation<typename model_t::value_type,model_t::dimension> spacetime_type;

	    /* Creates a forward s-type branching of the channel in. */

	    static branching_type* create_s_branching(particle_channel_type* in,particle_channel_type* out1,particle_channel_type* out2)
	    {
		bit_string_type bm=(out1->bitstring() & out2->bitstring());
		bit_string_type bp=(out1->bitstring() | out2->bitstring());
		if(bm.any() or bp!=in->bitstring())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-branching into invalid channels "<<out1->bitstring()<<","<<out2->bitstring()<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		return new s_branching<model_t,1,N_out,rng_t,spacetime_type>(in,out1,out2);
	    }
    };

    /* Phase space factory implementation for 2->N scattering processes. */

    template<class model_t,std::size_t N_out,class rng_t>class ps_branching_factory<model_t,2,N_out,rng_t,Minkowski_type>: public ps_channel_factory<model_t,2,N_out,rng_t,Minkowski_type>
    {
	typedef ps_channel_factory<model_t,2,N_out,rng_t,Minkowski_type> base_type;
	
	public:

	    /* Type definitions: */

	    typedef typename base_type::model_type model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_channel_type momentum_channel_type;
	    typedef typename base_type::particle_channel_type particle_channel_type;
	    typedef typename base_type::branching_type branching_type;
	    typedef typename base_type::particle_type particle_type;
	    typedef typename base_type::bit_string_type bit_string_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::size_type size_type;
	    
	    typedef typename Minkowski_type::template implementation<typename model_t::value_type,model_t::dimension> spacetime_type;

	    /* Creates a forward s-type branching of the channel in: */

	    static branching_type* create_s_branching(particle_channel_type* in,particle_channel_type* out1,particle_channel_type* out2)
	    {
		if(in->bitstring()[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"spacelike incoming channel "<<in->bitstring()<<" allows no s-branching--returning NULL"<<endlog;
		    return NULL;
		}
		bit_string_type bm=(out1->bitstring() & out2->bitstring());
		bit_string_type bp=(out1->bitstring() | out2->bitstring());
		if(bm.any() or bp!=in->bitstring())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-branching into invalid channels "<<out1->bitstring()<<","<<out2->bitstring()<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		return new s_branching<model_t,2,N_out,rng_t,spacetime_type>(in,out1,out2);
	    }

	    /* Creates a forward momentum sum branching: */

	    static branching_type* create_sum_branching(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type ss=out->bitstring();
		sa.flip();
		
		if(!(sb.count()==1 and sb[0]))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming channel "<<sb<<" is not the second beam--returning NULL"<<endlog;
		    return NULL;
		}
		return new sum_branching<model_t,2,N_out,rng_t>(in1,in2,out);
	    }

	    /* Creates a forward t-type branching of the channel in1: */
	    
	    static branching_type* create_t_branching(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out1,particle_channel_type* out2,particle_channel_type* out3,particle_channel_type* out4)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type s1=out1->bitstring();
		bit_string_type st=out2->bitstring();
		bit_string_type s2=out3->bitstring();
		bit_string_type stot=out4->bitstring();

		if(!sa[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"incoming channel "<<sa<<" cannot be part of a t-channel--returning NULL"<<endlog;
		    return NULL;
		}
		if(!(sb[0] and sb.count()==1))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming momentum "<<sb<<" is not the second beam momentum--returning NULL"<<endlog;
		    return NULL;
		}
		if(s1[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first outgoing channel "<<s1<<" is not timelike--returning NULL"<<endlog;
		    return NULL;
		}
		if((st & s1).any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		if((st | s1)!=sa)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" and incoming channel "<<sa<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		if((s1 | s2)!=stot)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"total outgoing channel "<<stot<<" is not the combination of the outgoing channels "<<s1<<" and "<<s2<<"--returning NULL"<<endlog;
		    return NULL;
		}
		bit_string_type bt=st;
		bt.reset(0);
		if(bt!=s2)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and second s-channel "<<s2<<" ignored--returning NULL"<<endlog;
		   return NULL;
		}
		return new t_branching<model_t,N_out,rng_t,spacetime_type>(in1,in2,out1,out2,out3,out4,false);
	    }
    };
}

#endif /*CAMGEN_PS_BRANCHING_FAC_H_*/

