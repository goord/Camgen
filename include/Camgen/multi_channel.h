//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file multi_channel.h
    \brief Muli-channel aggregation of Monte Carlo generators.
 */

#ifndef CAMGEN_MULTI_CHANNEL_H_
#define CAMGEN_MULTI_CHANNEL_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Multi-channel Monte Carlo generator, created from a vector of single-channel  *
 * generators                                                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <vector>
#include <algorithm>
#include <Camgen/MC_gen.h>
#include <Camgen/MC_int_base.h>

namespace Camgen
{
    /// Multichannel Monte Carlo generator class.

    template<class value_t,class rng_t>class multi_channel: public virtual MC_generator<value_t>, public MC_integrator_base<value_t>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef value_t value_type;
	    typedef std::size_t size_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef MC_generator<value_t> generator_type;

	    // Simple utility struct holding generator, multichannel weight and cumulant: 

	    struct channel_type
	    {
		generator_type* generator;
		value_type alpha;
		value_type W;
	    };

	    typedef typename std::vector<channel_type>::iterator channel_iterator;
	    typedef typename std::vector<channel_type>::const_iterator const_channel_iterator;

	    /* Public data members: */
	    /*----------------------*/

	    /// Adaption exponent parameter.

	    value_type adaptive_exponent;
	    
	    /// Multichannel threshold parameter. Channels with weights below
	    /// this value are discarded.
	    
	    value_type channel_threshold;

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Default onstructor.

	    multi_channel():adaptive_exponent(1),channel_threshold(0),update_flag(false),update_counter(0){}

	    /// Constructor with channels argument.

	    multi_channel(const std::vector<generator_type*>& generators):adaptive_exponent(1),channel_threshold(0),update_flag(false),update_counter(0)
	    {
		if(generators.size()==0)
		{
		    return;
		}
		value_type alpha=((value_type)1)/generators.size();
		for(int i=0;i<generators.size();++i)
		{
		    struct channel_type channel={generators[i],alpha,0};
		    channels.push_back(channel);
		}
		current_channel_iterator=generators.end();
	    }

	    /// Destructor.

	    virtual ~multi_channel(){}

	    /// Method that picks the generator channel according to the
	    /// relative weights.

	    generator_type* select_channel()
	    {
		if(channels.size()==0)
		{
		    current_channel_iterator=channels.end();
		    return NULL;
		}
		current_channel_iterator=channels.begin();
		if(channels.size()==1)
		{
		    return current_channel_iterator->generator;
		}
		value_type rho=rn_stream::throw_number();
		value_type r(0);
		while(r<=rho and current_channel_iterator!=channels.end())
		{
		    r+=(current_channel_iterator->alpha);
		    ++current_channel_iterator;
		}
		--current_channel_iterator;
		if(r>(value_type)1+numeric_configuration<value_type>::epsilon_abs)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"multichannel weights incorrectly normalized."<<r<<endlog;
		}
		return current_channel();
	    }

	    /// Resets the chosen sub-channel.

	    void reset_current_channel()
	    {
		current_channel_iterator=channels.end();
	    }

	    /// Returns the channel chosen by select_channel.

	    generator_type* current_channel()
	    {
		return (current_channel_iterator==channels.end())?NULL:current_channel_iterator->generator;
	    }

	    /// Generation method implementation.

	    bool generate()
	    {
		select_channel();
		if(channels.size()==0)
		{
		    this->weight()=(value_type)1;
		    return true;
		}
		generator_type* channel=current_channel();
		if(channel==NULL)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!channel->generate())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		return evaluate_weight();
	    }

	    /// Weight evaluation method.
	    
	    bool evaluate_weight()
	    {
		return evaluate_weight(false);
	    }

	    /// Weight evaluation method with recursion flag.

	    bool evaluate_weight(bool recurse)
	    {
		if(channels.size()==0)
		{
		    this->weight()==(value_type)1;
		    return true;
		}
		value_type g(0);
		bool zero_weight=false;
		for(channel_iterator it=channels.begin();it!=channels.end();++it)
		{
		    if(recurse)
		    {
			if(!it->generator->evaluate_weight())
			{
			    this->weight()=(value_type)0;
			    return false;
			}
		    }
		    if(it->generator->weight()==(value_type)0 and it->alpha>(value_type)0)
		    {
			zero_weight=true;
		    }
		    else if(it->generator->weight()!=(value_type)0)
		    {
			g+=(it->alpha/it->generator->weight());
		    }
		}
		this->weight()=zero_weight?(value_type)0:((value_type)1)/g;
		return true;
	    }

	    /// Updates internal multichannel weights.

	    void update()
	    {
		if(channels.size()<=1)
		{
		    return;
		}
		value_type w3=std::pow(this->integrand(),(int)2)*(this->weight());
		if(w3!=(value_type)0)
		{
		    update_flag=true;
		}
		for(channel_iterator it=channels.begin();it!=channels.end();++it)
		{
		    it->W+=(w3/it->generator->weight());
		}
		++update_counter;
	    }

	    /// Adapts multichannel weights.

	    void adapt()
	    {
		if(channels.size()<2 || !update_flag)
		{
		    return;
		}
		value_type beta=(value_type)0.5*adaptive_exponent;
		for(channel_iterator it=channels.begin();it!=channels.end();++it)
		{
		    it->W/=update_counter;
		    it->alpha*=std::pow(it->W,beta);
		}
		normalise_channel_weights();
		update_flag=false;
		update_counter=0;
	    }

	    /// Adds a sub-channel, resets all weights if necessary.

	    void add_generator(generator_type* generator)
	    {
		struct channel_type channel={generator,0,0};
		channels.push_back(channel);
		reset_channel_weights();
	    }

	    /// Removes a sub-channel, resets all weights if necessary.

	    void remove_generator(generator_type* generator)
	    {
		match_channels predicate(generator);
		channel_iterator it=std::remove_if(channels.begin(),channels.end(),predicate);
		channels.erase(it,channels.end());
		reset_channel_weights();
	    }

	    /// Replaces a sub-channel.

	    void replace_generator(generator_type* first,generator_type* second)
	    {
		match_channels predicate(first);
		channel_iterator it=std::find_if(channels.begin(),channels.end(),predicate);
		if(it==channels.end())
		{
		    return;
		}
		it->generator=second;

	    }

	    /// Queries whether the given generator is in the multichannel

	    bool contains_generator(generator_type* generator) const
	    {
		match_channels predicate(generator);
		return std::find_if(channels.begin(),channels.end(),predicate)!=channels.end();
	    }


	    /* Public readout functions: */
	    /*---------------------------*/

	    /// Returns the number of channels, whee the argument controls
	    /// whether to count the zero-weight channels.

	    size_type channel_count() const
	    {
		return channels.size();
	    }

	    /// Returns the i-th channel (no bound-checking).

	    const generator_type* generator(size_type i) const
	    {
		return channels[i].generator;
	    }

	    /// Returns the i-th multichannel weight (no bound-checking).

	    value_type alpha(size_type i) const
	    {
		return channels[i].alpha;
	    }

	private:

	    /* List sub-channels: */

	    std::vector<channel_type> channels;

	    /* Iterator pointing to the last chosen channel: */

	    channel_iterator current_channel_iterator;
	    
	    /* Update flag: */
	    
	    bool update_flag;

	    /* Update counter: */

	    size_type update_counter;

	    /* Normalises the channel weights: */

	    void normalise_channel_weights()
	    {
		if(channels.size()==0)
		{
		    return;
		}
		value_type norm(0);
		for(channel_iterator it=channels.begin();it!=channels.end();++it)
		{
		    norm+=(it->alpha);
		}
		value_type alpha_min=channel_threshold*norm/channels.size();
		value_type normbis(0);
		for(channel_iterator it=channels.begin();it!=channels.end();++it)
		{
		    if(it->alpha<alpha_min)
		    {
			it->alpha=(value_type)0;
		    }
		    else
		    {
			normbis+=it->alpha;
		    }
		}
		for(channel_iterator it=channels.begin();it!=channels.end();++it)
		{
		    it->alpha/=normbis;
		}
	    }

	    /* Resets all multichannel weights: */

	    void reset_channel_weights()
	    {
		value_type alpha=((value_type)1)/channels.size();
		for(channel_iterator it=channels.begin();it!=channels.end();++it)
		{
		    it->alpha=alpha;
		}
	    }

	    /* Channel matching predicate class: */

	    class match_channels
	    {
		public:

		    const generator_type* generator;

		    match_channels(const generator_type* generator_):generator(generator_){}

		    bool operator()(const channel_type& channel)
		    {
			return channel.generator==generator;
		    }
	    };
    };
}

#endif /*CAMGEN_MULTI_CHANNEL_H_*/

