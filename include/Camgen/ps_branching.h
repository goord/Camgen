//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PS_BRANCHING_H_
#define CAMGEN_PS_BRANCHING_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the phase space branching base class, which *
 * serves as base for all objects splitting a channel into one or more other *
 * phase space channels.                                                     *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/vector.h>
#include <Camgen/bit_string.h>
#include <Camgen/type_holders.h>
#include <Camgen/rn_strm.h>
#include <Camgen/MC_gen.h>
#include <Camgen/val_gen.h>

namespace Camgen
{
    /* Base class template for phase space channel branching algorithms. */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_branching: public MC_generator<typename model_t::value_type>, 
    											      public MC_integrator_base<typename model_t::value_type>
    {
	/* Private type definitions: */

	typedef MC_generator<typename model_t::value_type> base_type;

	public:

	    /* Public type definitions: */

	    typedef model_t model_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef typename model_t::value_type value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef std::size_t size_type;
	    typedef momentum_channel<model_t,N_in,N_out,rng_t> ps_channel_type;
	    typedef particle_channel<model_t,N_in,N_out,rng_t> channel_type;
	    typedef value_generator<value_type,rng_t> s_generator_type;
	    typedef typename std::vector<channel_type*>::iterator channel_iterator;
	    typedef typename std::vector<channel_type*>::const_iterator const_channel_iterator;

	    static const std::size_t N_bits=N_in+N_out-1;
	    static const std::size_t beam_direction=model_t::beam_direction;
	    
	    typedef bit_string<N_bits> bit_string_type;

	    /* Inupt phase space channel: */
	    
	    channel_type* const incoming_channel;

	    /* Controls the direction of s-sampling: */

	    bool backward_s_sampling;

	    /* Denotes whether this branching can sample the incoming s-hat: */

	    bool shat_sampling;

	    /* Utility flag for generation: */

	    bool generated;

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Virtual destructor: */

	    virtual ~ps_branching(){}

	    /* Samples incoming positive invariants (no recursion): */
	    
	    virtual bool generate_s()=0;

	    /* Generates negative invariants (no recursion): */

	    virtual bool generate_t()=0;

	    /* Generates momenta (no recursion): */

	    virtual bool generate_p()=0;

	    /* Recursively generates momenta */
	    
	    bool generate()
	    {
		if(!generate_s())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!generate_t())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!generate_p())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		return true;
	    }

	    /* Computes the momentum branching weight. */

	    virtual bool evaluate_branching_weight()=0;

	    bool evaluate_weight()
	    {
		return evaluate_weight(true);
	    }

	    /* Computes the recursive channel branching weight. */
	    
	    bool evaluate_weight(bool eval_bw)
	    {
		for(size_type i=0;i<channels.size();++i)
		{
		    if(channels[i]->get_status()!=ps_channel_type::p_set)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"evaluating weight of incomplete event..."<<endlog;
		    }
		}
		if(eval_bw and !evaluate_branching_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->weight()=branching_weight;
		for(size_type i=0;i<channels.size();++i)
		{
		    this->weight()*=channels[i]->weight();
		}
		return true;
	    }

	    /* Trivial update method: */

	    virtual void update(){}

	    /* Adapts grids: */

	    virtual void adapt_grids(){}

	    /* Adapts multichannel weights: */

	    virtual void adapt_channels(){}

	    /* Full adaptation method: */

	    void adapt()
	    {
		adapt_grids();
		adapt_channels();
	    }
	    
	    /* Public readout functions: */
	    /*---------------------------*/

	    /* Returns the starting iterator of the outgoing channels: */
	    
	    channel_iterator begin_channels_out()
	    {
		return channels.begin();
	    }

	    /* Returns the const starting iterator of the outgoing channels: */
	    
	    const_channel_iterator begin_channels_out() const
	    {
		return channels.begin();
	    }

	    /* Returns the end iterator of the outgoing channels: */

	    channel_iterator end_channels_out()
	    {
		return channels.end();
	    }

	    /* Returns the const end iterator of the outgoing channels: */

	    const_channel_iterator end_channels_out() const
	    {
		return channels.end();
	    }

	    /* Returns the i-th outgoing particle channel: */

	    channel_type* channel(size_type i)
	    {
		return channels[i];
	    }

	    /* Returns the i-th const outgoing particle channel: */

	    const channel_type* channel(size_type i) const
	    {
		return channels[i];
	    }

	    /* Number of outgoing phase space channels. */

	    size_type out_rank() const
	    {
		return channels.size();
	    }

	    virtual bool t_type() const
	    {
		return false;
	    }

	    virtual bool s_type() const
	    {
		return false;
	    }

	    /* Returns a const reference to the incoming momentum. */
	    
	    const momentum_type& p_in() const
	    {
		return incoming_channel->p();
	    }

	    virtual momentum_type evaluate_p_in() const
	    {
		return p_out();
	    }

	    /* Returns the total outgoing momentum. */

	    momentum_type p_out() const
	    {
		momentum_type p(channels[0]->p());
		for(size_type i=1;i<channels.size();++i)
		{
		    p+=channels[i]->p();
		}
		return p;
	    }

	    /* Returns a reference to the i-th outgoing momentum. */

	    momentum_type& p_out(size_type i)
	    {
		return channels[i]->p();
	    }

	    /* Returns a const reference to the i-th outgoing momentum. */

	    const momentum_type& p_out(size_type i) const
	    {
		return channels[i]->p();
	    }
	    
	    /* Returns a const reference to the incoming invariant
	    mass-squared. */
	    
	    const value_type& s_in() const
	    {
		return incoming_channel->s();
	    }

	    /* Returns a reference to the i-th outgoing invariant mass-squared.
	     * */

	    value_type& s_out(size_type i)
	    {
		return channels[i]->s();
	    }

	    /* Returns a const reference to the i-th outgoing invariant
	     * mass-squared. */

	    const value_type& s_out(size_type i) const
	    {
		return channels[i]->s();
	    }
	    
	    /* Returns a const reference to the incoming invariant mass. */
	    
	    value_type m_in() const
	    {
		return incoming_channel->m();
	    }

	    /* Returns a reference to the i-th outgoing invariant mass. */

	    value_type m_out(size_type i) const
	    {
		return channels[i]->m();
	    }

	    /* Returns the type identifyer for printing. */

	    virtual std::string type() const=0;

	    /* Returns equivalence of branchings: */

	    virtual bool equiv(const ps_branching<model_t,N_in,N_out,rng_t>* other) const
	    {
		return false;
	    }

	    value_type branch_weight() const
	    {
		return branching_weight;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Printing method. */

	    virtual std::ostream& print(std::ostream& os) const
	    {
		std::stringstream ss1;
		size_type n=out_rank()-1;
		for(size_type i=0;i<n;++i)
		{
		    ss1<<channels[i]->name<<",";
		}
		ss1<<channels.back()->name;
		os<<std::setw(2*(N_in+N_out)+20)<<std::left<<ss1.str();
		os<<std::setw(15)<<std::left<<type();
		os<<std::setw(15)<<std::left<<this->weight();
		return os;
	    }

	protected:

	    /* Constructor 1->1 branching. */
	    
	    ps_branching(channel_type* incoming_channel_,channel_type* channel):incoming_channel(incoming_channel_),backward_s_sampling(false),shat_sampling(false),generated(false)
	    {
		channels.push_back(channel);
	    }

	    /* Constructor 1->2 branching. */
	    
	    ps_branching(channel_type* incoming_channel_,channel_type* channel1,channel_type* channel2):incoming_channel(incoming_channel_),backward_s_sampling(false),shat_sampling(false),generated(false)
	    {
		channels.push_back(channel1);
		channels.push_back(channel2);
	    }

	    /* Returns the i-th invariant mass generator: */

	    value_generator<value_type,rng_t>* s_generator(size_type i)
	    {
		return channels[i]->s_gen;
	    }

	    static value_generator<value_type,rng_t>* s_generator(channel_type* channel)
	    {
		return channel->s_gen;
	    }

	    /* Returns whether the argument branching has the same channels: */

	    bool same_channels(const ps_branching<model_t,N_in,N_out,rng_t>* other,bool sorted=true) const
	    {
		if(incoming_channel!=other->incoming_channel)
		{
		    return false;
		}
		if(sorted)
		{
		    return (channels==other->channels);
		}
		else
		{
		    std::set<channel_type*>s1,s2;
		    s1.insert(channels.begin(),channels.end());
		    s2.insert(other->channels.begin(),other->channels.end());
		    return (s1==s2);
		}
	    }

	    /* Returns whether the incoming channel has the same incoming
	     * channel: */

	    bool same_incoming_channel(const ps_branching<model_t,N_in,N_out,rng_t>* other) const
	    {
		return (incoming_channel==other->incoming_channel);
	    }

	    /* Branching weight: */

	    value_type branching_weight;

	private:

	    /* Outgoing phase space channels: */

	    std::vector<channel_type*>channels;
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t ps_branching<model_t,N_in,N_out,rng_t>::beam_direction;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t ps_branching<model_t,N_in,N_out,rng_t>::N_bits;
}

#endif /*CAMGEN_PS_BRANCHING_H_*/

