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
#include <Camgen/ps_decl.h>
#include <Camgen/MC_gen.h>
#include <Camgen/s_gen.h>

namespace Camgen
{
    /* Base class template for phase space channel branching algorithms. */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_branching: public MC_generator<typename model_t::value_type>
    {
	/* Friend classes: */

	friend class particle_channel<model_t,N_in,N_out,rng_t>;
	friend class ps_tree<model_t,N_in,N_out,rng_t>;
	friend class ps_factory<model_t,N_in,N_out,rng_t,typename model_t::spacetime_type>;

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
	    typedef s_generator<value_type,rng_t> s_generator_type;

	    static const std::size_t N_bits=N_in+N_out-1;
	    static const std::size_t beam_direction=model_t::beam_direction;
	    
	    typedef bit_string<N_bits> bit_string_type;

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Recursively sets the redundant flag: */
	    
	    void set_redundant()
	    {
		redundant=true;
	    }

	    /* Recursively unsets the redundant flag: */
	    
	    virtual void unset_redundant()
	    {
		redundant=false;
		for(size_type i=0;i<channels.size();++i)
		{
		    channels[i]->unset_redundant();
		}
	    }

	    /* Sets the generation flag to false: */

	    void reset_generation_flags()
	    {
		generation_flag=false;
		for(size_type i=0;i<channels.size();++i)
		{
		    channels[i]->reset_generation_flags();
		}
	    }

	    /* Passes the channel picking procedure to the produced momentum
	     * channels: */

	    bool choose_channel(std::vector<ps_branching<model_t,N_in,N_out,rng_t>*>& b)
	    {
		generation_flag=true;
		++callcount;
		for(size_type i=0;i<channels.size();++i)
		{
		    if(!channels[i]->choose_branching(b))
		    {
			return false;
		    }
		}
		return true;
	    }

	    /* Generates positive invariants: */

	    virtual bool generate_s()
	    {
		return true;
	    }

	    /* Generates negative invariants: */

	    virtual bool generate_t()
	    {
		return true;
	    }

	    /* Generates momenta: */

	    virtual bool generate_p()=0;
	    
	    /* Recursively generates momenta. */
	    
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
		for(size_type i=0;i<channels.size();++i)
		{
		    if(!channels[i]->generate())
		    {
			return false;
		    }
		}
		return true;
	    }

	    /* Computes the momentum branching weight. */

	    virtual bool evaluate_branching_weight()=0;

	    /* Computes the recursive channel branching weight. */

	    bool evaluate_weight()
	    {
		value_type channel_weights=(value_type)1;
		for(size_type i=0;i<channels.size();++i)
		{
		    if(channels[i]->weight()==(value_type)0)
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    channel_weights*=channels[i]->weight();
		}
		if(!generation_flag)
		{
		    if(!evaluate_branching_weight())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		}
		else if(this->weight()==(value_type)0)
		{
		    return false;
		}
		this->weight()*=channel_weights;
		return true;
	    }

	    /* Resets cross sections, adaptive grids and multichannels weights:
	     * */

	    virtual void reset()
	    {
		this->base_type::reset();
		callcount=0;
		generation_flag=false;
	    }

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

	    /* Number of outgoing phase space channels. */

	    size_type out_rank() const
	    {
		return channels.size();
	    }

	    bool t_type() const
	    {
		return channels[0]->spacelike();
	    }

	    /* Returns a const reference to the incoming momentum. */
	    
	    const momentum_type& p_in() const
	    {
		return incoming_channel->p();
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
	    
	    const value_type& m_in() const
	    {
		return incoming_channel->m();
	    }

	    /* Returns a reference to the i-th outgoing invariant mass. */

	    const value_type& m_out(size_type i) const
	    {
		return channels[i]->m();
	    }

	    /* Returns the type identifyer for printing. */

	    virtual std::string type() const=0;

	    /* Returns the redundancy flag. */

	    bool is_redundant() const
	    {
		return redundant;
	    }

	    /* Returns equivalence of branchings: */

	    virtual bool equiv(const ps_branching<model_t,N_in,N_out,rng_t>* other) const
	    {
		return false;
	    }

	    /* Returns the maximally allowed s-pair throws per event: */

	    size_type max_s_pairs() const
	    {
		return (spairs==NULL)?1:(*spairs);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Printing method. */

	    std::ostream& print(std::ostream& os) const
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
		os<<std::setw(15)<<std::left<<alpha;
		os<<std::setw(15)<<std::left<<callcount;
		os<<std::setw(15)<<std::left<<W/callcount;
		os<<std::setw(10)<<std::left<<generation_flag;
		return os;
	    }

	    /* Overridden loading method: */

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		safe_read(is,alpha);
		safe_read(is,W);
		is>>callcount;
		load_data(is);
		return is;
	    }

	    /* Virtual method for derived classes to load their subclass-data:
	     * */

	    virtual std::istream& load_data(std::istream& is)
	    {
		return is;
	    }

	    /* Overridden saving method: */

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<branching>"<<std::endl;
		os<<type()<<std::endl;
		save_channels(os);
		this->base_type::save(os);
		safe_write(os,alpha);
		os<<"\t";
		safe_write(os,W);
		os<<"\t"<<callcount<<std::endl;
		save_data(os);
		os<<"</branching>"<<std::endl;
		return os;
	    }

	    /* Virtual method for derived classes to save their subclass-data:
	     * */

	    virtual std::ostream& save_data(std::ostream& os) const
	    {
		return os;
	    }

	    virtual std::ostream& save_channels(std::ostream& os) const
	    {
		os<<incoming_channel->name;
		os<<std::endl;
		for(size_type i=0;i<channels.size();++i)
		{
		    os<<channels[i]->name<<"\t";
		}
		os<<std::endl;
		return os;
	    }

	protected:

	    /* Inupt phase space channel: */
	    
	    channel_type* incoming_channel;

	    /* Outgoing phase space channels: */

	    std::vector<channel_type*>channels;

	    /* Constructor. */
	    
	    ps_branching(channel_type* incoming_channel_,size_type n_out):incoming_channel(incoming_channel_),generation_flag(false),W(0),redundant(true),callcount(0),spairs(NULL)
	    {
		if(n_out<=0)
		{
		    log(log_level::error)<<CAMGEN_STREAMLOC<<"attempt to construct invalid branching of rank "<<n_out<<endlog;
		}
		else
		{
		    channels.resize(n_out,NULL);
		}
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

	private:

	    /* Private data members: */
	    /*-----------------------*/

	    /* Generation flag: */

	    bool generation_flag;

	    /* Multichannel weight: */

	    value_type alpha;
	    
	    /* Multichannel variance measure: */
	    
	    value_type W;

	    /* Redundancy flag: */

	    bool redundant;

	    /* Call counter utility: */

	    size_type callcount;

	    /* Maximal number of s-pair generations: */

	    const size_type* spairs;

	    /* Private default constructor: */

	    ps_branching():spairs(NULL){}
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t ps_branching<model_t,N_in,N_out,rng_t>::beam_direction;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t ps_branching<model_t,N_in,N_out,rng_t>::N_bits;
    
    /* Orders two branching instances by decreasing multichannel weight. */

    template<class spacetime_t,std::size_t N_in,std::size_t N_out,class rng_t>bool multichannel_ordering(const ps_branching<spacetime_t,N_in,N_out,rng_t>* b1,const ps_branching<spacetime_t,N_in,N_out,rng_t>* b2)
    {
	return b1->alpha>b2->alpha;
    }
}

#endif /*CAMGEN_PS_BRANCHING_H_*/

