//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_P_CHANNEL_H_
#define CAMGEN_P_CHANNEL_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Particle channel class definition. The particle channel is a sub-channel of a *
 * momentum channel, sampling the invariant mass acoording to the type of        *
 * particle propagating through the channel. If the particle type is not defined *
 * (multi-particle invariant masses), the sampling is performed by a power-law.  *
 * Adaptive grids may be overlayed to optimise the propagator sampling.          *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <list>
#include <limits>
#include <Camgen/bit_string.h>
#include <Camgen/particle.h>
#include <Camgen/flav_comp.h>
#include <Camgen/multi_channel.h>
#include <Camgen/uni_val_gen.h>
#include <Camgen/Dirac_delta.h>

namespace Camgen
{
    /* Phase space particle channel class: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class particle_channel: public MC_generator<typename model_t::value_type>, public MC_integrator_base<typename model_t::value_type>
    {
	friend class ps_branching<model_t,N_in,N_out,rng_t>;
	typedef MC_generator<typename model_t::value_type> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef typename model_t::value_type value_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef std::size_t size_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef momentum_channel<model_t,N_in,N_out,rng_t> momentum_channel_type;
	    typedef ps_branching<model_t,N_in,N_out,rng_t> branching_type;
	    typedef value_generator<value_type,rng_t> s_generator_type;
	    typedef typename momentum_channel_type::status_type status_type;

	    static const std::size_t N_bits=N_in+N_out-1;
	    typedef bit_string<N_bits> bit_string_type;

	    /* Public data members: */
	    /*----------------------*/

	    /* Particle type: */

	    const particle<model_t>* const particle_type;

	    /* Channel name: */

	    const std::string name;

	    /* Invariant mass sampling exponent: */

	    value_type s_sampling_exponent;

	    /* Constructors: */
	    /* ------------- */

	    particle_channel(momentum_channel_type* ps_channel_,const particle<model_t>* particle_type_=NULL):particle_type(particle_type_),name(make_channel_name(particle_type_,ps_channel_->bitstring)),s_sampling_exponent(0),ps_channel(ps_channel_),s_gen(new uniform_value_generator<value_type,rng_t>())
	    {
		s_gen->set_value(&(ps_channel->s()));
	    }

	    /* Destructor: */
	    /*-------------*/

	    virtual ~particle_channel()
	    {
		delete s_gen;
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Adds a branching instance to the multichannel: */

	    bool insert_branching(branching_type* br)
	    {
		if(ps_channel->bitstring.count()==1)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to attach phase space branching to external channel ignored--returning NULL"<<endlog;
		    return false;
		}
		std::vector<const base_type*>channel_generators(branching_multichannel.channel_count());
		for(size_type i=0;i<channel_generators.size();++i)
		{
		    channel_generators[i]=branching_multichannel.generator(i);
		}
		same_branching pred(br);
		typename std::vector<const base_type*>::iterator it=std::find_if(channel_generators.begin(),channel_generators.end(),pred);
		if(it!=channel_generators.end())
		{
		    return false;
		}
		branching_multichannel.add_generator(br);
		return true;
	    }

	    /* Removes a branching instance to the multichannel: */

	    bool remove_branching(branching_type* br)
	    {
		if(ps_channel->bitstring.count()==1)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to remove phase space branching to external channel ignored"<<endlog;
		    return false;
		}
		if(!branching_multichannel.contains_generator(br))
		{
		    return false;
		}
		branching_multichannel.remove_generator(br);
		return true;
	    }

	    /* Replaces the first branching by the second: */

	    bool replace_branching(branching_type* first,branching_type* second)
	    {
		if(!branching_multichannel.contains_generator(first))
		{
		    return false;
		}
		branching_multichannel.replace_generator(first,second);
		return true;
	    }

	    /* Sets the invariant mass generator object: */

	    bool set_s_generator(s_generator_type* s_gen_)
	    {
		if(s_gen_==NULL)
		{
		    return false;
		}
		delete s_gen;
		s_gen=s_gen_;
		s_gen->set_value(&(ps_channel->s()));
		
		Dirac_delta<value_type,rng_t>* dirac_delta=dynamic_cast<Dirac_delta<value_type,rng_t>*>(s_gen);

		if(dirac_delta!=NULL)
		{
		    const value_type* mass=dirac_delta->m;
		    set_m_min(mass==NULL?0:(*mass));
		    set_m_max(mass==NULL?0:(*mass));
		}
	    }

	    /* Sets the minimal invariant mass-squared. */
	    
	    bool set_s_min(const value_type& smin)
	    {
		return s_gen->set_lower_bound(smin);
	    }

	    /* Sets the maximal invariant mass-squared. */
	    
	    bool set_s_max(const value_type& smax)
	    {
		return s_gen->set_upper_bound(smax);
	    }

	    /* Sets the minimal signed invariant mass. */
	    
	    bool set_m_min(const value_type& mmin)
	    {
		return set_s_min(sgn_sq(mmin));
	    }

	    /* Sets the maximal signed invariant mass. */
	    
	    bool set_m_max(const value_type& mmax)
	    {
		return set_s_max(sgn_sq(mmax));
	    }

	    /* Sets the minimal minimal invariant mass-squared. */
	    
	    bool set_s_min_min(const value_type& sminmin)
	    {
		adaptive_value_generator<value_type,rng_t>* s_gen_grid=dynamic_cast<adaptive_value_generator<value_type,rng_t>*>(s_gen);
		if(s_gen_grid!=NULL)
		{
		    return s_gen_grid->set_mapping_lower_bound(sminmin);
		}
		else
		{
		    return s_gen->set_lower_bound(sminmin);
		}
	    }

	    /* Sets the maximal maximal invariant mass-squared. */
	    
	    bool set_s_max_max(const value_type& smaxmax)
	    {
		adaptive_value_generator<value_type,rng_t>* s_gen_grid=dynamic_cast<adaptive_value_generator<value_type,rng_t>*>(s_gen);
		if(s_gen_grid!=NULL)
		{
		    return s_gen_grid->set_mapping_upper_bound(smaxmax);
		}
		else
		{
		    return s_gen->set_upper_bound(smaxmax);
		}
	    }

	    /* Sets the invariant mass-squared range. */
	    
	    bool set_s_range(const value_type& smin,const value_type& smax)
	    {
		return s_gen->set_bounds(std::max(s_min_min(),smin),std::min(s_max_max(),smax));
	    }

	    /* Sets the invariant mass range. */
	    
	    bool set_m_range(const value_type& mmin,const value_type& mmax)
	    {
		return set_s_range(sgn_sq(mmin),sgn_sq(mmax));
	    }

	    /* Selects a branching: */

	    branching_type* select_branching()
	    {
		if(branching_multichannel.channel_count()==0)
		{
		    return NULL;
		}
		return static_cast<branching_type*>(branching_multichannel.select_channel());
	    }

	    /* Recursive momentum generation method. */

	    bool generate()
	    {
		if(branching_multichannel.channel_count()==0)
		{
		    this->weight()=(value_type)1;
		    return true;
		}
		base_type* selected_branching=select_branching();
		if(selected_branching==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Failed to select ps branching in particle channel "<<bitstring()<<endlog;
		    this->weight()=(value_type)0;
		    return false;
		}
		return selected_branching->generate();
	    }

	    /* Recursive weight evaluation method. */

	    bool evaluate_weight()
	    {
		bool q=branching_multichannel.evaluate_weight(false);
		this->weight()=branching_multichannel.weight();
		return q;
	    }

	    /* Generates the incoming invariant mass. */

	    bool generate_s()
	    {
		return s_gen->generate();
	    }

	    /* Evaluates the weight of the invariant mass generation. */

	    bool evaluate_s_weight()
	    {
		return s_gen->evaluate_weight();
	    }

	    /* Returns the weight corresponding to the generated incoming
	     * invariant mass. */

	    const value_type& s_weight() const
	    {
		return s_gen->weight();
	    }

	    /* Updates weights in the tree. */

	    void update()
	    {
		branching_multichannel.integrand()=this->integrand();
		branching_multichannel.update();
		s_gen->integrand()=this->integrand();
		s_gen->update();
	    }

	    /* Adapts the multichannel weights. */

	    void adapt_channels()
	    {
		branching_multichannel.adapt();
	    }

	    /* Adapts the vegas grids: */

	    void adapt_grids()
	    {
		s_gen->adapt();
	    }

	    /* Adapts multichannel weights and vegas grids: */

	    void adapt()
	    {
		branching_multichannel.adapt();
		s_gen->adapt();
	    }

	    /* Resets cross sections, adaptive grids and multichannel weights:
	     * */

	    void reset()
	    {
		s_gen->reset();
		branching_multichannel.reset();
	    }
	    
	    /* Throws out all branchings that have zero multichannel weight. */

	    void clean()
	    {
		/* TODO: recursively clean up... */
	    }

	    /* Evaluates the invariant mass from the momentum. */

	    value_type evaluate_s()
	    {
		return ps_channel->evaluate_s();
	    }

	    /* Evaluates the invariant mass from the momentum. */

	    value_type evaluate_m()
	    {
		return ps_channel->evaluate_m();
	    }

	    /* Refreshes internal generation parameters. */

	    bool refresh_params()
	    {
		return s_gen->refresh_params();
	    }

	    /* Set the status to s-generated: */

	    void set_status_s_generated()
	    {
		ps_channel->set_status_s_generated();
	    }

	    /* Set the status to p-generated: */

	    void set_status_p_generated()
	    {
		ps_channel->set_status_p_generated();
	    }

	    /* Reset the status: */

	    void reset_status()
	    {
		ps_channel->reset_status();
	    }

	    /* Public readout functions: */
	    /*---------------------------*/

	    /* Momentum channel readout: */
	    
	    bit_string_type bitstring() const
	    {
		return ps_channel->bitstring;
	    }

	    /* Returns a reference to the momentum flowing through the channel.
	     * */

	    momentum_type& p()
	    {
		return ps_channel->p();
	    }

	    /* Returns a const reference to the momentum flowing through the
	     * channel. */

	    const momentum_type& p() const
	    {
		return ps_channel->p();
	    }

	    /* Returns a reference to the i-th momentum component. */

	    value_type& p(size_type i)
	    {
		return ps_channel->p(i);
	    }

	    /* Returns a const reference to the i-th momentum component. */

	    const value_type& p(size_type i) const
	    {
		return ps_channel->p(i);
	    }

	    /* Returns a reference to the invariant mass-squared. */

	    value_type& s()
	    {
		return ps_channel->s();
	    }

	    /* Returns a const reference to the invariant mass-squared. */

	    const value_type& s() const
	    {
		return ps_channel->s();
	    }

	    /* Returns the signed invariant mass. */

	    value_type m() const
	    {
		return ps_channel->m();
	    }

	    /* Returns a const reference to the minimal invariant mass-squared.
	     * */

	    const value_type& s_min() const
	    {
		return s_gen->lower_bound();
	    }

	    /* Returns the signed minimal invariant mass. */

	    value_type m_min() const
	    {
		return sgn_sqrt(s_min());
	    }

	    /* Returns a const reference to the minimal minimal invariant
	     * mass-squared.  */

	    const value_type& s_min_min() const
	    {
		return ps_channel->s_min_min();
	    }

	    /* Returns the minimal minimal signed invariant mass. */

	    value_type m_min_min() const
	    {
		return ps_channel->m_min_min();
	    }

	    /* Returns a const reference to the maximal invariant mass-squared.
	     * */

	    const value_type& s_max() const
	    {
		return s_gen->upper_bound();
	    }

	    /* Returns a const reference to the maximal signed invariant mass.
	     * */

	    value_type m_max() const
	    {
		return sgn_sqrt(s_max());
	    }

	    /* Returns a const reference to the maximal maximal invariant
	     * mass-squared.  */

	    const value_type& s_max_max() const
	    {
		return ps_channel->s_max_max();
	    }

	    /* Returns the maximal maximal signed invariant mass.  */

	    value_type m_max_max() const
	    {
		return ps_channel->m_max_max();
	    }

	    /* Returns the minimal pT: */

	    value_type pT_min() const
	    {
		return ps_channel->pT_min();
	    }

	    /* Returns the maximal cosine theta: */

	    value_type ct_max() const
	    {
		return ps_channel->ct_max();
	    }

	    /* Returns whether the channel is external. */

	    bool on_shell() const
	    {
		return ps_channel->on_shell();
	    }

	    /* Returns whether the channel is internal. */

	    bool off_shell() const
	    {
		return ps_channel->off_shell();
	    }

	    /* Returns whether we have a t-type channel.*/

	    bool spacelike() const
	    {
		return ps_channel->spacelike();
	    }

	    /* Returns whether we have an s-type channel. */

	    bool timelike() const
	    {
		return ps_channel->timelike();
	    }

	    /* Returns the propagator sampling type */

	    std::string type() const
	    {
		return s_gen->type();
	    }

	    /* Returns the number of branchings */

	    size_type branching_count() const
	    {
		return branching_multichannel.channel_count();
	    }

	    /* Returns the i-th branching instance: */

	    const branching_type* branching(size_type i) const
	    {
		return static_cast<const branching_type*>(branching_multichannel.generator(i));
	    }

	    /* Returns the channel generation status for readout: */

	    status_type get_status() const
	    {
		return ps_channel->get_status();
	    }

	    /* Printing method. */
	    
	    std::ostream& print(std::ostream& os) const
	    {
		if(particle_type!=NULL)
		{
		    os<<std::setw(10)<<particle_type->get_name();
		}
		else
		{
		    os<<std::setw(10)<<"X";
		}
		if(branching_multichannel.channel_count()!=0)
		{
		    size_type i=0;
		    os<<"  --->  ";
		    branching(i)->print(os);
		    os<<"\t"<<branching_multichannel.alpha(i);
		    os<<std::endl;
		    ++i;
		    for(;i<branching_count();++i)
		    {
			os<<std::setw(10)<<' ';
			os<<"  --->  ";
			branching(i)->print(os);
			os<<"\t"<<branching_multichannel.alpha(i);
			os<<std::endl;
		    }
		    os<<std::setw(12)<<' ';
		}
		else
		{
		    os<<std::setw(2)<<' ';
		}
		os<<"recursive weight:"<<this->weight()<<", s generator type: "<<s_gen->detailed_type()<<std::endl;
		return os;
	    }

	    /* Other printing method: */

	    std::ostream& shortprint(std::ostream& os) const
	    {
		os<<name<<", p = "<<p()<<", m = "<<m()<<", [m--,m++] = ["<<this->m_min_min()<<','<<this->m_max_max()<<"], w = "<<this->weight();
		return os;
	    }

	    std::ostream& print_s_generator(std::ostream& os) const
	    {
		os<<s_gen->detailed_type()<<", s = "<<m()<<" [m--,m++] = ["<<this->m_min_min()<<','<<this->m_max_max()<<"]"<<std::endl;
		return os;
	    }

	    /* Warning utility: */

	    void s_warning() const
	    {
		if(!s_gen->normalisable())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-sampler "<<s_gen->detailed_type()<<" not normalizable within ["<<m_min()<<','<<m_max()<<"] in channel\n"<<*this<<endlog;
		}
	    }

	private:

	    /* Private static functions: */
	    /*---------------------------*/

	    /* Method that determines the id string: */

	    static std::string make_channel_name(const particle<model_t>* phi,const bit_string_type p)
	    {
		std::stringstream ss;
		if(phi==NULL)
		{
		    ss<<"X("<<p<<')';
		}
		else
		{
		    ss<<phi->get_name()<<'('<<p<<')';
		}
		return ss.str();
	    }

	    /* Private data members: */
	    /*-----------------------*/

	    /* Momentum channel: */

	    momentum_channel_type* const ps_channel;

	    /* Invariant mass generator: */

	    s_generator_type* s_gen;

	    /* Branchings multi-channel: */

	    multi_channel<value_type,rn_engine> branching_multichannel;

	    class same_branching
	    {
		public:

		    const branching_type* const b;

		    same_branching(const branching_type* b_):b(b_){}
		    bool operator()(const base_type* other) const
		    {
			return b->equiv(static_cast<const branching_type*>(other));
		    }
	    };
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t particle_channel<model_t,N_in,N_out,rng_t>::N_bits;

    /* Output stream overload: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>std::ostream& operator << (std::ostream& os,const particle_channel<model_t,N_in,N_out,rng_t>& c)
    {
	return c.shortprint(os);
    }

    /* Particle channel comparison class implementation: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class particle_channel_comp
    {
	public:

	    typedef particle_channel<model_t,N_in,N_out,rng_t> value_type;

	    flavour_comp< particle<model_t> > particle_comp;

	    particle_channel_comp(){}
	    bool operator ()(const value_type& a,const value_type& b) const
	    {
		return (a.bitstring()==b.bitstring())?(particle_comp(a.particle_type,b.particle_type)):(a.bitstring()<b.bitstring());
	    }
	    bool operator ()(const value_type* a,const value_type* b) const
	    {
		return (a->bitstring()==b->bitstring())?(particle_comp(a->particle_type,b->particle_type)):(a->bitstring()<b->bitstring());
	    }
    };

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>typename model_t::value_type integrate(const particle_channel<model_t,N_in,N_out,rng_t>* first,const particle_channel<model_t,N_in,N_out,rng_t>* second,const typename model_t::value_type& sqrts)
    {
	return integrate<typename model_t::value_type,rng_t>(first->s_gen,second->s_gen,sqrts);
    }
}

#endif /*CAMGEN_P_CHANNEL_H_*/

