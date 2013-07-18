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
#include <Camgen/MC_config.h>
#include <Camgen/ps_decl.h>
#include <Camgen/MC_gen.h>
#include <Camgen/sgen_grid.h>

namespace Camgen
{

    /* Integrates the s-densities within the allowed phase space: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>typename model_t::value_type integrate(const particle_channel<model_t,N_in,N_out,rng_t>*,const particle_channel<model_t,N_in,N_out,rng_t>*,const typename model_t::value_type&);
    
    /* Phase space particle channel class: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class particle_channel: public MC_generator<typename model_t::value_type>
    {
	friend typename model_t::value_type integrate<model_t,N_in,N_out,rng_t>(const particle_channel<model_t,N_in,N_out,rng_t>*,const particle_channel<model_t,N_in,N_out,rng_t>*,const typename model_t::value_type&);
	friend class momentum_channel<model_t,N_in,N_out,rng_t>;
	friend class ps_factory<model_t,N_in,N_out,rng_t,typename model_t::spacetime_type>;
	friend class ps_tree<model_t,N_in,N_out,rng_t>;
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
	    typedef s_generator<value_type,rng_t> s_generator_type;
	    typedef BW_s_generator<value_type,rng_t> Breit_Wigner_type;
	    typedef pl_s_generator<value_type,rng_t> power_law_type;
	    typedef uni_s_generator<value_type,rng_t> uniform_type;
	    typedef Dd_s_generator<value_type,rng_t> Dirac_delta_type;
	    static const std::size_t N_bits=N_in+N_out-1;	    
	    typedef bit_string<N_bits> bit_string_type;
	    
	    /* Public data members: */
	    /*----------------------*/

	    /* Particle type: */

	    const particle<model_t>* const particle_type;

	    /* Channel name: */

	    const std::string name;

	    /* Dirac delta generation flag: */

	    const bool Dirac_delta;

	    /* Destructor: */
	    /*-------------*/

	    virtual ~particle_channel()
	    {
		delete s_gen;
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Sets the minimal invariant mass-squared. */
	    
	    bool set_s_min(const value_type& smin)
	    {
		return s_gen->set_s_min(std::max(ps_channel->s_min_min(),smin));
	    }

	    /* Sets the minimal minimal invariant mass-squared. */
	    
	    bool set_s_min_min(const value_type& sminmin)
	    {
		bool q=s_gen->set_s_min_min(sminmin);
		if(sminmin>this->s_min())
		{
		    q&=(s_gen->set_s_min(sminmin));
		}
		return q;
	    }

	    /* Sets the minimal signed invariant mass. */
	    
	    bool set_m_min(const value_type& mmin)
	    {
		return s_gen->set_m_min(std::max(ps_channel->m_min_min(),mmin));
	    }

	    /* Sets the maximal invariant mass-squared. */
	    
	    bool set_s_max(const value_type& smax)
	    {
		return s_gen->set_s_max(std::min(ps_channel->s_max_max(),smax));
	    }

	    /* Sets the maximal maximal invariant mass-squared. */
	    
	    bool set_s_max_max(const value_type& smaxmax)
	    {
		bool q=s_gen->set_s_max_max(smaxmax);
		if(smaxmax<this->s_max())
		{
		    q&=s_gen->set_s_max(smaxmax);
		}
		return q;
	    }

	    /* Sets the maximal signed invariant mass. */
	    
	    bool set_m_max(const value_type& mmax)
	    {
		return s_gen->set_m_max(std::min(ps_channel->m_max_max(),mmax));
	    }

	    /* Sets the invariant mass-squared range. */
	    
	    bool set_s_range(const value_type& smin,const value_type& smax)
	    {
		if(s_gen->set_s_range(std::max(ps_channel->s_min_min(),smin),std::min(ps_channel->s_max_max(),smax)))
		{
		    return true;
		}
		return false;
	    }

	    /* Sets the invariant mass range. */
	    
	    bool set_m_range(const value_type& mmin,const value_type& mmax)
	    {
		return s_gen->set_m_range(std::max(ps_channel->m_min_min(),mmin),std::min(ps_channel->m_max_max(),mmax));
	    }

	    /* Sets the multichannel adaptivity and threshold to the argument
	     * pair: */

	    void set_multichannel_params(const std::pair<value_type,value_type>& x)
	    {
		adaptivity=&(x.first);
		threshold=&(x.second);
	    }

	    /* Recursively picks a channel according to the multichannel weight of the
	     * branchings, stores the selected branching in the argument vector,
	     * and proceeds down the tree. Returns false upon failure. */

	    bool choose_branching(std::vector<branching_type*>& b)
	    {
		generation_flag=true;
		if(branchings.size()==0)
		{
		    return true;
		}
		if(branchings.size()==1)
		{
		    branching_iterator=branchings.begin();
		    b.push_back(*branching_iterator);
		    (*branching_iterator)->choose_channel(b);
		    return true;
		}
		else
		{
		    branching_iterator=branchings.begin();
		    value_type rho=rn_stream::throw_number();
		    value_type r(0);
		    while(r<rho and branching_iterator!=branchings.end())
		    {
			r+=((*branching_iterator)->alpha);
			++branching_iterator;
		    }
		    --branching_iterator;
		}
		if(branching_iterator==branchings.end())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"multichannel weights incorrectly normalized: channel iterator overflow encountered"<<endlog;
		    b.clear();
		    return false;
		}
		b.push_back(*branching_iterator);
		return (*branching_iterator)->choose_channel(b);
	    }

	    /* Calls the selected branching generation method. */

	    bool generate()
	    {
		if(branchings.size()==0)
		{
		    return true;
		}
		return (*branching_iterator)->generate();
	    }

	    /* Generates the incoming invariant mass. */

	    bool generate_s()
	    {
		if(s_gen->generate())
		{
		    ps_channel->evaluate_m();
		    return true;
		}
		return false;
	    }

	    /* Evaluates the weight of the invariant mass generation. */

	    bool evaluate_s_weight()
	    {
		return s_gen->evaluate_weight();
	    }

	    /* Retruns whether the s-generation is normalisable: */

	    bool s_normalisable()
	    {
		return s_gen->normalisable();
	    }

	    /* Normalisation checking utility: */

	    std::ostream& check_s_norm(std::ostream& os) const
	    {
		if(!s_gen->normalisable())
		{
		    os<<"in particle channel "<<name<<": s-gen "<<s_gen->detailed_type();
		    os<<" not normalisable within ["<<s_gen->s_min()<<','<<s_gen->s_max()<<"]..."<<std::endl;
		}
		return os;
	    }

	    /* Returns the weight corresponding to the generated incoming
	     * invariant mass. */

	    const value_type& s_weight() const
	    {
		return s_gen->weight();
	    }

	    /* Weight evaluation method. */

	    bool evaluate_weight()
	    {
		if(branchings.size()==0)
		{
		    this->weight()=(value_type)1;
		    return true;
		}
		value_type g(0);
		for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		{
		    if((*it)->alpha!=(value_type)0)
		    {
			if(!(*it)->evaluate_weight())
			{
			    this->weight()=(value_type)0;
			    return false;
			}
			g+=((*it)->alpha/(*it)->weight());
		    }
		}
		this->weight()=(value_type)1/g;
		return true;
	    }

	    /* Resets the generation flag. */

	    void reset_generation_flags()
	    {
		generation_flag=false;
		if(branching_iterator!=branchings.end())
		{
		    (*branching_iterator)->reset_generation_flags();
		}
	    }

	    /* Refreshes the estimators of the individual channels. */

	    void update()
	    {
		if(this->weight()==(value_type)0)
		{
		    return;
		}
		if(branchings.size()>1)
		{
		    value_type w3=std::pow(this->integrand(),(int)2)*(this->weight());
		    if(w3!=(value_type)0)
		    {
			multichannel_flag=true;
		    }
		    for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		    {
			if((*it)->alpha!=(value_type)0 and (*it)->weight()!=(value_type)0)
			{
			    (*it)->W+=(w3/(*it)->weight());
			}
		    }
		}
		++update_counter;
		s_gen->integrand()=this->integrand();
		s_gen->update_weight();
	    }

	    /* Adapts the multichannel weights. */

	    void adapt_channels()
	    {
		if(branchings.size()>1 and multichannel_flag)
		{
		    value_type D(0);
		    for(typename std::list<branching_type*>::iterator it1=branchings.begin();it1!=branchings.end();++it1)
		    {
			if((*it1)->alpha!=(value_type)0)
			{
			    (*it1)->W/=update_counter;
			    for(typename std::list<branching_type*>::iterator it2=branchings.begin();it2!=it1;++it2)
			    {
				if((*it2)->alpha!=(value_type)0)
				{
				    D=std::max(D,std::abs((*it1)->W-(*it2)->W));
				}
			    }
			}
		    }
		    value_type norm(0),beta((value_type)0.5*multichannel_adaptivity());
		    for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		    {
			if((*it)->alpha!=(value_type)0)
			{
			    (*it)->alpha*=std::pow((*it)->W,beta);
			    norm+=(*it)->alpha;
			}
		    }
		    value_type thresh=multichannel_threshold()*norm;
		    value_type normbis(0);
		    for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		    {
			if((*it)->alpha<thresh/branchings.size())
			{
			    (*it)->alpha=(value_type)0;
			}
			normbis+=(*it)->alpha;
		    }
		    for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		    {
			if((*it)->alpha!=(value_type)0)
			{
			    (*it)->alpha/=normbis;
			}
			(*it)->W=(value_type)0;
		    }
		    branchings.sort(branching_ordering);
		}
		update_counter=0;
		multichannel_flag=false;
	    }

	    /* Sorts the branchings descending in multichannel weight: */

	    void sort_branchings()
	    {
		value_type norm(0);
		for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		{
		    norm+=(*it)->alpha;
		}
		if(!(norm>(value_type)0))
		{
		    value_type a=(value_type)1/branchings.size();
		    for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		    {
			(*it)->alpha=a;
		    }
		}
		else
		{
		    for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		    {
			(*it)->alpha/=norm;
		    }
		    branchings.sort(branching_ordering);
		}
	    }

	    /* Adapts the vegas grids: */

	    void adapt_grids()
	    {
		s_gen->adapt();
	    }

	    /* Adapts multichannel weights and vegas grids: */

	    void adapt()
	    {
		adapt_grids();
		adapt_channels();
	    }

	    /* Resets cross sections, adaptive grids and multichannel weights:
	     * */

	    void reset()
	    {
		this->base_type::reset();
		for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		{
		    (*it)->W=(value_type)0;
		    (*it)->alpha=(value_type)1/branchings.size();
		    (*it)->reset();
		}
		generation_flag=false;
		update_counter=0;
		multichannel_flag=false;
		s_gen->reset();
	    }

	    /* Resets cross sections of branchings and s-generator: */

	    void reset_cross_section()
	    {
		this->base_type::reset_cross_section();
		for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		{
		    (*it)->reset_cross_section();
		}
		s_gen->reset_cross_section();
	    }
	    
	    static bool redundant_branching(const branching_type* c)
	    {
		return c->is_redundant();
	    }

	    /* Throws out all branchings that have zero multichannel weight. */

	    void clean()
	    {
		if(generation_flag)
		{
		    reset_generation_flags();
		}
	    	branchings.remove_if(redundant_branching);
		branching_iterator=branchings.begin();
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
		return *(ps_channel->momentum);
	    }

	    /* Returns a const reference to the momentum flowing through the
	     * channel. */

	    const momentum_type& p() const
	    {
		return *(ps_channel->momentum);
	    }

	    /* Returns a reference to the i-th momentum component. */

	    value_type& p(size_type i)
	    {
		return (*(ps_channel->momentum))[i];
	    }

	    /* Returns a const reference to the i-th momentum component. */

	    const value_type& p(size_type i) const
	    {
		return (*(ps_channel->momentum))[i];
	    }

	    /* Returns a reference to the invariant mass-squared. */

	    value_type& s()
	    {
		return ps_channel->invariant_mass;
	    }

	    /* Returns a const reference to the invariant mass-squared. */

	    const value_type& s() const
	    {
		return ps_channel->invariant_mass;
	    }

	    /* Returns a reference to the signed invariant mass. */

	    const value_type& m() const
	    {
		return ps_channel->signed_mass;
	    }

	    /* Returns a const reference to the minimal invariant mass-squared.
	     * */

	    const value_type& s_min() const
	    {
		return s_gen->s_min();
	    }

	    /* Returns a const reference to the minimal minimal invariant
	     * mass-squared.  */

	    const value_type& s_min_min() const
	    {
		return ps_channel->sminmin;
	    }

	    /* Returns a const reference to the minimal signed invariant mass.
	     * */

	    const value_type& m_min() const
	    {
		return s_gen->m_min();
	    }

	    /* Returns a const reference to the minimal minimal signed invariant
	     * mass.  */

	    const value_type& m_min_min() const
	    {
		return ps_channel->mminmin;
	    }

	    /* Returns a const reference to the maximal invariant mass-squared.
	     * */

	    const value_type& s_max() const
	    {
		return s_gen->s_max();
	    }

	    /* Returns a const reference to the maximal maximal invariant
	     * mass-squared.  */

	    const value_type& s_max_max() const
	    {
		return ps_channel->smaxmax;
	    }

	    /* Returns a const reference to the maximal signed invariant mass.
	     * */

	    const value_type& m_max() const
	    {
		return s_gen->m_max();
	    }

	    /* Returns a const reference to the maximal maximal signed invariant
	     * mass.  */

	    const value_type& m_max_max() const
	    {
		return ps_channel->mmaxmax;
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

	    /* Multichannel adaptivity: */

	    value_type multichannel_adaptivity() const
	    {
		return (adaptivity==NULL)?(value_type)1:(*adaptivity);
	    }

	    /* Multichannel threshold: */

	    value_type multichannel_threshold() const
	    {
		if(threshold==NULL)
		{
		    return (value_type)0;
		}
		if(*threshold<(value_type)0)
		{
		    return (value_type)0;
		}
		if(*threshold>(value_type)1)
		{
		    return (value_type)1;
		}
		return *threshold;
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
		if(branchings.size()!=0)
		{
		    typename std::list<branching_type*>::const_iterator it=branchings.begin();
		    os<<"  --->  ";
		    (*it)->print(os);
		    os<<std::endl;
		    ++it;
		    for(;it!=branchings.end();++it)
		    {
			os<<std::setw(10)<<' ';
			os<<"  --->  ";
			(*it)->print(os);
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

	    /* Returns the redundancy flag. */

	    bool is_redundant() const
	    {
		return redundant;
	    }

	    /* Recursively sets the redundant flag: */
	    
	    void set_redundant()
	    {
		redundant=true;
	    }

	    /* Recursively unsets the redundant flag: */
	    
	    void unset_redundant()
	    {
		redundant=false;
		for(typename std::list<branching_type*>::iterator it=branchings.begin();it!=branchings.end();++it)
		{
		    if((*it)->alpha>(value_type)0)
		    {
			(*it)->unset_redundant();
		    }
		}
	    }

	    /* Warning utility: */

	    void s_warning() const
	    {
		if(!s_gen->normalisable())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-sampler "<<s_gen->detailed_type()<<" not normalizable within ["<<s_gen->m_min()<<','<<s_gen->m_max()<<"] in channel\n"<<*this<<endlog;
		}
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Loads the instance from input stream: */

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		is>>update_counter>>multichannel_flag;
		return is;
	    }

	    /* Writes the instance to the output stream: */

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<particle>"<<std::endl;
		std::string str=(particle_type==NULL)?"X":particle_type->get_name();
		os<<str<<"\t"<<nu<<"\t"<<redundant<<std::endl;
		s_gen->save(os);
		this->base_type::save(os);
		os<<update_counter<<"\t"<<multichannel_flag<<std::endl;
		os<<"</particle>"<<std::endl;
		return os;
	    }

	private:

	    /* Private types: */
	    /*----- ----------*/

	    class same_branching
	    {
		public:

		    const branching_type* const b;

		    same_branching(const branching_type* b_):b(b_){}
		    bool operator()(const branching_type* other) const
		    {
			return b->equiv(other);
		    }
	    };

	    /* Private static functions: */
	    /*---------------------------*/

	    /* Inspects whether the momentum channel induces a Dirac-delta-like particle
	     * channel: */

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
	    static bool is_Dirac_delta(const momentum_channel_type* c,const value_type* E)
	    {
		return (c->on_shell() or (c->timelike() and (c->bitstring.count()==N_out) and (E!=NULL)));
	    }
	    static bool branching_ordering(const branching_type* first,const branching_type* second)
	    {
		return (first->alpha>second->alpha);
	    }

	    /* Private data members: */
	    /*-----------------------*/

	    /* Momentum channel: */

	    momentum_channel_type* const ps_channel;

	    /* Invariant mass generator: */

	    s_generator_type* s_gen;

	    /* Power-law exponent: */

	    value_type nu;

	    /* Branchings starting at this channel: */
	    
	    std::list<branching_type*>branchings;

	    /* Sampled branching by the multichannel: */

	    typename std::list<branching_type*>::iterator branching_iterator;
	    
	    /* Update counters: */
	    
	    size_type update_counter;

	    /* Flag setting nonzero multichannel weight estimators: */

	    bool multichannel_flag;

	    /* Generation flag: */

	    bool& generation_flag;

	    /* Cleaning flag: */

	    bool redundant;

	    /* Multichannel adaptivity: */

	    const value_type* adaptivity;

	    /* Multichannel threshold: */

	    const value_type* threshold;

	    /* Private constructors: */
	    /*-----------------------*/

	    /* Momentum-allocating constructor. */

	    particle_channel(momentum_channel_type* ps_channel_,const particle<model_t>* particle_type_,const value_type* Ecmhat=NULL):particle_type(particle_type_),name(make_channel_name(particle_type_,ps_channel_->bitstring)),Dirac_delta(is_Dirac_delta(ps_channel_,Ecmhat)),ps_channel(ps_channel_),s_gen(NULL),nu(0),branching_iterator(branchings.end()),update_counter(0),multichannel_flag(false),generation_flag(ps_channel_->generation_flag),redundant(true),adaptivity(NULL),threshold(NULL)
	    {
		value_type* s=&(ps_channel->invariant_mass);
		bool adaptive=false;
		bool total_s=(ps_channel_->timelike() and ps_channel_->bitstring.count()==N_out);
		s_generator<value_type,rng_t>* tempgen;
		const value_type* mass=(particle_type==NULL)?(NULL):(particle_type->get_mass_address());
		if(on_shell())
		{
		    tempgen=new Dirac_delta_type(s,mass);
		    adaptive=false;
		}
		else if(total_s and Ecmhat!=NULL)
		{
		    tempgen=new Dirac_delta_type(s,Ecmhat);
		    adaptive=false;
		}
		else if(particle_type==NULL)
		{
		    adaptive=(ps_channel->timelike())?adaptive_s_sampling():adaptive_t_sampling();
		    if(total_s)
		    {
			nu=shat_exponent();
		    }
		    else
		    {
			nu=auxiliary_exponent();
		    }
		    if(nu==(value_type)0)
		    {
			tempgen=new uniform_type(s);
		    }
		    else
		    {
			tempgen=new power_law_type(s,NULL,&nu);
		    }
		}
		else
		{
		    adaptive=(ps_channel->timelike())?adaptive_s_sampling():adaptive_t_sampling();
		    const value_type* mass=particle_type->get_mass_address();
		    const value_type* width=particle_type->get_width_address();
		    if(particle_type->is_auxiliary())
		    {
			tempgen=new uniform_type(s);
		    }
		    else if(spacelike() or (width==NULL))
		    {
			nu=spacelike()?spacelike_exponent(particle_type->get_name()):timelike_exponent(particle_type->get_name());
			tempgen=new power_law_type(s,mass,&nu);
		    }
		    else
		    {
			tempgen=new Breit_Wigner_type(s,mass,width);
		    }
		}
		if(adaptive)
		{
		    s_gen=new adaptive_s_generator<value_type,rng_t>(s,tempgen,grid_bins(),grid_mode());
		}
		else
		{
		    s_gen=tempgen;
		    if(Dirac_delta)
		    {
			set_m_min((mass==NULL)?0:(*mass));
			set_m_max((mass==NULL)?0:(*mass));
		    }
		}
	    }

	    /* Private modifiers: */
	    /*--------------------*/

	    /* Adds a branching instance to the multichannel: */

	    branching_type* insert_branching(branching_type* br)
	    {
		if(ps_channel->bitstring.count()==1)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to attach phase space branching to external channel ignored--returning NULL"<<endlog;
		    return NULL;
		}
		if(branchings.size()==0)
		{
		    branchings.push_back(br);
		    br->alpha=(value_type)1;
		    return br;
		}
		same_branching pred(br);
		typename std::list<branching_type*>::iterator it=std::find_if(branchings.begin(),branchings.end(),pred);
		if(it==branchings.end())
		{
		    branchings.push_back(br);
		    value_type alpha=(value_type)1/branchings.size();
		    for(typename std::list<branching_type*>::iterator it2=branchings.begin();it2!=branchings.end();++it2)
		    {
			(*it2)->alpha=alpha;
		    }
		    return br;
		}
		delete br;
		return *it;
	    }
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

