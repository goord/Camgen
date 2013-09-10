//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_tree.h
    \brief Phase space generator classes using the recursive phase space decomposition.
 */

#ifndef CAMGEN_PS_TREE_H_
#define CAMGEN_PS_TREE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Phase space multichannel tree class definition. Creates a recursive tree of *
 * phase space branchings and channels from a Camgen off-shell current tree.   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <Camgen/rn_strm.h>
#include <Camgen/ps_fac.h>
#include <Camgen/ps_gen.h>
#include <Camgen/rambo.h>

namespace Camgen
{
    /* Phase space decomposition tree generator for decay processes: */

    template<class model_t,std::size_t N_out,class rng_t>class ps_tree<model_t,1,N_out,rng_t>: public ps_generator<model_t,1,N_out>
    {
	typedef ps_generator<model_t,1,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::init_state_type init_state_type;

	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    typedef momentum_channel<model_t,1,N_out,rng_t> momentum_channel_type;
	    typedef particle_channel<model_t,1,N_out,rng_t> particle_channel_type;
	    typedef ps_branching<model_t,1,N_out,rng_t> branching_type;
	    typedef typename std::list<momentum_channel_type*> momentum_channel_container;
	    typedef typename std::list<particle_channel_type*> particle_channel_container;
	    typedef typename std::list<branching_type*> branching_container;
	    typedef ps_factory<model_t,1,N_out,rng_t,typename model_t::spacetime_type> factory_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;

	    /* Static utility functions: */

	    static bool necessary_branching(const branching_type* c)
	    {
		return !(c->is_redundant());
	    }
	    static bool non_empty_channel(const momentum_channel_type* c)
	    {
		return c->any();
	    }

	    /* Utility classes: */

	    class bitstring_eq
	    {
		public:

		    const bit_string_type bitstring;

		    bitstring_eq(const bit_string_type& bitstring_):bitstring(bitstring_){}

		    bool operator()(const momentum_channel_type* p)
		    {
			return (p->bitstring==bitstring);
		    }
	    };
	    
	    /* Public static methods */
	    /*-----------------------*/

	    /* Factory method. */

	    static ps_tree<model_t,1,N_out,rng_t>* create_instance(typename CM_algorithm<model_t,1,N_out>::tree_iterator it,init_state_type* is_)
	    {
		typedef typename CM_algorithm<model_t,1,N_out>::tree_type CM_tree_type;
		typedef typename CM_tree_type::current_type current_type;
		typedef typename current_type::particle_type particle_type;
		typedef typename CM_tree_type::current_tree_type current_tree_type;
		
		if(!is_->s_hat_sampling)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid request to let 1 -> "<<N_out<<" ps-tree generate partonic invariant mass--returning NULL"<<endlog;
		    return NULL;
		}
		if(current_tree_type::final_current()!=0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"phase space tree generator is only implemented for final current being zero--returning NULL"<<endlog;
		    return NULL;
		}
		ps_tree<model_t,1,N_out,rng_t>* result=new ps_tree<model_t,1,N_out,rng_t>(is_);
		if(!result->set_tree(it))
		{
		    return NULL;
		}
		bit_string_type b;
		b.set(0);
		const value_type* sqrtshat=(result->s_hat_sampling)?(NULL):(result->get_Ecm_hat());
		for(size_type i=0;i<N_out;++i)
		{
		    result->outgoing_momentum_channels[i]=new momentum_channel_type(b,&(result->p_out(i)));
		    result->outgoing_momentum_channels[i]->assign_mmin_components(result->eff_mmin);
		    result->outgoing_particle_channels[i]=result->outgoing_momentum_channels[i]->add_particle_channel(it->get_phase_space(i+1)->particle_type,sqrtshat);
		    result->outgoing_particle_channels[i]->set_multichannel_params(result->multichannel_params);
		    b>>=1;
		}
		b.set();
		result->incoming_momentum_channel=new momentum_channel_type(b,&(result->p_in(0)));
		result->incoming_momentum_channel->assign_mmin_components(result->eff_mmin);
		result->incoming_particle_channel=result->incoming_momentum_channel->add_particle_channel(it->get_phase_space(0)->particle_type,sqrtshat);
		result->incoming_particle_channel->set_multichannel_params(result->multichannel_params);
		for(typename CM_tree_type::const_interaction_iterator v_it=it->interactions_begin();v_it!=it->interactions_end();++v_it)
		{
		    if(v_it->is_coupled())
		    {
			bit_string<N_out>bs_in(v_it->get_produced_bit_string());
			std::vector< bit_string<N_out> >bs_out(v_it->get_incoming_bit_strings());
			const particle_type* phi_in=v_it->get_produced_particle();
			std::vector<const particle_type*>phi_out(v_it->get_incoming_particles());
			if(v_it->get_rank()==3)
			{
			    result->insert_branching(bs_in,phi_in->get_anti_particle(),bs_out[0],phi_out[0]->get_anti_particle(),bs_out[1],phi_out[1]->get_anti_particle());
			}
			else if(v_it->get_rank()==4)
			{
			    result->insert_branching((bs_out[0] | bs_out[1]),NULL,bs_out[0],phi_out[0]->get_anti_particle(),bs_out[1],phi_out[1]->get_anti_particle());
			    result->insert_branching(bs_in,phi_in->get_anti_particle(),(bs_out[0] | bs_out[1]),NULL,bs_out[2],phi_out[2]->get_anti_particle());
			    result->insert_branching((bs_out[0] | bs_out[2]),NULL,bs_out[0],phi_out[0]->get_anti_particle(),bs_out[2],phi_out[2]->get_anti_particle());
			    result->insert_branching(bs_in,phi_in->get_anti_particle(),(bs_out[0] | bs_out[2]),NULL,bs_out[1],phi_out[1]->get_anti_particle());
			    result->insert_branching((bs_out[1] | bs_out[2]),NULL,bs_out[1],phi_out[1]->get_anti_particle(),bs_out[2],phi_out[2]->get_anti_particle());
			    result->insert_branching(bs_in,phi_in->get_anti_particle(),(bs_out[1] | bs_out[2]),NULL,bs_out[0],phi_out[0]->get_anti_particle());
			}
			else
			{
			    std::cout<<"Phase space branching for "<<v_it->get_rank()<<"-rank vertex not (yet) implemented..."<<std::endl;
			}
		    }
		}
		result->clean();
		return result;
	    }
	    
	    /* Public data members */
	    /*---------------------*/

	    /* Const flags. */
	    
	    const bool backward_s_gen;
	    const bool adaptive_s;
	    const bool adaptive_theta;
	    const size_type max_grid_bins;

	    /* Maximal number of s-pair generations: */

	    size_type max_s_pairs;

	    /* Multichannel adaptivity and threshold: */

	    std::pair<value_type,value_type> multichannel_params;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Constructor. */

	    ps_tree(init_state_type* is_):base_type(is_),backward_s_gen(backward_s_sampling()),adaptive_s(adaptive_s_sampling()),adaptive_theta(adaptive_angles()),max_grid_bins(grid_bins()),max_s_pairs(Camgen::max_s_pairs()),multichannel_params(std::pair<value_type,value_type>(multichannel_adaptivity(),multichannel_threshold())),incoming_momentum_channel(NULL),incoming_particle_channel(NULL)
	    {
		outgoing_momentum_channels.assign(NULL);
		outgoing_particle_channels.assign(NULL);
	    }

	    /* Configurable constructor. */

	    ps_tree(init_state_type* is_,bool q1,bool q2,bool q3,size_type b):base_type(is_),backward_s_gen(q1),adaptive_s(q2),adaptive_theta(q3),max_grid_bins(b),max_s_pairs(Camgen::max_s_pairs()),multichannel_params(std::pair<value_type,value_type>(multichannel_adaptivity(),multichannel_threshold())),incoming_momentum_channel(NULL),incoming_particle_channel(NULL)
	    {
		outgoing_momentum_channels.assign(NULL);
		outgoing_particle_channels.assign(NULL);
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /* Destructor: */

	    ~ps_tree()
	    {
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    delete *it;
		}
		if(incoming_momentum_channel!=NULL)
		{
		    delete incoming_momentum_channel;
		}
		for(size_type i=0;i<N_out;++i)
		{
		    if(outgoing_momentum_channels[i]!=NULL)
		    {
			delete outgoing_momentum_channels[i];
		    }
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    delete *it;
		}
	    }

	    /* Public modifiers */
	    /*------------------*/

	    /* Channel choosing method: */

	    bool choose_channel()
	    {
		branching_sequence.clear();
		branching_sequence.reserve(2*(1+N_out));
		incoming_particle_channel->reset_generation_flags();
		return incoming_particle_channel->choose_branching(branching_sequence);
	    }

	    /* Partonic CM-energy generation: */

	    bool generate_s_hat(value_type& s)
	    {
		if(backward_s_gen)
		{
		    if(choose_channel())
		    {
			for(typename std::vector<branching_type*>::reverse_iterator it=branching_sequence.rbegin();it!=branching_sequence.rend();++it)
			{
			    if(!(*it)->generate_s())
			    {
				return false;
			    }
			}
			s=incoming_particle_channel->s();
			return true;
		    }
		    return false;
		}
		return true;
	    }

	    /* Final-state generation: */

	    bool generate_fs()
	    {
		if(!choose_channel())
		{
		    this->fsw=(value_type)0;
		    return false;
		}
		if(backward_s_gen)
		{
		    for(typename std::vector<branching_type*>::reverse_iterator it=branching_sequence.rbegin();it!=branching_sequence.rend();++it)
		    {
			if(!(*it)->generate_s())
			{
			    this->fsw=(value_type)0;
			    return false;
			}
		    }
		    for(typename std::vector<branching_type*>::iterator it=branching_sequence.begin();it!=branching_sequence.end();++it)
		    {
			if(!(*it)->generate_p())
			{
			    this->fsw=(value_type)0;
			    return false;
			}
		    }
		    return evaluate_fs_weight();
		}
		else if(incoming_particle_channel->evaluate_s() and incoming_particle_channel->generate())
		{
		    return evaluate_fs_weight();
		}
		this->fsw=(value_type)0;
		return false;
	    }

	    /* Weight evaluation method. */

	    bool evaluate_fs_weight()
	    {
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if(!(*it)->on_shell() and !(*it)->is_generated())
		    {
			(*it)->p().assign(0);
			for(size_type i=0;i<N_out;++i)
			{
			    if((*it)->bitstring[i])
			    {
				(*it)->p()+=(outgoing_momentum_channels[i]->p());
			    }
			}
			(*it)->evaluate_s();
		    }
		}
		bool q=true;
		for(size_type i=0;i<N_out;++i)
		{
		    q&=(outgoing_particle_channels[i]->evaluate_weight());
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if(!(*it)->on_shell())
		    {
			q&=((*it)->evaluate_weight());
		    }
		}
		q&=(incoming_particle_channel->evaluate_weight());
		this->fsw=q?(incoming_particle_channel->weight()):(value_type)0;
		return q;
	    }

	    /* Refreshes the minimal invariant masses: */

	    bool refresh_m_min()
	    {
		this->base_type::refresh_m_min();
		bool q=true;
		for(size_type i=0;i<N_out;++i)
		{
		    q&=(outgoing_momentum_channels[i]->refresh_m_min(this->s_tot(),this->s_tot(),this->s_tot()));
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    q&=((*it)->refresh_m_min(this->s_tot(),this->s_tot(),this->s_tot()));
		}
		q&=incoming_momentum_channel->refresh_m_min(this->s_tot(),this->s_tot(),this->s_tot());
		if(q)
		{
		    this->init_state()->set_m_hat_min(incoming_momentum_channel->m_min_min());
		}
		return q;
	    }
	    
	    /* Refreshes the total invariant mass and all depending objects: */
	    
	    bool refresh_Ecm()
	    {
		if(!this->init_state()->refresh_Ecm())
		{
		    return false;
		}
		bool q=true;
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    q&=((*it)->set_m_max_max(this->Ecm_hat()));
		}
		return q;
	    }

	    /* Refreshes the internal parameters: */

	    bool refresh_params()
	    {
		bool q=this->base_type::refresh_params();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    q&=((*it)->refresh_params());
		}
		return q;
	    }

	    /* Sets the minimal invariant mass i,j and refreshes depending
	     * channels: */

	    bool set_m_min(int i,int j,const value_type* sqrts)
	    {
		if(this->base_type::set_m_min(i,j,sqrts))
		{
		    bool q=true;
		    for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		    {
			if((*it)->bitstring[i-1] and (*it)->bitstring[j-1])
			{
			    q&=((*it)->refresh_m_min(this->s_tot(),this->s_tot(),this->s_tot()));
			}
		    }
		    q&=incoming_momentum_channel->refresh_m_min();
		    if(q)
		    {
			this->init_state()->set_m_hat_min(incoming_momentum_channel->m_min_min());
		    }
		    return q;
		}
		return false;
	    }

	    /* Overrides the update method: */

	    void update_fs()
	    {
		value_type f=(this->integrand())*(this->weight());
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    if((*it)->alpha!=(value_type)0)
		    {
			(*it)->integrand()=f;
			(*it)->update();
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->update(f);
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->update(f);
		}
		incoming_momentum_channel->update(f);
	    }

	    /* Overrides the vegas-grid adaptation method: */

	    void adapt_grids()
	    {
		this->base_type::adapt_grids();
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    if((*it)->alpha!=(value_type)0)
		    {
			(*it)->adapt_grids();
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->adapt_grids();
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->adapt_grids();
		}
		incoming_momentum_channel->adapt_grids();
	    }
	    
	    /* Overrides the multichannel adaptation method: */

	    void adapt_channels()
	    {
		this->base_type::adapt_channels();
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    if((*it)->alpha!=(value_type)0)
		    {
			(*it)->adapt_channels();
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->adapt_channels();
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->adapt_channels();
		}
		incoming_momentum_channel->adapt_channels();
		clean();
	    }

	    /* Resets cross-sections, multichannel weights and adaptive grids:
	     * */

	    void reset()
	    {
		this->base_type::reset();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->reset();
		}
	    }

	    /* Resets cross-sections: * */

	    void reset_cross_section()
	    {
		this->base_type::reset_cross_section();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->reset_cross_section();
		}
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Type identifier: */

	    std::string type() const
	    {
		if(backward_s_gen)
		{
		    return "pstree~";
		}
		else
		{
		    return "pstree";
		}
	    }
	    
	    /* Printing method. */
	    
	    std::ostream& print(std::ostream& os) const
	    {
		incoming_momentum_channel->print(os);
		for(typename momentum_channel_container::const_reverse_iterator it=momentum_channels.rbegin();it!=momentum_channels.rend();++it)
		{
		    (*it)->print(os);
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->print(os);
		}
		return os;
	    }

	    /// Virtual method printing the fs-generator settings.

	    std::ostream& print_fs_settings(std::ostream& os) const
	    {
		os<<std::setw(30)<<std::left<<"s sampling dir:"<<(backward_s_gen?"backward":"forward")<<std::endl;
		os<<std::setw(30)<<std::left<<"adaptive s-sampling:"<<(adaptive_s?"yes":"no")<<std::endl;
		os<<std::setw(30)<<std::left<<"adaptive angles:"<<(adaptive_theta?"yes":"no")<<std::endl;
		if(adaptive_s or adaptive_theta)
		{
		    os<<std::setw(30)<<std::left<<"max nr bins/grid:"<<max_grid_bins<<std::endl;
		}
		os<<std::setw(30)<<std::left<<"channel adaptivity:"<<multichannel_adaptivity()<<std::endl;
		os<<std::setw(30)<<std::left<<"channel threshold:"<<multichannel_threshold()<<std::endl;
		return os;
	    }

	    /* Generation channel printing method. */

	    std::ostream& print_channel(std::ostream& os) const
	    {
		if(branching_sequence.size()>0)
		{
		    for(size_type i=0;i<branching_sequence.size();++i)
		    {
			branching_sequence[i]->print(os);
			os<<std::endl;
		    }
		}
		else
		{
		    os<<"No branching sequence generated."<<std::endl;
		}
		return os;
	    }

	    /* Overridden loading method: */

	    std::istream& load_data(std::istream& is)
	    {
		is>>max_s_pairs;
		safe_read(is,multichannel_params.first);
		safe_read(is,multichannel_params.second);
		std::string flag;
		while(flag!="<channel>" and flag!="</fsgen>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(flag!="<channel>")
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"zero momentum channels read from file"<<endlog;
		    return is;
		}
		if(momentum_channels.size()!=0)
		{
		    for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		    {
			delete (*it);
		    }
		    momentum_channels.clear();
		}
		const value_type* sqrts=(this->s_hat_sampling)?NULL:(this->get_Ecm_hat());
		momentum_channel_type* c;
		while(flag=="<channel>" and !is.eof())
		{
		    c=factory_type::create_ps_channel(is,sqrts);
		    if(c==NULL)
		    {
			continue;
		    }
		    c->assign_mmin_components(this->eff_mmin);
		    c->set_multichannel_params(multichannel_params);
		    size_type i=0;
		    switch(c->bitstring.count())
		    {
			case N_out:
			    if(incoming_momentum_channel!=NULL)
			    {
				delete incoming_momentum_channel;
			    }
			    c->set_momentum(&(this->p_in(0)));
			    incoming_momentum_channel=c;
			    incoming_particle_channel=*(c->particle_channels.begin());
			    break;
			case 1:
			    while(!(c->bitstring[i]))
			    {
				++i;
			    }
			    c->set_momentum(&(this->p_out(i)));
			    if(outgoing_momentum_channels[i]!=NULL)
			    {
				delete outgoing_momentum_channels[i];
			    }
			    outgoing_momentum_channels[i]=c;
			    outgoing_particle_channels[i]=*(c->particle_channels.begin());
			    break;
			default:
			    momentum_channels.push_back(c);
		    }
		    std::getline(is,flag);
		}
		if(ps_branchings.size()!=0)
		{
		    for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		    {
			delete (*it);
		    }
		    ps_branchings.clear();
		}
		while(flag=="<branching>" and !is.eof())
		{
		    branching_type* br=factory_type::create_branching(is,incoming_particle_channel,momentum_channels,outgoing_particle_channels);
		    if(br!=NULL)
		    {
			br->spairs=&max_s_pairs;
			ps_branchings.push_back(br);
		    }
		    std::getline(is,flag);
		}
		incoming_particle_channel->sort_branchings();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    for(typename momentum_channel_type::particle_channel_list::iterator it2=(*it)->particle_channels.begin();it2!=(*it)->particle_channels.end();++it2)
		    {
			(*it2)->sort_branchings();
		    }
		}
		refresh_Ecm();
		refresh_m_min();
		return is;
	    }

	    /* Overridden saving method: */

	    std::ostream& save_data(std::ostream& os) const
	    {
		os<<adaptive_s<<"\t"<<adaptive_theta<<"\t"<<max_grid_bins<<std::endl;
		os<<max_s_pairs<<"\t";
		safe_write(os,multichannel_params.first);
		os<<"\t";
		safe_write(os,multichannel_params.second);
		os<<std::endl;
		for(typename momentum_channel_container::const_iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->save(os);
		}
		incoming_momentum_channel->save(os);
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->save(os);
		}
		for(typename branching_container::const_iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    (*it)->save(os);
		}
		return os;
	    }

	private:

	    /* Private modifiers */
	    /*-------------------*/

	    /* Adds a momentum channel to the list of no such with bit string bitstr
	     * exists, otherwise returns the respective list entry: */

	    momentum_channel_type* get_momentum_channel(bit_string_type bitstr)
	    {
		if(bitstr.count()==N_out)
		{
		    return incoming_momentum_channel;
		}
		if(bitstr.count()==1)
		{
		    return outgoing_momentum_channels[bitstr.to_integer()-1];
		}
		bitstring_eq pred(bitstr);
		typename momentum_channel_container::iterator it=std::find_if(momentum_channels.begin(),momentum_channels.end(),pred);
		if(it!=momentum_channels.end())
		{
		    return *it;
		}
		momentum_channel_type* c=new momentum_channel_type(bitstr);
		momentum_channel_comp<model_t,1,N_out,rng_t>comp;
		it=std::lower_bound(momentum_channels.begin(),momentum_channels.end(),c,comp);
		momentum_channels.insert(it,c);
		c->assign_mmin_components(this->eff_mmin);
		return c;
	    }

	    /* Adds a particle channel to the list of no such with bit string bitstr
	     * exists, otherwise returns the respective list entry: */

	    particle_channel_type* get_particle_channel(bit_string_type bitstr,const particle<model_t>* phi)
	    {
		if(bitstr.count()==N_out)
		{
		    return incoming_momentum_channel->find_particle_channel(phi);
		}
		if(bitstr.count()==1)
		{
		    return outgoing_momentum_channels[bitstr.to_integer()-1]->find_particle_channel(phi);
		}
		const value_type* sqrtshat=(this->s_hat_sampling)?(NULL):(this->get_Ecm_hat());
		particle_channel_type* result=get_momentum_channel(bitstr)->add_particle_channel(phi,sqrtshat);
		if(result!=NULL)
		{
		    result->set_multichannel_params(multichannel_params);
		}
		return result;
	    }

	    /* Returns the momentum channel with bitstring bitstr, returns NULL if not
	     * found. */

	    momentum_channel_type* find_momentum_channel(bit_string<N_out> bitstr)
	    {
		if(bitstr.count()==N_out)
		{
		    return incoming_momentum_channel;
		}
		if(bitstr.count()==1)
		{
		    return outgoing_momentum_channels[bitstr.to_integer()-1];
		}
		bitstring_eq pred(bitstr);
		typename momentum_channel_container::iterator it=std::find_if(momentum_channels.begin(),momentum_channels.end(),pred);
		if(it!=momentum_channels.end())
		{
		    return *it;
		}
		return NULL;
	    }

	    /* Returns the particle channel with bitstring bitstr, returns NULL if not
	     * found. */

	    particle_channel_type* find_particle_channel(bit_string<N_out> bitstr,const particle<model_t>* phi)
	    {
		if(bitstr.count()==N_out)
		{
		    return incoming_momentum_channel->find_particle_channel(phi);
		}
		if(bitstr.count()==1)
		{
		    return outgoing_momentum_channels[bitstr.to_integer()-1]->find_particle_channel(phi);
		}
		momentum_channel_type* pch=find_momentum_channel(bitstr);
		return (pch==NULL)?NULL:(pch->find_particle_channel(phi));
	    }

	    void set_redundant()
	    {
		incoming_particle_channel->set_redundant();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->set_redundant();
		}
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    (*it)->set_redundant();
		}
	    }

	    /* Cleans the phase space branching tree: */

	    void clean()
	    {
		set_redundant();
		incoming_particle_channel->unset_redundant();
		incoming_particle_channel->clean();
		
		for(typename momentum_channel_container::iterator it3=momentum_channels.begin();it3!=momentum_channels.end();++it3)
		{
		    (*it3)->clean();
		}
		typename std::list<branching_type*>::iterator it=std::stable_partition(ps_branchings.begin(),ps_branchings.end(),necessary_branching);
		for(typename std::list<branching_type*>::iterator it2=it;it2!=ps_branchings.end();++it2)
		{
		    delete *it2;
		}
		ps_branchings.erase(it,ps_branchings.end());
		typename std::list<momentum_channel_type*>::iterator it4=std::stable_partition(momentum_channels.begin(),momentum_channels.end(),non_empty_channel);
		for(typename std::list<momentum_channel_type*>::iterator it5=it4;it5!=momentum_channels.end();++it5)
		{
		    delete *it5;
		}
		momentum_channels.erase(it4,momentum_channels.end());
	    }

	    /* Inserts a two-fold branching in the generation tree: */

	    void insert_branching(bit_string_type b_in,const particle<model_t>* phi_in,
		    		  bit_string_type b_out1,const particle<model_t>* phi_out1,
			   	  bit_string_type b_out2,const particle<model_t>* phi_out2)
	    {
		if((b_in!=(b_out1 | b_out2)) or (b_out1 & b_out2).any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid channels "<<b_in<<"--->"<<b_out1<<','<<b_out2<<" encountered in branching construction--returning NULL"<<endlog;
		    return;
		}
		particle_channel_type* ch_in=get_particle_channel(b_in,phi_in);
		particle_channel_type* ch_out1=get_particle_channel(b_out1,phi_out1);
		particle_channel_type* ch_out2=get_particle_channel(b_out2,phi_out2);
		if(ch_in==NULL or ch_out1==NULL or ch_out2==NULL)
		{
		    return;
		}
		branching_type* br;
		if(backward_s_gen)
		{
		    br=factory_type::backward_s_branch(ch_in,ch_out1,ch_out2);
		}
		else
		{
		    br=factory_type::s_branch(ch_in,ch_out1,ch_out2);
		}
		if(br!=NULL)
		{
		    branching_type* brcpy=ch_in->insert_branching(br);
		    if(brcpy==br)
		    {
			insert_branching(brcpy);
		    }
		}
	    }

	    /* Inserts the branching in the list, provided it isn't already
	     * there: */

	    bool insert_branching(branching_type* br)
	    {
		if(std::find(ps_branchings.begin(),ps_branchings.end(),br)==ps_branchings.end())
		{
		    br->spairs=&max_s_pairs;
		    ps_branchings.push_back(br);
		    return true;
		}
		return false;
	    }

	    /* Private data members */
	    /*----------------------*/

	    /* Incoming momentum phase space channel: */

	    momentum_channel_type* incoming_momentum_channel;

	    /* Incoming particle channel: */

	    particle_channel_type* incoming_particle_channel;
	    
	    /* Internal momentum phase space channels: */
	    
	    momentum_channel_container momentum_channels;

	    /* Final-state phase space momentum channels: */

	    vector<momentum_channel_type*,N_out> outgoing_momentum_channels;

	    /* Final-state particle channels: */

	    vector<particle_channel_type*,N_out> outgoing_particle_channels;
	    
	    /* List of phase space branchings: */
	    
	    branching_container ps_branchings;

	    /* Currently chosen branchings: */

	    std::vector<branching_type*> branching_sequence;
    };

    /* Phase space decomposition tree generator for scattering processes: */

    template<class model_t,std::size_t N_out,class rng_t>class ps_tree<model_t,2,N_out,rng_t>: public ps_generator<model_t,2,N_out>
    {
	typedef ps_generator<model_t,2,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::init_state_type init_state_type;

	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    typedef momentum_channel<model_t,2,N_out,rng_t> momentum_channel_type;
	    typedef particle_channel<model_t,2,N_out,rng_t> particle_channel_type;
	    typedef ps_branching<model_t,2,N_out,rng_t> branching_type;
	    typedef typename std::list<momentum_channel_type*> momentum_channel_container;
	    typedef typename std::list<particle_channel_type*> particle_channel_container;
	    typedef typename std::list<branching_type*> branching_container;

	    typedef ps_factory<model_t,2,N_out,rng_t,typename model_t::spacetime_type> factory_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;

	    /* Static utility functions: */

	    static bool necessary_branching(const branching_type* c)
	    {
		return !(c->is_redundant());
	    }
	    static bool non_empty_channel(const momentum_channel_type* c)
	    {
		return c->any();
	    }
	    static const particle<model_t>* pmap(const bit_string_type b,const particle<model_t>* phi)
	    {
		if(b.count()==1 and b[0])
		{
		    return phi;
		}
		return (phi==NULL)?NULL:phi->get_anti_particle();
	    }

	    /* Utility classes: */

	    class bitstring_eq
	    {
		public:

		    const bit_string_type bitstring;

		    bitstring_eq(const bit_string_type& bitstring_):bitstring(bitstring_){}

		    bool operator()(const momentum_channel_type* p)
		    {
			return (p->bitstring==bitstring);
		    }
	    };
	    
	    /* Public static methods */
	    /*-----------------------*/

	    /* Factory method: */

	    static ps_tree<model_t,2,N_out,rng_t>* create_instance(typename CM_algorithm<model_t,2,N_out>::tree_iterator it,init_state_type* is_)
	    {
		typedef typename CM_algorithm<model_t,2,N_out>::tree_type CM_tree_type;
		typedef typename CM_tree_type::current_type current_type;
		typedef typename current_type::particle_type particle_type;
		typedef typename CM_tree_type::current_tree_type current_tree_type;

		if(current_tree_type::final_current()!=0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"phase space tree generator is only implemented for final current being zero--returning NULL"<<endlog;
		    return NULL;
		}
		ps_tree<model_t,2,N_out,rng_t>* result=new ps_tree<model_t,2,N_out,rng_t>(is_);
		if(!result->set_tree(it))
		{
		    return NULL;
		}
		const value_type* sqrts=(result->s_hat_sampling)?(NULL):(result->get_Ecm_hat());
		bit_string_type b;
		b.set(0);
		result->incoming_momentum_channels[1]=new momentum_channel_type(b,&(result->p_in(1)));
		result->incoming_particle_channels[1]=result->incoming_momentum_channels[1]->add_particle_channel(it->get_phase_space(1)->particle_type,sqrts);
		result->incoming_particle_channels[1]->set_multichannel_params(result->multichannel_params);
		for(size_type i=0;i<N_out;++i)
		{
		    b>>=1;
		    result->outgoing_momentum_channels[i]=new momentum_channel_type(b,&(result->p_out(i)));
		    result->outgoing_momentum_channels[i]->assign_mmin_components(result->eff_mmin);
		    result->outgoing_particle_channels[i]=result->outgoing_momentum_channels[i]->add_particle_channel(it->get_phase_space(i+2)->particle_type,sqrts);
		    result->outgoing_particle_channels[i]->set_multichannel_params(result->multichannel_params);
		}
		b.set();
		result->incoming_momentum_channels[0]=new momentum_channel_type(b,&(result->p_in(0)));
		result->incoming_particle_channels[0]=result->incoming_momentum_channels[0]->add_particle_channel(it->get_phase_space(0)->particle_type,sqrts);
		result->incoming_particle_channels[0]->set_multichannel_params(result->multichannel_params);

		for(typename CM_tree_type::const_interaction_iterator v_it=it->interactions_begin();v_it!=it->interactions_end();++v_it)
		{
		    if(v_it->is_coupled())
		    {
			bit_string<N_out+1>bs_in(v_it->get_produced_bit_string());
			std::vector< bit_string<N_out+1> >bs_out(v_it->get_incoming_bit_strings());
			const particle_type* phi_in=v_it->get_produced_particle();
			std::vector<const particle_type*>phi_out(v_it->get_incoming_particles());
			
			if(v_it->get_rank()==3)
			{
			    result->insert_branching(bs_in,pmap(bs_in,phi_in),bs_out[0],pmap(bs_out[0],phi_out[0]),bs_out[1],pmap(bs_out[1],phi_out[1]));
			}
			else if(v_it->get_rank()==4)
			{
			    result->insert_branching((bs_out[0] | bs_out[1]),NULL,bs_out[0],pmap(bs_out[0],phi_out[0]),bs_out[1],pmap(bs_out[1],phi_out[1]));
			    result->insert_branching(bs_in,pmap(bs_in,phi_in),(bs_out[0] | bs_out[1]),NULL,bs_out[2],pmap(bs_out[2],phi_out[2]));
			    result->insert_branching((bs_out[0] | bs_out[2]),NULL,bs_out[0],pmap(bs_out[0],phi_out[0]),bs_out[2],pmap(bs_out[2],phi_out[2]));
			    result->insert_branching(bs_in,pmap(bs_in,phi_in),(bs_out[0] | bs_out[2]),NULL,bs_out[1],pmap(bs_out[1],phi_out[1]));
			    result->insert_branching((bs_out[1] | bs_out[2]),NULL,bs_out[1],pmap(bs_out[1],phi_out[1]),bs_out[2],pmap(bs_out[2],phi_out[2]));
			    result->insert_branching(bs_in,pmap(bs_in,phi_in),(bs_out[1] | bs_out[2]),NULL,bs_out[0],pmap(bs_out[0],phi_out[0]));
			}
			else
			{
			    std::cout<<"Phase space branching for "<<v_it->get_rank()<<"-rank vertex not (yet) implemented..."<<std::endl;
			}
		    }
		}
		factory_type::flush_t_channels();
		result->clean();
		return result;
	    }

	    /* Public data members */
	    /*---------------------*/

	    /* Const flags. */
	    
	    const bool backward_s_gen;
	    const bool adaptive_s;
	    const bool adaptive_t;
	    const bool adaptive_theta;
	    const size_type max_grid_bins;

	    /* Maximal number of s-pair generations: */

	    size_type max_s_pairs;

	    /* Multichannel adaptivity and threshold: */

	    std::pair<value_type,value_type> multichannel_params;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Default constructor. */

	    ps_tree(init_state_type* is_):base_type(is_),backward_s_gen(backward_s_sampling()),adaptive_s(adaptive_s_sampling()),adaptive_t(adaptive_t_sampling()),adaptive_theta(adaptive_angles()),max_grid_bins(grid_bins()),max_s_pairs(Camgen::max_s_pairs()),multichannel_params(std::pair<value_type,value_type>(multichannel_adaptivity(),multichannel_threshold())),stot_channel(NULL)
	    {
		incoming_momentum_channels.assign(NULL);
		incoming_particle_channels.assign(NULL);
		outgoing_momentum_channels.assign(NULL);
		outgoing_particle_channels.assign(NULL);
	    }

	    /* Configurable constructor. */

	    ps_tree(init_state_type* is_,bool q1,bool q2,bool q3,bool q4,size_type b):base_type(is_),backward_s_gen(q1),adaptive_s(q2),adaptive_t(q3),adaptive_theta(q4),max_grid_bins(b),max_s_pairs(Camgen::max_s_pairs()),multichannel_params(std::pair<value_type,value_type>(multichannel_adaptivity(),multichannel_threshold())),stot_channel(NULL)
	    {
		incoming_momentum_channels.assign(NULL);
		incoming_particle_channels.assign(NULL);
		outgoing_momentum_channels.assign(NULL);
		outgoing_particle_channels.assign(NULL);
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /* Destructor. */

	    ~ps_tree()
	    {
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    delete *it;
		}
		if(incoming_momentum_channels[0]!=NULL)
		{
		    delete incoming_momentum_channels[0];
		}
		if(incoming_momentum_channels[1]!=NULL)
		{
		    delete incoming_momentum_channels[1];
		}
		for(size_type i=0;i<N_out;++i)
		{
		    if(outgoing_momentum_channels[i]!=NULL)
		    {
			delete outgoing_momentum_channels[i];
		    }
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    delete *it;
		}
	    }

	    /* Public modifiers */
	    /*------------------*/

	    /* Currently chosen branching sequence: */

	    std::vector<branching_type*> branching_sequence;

	    /* Topology generation method: */

	    bool choose_channel()
	    {
		branching_sequence.clear();
		branching_sequence.reserve(2*(2+N_out));
		incoming_particle_channels[0]->reset_generation_flags();
		bool q=incoming_particle_channels[0]->choose_branching(branching_sequence);
		return q;
	    }

	    /* Partonic CM-energy generation: */

	    bool generate_s_hat(value_type& s)
	    {
		if(backward_s_gen)
		{
		    if(stot_channel==NULL)
		    {
			s=(value_type)0;
			return false;
		    }
		    if(choose_channel())
		    {
			for(typename std::vector<branching_type*>::reverse_iterator it=branching_sequence.rbegin();it!=branching_sequence.rend();++it)
			{
			    if(!(*it)->generate_s())
			    {
				return false;
			    }
			}
			s=stot_channel->s();
			return true;
		    }
		    return false;
		}
		return true;
	    }

	    /* Final state generation method: */

	    bool generate_fs()
	    {
		incoming_particle_channels[0]->evaluate_s();
		incoming_particle_channels[1]->evaluate_s();

		if(this->s_hat_sampling)
		{
		    for(typename std::vector<branching_type*>::iterator it=branching_sequence.begin();it!=branching_sequence.end();++it)
		    {
			if(!((*it)->generate_t() and (*it)->generate_p()))
			{
			    this->fsw=(value_type)0;
			    return false;
			}
		    }
		    return evaluate_fs_weight();
		}
		else
		{
		    if(!choose_channel())
		    {
			this->fsw=(value_type)0;
			return false;
		    }
		    if(backward_s_gen)
		    {
			if(this->init_state()->hadronic)
			{
			    if(!refresh_Ecm_hat())
			    {
				this->fsw=(value_type)0;
				return false;
			    }
			}
			for(typename std::vector<branching_type*>::reverse_iterator it=branching_sequence.rbegin();it!=branching_sequence.rend();++it)
			{
			    if(!(*it)->generate_s())
			    {
				this->fsw=(value_type)0;
				return false;
			    }
			}
			for(typename std::vector<branching_type*>::iterator it=branching_sequence.begin();it!=branching_sequence.end();++it)
			{
			    if(!((*it)->generate_t() and (*it)->generate_p()))
			    {
				this->fsw=(value_type)0;
				return false;
			    }
			}
			return evaluate_fs_weight();
		    }
		    else
		    {
			if(!incoming_particle_channels[0]->generate())
			{
			    this->fsw=(value_type)0;
			    return false;
			}
			return evaluate_fs_weight();
		    }
		}
		this->fsw=(value_type)0;
		return false;
	    }

	    /* Weight evaluation method: */

	    bool evaluate_fs_weight()
	    {
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if(!(*it)->on_shell() and *it!=stot_channel)
		    {
			(*it)->p().assign(0);
			for(size_type i=1;i<N_out+1;++i)
			{
			    if((*it)->bitstring[i])
			    {
				(*it)->p()+=(outgoing_momentum_channels[i-1]->p());
			    }
			}
			if((*it)->bitstring[0])
			{
			    (*it)->p()-=incoming_momentum_channels[1]->p();
			}
			(*it)->evaluate_s();
		    }
		}
		bool q=true;
		for(size_type i=0;i<N_out;++i)
		{
		    q&=(outgoing_particle_channels[i]->evaluate_weight());
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if(!(*it)->on_shell())
		    {
			q&=((*it)->evaluate_weight());
		    }
		}
		q&=(incoming_particle_channels[0]->evaluate_weight());
		this->fsw=q?(incoming_particle_channels[0]->weight()):(value_type)0;
		return q;
	    }

	    /* Refreshes the minimal invariant masses: */

	    bool refresh_m_min()
	    {
		bool q=this->base_type::refresh_m_min();
		value_type sa=(this->m_in(0))*(this->m_in(0));
		value_type sb=(this->m_in(1))*(this->m_in(1));
		q&=incoming_momentum_channels[1]->refresh_m_min(this->s_tot(),sa,sb);
		for(size_type i=0;i<N_out;++i)
		{
		    q&=(outgoing_momentum_channels[i]->refresh_m_min(this->s_tot(),sa,sb));
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    q&=((*it)->refresh_m_min(this->s_tot(),sa,sb));
		}
		q&=incoming_momentum_channels[0]->refresh_m_min(this->s_tot(),sa,sb);
		if(q)
		{
		    this->init_state()->set_m_hat_min(incoming_momentum_channels[0]->m_min_min());
		}
		return q;
	    }

	    /* Refreshes the hadronic invariant mass: */

	    bool refresh_Ecm()
	    {
		if(!this->init_state()->refresh_Ecm())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"initial state does not allow proposed beam energies"<<endlog;
		    return false;
		}
		bool q=true;
		if(this->s_hat_sampling)
		{
		    for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		    {
			if(!(*it)->bitstring[0])
			{
			    q&=((*it)->set_m_max_max(this->Ecm()));
			    q&=((*it)->set_m_max(this->Ecm()));
			}
			else
			{
			    q&=((*it)->set_m_min_min(-this->Ecm()));
			    q&=((*it)->set_m_min(-this->Ecm()));
			}
		    }
		}
		else
		{
		    for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		    {
			if(!(*it)->bitstring[0])
			{
			    q&=((*it)->set_m_max_max(this->Ecm()));
			}
			else
			{
			    q&=((*it)->set_m_min_min(-this->Ecm()));
			}
		    }
		}
		if(!this->init_state()->hadronic)
		{
		    q&=refresh_Ecm_hat();
		}
		return q;
	    }

	    /* Refreshes the hadronic invariant mass: */

	    bool refresh_Ecm_hat()
	    {
		if(stot_channel!=NULL)
		{
		    stot_channel->refresh_params();
		}
		bool q=true;
		if(backward_s_gen and !this->s_hat_sampling)
		{
		    for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		    {
			if(!(*it)->bitstring[0] and *it!=stot_channel)
			{
			    q&=((*it)->set_m_max(this->Ecm_hat()));
			}
		    }
		}
		return q;
	    }

	    /* Refreshes the internal parameters: */

	    bool refresh_params()
	    {
		bool q=this->base_type::refresh_params();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    q&=((*it)->refresh_params());
		}
		return q;
	    }

	    /* Sets the minimal invariant mass i,j and refreshes depending
	     * channels: */

	    bool set_m_min(int i,int j,const value_type* sqrts)
	    {
		if(this->base_type::set_m_min(i,j,sqrts))
		{
		    bool q1=true,q2=false;
		    value_type mhatmin(0);
		    value_type sa=(this->m_in(0))*(this->m_in(0));
		    value_type sb=(this->m_in(1))*(this->m_in(1));
		    for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		    {
			if((*it)->bitstring[i] and (*it)->bitstring[j])
			{
			    q1&=((*it)->refresh_m_min(this->s_tot(),sa,sb));
			    if(!(*it)->bitstring[0] and (*it)->bitstring.count()==N_out)
			    {
				mhatmin=(*it)->m_min_min();
				q2=true;
			    }
			}
		    }
		    if(q2)
		    {
			this->init_state()->set_m_hat_min(mhatmin);
		    }
		    else
		    {
			q1&=(this->init_state()->refresh_m_min());
		    }
		    return q1;
		}
		return false;
	    }

	    /* Overrides the update method: */

	    void update_fs()
	    {
		value_type f=this->integrand()*this->weight();
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    if((*it)->alpha!=(value_type)0)
		    {
			(*it)->integrand()=f;
			(*it)->update();
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->update(f);
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->update(f);
		}
		incoming_momentum_channels[0]->update(f);
		incoming_momentum_channels[1]->update(f);
	    }

	    /* Overrides the vegas-grid adaptation method: */

	    void adapt_grids()
	    {
		this->base_type::adapt_grids();
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    if((*it)->alpha!=(value_type)0)
		    {
			(*it)->adapt_grids();
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->adapt_grids();
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->adapt_grids();
		}
		incoming_momentum_channels[0]->adapt_grids();
		incoming_momentum_channels[1]->adapt_grids();
	    }
	    
	    /* Overrides the multichannel adaptation method: */

	    void adapt_channels()
	    {
		this->base_type::adapt_channels();
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    if((*it)->alpha!=(value_type)0)
		    {
			(*it)->adapt_channels();
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->adapt_channels();
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->adapt_channels();
		}
		incoming_momentum_channels[0]->adapt_channels();
		incoming_momentum_channels[1]->adapt_channels();
		clean();
	    }

	    /* Resets cross-sections, multichannel weights and adaptive grids:
	     * */

	    void reset()
	    {
		this->base_type::reset();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->reset();
		}
		incoming_momentum_channels[0]->reset();
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Type identifier: */

	    std::string type() const
	    {
		if(backward_s_gen)
		{
		    return "pstree~";
		}
		else
		{
		    return "pstree";
		}
	    }

	    /* Printing method. */
	    
	    std::ostream& print(std::ostream& os) const
	    {
		incoming_momentum_channels[0]->print(os);
		incoming_momentum_channels[1]->print(os);
		for(typename momentum_channel_container::const_reverse_iterator it=momentum_channels.rbegin();it!=momentum_channels.rend();++it)
		{
		    (*it)->print(os);
		}
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->print(os);
		}
		return os;
	    }

	    /// Virtual method printing the fs-generator settings.

	    std::ostream& print_fs_settings(std::ostream& os) const
	    {
		os<<std::setw(30)<<std::left<<"s sampling dir:"<<(backward_s_gen?"backward":"forward")<<std::endl;
		os<<std::setw(30)<<std::left<<"adaptive s-sampling:"<<(adaptive_s?"yes":"no")<<std::endl;
		os<<std::setw(30)<<std::left<<"adaptive t-sampling:"<<(adaptive_t?"yes":"no")<<std::endl;
		os<<std::setw(30)<<std::left<<"adaptive angles:"<<(adaptive_theta?"yes":"no")<<std::endl;
		if(adaptive_s or adaptive_t or adaptive_theta)
		{
		    os<<std::setw(30)<<std::left<<"max nr bins/grid:"<<max_grid_bins<<std::endl;
		}
		os<<std::setw(30)<<std::left<<"channel adaptivity:"<<multichannel_adaptivity()<<std::endl;
		os<<std::setw(30)<<std::left<<"channel threshold:"<<multichannel_threshold()<<std::endl;
		return os;
	    }

	    /* Generation channel printing method. */

	    std::ostream& print_channel(std::ostream& os) const
	    {
		if(branching_sequence.size()>0)
		{
		    for(size_type i=0;i<branching_sequence.size();++i)
		    {
			branching_sequence[i]->print(os);
			os<<std::endl;
		    }
		}
		else
		{
		    os<<"No branching sequence generated."<<std::endl;
		}
		return os;
	    }

	    /* Overridden loading method: */

	    std::istream& load_data(std::istream& is)
	    {
		is>>max_s_pairs;
		safe_read(is,multichannel_params.first);
		safe_read(is,multichannel_params.second);
		std::string flag;
		while(flag!="<channel>" and flag!="</fsgen>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(flag!="<channel>")
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"zero momentum channels read from file"<<endlog;
		    return is;
		}
		if(momentum_channels.size()!=0)
		{
		    for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		    {
			delete (*it);
		    }
		    momentum_channels.clear();
		}
		const value_type* sqrts=(this->s_hat_sampling)?NULL:(this->get_Ecm_hat());
		momentum_channel_type* c;
		while(flag=="<channel>" and !is.eof())
		{
		    c=factory_type::create_ps_channel(is,sqrts);
		    if(c==NULL)
		    {
			continue;
		    }
		    c->assign_mmin_components(this->eff_mmin);
		    c->set_multichannel_params(multichannel_params);
		    switch(c->bitstring.count())
		    {
			case (N_out+1):
			    if(incoming_momentum_channels[0]!=NULL)
			    {
				delete incoming_momentum_channels[0];
			    }
			    c->set_momentum(&(this->p_in(0)));
			    incoming_momentum_channels[0]=c;
			    incoming_particle_channels[0]=*(c->particle_channels.begin());
			    break;
			case 1:
			    if(c->bitstring[0])
			    {
				if(incoming_momentum_channels[1]!=NULL)
				{
				    delete incoming_momentum_channels[1];
				}
				c->set_momentum(&(this->p_in(1)));
				incoming_momentum_channels[1]=c;
				incoming_particle_channels[1]=*(c->particle_channels.begin());
			    }
			    else
			    {
				size_type i=0;
				while(!(c->bitstring[i+1]))
				{
				    ++i;
				}
				c->set_momentum(&(this->p_out(i)));
				if(outgoing_momentum_channels[i]!=NULL)
				{
				    delete outgoing_momentum_channels[i];
				}
				outgoing_momentum_channels[i]=c;
				outgoing_particle_channels[i]=*(c->particle_channels.begin());
			    }
			    break;
			default:
			    if(c->timelike() and c->bitstring.count()==N_out)
			    {
				stot_channel=c;
			    }
			    momentum_channels.push_back(c);
		    }
		    std::getline(is,flag);
		}
		if(ps_branchings.size()!=0)
		{
		    for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		    {
			delete (*it);
		    }
		    ps_branchings.clear();
		}
		while(flag=="<branching>" and !is.eof())
		{
		    branching_type* br=factory_type::create_branching(is,incoming_particle_channels,momentum_channels,outgoing_particle_channels);
		    if(br!=NULL)
		    {
			br->spairs=&max_s_pairs;
			ps_branchings.push_back(br);
		    }
		    std::getline(is,flag);
		}
		clean();
		incoming_particle_channels[0]->sort_branchings();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    for(typename momentum_channel_type::particle_channel_list::iterator it2=(*it)->particle_channels.begin();it2!=(*it)->particle_channels.end();++it2)
		    {
			(*it2)->sort_branchings();
		    }
		}
		return is;
	    }

	    /* Overridden saving method: */

	    std::ostream& save_data(std::ostream& os) const
	    {
		os<<max_s_pairs<<"\t";
		safe_write(os,multichannel_params.first);
		os<<"\t";
		safe_write(os,multichannel_params.second);
		os<<std::endl;
		incoming_momentum_channels[1]->save(os);
		for(typename momentum_channel_container::const_iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->save(os);
		}
		incoming_momentum_channels[0]->save(os);
		for(size_type i=0;i<N_out;++i)
		{
		    outgoing_momentum_channels[i]->save(os);
		}
		for(typename branching_container::const_iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    (*it)->save(os);
		}
		return os;
	    }
	    
	    /* Outputs the const data: */

	    std::ostream& save_preamble(std::ostream& os) const
	    {
		os<<adaptive_s<<"\t"<<adaptive_t<<"\t"<<adaptive_theta<<"\t"<<max_grid_bins<<std::endl;
		return os;
	    }

	private:

	    /* Private modifiers */
	    /*-------------------*/

	    /* Adds a phase space channel to the list of no such with bit string bitstr
	     * exists, otherwise returns the respective list entry: */

	    momentum_channel_type* get_momentum_channel(bit_string_type bitstr)
	    {
		if(bitstr.count()==N_out+1)
		{
		    return incoming_momentum_channels[0];
		}
		if(bitstr.count()==1)
		{
		    if(bitstr[0])
		    {
			return incoming_momentum_channels[1];
		    }
		    return outgoing_momentum_channels[bitstr.to_integer()-2];
		}
		bitstring_eq pred(bitstr);
		typename momentum_channel_container::iterator it=std::find_if(momentum_channels.begin(),momentum_channels.end(),pred);
		if(it!=momentum_channels.end())
		{
		    return *it;
		}
		momentum_channel_type* c=new momentum_channel_type(bitstr);
		if(!bitstr[0] and bitstr.count()==N_out)
		{
		    stot_channel=c;
		}
		momentum_channel_comp<model_t,2,N_out,rng_t>comp;
		it=std::lower_bound(momentum_channels.begin(),momentum_channels.end(),c,comp);
		momentum_channels.insert(it,c);
		c->assign_mmin_components(this->eff_mmin);
		return c;
	    }

	    /* Adds a phase space momentum channel and particle channel with
	     * momentum bit string bitstr and flavour phi. Returns the instance
	     * if such object already exist: */

	    particle_channel_type* get_particle_channel(bit_string_type bitstr,const particle<model_t>* phi)
	    {
		if(bitstr.count()==N_out+1)
		{
		    return incoming_momentum_channels[0]->find_particle_channel(phi);
		}
		if(bitstr.count()==1)
		{
		    if(bitstr[0])
		    {
			return incoming_momentum_channels[1]->find_particle_channel(phi);
		    }
		    return outgoing_momentum_channels[bitstr.to_integer()-2]->find_particle_channel(phi);
		}
		const value_type* sqrtshat=(this->s_hat_sampling)?NULL:(this->get_Ecm_hat());
		particle_channel_type* result=get_momentum_channel(bitstr)->add_particle_channel(phi,sqrtshat);
		if(result!=NULL)
		{
		    result->set_multichannel_params(multichannel_params);
		}
		return result;
	    }

	    /* Returns the channel with bitstring bitstr, returns NULL if not
	     * found. */

	    momentum_channel_type* find_momentum_channel(bit_string_type bitstr)
	    {
		if(bitstr.count()==N_out+1)
		{
		    return incoming_momentum_channels[0];
		}
		if(bitstr.count()==1)
		{
		    if(bitstr[0])
		    {
			return incoming_momentum_channels[1];
		    }
		    return outgoing_momentum_channels[bitstr.to_integer()-2];
		}
		bitstring_eq pred(bitstr);
		typename momentum_channel_container::iterator it=std::find_if(momentum_channels.begin(),momentum_channels.end(),pred);
		if(it!=momentum_channels.end())
		{
		    return *it;
		}
		return NULL;
	    }

	    /* Finds a phase space particle channel with momentum bit string
	     * bitstr and flavour phi. Returns the null pointer if no such
	     * object was found. */

	    particle_channel_type* find_particle_channel(bit_string_type bitstr,const particle<model_t>* phi)
	    {
		if(bitstr.count()==N_out+1)
		{
		    return incoming_momentum_channels[0]->find_particle_channel(phi);
		}
		if(bitstr.count()==1)
		{
		    if(bitstr[0])
		    {
			return incoming_momentum_channels[1]->find_particle_channel(phi);
		    }
		    return outgoing_momentum_channels[bitstr.to_integer()-2]->find_particle_channel(phi);
		}
		momentum_channel_type* pch=find_momentum_channel(bitstr);
		return (pch==NULL)?NULL:(pch->find_particle_channel(phi));
	    }

	    void set_redundant()
	    {
		incoming_particle_channels[0]->set_redundant();
		incoming_particle_channels[1]->set_redundant();
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->set_redundant();
		}
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    (*it)->set_redundant();
		}
	    }

	    /* Cleans the phase space branching tree: */

	    void clean()
	    {
		set_redundant();
		incoming_particle_channels[0]->unset_redundant();
		incoming_particle_channels[0]->clean();

		for(typename momentum_channel_container::iterator it3=momentum_channels.begin();it3!=momentum_channels.end();++it3)
		{
		    (*it3)->clean();
		}
		typename std::list<branching_type*>::iterator it=std::stable_partition(ps_branchings.begin(),ps_branchings.end(),necessary_branching);
		for(typename std::list<branching_type*>::iterator it2=it;it2!=ps_branchings.end();++it2)
		{
		    delete *it2;
		}
		ps_branchings.erase(it,ps_branchings.end());
		if(stot_channel!=NULL)
		{
		    if(stot_channel->empty())
		    {
			stot_channel=NULL;
		    }
		}
		typename std::list<momentum_channel_type*>::iterator it4=std::stable_partition(momentum_channels.begin(),momentum_channels.end(),non_empty_channel);
		for(typename std::list<momentum_channel_type*>::iterator it5=it4;it5!=momentum_channels.end();++it5)
		{
		    delete *it5;
		}
		momentum_channels.erase(it4,momentum_channels.end());
	    }

	    /* Inserts a two-fold branching in the generation tree: */

	    void insert_branching(bit_string_type b_in,const particle<model_t>* phi_in,
		    		  bit_string_type b_out1,const particle<model_t>* phi_out1,
				  bit_string_type b_out2,const particle<model_t>* phi_out2)
	    {
		if((b_in!=(b_out1 | b_out2)) or (b_out1 & b_out2).any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid channels "<<b_in<<"--->"<<b_out1<<','<<b_out2<<" encountered in branching construction--returning NULL"<<endlog;
		    return;
		}
		particle_channel_type* ch_in=get_particle_channel(b_in,phi_in);
		particle_channel_type* ch_out1=get_particle_channel(b_out1,phi_out1);
		particle_channel_type* ch_out2=get_particle_channel(b_out2,phi_out2);
		if(ch_in==NULL or ch_out1==NULL or ch_out2==NULL)
		{
		    return;
		}
		if(b_in[0])
		{
		    if(b_out1[0])
		    {
			std::swap(ch_out1,ch_out2);
		    }
		    if(ch_out2->bitstring().count()==1)
		    {
			if(ch_out1->bitstring().count()==N_out)
			{
			    branching_type* br;
			    if(backward_s_gen)
			    {
				br=factory_type::backward_sum_branch(ch_in,ch_out2,ch_out1);
			    }
			    else
			    {
				br=factory_type::sum_branch(ch_in,ch_out2,ch_out1);
			    }
			    if(br!=NULL)
			    {
				branching_type* brcpy=ch_in->insert_branching(br);
				if(brcpy==br)
				{
				    insert_branching(brcpy);
				}
			    }
			}
			else
			{
			    factory_type::register_t_channel(ch_in,ch_out1,true);
			}
		    }
		    else
		    {
			std::vector<branching_type*>brvec;
			if(backward_s_gen)
			{
			    brvec=factory_type::backward_t_branch(ch_in,incoming_particle_channels[1],ch_out1,ch_out2);
			}
			else
			{
			    brvec=factory_type::t_branch(ch_in,incoming_particle_channels[1],ch_out1,ch_out2);
			}
			for(size_type i=0;i<brvec.size();++i)
			{
			    if(backward_s_gen)
			    {
				bit_string_type s12=(b_out1 | b_out2);
				s12.reset(0);
				static_cast<backward_t_branching<model_t,N_out,rng_t,spacetime_type>*>(brvec[i])->set_s12_channel(get_particle_channel(s12,NULL));
			    }
			    branching_type* brcpy=ch_in->insert_branching(brvec[i]);
			    if(brcpy==brvec[i])
			    {
				insert_branching(brcpy);
			    }
			}
			bit_string_type bstot=(b_out1 | b_out2);
			bstot[0]=false;
			factory_type::register_t_channel(ch_in,get_particle_channel(bstot,NULL),false);
		    }
		}
		else
		{
		    branching_type* br;
		    if(backward_s_gen)
		    {
			br=factory_type::backward_s_branch(ch_in,ch_out1,ch_out2);
		    }
		    else
		    {
			br=factory_type::s_branch(ch_in,ch_out1,ch_out2);
		    }
		    if(br!=NULL)
		    {
			branching_type* brcpy=ch_in->insert_branching(br);
			if(brcpy==br)
			{
			    insert_branching(brcpy);
			}
		    }
		}
	    }

	    /* Inserts the branching in the list, provided it isn't already
	     * there: */

	    bool insert_branching(branching_type* br)
	    {
		if(std::find(ps_branchings.begin(),ps_branchings.end(),br)==ps_branchings.end())
		{
		    br->spairs=&max_s_pairs;
		    ps_branchings.push_back(br);
		    return true;
		}
		return false;
	    }

	    /* Private data members */
	    /*----------------------*/

	    /* Incoming phase space momentum channels: */

	    vector<momentum_channel_type*,2> incoming_momentum_channels;

	    /* Incoming phase space particle channels: */

	    vector<particle_channel_type*,2> incoming_particle_channels;
	    
	    /* Internal phase space momentum channels: */
	    
	    momentum_channel_container momentum_channels;

	    /* Final-state phase space momentum channels: */

	    vector<momentum_channel_type*,N_out> outgoing_momentum_channels;

	    /* Total outgoing momentum channel: */

	    momentum_channel_type* stot_channel;

	    /* Final-state phase space momentum channels: */

	    vector<particle_channel_type*,N_out> outgoing_particle_channels;
	    
	    /* List of phase space branchings: */
	    
	    branching_container ps_branchings;
    };
}

#endif /*CAMGEN_PS_TREE_H_*/

