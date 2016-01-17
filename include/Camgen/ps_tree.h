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
#include <Camgen/ps_branching_fac.h>
#include <Camgen/ps_tree_base.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_tree;

    /// Phase space decomposition tree generator for decay processes.

    template<class model_t,std::size_t N_out,class rng_t>class ps_tree<model_t,1,N_out,rng_t>: public ps_tree_base<model_t,1,N_out,rng_t>
    {
	typedef ps_tree_base<model_t,1,N_out,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::model_type model_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::init_state_type init_state_type;
	    typedef typename base_type::rn_engine rn_engine;
	    typedef typename base_type::rn_stream rn_stream;
	    typedef typename base_type::momentum_channel_type momentum_channel_type;
	    typedef typename base_type::particle_channel_type particle_channel_type;
	    typedef typename base_type::branching_type branching_type;
	    typedef typename base_type::momentum_channel_container momentum_channel_container;
	    typedef typename base_type::particle_channel_container particle_channel_container;
	    typedef typename base_type::branching_container branching_container;
	    typedef typename base_type::branching_factory_type branching_factory_type;
	    typedef typename base_type::bit_string_type bit_string_type;
	    
	    /* Public constructors: */
	    /*----------------------*/

	    /// Constructor without recursive amplitude.

	    ps_tree(init_state_type* is_,bool backward_s_gen_=false):base_type(is_,backward_s_gen_){}

	    /* Public modifiers */
	    /*------------------*/

	    /// Overriden amplitude insertion method.

	    bool set_amplitude(typename CM_algorithm<model_t,1,N_out>::tree_iterator it)
	    {
		typedef CM_algorithm<model_t,1,N_out> CM_algorithm_type;
		typedef typename CM_algorithm_type::tree_type current_tree_type;
		typedef typename current_tree_type::const_interaction_iterator const_interaction_iterator;
		typedef typename current_tree_type::particle_type particle_type;

		this->delete_tree();
		if(!this->base_type::set_amplitude(it))
		{
		    return false;
		}

		this->set_external_particles(it);
		
		for(const_interaction_iterator v_it=it->interactions_begin();v_it!=it->interactions_end();++v_it)
		{
		    bit_string_type bs_in(v_it->get_produced_bit_string());
		    std::vector< bit_string_type >bs_out(v_it->get_incoming_bit_strings());
		    const particle_type* phi_in=v_it->get_produced_particle();
		    std::vector<const particle_type*>phi_out(v_it->get_incoming_particles());
		    
		    if(v_it->get_rank()==3)
		    {
			const particle_type* phi_bar_out[2]={phi_out[0]->get_anti_particle(),phi_out[1]->get_anti_particle()};
			
			this->insert_branching(bs_in,phi_in->get_anti_particle(),bs_out[0],phi_bar_out[0],bs_out[1],phi_bar_out[1]);
		    }
		    else if(v_it->get_rank()==4)
		    {
			const particle_type* phi_bar_out[3]={phi_out[0]->get_anti_particle(),
			    				     phi_out[1]->get_anti_particle(),
							     phi_out[2]->get_anti_particle()};

			particle_type* rho=this->create_auxiliary_particle();

			this->insert_branching((bs_out[0] | bs_out[1]),rho,bs_out[0],phi_bar_out[0],bs_out[1],phi_bar_out[1]);
			this->insert_branching(bs_in,phi_in->get_anti_particle(),(bs_out[0] | bs_out[1]),rho,bs_out[2],phi_bar_out[2]);
			this->insert_branching((bs_out[0] | bs_out[2]),rho,bs_out[0],phi_bar_out[0],bs_out[2],phi_bar_out[2]);
			this->insert_branching(bs_in,phi_in->get_anti_particle(),(bs_out[0] | bs_out[2]),rho,bs_out[1],phi_bar_out[1]);
			this->insert_branching((bs_out[1] | bs_out[2]),rho,bs_out[1],phi_bar_out[1],bs_out[2],phi_bar_out[2]);
			this->insert_branching(bs_in,phi_in->get_anti_particle(),(bs_out[1] | bs_out[2]),rho,bs_out[0],phi_bar_out[0]);
		    }
		    else
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"Phase space branching for "<<v_it->get_rank()<<"-rank vertex not (yet) implemented, vertex ignored."<<endlog;
		    }
		}
		if(this->empty())
		{
		    std::stringstream ss;
		    it->print_process(ss);
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"The phase space generator for "<<ss.str()<<" contains no branchings."<<endlog;
		}
		return true;
	    }

	    bool generate()
	    {
		this->initialise_channels();

		for(typename branching_container::iterator it=this->ps_branchings.begin();it!=this->ps_branchings.end();++it)
		{
		    (*it)->weight()=(value_type)0;
		}

		branching_container selected_branchings=this->select_branchings();

		if(selected_branchings.size()==0)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!this->generate_is())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!this->check_sufficient_shat())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->incoming_particle_channels[0]->set_status_p_generated();
		if(!this->refresh_Ecm_hat())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(this->backward_s_sampling())
		{
		    if(!generate_s_back(selected_branchings))
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		}
		else
		{
		    if(!base_type::generate_s(selected_branchings))
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		}
		if(!base_type::generate_p(selected_branchings))
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		base_type::set_generated(selected_branchings);
		bool q=this->evaluate_weight();
		base_type::reset_generated(selected_branchings);
		return q;
	    }

	private:


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
		particle_channel_type* ch_in=this->find_or_add_particle_channel(b_in,phi_in);
		particle_channel_type* ch_out1=this->find_or_add_particle_channel(b_out1,phi_out1);
		particle_channel_type* ch_out2=this->find_or_add_particle_channel(b_out2,phi_out2);
		if(ch_in==NULL or ch_out1==NULL or ch_out2==NULL)
		{
		    return;
		}
		branching_type* s_branching=branching_factory_type::create_s_branching(ch_in,ch_out1,ch_out2);
		this->try_insert_branching_or_delete(ch_in,s_branching);
	    }

	    /* Evaluates internal momentum: */

	    momentum_type p_internal(bit_string_type bs)
	    {
		momentum_type result;
		result.assign((value_type)0);
		for(size_type i=0;i<N_out;++i)
		{
		    if(!bs[i]) continue;
		    result+=(this->outgoing_particle_channels[i]->p());
		}
		return result;
	    }

	    /* Adds a momentum channel to the list if no such with bit string bs
	     * exists, otherwise returns the respective channel: */

	    momentum_channel_type* find_or_add_momentum_channel(bit_string_type bs)
	    {
		bitstring_compare<model_t,1,N_out,rng_t> pred(bs);
		typename momentum_channel_container::iterator it=std::find_if(this->momentum_channels.begin(),this->momentum_channels.end(),pred);
		if(it!=this->momentum_channels.end() and (*it)->bitstring==bs)
		{
		    return *it;
		}
		momentum_channel_type* c=NULL;
		if(bs.count()==1)
		{
		    size_type i=0;
		    while(!bs[i])
		    {
			++i;
		    }
		    c=new momentum_channel_type(bs,&(this->p_out(i)));
		}
		else if(bs.count()==N_out)
		{
		    c=new momentum_channel_type(bs,&(this->p_in(0)));
		}
		else
		{
		    c=new momentum_channel_type(bs);
		}
		c->assign_mmin_components(this->eff_mmin);
		this->momentum_channels.insert(it,c);
		return c;
	    }

	    /* Generates positive invariants for the given sequence of branchings backward: */

	    static bool generate_s_back(branching_container& selected_branchings)
	    {
		for(typename branching_container::reverse_iterator it=selected_branchings.rbegin();it!=selected_branchings.rend();++it)
		{
		    if(!(*it)->generate_s())
		    {
			return false;
		    }
		}
		return true;
	    }
    };

    /* Phase space decomposition tree generator for scattering processes: */

    template<class model_t,std::size_t N_out,class rng_t>class ps_tree<model_t,2,N_out,rng_t>: public ps_tree_base<model_t,2,N_out,rng_t>
    {
	typedef ps_tree_base<model_t,2,N_out,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::model_type model_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::init_state_type init_state_type;
	    typedef typename base_type::rn_engine rn_engine;
	    typedef typename base_type::rn_stream rn_stream;
	    typedef typename base_type::momentum_channel_type momentum_channel_type;
	    typedef typename base_type::particle_channel_type particle_channel_type;
	    typedef typename base_type::branching_type branching_type;
	    typedef typename base_type::momentum_channel_container momentum_channel_container;
	    typedef typename base_type::particle_channel_container particle_channel_container;
	    typedef typename base_type::branching_container branching_container;
	    typedef typename base_type::branching_factory_type branching_factory_type;
	    typedef typename base_type::bit_string_type bit_string_type;

	    /* Static utility functions: */

	    static const particle<model_t>* get_anti_particle(const bit_string_type b,const particle<model_t>* phi)
	    {
		if(b.count()==1 and b[0])
		{
		    return phi;
		}
		return (phi==NULL)?NULL:phi->get_anti_particle();
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /// Constructor without recursive amplitude.

	    ps_tree(init_state_type* is_,bool backward_s_gen_=false):base_type(is_,backward_s_gen_){}

	    /* Public modifiers */
	    /*------------------*/

	    /// Overriden amplitude insertion method.

	    bool set_amplitude(typename CM_algorithm<model_t,2,N_out>::tree_iterator it)
	    {
		typedef CM_algorithm<model_t,2,N_out> CM_algorithm_type;
		typedef typename CM_algorithm_type::tree_type current_tree_type;
		typedef typename current_tree_type::const_interaction_iterator const_interaction_iterator;
		typedef typename current_tree_type::particle_type particle_type;

		this->delete_tree();
		if(!this->base_type::set_amplitude(it))
		{
		    return false;
		}
		
		this->set_external_particles(it);

		for(const_interaction_iterator v_it=it->interactions_begin();v_it!=it->interactions_end();++v_it)
		{
		    if(!v_it->is_coupled())
		    {
			continue;
		    }

		    bit_string<N_out+1>bs_in(v_it->get_produced_bit_string());
		    std::vector< bit_string<N_out+1> >bs_out(v_it->get_incoming_bit_strings());
		    const particle_type* phi_in=v_it->get_produced_particle();
		    const particle_type* phi_in_bar=get_anti_particle(bs_in,phi_in);
		    std::vector<const particle_type*>phi_out(v_it->get_incoming_particles());

		    if(v_it->get_rank()==3)
		    {
			const particle_type* phi_out_bar[2]={get_anti_particle(bs_out[0],phi_out[0]),
			    				     get_anti_particle(bs_out[1],phi_out[1])};

			this->insert_branching(bs_in,phi_in_bar,bs_out[0],phi_out_bar[0],bs_out[1],phi_out_bar[1]);
		    }
		    else if(v_it->get_rank()==4)
		    {
			const particle_type* phi_out_bar[3]={get_anti_particle(bs_out[0],phi_out[0]),
			    				     get_anti_particle(bs_out[1],phi_out[1]),
							     get_anti_particle(bs_out[2],phi_out[2])};

			particle_type* rho=this->create_auxiliary_particle();

			this->insert_branching((bs_out[0] | bs_out[1]),rho,bs_out[0],phi_out_bar[0],bs_out[1],phi_out_bar[1]);
			this->insert_branching(bs_in,phi_in_bar,(bs_out[0] | bs_out[1]),rho,bs_out[2],phi_out_bar[2]);
			this->insert_branching((bs_out[0] | bs_out[2]),rho,bs_out[0],phi_out_bar[0],bs_out[2],phi_out_bar[2]);
			this->insert_branching(bs_in,phi_in_bar,(bs_out[0] | bs_out[2]),rho,bs_out[1],phi_out_bar[1]);
			this->insert_branching((bs_out[1] | bs_out[2]),rho,bs_out[1],phi_out_bar[1],bs_out[2],phi_out_bar[2]);
			this->insert_branching(bs_in,phi_in_bar,(bs_out[1] | bs_out[2]),rho,bs_out[0],phi_out_bar[0]);
		    }
		    else
		    {
			std::cout<<"Phase space branching for "<<v_it->get_rank()<<"-rank vertex not (yet) implemented..."<<std::endl;
		    }
		}

		replace_obsolete_sum_branchings();

		this->clean_tree();

		return true;
	    }

	    /// Initial state generation method.

	    bool generate_is()
	    {
		if(!this->base_type::generate_is())
		{
		    return false;
		}
		this->incoming_particle_channels[0]->evaluate_s();
		this->incoming_particle_channels[0]->set_status_p_generated();
		this->incoming_particle_channels[1]->evaluate_s();
		this->incoming_particle_channels[1]->set_status_p_generated();
		shat_channel->p()=this->incoming_particle_channels[0]->p()+this->incoming_particle_channels[1]->p();
		shat_channel->evaluate_s();
		shat_channel->set_status_p_generated();
		return true;
	    }

	    /// Generation method.

	    bool generate()
	    {
		this->initialise_channels();

		for(typename branching_container::iterator it=this->ps_branchings.begin();it!=this->ps_branchings.end();++it)
		{
		    (*it)->weight()=(value_type)0;
		}

		branching_container selected_branchings=this->select_branchings();

		if(selected_branchings.size()==0)
		{
		    this->weight()=(value_type)0;
		    return false;
		}

		if(this->backward_shat_sampling())
		{
		    if(!this->check_sufficient_s())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    if(!generate_s_back(selected_branchings))
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    if(!base_type::generate_t(selected_branchings))
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    this->is_generator()->set_s_hat(this->shat_channel->s());
		    if(!generate_is())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    if(!this->refresh_Ecm_hat())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    if(!base_type::generate_p(selected_branchings))
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		}
		else
		{
		    if(!generate_is())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    if(!this->refresh_Ecm_hat())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    if(this->backward_s_sampling())
		    {
			if(!generate_s_back(selected_branchings))
			{
			    this->weight()=(value_type)0;
			    return false;
			}
		    }
		    else
		    {
			if(!base_type::generate_s(selected_branchings))
			{
			    this->weight()=(value_type)0;
			    return false;
			}
		    }
		    if(!base_type::generate_t(selected_branchings))
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    if(!base_type::generate_p(selected_branchings))
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		}
		base_type::set_generated(selected_branchings);
		bool q=this->evaluate_weight();
		base_type::reset_generated(selected_branchings);
		return q;
	    }

	private:

	    /* Private modifiers */
	    /*-------------------*/

	    /* Adds a momentum channel to the list of no such with bit string bs
	     * exists, otherwise returns the respective list entry: */

	    momentum_channel_type* find_or_add_momentum_channel(bit_string_type bs)
	    {
		bitstring_compare<model_t,2,N_out,rng_t> pred(bs);
		typename momentum_channel_container::iterator it=std::find_if(this->momentum_channels.begin(),this->momentum_channels.end(),pred);
		if(it!=this->momentum_channels.end() and (*it)->bitstring==bs)
		{
		    return *it;
		}
		momentum_channel_type* c=NULL;
		if(bs.count()==1)
		{
		    size_type i=0;
		    while(!bs[i])
		    {
			++i;
		    }
		    c=new momentum_channel_type(bs,i==0?&(this->p_in(1)):&(this->p_out(i-1)));
		}
		else if(bs.count()==N_out+1)
		{
		    c=new momentum_channel_type(bs,&(this->p_in(0)));
		}
		else
		{
		    c=new momentum_channel_type(bs);
		}
		this->momentum_channels.insert(it,c);
		c->assign_mmin_components(this->eff_mmin);
		if(!bs[0] and bs.count()==N_out)
		{
		    shat_channel=c;
		}
		return c;
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
		particle_channel_type* ch_in=this->find_or_add_particle_channel(b_in,phi_in);
		particle_channel_type* ch_out1=this->find_or_add_particle_channel(b_out1,phi_out1);
		particle_channel_type* ch_out2=this->find_or_add_particle_channel(b_out2,phi_out2);
		
		if(ch_in==NULL or ch_out1==NULL or ch_out2==NULL)
		{
		    return;
		}
		if(ch_in->bitstring()[0])
		{
		    if(ch_out1->bitstring()[0])
		    {
			// swap channels such that ch_out1 is the timelike outgoing momentum
			std::swap(ch_out1,ch_out2);
		    }
		    if(ch_out2==this->incoming_particle_channels[1])
		    {
			// sum branching
			branching_type* sum_branching=branching_factory_type::create_sum_branching(ch_in,this->incoming_particle_channels[1],ch_out1);
			this->try_insert_branching_or_delete(ch_in,sum_branching);
		    }
		    else
		    {
			// t-branching
			// find or add the second outgoing s-channel
			bit_string_type b_out3=ch_out2->bitstring();
			b_out3.reset(0);
			momentum_channel_type* p_ch_out3=find_or_add_momentum_channel(b_out3);
			particle_channel_type* ch_out3;
			if(p_ch_out3->on_shell())
			{
			    ch_out3=*(p_ch_out3->begin_particle_channels());
			}
			else
			{
			    ch_out3=this->find_or_add_particle_channel(b_out3,NULL);
			}

			// find or add the total outgoing s-channel
			bit_string_type b_out4=ch_out1->bitstring()|ch_out3->bitstring();
			particle_channel_type* ch_out4=this->find_or_add_particle_channel(b_out4,NULL);

			branching_type* t_branching=branching_factory_type::create_t_branching(ch_in,this->incoming_particle_channels[1],ch_out1,ch_out2,ch_out3,ch_out4);
			this->try_insert_branching_or_delete(ch_in,t_branching);
		    }
		}
		else
		{
		    // s-branching
		    branching_type* s_branching=branching_factory_type::create_s_branching(ch_in,ch_out1,ch_out2);
		    this->try_insert_branching_or_delete(ch_in,s_branching);
		}
	    }

	    /* Checks whether the branching is at the end of a channel: */

	    bool static final_t_channel_sum(branching_type* branching)
	    {
		return !branching->s_type() and !branching->t_type() and !branching->incoming_channel->on_shell();
	    }

	    void replace_obsolete_sum_branchings()
	    {
		branching_container obsolete_sum_branchings;
		std::map<branching_type*,branching_container> t_branching_clones;

		for(typename branching_container::iterator it=this->ps_branchings.begin();it!=this->ps_branchings.end();++it)
		{
		    if(!final_t_channel_sum(*it)) continue;

		    create_final_t_channels(t_branching_clones,(*it)->incoming_channel,(*it)->channel(0));
		    obsolete_sum_branchings.push_back(*it);
		}
		for(typename branching_container::iterator it=obsolete_sum_branchings.begin();it!=obsolete_sum_branchings.end();++it)
		{
		    this->try_remove_branching_and_delete(*it);
		}

		for(typename std::map<branching_type*,branching_container>::iterator it=t_branching_clones.begin();it!=t_branching_clones.end();++it)
		{
		    bool replacing=true;
		    for(typename branching_container::iterator it2=it->second.begin();it2!=it->second.end();++it2)
		    {
			if(*it2==it->first)
			{
			    replacing=false;
			    continue;
			}
			else if(replacing)
			{
			    this->try_replace_branching_or_delete(it->first,*it2);
			    break;
			}
			else
			{
			    typename branching_container::iterator it3;
			    it3=std::find(this->ps_branchings.begin(),this->ps_branchings.end(),it->first);
			    this->try_insert_branching_or_delete((*it2)->incoming_channel,it3+1,*it2);
			}
		    }
		}
	    }

	    std::map<branching_type*,branching_container>& create_final_t_channels(std::map<branching_type*,branching_container>& mapping, particle_channel_type* t_channel,particle_channel_type* s2_channel)
	    {
		for(typename branching_container::iterator it=this->ps_branchings.begin();it!=this->ps_branchings.end();++it)
		{
		    if(!((*it)->t_type() and (*it)->channel(1)==t_channel)) continue;

		    t_branching<model_type,N_out,rng_t,spacetime_type>* t_branching_=static_cast<t_branching<model_type,N_out,rng_t,spacetime_type>*>(*it);
		    branching_container branchings;

		    if(t_channel->branching_count()>1)
		    {
			branchings.push_back(*it);
		    }
		    branchings.push_back(t_branching_->copy_to_final_branching(s2_channel));

		    mapping[t_branching_]=branchings;
		}
		return mapping;
	    }

	    /* Evaluates internal momentum: */

	    momentum_type p_internal(bit_string_type bs)
	    {
		momentum_type result;
		result.assign((value_type)0);
		for(size_type i=1;i<=N_out;++i)
		{
		    if(!bs[i]) continue;
		    result+=(this->outgoing_particle_channels[i-1]->p());
		}
		if(bs[0])
		{
		    result-=(this->incoming_particle_channels[1]->p());
		}
		return result;
	    }

	    /* Generates positive invariants for the given sequence of branchings: */

	    static bool generate_s_back(branching_container& selected_branchings)
	    {
		for(typename branching_container::reverse_iterator it=selected_branchings.rbegin();it!=selected_branchings.rend();++it)
		{
		    if((*it)->s_type())
		    {
			if(!(*it)->generate_s())
			{
			    return false;
			}
		    }
		}
		for(typename branching_container::iterator it=selected_branchings.begin();it!=selected_branchings.end();++it)
		{
		    if((*it)->t_type())
		    {
			if(!(*it)->generate_s())
			{
			    return false;
			}
		    }
		}
		return true;
	    }

	    /* Private data members */
	    /*----------------------*/

	    /* total momentum channel: */

	    momentum_channel_type* shat_channel;
    };
}

#endif /*CAMGEN_PS_TREE_H_*/

