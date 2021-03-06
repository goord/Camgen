//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_tree_base.h
    \brief Base class for phase space decomposition tree.
 */

#ifndef CAMGEN_PS_TREE_BASE_H_
#define CAMGEN_PS_TREE_BASE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility base class for phase space decomposition generators.*
 *                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <Camgen/rn_strm.h>
#include <Camgen/ps_branching_fac.h>
#include <Camgen/ps_gen.h>

namespace Camgen
{
    /* Comparison functor for momentum channel ordering: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class bitstring_compare
    {
	typedef bit_string<N_in+N_out-1> bit_string_type;
    };

    /* Comparison functor for momentum channel ordering in 1->N_out processes: */

    template<class model_t,std::size_t N_out,class rng_t>class bitstring_compare<model_t,1,N_out,rng_t>
    {
	public:

	    typedef momentum_channel<model_t,1,N_out,rng_t> momentum_channel_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;
	    typedef std::size_t size_type;

	    const bit_string_type bitstring;

	    bitstring_compare(const bit_string_type& bitstring_):bitstring(bitstring_){}

	    bool operator()(const momentum_channel_type* p) const
	    {
		return (bitstring>=p->bitstring);
	    }
    };

    /* Comparison functor for momentum channel ordering in 2->N_out processes: */

    template<class model_t,std::size_t N_out,class rng_t>class bitstring_compare<model_t,2,N_out,rng_t>
    {
	public:

	    typedef momentum_channel<model_t,2,N_out,rng_t> momentum_channel_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;
	    typedef std::size_t size_type;

	    const bit_string_type bitstring;

	    bitstring_compare(const bit_string_type& bitstring_):bitstring(bitstring_){}

	    bool operator()(const momentum_channel_type* p) const
	    {
		return less(bitstring,p->bitstring);
	    }

	    static bool less(const bit_string_type& b1,const bit_string_type& b2)
	    {
		if(b1==b2) return true;

		size_type n=b1.count();
		size_type m=b2.count();

		if(n>N_out)
		{
		    return true;
		}
		if(m>N_out)
		{
		    return false;
		}
		if(n==1 && b1[0])
		{
		    return true;
		}
		if(m==1 && b2[0])
		{
		    return false;
		}
		return b1>=b2;
	    }
    };

    /// Phase space decomposition tree base class.

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_tree_base: public ps_generator<model_t,N_in,N_out>
    {
	typedef ps_generator<model_t,N_in,N_out> base_type;

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
	    typedef momentum_channel<model_t,N_in,N_out,rng_t> momentum_channel_type;
	    typedef particle_channel<model_t,N_in,N_out,rng_t> particle_channel_type;
	    typedef ps_branching<model_t,N_in,N_out,rng_t> branching_type;
	    typedef typename std::vector<momentum_channel_type*> momentum_channel_container;
	    typedef typename std::vector<particle_channel_type*> particle_channel_container;
	    typedef typename std::vector<branching_type*> branching_container;
	    typedef ps_branching_factory<model_t,N_in,N_out,rng_t,typename model_t::spacetime_type> branching_factory_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;

	    typedef particle<model_t> particle_type;

	    /// Determines whether invariant masses are generated in backward
	    /// direction.
	    
	    const bool backward_s_gen;

	    /// Constructor:

	    ps_tree_base(init_state_type* is,bool backward_s_gen_):base_type(is),backward_s_gen(backward_s_gen_){}

	    /// Destructor:
	    
	    virtual ~ps_tree_base()
	    {
		delete_tree();
		for(typename std::vector<const particle_type*>::iterator it=auxiliary_particles.begin();it!=auxiliary_particles.end();++it)
		{
		    delete *it;
		}
	    }

	    /// Evaluates the internal momentum

	    virtual momentum_type p_internal(bit_string_type bs)=0;

	    /// Returns whether the generator will sample invariants in backward order.

	    bool backward_s_sampling() const
	    {
		return backward_s_gen;
	    }

	    /// Returns whether the generator will sample the total partonic invariant mass.

	    bool backward_shat_sampling() const
	    {
		return !this->is_generator()->s_hat_sampling;
	    }

	    /// Weight evaluation method.

	    bool evaluate_fs_weight()
	    {
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if((*it)->get_status()==momentum_channel_type::p_set) continue;

		    (*it)->p()=p_internal((*it)->bitstring);
		    (*it)->evaluate_s();
		    (*it)->set_status_p_generated();
		}
		typename branching_container::iterator it=ps_branchings.begin();
		const particle_channel_type* incoming_channel=(*it)->incoming_channel;
		while(it!=ps_branchings.end())
		{
		    if((*it)->incoming_channel!=incoming_channel)
		    {
			if(!const_cast<particle_channel_type*>(incoming_channel)->evaluate_weight())
			{
			    this->fs_weight()=(value_type)0;
			    return false;
			}
			incoming_channel=(*it)->incoming_channel;
		    }
		    if(!(*it)->evaluate_weight(!(*it)->generated))
		    {
			this->fs_weight()=(value_type)0;
			return false;
		    }
		    ++it;
		}
		this->incoming_particle_channels[0]->evaluate_weight();
		this->fs_weight()=this->incoming_particle_channels[0]->weight();
		return true;
	    }

            /// Overrides the MC_generator generate_unweighted method.

            bool generate_unweighted()
            {
                this->MC_integrator<value_type>::template generate_unweighted<rng_t>();
                return true;
            }

	    /// Overrides the update method.

	    void update()
	    {
		this->base_type::update();
		value_type f=(this->integrand())*(this->weight());
		particle_channel_update update_callback(f);
		apply_to_particle_channels(update_callback);
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    static_cast<MC_integrator_base<value_type>*>(*it)->update(f);
		}
	    }

	    /// Overrides the parni-grid adaptation method.

	    void adapt_grids()
	    {
		this->base_type::adapt_grids();
		apply_to_particle_channels(particle_channel_adapt_grids);
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    (*it)->adapt_grids();
		}
	    }
	    
	    /// Overrides the multichannel adaptation method.

	    void adapt_channels()
	    {
		this->base_type::adapt_channels();
		apply_to_particle_channels(particle_channel_adapt_channels);
		particle_channel_clean_branchings cleanup_callback;
		apply_to_particle_channels(cleanup_callback);

		typename branching_container::iterator end_it=ps_branchings.end();
		for(typename branching_container::iterator it=cleanup_callback.branchings.begin();it!=cleanup_callback.branchings.end();++it)
		{
		    end_it=std::remove(ps_branchings.begin(),end_it,*it);
		}
		ps_branchings.erase(end_it,ps_branchings.end());
	    }

	    /// Refreshes the internal parameters.

	    bool refresh_params()
	    {
		particle_channel_refresh_params refresh_callback;
		apply_to_particle_channels(refresh_callback);
		return refresh_callback.result;
	    }

	    /// Resets cross-sections, multichannel weights and adaptive grids.

	    void reset()
	    {
		this->base_type::reset();
		apply_to_particle_channels(particle_channel_reset);
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    (*it)->reset();
		}
	    }

	    /// Refreshes the minimal invariant masses.

	    bool refresh_m_min()
	    {
		if(!this->base_type::refresh_m_min())
		{
		    return false;
		}
                value_type min0=this->M_in(0);
                value_type min1=N_in>1?this->M_in(1):min0;

		value_type sa=min0*min0;
		value_type sb=min1*min1;

		bool success=true;
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if((*it)->on_shell())
		    {
			success&=((*it)->set_s_min_min(0));
			success&=((*it)->set_s_min(0));
		    }
		    else
		    {
			success&=((*it)->refresh_s_min(this->s_tot(),sa,sb));
		    }
		}
		return success;
	    }

	    /// Refreshes the total hadronic invariant mass.

	    bool refresh_Ecm()
	    {
		if(!this->base_type::refresh_Ecm())
		{
		    return false;
		}
		if(N_in==2 and !this->check_sufficient_s())
		{
		    return false;
		}
		value_type s=this->Ecm()*this->Ecm();
		bool success=true;
		for(size_type i=0;i<N_in;++i)
		{
		    success&=incoming_particle_channels[i]->set_s_max_max(s);
		    success&=incoming_particle_channels[i]->refresh_params();
		}
		for(size_type i=0;i<N_out;++i)
		{
		    success&=outgoing_particle_channels[i]->set_s_max_max(s);
		    success&=outgoing_particle_channels[i]->refresh_params();
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if((*it)->shat_channel())
		    {
			if(this->backward_shat_sampling())
			{
                            success&=((*it)->set_max_s_range((value_type)0,s));
			}
			continue;
		    }
		    if((*it)->on_shell())
		    {
			continue;
		    }
		    if((*it)->spacelike())
		    {
                        success&=((*it)->set_max_s_range(-s,(value_type)0));
		    }
		    else
		    {
                        success&=((*it)->set_max_s_range((value_type)0,s));
		    }
		}
		return success;
	    }

	    /// Refreshes the partonic invariant mass.

	    bool refresh_Ecm_hat()
	    {
		if(!this->base_type::refresh_Ecm_hat())
		{
		    return false;
		}
		value_type shat=this->Ecm_hat()*this->Ecm_hat();
		bool success=true;
		if(N_in==1)
		{
		    incoming_particle_channels[0]->s()=this->is_generator()->s_hat();
		}
		else
		{
		    for(size_type i=0;i<N_in;++i)
		    {
			success&=incoming_particle_channels[i]->set_s_max(shat);
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    success&=outgoing_particle_channels[i]->set_s_max(shat);
		    if(backward_s_sampling())
		    {
			success&=outgoing_particle_channels[i]->generate_s();
		    }
		    outgoing_particle_channels[i]->weight()=(value_type)1;
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if((*it)->shat_channel())
		    {
			success&=(*it)->refresh_params();
			(*it)->s()=shat;
			continue;
		    }
		    if((*it)->on_shell())
		    {
			continue;
		    }
		    if((*it)->spacelike())
		    {
			success&=((*it)->set_s_min(-shat));
		    }
		    else
		    {
			success&=((*it)->set_s_max(shat));
		    }
		}
		return success;
	    }

	    /// Sets the minimal invariant mass i,j and refreshes depending channels.

	    bool set_m_min(int i,int j,const value_type* sqrts)
	    {
		if(!this->base_type::set_m_min(i,j,sqrts)) return false;

                value_type min0=this->get_event().M_in(0);
                value_type min1=N_in>1?this->get_event().M_in(1):min0;

		value_type sa=min0*min0;
		value_type sb=min1*min1;
		
                bool success=true;
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    if((*it)->on_shell())
		    {
			continue;
		    }
		    bit_string_type bs=(*it)->bitstring;
		    if(bs[i+N_in-2] and bs[j+N_in-2])
		    {
			success&=((*it)->refresh_s_min(this->s_tot(),sa,sb));
		    }
		}
		return success;
	    }

	    /// Returns whether we have branchings or not.

	    bool empty() const
	    {
		return ps_branchings.size()==0;
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
	    
	    /// Printing method.
	    
	    std::ostream& print(std::ostream& os) const
	    {
		for(typename momentum_channel_container::const_iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->print(os);
		}
		return os;
	    }

	    /// Generation channel printing method.

	    std::ostream& print_channel(std::ostream& os) const
	    {
		return incoming_particle_channels[0]->print_channel(os);
	    }


	protected:

	    /* Protected modifiers */
	    /*---------------------*/

	    /* Adds momentum and particle channels for the external particles: */

	    void set_external_particles(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		if(N_in==0) return;

		bit_string_type b_in0;
		b_in0.set();
		const particle_type* phi_in0=it->get_phase_space(0)->particle_type;
		this->incoming_particle_channels[0]=this->find_or_add_particle_channel(b_in0,phi_in0);
		
		for(size_type i=1;i<N_in;++i)
		{
		    bit_string_type b_in;
		    b_in.set(i-1);
		    const particle_type* phi_in=it->get_phase_space(i)->particle_type;
		    this->incoming_particle_channels[i]=this->find_or_add_particle_channel(b_in,phi_in);
		}

		for(size_type i=0;i<N_out;++i)
		{
		    bit_string_type b_out;
		    b_out.set(N_in+i-1);
		    const particle_type* phi_out=it->get_phase_space(N_in+i)->particle_type;
		    this->outgoing_particle_channels[i]=this->find_or_add_particle_channel(b_out,phi_out);
		}
	    }

	    /* Returns the momentum channel with bitstring bs, returns NULL if not
	     * found. */

	    momentum_channel_type* find_momentum_channel(bit_string_type bs)
	    {
		bitstring_compare<model_t,N_in,N_out,rng_t>pred(bs);
		typename momentum_channel_container::iterator it=std::find_if(momentum_channels.begin(),momentum_channels.end(),pred);
		return it==momentum_channels.end()?NULL:((*it)->bitstring!=bs?NULL:(*it));;
	    }

	    /* Returns the particle channel with bitstring bitstr, returns NULL if not
	     * found. */

	    particle_channel_type* find_particle_channel(bit_string_type bs,const particle_type* phi)
	    {
		momentum_channel_type* pch=find_momentum_channel(bs);
		return (pch==NULL)?NULL:(pch->find_particle_channel(phi));
	    }

	    /* Adds a particle channel to the list of no such with bit string bs
	     * exists, otherwise returns the respective list entry: */

	    particle_channel_type* find_or_add_particle_channel(bit_string_type bs,const particle_type* phi)
	    {
		const value_type* sqrtshat=this->is_generator()->s_hat_sampling?this->get_Ecm_hat():NULL;
		momentum_channel_type* c=find_or_add_momentum_channel(bs);
		particle_channel_type* result=branching_factory_type::find_or_add_particle_channel(c,phi,sqrtshat);
		if(bs.count()==N_in+N_out-1)
		{
		    incoming_particle_channels[0]=result;
		}
		else if(bs.count()==1)
		{
		    size_type i=0;
		    for(;i<bs.size();++i)
		    {
			if(bs[i])
			{
			    break;
			}
		    }
		    if(N_in==1)
		    {
			outgoing_particle_channels[i]=result;
		    }
		    if(N_in==2)
		    {
			if(i==0)
			{
			    incoming_particle_channels[1]=result;
			}
			else
			{
			    outgoing_particle_channels[i-1]=result;
			}
		    }
		}
		return result;
	    }

	    virtual momentum_channel_type* find_or_add_momentum_channel(bit_string_type bs)=0;

	    /* Recursively selects a branching sequence: */

	    branching_container select_branchings()
	    {
		branching_container result;
		if(!select_branchings(incoming_particle_channels[0],result))
		{
		    result.clear();
		}
		return result;
	    }

	    /* Inserts the branching in the list, provided it isn't already
	     * there. If so deletes the instance. */

	    bool try_insert_branching_or_delete(particle_channel_type* channel,branching_type* branching)
	    {
		if(branching==NULL)
		{
		    return false;
		}
		if(channel->insert_branching(branching) and std::find(ps_branchings.begin(),ps_branchings.end(),branching)==ps_branchings.end())
		{
		    ps_branchings.push_back(configure_branching(branching));
		    return true;
		}
		else
		{
		    delete branching;
		}
		return false;
	    }

	    /* Inserts the branching in the list, provided it isn't already
	     * there. If so deletes the instance. */

	    bool try_insert_branching_or_delete(particle_channel_type* channel,typename branching_container::iterator it,branching_type* branching)
	    {
		if(branching==NULL)
		{
		    return false;
		}
		if(channel->insert_branching(branching) and std::find(ps_branchings.begin(),ps_branchings.end(),branching)==ps_branchings.end())
		{
		    ps_branchings.insert(it,configure_branching(branching));
		    return true;
		}
		else
		{
		    delete branching;
		}
		return false;
	    }

	    /* Tries to remove the branching everywhere in the tree and deletes the instance: */

	    bool try_remove_branching_and_delete(branching_type* branching)
	    {
		typename branching_container::iterator it=std::find(ps_branchings.begin(),ps_branchings.end(),branching);
		if(it==ps_branchings.end())
		{
		    return false;
		}
		particle_channel_type* channel=branching->incoming_channel;
		channel->remove_branching(branching);
		delete branching;
		*it=NULL;
		ps_branchings.erase(it);
		return true;
	    }

	    /* Tries to replace the first branching with the second. If succesful, deletes the first branching.
	     * Otherwise, deletes the second branching: */

	    bool try_replace_branching_or_delete(branching_type* first,branching_type* second)
	    {
		typename branching_container::iterator it=std::find(ps_branchings.begin(),ps_branchings.end(),first);
		if(it==ps_branchings.end())
		{
		    delete second;
		    return false;
		}
		particle_channel_type* channel=first->incoming_channel;
		
		if(channel->replace_branching(first,second))
		{
		    *it=configure_branching(second);
		    delete first;
		    return true;
		}
		return false;
	    }

	    /* De-allocates the complete tree: */

	    void delete_tree()
	    {
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    delete *it;
		}
		momentum_channels.clear();
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    delete *it;
		}
		ps_branchings.clear();
		incoming_particle_channels.assign(NULL);
		outgoing_particle_channels.assign(NULL);
	    }

	    /* Cleans disjoint pieces from the tree: */

	    void clean_tree()
	    {
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    if(dead_end_branching(*it))
		    {
			particle_channel_type* channel=(*it)->incoming_channel;
			channel->remove_branching(*it);
			delete *it;
			*it=NULL;
		    }
		}
		typename branching_container::iterator it_end=std::remove(ps_branchings.begin(),ps_branchings.end(),static_cast<branching_type*>(NULL));
		ps_branchings.erase(it_end,ps_branchings.end());
		for(typename branching_container::iterator it=ps_branchings.begin();it!=ps_branchings.end();++it)
		{
		    particle_channel_type* channel=(*it)->incoming_channel;
		    if(channel->on_shell())
		    {
			continue;
		    }
		    bool disconnected=true;
		    for(typename branching_container::iterator it2=ps_branchings.begin();it2!=ps_branchings.end();++it2)
		    {
			if(*it2!=NULL and std::find((*it2)->begin_channels_out(),(*it2)->end_channels_out(),channel)!=(*it2)->end_channels_out())
			{
			    disconnected=false;
			    break;
			}
		    }
		    if(disconnected)
		    {
			channel->remove_branching(*it);
			delete *it;
			*it=NULL;
		    }
		}
		it_end=std::remove(ps_branchings.begin(),ps_branchings.end(),static_cast<branching_type*>(NULL));
		ps_branchings.erase(it_end,ps_branchings.end());
	    }

	    /* Resets all generation flags: */

	    void initialise_channels()
	    {
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    (*it)->p().assign((value_type)0);
		    (*it)->s()=(value_type)0;
		    (*it)->reset_status();
		}
	    }

	    /* Creates a new auxiliary particle instance: */

	    particle_type* create_auxiliary_particle()
	    {
		std::string name("rho");
		name.append(make_string(auxiliary_particles.size()));
		particle_type* rho=particle_type::create_auxiliary(name);
		auxiliary_particles.push_back(rho);
		return rho;
	    }

	    /* Stets correct flags to the given branching: */

	    branching_type* configure_branching(branching_type* branching) const
	    {
		branching->backward_s_sampling=this->backward_s_sampling();
		branching->shat_sampling=backward_shat_sampling();
		return branching;
	    }

	    /* Utility method, applies the functional to all particle channels: */

	    template<class T>void apply_to_particle_channels(T& func)
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    func(incoming_particle_channels[i]);
		}
		for(typename momentum_channel_container::iterator it=momentum_channels.begin();it!=momentum_channels.end();++it)
		{
		    for(typename momentum_channel_type::particle_channel_iterator p_it=(*it)->begin_particle_channels();p_it!=(*it)->end_particle_channels();++p_it)
		    {
			func(*p_it);
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    func(outgoing_particle_channels[i]);
		}
	    }

	    /* Protected data members */
	    /*----------------------*/

	    /* Incoming particle channel: */

	    vector<particle_channel_type*,N_in> incoming_particle_channels;

	    /* Outgoing particle channels: */

	    vector<particle_channel_type*,N_out> outgoing_particle_channels;
	    
	    /* Internal momentum phase space channels: */
	    
	    momentum_channel_container momentum_channels;
	    
	    /* List of phase space branchings: */
	    
	    branching_container ps_branchings;

	    /* Extra auxiliary particles: */

	    std::vector<const particle_type*> auxiliary_particles;

	    /* Protected static methods: */
	    /*---------------------------*/

	    /* Recursively selects a sequence of branchings for the given channel: */

	    static bool select_branchings(particle_channel_type* channel,branching_container& branchings)
	    {
		if(channel->branching_count()==0)
		{
		    return true;
		}
		branching_type* selected_branching=channel->select_branching();
		if(selected_branching==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Failed to select branching in channel "<<channel->name<<endlog;
		    return false;
		}
		branchings.push_back(selected_branching);
		for(typename branching_type::channel_iterator it=selected_branching->begin_channels_out();it!=selected_branching->end_channels_out();++it)
		{
		    if(!select_branchings(*it,branchings))
		    {
			return false;
		    }
		}

		return true;
	    }

	    /* Generates positive invariants for the given sequence of branchings: */

	    static bool generate_s(branching_container& selected_branchings)
	    {
		for(typename branching_container::iterator it=selected_branchings.begin();it!=selected_branchings.end();++it)
		{
		    if(!(*it)->generate_s())
		    {
			return false;
		    }
		}
		return true;
	    }

	    /* Generates negative invariants for the given sequence of branchings: */

	    static bool generate_t(branching_container& selected_branchings)
	    {
		for(typename branching_container::iterator it=selected_branchings.begin();it!=selected_branchings.end();++it)
		{
		    if(!(*it)->generate_t())
		    {
			return false;
		    }
		}
		return true;
	    }

	    /* Generates momentum for the given sequence of branchings: */

	    static bool generate_p(branching_container& selected_branchings)
	    {
		for(typename branching_container::iterator it=selected_branchings.begin();it!=selected_branchings.end();++it)
		{
		    if(!(*it)->generate_p())
		    {
			return false;
		    }
		}
		return true;
	    }

	    /* Sets the generation flag of the branchings: */

	    static void set_generated(branching_container& selected_branchings)
	    {
		for(typename branching_container::iterator it=selected_branchings.begin();it!=selected_branchings.end();++it)
		{
		    (*it)->generated=true;
		}
	    }

	    /* Resets the generation flag of the branchings: */

	    static void reset_generated(branching_container& selected_branchings)
	    {
		for(typename branching_container::iterator it=selected_branchings.begin();it!=selected_branchings.end();++it)
		{
		    (*it)->generated=false;
		}
	    }

	    /* Checks whether the branching does not continue: */

	    static bool dead_end_branching(const branching_type* branching)
	    {
		for(typename branching_type::const_channel_iterator it=branching->begin_channels_out();it!=branching->end_channels_out();++it)
		{
		    if(!(*it)->on_shell() and (*it)->branching_count()==0)
		    {
			return true;
		    }
		}
		return false;
	    }

	    /* Callback for grid adaptation: */

	    static void particle_channel_adapt_grids(particle_channel_type* p)
	    {
		p->adapt_grids();
	    }

	    /* Callback for channel adaptation: */

	    static void particle_channel_adapt_channels(particle_channel_type* p)
	    {
		p->adapt_channels();
	    }

	    /* Callback for reset: */

	    static void particle_channel_reset(particle_channel_type* p)
	    {
		p->reset();
	    }

	    /* Functor for refreshing parameters: */
	    
	    class particle_channel_refresh_params
	    {
		public:

		    bool result;
		    
		    particle_channel_refresh_params():result(true){}
		    
		    void operator()(particle_channel_type* p)
		    {
			result&=p->refresh_params();
		    }
	    };

	    /* Functor for updating weight sums: */

	    class particle_channel_update
	    {
		public:

		    value_type w;

		    particle_channel_update(const value_type& w_):w(w_){}

		    void operator()(particle_channel_type* p)
		    {
			static_cast<MC_integrator_base<value_type>*>(p)->update(w);
		    }
	    };

	    /* Functor for cleaning branchings: */

	    class particle_channel_clean_branchings
	    {
		public:

		    std::vector<branching_type*> branchings;

		    void operator()(particle_channel_type* p)
		    {
			std::vector<branching_type*>b=p->clean_branchings();
			branchings.insert(branchings.end(),b.begin(),b.end());
		    }
	    };
    };
}

#endif /*CAMGEN_PS_TREE_BASE_H_*/

