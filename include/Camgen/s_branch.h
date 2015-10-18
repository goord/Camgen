//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_S_BRANCH_H_
#define CAMGEN_S_BRANCH_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the s-type momentum branching class. The        *
 * algorithm samples the 2 outgoing invariant masses and creates the 2 outgoing  *
 * momenta isotropically back-to-back in the incoming momentum rest-frame. Then  *
 * the output channels are boosted to the lab frame.                             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/Minkowski.h>
#include <Camgen/MC_config.h>
#include <Camgen/ps_branching.h>
#include <Camgen/uni_sphere.h>
#include <Camgen/Lorentz.h>
#include <Camgen/ps_vol.h>
#include <Camgen/ss_gen_fac.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t>class s_branching;

    /* Specialisation for Minkowski spacetimes: */
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,std::size_t D>class s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,D> >: public ps_branching<model_t,N_in,N_out,rng_t>
    {
	typedef ps_branching<model_t,N_in,N_out,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_t::value_type value_type;
	    typedef typename Minkowski_type::template implementation<value_type,model_t::dimension> spacetime_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::channel_type channel_type;
	    typedef typename base_type::ps_channel_type ps_channel_type;
	    typedef uniform_sphere<value_type,model_t::dimension-2,rn_engine> sphere_generator_type;
	    typedef s_pair_generator<value_type,rng_t> s_pair_generator_type;

	    /* Constructor: */

	    s_branching(const channel_type* in,channel_type* out1,channel_type* out2):base_type(in,out1,out2),ssgen(s_pair_generator_factory<value_type,rng_t>::create(this->s_generator(out1),this->s_generator(out2),in->s())),sphere_generator(){}

	    /* Destructor: */

	    virtual ~s_branching()
	    {
		delete ssgen;
	    }

	    /* Generates the positive invariant masses: */

	    bool generate_s()
	    {
		if(this->backward_s_sampling)
		{
		    if(this->incoming_channel->stot_channel() and !this->shat_sampling)
		    {
			sweight=(value_type)1;
			return true;
		    }
		    value_type mmin=this->m_out(0)+this->m_out(1);
		    if(this->incoming_channel->m_max()<mmin)
		    {
			sweight=(value_type)0;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    this->incoming_channel->set_m_min(mmin);
		    if(!this->incoming_channel->generate_s())
		    {
			sweight=(value_type)0;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    this->incoming_channel->set_status_s_generated();
		    sweight=this->incoming_channel->s_weight();
		}
		else
		{
		    if(!ssgen->generate())
		    {
			sweight=(value_type)0;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    this->channel(0)->set_status_s_generated();
		    this->channel(1)->set_status_s_generated();
		    sweight=ssgen->weight();
		}
		return true;
	    }

	    bool generate_t()
	    {
		return true;
	    }

	    /* Generates momenta with current invariant masses: */

	    bool generate_p()
	    {
		value_type s=this->s_in();
		value_type sqrts=this->m_in();
		value_type s1=this->s_out(0);
		value_type s2=this->s_out(1);
		value_type k12=Kallen(s,s1,s2);
		if(sqrts<=(value_type)0 or k12<(value_type)0)
		{
		    this->branching_weight=(value_type)0;
		    return false;
		}
		value_type lambda=std::sqrt(k12);
		value_type prefactor=(value_type)0.5/sqrts;
		(this->p_out(0))[0]=prefactor*(s+s1-s2);
		value_type phat=prefactor*lambda;
		sphere_generator.generate();
		for(size_type mu=1;mu<model_t::dimension;++mu)
		{
		    (this->p_out(0))[mu]=phat*(sphere_generator.object())[mu-1];
		}
		boost_from_restframe(this->p_out(0),this->p_in(),sqrts);
		this->channel(0)->evaluate_s();
		this->channel(0)->set_status_p_generated();

		this->p_out(1)=this->p_in()-this->p_out(0);
		this->channel(1)->evaluate_s();
		this->channel(1)->set_status_p_generated();

		this->branching_weight=massless_ps<value_type,2,model_t::dimension>::volume(sqrts)*std::pow(lambda/s,int(model_t::dimension-3))*sweight;
		
		
		return (this->branching_weight>(value_type)0);
	    }

	    /* Overridden weight evaluation method: */

            bool evaluate_branching_weight()
            {
		if(!evaluate_s_weight())
		{
		    this->branching_weight=(value_type)0;
		    return false;
		}
		value_type s=this->s_in();
		value_type k12=Kallen(s,this->s_out(0),this->s_out(1));
		if(s<=0 or k12<(value_type)0)
		{
		    this->branching_weight=(value_type)0;
		    return false;
		}
		value_type lambda=std::sqrt(k12);
		this->branching_weight=massless_ps<value_type,2,model_t::dimension>::volume(this->m_in())*std::pow(lambda/s,int(model_t::dimension-3))*sweight;
		return (this->branching_weight>(value_type)0);
            }

	    bool s_type() const
	    {
		return true;
	    }

	    /* Returns the branching topology type string: */

	    std::string type() const
	    {
		return "s";
	    }

	    /* Returns whether the branchings are equivalent: */

	    bool equiv(const ps_branching<model_t,N_in,N_out,rng_t>* other) const
	    {
		if(dynamic_cast<const s_branching<model_t,N_in,N_out,rng_t,spacetime_type>*>(other)!=NULL)
		{
		    return this->same_channels(other,false);
		}
		return false;
	    }

	private:

	    /* s-generation weight: */

	    value_type sweight;

	    /* S-pair generator: */

	    s_pair_generator_type* ssgen;

	    /* Solid angle generator: */

	    sphere_generator_type sphere_generator;

	    /* Evaluates sweight for backward s-sampling: */

	    bool evaluate_s_weight()
	    {
		if(this->backward_s_sampling)
		{
		    if(this->incoming_channel->shat_channel() and !this->shat_sampling)
		    {
			sweight=(value_type)1;
			return true;
		    }
		    value_type mmin=this->m_out(0)+this->m_out(1);
		    if(this->incoming_channel->m_max()<mmin)
		    {
			sweight=(value_type)0;
			return false;
		    }
		    this->incoming_channel->set_m_min(mmin);
		    if(!this->incoming_channel->evaluate_s_weight())
		    {
			sweight=(value_type)0;
			return false;
		    }
		    sweight=this->incoming_channel->s_weight();
		}
		else
		{
		    if(!ssgen->evaluate_weight())
		    {
			sweight=(value_type)0;
			return false;
		    }
		    sweight=ssgen->weight();
		}
		return true;
	    }
    };

    /* Specialisation for 4d-Minkowski spacetimes, adapting the polar decay
     * angle as well: */
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >: public ps_branching<model_t,N_in,N_out,rng_t>
    {
	typedef ps_branching<model_t,N_in,N_out,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef typename model_t::value_type value_type;
	    typedef typename Minkowski_type::template implementation<typename model_t::value_type,4> spacetime_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::channel_type channel_type;
	    typedef typename base_type::ps_channel_type ps_channel_type;
	    typedef s_pair_generator<value_type,rng_t> s_pair_generator_type;

	    /* Static data/methods: */

	    static const value_type twopi;
	    static const std::size_t long_dir=base_type::beam_direction;
	    static const std::size_t trans_dir1=(base_type::beam_direction==1)?2:1;
	    static const std::size_t trans_dir2=(base_type::beam_direction==3)?2:3;

	    /* Constructor: */

	    s_branching(channel_type* in,channel_type* out1,channel_type* out2):base_type(in,out1,out2),ssgen(s_pair_generator_factory<value_type,rng_t>::create(this->s_generator(out1),this->s_generator(out2),in->s())),theta_grid(NULL)
	    {
		if(adaptive_angles())
		{
		    theta_grid=new parni_generator<value_type,1,rng_t>(0,1,grid_bins(),grid_mode());
		}
	    }

	    /* Destructor: */

	    ~s_branching()
	    {
		delete ssgen;
		if(theta_grid!=NULL)
		{
		    delete theta_grid;
		}
	    }

	    /* Generates the positive invariant masses: */

	    bool generate_s()
	    {
		if(this->backward_s_sampling)
		{
		    if(this->incoming_channel->shat_channel() and !this->shat_sampling)
		    {
			sweight=(value_type)1;
			return true;
		    }
		    value_type mmin=std::max(this->m_out(0)+this->m_out(1),this->incoming_channel->m_min_min());
		    if(this->incoming_channel->m_max()<mmin)
		    {
			sweight=(value_type)0;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    this->incoming_channel->set_m_min(mmin);
		    if(!this->incoming_channel->generate_s())
		    {
			sweight=(value_type)0;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    this->incoming_channel->set_status_s_generated();
		    sweight=this->incoming_channel->s_weight();
		    return true;
		}
		else
		{
		    if(!ssgen->generate())
		    {
			sweight=(value_type)0;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    this->channel(0)->set_status_s_generated();
		    this->channel(1)->set_status_s_generated();
		    sweight=ssgen->weight();
		}
		return true;
	    }

	    bool generate_t()
	    {
		return true;
	    }

	    /* Generation implementation: */

	    bool generate_p()
	    {
		value_type s=this->s_in();
		value_type sqrts=this->m_in();
		value_type s1=this->s_out(0);
		value_type s2=this->s_out(1);
		value_type k12=Kallen(s,s1,s2);
		if(sqrts<=(value_type)0 or k12<(value_type)0)
		{
		    this->branching_weight=(value_type)0;
		    return false;
		}
		value_type lambda=std::sqrt(k12);
		value_type prefactor=(value_type)0.5/sqrts;
		(this->p_out(0))[0]=prefactor*(s+s1-s2);
		value_type phat=prefactor*lambda;
		value_type thetaweight(1);
		value_type r0(0);
		if(theta_grid!=NULL)
		{
		    if(!theta_grid->generate())
		    {
			this->branching_weight=(value_type)0;
			return false;
		    }
		    thetaweight=theta_grid->weight();
		    r0=theta_grid->object();
		}
		else
		{
		    r0=rn_stream::throw_number();
		}
		value_type costheta=(value_type)2*r0-(value_type)1;
		value_type sintheta=std::sqrt((value_type)1-costheta*costheta);
		value_type phi=rn_stream::throw_number(0,twopi);

		(this->p_out(0))[long_dir]=phat*costheta;
		(this->p_out(0))[trans_dir1]=phat*sintheta*std::cos(phi);
		(this->p_out(0))[trans_dir2]=phat*sintheta*std::sin(phi);
		
		boost_from_restframe(this->p_out(0),this->p_in(),sqrts);
		this->channel(0)->evaluate_s();
		this->channel(0)->set_status_p_generated();

		this->p_out(1)=this->p_in()-this->p_out(0);
		this->channel(1)->evaluate_s();
		this->channel(1)->set_status_p_generated();

		this->branching_weight=massless_ps<value_type,2,4>::volume(sqrts)*(lambda/s)*sweight*thetaweight;
		return (this->branching_weight>(value_type)0);
	    }

	    /* Overridden weight evaluation method: */

            bool evaluate_branching_weight()
            {
		if(!evaluate_s_weight())
		{
		    this->branching_weight=(value_type)0;
		    return false;
		}

		value_type s=this->s_in();
		value_type sqrts=this->m_in();
		value_type s1=this->s_out(0);
		value_type s2=this->s_out(1);

		value_type k12=Kallen(s,s1,s2);
		if(sqrts<=(value_type)0 or k12<(value_type)0)
		{
		    this->branching_weight=(value_type)0;
		    return false;
		}

		value_type lambda=std::sqrt(k12);
		value_type phat=(value_type)0.5*lambda/sqrts;
		value_type thetaweight(1);
		if(theta_grid!=NULL)
		{
		    value_type costheta=copy_boost_to_restframe(this->p_out(0),this->p_in(),sqrts)[long_dir]/phat;
		    value_type r0=(value_type)0.5*(costheta+(value_type)1);
		    if(r0<(value_type)0 or r0>(value_type)1)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"cos polar angle "<<costheta<<" out of range--returning weight 0."<<endlog;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    theta_grid->object()=r0;
		    if(!theta_grid->evaluate_weight())
		    {
			this->branching_weight=(value_type)0;
			return false;
		    }
		    thetaweight=(theta_grid->weight());
		}
		this->branching_weight=massless_ps<value_type,2,4>::volume(sqrts)*(lambda/s)*sweight*thetaweight;
		return (this->branching_weight>(value_type)0);
            }

	    bool s_type() const
	    {
		return true;
	    }

	    /* Returns the branching topology type: */

	    std::string type() const
	    {
		return "s";
	    }

	    /* Update method implementation: */

	    void update()
	    {
		if(theta_grid!=NULL)
		{
		    theta_grid->integrand()=this->integrand();
		    theta_grid->update();
		}
	    }

	    /* Grid adaption method implementation: */

	    void adapt_grids()
	    {
		if(theta_grid!=NULL)
		{
		    theta_grid->adapt();
		}
	    }

	    /* Returns whether the branchings are equivalent: */

	    bool equiv(const ps_branching<model_t,N_in,N_out,rng_t>* other) const
	    {
		if(dynamic_cast<const s_branching<model_t,N_in,N_out,rng_t,spacetime_type>*>(other)!=NULL)
		{
		    return this->same_channels(other,false);
		}
		return false;
	    }

	private:

	    /* Weight of the s-pair generation: */

	    value_type sweight;

	    /* s-pair generator: */

	    s_pair_generator_type* ssgen;

	    /* Polar angle adaptive grid: */

	    parni_generator<value_type,1,rng_t>* theta_grid;

	    /* Evaluates sweight for backward s-sampling: */

	    bool evaluate_s_weight()
	    {
		if(this->backward_s_sampling)
		{
		    if(this->incoming_channel->shat_channel() and !this->shat_sampling)
		    {
			sweight=(value_type)1;
			return true;
		    }
		    value_type mmin=this->m_out(0)+this->m_out(1);
		    if(this->incoming_channel->m_max()<mmin)
		    {
			sweight=(value_type)0;
			return false;
		    }
		    this->incoming_channel->set_m_min(mmin);
		    if(!this->incoming_channel->evaluate_s_weight())
		    {
			sweight=(value_type)0;
			return false;
		    }
		    sweight=this->incoming_channel->s_weight();
		}
		else
		{
		    if(!ssgen->evaluate_weight())
		    {
			sweight=(value_type)0;
			return false;
		    }
		    sweight=ssgen->weight();
		}
		return true;
	    }
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const typename model_t::value_type s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::twopi=(typename model_t::value_type)2*std::acos(-(typename model_t::value_type)1);
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::long_dir;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::trans_dir1;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::trans_dir2;
}

#endif /*CAMGEN_S_BRANCH_H_*/

