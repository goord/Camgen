//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_S_BRANCH_BACK_H_
#define CAMGEN_S_BRANCH_BACK_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the backward s-type momentum branching class.   *
 * The algorithm samples the incoming invariant mass and creates the 2 outgoig   *
 * momenta isotropically back-to-back in the incoming momentum rest-frame. Then  *
 * the output channels are boosted to the lab frame.                             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/Minkowski.h>
#include <Camgen/branching.h>
#include <Camgen/uni_sphere.h>
#include <Camgen/Lorentz.h>
#include <Camgen/ps_vol.h>
#include <Camgen/parni.h>

namespace Camgen
{
    /* Specialisation for Minkowski spacetimes: */
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,std::size_t D>class backward_s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,D> >: public ps_branching<model_t,N_in,N_out,rng_t>
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
	    typedef typename base_type::s_generator_type s_gen_type;
	    typedef uniform_sphere<value_type,model_t::dimension-2,rn_engine> sphere_gen_type;

	    /* On-shell flags: */

	    const bool s1_on_shell;
	    const bool s2_on_shell;

	    /* Constructor: */

	    backward_s_branching(channel_type* in,channel_type* out1,channel_type* out2):base_type(in,2),s1_on_shell(out1->on_shell()),s2_on_shell(out2->on_shell()),sphere_generator(new sphere_gen_type(&phat))
	    {
		this->channels[0]=out1;
		this->channels[1]=out2;
	    }

	    /* Destructor: */

	    ~backward_s_branching()
	    {
		delete sphere_generator;
	    }

	    /* Generates the incoming invariant mass: */

	    bool generate_s()
	    {
		if(s1_on_shell)
		{
		    this->channel(0)->generate_s();
		}
		if(s2_on_shell)
		{
		    this->channel(1)->generate_s();
		}
		if(!refresh_s_bounds())
		{
		    sweight=(value_type)0;
		    return false;
		}
		if(!this->incoming_channel->generate_s())
		{
		    sweight=(value_type)0;
		    return false;
		}
		sweight=this->incoming_channel->s_weight();
		return (sweight>(value_type)0);
	    }

	    /* Generates momenta with current invariant masses: */

	    bool generate_p()
	    {
		value_type s=this->s_in();
		value_type sqrts=this->m_in();
		value_type s1=this->channel(0)->s();
		value_type s2=this->channel(1)->s();
		value_type lambda=std::sqrt(Kallen(s,s1,s2));
		value_type prefactor=(value_type)0.5/sqrts;
		(this->p_out(0))[0]=prefactor*(s+s1-s2);
		phat=prefactor*lambda;
		sphere_generator->generate();
		for(size_type mu=1;mu<model_t::dimension;++mu)
		{
		    (this->p_out(0))[mu]=(sphere_generator->object())[mu-1];
		}
		boost_from_restframe(this->p_out(0),this->p_in(),sqrts);
		this->p_out(1)=this->p_in()-this->p_out(0);
		this->weight()=massless_ps<value_type,2,model_t::dimension>::volume(sqrts)*std::pow(lambda/s,int(model_t::dimension-3))*sweight;
		return (this->weight()>(value_type)0);
	    }

	    /* Weight evaluation method: */

	    bool evaluate_branching_weight()
	    {
		if(!refresh_s_bounds())
		{
		    sweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type sqrts=this->m_in();
		if((this->m_out(0)+this->m_out(1)>sqrts) or !this->incoming_channel->evaluate_s_weight())
		{
		    sweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		sweight=this->incoming_channel->s_weight();
		value_type s=this->s_in();
		value_type lambda=std::sqrt(Kallen(s,this->s_out(0),this->s_out(1)));
		this->weight()=massless_ps<value_type,2,model_t::dimension>::volume(sqrts)*std::pow(lambda/s,int(model_t::dimension-3))*sweight;
		return (this->weight()>(value_type)0);
	    }

	    /* Returns the branching topology type string: */

	    std::string type() const
	    {
		return "s~";
	    }

	    /* Returns whether the branchings are equivalent: */

	    bool equiv(const ps_branching<model_t,N_in,N_out,rng_t>* other) const
	    {
		if(dynamic_cast<const backward_s_branching<model_t,N_in,N_out,rng_t,spacetime_type>*>(other)!=NULL)
		{
		    return this->same_channels(other,false);
		}
		return false;
	    }

	protected:

	    /* Updates minimal invariant masses backward: */

	    bool refresh_s_bounds()
	    {
		value_type mmin=std::max(this->m_out(0)+this->m_out(1),this->incoming_channel->m_min());
		if(mmin>this->incoming_channel->m_max())
		{
		    return false;
		}
		if(!(this->incoming_channel->set_m_min(mmin)))
		{
		    this->incoming_channel->s_warning();
		    return false;
		}
		return true;
	    }

	private:

	    /* s-generation weight: */

	    value_type sweight;

	    /* Solid angle generator: */

	    sphere_gen_type* sphere_generator;

	    /* Rest-frame momentum of decay product: */

	    value_type phat;
    };

    /* Specialisation for 4d-Minkowski spacetimes, adapting the polar decay
     * angle as well: */
    
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class backward_s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >: public ps_branching<model_t,N_in,N_out,rng_t>
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
	    typedef typename base_type::s_generator_type s_gen_type;

	    /* Static data/methods: */

	    static const value_type twopi;
	    static const std::size_t long_dir=base_type::beam_direction;
	    static const std::size_t trans_dir1=(base_type::beam_direction==1)?2:1;
	    static const std::size_t trans_dir2=(base_type::beam_direction==3)?2:3;

	    /* On-shell flags: */

	    const bool s1_on_shell;
	    const bool s2_on_shell;

	    /* Constructor: */

	    backward_s_branching(channel_type* in,channel_type* out1,channel_type* out2):base_type(in,2),s1_on_shell(out1->on_shell()),s2_on_shell(out2->on_shell()),theta_grid(NULL)
	    {
		this->channels[0]=out1;
		this->channels[1]=out2;
		if(adaptive_angles())
		{
		    theta_grid=new parni<value_type,1,rng_t>(&r0,0,1,grid_bins(),grid_mode());
		}
	    }

	    /* Destructor: */

	    ~backward_s_branching()
	    {
		if(theta_grid!=NULL)
		{
		    delete theta_grid;
		}
	    }

	    /* Generates the incoming invariant mass: */

	    bool generate_s()
	    {
		if(s1_on_shell)
		{
		    this->channel(0)->generate_s();
		}
		if(s2_on_shell)
		{
		    this->channel(1)->generate_s();
		}
		if(!refresh_s_bounds())
		{
		    sweight=(value_type)0;
		    return false;
		}
		if(!this->incoming_channel->generate_s())
		{
		    sweight=(value_type)0;
		    return false;
		}
		sweight=this->incoming_channel->s_weight();
		return (sweight>(value_type)0);
	    }

	    /* Generation implementation: */

	    bool generate_p()
	    {
		value_type s=this->s_in();
		value_type sqrts=this->m_in();
		value_type s1=this->s_out(0);
		value_type s2=this->s_out(1);
		value_type lambda=std::sqrt(Kallen(s,s1,s2));
		value_type prefactor=(value_type)0.5/sqrts;
		(this->p_out(0))[0]=prefactor*(s+s1-s2);
		phat=prefactor*lambda;
		value_type thetaweight(1);
		if(theta_grid!=NULL)
		{
		    if(!theta_grid->generate())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    thetaweight=theta_grid->weight();
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
		this->p_out(1)=this->p_in()-this->p_out(0);
		this->weight()=massless_ps<value_type,2,4>::volume(sqrts)*(lambda/s)*sweight*thetaweight;
		return (this->weight()>(value_type)0);
	    }

	    /* Overridden weight evaluation method: */

            bool evaluate_branching_weight()
            {
		if(!refresh_s_bounds())
		{
		    sweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type sqrts=this->m_in();
		if((this->m_out(0)+this->m_out(1)>sqrts) or !this->incoming_channel->evaluate_s_weight())
		{
		    sweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		sweight=this->incoming_channel->s_weight();
		value_type s=this->s_in();
		value_type lambda=std::sqrt(Kallen(s,this->s_out(0),this->s_out(1)));
		phat=(value_type)0.5*lambda/sqrts;
		value_type rweight(1);
		if(theta_grid!=NULL)
		{
		    value_type costheta=copy_boost_to_restframe(this->p_out(0),this->p_in(),sqrts)[long_dir]/phat;
		    r0=(value_type)0.5*(costheta+(value_type)1);
		    r0=std::max(r0,(value_type)0);
		    r0=std::min(r0,(value_type)1);
		    if(r0<(value_type)0 or r0>(value_type)1)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"cos polar angle "<<costheta<<" out of range--returning weight 0."<<endlog;
			this->weight()=(value_type)0;
			return false;
		    }
		    if(!theta_grid->evaluate_weight())
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    rweight=(theta_grid->weight());
		}
		this->weight()=massless_ps<value_type,2,model_t::dimension>::volume(sqrts)*std::pow(lambda/s,int(model_t::dimension-3))*sweight*rweight;
		return (this->weight()>(value_type)0);
            }

	    /* Returns the branching topology type: */

	    std::string type() const
	    {
		return "s~";
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
		if(dynamic_cast<const backward_s_branching<model_t,N_in,N_out,rng_t,spacetime_type>*>(other)!=NULL)
		{
		    return this->same_channels(other,false);
		}
		return false;
	    }

	    /* Subclass loading method: */

	    std::istream& load_data(std::istream& is)
	    {
		if(theta_grid!=NULL)
		{
		    delete theta_grid;
		}
		is>>std::ws;
		if(is.peek()=='N')
		{
		    is.ignore(1);
		    theta_grid=NULL;
		}
		else
		{
		    theta_grid=parni<value_type,1,rng_t>::create_instance(&r0,is);
		}
		return is;
	    }
	    
	    /* Subclass loading method: */

	    std::ostream& save_data(std::ostream& os) const
	    {
		if(theta_grid!=NULL)
		{
		    theta_grid->save(os);
		}
		else
		{
		    os<<'N'<<std::endl;
		}
		return os;
	    }

	protected:

	    /* Updates minimal invariant masses backward: */

	    bool refresh_s_bounds()
	    {
		value_type mmin=std::max(this->m_out(0)+this->m_out(1),this->incoming_channel->m_min_min());
		if(mmin>this->incoming_channel->m_max())
		{
		    return false;
		}
		if(!(this->incoming_channel->set_m_min(mmin)))
		{
		    this->incoming_channel->s_warning();
		    return false;
		}
		return true;
	    }

	private:

	    /* Weight of the s-pair generation: */

	    value_type sweight;

	    /* Rest-frame momentum of decay product: */

	    value_type phat;

	    /* Polar angle adaptive grid: */

	    parni<value_type,1,rng_t>* theta_grid;

	    /* Polar angle random number: */

	    value_type r0;
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const typename model_t::value_type backward_s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::twopi=(typename model_t::value_type)2*std::acos(-(typename model_t::value_type)1);
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t backward_s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::long_dir;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t backward_s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::trans_dir1;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t backward_s_branching<model_t,N_in,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,4> >::trans_dir2;
}

#endif /*CAMGEN_S_BRANCH_BACK_H_*/

