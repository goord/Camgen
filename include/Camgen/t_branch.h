//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_T_BRANCH_H_
#define CAMGEN_T_BRANCH_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Phase space branching module splitting a (possibly) spacelike momentum, with  *
 * the help of the second incoming momentum into a spacelike and 2 timelike      *
 * momenta, one of which is considered to be auxiliary.                          *
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
    template<class model_t,std::size_t N_out,class rng_t,class spacetime_t>class t_branching;

    /* t-branching specialisation for Minkowski spacetimes: */
    
    template<class model_t,std::size_t N_out,class rng_t>class t_branching<model_t,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,model_t::dimension> >: public ps_branching<model_t,2,N_out,rng_t>
    {
	typedef ps_branching<model_t,2,N_out,rng_t> base_type;
        
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

	    /* Second beam: */
	    
	    const channel_type* const sb_channel;

	    /* Total output channel: */

	    channel_type* const stot_channel;

	    /* Flag denoting whether second outgoing branching is timelike (final t-channel): */

	    const bool is_final;

	    /* Flag denoting whether the branching should generate s2_channel: */

	    /* Constructor. First argument is the incoming spacelike channel or
	     * first beam, second argument is the second beam, third argument is
	     * the first outgoing timelike channel, third argument is the
	     * spacelike outgoing channel and the last argument denotes the
	     * second timelike channel. The last boolean denotes whether the
	     * branching should allow more t-channels to be chained. */

	    t_branching(channel_type* in1_,
		        const channel_type* in2_,
			channel_type* out1_,
			channel_type* out2_,
			channel_type* out3_,
			channel_type* out13_,
			bool is_final_):base_type(in1_,out1_,is_final_?out3_:out2_),
					sb_channel(in2_),
					s1_channel(out1_),
					t_channel(out2_),
					s2_channel(out3_),
					stot_channel(out13_),
					is_final(is_final_),
					ssgen(s_pair_generator_factory<value_type,rng_t>::create(this->s_generator(out1_),this->s_generator(out3_),stot_channel->s())){}

	    /* Destructor: */

	    ~t_branching()
	    {
		delete ssgen;
	    }

	    bool s1_sampling() const
	    {
		return !this->backward_s_sampling or (s1_channel->off_shell() and s1_channel->branching_count()==0);
	    }

	    bool s2_sampling() const
	    {
		return !this->backward_s_sampling or (s2_channel->off_shell() and s2_channel->branching_count()==0);
	    }

	    bool stot_sampling() const
	    {
		return this->backward_s_sampling and this->shat_sampling and stot_channel->shat_channel();
	    }

	    /* Generates positive invariants: */

	    bool generate_s()
	    {
		sweight=(value_type)1;

		if(this->backward_s_sampling)
		{
		    bool q1=s1_sampling();
		    bool q2=s2_sampling();
		    bool q12=stot_sampling();

		    if(q12)
		    {
			value_type mmin1=q1?s1_channel->m():s1_channel->m_min();
			value_type mmin2=q2?s2_channel->m():s2_channel->m_min();
			stot_channel->set_m_min(mmin1+mmin2);
			if(!stot_channel->generate_s())
			{
			    sweight=(value_type)0;
			    this->branching_weight=(value_type)0;
			    return false;
			}
			stot_channel->set_status_s_generated();
			sweight*=stot_channel->s_weight();
		    }
		    else if(stot_channel->get_status()==ps_channel_type::reset)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"sampling invariant mass without total invariant mass generated..."<<endlog;
		    }
		    if(q1 and q2)
		    {
			log(log_level::error)<<CAMGEN_STREAMLOC<<"sampling 2 auxiliary invariant masses in t-branching is not supported yet"<<endlog;
		    }
		    if(q1 or q2)
		    {
			channel_type* q1_channel=q1?s1_channel:s2_channel;
			channel_type* q2_channel=q2?s2_channel:s1_channel;
			q1_channel->set_m_max(stot_channel->m()-q2_channel->m());
			if(!q1_channel->generate_s())
			{
			    sweight=(value_type)0;
			    this->branching_weight=(value_type)0;
			    return false;
			}
			q1_channel->set_status_s_generated();
			sweight*=q1_channel->s_weight();
		    }
		}
		else
		{
		    if(!ssgen->generate())
		    {
			sweight=(value_type)0;
			this->branching_weight=(value_type)0;
			return false;
		    }
		    s1_channel->set_status_s_generated();
		    s2_channel->set_status_s_generated();
		    sweight=ssgen->weight();
		}
		return true;
	    }

	    /* Generates the spacelike invariant: */

	    bool generate_t()
	    {
		value_type sa=this->s_in();
		value_type sb=sb_channel->s();
		value_type s1=s1_channel->s();
		value_type s2=s2_channel->s();
		value_type s=stot_channel->s();

		value_type A=0.5*(sa+sb+s1+s2-s-(sa-sb)*(s1-s2)/s);
		value_type l12=Kallen(s,s1,s2);
		if(l12<(value_type)0)
		{
		    tweight=(value_type)0;
		    this->branching_weight=(value_type)0;
		    return false;
		}
		value_type B=0.5*std::sqrt(Kallen(s,sa,sb)*l12)/s;
		value_type tmin=A-B;
		value_type tmax=A+B;

		if(!(t_channel->set_s_range(tmin,tmax)))
		{
		    tweight=(value_type)0;
		    this->branching_weight=(value_type)0;
		    return false;
		}
		if(!t_channel->generate_s())
		{
		    tweight=(value_type)0;
		    this->branching_weight=(value_type)0;
		    return false;
		}
		t_channel->set_status_s_generated();
		tweight=t_channel->s_weight();
		return (tweight>(value_type)0);
	    }

	    /* Generates momenta with current invariant masses: */

	    bool generate_p()
	    {
		value_type s1=s1_channel->s();
		value_type s2=s2_channel->s();
		value_type s=stot_channel->s();

		value_type m=stot_channel->m();
		momentum_type pahat=this->p_in();
		momentum_type Ptot=pahat+sb_channel->p();

		boost_to_restframe(pahat,Ptot,m);
		sphere_generator.generate();
		value_type pdots(0);
		for(size_type i=0;i<model_t::dimension-1;++i)
		{
		    pdots+=(pahat[i+1]*sphere_generator.object()[i]);
		}
		value_type pahatabs=std::sqrt(spacetime_type::space_dot(pahat,pahat));
		value_type k12=Kallen(s,s1,s2);
		if(k12<(value_type)0)
		{
		    this->branching_weight=(value_type)0;
		    return false;
		}
		value_type l12=std::sqrt(k12);
		value_type qhatabs=(value_type)0.5*l12/m;

		value_type A=0.5*(t_channel->s_max()+t_channel->s_min());
		value_type B=0.5*(t_channel->s_max()-t_channel->s_min());
		value_type cos_theta=(t_channel->s()-A)/B;
		value_type sin_theta=std::sqrt((1-cos_theta)*(1+cos_theta));

		value_type c1,c2;
		if(equals(pdots,(value_type)0))
		{
		    c1=qhatabs*cos_theta/pahatabs;
		    c2=qhatabs*sin_theta;
		}
		else
		{
		    value_type alpha=pahatabs/pdots;
		    value_type y=sin_theta/std::sqrt((alpha-(value_type)1)*(alpha+(value_type)1));
		    c1=qhatabs*(cos_theta+y)/pahatabs;
		    c2=qhatabs*alpha*y;
		}
		(s1_channel->p())[0]=(value_type)0.5*(s+s1-s2)/m;


		for(size_type i=1;i<model_t::dimension;++i)
		{
		    (s1_channel->p())[i]=c1*pahat[i]-c2*sphere_generator.object()[i-1];
		}
		boost_from_restframe(s1_channel->p(),Ptot,m);
		s1_channel->set_status_p_generated();
		s1_channel->evaluate_s();
		s2_channel->p()=Ptot-s1_channel->p();
		s2_channel->set_status_p_generated();
		t_channel->p()=this->p_in()-s1_channel->p();
		t_channel->set_status_p_generated();

		this->branching_weight=massless_ps<value_type,2,model_t::dimension>::volume(m)*std::pow(l12/s,int(model_t::dimension-4))*sweight*tweight/((value_type)2*m*pahatabs);
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

		value_type sa=this->s_in();
		value_type sb=sb_channel->s();
		value_type s1=s1_channel->s();

		value_type s=stot_channel->s();
		value_type m=stot_channel->m();
		value_type s2=s2_channel->s();

		value_type A=0.5*(sa+sb+s1+s2-s-(sa-sb)*(s1-s2)/s);
		value_type lab=std::sqrt(Kallen(s,sa,sb));
		value_type l12=std::sqrt(Kallen(s,s1,s2));
		value_type B=0.5*lab*l12/s;
		value_type tmin=A-B;
		value_type tmax=A+B;

		if(!t_channel->set_s_range(tmin,tmax))
		{
		    tweight=(value_type)0;
		    this->branching_weight=(value_type)0;
		    return false;
		}
		if(!t_channel->evaluate_s_weight())
		{
		    tweight=(value_type)0;
		    this->branching_weight=(value_type)0;
		    return false;
		}
		tweight=t_channel->s_weight();
		this->branching_weight=massless_ps<value_type,2,model_t::dimension>::volume(m)*std::pow(l12/s,int(model_t::dimension-4))*sweight*tweight/lab;
		return (this->branching_weight>(value_type)0);
            }

	    /* Returns whether the branchings are equivalent: */

	    bool equiv(const ps_branching<model_t,2,N_out,rng_t>* other) const
	    {
		const t_branching<model_t,N_out,rng_t,spacetime_type>* othercast=dynamic_cast<const t_branching<model_t,N_out,rng_t,spacetime_type>*>(other);
		if(othercast!=NULL)
		{
		    return (this->same_incoming_channel(other) and sb_channel==othercast->sb_channel 
			    				       and s1_channel==othercast->s1_channel
							       and t_channel==othercast->t_channel
							       and s2_channel==othercast->s2_channel);
		}
		return false;
	    }

	    bool t_type() const
	    {
		return true;
	    }

	    /* Returns the t-type for printing: */

	    std::string type() const
	    {
		return is_final?"t.":"t";
	    }

	    const channel_type* t_channel_out() const
	    {
		return t_channel;
	    }

	    const channel_type* s2_channel_out() const
	    {
		return s2_channel;
	    }

	    t_branching<model_type,N_out,rn_engine,spacetime_type>* copy_to_final_branching(channel_type* s2_channel_) const
	    {
		return new t_branching<model_type,N_out,rn_engine,spacetime_type>(this->incoming_channel,sb_channel,s1_channel,t_channel,s2_channel_,stot_channel,true);
	    }

	    std::ostream& print(std::ostream& os) const
	    {
		this->base_type::print(os);
		os<<std::setw(15)<<std::left<<(is_final?t_channel->name:s2_channel->name);
		os<<std::setw(15)<<std::left<<stot_channel->name;
		return os;
	    }

	    
	    momentum_type evaluate_p_in() const
	    {
		return s1_channel->p()+s2_channel->p()-sb_channel->p();
	    }

	private:

	    /* Channels pointer copies: */

	    channel_type* s1_channel;
	    channel_type* t_channel;
	    channel_type* s2_channel;

	    /* Timelike and spacelike weights: */

	    value_type sweight,tweight;

	    /* Invariant mass pair generator: */

	    s_pair_generator_type* ssgen;

	    /* Azimuthal angle generator: */

	    sphere_generator_type sphere_generator;

	    bool evaluate_s_weight()
	    {
		sweight=(value_type)1;
		if(this->backward_s_sampling)
		{
		    bool q1=s1_sampling();
		    bool q2=s2_sampling();
		    bool q12=stot_sampling();

		    if(q12)
		    {
			value_type mmin1=q1?s1_channel->m():s1_channel->m_min();
			value_type mmin2=q2?s2_channel->m():s2_channel->m_min();
			stot_channel->set_m_min(mmin1+mmin2);
			if(!stot_channel->evaluate_s_weight())
			{
			    sweight=(value_type)0;
			    return false;
			}
			sweight*=stot_channel->s_weight();
		    }
		    else if(stot_channel->get_status()==ps_channel_type::reset)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"sampling invariant mass without total invariant mass generated..."<<endlog;
		    }
		    if(q1 and q2)
		    {
			log(log_level::error)<<CAMGEN_STREAMLOC<<"sampling 2 auxiliary invariant masses in t-branching is not supported yet"<<endlog;
		    }
		    if(q1 or q2)
		    {
			channel_type* q1_channel=q1?s1_channel:s2_channel;
			channel_type* q2_channel=q2?s1_channel:s2_channel;
			q1_channel->set_m_max(stot_channel->m()-q2_channel->m());
			if(!q1_channel->evaluate_s_weight())
			{
			    sweight=(value_type)0;
			    return false;
			}
			sweight*=q1_channel->s_weight();
		    }
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
}

#endif /*CAMGEN_T_BRANCH_H_*/

