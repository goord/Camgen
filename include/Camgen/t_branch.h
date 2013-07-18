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
#include <Camgen/branching.h>
#include <Camgen/uni_sphere.h>
#include <Camgen/Lorentz.h>
#include <Camgen/ps_vol.h>
#include <Camgen/parni.h>

namespace Camgen
{
    /* t-branching specialisation for Minkowski spacetimes: */
    
    template<class model_t,std::size_t N_out,class rng_t>class t_branching<model_t,N_out,rng_t,typename Minkowski_type::template implementation<typename model_t::value_type,model_t::dimension> >: public ps_branching<model_t,2,N_out,rng_t>
    {
	typedef ps_branching<model_t,2,N_out,rng_t> base_type;
        
	public:

	    /* Type definitions: */

	    typedef typename model_t::value_type value_type;
	    typedef typename Minkowski_type::template implementation<value_type,model_t::dimension> spacetime_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::channel_type channel_type;
	    typedef uniform_sphere<value_type,model_t::dimension-2,rn_engine> sphere_gen_type;

	    /* Initial t-channel branching flag: */

	    const bool init_t_branching;

	    /* Final t-channel branching flag: */

	    const bool final_t_branching;

	    /* Second beam: */
	    
	    const channel_type* aux_in_channel;

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
			bool final=false):base_type(in1_,2),
	    				  init_t_branching(in1_->on_shell()),
					  final_t_branching(final),
					  aux_in_channel(in2_),
					  s1_channel(out1_),
					  t_channel(out2_),
					  s2_channel(out3_),
					  sphere_generator(new sphere_gen_type)
	    {
		/* t-channel goes first...*/

		this->channels[0]=final?out1_:out2_;
		this->channels[1]=final?out3_:out1_;
	    }

	    /* Destructor: */

	    ~t_branching()
	    {
		delete sphere_generator;
	    }

	    /* Generates positive invariants: */

	    bool generate_s()
	    {
		if(!refresh_s_bounds())
		{
		    sweight=(value_type)0;
		    return false;
		}
		const value_type& sqrts1=s1_channel->m();
		const value_type& sqrts2=s2_channel->m();
		if(this->max_s_pairs()>1)
		{
		    std::size_t n=0;
		    bool q,q1,q2;
		    do
		    {
			q1=s1_channel->generate_s();
			q2=s2_channel->generate_s();
			q=(q1 and q2 and (sqrts1+sqrts2<=sqrts));
			++n;
		    }
		    while(!q and (n<this->max_s_pairs()));
		    if(!q)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"no succesful s-pair generated after "<<n<<" throws--returning weight 0"<<endlog;
			sweight=(value_type)0;
			return false;
		    }
		    sweight=integrate(s1_channel,s2_channel,sqrts)*(s1_channel->s_weight())*(s2_channel->s_weight());
		}
		else
		{
		    bool q1=s1_channel->generate_s();
		    bool q2=s2_channel->generate_s();
		    if(!(q1 and q2) or ((sqrts1+sqrts2)>sqrts))
		    {
			sweight=(value_type)0;
			return false;
		    }
		    sweight=(s1_channel->s_weight())*(s2_channel->s_weight());
		}
		return (sweight>(value_type)0);
	    }

	    /* Generates the spacelike invariant: */

	    bool generate_t()
	    {
		value_type sa=this->s_in();
		value_type sb=aux_in_channel->s();
		value_type s1=s1_channel->s();
		value_type s2=s2_channel->s();
		
		value_type A=0.5*(sa+sb+s1+s2-s-(sa-sb)*(s1-s2)/s);
		value_type B=0.5*std::sqrt(Kallen(s,sa,sb)*Kallen(s,s1,s2))/s;
		value_type tmin=A-B;
		value_type tmax=A+B;

		if(!(t_channel->set_s_range(tmin,tmax)))
		{
		    t_channel->s_warning();
		    tweight=(value_type)0;
		    return false;
		}
		if(!t_channel->generate_s())
		{
		    tweight=(value_type)0;
		    return false;
		}
		tweight=t_channel->s_weight();
		cos_theta=(t_channel->s()-A)/B;
		return (tweight>(value_type)0);
	    }

	    /* Generates momenta with current invariant masses: */

	    bool generate_p()
	    {
		momentum_type pahat=this->p_in();
		boost_to_restframe(pahat,Ptot,sqrts);
		sphere_generator->generate();
		value_type pdots(0);
		for(size_type i=0;i<model_t::dimension-1;++i)
		{
		    pdots+=(pahat[i+1]*sphere_generator->object()[i]);
		}
		value_type pahatabs=std::sqrt(spacetime_type::space_dot(pahat,pahat));
		value_type l12=std::sqrt(Kallen(s,s1_channel->s(),s2_channel->s()));
		value_type qhatabs=(value_type)0.5*l12/sqrts;
		value_type c1,c2,sin_theta=std::sqrt((1-cos_theta)*(1+cos_theta));
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
		(s1_channel->p())[0]=(value_type)0.5*(s+s1_channel->s()-s2_channel->s())/sqrts;
		for(size_type i=1;i<model_t::dimension;++i)
		{
		    (s1_channel->p())[i]=c1*pahat[i]-c2*sphere_generator->object()[i-1];
		}
		boost_from_restframe(s1_channel->p(),Ptot,sqrts);
		s2_channel->p()=Ptot-s1_channel->p();
		t_channel->p()=this->p_in()-s1_channel->p();
		this->weight()=massless_ps<value_type,2,model_t::dimension>::volume(sqrts)*std::pow(l12/s,int(model_t::dimension-4))*sweight*tweight/((value_type)2*sqrts*pahatabs);
		return (this->weight()>(value_type)0);
	    }

	    /* Overridden weight evaluation method: */
            
	    bool evaluate_branching_weight()
            {
		if(!refresh_s_bounds())
		{
		    sweight=(value_type)0;
		    tweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		if((s1_channel->m()+s2_channel->m()>sqrts) or !s1_channel->evaluate_s_weight() or !s2_channel->evaluate_s_weight())
		{
		    sweight=(value_type)0;
		    tweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		if(this->max_s_pairs()>1)
		{
		    sweight=integrate(s1_channel,s2_channel,sqrts)*(s1_channel->s_weight())*(s2_channel->s_weight());
		}
		else
		{
		    sweight=(s1_channel->s_weight())*(s2_channel->s_weight());
		}
		value_type sa=this->s_in();
		value_type sb=aux_in_channel->s();
		value_type s1=s1_channel->s();
		value_type s2=s2_channel->s();

		value_type A=0.5*(sa+sb+s1+s2-s-(sa-sb)*(s1-s2)/s);
		value_type lab=std::sqrt(Kallen(s,sa,sb));
		value_type l12=std::sqrt(Kallen(s,s1,s2));
		value_type B=0.5*lab*l12/s;
		value_type tmin=A-B;
		value_type tmax=A+B;
		
		if(!t_channel->set_s_range(tmin,tmax))
		{
		    t_channel->s_warning();
		    tweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!t_channel->evaluate_s_weight())
		{
		    tweight=(value_type)0;
		    this->weight()=(value_type)0;
		    return false;
		}
		tweight=t_channel->s_weight();
		this->weight()=massless_ps<value_type,2,model_t::dimension>::volume(sqrts)*std::pow(l12/s,int(model_t::dimension-4))*sweight*tweight/lab;
		return (this->weight()>(value_type)0);
            }

	    /* Returns whether the branchings are equivalent: */

	    bool equiv(const ps_branching<model_t,2,N_out,rng_t>* other) const
	    {
		const t_branching<model_t,N_out,rng_t,spacetime_type>* othercast=dynamic_cast<const t_branching<model_t,N_out,rng_t,spacetime_type>*>(other);
		if(othercast!=NULL)
		{
		    return (this->same_incoming_channel(other) and aux_in_channel==othercast->aux_in_channel 
			    				       and s1_channel==othercast->s1_channel
							       and t_channel==othercast->t_channel
							       and s2_channel==othercast->s2_channel);
		}
		return false;
	    }

	    /* Returns the t-type for printing: */

	    std::string type() const
	    {
		if(final_t_branching)
		{
		    return "t.";
		}
		return "t";
	    }

	    /* Serialization helper: */

	    std::ostream& save_channels(std::ostream& os) const
	    {
		os<<this->incoming_channel->name<<"\t"<<aux_in_channel->name<<std::endl;
		os<<s1_channel->name<<"\t"<<t_channel->name<<"\t"<<s2_channel->name<<std::endl;
		return os;
	    }

	protected:

	    /* Updates maximal invariant masses: */
	    
	    bool refresh_s_bounds()
	    {
		Ptot=this->p_in()+aux_in_channel->p();
		s=spacetime_type::dot(Ptot,Ptot);
		if(s<(value_type)0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"negative total incoming invariant mass encountered in t-branching"<<endlog;
		    return false;
		}
		sqrts=std::sqrt(s);
		if(sqrts<s1_channel->m_min()+s2_channel->m_min())
		{
		    return false;
		}
		if(!(s1_channel->set_m_max(sqrts-s2_channel->m_min())))
		{
		    s1_channel->s_warning();
		    return false;
		}
		if(!(s2_channel->set_m_max(sqrts-s1_channel->m_min())))
		{
		    s2_channel->s_warning();
		    return false;
		}
		return true;
	    }

	    /* Recursively unsets the redundant flag: */
	    
	    void unset_redundant()
	    {
		this->base_type::unset_redundant();
		if(final_t_branching)
		{
		    t_channel->unset_redundant();
		}
		else
		{
		    s2_channel->unset_redundant();
		}
	    }

	private:

	    /* Channels pointer copies: */

	    channel_type* s1_channel;
	    channel_type* t_channel;
	    channel_type* s2_channel;

	    /* Timelike and spacelike weights: */

	    value_type sweight,tweight;

	    /* Azimuthal angle generator: */

	    sphere_gen_type* sphere_generator;
	    
	    /* Total incoming momentum: */
	    
	    momentum_type Ptot;

	    /* Total incoming invariant mass-squared: */
	    
	    value_type s;

	    /* Total incoming invariant mass: */

	    value_type sqrts;

	    /* Cosine of the polar scattering angle: */

	    value_type cos_theta;
    };
}

#endif /*CAMGEN_T_BRANCH_H_*/

