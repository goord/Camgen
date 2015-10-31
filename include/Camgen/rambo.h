//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file rambo.h
    \brief Uniform phase space generator class using the RAMBO algorithm.
 */

#ifndef CAMGEN_RAMBO_H_
#define CAMGEN_RAMBO_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Phase space generator class using the RAMBO algorithm to generate uniformly *
 * distributed momenta. In the case of a massive final state, the massless     *
 * momenta are rescaled by a factor which is computed with the Newton-Raphson  *
 * method.                                                                     *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/Minkowski.h>
#include <Camgen/MC_config.h>
#include <Camgen/Lorentz.h>
#include <Camgen/rn_strm.h>
#include <Camgen/uni_sphere.h>
#include <Camgen/ps_gen.h>
#include <Camgen/ps_vol.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t=typename model_t::spacetime_type>class rambo;

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class rambo<model_t,N_in,N_out,rng_t,Minkowski_type>: public ps_generator<model_t,N_in,N_out>
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

	    static const size_type D=model_t::dimension;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Default constructor: */

	    rambo(init_state_type* is_):base_type(is_),spheregen(new uniform_sphere<value_type,D-2,rng_t>(&phat)){}

	    /* Copy constructor: */

	    rambo(const rambo<model_t,N_in,N_out,rng_t>& other):base_type(other),spheregen(new uniform_sphere<value_type,D-2,rng_t>(&phat)){}

	    /* Destructor: */
	    
	    ~rambo()
	    {
		delete spheregen;
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generates the final state. */

	    bool generate_fs()
	    {
		momentum_type P=this->P_in();
		value_type stot=spacetime_type::dot(P,P);
		if(stot<(value_type)0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"negative total invariant mass-squared encountered: P = "<<P<<endlog;
		    this->fs_weight()=(value_type)0;
		    return false;
		}
		value_type Etot=std::sqrt(stot);
		if(Etot<this->ps_generator_base<model_type>::M_out_sum())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"total energy "<<Etot<<" too small to create final state with M = "<<this->ps_generator_base<model_type>::M_out_sum()<<endlog;
		    this->fs_weight()=(value_type)0;
		    return false;
		}

		/* CM-frame massless momenta generation: */

		momentum_type k;
		k.assign((value_type)0);
		for(size_type i=0;i<N_out;++i)
		{
		    value_type rho=rn_stream::throw_number();
		    for(size_type n=0;n<D-3;++n)
		    {
			rho*=rn_stream::throw_number();
		    }
		    phat=-std::log(rho);
		    (this->p_out(i))[0]=phat;
		    spheregen->generate();
		    for(size_type mu=1;mu<D;++mu)
		    {
			this->p_out(i)[mu]=spheregen->object()[mu-1];
		    }
		    k+=(this->p_out(i));
		}

		/* Rescaling: */

		value_type ksq=spacetime_type::dot(k,k);
		value_type scale_factor=Etot/ksq;
		value_type u=std::sqrt(std::abs(ksq));

		value_type kq,w;
		for(size_type i=0;i<N_out;++i)
		{
		    kq=spacetime_type::dot(k,this->p_out(i));
		    value_type& e=this->p_out(i)[0];
		    w=(kq+u*e)/(u+k[0]);
		    e=scale_factor*kq;
		    for(size_type mu=1;mu<D;++mu)
		    {
			this->p_out(i)[mu]=scale_factor*(u*(this->p_out(i)[mu])-w*k[mu]);
		    }
		}
		if(!this->massive_fs())
		{
		    this->fs_weight()=massless_ps<value_type,N_out,D>::volume(Etot);
		}
		else
		{
		    /* Massive case, with extra constraction: */

		    value_type xi=(value_type)1;
		    value_type f,df,d;
		    vector<value_type,N_out>Esq;
		    for(size_type i=0;i<N_out;++i)
		    {
			Esq[i]=(this->p_out(i)[0])*(this->p_out(i)[0]);
		    }
		    for(size_type i=0;i<NR_iterations();++i)
		    {
			f=-Etot;
			df=(value_type)0;
			for(size_type j=0;j<N_out;++j)
			{
			    d=std::sqrt(xi*xi*Esq[j]+this->M2_out(j));
			    f+=d;
			    df+=Esq[j]/d;
			}
			xi-=f/(xi*df);
		    }
		    for(size_type i=0;i<N_out;++i)
		    {
			this->p_out(i)[0]=std::sqrt(xi*xi*Esq[i]+this->M2_out(i));
			for(size_type mu=1;mu<D;++mu)
			{
			    this->p_out(i)[mu]*=xi;
			}
		    }

		    /* Extra weight due to phase space contraction: */

		    (this->fs_weight())=massless_ps<value_type,N_out,D>::volume(Etot)*std::pow(xi,int((D-2)*N_out-D+1))*Etot;
		    value_type factor=(value_type)1;
		    value_type denom=(value_type)0;
		    for(size_type i=0;i<N_out;++i)
		    {
			value_type e=this->p_out(i)[0];
			value_type psq=e*e-this->M2_out(i);
			factor*=(std::sqrt(psq)/e);
			denom+=(psq/e);
		    }
		    (this->fs_weight())*=(factor/denom);
		}
		
		/* Boosting CM-momenta to the lab frame: */

		for(size_type i=0;i<N_out;++i)
		{
		    boost_from_restframe(this->p_out(i),P);
		}

		return true;
	    }

	    /* Evaluates the weight for the given final state. */

	    bool evaluate_fs_weight()
	    {
		value_type stot=this->s_hat();
		if(stot<(value_type)0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"negative total invariant mass-squared encountered"<<endlog;
		    this->fs_weight()=(value_type)0;
		    return false;
		}
		value_type Etot=std::sqrt(stot);
		if(Etot<this->ps_generator_base<model_type>::M_out_sum())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"total energy too small to create final state"<<endlog;
		    this->fs_weight()=(value_type)0;
		    return false;
		}
		this->fs_weight()=massless_ps<value_type,N_out,D>::volume(Etot);
		if(this->massive_fs())
		{
		    value_type factor=(value_type)1;
		    value_type denom=(value_type)0;
		    value_type num=(value_type)0;
		    momentum_type P=this->P_in();
		    for(size_type i=0;i<N_out;++i)
		    {
			momentum_type p=this->p_out(i);
			boost_to_restframe(p,P,Etot);
			value_type E=p[0];
			value_type pvec=std::sqrt(E*E-this->M2_out(i));
			factor*=(pvec/E);
			num+=pvec;
			denom+=(pvec*pvec/E);
		    }
		    this->fs_weight()*=(Etot*std::pow(num/Etot,int((D-2)*N_out-D+1))*factor/denom);
		}
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method implementation. */

	    virtual rambo<model_t,N_in,N_out,rng_t,Minkowski_type>* clone() const
	    {
		return new rambo<model_t,N_in,N_out,rng_t,Minkowski_type>(*this);
	    }

	    std::string type() const
	    {
		return "rambo";
	    }

	private:

	    /* Private data: */

	    value_type phat;
	    uniform_sphere<value_type,D-2,rng_t>* spheregen;
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const typename rambo<model_t,N_in,N_out,rng_t,Minkowski_type>::size_type rambo<model_t,N_in,N_out,rng_t,Minkowski_type>::D;
}

#endif /*CAMGEN_RAMBO_H_*/

