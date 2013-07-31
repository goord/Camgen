//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_gen_base.h
    \brief Abstract interface base class for phase space generators.
 */

#ifndef PS_GEN_BASE_H_
#define PS_GEN_BASE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for phase space generators providing a common interface   *
 * to use in phase space cut objects and output interfaces. The class requires   *
 * subclasses to implement const access to incoming and outgoing momenta and the *
 * hadronic and partonic invariant masses. It provides some kinematic functions. *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/vector.h>
#include <Camgen/type_holders.h>
#include <Camgen/MC_int.h>

namespace Camgen
{
    /// Abstract interface for phase space generators, used by output interfaces
    /// and cut objects.
    
    template<class model_t>class ps_generator_base
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;

	    /* Static const integers: */

	    static const size_type k0=spacetime_type::timelike_direction;
	    static const size_type kL=(model_type::beam_direction<0)?(-model_type::beam_direction):(model_type::beam_direction);
	    static const size_type kT1=(kL==1)?2:1;
	    static const size_type kT2=(kL==3)?2:3;
	    
	    /* Static utility methods: */
	    /*-------------------------*/

	    /// Contracts the argument Lorentz vectors.

	    static value_type dot(const momentum_type& v1,const momentum_type& v2)
	    {
		return spacetime_type::dot(v1,v2);
	    }

	    /// Takes the (Euclidean) inner product of the argument vectors.

	    static value_type space_dot(const momentum_type& v1,const momentum_type& v2)
	    {
		return spacetime_type::space_dot(v1,v2);
	    }

	    /// Takes the invariant mass-squared of the argument vector.

	    static value_type s(const momentum_type& v)
	    {
		return dot(v,v);
	    }

	    /// Takes the invariant mass-squared of the sum of the argument vectors.

	    static value_type s(const momentum_type& v1,const momentum_type& v2)
	    {
		return s(v1+v2);
	    }

	    /// Takes the (signed) invariant mass of the argument vector.

	    static value_type m(const momentum_type& v)
	    {
		value_type sq(s(v));
		return (sq<(value_type)0)?(-std::sqrt(-sq)):(std::sqrt(sq));
	    }

	    /// Takes the (signed) invariant mass of the sum of the argument vectors.

	    static value_type m(const momentum_type& v1,const momentum_type& v2)
	    {
		return m(v1+v2);
	    }

	    /// Takes the energy component of the argument vector.

	    static const value_type& E(const momentum_type& v)
	    {
		return v[k0];
	    }

	    /// Projects out the energy component of the argument vector.

	    static momentum_type vec(const momentum_type& v)
	    {
		momentum_type result(v);
		result[kL]=(value_type)0;
		return v;
	    }

	    /// Returns the Euclidean momentum-squared of the argument vector.

	    static value_type P2(const momentum_type& v)
	    {
		return space_dot(v,v);
	    }

	    /// Returns the Euclidean momentum of the argument vector.

	    static value_type P(const momentum_type& v)
	    {
		return std::sqrt(space_dot(v,v));
	    }

	    /// Projects out the energy and longitudinal component of the argument vector.

	    static momentum_type vecT(const momentum_type& v)
	    {
		momentum_type result(v);
		result[k0]=(value_type)0;
		result[kL]=(value_type)0;
		return v;
	    }

	    /// Returns the transverse momentum-squared of the argument vector.
	    
	    static value_type pT2(const momentum_type& v)
	    {
		momentum_type result(v);
		result[kL]=(value_type)0;
		return P2(result);
	    }

	    /// Returns the transverse momentum of the argument vector.
	    
	    static value_type pT(const momentum_type& v)
	    {
		momentum_type result(v);
		result[kL]=(value_type)0;
		return P(result);
	    }

	    /// Returns the rapidity of the argument vector.

	    static value_type y(const momentum_type& v)
	    {
		return (value_type)0.5*std::log((v[k0]+v[kL])/(v[k0]-v[kL]));
	    }

	    /// Returns the rapidity difference of the agument vectors.

	    static value_type d_y(const momentum_type& v1,const momentum_type& v2)
	    {
		return (value_type)0.5*std::log((v1[k0]+v1[kL])*(v2[k0]-v2[kL])/((v1[k0]-v1[kL])*(v2[k0]+v2[kL])));
	    }

	    /// Returns the pseudo-rapidity of the argument vector.

	    static value_type eta(const momentum_type& v)
	    {
		value_type q(P(v));
		return (value_type)0.5*std::log((q+v[kL])/(q-v[kL]));
	    }

	    /// Returns the pseudo-rapidity difference of the agument vectors.

	    static value_type d_eta(const momentum_type& v1,const momentum_type& v2)
	    {
		value_type q1(P(v1));
		value_type q2(P(v2));
		return (value_type)0.5*std::log((q1+v1[kL])*(q2-v2[kL])/((q1-v1[kL])*(q2+v2[kL])));
	    }

	    /// Returns the cosine of the angle w.r.t. the beam axis of the argument
	    /// vector.

	    static value_type cos_theta(const momentum_type& v)
	    {
		return v[kL]/P(v);
	    }

	    /// Returns the polar angle w.r.t. the beam axis of the argument vector.

	    static value_type theta(const momentum_type& v)
	    {
		return std::acos(v[kL]/P(v));
	    }

	    /// Returns the azimuthal angle w.r.t. the beam axis of the argument vector.

	    static value_type phi(const momentum_type& v)
	    {
		return std::atan2(v[kT1],v[kT2]);
	    }

	    /// Returns the separation in eta-phi space between the argument vectors.

	    static value_type d_R(const momentum_type& v1,const momentum_type& v2)
	    {
		value_type dphi=phi(v1)-phi(v2);
		value_type deta=d_eta(v1,v2);
		return std::sqrt(deta*deta+dphi*dphi);
	    }

	    /// Returns the cosine of the angle between the argument momenta.

	    static value_type cos_alpha(const momentum_type& v1,const momentum_type& v2)
	    {
		return space_dot(v1,v2)/(P(v1)*P(v2));
	    }

	    /// Returns the angle between the argument momenta.

	    static value_type alpha(const momentum_type& v1,const momentum_type& v2)
	    {
		return std::acos(space_dot(v1,v2)/(P(v1)*P(v2)));
	    }

	    /* Default constructor: */

	    ps_generator_base(){}

	    /* Destructor: */

	    virtual ~ps_generator_base(){}

	    /// Generates new event according to the given strategy argument:

	    virtual bool next_event(int strategy=0)=0;

	    /// Abstract method returning the number of incoming particles.

	    virtual size_type n_in() const=0;

	    /// Abstract method returning the number of outgoing particles.

	    virtual size_type n_out() const=0;

	    /// Virtual method returning the number of subprocesses.

	    virtual size_type processes() const
	    {
		return 1;
	    }

	    /// Virtual method returning the subprocess counter.

	    virtual size_type process_id() const
	    {
		return 0;
	    }

	    /// Virtual method returning the subprocess id of the i-th subprocess.

	    virtual size_type process_id(size_type i) const
	    {
		return 0;
	    }

	    /// Abstract method returning the total cross section.

	    virtual MC_integral<value_type> xsec() const=0;

	    /// Virtual method returning the i-th subprocess cross section.

	    virtual MC_integral<value_type> xsec(size_type i) const
	    {
		if(i==process_id())
		{
		    return xsec();
		}
		MC_integral<value_type>result;
		return result;
	    }

	    /// Virtual method returning the current subprocess cross section.

	    virtual MC_integral<value_type> sub_xsec() const
	    {
		return xsec();
	    }

	    /// Abstract method reading out an incoming momentum.

	    virtual const momentum_type& p_in(size_type) const=0;

	    /// Abstract method reading out an outgoing momentum.
	    
	    virtual const momentum_type& p_out(size_type) const=0;

	    /// Abstract method returning the beam energies.

	    virtual value_type beam_energy(int) const=0;

	    /// Abstract method returning the total (hadronic) invariant mass.

	    virtual value_type Ecm() const=0;

	    /// Abstract method reading out an incoming mass.

	    virtual value_type M_in(size_type) const=0;

	    /// Abstract method reading out an outgoing mass.
	    
	    virtual value_type M_out(size_type) const=0;

	    /// Virtual method returning the i-th incoming invariant mass-squared.

	    virtual value_type s_in(size_type i) const
	    {
		return spacetime_type::dot(p_in(i),p_in(i));
	    }

	    /// Virtual method returning the i-th outgoing invariant mass-squared.

	    virtual value_type s_out(size_type i) const
	    {
		return spacetime_type::dot(p_out(i),p_out(i));
	    }

	    /// Virtual method returning the i-th incoming invariant mass.

	    virtual value_type m_in(size_type i) const
	    {
		return std::sqrt(s_in(i));
	    }

	    /// Virtual method returning the i-th outgoing invariant mass.

	    virtual value_type m_out(size_type i) const
	    {
		return std::sqrt(s_out(i));
	    }

	    /// Virtual method returning the i-th incoming particle id (returns
	    /// zero by default).

	    virtual int id_in(size_type i) const
	    {
		return 0;
	    }

	    /// Virtual method returning the i-th outgoing particle id (returns
	    /// zero by default).

	    virtual int id_out(size_type i) const
	    {
		return 0;
	    }

	    /// Virtual method returning the total (hadronic) invariant mass-squared.

	    virtual value_type s_tot() const
	    {
		return Ecm()*Ecm();
	    }

	    /// Virtual method returning the effective (partonic) invariant
	    /// mass.

	    virtual value_type Ecm_hat() const
	    {
		momentum_type q(P_in());
		return m(q);
	    }

	    /// Virtual method returning the effective (partonic) invariant mass-squared.

	    virtual value_type s_hat() const
	    {
		momentum_type q(P_in());
		return s(q);
	    }

	    /// Virtual method returning the total event weight (returns one by default).

	    virtual value_type w() const
	    {
		return (value_type)1;
	    }

	    /// Virtual method returning the maximal event weight (returns one by default).

	    virtual value_type max_w() const
	    {
		return (value_type)1;
	    }

	    /// Virtual method returning the event weight without flux and
	    /// symmetry factors etc. The first argument denotes whether
	    /// helicities are included, the second whether colours are
	    /// (returns one by default).

	    virtual value_type ps_weights(bool with_hels,bool with_cols) const
	    {
		return (value_type)1;
	    }

	    /// Virtual method returning the flux, symmetry and conversion
	    /// factors. The first argument denotes whther helicities are
	    /// included, the second whether colours are (returns one by
	    /// default).

	    virtual value_type ps_factors(bool with_hels,bool with_cols) const
	    {
		return (value_type)1;
	    }

	    /// Returns the total incoming momentum.

	    virtual momentum_type P_in() const
	    {
		momentum_type q(p_in(0));
		size_type n(n_in());
		for(size_type i=1;i<n;++i)
		{
		    q+=p_in(i);
		}
		return q;
	    }

	    /// Returns the total outgoing momentum.

	    virtual momentum_type P_out() const
	    {
		momentum_type q(p_out(0));
		size_type n(n_out());
		for(size_type i=1;i<n;++i)
		{
		    q+=p_out(i);
		}
		return q;
	    }

	    /// Returns the mu-th component of the total incoming momentum.

	    virtual value_type P_in(size_type mu) const
	    {
		value_type q(p_in(0,mu));
		size_type n(n_in());
		for(size_type i=1;i<n;++i)
		{
		    q+=p_in(i,mu);
		}
		return q;
	    }

	    /// Returns the mu-th component of the total outgoing momentum.

	    virtual value_type P_out(size_type mu) const
	    {
		value_type q(p_out(0,mu));
		size_type n(n_out());
		for(size_type i=1;i<n;++i)
		{
		    q+=p_out(i,mu);
		}
		return q;
	    }

	    /// Returns the total partonic incoming energy.

	    virtual value_type E_in() const
	    {
		return P_in(k0);
	    }

	    /// Returns the total partonic outgoing energy.

	    virtual value_type E_out() const
	    {
		return P_out(k0);
	    }

	    /// Returns the sum of masses of incoming particles.

	    virtual value_type M_in() const
	    {
		value_type q(m_in(0));
		size_type n(n_in());
		for(size_type i=1;i<n;++i)
		{
		    q+=m_in(i);
		}
		return q;
	    }

	    /// Returns the sum of masses of outgoing particles.

	    virtual value_type M_out() const
	    {
		value_type q(m_out(0));
		size_type n(n_out());
		for(size_type i=1;i<n;++i)
		{
		    q+=m_out(i);
		}
		return q;
	    }

	    /// Returns the i-th momentum, where i<0 means incoming momenta and
	    /// i>0 outgoing ones. If i==0,i<-N_in or i>N_out, an error will
	    /// occur.

	    virtual const momentum_type& p(int i) const
	    {
		return (i<0)?p_in(-i-1):p_out(i-1);
	    }

	    /// Returns the i-th mass, where i<0 means incoming masses and
	    /// i>0 outgoing ones. If i==0,i<-N_in or i>N_out, an error will
	    /// occur.

	    virtual value_type M(int i) const
	    {
		return (i<0)?M_in(-i-1):M_out(i-1);
	    }
	   

	    /// Reads the i-th incoming momentum's mu-th component.
	    
	    virtual const value_type& p_in(size_type i,size_type mu) const
	    {
		return p_in(i)[mu];
	    }

	    /// Reads the i-th outgoing momentum's mu-th component.
	    
	    virtual const value_type& p_out(size_type i,size_type mu) const
	    {
		return p_out(i)[mu];
	    }

	    /// Returns the i-th momentum's mu-th component.

	    virtual const value_type& p(int i,size_type mu) const
	    {
		return (i<0)?p_in(-i-1,mu):p_out(i-1,mu);
	    }

	    /// Method returning the i-th invariant mass-squared.

	    virtual value_type s(int i) const
	    {
		return (i<0)?s_in(-i-1):s_out(i-1);
	    }

	    /// Method returning the i-th invariant mass.

	    virtual value_type m(int i) const
	    {
		return (i<0)?m_in(-i-1):m_out(i-1);
	    }

	    /// Virtual method returning the i-th incoming beam id.

	    virtual int beam_id(int i) const
	    {
		return 0;
	    }

	    /// Returns the cernlib pdf group number for the i-th incoming beam.

	    virtual int pdfg(int i) const
	    {
		return -1;
	    }

	    /// Returns the cernlib pdf set number for the i-th incoming beam.

	    virtual int pdfs(int i) const
	    {
		return -1;
	    }

	    /// Virtual method returning the factorisation scale.

	    virtual value_type mu_F() const
	    {
		return 0;
	    }

	    /// Method returning the i-th particle id.

	    virtual int id(int i) const
	    {
		return (i<0)?id_in(-i-1):id_out(i-1);
	    }

	    /// Method determining the colour connection for the event.

	    virtual void fill_colours(std::vector<int>& c,std::vector<int>& cbar) const
	    {
		return;
	    }

	    /// Contracts the argument external momenta.

	    value_type dot(int i1,int i2) const
	    {
		return dot(p(i1),p(i2));
	    }

	    /// Takes the (Euclidean) inner product of the argument external momenta.

	    value_type space_dot(int i1,int i2) const
	    {
		return space_dot(p(i1),p(i2));
	    }

	    /// Takes the invariant mass-squared of the sum of the argument external
	    /// momenta.

	    value_type s(int i,int j) const
	    {
		return s(p(i),p(j));
	    }

	    /// Takes the (signed) invariant mass of the sum of the argument external
	    /// momenta.

	    value_type m(int i,int j) const
	    {
		return m(p(i),p(j));
	    }

	    /// Takes the energy component of the argument external momentum.

	    const value_type& E(int i) const
	    {
		return p(i,k0);
	    }

	    /// Projects out the energy component of the argument external momentum.

	    momentum_type vec(int i) const
	    {
		return vec(p(i));
	    }

	    /// Returns the Euclidean momentum-squared of the argument external momentum.

	    value_type P2(int i) const
	    {
		return P2(p(i));
	    }

	    /// Returns the Euclidean momentum of the argument external momentum.

	    value_type P(int i) const
	    {
		return P(p(i));
	    }

	    /// Projects out the energy and longitudinal component of the argument
	    /// external momentum.

	    momentum_type vecT(int i) const
	    {
		return vecT(p(i));
	    }

	    /// Returns the transverse momentum-squared of the argument external momentum.
	    
	    value_type pT2(int i) const
	    {
		return pT2(p(i));
	    }

	    /// Returns the transverse momentum of the argument external momentum.
	    
	    value_type pT(int i) const
	    {
		return pT(p(i));
	    }

	    /// Returns the rapidity of the argument external momentum.

	    value_type y(int i) const
	    {
		return y(p(i));
	    }

	    /// Returns the rapidity difference of the agument external momenta.

	    value_type d_y(int i1,int i2) const
	    {
		return d_y(p(i1),p(i2));
	    }

	    /// Returns the pseudo-rapidity of the argument external momentum.

	    value_type eta(int i) const
	    {
		return eta(p(i));
	    }

	    /// Returns the pseudo-rapidity difference of the agument external momenta.

	    value_type d_eta(int i1,int i2)
	    {
		return d_eta(p(i1),p(i2));
	    }

	    /// Returns the cosine of the angle w.r.t. the beam axis of the argument
	    /// external momentum.

	    value_type cos_theta(int i) const
	    {
		return cos_theta(p(i));
	    }

	    /// Returns the polar angle w.r.t. the beam axis of the argument external
	    /// momentum.

	    value_type theta(int i) const
	    {
		return theta(p(i));
	    }

	    /// Returns the azimuthal angle w.r.t. the beam axis of the argument external
	    /// momentum.

	    value_type phi(int i) const
	    {
		return phi(p(i));
	    }

	    /// Returns the separation in eta-phi space between the argument external momenta.

	    value_type d_R(int i1,int i2) const
	    {
		return d_R(p(i1),p(i2));
	    }

	    /// Returns the cosine of the angle between the argument external momenta.

	    value_type cos_alpha(int i1,int i2) const
	    {
		return cos_alpha(p(i1),p(i2));
	    }

	    /// Returns the angle between the argument external momenta.

	    value_type alpha(int i1,int i2) const
	    {
		return alpha(p(i1),p(i2));
	    }
    };
    template<class model_t>const typename ps_generator_base<model_t>::size_type ps_generator_base<model_t>::k0;
    template<class model_t>const typename ps_generator_base<model_t>::size_type ps_generator_base<model_t>::kL;
    template<class model_t>const typename ps_generator_base<model_t>::size_type ps_generator_base<model_t>::kT1;
    template<class model_t>const typename ps_generator_base<model_t>::size_type ps_generator_base<model_t>::kT2;

    /// Viewer class for phase space generators.

    template<class model_t>class ps_generator_viewer: public ps_generator_base<model_t>
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;

	    /* Public constructors: */
	    /*----------------------*/

	    /* Default constructor. */

	    ps_generator_viewer():gen(NULL),lock(false){}

	    /* Constructor with generator instance. */

	    ps_generator_viewer(const ps_generator_base<model_t>* gen_):gen(gen_),lock(false){}

	    /* Public destructors: */
	    /*---------------------*/

	    /* Destructor: */

	    virtual ~ps_generator_viewer(){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generates new event according to the given strategy argument: */

	    bool next_event(int strategy)
	    {
		return false;
	    }

	    /* Sets the generator instance. */

	    bool set_generator(const ps_generator_base<model_t>* gen_)
	    {
		if(lock)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"request to set generator instance discarded by lock"<<endlog;
		    return false;
		}
		gen=gen_;
		return (gen!=NULL);
	    }

	    /* Sets the generator, foor good: */

	    bool lock_generator(const ps_generator_base<model_t>* gen_)
	    {
		if(lock)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"request to set generator instance discarded by lock"<<endlog;
		    return false;
		}
		if(gen_!=NULL)
		{
		    gen=gen_;
		    lock=true;
		    return true;
		}
		return false;
	    }

	    bool locked() const
	    {
		return lock;
	    }

	    /* Public readout functions: */
	    /*---------------------------*/

	    /* Returns whether the state of the generator is ready to read out
	    kinematic variables. */
	    
	    bool valid() const
	    {
		return (gen!=NULL);
	    }

	    /* Returns the number of incoming particles. */

	    size_type n_in() const
	    {
		return gen->n_in();
	    }

	    /* Returns the number of outgoing particles. */

	    size_type n_out() const
	    {
		return gen->n_out();
	    }

	    /* Returning the number of subprocesses. */

	    size_type processes() const
	    {
		return gen->processes();
	    }

	    /* Returns the subprocess counter. */

	    size_type process_id() const
	    {
		return gen->process_id();
	    }

	    /* Returns the subprocess id of the i-th subprocess. */

	    size_type process_id(size_type i) const
	    {
		return gen->process_id(i);
	    }

	    /* Returns the total cross section. */

	    MC_integral<value_type> xsec() const
	    {
		return gen->xsec();
	    }

	    /* Returns the i-th subprocess cross section. */

	    MC_integral<value_type> xsec(size_type i) const
	    {
		return gen->xsec(i);
	    }

	    /* Reads out an incoming momentum. */

	    const momentum_type& p_in(size_type i) const
	    {
		return gen->p_in(i);
	    }

	    /* Reads out an outgoing momentum. */
	    
	    const momentum_type& p_out(size_type i) const
	    {
		return gen->p_out(i);
	    }

	    /* Reads out an incoming momentum. */

	    value_type M_in(size_type i) const
	    {
		return gen->M_in(i);
	    }

	    /* Reads out an outgoing momentum. */
	    
	    value_type M_out(size_type i) const
	    {
		return gen->M_out(i);
	    }

	    /* Abstract method returning the beam energies. */

	    value_type beam_energy(int i) const
	    {
		return gen->beam_energy(i);
	    }

	    /* Returns the total (hadronic) invariant mass. */

	    value_type Ecm() const
	    {
		return gen->Ecm();
	    }

	    /* Returns the i-th incoming invariant mass-squared. */

	    value_type s_in(size_type i) const
	    {
		return gen->s_in(i);
	    }

	    /* Returns the i-th outgoing invariant mass-squared. */

	    value_type s_out(size_type i) const
	    {
		return gen->s_out(i);
	    }

	    /* Returns the i-th incoming invariant mass. */

	    value_type m_in(size_type i) const
	    {
		return gen->m_in(i);
	    }

	    /* Returns the i-th outgoing invariant mass. */

	    value_type m_out(size_type i) const
	    {
		return gen->m_out(i);
	    }

	    /* Virtual method returning the i-th incoming beam id. */

	    int beam_id(int i) const
	    {
		return gen->beam_id(i);
	    }

	    /// Returns the cernlib pdf group number for the i-th incoming beam.

	    int pdfg(int i) const
	    {
		return gen->pdfg(i);
	    }

	    /// Returns the cernlib pdf set number for the i-th incoming beam.

	    int pdfs(int i) const
	    {
		return gen->pdfs(i);
	    }

	    /// Returns the factorisation scale.

	    value_type mu_F() const
	    {
		return gen->mu_F();
	    }

	    /* Returns the i-th incoming particle id. */

	    int id_in(size_type i) const
	    {
		return gen->id_in(i);
	    }

	    /* Returns the i-th outgoing particle id. */

	    int id_out(size_type i) const
	    {
		return gen->id_out(i);
	    }

	    /// Method determining the colour connection for the event.

	    virtual void fill_colours(std::vector<int>& c,std::vector<int>& cbar) const
	    {
		gen->fill_colours(c,cbar);
	    }

	    /* Returns the total (hadronic) invariant mass-squared. */

	    value_type s_tot() const
	    {
		return gen->s_tot();
	    }

	    /* Returns the effective (partonic) invariant mass. */

	    value_type Ecm_hat() const
	    {
		return gen->Ecm_hat();
	    }

	    /* Returns the effective (partonic) invariant mass-squared. */

	    value_type s_hat() const
	    {
		return gen->s_hat();
	    }

	    /* Returns the total event weight. */

	    virtual value_type w() const
	    {
		return gen->w();
	    }

	    /* Returns the maximal event weight. */

	    virtual value_type max_w() const
	    {
		return gen->max_w();
	    }

	    /* Returns the partial event weights. */

	    value_type ps_weights(bool with_hels,bool with_cols) const
	    {
		return gen->ps_weights(with_hels,with_cols);
	    }

	    /* Returns the partial event factors. */

	    value_type ps_factors(bool with_hels,bool with_cols) const
	    {
		return gen->ps_factors(with_hels,with_cols);
	    }

	protected:

	    const ps_generator_base<model_t>* gen;

	private:

	    bool lock;
    };

    template<class model_t>class test_generator: public ps_generator_base<model_t>
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;

	    const size_type N_in;
	    const size_type N_out;

	    /* Constructor: */

	    test_generator(size_type N_in_,size_type N_out_):N_in(N_in_),N_out(N_out_),dummy_s(0),dummy_n(0),validity_flag(true)
	    {
		dummy_p.assign(0);
	    }

	    /* Returns the validity flag: */

	    bool valid() const
	    {
		return validity_flag;
	    }

	    /* Generates new event according to the given strategy argument: */

	    bool next_event(int strategy)
	    {
		return true;
	    }

	    /* Returns the nr incoming of particles. */

	    size_type n_in() const
	    {
		return N_in;
	    }

	    /* Returns the nr of outgoing particles. */

	    size_type n_out() const
	    {
		return N_out;
	    }

	    /* Returns the total cross section. */

	    MC_integral<value_type> xsec() const
	    {
		MC_integral<value_type>result;
		return result;
	    }

	    /* Returns the i-th incoming momentum. */

	    const momentum_type& p_in(size_type i) const
	    {
		if(i>=N_in)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested incoming momentum "<<i<<" out of range for "<<N_in<<" incoming particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_p;
	    }

	    /* Returns the i-th outgoing momentum. */
	    
	    const momentum_type& p_out(size_type i) const
	    {
		if(i>=N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested outgoing momentum "<<i<<" out of range for "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_p;
	    }

	    /* Returns the i-th incoming mass. */

	    value_type M_in(size_type i) const
	    {
		if(i>=N_in)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested incoming mass "<<i<<" out of range for "<<N_in<<" incoming particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th outgoing momentum. */
	    
	    value_type M_out(size_type i) const
	    {
		if(i>=N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested outgoing mass "<<i<<" out of range for "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th incoming momentum's mu-th component. */
	    
	    const value_type& p_in(size_type i,size_type mu) const
	    {
		if(mu<dummy_p.size())
		{
		    return p_in(i)[mu];
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested momentum component "<<mu<<" out of range for "<<dummy_p.size()<<" spacetime dimensions"<<endlog;
		validity_flag=false;
		return dummy_s;
	    }

	    /* Returns the i-th outgoing momentum's mu-th component. */
	    
	    const value_type& p_out(size_type i,size_type mu) const
	    {
		if(mu<dummy_p.size())
		{
		    return p_out(i)[mu];
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested momentum component "<<mu<<" out of range for "<<dummy_p.size()<<" spacetime dimensions"<<endlog;
		validity_flag=false;
		return dummy_s;
	    }

	    /* Returns the total (hadronic) invariant mass. */

	    value_type Ecm() const
	    {
		return dummy_s;
	    }

	    /* Returns the beam energies. */

	    value_type beam_energy(int i) const
	    {
		return 0.5*dummy_s;
	    }

	    /* Returns the i-th incoming invariant mass-squared. */

	    value_type s_in(size_type i) const
	    {
		if(i>=N_in)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested incoming inv. mass^2 "<<i<<" out of range for "<<N_in<<" incoming particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th outgoing invariant mass-squared. */

	    value_type s_out(size_type i) const
	    {
		if(i>=N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested outgoing inv. mass^2 "<<i<<" out of range for "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th incoming invariant mass. */

	    value_type m_in(size_type i) const
	    {
		if(i>=N_in)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested incoming inv. mass "<<i<<" out of range for "<<N_in<<" incoming particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th outgoing invariant mass. */

	    value_type m_out(size_type i) const
	    {
		if(i>=N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested outgoing inv. mass "<<i<<" out of range for "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th incoming particle id. */

	    int id_in(size_type i) const
	    {
		if(i>=N_in)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested incoming id "<<i<<" out of range for "<<N_in<<" incoming particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_n;
	    }

	    /* Returns the i-th outgoing particle id. */

	    int id_out(size_type i) const
	    {
		if(i>=N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested outgoing id "<<i<<" out of range for "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_n;
	    }

	    /* Returns the i-th momentum */

	    const momentum_type& p(int i) const
	    {
		if(i<-(int)N_in or i==0 or i<N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested momentum "<<i<<" out of range for "<<N_in<<" incoming and "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_p;
	    }

	    /* Returns the i-th mass */

	    value_type M(int i) const
	    {
		if(i<-(int)N_in or i==0 or i<N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested mass "<<i<<" out of range for "<<N_in<<" incoming and "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return 0;
	    }

	    /* Returns the i-th momentum's mu-th component. */

	    const value_type& p(int i,size_type mu) const
	    {
		if(mu<dummy_p.size())
		{
		    return p(i)[mu];
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested component "<<mu<<" out of range in "<<dummy_p.size()<<" spacetime dimensions"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th invariant mass-squared. */

	    value_type s(int i) const
	    {
		if(i<-(int)N_in or i==0 or i<N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Requested inv. mass^2 "<<i<<" out of range for "<<N_in<<" incoming and "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th invariant mass. */

	    value_type m(int i) const
	    {
		if(i<-(int)N_in or i==0 or i<N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested inv. mass "<<i<<" out of range for "<<N_in<<" incoming and "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_s;
	    }

	    /* Returns the i-th particle id. */

	    int id(int i) const
	    {
		if(i<-(int)N_in or i==0 or i<N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requested id "<<i<<" out of range for "<<N_in<<" incoming and "<<N_out<<" outgoing particles"<<endlog;
		    validity_flag=false;
		}
		return dummy_n;
	    }

	private:

	    momentum_type dummy_p;
	    value_type dummy_s;
	    int dummy_n;
	    mutable bool validity_flag;
    };
}

#endif /*PS_GEN_BASE_H_*/

