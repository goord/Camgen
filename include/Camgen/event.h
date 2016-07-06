//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file event.h
    \brief Abstract base class for events in Camgen
 */

#ifndef CAMGEN_EVENT_H_
#define CAMGEN_EVENT_H_

#include <Camgen/sub_proc.h>
#include <Camgen/type_holders.h>
#include <Camgen/MC_integral.h>

namespace Camgen
{
    /// Event interface class in Camgen. Provides access to phase space data.

    template<class model_t,std::size_t N_in,std::size_t N_out>class event
    {
        public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
            typedef particle<model_t> particle_type;
	    
	    static const size_type k0=spacetime_type::timelike_direction;

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
		result[0]=(value_type)0;
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

	    event():sub_proc(NULL),procid(1){}

	    /* Destructor: */

	    virtual ~event()
            {
                if(sub_proc!=NULL)
                {
                    delete sub_proc;
                }
            }

            /// Abstract method returning incoming momenta.

            virtual momentum_type p_in(size_type i) const=0;

            /// Abstract method returning outgoing momenta.

            virtual momentum_type p_out(size_type i) const=0;

            /// Virtual method returning the process id.

            virtual int process_id() const
            {
                return procid;
            }

            /// Virtual method returning the event weight.

            virtual value_type w() const
            {
                return (value_type)1;
            }

            /// Virtual method returning the total cross sesction.

            virtual MC_integral<value_type> xsec() const
            {
                return MC_integral<value_type>(0);
            }

            /// Virtual method returning the subprocess cross sesction.

            virtual MC_integral<value_type> process_xsec() const
            {
                return MC_integral<value_type>(0);
            }

	    /// Returns the i-th momentum, where i<0 means incoming momenta and
	    /// i>0 outgoing ones. If i==0,i<-N_in or i>N_out, an error will
	    /// occur.

	    momentum_type p(int i) const
	    {
		return (i<0)?p_in(-i-1):p_out(i-1);
	    }

	    /// Projects out the energy component of the argument external momentum.

	    momentum_type vec(int i) const
	    {
		return vec(p(i));
	    }

	    /// Reads the i-th incoming momentum's mu-th component.
	    
	    value_type p_in(size_type i,size_type mu) const
	    {
		return p_in(i)[mu];
	    }

	    /// Reads the i-th outgoing momentum's mu-th component.
	    
	    value_type p_out(size_type i,size_type mu) const
	    {
		return p_out(i)[mu];
	    }

	    /// Returns the i-th momentum's mu-th component.

	    value_type p(int i,size_type mu) const
	    {
		return (i<0)?p_in(-i-1,mu):p_out(i-1,mu);
	    }

	    /// Returns the Euclidean momentum-squared of the i-th incoming momentum.

	    value_type P2_in(size_type i) const
	    {
		return P2(p_in(i));
	    }

	    /// Returns the Euclidean momentum-squared of the i-th outgoing momentum.

	    value_type P2_out(size_type i) const
	    {
		return P2(p_out(i));
	    }

	    /// Returns the Euclidean momentum-squared of the argument external momentum.

	    value_type P2(int i) const
	    {
		return P2(p(i));
	    }

	    /// Returns the Euclidean momentum of the i-th incoming momentum.

	    value_type P_in(size_type i) const
	    {
		return P(p_in(i));
	    }

	    /// Returns the Euclidean momentum of the i-th outgoing momentum.

	    value_type P_out(size_type i) const
	    {
		return P(p_out(i));
	    }

	    /// Returns the Euclidean momentum of the argument external momentum.

	    value_type P(int i) const
	    {
		return P(p(i));
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

	    /// Returns the total incoming momentum.

	    momentum_type p_tot_in() const
	    {
		momentum_type q;
                q.assign(0);
		for(size_type i=0;i<N_in;++i)
		{
		    q+=p_in(i);
		}
		return q;
	    }

	    /// Returns the total outgoing momentum.

	    momentum_type p_tot_out() const
	    {
		momentum_type q;
                q.assign(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=p_out(i);
		}
		return q;
	    }

            /// Returns the total incoming momentum.

            momentum_type p_tot() const
            {
                return p_tot_in();
            }

	    /// Returns the mu-th component of the total incoming momentum.

	    value_type p_tot_in(size_type mu) const
	    {
		value_type q(0);
		for(size_type i=0;i<N_in;++i)
		{
		    q+=p_in(i,mu);
		}
		return q;
	    }

	    /// Returns the mu-th component of the total outgoing momentum.

	    value_type p_tot_out(size_type mu) const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=p_out(i,mu);
		}
		return q;
	    }

	    /// Returns the mu-th component of the total incoming momentum.

            momentum_type p_tot(size_type mu) const
            {
                return p_tot_in(mu);
            }

	    /// Returns the energy component of the i-th incoming momentum

	    value_type E_in(size_type i) const
	    {
		return p_in(i,k0);
	    }

	    /// Returns the energy component of the i-th outgoing momentum

	    value_type E_out(size_type i) const
	    {
		return p_out(i,k0);
	    }

	    /// Returns the energy component of the argument external momentum.

	    value_type E(int i) const
	    {
		return p(i,k0);
	    }

	    /// Returns the total partonic incoming energy.

	    value_type E_tot_in() const
	    {
		return p_tot_in(k0);
	    }

	    /// Returns the total partonic outgoing energy.

	    value_type E_tot_out() const
	    {
		return p_tot_out(k0);
	    }

	    /// Returns the total partonic incoming energy.

	    value_type E_tot() const
	    {
		return E_tot_in();
	    }

            /// Returns the i-th incoming invariant mass-squared.

            virtual value_type s_in(size_type i) const
            {
                return s(p_in(i));
            }

            /// Returns the i-th outgoing invariant mass-squared.

            virtual value_type s_out(size_type i) const
            {
                return s(p_out(i));
            }

	    /// Method returning the i-th invariant mass-squared.

	    value_type s(int i) const
	    {
		return s(p(i));
	    }

            /// Returns the total incoming partonic invariant mass.

            virtual value_type s_tot_in() const
            {
                return s(p_tot_in());
            }

            /// Returns the total outgoing partonic invariant mass.

            virtual value_type s_tot_out() const
            {
                return s(p_tot_out());
            }

            /// Returns the total incoming partonic invariant mass (s-hat).

            value_type s_tot() const
            {
                return s_tot_in();
            }

            /// Returns the total incoming center-of-mass energy.

            virtual value_type Ecm_hat() const
            {
                return std::sqrt(s_tot());
            }

	    /// Takes the invariant mass-squared of the sum of the argument external
	    /// momenta.

	    value_type s(int i,int j) const
	    {
		return s(p(i),p(j));
	    }

	    /// Returns the sum of invariant masses squared of incoming particles.

	    value_type s_in_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_in;++i)
		{
		    q+=s_in(i);
		}
		return q;
	    }

	    /// Returns the sum of invariant masses squared of outgoing particles.

	    value_type s_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=s_out(i);
		}
		return q;
	    }

            /// Returns the i-th incoming invariant mass.

            virtual value_type m_in(size_type i) const
            {
                return std::sqrt(s_in(i));
            }

            /// Returns the i-th outgoing invariant mass.

            virtual value_type m_out(size_type i) const
            {
                return std::sqrt(s_out(i));
            }

	    /// Method returning the i-th invariant mass.

	    value_type m(int i) const
	    {
		return (i<0)?m_in(-i-1):m_out(i-1);
	    }

	    /// Takes the (signed) invariant mass of the sum of the argument external
	    /// momenta.

	    value_type m(int i,int j) const
	    {
		return m(p(i),p(j));
	    }

	    /// Returns the sum of invariant masses of incoming particles.

	    value_type m_in_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_in;++i)
		{
		    q+=m_in(i);
		}
		return q;
	    }

	    /// Returns the sum of invariant masses of outgoing particles.

	    value_type m_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=m_out(i);
		}
		return q;
	    }

            /// Returns the i-th incoming particle.

            const particle_type* get_particle_in(size_type i) const
            {
                return sub_proc->get_particle_in(i);
            }

            /// Returns the i-th outgoing particle.

            const particle_type* get_particle_out(size_type i) const
            {
                return sub_proc->get_particle_out(i);
            }

            /// Returns the i-th particle.

            const particle_type* get_particle(int i) const
            {
                return sub_proc->get_particle(i);
            }

            /// Returns the i-th incoming particle pdg id.

            int id_in(size_type i) const
            {
                return sub_proc->id_in(i);
            }

            /// Returns the i-th outgoing particle pdg id.

            int id_out(size_type i) const
            {
                return sub_proc->id_out(i);
            }

            /// Returns the i-th particle pdg id.

            int id(int i) const 
            {
                return sub_proc->id(i);
            }

            /// Returns the i-th incoming particle name.

            const std::string& name_in(size_type i) const
            {
                return sub_proc->name_in(i);
            }

            /// Returns the i-th outgoing particle name.

            const std::string& name_out(size_type i) const
            {
                return sub_proc->name_out(i);
            }

            /// Returns the i-th particle name.

            const std::string& name(int i) const 
            {
                return sub_proc->name(i);
            }

            /// Returns the i-th incoming particle mass.

            value_type M_in(size_type i) const
            {
                return sub_proc->m_in(i);
            }

            /// Returns the i-th outgoing particle mass.

            value_type M_out(size_type i) const
            {
                return sub_proc->m_out(i);
            }

	    /// Returns the i-th mass, where i<0 means incoming masses and
	    /// i>0 outgoing ones. If i==0,i<-N_in or i>N_out, an error will
	    /// occur.

	    value_type M(int i) const
	    {
		return (i<0)?M_in(-i-1):M_out(i-1);
	    }

	    /// Returns the sum of masses of incoming particles.

	    value_type M_in_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_in;++i)
		{
		    q+=M_in(i);
		}
		return q;
	    }

	    /// Returns the sum of masses of outgoing particles.

	    value_type M_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=M_out(i);
		}
		return q;
	    }

            /// Returns the i-th incoming particle mass-squared.

            value_type M2_in(size_type i) const
            {
                return M_in(i)*M_in(i);
            }

            /// Returns the i-th outgoing particle mass-squared.

            value_type M2_out(size_type i) const
            {
                return M_out(i)*M_out(i);
            }

	    /// Returns the i-th mass-squared, where i<0 means incoming masses and
	    /// i>0 outgoing ones. If i==0,i<-N_in or i>N_out, an error will
	    /// occur.

	    value_type M2(int i) const
	    {
		return (i<0)?M2_in(-i-1):M2_out(i-1);
	    }

	    /// Returns the sum of masses squared of incoming particles.

	    value_type M2_in_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_in;++i)
		{
		    q+=M2_in(i);
		}
		return q;
	    }

	    /// Returns the sum of masses squared of outgoing particles.

	    value_type M2_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=M2_out(i);
		}
		return q;
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

            /// Returns the sub-process.

            const sub_process<model_type,N_in,N_out>& get_process() const
            {
                return *sub_proc;
            } 
            
            /* Predefined event validation functions. */
            
	    /// Checks whether momentum is conserved.

	    bool check_p_conservation() const
	    {
		momentum_type Pin=p_tot_in();
		momentum_type Pout=p_tot_out();
		bool q=true;
		for(size_type mu=0;mu!=model_t::dimension;++mu)
		{
		    q&=equals(Pin[mu]/Ecm_hat(),Pout[mu]/Ecm_hat());
		}
		if(!q)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"momentum conservation violation detected "<<Pin<<" not equal to "<<Pout<<endlog;
		}
		return q;
	    }

	    /// Checks whether momenta are on-shell.

	    bool check_p_on_shell() const
	    {
		bool q=true;
		for(size_type i=0;i<N_in;++i)
		{
		    value_type s=s_in(i);
		    if(!equals(s/Ecm_hat(),M2_in(i)/Ecm_hat()))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"incoming momentum "<<i<<": "<<p_in(i)<<" with invariant mass "<<s<<" not equal to "<<M2_in(i)<<" detected"<<endlog;
			q=false;
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    value_type s=s_out(i);
		    if(!equals(s/Ecm_hat(),M2_out(i)/Ecm_hat()))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"outgoing momentum "<<i<<": "<<p_out(i)<<" with mass-squared "<<s<<" not equal to "<<M2_out(i)<<" detected"<<endlog;
			q=false;
		    }
		}
		return q;
	    }

	    /// Checks parton level momentum conservation condition

	    bool check_sufficient_shat() const
	    {
                value_type msum=M_out_sum();
		if(Ecm_hat()<=msum)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"total CM-energy "<<Ecm_hat()<<" insufficient to accomodate outgoing particles with total mass "<<msum<<endlog;
		    return false;
		}
		return true;
	    }

	    /// Checks hadron level momentum conservation condition

	    bool check_sufficient_s() const
	    {
                return check_sufficient_shat();
	    }

        protected:

            const sub_process<model_type,N_in,N_out>* sub_proc;
            int procid;
    };

    /* Helper template class for transverse directions of 2-particle initial states in 3d: */

    template<std::size_t N>class directions_3d;

    template<>class directions_3d<1>
    {
	public:
	    static const std::size_t kT1=2;
	    static const std::size_t kT2=3;
    };
    template<>class directions_3d<2>
    {
	public:
	    static const std::size_t kT1=3;
	    static const std::size_t kT2=1;
    };
    template<>class directions_3d<3>
    {
	public:
	    static const std::size_t kT1=1;
	    static const std::size_t kT2=2;
    };

    /// Event interface specialization for 2-particle initial states.

    template<class model_t,std::size_t N_out>class event<model_t,2,N_out>
    {
        public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
            typedef particle<model_t> particle_type;
	    
	    static const size_type k0=spacetime_type::timelike_direction;
	    static const size_type kL=(model_type::beam_direction<0)?(-model_type::beam_direction):(model_type::beam_direction);
	    static const size_type kT1=directions_3d<kL>::kT1;
	    static const size_type kT2=directions_3d<kL>::kT2;

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
		result[0]=(value_type)0;
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
		return std::atan2(v[kT2],v[kT1]);
	    }

	    /// Returns the separation in eta-phi space between the argument vectors.

	    static value_type d_R(const momentum_type& v1,const momentum_type& v2)
	    {
		value_type dphi=phi(v1)-phi(v2);
		value_type deta=d_eta(v1,v2);
		return std::sqrt(deta*deta+dphi*dphi);
	    }

	    /* Default constructor: */

	    event():sub_proc(NULL),procid(1){}

	    /* Destructor: */

	    virtual ~event()
            {
                if(sub_proc!=NULL)
                {
                    delete sub_proc;
                }
            }

            /// Virtual method returning the process id.

            virtual int process_id() const
            {
                return procid;
            }

            /// Virtual method returning the event weight.

            virtual value_type w() const
            {
                return (value_type)1;
            }

            /// Virtual method returning the process cross sesction.

            virtual MC_integral<value_type> xsec() const
            {
                return MC_integral<value_type>(0);
            }

            /// Abstract method returning incoming momenta.

            virtual momentum_type p_in(size_type i) const=0;

            /// Abstract method returning outgoing momenta.

            virtual momentum_type p_out(size_type i) const=0;

	    /// Returns the i-th momentum, where i<0 means incoming momenta and
	    /// i>0 outgoing ones. If i==0,i<-2 or i>N_out, an error will
	    /// occur.

	    momentum_type p(int i) const
	    {
		return (i<0)?p_in(-i-1):p_out(i-1);
	    }

	    /// Projects out the energy component of the argument external momentum.

	    momentum_type vec(int i) const
	    {
		return vec(p(i));
	    }

	    /// Reads the i-th incoming momentum's mu-th component.
	    
	    value_type p_in(size_type i,size_type mu) const
	    {
		return p_in(i)[mu];
	    }

	    /// Reads the i-th outgoing momentum's mu-th component.
	    
	    value_type p_out(size_type i,size_type mu) const
	    {
		return p_out(i)[mu];
	    }

	    /// Returns the i-th momentum's mu-th component.

	    value_type p(int i,size_type mu) const
	    {
		return (i<0)?p_in(-i-1,mu):p_out(i-1,mu);
	    }

	    /// Returns the Euclidean momentum-squared of the i-th incoming momentum.

	    value_type P2_in(size_type i) const
	    {
		return P2(p_in(i));
	    }

	    /// Returns the Euclidean momentum-squared of the i-th outgoing momentum.

	    value_type P2_out(size_type i) const
	    {
		return P2(p_out(i));
	    }

	    /// Returns the Euclidean momentum-squared of the argument external momentum.

	    value_type P2(int i) const
	    {
		return P2(p(i));
	    }

	    /// Returns the Euclidean momentum of the i-th incoming momentum.

	    value_type P_in(size_type i) const
	    {
		return P(p_in(i));
	    }

	    /// Returns the Euclidean momentum of the i-th outgoing momentum.

	    value_type P_out(size_type i) const
	    {
		return P(p_out(i));
	    }

	    /// Returns the Euclidean momentum of the argument external momentum.

	    value_type P(int i) const
	    {
		return P(p(i));
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

	    /// Returns the total incoming momentum.

	    momentum_type p_tot_in() const
	    {
		return p_in(0)+p_in(1);
	    }

	    /// Returns the total outgoing momentum.

	    momentum_type p_tot_out() const
	    {
		momentum_type q;
                q.assign(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=p_out(i);
		}
		return q;
	    }

            /// Returns the total incoming momentum.

            momentum_type p_tot() const
            {
                return p_tot_in();
            }

	    /// Returns the mu-th component of the total incoming momentum.

	    value_type p_tot_in(size_type mu) const
	    {
		return p_in(0,mu)+p_in(1,mu);
	    }

	    /// Returns the mu-th component of the total outgoing momentum.

	    value_type p_tot_out(size_type mu) const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=p_out(i,mu);
		}
		return q;
	    }

	    /// Returns the mu-th component of the total incoming momentum.

            momentum_type p_tot(size_type mu) const
            {
                return p_tot_in(mu);
            }

	    /// Returns the energy component of the i-th incoming momentum

	    value_type E_in(size_type i) const
	    {
		return p_in(i,k0);
	    }

	    /// Returns the energy component of the i-th outgoing momentum

	    value_type E_out(size_type i) const
	    {
		return p_out(i,k0);
	    }

	    /// Returns the energy component of the argument external momentum.

	    value_type E(int i) const
	    {
		return p(i,k0);
	    }

	    /// Returns the total partonic incoming energy.

	    value_type E_tot_in() const
	    {
		return p_tot_in(k0);
	    }

	    /// Returns the total partonic outgoing energy.

	    value_type E_tot_out() const
	    {
		return p_tot_out(k0);
	    }

	    /// Returns the total partonic incoming energy.

	    value_type E_tot() const
	    {
		return E_tot_in();
	    }

            /// Returns the i-th incoming invariant mass-squared.

            virtual value_type s_in(size_type i) const
            {
                return s(p_in(i));
            }

            /// Returns the i-th outgoing invariant mass-squared.

            virtual value_type s_out(size_type i) const
            {
                return s(p_out(i));
            }

	    /// Method returning the i-th invariant mass-squared.

	    value_type s(int i) const
	    {
		return s(p(i));
	    }

            /// Returns the total incoming partonic invariant mass.

            virtual value_type s_tot_in() const
            {
                return s(p_tot_in());
            }

            /// Returns the total outgoing partonic invariant mass.

            virtual value_type s_tot_out() const
            {
                return s(p_tot_out());
            }

            /// Returns the total partonic invariant mass (s-hat).

            virtual value_type s_hat() const
            {
                return s_tot_in();
            }

            /// Returns the total partonic center-of-mass energy.

            virtual value_type Ecm_hat() const
            {
                return std::sqrt(s_hat());
            }

            /// Returns the total hadronic invariant mass-squared.

            virtual value_type s_beams() const
            {
                return s_hat();
            }

            /// Returns the total hadronic center-of-mass energy.

            virtual value_type Ecm_beams() const
            {
                return std::sqrt(s_beams());
            }

	    /// Takes the invariant mass-squared of the sum of the argument external
	    /// momenta.

	    value_type s(int i,int j) const
	    {
		return s(p(i),p(j));
	    }

	    /// Returns the sum of invariant masses squared of incoming particles.

	    value_type s_in_sum() const
	    {
		return s_in(0)+s_in(1);
	    }

	    /// Returns the sum of invariant masses squared of outgoing particles.

	    value_type s_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=s_out(i);
		}
		return q;
	    }

            /// Returns the i-th incoming invariant mass.

            virtual value_type m_in(size_type i) const
            {
                return std::sqrt(s_in(i));
            }

            /// Returns the i-th outgoing invariant mass.

            virtual value_type m_out(size_type i) const
            {
                return std::sqrt(s_out(i));
            }

	    /// Method returning the i-th invariant mass.

	    value_type m(int i) const
	    {
		return (i<0)?m_in(-i-1):m_out(i-1);
	    }

	    /// Takes the (signed) invariant mass of the sum of the argument external
	    /// momenta.

	    value_type m(int i,int j) const
	    {
		return m(p(i),p(j));
	    }

	    /// Returns the sum of invariant masses of incoming particles.

	    value_type m_in_sum() const
	    {
		return m_in(0)+m_in(1);
	    }

	    /// Returns the sum of invariant masses of outgoing particles.

	    value_type m_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=m_out(i);
		}
		return q;
	    }

            /// Returns the i-th incoming particle.

            const particle_type* get_particle_in(size_type i) const
            {
                return sub_proc->get_particle_in(i);
            }

            /// Returns the i-th outgoing particle.

            const particle_type* get_particle_out(size_type i) const
            {
                return sub_proc->get_particle_out(i);
            }

            /// Returns the i-th particle.

            const particle_type* get_particle(int i) const
            {
                return sub_proc->get_particle(i);
            }

            /// Returns the i-th incoming particle pdg id.

            int id_in(size_type i) const
            {
                return sub_proc->id_in(i);
            }

            /// Returns the i-th outgoing particle pdg id.

            int id_out(size_type i) const
            {
                return sub_proc->id_out(i);
            }

            /// Returns the i-th particle pdg id.

            int id(int i) const 
            {
                return sub_proc->id(i);
            }

            /// Returns the i-th incoming particle name.

            const std::string& name_in(size_type i) const
            {
                return sub_proc->name_in(i);
            }

            /// Returns the i-th outgoing particle name.

            const std::string& name_out(size_type i) const
            {
                return sub_proc->name_out(i);
            }

            /// Returns the i-th particle name.

            const std::string& name(int i) const 
            {
                return sub_proc->name(i);
            }

            /// Returns the i-th incoming particle mass.

            value_type M_in(size_type i) const
            {
                return sub_proc->m_in(i);
            }

            /// Returns the i-th outgoing particle mass.

            value_type M_out(size_type i) const
            {
                return sub_proc->m_out(i);
            }

	    /// Returns the i-th mass, where i<0 means incoming masses and
	    /// i>0 outgoing ones. If i==0,i<-2 or i>N_out, an error will
	    /// occur.

	    value_type M(int i) const
	    {
		return (i<0)?M_in(-i-1):M_out(i-1);
	    }

	    /// Returns the sum of masses of incoming particles.

	    value_type M_in_sum() const
	    {
		return M_in(0)+M_in(1);
	    }

	    /// Returns the sum of masses of outgoing particles.

	    value_type M_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=M_out(i);
		}
		return q;
	    }

            /// Returns the i-th incoming particle mass-squared.

            value_type M2_in(size_type i) const
            {
                return M_in(i)*M_in(i);
            }

            /// Returns the i-th outgoing particle mass-squared.

            value_type M2_out(size_type i) const
            {
                return M_out(i)*M_out(i);
            }

	    /// Returns the i-th mass-squared, where i<0 means incoming masses and
	    /// i>0 outgoing ones. If i==0,i<-2 or i>N_out, an error will
	    /// occur.

	    value_type M2(int i) const
	    {
		return (i<0)?M2_in(-i-1):M2_out(i-1);
	    }

	    /// Returns the sum of masses squared of incoming particles.

	    value_type M2_in_sum() const
	    {
		return M2_in(0)+M2_in(1);
	    }

	    /// Returns the sum of masses squared of outgoing particles.

	    value_type M2_out_sum() const
	    {
		value_type q(0);
		for(size_type i=0;i<N_out;++i)
		{
		    q+=M2_out(i);
		}
		return q;
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

	    /// Projects out the energy and longitudinal component of the i-th incoming momentum.

	    momentum_type vecT_in(size_type i) const
	    {
		return vecT(p_in(i));
	    }

	    /// Projects out the energy and longitudinal component of the i-th outgoing momentum.

	    momentum_type vecT_out(size_type i) const
	    {
		return vecT(p_out(i));
	    }

	    /// Projects out the energy and longitudinal component of the argument
	    /// external momentum.

	    momentum_type vecT(int i) const
	    {
		return vecT(p(i));
	    }

	    /// Returns the transverse momentum-squared of the i-th incoming momentum.
	    
	    value_type pT2_in(size_type i) const
	    {
		return pT2(p_in(i));
	    }

	    /// Returns the transverse momentum-squared of the i-th outgoing momentum.
	    
	    value_type pT2_out(size_type i) const
	    {
		return pT2(p_out(i));
	    }

	    /// Returns the transverse momentum-squared of the argument external momentum.
	    
	    value_type pT2(int i) const
	    {
		return pT2(p(i));
	    }

	    /// Returns the transverse momentum of the i-th incoming momentum.
	    
	    value_type pT_in(int i) const
	    {
		return pT(p_in(i));
	    }

	    /// Returns the transverse momentum of the i-th outgoing momentum.
	    
	    value_type pT_out(int i) const
	    {
		return pT(p_out(i));
	    }

	    /// Returns the transverse momentum of the argument external momentum.
	    
	    value_type pT(int i) const
	    {
		return pT(p(i));
	    }

	    /// Returns the rapidity of the i-th incoming momentum.

	    value_type y_in(size_type i) const
	    {
		return y(p_in(i));
	    }

	    /// Returns the rapidity of the i-th outgoing momentum.

	    value_type y_out(size_type i) const
	    {
		return y(p_out(i));
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

	    /// Returns the pseudo-rapidity of the i-th incoming momentum.

	    value_type eta_in(size_type i) const
	    {
		return eta(p_in(i));
	    }

	    /// Returns the pseudo-rapidity of the i-th outgoing momentum.

	    value_type eta_out(size_type i) const
	    {
		return eta(p_out(i));
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

	    /// Returns the cosine of the angle w.r.t. the beam axis of the i-th incoming momentum.

	    value_type cos_theta_in(size_type i) const
	    {
		return cos_theta(p_in(i));
	    }

	    /// Returns the cosine of the angle w.r.t. the beam axis of the i-th outgoing momentum.

	    value_type cos_theta_out(size_type i) const
	    {
		return cos_theta(p_out(i));
	    }

	    /// Returns the cosine of the angle w.r.t. the beam axis of the argument
	    /// external momentum.

	    value_type cos_theta(int i) const
	    {
		return cos_theta(p(i));
	    }

	    /// Returns the polar angle w.r.t. the beam axis of the i-th incoming momentum.

	    value_type theta_in(size_type i) const
	    {
		return theta(p_in(i));
	    }

	    /// Returns the polar angle w.r.t. the beam axis of the i-th outgoing momentum.

	    value_type theta_out(size_type i) const
	    {
		return theta(p_out(i));
	    }

	    /// Returns the polar angle w.r.t. the beam axis of the argument external
	    /// momentum.

	    value_type theta(int i) const
	    {
		return theta(p(i));
	    }

	    /// Returns the azimuthal angle w.r.t. the beam axis of the i-th incoming momentum.

	    value_type phi_in(size_type i) const
	    {
		return phi(p_in(i));
	    }

	    /// Returns the azimuthal angle w.r.t. the beam axis of the i-th outgoing momentum.

	    value_type phi_out(size_type i) const
	    {
		return phi(p_out(i));
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

            /// Returns the sub-process.

            const sub_process<model_type,2,N_out>& get_process() const
            {
                return *sub_proc;
            } 
            
            /* Predefined event validation functions. */
            
	    /// Checks whether momentum is conserved.

	    bool check_p_conservation() const
	    {
		momentum_type Pin=p_tot_in();
		momentum_type Pout=p_tot_out();
		bool q=true;
		for(size_type mu=0;mu!=model_t::dimension;++mu)
		{
		    q&=equals(Pin[mu]/Ecm_hat(),Pout[mu]/Ecm_hat());
		}
		if(!q)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"momentum conservation violation detected "<<Pin<<" not equal to "<<Pout<<endlog;
		}
		return q;
	    }

	    /// Checks whether momenta are on-shell.

	    bool check_p_on_shell() const
	    {
		bool q=true;
		for(size_type i=0;i<2;++i)
		{
		    value_type s=s_in(i);
		    if(!equals(s/Ecm_hat(),M2_in(i)/Ecm_hat()))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"incoming momentum "<<i<<": "<<p_in(i)<<" with invariant mass "<<s<<" not equal to "<<M2_in(i)<<" detected"<<endlog;
			q=false;
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    value_type s=s_out(i);
		    if(!equals(s/Ecm_hat(),M2_out(i)/Ecm_hat()))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"outgoing momentum "<<i<<": "<<p_out(i)<<" with mass-squared "<<s<<" not equal to "<<M2_out(i)<<" detected"<<endlog;
			q=false;
		    }
		}
		return q;
	    }

	    /// Checks parton level momentum conservation condition

	    bool check_sufficient_shat() const
	    {
                value_type msum=M_out_sum();
		if(Ecm_hat()<=msum)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"total CM-energy "<<Ecm_hat()<<" insufficient to accomodate outgoing particles with total mass "<<msum<<endlog;
		    return false;
		}
		return true;
	    }

	    /// Checks hadron level momentum conservation condition

	    bool check_sufficient_s() const
	    {
                value_type msum=M_out_sum();
		if(Ecm_beams()<=msum)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"total CM-energy "<<Ecm_beams()<<" insufficient to accomodate outgoing particles with total mass "<<msum<<endlog;
		    return false;
		}
		return true;
	    }

        protected:

            const sub_process<model_type,2,N_out>* sub_proc;
            int procid;
    };
}

#endif /*CAMGEN_EVENT_H_*/

