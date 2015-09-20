//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PS_CHANNEL_H_
#define CAMGEN_PS_CHANNEL_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Momentum channel class definition. The momentum channel holds momenta,      *
 * invariant masses and the extremal invariant mass bounds for a particular    *
 * channel in the recursive phase space decomposition. It also contains a list *
 * of particle channels, representing the off-shell invariant mass sampling    *
 * channels within.                                                            *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <list>
#include <limits>
#include <Camgen/bit_string.h>
#include <Camgen/particle.h>
#include <Camgen/flav_comp.h>
#include <Camgen/bipart.h>
#include <Camgen/ps_vol.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class particle_channel;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_branching;

    /* Phase space momentum channel class: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class momentum_channel
    {
	public:

	    enum generation_status
	    {
		reset,
		s_set,
		p_set
	    };

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef std::size_t size_type;
	    typedef typename model_t::value_type value_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef particle_channel<model_t,N_in,N_out,rng_t> particle_channel_type;
	    typedef std::list<particle_channel_type*> particle_channel_list;
	    typedef typename particle_channel_list::iterator particle_channel_iterator;
	    typedef typename particle_channel_list::const_iterator const_particle_channel_iterator;
	    static const std::size_t N_bits=N_in+N_out-1;
	    typedef bit_string<N_bits> bit_string_type;
	    typedef generation_status status_type;

	    /* Static members: */
	    /*-----------------*/

	    /* Particle channel comparator: */

	    static bool particle_comp(const particle_channel_type* first,const particle_channel_type* second)
	    {
		if(first->particle_type==NULL)
		{
		    if(second->particle_type==NULL)
		    {
			return false;
		    }
		    else
		    {
			return true;
		    }
		}
		if(second->particle_type==NULL)
		{
		    return false;
		}
		flavour_comp< particle<model_t> >comp;
		return comp(first->particle_type,second->particle_type);
	    }

	    /* Public data members */
	    /*---------------------*/

	    /* Channel encoding bit string. */

	    const bit_string_type bitstring;

	    /* Public constructors */
	    /*---------------------*/

	    /* Momentum-allocating constructor. */

	    momentum_channel(bit_string_type bitstring_):bitstring(bitstring_),momentum(new momentum_type),alloc_momentum(true),invariant_mass(0),status(reset)
	    {
		if(timelike())
		{
		    sminmin=(value_type)0;
		    smaxmax=std::numeric_limits<value_type>::infinity();
		}
		if(spacelike())
		{
		    sminmin=-std::numeric_limits<value_type>::infinity();
		    smaxmax=(value_type)0;
		}
		momentum->assign((value_type)0);
	    }

	    /* Momentum reference-copying constructor. */

	    momentum_channel(bit_string_type bitstring_,momentum_type* momentum_):bitstring(bitstring_),momentum(momentum_),alloc_momentum(false),invariant_mass(0),status(reset)
	    {
		if(timelike())
		{
		    sminmin=(value_type)0;
		    smaxmax=std::numeric_limits<value_type>::infinity();
		}
		if(spacelike())
		{
		    sminmin=-std::numeric_limits<value_type>::infinity();
		    smaxmax=(value_type)0;
		}
		momentum->assign((value_type)0);
	    }

	    /* Destructor */
	    /*------------*/

	    virtual ~momentum_channel()
	    {
		if(alloc_momentum)
		{
		    delete momentum;
		}
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    delete *it;
		}
	    }

	    /* Public modifyers */
	    /*------------------*/

	    void add_particle_channel(particle_channel_type* channel)
	    {
		if(on_shell())
		{
		    particle_channel_iterator it=particle_channels.begin();
		    if(it!=particle_channels.end())
		    {
			delete *it;
			*it=channel;
		    }
		    else
		    {
			particle_channels.insert(it,channel);
		    }
		}
		else
		{
		    particle_channel_iterator it=std::lower_bound(particle_channels.begin(),particle_channels.end(),channel,particle_comp);

		    if(it!=particle_channels.end() and (*it)->particle_type==channel->particle_type)
		    {
			delete *it;
			*it=channel;
		    }
		    else
		    {
			particle_channels.insert(it,channel);
		    }
		}
	    }

	    /* Evaluates all subchannel weights: */

	    bool evaluate_weight()
	    {
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=(*it)->evaluate_weight();
		}
		return q;
	    }

	    /* Returns the particle channel of flavour phi. Returns the null
	     * pointer if not found: */

	    particle_channel_type* find_particle_channel(const particle<model_t>* phi)
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    if((*it)->particle_type==phi)
		    {
			return *it;
		    }
		}
		return NULL;
	    }

	    /* Returns the particle channel of particle name phi. Returns the null
	     * pointer if not found: */

	    particle_channel_type* find_particle_channel(const std::string& phi)
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    if((*it)->particle_type->get_name()==phi)
		    {
			return *it;
		    }
		}
		return NULL;
	    }

	    /* Evaluates the invariant mass from the momentum. */

	    value_type evaluate_s()
	    {
		invariant_mass=spacetime_type::dot(*momentum,*momentum);
		return invariant_mass;
	    }

	    /* Sets the momentum pointer: */

	    void set_momentum(momentum_type* p_)
	    {
		if(alloc_momentum)
		{
		    delete momentum;
		}
		momentum=p_;
		alloc_momentum=false;
	    }

	    /* Sets the minimal minimal invariant mass-squared and updates all particle
	     * channels: */
	    
	    bool set_s_min_min(const value_type& sminmin_)
	    {
		sminmin=sminmin_;
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_min_min(sminmin));
		}
		return q;
	    }

	    /* Sets the maximal invariant mass-squared for all channels: */
	    
	    bool set_s_max_max(const value_type& smaxmax_)
	    {
		smaxmax=smaxmax_;
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_max_max(smaxmax));
		}
		return q;
	    }

	    /* Sets the minimal minimal invariant mass and updates all particle
	     * channels: */
	    
	    bool set_m_min_min(const value_type& mminmin_)
	    {
		return set_s_min_min(sgn_sq(mminmin_));
	    }

	    /* Sets the maximal invariant mass for all channels: */
	    
	    bool set_m_max_max(const value_type& mmaxmax_)
	    {
		return set_s_max_max(sgn_sq(mmaxmax_));
	    }

	    /* Sets the minimal invariant mass-squared for all particle
	     * channels: */
	    
	    bool set_s_min(const value_type& smin_)
	    {
		bool q=true;
		value_type sm=std::max(sminmin,smin_);
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_min(sm));
		}
		return q;
	    }

	    /* Sets the maximal invariant mass-squared for all particle
	     * channels: */
	    
	    bool set_s_max(const value_type& smax_)
	    {
		bool q=true;
		value_type sm=std::min(smaxmax,smax_);
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_max(sm));
		}
		return q;
	    }

	    /* Sets the minimal invariant mass for all particle channels: */
	    
	    bool set_m_min(const value_type& mmin_)
	    {
		return set_s_min(sgn_sq(mmin_));
	    }

	    /* Sets the maximal invariant mass for all channels: */
	    
	    bool set_m_max(const value_type& mmax_)
	    {
		return set_s_max(sgn_sq(mmax_));
	    }

	    /* Function filling the mmin components for timelike channels: */

	    void assign_mmin_components(value_type mminout[N_out][N_out])
	    {
		mmin_decomps.clear();
		mmin_decomps_bar.clear();
		std::set<size_type>s,sbar;
		for(size_type i=0;i<N_out;++i)
		{
		    if(bitstring[N_in-1+i])
		    {
			s.insert(i);
		    }
		    else
		    {
			sbar.insert(i);
		    }
		}
		bi_partition part(s);
		std::vector<const value_type*>row;
		for(size_type i=0;i<part.size();++i)
		{
		    row.clear();
		    for(size_type j=0;j<part[i].size();++j)
		    {
			row.push_back(&mminout[part[i][j].first][part[i][j].second]);
		    }
		    mmin_decomps.push_back(row);
		}
		bi_partition partbar(sbar);
		for(size_type i=0;i<partbar.size();++i)
		{
		    row.clear();
		    for(size_type j=0;j<partbar[i].size();++j)
		    {
			row.push_back(&mminout[partbar[i][j].first][partbar[i][j].second]);
		    }
		    mmin_decomps_bar.push_back(row);
		}
	    }

	    /* Re-evaluates mminmin from branchings. */

	    bool refresh_s_min(const value_type& s,const value_type& sa,const value_type& sb)
	    {
		if(on_shell())
		{
		    return true;
		}
		if(timelike())
		{
		    value_type result(0);
		    for(size_type i=0;i<mmin_decomps.size();++i)
		    {
			value_type m(0);
			for(size_type j=0;j<mmin_decomps[i].size();++j)
			{
			    m+=(*mmin_decomps[i][j]);
			}
			result=std::max(result,m);
		    }
		    return set_s_min_min(result*result);
		}
		else if(spacelike())
		{
		    value_type s1(0);
		    for(size_type i=0;i<mmin_decomps.size();++i)
		    {
			value_type s(0);
			for(size_type j=0;j<mmin_decomps[i].size();++j)
			{
			    s+=(*mmin_decomps[i][j]);
			}
			s1=std::max(s1,s);
		    }
		    value_type s2(0);
		    for(size_type i=0;i<mmin_decomps_bar.size();++i)
		    {
			value_type s(0);
			for(size_type j=0;j<mmin_decomps_bar[i].size();++j)
			{
			    s+=(*mmin_decomps_bar[i][j]);
			}
			s2=std::max(s2,s);
		    }
		    value_type lambda=std::sqrt(Kallen(s,sa,sb)*Kallen(s,s1,s2));

		    value_type t1min=sb+s1-0.5*((s+sb-sa)*(s+s1-s2)+lambda)/s;
		    value_type t1max=sb+s1-0.5*((s+sb-sa)*(s+s1-s2)-lambda)/s;

		    value_type t2min=sa+s2-0.5*((s+sa-sb)*(s+s2-s1)+lambda)/s;
		    value_type t2max=sa+s2-0.5*((s+sa-sb)*(s+s2-s1)-lambda)/s;

		    bool q1,q2;
		    if(t1max<t2max)
		    {
			q1=set_s_min_min(t1min);
			q2=set_s_max_max(t1max);
		    }
		    else
		    {
			q1=set_s_min_min(t2min);
			q2=set_s_max_max(t2max);
		    }
		    return (q1 and q2);

		    return true;
		}
		return false;
	    }

	    /* Refreshes internal parameters: */

	    void refresh_params()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->refresh_params();
		}
	    }

	    /* Set the status to s-generated: */

	    void set_status_s_generated()
	    {
		if(status==s_set)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invariant mass regenerated..."<<endlog;
		}
		if(status==p_set)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invariant mass regenerated after momentum was set..."<<endlog;
		}
		status=s_set;
	    }

	    /* Set the status to p-generated: */

	    void set_status_p_generated()
	    {
		if(status==p_set)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"momentum regenerated after momentum was set..."<<endlog;
		}
		status=p_set;
	    }

	    /* Reset the status: */

	    void reset_status()
	    {
		status=reset;
	    }


	    /* Public readout functions */
	    /*--------------------------*/

	    /* Returns the starting iterator to the particle channels: */

	    particle_channel_iterator begin_particle_channels()
	    {
		return particle_channels.begin();
	    }

	    /* Returns the end iterator of the particle channels: */

	    particle_channel_iterator end_particle_channels()
	    {
		return particle_channels.end();
	    }

	    /* Returns the starting const iterator to the particle channels: */

	    const_particle_channel_iterator begin_particle_channels() const
	    {
		return particle_channels.begin();
	    }

	    /* Returns the end const iterator of the particle channels: */

	    const_particle_channel_iterator end_particle_channels() const
	    {
		return particle_channels.end();
	    }

	    /* Returns a reference to the momentum flowing through the channel.
	     * */

	    momentum_type& p()
	    {
		return *momentum;
	    }

	    /* Returns a const reference to the momentum flowing through the
	     * channel. */

	    const momentum_type& p() const
	    {
		return *momentum;
	    }

	    /* Returns a reference to the i-th momentum component. */

	    value_type& p(size_type i)
	    {
		return (*momentum)[i];
	    }

	    /* Returns a const reference to the i-th momentum component. */

	    const value_type& p(size_type i) const
	    {
		return (*momentum)[i];
	    }

	    /* Returns a reference to the invariant mass-squared. */

	    value_type& s()
	    {
		return invariant_mass;
	    }

	    /* Returns a const reference to the invariant mass-squared. */

	    const value_type& s() const
	    {
		return invariant_mass;
	    }

	    /* Returns a reference to the signed invariant mass. */

	    value_type m() const
	    {
		return sgn_sqrt(invariant_mass);
	    }

	    /* Returns a const reference to the minimal invariant mass-squared.
	     * */

	    const value_type& s_min_min() const
	    {
		return sminmin;
	    }

	    /* Returns a const reference to the minimal signed invariant mass.
	     * */

	    value_type m_min_min() const
	    {
		return sgn_sqrt(sminmin);
	    }

	    /* Returns a const reference to the minimal invariant mass-squared.
	     * */

	    const value_type& s_max_max() const
	    {
		return smaxmax;
	    }

	    /* Returns a const reference to the minimal signed invariant mass.
	     * */

	    value_type m_max_max() const
	    {
		return sgn_sqrt(smaxmax);
	    }
	    
	    /* Returns whether there are no particles propagating: */

	    bool empty() const
	    {
		return (particle_channels.size()==0);
	    }

	    /* Returns whther there are particles propagating: */

	    bool any() const
	    {
		return (particle_channels.size()!=0);
	    }

	    /* Returns the number of particles propagating: */

	    size_type count_particle_channels()
	    {
		return particle_channels.size();
	    }

	    /* Returns whether the channel is external. */

	    bool on_shell() const
	    {
		size_type c=bitstring.count();
		return (c==0 or c==1 or c==(N_in+N_out-1));
	    }

	    /* Returns whether the channel is internal. */

	    bool off_shell() const
	    {
		return !on_shell();
	    }

	    /* Returns whether we have a t-type channel.*/

	    bool spacelike() const
	    {
		return (N_in==1)?false:(bitstring.count()<N_bits and bitstring.count()!=1 and bitstring[0]);
	    }

	    /* Returns whether we have an s-type channel. */

	    bool timelike() const
	    {
		return !spacelike();
	    }

	    /* Returns whether this is the s-hat channel.*/

	    bool shat_channel() const
	    {
		return (timelike() and bitstring.count()==N_out);
	    }

	    /* Returns the generation status: */

	    status_type get_status() const
	    {
		return status;
	    }

	    /* Printing method. */
	    
	    std::ostream& print(std::ostream& os) const
	    {
		std::stringstream sstr1,sstr2;
		sstr1<<bitstring;
		int n=N_in+N_out+10;
		os<<std::endl<<std::setw(N_in+N_out+10)<<std::left<<sstr1.str();
		n+=6*model_t::dimension+15;
		sstr2<<*momentum;
		os<<std::setw(6*model_t::dimension+20)<<std::left<<sstr2.str();
		os<<std::setw(23)<<std::left<<m();
		os<<std::setw(23)<<std::left<<m_min_min();
		os<<std::setw(23)<<std::left<<m_max_max()<<std::endl;
		os<<std::setw(n+70)<<std::setfill('.')<<'.'<<std::endl;
		os<<std::setfill(' ');
		os<<std::endl;
		for(const_particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->print(os);
		    os<<std::endl;
		}
		os<<std::setw(n+70)<<std::setfill('-')<<'-'<<std::endl;
		os<<std::setfill(' ');
		return os;
	    }

	    /* Other printing method: */

	    std::ostream& shortprint(std::ostream& os) const
	    {
		os<<bitstring<<"\t, p = "<<*momentum<<", m = "<<m()<<", [m--,m++] = ["<<m_min_min()<<','<<m_max_max()<<']';
		return os;
	    }

	private:

	    /* Private data members: */
	    /*-----------------------*/

	    /* Particles propagating in the momentum channel: */

	    particle_channel_list particle_channels;

	    /* Momentum: */

	    momentum_type* momentum;

	    /* Momentum allocator flag: */

	    bool alloc_momentum;

	    /* Invariant mass-squared: */

	    value_type invariant_mass;

	    /* Extremal invariant mass bounds: */

	    value_type sminmin,smaxmax;

	    /* Vector of references for computing smin: */

	    std::vector<std::vector<const value_type*> > mmin_decomps;

	    /* Vector of references for computing the adjoint smin: */

	    std::vector<std::vector<const value_type*> > mmin_decomps_bar;

	    /* Channel generation status: */

	    status_type status;

	    /* Private modifiers: */
	    /*--------------------*/

    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t momentum_channel<model_t,N_in,N_out,rng_t>::N_bits;

    /* Output stream overload: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>std::ostream& operator << (std::ostream& os,const momentum_channel<model_t,N_in,N_out,rng_t>& c)
    {
	return c.shortprint(os);
    }

    /* Phase space momentum channel comparison class, ordering according to the
     * Caravaglios-Moretti algorithm, */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class momentum_channel_comp
    {
	public:

	    typedef momentum_channel<model_t,N_in,N_out,rng_t> value_type;

	    bool operator ()(const value_type& a,const value_type& b) const
	    {
		return a.bitstring<b.bitstring;
	    }
	    bool operator ()(const value_type* a,const value_type* b) const
	    {
		return a->bitstring<b->bitstring;
	    }
    };
}

#endif /*CAMGEN_PS_CHANNEL_H_*/

