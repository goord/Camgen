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
#include <Camgen/ps_decl.h>
#include <Camgen/ps_vol.h>

namespace Camgen
{
    /* Phase space momentum channel class: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class momentum_channel
    {
	friend class particle_channel<model_t,N_in,N_out,rng_t>;
	friend class ps_factory<model_t,N_in,N_out,rng_t,typename model_t::spacetime_type>;
	friend class ps_tree<model_t,N_in,N_out,rng_t>;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef typename model_t::value_type value_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef std::size_t size_type;
	    typedef particle_channel<model_t,N_in,N_out,rng_t> particle_channel_type;
	    typedef std::list<particle_channel_type*> particle_channel_list;
	    typedef typename particle_channel_list::iterator particle_channel_iterator;
	    typedef typename particle_channel_list::const_iterator const_particle_channel_iterator;
	    static const std::size_t N_bits=N_in+N_out-1;
	    typedef bit_string<N_bits> bit_string_type;

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

	    /* Particle channel redundancy tester: */

	    static bool necessary_particle(const particle_channel_type* p)
	    {
		return !(p->is_redundant());
	    }

	    /* Public data members */
	    /*---------------------*/

	    /* Channel encoding bit string. */

	    const bit_string_type bitstring;

	    /* Public constructors */
	    /*---------------------*/

	    /* Momentum-allocating constructor. */

	    momentum_channel(bit_string_type bitstring_):bitstring(bitstring_),momentum(new momentum_type),alloc_momentum(true),invariant_mass(0),signed_mass(0),generation_flag(false)
	    {
		if(timelike())
		{
		    sminmin=(value_type)0;
		    mminmin=(value_type)0;
		    smaxmax=std::numeric_limits<value_type>::infinity();
		    mmaxmax=std::numeric_limits<value_type>::infinity();
		}
		if(spacelike())
		{
		    sminmin=-std::numeric_limits<value_type>::infinity();
		    mminmin=-std::numeric_limits<value_type>::infinity();
		    smaxmax=(value_type)0;
		    mmaxmax=(value_type)0;
		}
		momentum->assign((value_type)0);
	    }

	    /* Momentum reference-copying constructor. */

	    momentum_channel(bit_string_type bitstring_,momentum_type* momentum_):bitstring(bitstring_),momentum(momentum_),alloc_momentum(false),invariant_mass(0),signed_mass(0),generation_flag(false)
	    {
		if(timelike())
		{
		    sminmin=(value_type)0;
		    mminmin=(value_type)0;
		    smaxmax=std::numeric_limits<value_type>::infinity();
		    mmaxmax=std::numeric_limits<value_type>::infinity();
		}
		if(spacelike())
		{
		    sminmin=-std::numeric_limits<value_type>::infinity();
		    mminmin=-std::numeric_limits<value_type>::infinity();
		    smaxmax=(value_type)0;
		    mmaxmax=(value_type)0;
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

	    /* Adds a particle channel: */

	    particle_channel_type* add_particle_channel(const particle<model_t>* phi,const value_type* Ecmhat=NULL)
	    {
		/* If the channel is external, replace the particle channels: */
		
		if(on_shell())
		{
		    if(any())
		    {
			for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
			{
			    delete *it;
			}
			particle_channels.clear();
		    }
		    particle_channel_type* p=new particle_channel_type(this,phi,Ecmhat);
		    particle_channel_iterator it=std::lower_bound(particle_channels.begin(),particle_channels.end(),p,particle_comp);
		    particle_channels.insert(it,p);
		    return p;
		}

		/* If the particle channel is already constructed, return the
		 * instance: */

		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    if((*it)->particle_type==phi)
		    {
			return *it;
		    }
		}

		/* If the particle is pure gauge, do nothing and return the null
		 * pointer: */

		if(phi!=NULL)
		{
		    if(phi->get_pdg_id()==0 and !phi->is_auxiliary())
		    {
			return NULL;
		    }
		}
		particle_channel_type* p=new particle_channel_type(this,phi,Ecmhat);
		particle_channel_iterator it=std::lower_bound(particle_channels.begin(),particle_channels.end(),p,particle_comp);
		particle_channels.insert(it,p);
		return p;
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
		signed_mass=(invariant_mass<0)?(-std::sqrt(-invariant_mass)):std::sqrt(invariant_mass);
		return invariant_mass;
	    }

	    /* Evaluates the invariant mass from the invariant mass-squared. */

	    value_type evaluate_m()
	    {
		signed_mass=(invariant_mass<0)?(-std::sqrt(-invariant_mass)):std::sqrt(invariant_mass);
		return signed_mass;
	    }

	    /* Resets the generation flags. */

	    void reset_generation_flags()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->reset_generation_flags();
		}
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
		mminmin=sgn_sqrt(sminmin);
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_min_min(sminmin));
		}
		return q;
	    }

	    /* Sets the minimal minimal signed invariant mass. Updates all propagator
	    samplers invariant mass limits. */
	    
	    bool set_m_min_min(const value_type& mminmin_)
	    {
		mminmin=mminmin_;
		sminmin=sgn_sq(mminmin);
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_min_min(sminmin));
		}
		return q;
	    }

	    /* Sets the minimal invariant mass-squared for all particle
	     * channels: */
	    
	    bool set_s_min(const value_type& smin_)
	    {
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_min(smin_));
		}
		return q;
	    }

	    /* Sets the minimal invariant mass for all particle channels: */
	    
	    bool set_m_min(const value_type& mmin_)
	    {
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_m_min(mmin_));
		}
		return q;
	    }

	    /* Sets the maximal invariant mass-squared for all channels: */
	    
	    bool set_s_max_max(const value_type& smaxmax_)
	    {
		smaxmax=smaxmax_;
		mmaxmax=sgn_sqrt(smaxmax);
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_max_max(smaxmax));
		}
		return q;
	    }

	    /* Sets the minimal signed invariant mass. Updates all propagator
	    samplers invariant mass limits. */
	    
	    bool set_m_max_max(const value_type& mmaxmax_)
	    {
		mmaxmax=mmaxmax_;
		smaxmax=sgn_sq(mmaxmax);
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_max_max(smaxmax));
		}
		return q;
	    }

	    /* Sets the maximal invariant mass-squared for all particle
	     * channels: */
	    
	    bool set_s_max(const value_type& smax_)
	    {
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_s_max(smax_));
		}
		return q;
	    }

	    /* Sets the maximal invariant mass for all particle channels: */
	    
	    bool set_m_max(const value_type& mmax_)
	    {
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->set_m_max(mmax_));
		}
		return q;
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

	    bool refresh_m_min(const value_type& s,const value_type& sa,const value_type& sb)
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
			value_type s(0);
			for(size_type j=0;j<mmin_decomps[i].size();++j)
			{
			    s+=(*mmin_decomps[i][j]);
			}
			result=std::max(result,s);
		    }
		    return set_m_min_min(result);
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

	    /* Refreshes internal generation parameters. */

	    bool refresh_params()
	    {
		bool q=true;
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    q&=((*it)->refresh_params());
		}
		return q;
	    }

	    /* Updates all subchannels with the argument value: */

	    void update(const value_type& integrand)
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->integrand()=integrand;
		    (*it)->update();
		}
	    }

	    /* Adapts all subchannels: */

	    void adapt_grids()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->adapt_grids();
		}
	    }

	    /* Adapts all subchannels: */

	    void adapt_channels()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->adapt_channels();
		}
	    }

	    /* Adapts all subchannels: */

	    void adapt()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->adapt();
		}
	    }

	    /* Resets all subchannel multichannel configurations:
	     * */

	    void reset()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->reset();
		}
	    }

	    /* Resets all subchannel cross sections: * */

	    void reset_cross_section()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->reset_cross_section();
		}
	    }

	    /* Sets all particle channel redundancy flags: */
	    
	    void set_redundant()
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->set_redundant();
		}
	    }

	    /* Cleans the channel from all redundant subchannels: */

	    void clean()
	    {
		particle_channel_iterator it=std::stable_partition(particle_channels.begin(),particle_channels.end(),necessary_particle);
		for(particle_channel_iterator it2=particle_channels.begin();it2!=it;++it2)
		{
		    (*it2)->clean();
		}
		for(particle_channel_iterator it2=it;it2!=particle_channels.end();++it2)
		{
		    delete *it2;
		}
		particle_channels.erase(it,particle_channels.end());
	    }

	    /* Sets the multichannel adaptivity and threshold to the argument
	     * pair: */

	    void set_multichannel_params(const std::pair<value_type,value_type>& x)
	    {
		for(particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->set_multichannel_params(x);
		}
	    }

	    /* Public readout functions */
	    /*--------------------------*/

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

	    const value_type& m() const
	    {
		return signed_mass;
	    }

	    /* Returns a const reference to the minimal invariant mass-squared.
	     * */

	    const value_type& s_min_min() const
	    {
		return sminmin;
	    }

	    /* Returns a const reference to the minimal signed invariant mass.
	     * */

	    const value_type& m_min_min() const
	    {
		return mminmin;
	    }

	    /* Returns a const reference to the minimal invariant mass-squared.
	     * */

	    const value_type& s_max_max() const
	    {
		return smaxmax;
	    }

	    /* Returns a const reference to the minimal signed invariant mass.
	     * */

	    const value_type& m_max_max() const
	    {
		return mmaxmax;
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

	    /* Returns whether the momentum channel is redundant: */

	    bool is_redundant() const
	    {
		for(const_particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    if(!((*it)->is_redundant()))
		    {
			return false;
		    }
		}
		return true;
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

	    /* Returns whether the channel has just been generated: */

	    bool is_generated() const
	    {
		return generation_flag;
	    }

	    /* Normalisation checking utility: */

	    std::ostream& check_s_norm(std::ostream& os) const
	    {
		if(on_shell())
		{
		    return os;
		}
		for(const_particle_channel_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->check_s_norm(os);
		}
		return os;
	    }

	    /* Serialization */
	    /*---------------*/

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
		os<<std::setw(23)<<std::left<<signed_mass;
		os<<std::setw(23)<<std::left<<mminmin;
		os<<std::setw(23)<<std::left<<mmaxmax<<std::endl;
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
		os<<bitstring<<"\t, p = "<<*momentum<<", m = "<<signed_mass<<", [m--,m++] = ["<<mminmin<<','<<mmaxmax<<']';
		return os;
	    }

	    /* Overridden loading function: */

	    std::istream& load(std::istream& is)
	    {
		safe_read(is,signed_mass);
		safe_read(is,mminmin);
		safe_read(is,mmaxmax);
		
		invariant_mass=sgn_sq(signed_mass);
		sminmin=sgn_sq(mminmin);
		smaxmax=sgn_sq(mmaxmax);
		std::string initflag;
		std::getline(is,initflag);
		return is;
	    }

	    /* Overridden saving function: */

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<channel>"<<std::endl;
		os<<bitstring<<std::endl;
		safe_write(os,signed_mass);
		os<<"\t";
		safe_write(os,mminmin);
		os<<"\t";
		safe_write(os,mmaxmax);
		os<<std::endl;
		for(typename particle_channel_list::const_iterator it=particle_channels.begin();it!=particle_channels.end();++it)
		{
		    (*it)->save(os);
		}
		os<<"</channel>"<<std::endl;
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
	    
	    /* Signed invariant mass: */
	    
	    value_type signed_mass;

	    /* Extremal invariant mass bounds: */

	    value_type sminmin,mminmin,smaxmax,mmaxmax;

	    /* Vector of references for computing smin: */

	    std::vector<std::vector<const value_type*> > mmin_decomps;

	    /* Vector of references for computing the adjoint smin: */

	    std::vector<std::vector<const value_type*> > mmin_decomps_bar;

	    /* Generation flag: */

	    bool generation_flag;
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

