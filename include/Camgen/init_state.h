//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file init_state.h
    \brief Base class for initial state generators.
 */

#ifndef CAMGEN_INIT_STATE_H_
#define CAMGEN_INIT_STATE_H_

#include <bitset>
#include <Camgen/particle.h>
#include <Camgen/rn_strm.h>
#include <Camgen/bipart.h>
#include <Camgen/MC_gen.h>

namespace Camgen
{
    template<class model_t,std::size_t N,class spacetime_t=typename model_t::spacetime_type>class partonic_is;
    template<class model_t,class rng_t,class spacetime_t=typename model_t::spacetime_type>class hadronic_is_xx;
    template<class model_t,class rng_t,class spacetime_t=typename model_t::spacetime_type>class hadronic_is_sy;
    template<class model_t,class rng_t,class spacetime_t=typename model_t::spacetime_type>class hadronic_is_y;
    
    /// Initial-state Monte Carlo generator base class template.
    
    template<class model_t,std::size_t N>class initial_state: public MC_generator<typename model_t::value_type>
    {
	typedef MC_generator<typename model_t::value_type> base_type;

	public:

	    /* Type definitions: */

	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef typename model_t::value_type value_type;
	    typedef vector<value_type,model_t::dimension> momentum_type;
	    typedef typename momentum_type::size_type size_type;

	    /* Public data members */
	    /*---------------------*/

	    /// Boolean denoting whether the generator is hadronic (no fixed
	    /// partonic invariant mass).

	    const bool hadronic;

	    /// Boolean denoting whether the generator is adaptive.

	    const bool adaptive;

	    /// Boolean denoting whether the generator samples s-hat or not.

	    const bool s_hat_sampling;

	    /// Factorisation scale (default value MZ=91.1876).

	    value_type mu_F;

	    /* Public constructors */
	    /*---------------------*/

	    /// Copy constructor.

	    initial_state(const initial_state<model_t,N>& other):hadronic(other.hadronic),adaptive(other.adaptive),s_hat_sampling(other.s_hat_sampling),mu_F(other.mu_F),E_beams(other.E_beams),ecm(other.ecm),stot(other.stot),ecmhat(other.ecmhat),shat(other.shat),alloc_momenta(other.alloc_momenta),masses(other.masses),partons(other.partons),mhatmin(other.mhatmin)
	    {
		for(size_type i=0;i<N;++i)
		{
		    if(alloc_momenta[i])
		    {
			momenta[i]=new momentum_type(*(other.momenta[i]));
		    }
		    else
		    {
			momenta[i]=other.momenta[i];
		    }
		}
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /// Destructor.

	    virtual ~initial_state()
	    {
		for(size_type i=0;i<N;++i)
		{
		    if(alloc_momenta[i])
		    {
			delete momenta[i];
		    }
		}
	    }

	    /* Public modifiers */
	    /*------------------*/

	    /// Sets the i-th beam energy

	    bool set_beam_energy(size_type i,const value_type& E)
	    {
		if(i>=N)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"beam label "<<i<<" out of range"<<endlog;
		    return false;
		}
		if(E<(value_type)0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"assigned beam energy "<<E<<" invalid"<<endlog;
		    return false;
		}
		E_beams[i]=E;
		return true;
	    }

	    /// Sets the partonic invariant mass.
	    
	    bool set_Ecm_hat(const value_type& sqrts)
	    {
		if(!(sqrts<mhatmin))
		{
		    ecmhat=sqrts;
		    shat=sqrts*sqrts;
		    return true;
		}
		return false;
	    }

	    /// Sets the partonic invariant mass-squared.
	    
	    bool set_s_hat(const value_type& s)
	    {
		if(!(s<shatmin))
		{
		    shat=s;
		    ecmhat=std::sqrt(s);
		    return true;
		}
		return false;
	    }

	    /// Returns the i-th parton momentum.

	    momentum_type& p(size_type i)
	    {
		return *(momenta[i]);
	    }

	    /// Sets the i-th momentum address.

	    virtual void set_p(size_type i,momentum_type* p)
	    {
		if(alloc_momenta[i])
		{
		    delete momenta[i];
		    alloc_momenta.reset(i);
		}
		momenta[i]=p;
	    }

	    /// Sets the i-th mass address.

	    virtual void set_m(size_type i,const value_type* m=NULL)
	    {
		masses[i]=m;
	    }

	    /// Sets the i-th initial state parton flavour.
	    
	    virtual void set_parton(size_type i,int j=0)
	    {
		partons[i]=j;
	    }

	    /// Sets the minimal partonic invariant mass to the argument.
	    
	    virtual bool set_m_hat_min(const value_type& m)
	    {
		if(m>(value_type)0 and m<ecm)
		{
		    mhatmin=m;
		    shatmin=m*m;
		    return true;
		}
		return false;
	    }

	    /// Sets the minimal partonic invariant mass-squared to the argument.
	    
	    virtual bool set_s_hat_min(const value_type& s)
	    {
		if(s>(value_type)0 and s<stot)
		{
		    shatmin=s;
		    mhatmin=std::sqrt(s);
		    return true;
		}
		return false;
	    }

	    /// Re-evaluates mhat_min from branchings.

	    bool refresh_m_min()
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
		value_type msum=(value_type)0;
		for(size_type i=0;i<N;++i)
		{
		    msum+=m(i);
		}
		return set_m_hat_min(std::max(result,msum));
	    }

	    /// Re-evaluates the total hadronic CM-energy.

	    virtual bool refresh_Ecm()=0;

	    /// Refreshes internal parameters.

	    virtual bool refresh_params()
	    {
		return true;
	    }

	    /* Function filling the mmin components for timelike channels: */

	    template<std::size_t N_out>void assign_mmin_components(value_type mminout[N_out][N_out])
	    {
		std::set<size_type>s;
		for(size_type i=0;i<N_out;++i)
		{
		    s.insert(i);
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
	    }

	    /// VEGAS-grid adaptation method.

	    virtual void adapt_grids(){}

	    /// Multichannel weight adaptation method.

	    virtual void adapt_channels(){}

	    /// Overridden adaption method.

	    void adapt()
	    {
		adapt_grids();
		adapt_channels();
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /// Virtual cloning method.

	    virtual initial_state<model_t,N>* clone() const=0;
	    
	    /// Returns the i-th beam energy (if the pointer equals zero,
	    /// returns the parton mass).

	    value_type beam_energy(size_type i) const
	    {
		return E_beams[i];
	    }

	    /// Returns the i-th beam particle type.

	    virtual int beam_id(size_type i) const
	    {
		return 0;
	    }

	    /// Returns the cernlib pdf group number.

	    virtual int pdfg(size_type i) const
	    {
		return 0;
	    }

	    /// Returns the cernlib pdf set number.

	    virtual int pdfs(size_type i) const
	    {
		return 0;
	    }

	    /// Returns the i-th parton momentum.

	    const momentum_type& p(size_type i) const
	    {
		return *(momenta[i]);
	    }

	    /// Returns the i-th parton mass.

	    value_type m(size_type i) const
	    {
		return (masses[i]==NULL)?(value_type)0:(*(masses[i]));
	    }

	    /// Returns the i-th parton mass-squared.

	    value_type m2(size_type i) const
	    {
		return (masses[i]==NULL)?(value_type)0:((*(masses[i]))*(*(masses[i])));
	    }

	    /// Returns the i-th parton pdg id.

	    int parton_id(size_type i) const
	    {
		return partons[i];
	    }

	    /// Returns the collider CM-energy.
	    
	    const value_type& Ecm() const
	    {
		return ecm;
	    }

	    /// Returns the collider energy-squared.

	    const value_type& s() const
	    {
		return stot;
	    }

	    /// Returns the partonic invariant mass.

	    const value_type& Ecm_hat() const
	    {
		return ecmhat;
	    }

	    /// Returns the partonic invariant mass-squared.

	    const value_type& s_hat() const
	    {
		return shat;
	    }

	    /// Returns the minimal partonic invariant mass const reference.

	    const value_type& m_hat_min() const
	    {
		return mhatmin;
	    }

	    /// Returns the minimal partonic invariant mass-squared.

	    const value_type& s_hat_min() const
	    {
		return shatmin;
	    }

	    /// Returns the flux factor.

	    virtual value_type flux_factor() const=0;

	    /// Returns the strong coupling.
	    
	    virtual value_type alpha_s() const
	    {
		return value_type(-1);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Derived type identifier.

	    virtual std::string type() const=0;

	    /// Prints settings.
	    
	    virtual std::ostream& print_settings(std::ostream& os) const
	    {
		os<<std::setw(30)<<std::left<<"hadronic:"<<(hadronic?"yes":"no")<<std::endl;
		os<<std::setw(30)<<std::left<<"adaptive:"<<(adaptive?"yes":"no")<<std::endl;
		for(size_type i=0;i<N;++i)
		{
		    std::stringstream ss;
		    ss<<"beam energy -"<<(i+1)<<':';
		    os<<std::setw(30)<<ss.str()<<E_beams[i]<<std::endl;
		}
		os<<std::setw(30)<<std::left<<"Ecm:"<<ecm<<std::endl;
		if(hadronic)
		{
		    os<<std::setw(30)<<std::left<<"mu_F:"<<mu_F<<std::endl;
		    os<<std::setw(30)<<std::left<<"mhatmin:"<<m_hat_min()<<std::endl;
		}
		return os;
	    }

	    static initial_state<model_t,N>* create_instance(std::istream& is)
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<isgen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read"<<endlog;
		    return NULL;
		}
		initial_state<model_t,N>* result;
		std::string br_id;
		is>>br_id;
	    }

	    /// Loads derived type data.

	    virtual std::istream& load_data(std::istream& is)
	    {
		return is;
	    }

	    /// Loading method implementation.

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		for(size_type i=0;i<N;++i)
		{
		    is>>E_beams[i];
		}
		safe_read(is,mhatmin);
		safe_read(is,mu_F);
		stot=ecm*ecm;
		shat=ecmhat*ecmhat;
		shatmin=mhatmin*mhatmin;
		for(size_type i=0;i<N;++i)
		{
		    is>>partons[i];
		}
		return load_data(is);
	    }

	    /// Outputs derived type data.

	    virtual std::ostream& save_data(std::ostream& os) const
	    {
		return os;
	    }

	    /// Serializing function.

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<isgen>"<<std::endl;
		os<<type()<<std::endl;
		this->base_type::save(os);
		for(size_type i=0;i<N;++i)
		{
		    safe_write(os,E_beams[i]);
		    os<<"\t";
		}
		safe_write(os,mhatmin);
		os<<"\t";
		safe_write(os,mu_F);
		os<<std::endl;
		for(size_type i=0;i<N;++i)
		{
		    os<<partons[i]<<"\t"<<std::endl;
		}
		save_data(os);
		os<<"</isgen>"<<std::endl;
		return os;
	    }
	    
	protected:

	    /* Beam energies: */

	    vector<value_type,N>E_beams;

	    /* Protected constructors */
	    /*------------------------*/

	    /* Constructor with all beam energies NULL, i.e. all incoming
	     * particles at rest: */

	    initial_state(bool hadronic_=false,bool adaptive_=false,bool s_hat_sampling_=true):hadronic(hadronic_),adaptive(adaptive_),s_hat_sampling(s_hat_sampling_),mu_F(91.1876),mhatmin(0),shatmin(0)
	    {
		E_beams.assign(0);
		masses.assign(NULL);
		partons.assign(0);
		for(size_type i=0;i<N;++i)
		{
		    momenta[i]=new momentum_type;
		}
		alloc_momenta.set();
	    }

	    /* Constructor from beam energies: */

	    initial_state(const vector<value_type,N>& E_beams_,bool hadronic_=false,bool adaptive_=false,bool s_hat_sampling_=true):hadronic(hadronic_),adaptive(adaptive_),s_hat_sampling(s_hat_sampling_),mu_F(91.1876),E_beams(E_beams_),mhatmin(0),shatmin(0)
	    {
		masses.assign(NULL);
		partons.assign(0);
		for(size_type i=0;i<N;++i)
		{
		    momenta[i]=new momentum_type;
		}
		alloc_momenta.set();
	    }

	    /* Protected methods */
	    /*-------------------*/

	    /// Sets the hadronic invariant mass.
	    
	    bool set_Ecm(const value_type& sqrts)
	    {
		if(!(sqrts<mhatmin))
		{
		    ecm=sqrts;
		    stot=sqrts*sqrts;
		    return true;
		}
		return false;
	    }

	    /// Sets the hadronic invariant mass-squared.
	    
	    bool set_s(const value_type& s)
	    {
		if(s<shatmin)
		{
		    return false;
		}
		stot=s;
		ecm=std::sqrt(s);
		return true;
	    }

	private:

	    /* Private data members */
	    /*----------------------*/

	    /* Hadronic invariant mass: */

	    value_type ecm;

	    /* Hadronic invariant mass-squared:*/
	    
	    value_type stot;

	    /* Partonic invariant mass: */

	    value_type ecmhat;

	    /* Partonic invariant mass-squared: */

	    value_type shat;

	    /* Incoming phase space channels, whose momenta are generated by the
	     * initial state generator: */

	    vector<momentum_type*,N> momenta;

	    /* Momenta allocation flags: */

	    std::bitset<N> alloc_momenta;

	    /* Parton masses: */

	    vector<const value_type*,N>masses;
	    
	    /* Parton flavours: */
	    
	    vector<int,N>partons;

	    /* Minimal partonic invariant mass: */

	    value_type mhatmin;

	    /* Minimal partonic invariant mass-squared: */

	    value_type shatmin;

	    /* Minimal invariant mass decompositions: */

	    std::vector< std::vector<const value_type*> >mmin_decomps;
    };
}

#endif /*CAMGEN_INIT_STATE_H_*/

