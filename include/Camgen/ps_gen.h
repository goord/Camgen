//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_gen.h
    \brief Phase space generator base class.
 */

#ifndef CAMGEN_PS_GEN_H_
#define CAMGEN_PS_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * *
 * Base class for momentum generators in Camgen.   *
 *                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/CM_algo.h>
#include <Camgen/init_state.h>
#include <Camgen/gen_conf.h>
#include <Camgen/MC_int.h>
#include <Camgen/evt_gen_base.h>

namespace Camgen
{
    /// Phase space generator base class template.

    template<class model_t,std::size_t N_in,std::size_t N_out>class ps_generator: public MC_integrator<typename model_t::value_type>,
                                                                                  public event_generator_base<model_t,N_in,N_out>
    {
	typedef MC_integrator<typename model_t::value_type> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
            typedef typename event_generator_base<model_t,N_in,N_out>::event_type event_type;
            typedef typename event_type::spacetime_type spacetime_type;
	    typedef typename event_type::momentum_type momentum_type;
	    typedef typename event_type::value_type value_type;
	    typedef typename event_type::size_type size_type;
	    typedef initial_state<model_t,N_in> init_state_type;

	    /* Public data members */
	    /*---------------------*/

	    /* Momentum-space integration factor: */

	    static value_type ps_factor;

	    /* Public constructors */
	    /*---------------------*/

	    /// Constructor with initial state argument.

	    ps_generator(init_state_type* is_):is(is_),isw(0),fsw(0)
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    pin[i]=new momentum_type;
		    alloc_pin.set(i);
		    is->set_p(i,pin[i]);
		    min[i]=NULL;
		}
		for(size_type i=0;i<N_out;++i)
		{
		    pout[i]=new momentum_type;
		    alloc_pout.set(i);
		    mout[i]=NULL;
		    for(size_type j=0;j<N_out;++j)
		    {
			mmin[i][j]=(value_type)0;
			eff_mmin[i][j]=(value_type)0;
		    }
		}
		is->template assign_mmin_components<N_out>(eff_mmin);
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /// Destructor.

	    virtual ~ps_generator()
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    if(alloc_pin[i])
		    {
			delete pin[i];
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    if(alloc_pout[i])
		    {
			delete pout[i];
		    }
		}
		if(is!=NULL)
		{
		    delete is;
		}
	    }

	    /* Public modifiers */
	    /*------------------*/

	    /// Sets the i-th beam energy pointer

	    bool set_beam_energy(int i,const value_type& E)
	    {
		if(E>(value_type)0 and i<0 and i>=-(int)N_in)
		{
		    is->set_beam_energy(-i-1,E);
		    return true;
		}
		return false;
	    }

	    /// Sets the i-th momentum (i<0 denotes incoming momenta) address to
	    /// the pointer p.

	    bool set_p(int i,momentum_type* p)
	    {
		if(p==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to assign null pointer to momentum address ignored"<<endlog;
		}
		if(i<(int)0 and i>=-(int)N_in)
		{
		    size_type n=-i-1;
		    if(alloc_pin[n])
		    {
			delete pin[n];
		    }
		    pin[n]=p;
		    alloc_pin[n]=false;
		    is->set_p(n,p);
		    return true;
		}
		if(i>(int)0 and i<=(int)N_out)
		{
		    size_type n=i-1;
		    pout[n]=p;
		    alloc_pout[n]=false;
		    return true;
		}
		return false;
	    }

	    /// Sets the i-th mass (i<0 denotes incoming momenta) address to
	    /// the pointer m.

	    bool set_m(int i,const value_type* m)
	    {
		if(i<(int)0 and i>=-(int)N_in)
		{
		    size_type n=-i-1;
		    min[n]=m;
		    is->set_m(n,m);
		    return true;
		}
		if(i>(int)0 and i<=(int)N_out)
		{
		    size_type n=i-1;
		    mout[n]=m;
		    return true;
		}
		return false;
	    }

	    /// Sets the minimal invariant dimass of outgoing particles i and j
	    /// to sqrts.

	    virtual bool set_m_min(int i,int j,const value_type& sqrts)
	    {
		if(i!=j and i>0 and j>0 and i<=(int)N_out and j<=(int)N_out)
		{
		    mmin[i-1][j-1]=sqrts;
		    mmin[j-1][i-1]=sqrts;
		    return true;
		}
		return false;
	    }

	    /// Refreshes internal parameters to lower dimass limits.

	    virtual bool refresh_m_min()
	    {
		for(size_type i=0;i<N_out;++i)
		{
		    for(size_type j=i;j<N_out;++j)
		    {
			eff_mmin[i][j]=eval_m_min(i,j);
			eff_mmin[j][i]=eff_mmin[i][j];
		    }
		}
		is->refresh_m_min();
		return true;
	    }

	    /// Refreshes the total hadronic invariant mass.

	    virtual bool refresh_Ecm()
	    {
		if(is->refresh_Ecm())
                {
                    for(size_type i=0;i<N_in;++i)
                    {
                        this->get_event_ptr()->set_beam_energy(-1-i,is->beam_energy(i));
                        this->get_event_ptr()->set_beam_id(-1-i,is->beam_id(i));
                    }
                    return true;
                }
                return false;
	    }

	    /// Refreshes internal parameters.

	    virtual bool refresh_params()
	    {
		return true;
	    }

	    /// Refreshes the partonic invariant mass.

	    virtual bool refresh_Ecm_hat()
	    {
                this->get_event_ptr()->set_Ecm_hat(Ecm_hat());
		return check_sufficient_shat();
            }

	    /// Generation method.

	    virtual bool generate()
	    {
                this->get_event_ptr()->reset();
		if(!generate_is())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!refresh_Ecm_hat())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!generate_fs())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
                copy_event_data();
		collect_integrand();
		collect_weight();
		return true;
	    }

	    /// Weight evaluation implementation.

	    virtual bool evaluate_weight()
	    {
		if(!evaluate_is_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(!evaluate_fs_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		collect_integrand();
		collect_weight();
		return true;
	    }

	    /// Overridden updating method.

	    virtual void update()
	    {
		is->integrand()=(this->integrand())*(this->weight());
		is->update();
	    }

	    /// VEGAS-grid adaptation method.

	    virtual void adapt_grids()
	    {
		is->adapt_grids();
	    }

	    /// Multichannel weight adaptation method.

	    virtual void adapt_channels()
	    {
		is->adapt_channels();
	    }

	    /// Overridden adaption method.

	    void adapt()
	    {
		adapt_grids();
		adapt_channels();
	    }

	    /// Replaces current multichannel weights by there optimal values,
	    /// collected over several adaptations.

	    virtual void set_optimal_multichannel(){}

	    /// Sets all momentum and mass addresses to the current tree's
	    /// external legs momenta and masses.

	    virtual bool set_amplitude(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
                if(this->get_event_ptr()!=NULL and this->allocated_event())
                {
                    this->get_event_ptr()->set_process(new sub_process<model_t,N_in,N_out>(it->get_phi_in(),it->get_phi_out()));
                }
		for(size_type i=0;i<N_in;++i)
		{
		    if(it->get_phase_space(i)==NULL)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"NULL phase space instance encountered in Camgen current tree--returning NULL"<<endlog;
			return false;
		    }
		    set_p_in(i,it->get_phase_space(i)->get_momentum());
		    set_m_in(i,it->get_phase_space(i)->particle_type->get_mass_address());
		    is->set_parton(i,it->get_phase_space(i)->particle_type->get_pdg_id());
		}
		for(size_type i=0;i<N_out;++i)
		{
		    size_type n=i+N_in;
		    if(it->get_phase_space(n)==NULL)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"NULL phase space instance encountered in Camgen current tree--returning NULL"<<endlog;
			return false;
		    }
		    set_p_out(i,it->get_phase_space(n)->get_momentum());
		    set_m_out(i,it->get_phase_space(n)->particle_type->get_mass_address());
		}
		for(size_type i=0;i<N_out;++i)
		{
		    for(size_type j=i;j<N_out;++j)
		    {
			eff_mmin[i][j]=eval_m_min(i,j);
			eff_mmin[j][i]=eff_mmin[i][j];
		    }
		}
		return true;
	    }

	    /* Public readout methods */
	    /*------------------------*/

            /// Subprocess cross section implementation: returns the current cross section if the argument is 1.

            MC_integral<value_type> process_xsec(int proc_id) const
            {
                return proc_id==1?(this->cross_section()):MC_integral<value_type>();
            }

            /// Subprocess maximum weight implementation: returns the current maximal weight if the argument is 1.

            value_type process_maxw(int proc_id) const
            {
                return proc_id==1?(this->max_weight()):(value_type)0;
            }

	    /// Returns the i-th incoming beam identifier.

	    int beam_id(int i) const
	    {
		return is->beam_id(-i-1);
	    }

	    /// Returns the cernlib pdf group number for the i-th incoming beam.

	    int pdfg(int i) const
	    {
		return is->pdfg(-i-1);
	    }

	    /// Returns the cernlib pdf set number for the i-th incoming beam.

	    int pdfs(int i) const
	    {
		return is->pdfs(-i-1);
	    }

	    /// Returns the const reference to the i-th incoming momentum (no
	    /// range checking on i).

	    const momentum_type& p_in(size_type i) const
	    {
		return *(pin[i]);
	    }

            /// Returns the total incoming momentum.

            momentum_type P_in() const
            {
                momentum_type p(p_in(0));
                for(size_type i=1;i<N_in;++i)
                {
                    p+=p_in(i);
                }
                return p;
            }

	    /// Returns the const reference to the i-th incoming momentum (no
	    /// range checking on i).

	    const momentum_type& p_out(size_type i) const
	    {
		return *(pout[i]);
	    }

            /// Returns the total outgoing momentum.

            momentum_type P_out() const
            {
                momentum_type p(p_out(0));
                for(size_type i=1;i<N_out;++i)
                {
                    p+=p_out(i);
                }
                return p;
            }
	    
	    /// Returns the i-th incoming mass (no range checking on i).

	    value_type M_in(size_type i) const
	    {
		return (min[i]==NULL)?((value_type)0):(*(min[i]));
	    }

	    /// Returns the i-th outgoing mass (no range checking on i).

	    value_type M_out(size_type i) const
	    {
		return (mout[i]==NULL)?((value_type)0):(*(mout[i]));
	    }

            /// Returns the incoming mass sum.

            value_type M_in() const
            {
                value_type m(0);
                for(size_type i=0;i<N_in;++i)
                {
                    m+=M_in(i);
                }
                return m;
            }

            /// Returns the outgoing mass sum.

            value_type M_out() const
            {
                value_type m(0);
                for(size_type i=0;i<N_out;++i)
                {
                    m+=M_out(i);
                }
                return m;
            }

	    /// Returns the total hadronic CM-frame energy.

	    value_type Ecm() const
	    {
		return is->Ecm();
	    }

	    /// Returns the total hadronic invariant mass-squared.

	    value_type s_tot() const
	    {
		return is->s();
	    }

	    /// Returns the total partonic CM-frame energy.

	    value_type Ecm_hat() const
	    {
		return is->Ecm_hat();
	    }

	    /// Returns the total partonic invariant mass-squared.

	    value_type s_hat() const
	    {
		return is->s_hat();
	    }

	    /// Returns the initial-state weight.

	    const value_type& is_weight() const
	    {
		return isw;
	    }

	    /// Returns the final-state weight.

	    const value_type& fs_weight() const
	    {
		return fsw;
	    }

	    /// Returns the i-th beam energy.

	    value_type beam_energy(int i) const
	    {
		return is->beam_energy(-i-1);
	    }

	    /// Returns whether s-hat is kept fixed and equal to s by the
	    /// initial state.

	    bool fixed_s_hat() const
	    {
		return !(is->hadronic);
	    }

	    /// Returns the minimal invariant dimass of particles i,j.

	    value_type m_min(int i,int j) const
	    {
		return eff_mmin[i-1][j-1];
	    }
	    
	    /// Incoming state flux factor.

	    value_type flux_factor()
	    {
		is->mu_F=this->event_generator_base<model_t,N_in,N_out>::F_scale();
		return is->flux_factor();
	    }

	    /// Returns the const pointer to the total hadronic invariant mass

	    const value_type* get_Ecm() const
	    {
		return &(is->Ecm());
	    }

	    /// Returns the const pointer to the total partonic invariant mass

	    const value_type* get_Ecm_hat() const
	    {
		return &(is->Ecm_hat());
	    }

	    value_type alpha_s() const
	    {
		return is->alpha_s();
	    }

	    /// Checks parton level momentum conservation condition

	    bool check_sufficient_shat() const
	    {
                return this->get_event().check_sufficient_shat();
	    }

	    /// Checks hadron level momentum conservation condition

	    bool check_sufficient_s() const
	    {
                return this->get_event().check_sufficient_s();
	    }

	    /// Checks momentum conservation and masses.

	    bool check() const
	    {
		return this->get_event().check_p_on_shell() and this->get_event().check_p_conservation();
	    }

	    /// Checks whether the minimal dimass constraints are fulfilled (no
	    /// logging).

	    bool pass_cuts()
	    {
		for(size_type i=0;i<N_out;++i)
		{
		    momentum_type q(p_out(i));
		    for(size_type j=i+1;j<N_out;++j)
		    {
			if(eff_mmin[i][j]>(value_type)0)
			{
			    if(event<model_t,N_in,N_out>::s(q+p_out(j))<eff_mmin[i][j]*eff_mmin[i][j])
			    {
				return false;
			    }
			}
		    }
		}
		return this->event_generator_base<model_t,N_in,N_out>::pass_cuts();
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Polymorphic type identifier.

	    virtual std::string type() const=0;

	    /// Printing method (trivial by default).

	    virtual std::ostream& print(std::ostream& os) const
	    {
		return os;
	    }

	    /// Prints the minimal dimass cuts.

	    std::ostream& print_m_min(std::ostream& os) const
	    {
		os<<"sqrt(s-) = ";
		for(size_type i=0;i<N_out;++i)
		{
		    if(i==0)
		    {
			os<<std::setw(15)<<'[';
		    }
		    else
		    {
			os<<"           "<<std::setw(15)<<'[';
		    }
		    for(size_type j=0;j<N_out;++j)
		    {
			os<<std::setw(15)<<eff_mmin[i][j];
		    }
		    os<<']'<<std::endl;
		}
		os<<std::endl;
		return os;
	    }

	    /// Virtual method printing the fs-generator settings.

	    virtual std::ostream& print_fs_settings(std::ostream& os) const
	    {
		return os;
	    }

	    /// Method printing the full generator settings.

	    std::ostream& print_settings(std::ostream& os) const
	    {
		os<<std::setw(30)<<std::left<<"s-hat sampling:"<<(is->s_hat_sampling?"fs":"is")<<std::endl<<std::endl;
		os<<std::setw(30)<<std::left<<"initial state:"<<is->type()<<std::endl;
		is->print_settings(os);
		os<<std::endl<<std::setw(30)<<std::left<<"final state:"<<type()<<std::endl;
		os<<std::setw(30)<<std::left<<"sqrt(s-):";
		for(size_type i=0;i<N_out;++i)
		{
		    if(i==0)
		    {
			os<<std::setw(15)<<'[';
		    }
		    else
		    {
			os<<std::setw(30)<<std::left<<" "<<std::setw(15)<<'[';
		    }
		    for(size_type j=0;j<N_out;++j)
		    {
			os<<std::setw(15)<<mmin[i][j];
		    }
		    os<<']'<<std::endl;
		}
		print_fs_settings(os);
		return os;
	    }

	    /// Prints the status of the generator.
	    
	    std::ostream& print_status(std::ostream& os=std::cout) const
	    {
		os<<"###############################################################################################"<<std::endl;
		os<<"Nr of events contributing to cross-section:        "<<std::scientific<<this->calls()<<std::endl;
		os<<"Monte Carlo efficiency (%):                        "<<std::scientific<<this->efficiency()<<std::endl;
		os<<"Cross section (pb):                                "<<std::scientific<<this->cross_section()<<std::endl;
		os<<"###############################################################################################"<<std::endl;
		return os;
	    }

	    /// Prints the ingredients of the returned weight.

	    std::ostream& print_weights(std::ostream& os=std::cout) const
	    {
		os<<"###############################################################################################"<<std::endl;
		os<<"Initial state weight:				"<<std::scientific<<is_weight()<<std::endl;
		os<<"Final state weight:				"<<std::scientific<<fs_weight()<<std::endl;
		os<<"Full integrand:                                    "<<std::scientific<<this->integrand()<<std::endl;
		os<<"Full weight:                                       "<<std::scientific<<this->weight()<<std::endl;
		os<<"###############################################################################################"<<std::endl;
		return os;
	    }

	protected:

	    /* Minimal invariant mass pointers: */

	    value_type mmin[N_out][N_out];

	    /* Effective minimal outgoing invariant dimasses: */

	    value_type eff_mmin[N_out][N_out];

	    /// Calls the initial state generator.

	    virtual bool generate_is()
	    {
		if(!(is->generate()))
		{
		    isw=(value_type)0;
		    return false;
		}
		isw=is->weight();
		return true;
	    }

	    /// Calls the initial state weight evaluation method.
	    
	    bool evaluate_is_weight()
	    {
		if(!is->evaluate_weight())
		{
		    isw=(value_type)0;
		    return false;
		}
		isw=is->weight();
		return true;
	    }

	    /// Virtual method generating the final state.

	    virtual bool generate_fs()
	    {
		return true;
	    }

	    /// Virtual method evaluating the final state weight.
	    
	    virtual bool evaluate_fs_weight()
	    {
		return true;
	    }

	    /// Sets the i-th incoming momentum (no range checking on i)

	    void set_p_in(size_type i,momentum_type* p)
	    {
		if(alloc_pin[i])
		{
		    delete pin[i];
		}
		alloc_pin[i]=false;
		pin[i]=p;
		is->set_p(i,p);
	    }

	    /// Sets the i-th outgoing momentum (no range checking on i)

	    void set_p_out(size_type i,momentum_type* p)
	    {
		if(alloc_pout[i])
		{
		    delete pout[i];
		}
		alloc_pout[i]=false;
		pout[i]=p;
	    }

	    /// Sets the i-th incoming mass pointer (no range checking on i).

	    void set_m_in(size_type i,const value_type* m)
	    {
		min[i]=m;
		is->set_m(i,m);
	    }

	    /// Sets the i-th outgoing mass pointer (no range checking on i).

	    void set_m_out(size_type i,const value_type* m)
	    {
		mout[i]=m;
	    }

	    /// Returns the reference to the i-th incoming momentum (no range
	    /// checking on i).

	    momentum_type& p_in(size_type i)
	    {
		return *(pin[i]);
	    }

	    /// Returns the reference to the i-th outgoing momentum (no range
	    /// checking on i).

	    momentum_type& p_out(size_type i)
	    {
		return *(pout[i]);
	    }

	    /// Returns the initial-state pointer.
	    
	    init_state_type* is_generator()
	    {
		return is;
	    }

	    /// Returns the initial-state const-pointer.
	    
	    const init_state_type* is_generator() const
	    {
		return is;
	    }

	    /// Returns the initial state weight reference.

	    value_type& is_weight()
	    {
		return isw;
	    }

	    /// Returns the final state weight reference.

	    value_type& fs_weight()
	    {
		return fsw;
	    }

	    /// Returns whether the outgoing state has massive particles.

	    bool massive_fs() const
	    {
		for(size_type i=0;i<N_out;++i)
		{
		    if(mout[i]!=NULL)
		    {
			if(*(mout[i])!=(value_type)0)
			{
			    return true;
			}
		    }
		}
		return false;
	    }

	    /* Evaluates the minimal s_{ij}: */

	    value_type eval_m_min(size_type i,size_type j) const
	    {
		if(i==j)
		{
		    return this->M_out(i);
		}
		return std::max(this->M_out(i)+this->M_out(j),mmin[i][j]);
	    }

	    /* Evaluates the generated integrand: */

	    void collect_integrand()
	    {
		this->integrand()=ps_factor*flux_factor();
	    }

	    /* Evaluates the generated weight: */

	    void collect_weight()
	    {
		this->weight()=pass_cuts()?(isw*fsw):(value_type)0;
	    }

            /* Copies momenta to the event instance: */

            void copy_event_data()
            {
                this->get_event_ptr()->set_w(this->weight());
                this->get_event_ptr()->set_max_w(this->max_weight());
                this->get_event_ptr()->set_xsec(this->cross_section());
                this->get_event_ptr()->set_process_xsec(this->cross_section());
                for(size_type i=0;i<N_in;++i)
                {
                    this->get_event_ptr()->set_p_in(i,*(pin[i]));
                }
                for(size_type i=0;i<N_out;++i)
                {
                    this->get_event_ptr()->set_p_out(i,*(pout[i]));
                }
            }
	    
	private:

	    /* Initial state generator: */
	    
	    init_state_type* is;

	    /* Initial-state weight: */

	    value_type isw;

	    /* Final-state weight: */

	    value_type fsw;

	    /* Incoming momenta: */

	    momentum_type* pin[N_in];
	    
	    /* Outgoing momenta: */
	    
	    momentum_type* pout[N_out];

	    /* Incoming momentum allocation flags: */

            // TODO: Use smart pointers for the momenta
	    std::bitset<N_in>alloc_pin;
	    
	    /* Outgoing momentum allocation flags: */
	    
	    std::bitset<N_out>alloc_pout;

	    /* Incoming masses: */

	    const value_type* min[N_in];
	    
	    /* Outgoing masses: */
	    
	    const value_type* mout[N_out];
    };
    template<class model_t,std::size_t N_in,std::size_t N_out>typename ps_generator<model_t,N_in,N_out>::value_type ps_generator<model_t,N_in,N_out>::ps_factor=std::pow((typename model_t::value_type)2*std::acos(-(typename model_t::value_type)1),(int)model_t::dimension-(int)(model_t::dimension-1)*(int)N_out);
}

#endif /*CAMGEN_PS_GEN_H_*/

