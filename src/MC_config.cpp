//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <cmath>
#include <iostream>
#include <config.h>
#include <Camgen/MC_config.h>

namespace Camgen
{
    /* Sets the default helicity generator. */

    void set_helicity_generator_type(helicity_generators::type gen)
    {
	MC_config::set_helicity_generator_type(gen);
    }
    void MC_config::set_helicity_generator_type(helicity_generators::type gen)
    {
	MC_config::helgen=gen;
    }
    
    /* Sets the default colour generator. */

    void set_colour_generator_type(colour_generators::type gen)
    {
	MC_config::set_colour_generator_type(gen);
    }
    void MC_config::set_colour_generator_type(colour_generators::type gen)
    {
	MC_config::colgen=gen;
    }

    /* Sets the default initial state. */

    void set_initial_state_type(initial_states::type gen)
    {
	MC_config::set_initial_state_type(gen);
    }
    void MC_config::set_initial_state_type(initial_states::type gen)
    {
	MC_config::isgen=gen;
    }

    /* Sets the default phase space generator. */

    void set_phase_space_generator_type(phase_space_generators::type gen)
    {
	MC_config::set_phase_space_generator_type(gen);
    }
    void MC_config::set_phase_space_generator_type(phase_space_generators::type gen)
    {
	if(gen==phase_space_generators::recursive_backward_shat and (MC_config::isgen==initial_states::partonic or MC_config::isgen==initial_states::eplus_eminus))
	{
	    MC_config::psgen=phase_space_generators::recursive_backward_s;
	}
	else
	{
	    MC_config::psgen=gen;
	}
    }

    /* Sets the i-th beam energy, where i starts from 1,... */

    bool set_beam_energy(int i,double e)
    {
	return MC_config::set_beam_energy(i,e);
    }
    bool MC_config::set_beam_energy(int i,double e)
    {
	if(e>(double)0)
	{
	    if(i==1 or i==-1)
	    {
		MC_config::E1=e;
		return true;
	    }
	    if(i==2 or i==-2)
	    {
		MC_config::E2=e;
		return true;
	    }
	}
	return false;
    }

    /* Sets the first beam energy. Returns false at failure (negative argument).
     * */

    bool set_first_beam_energy(double e)
    {
	return MC_config::set_first_beam_energy(e);
    }
    bool MC_config::set_first_beam_energy(double e)
    {
	if(e>(double)0)
	{
	    MC_config::E1=e;
	    return true;
	}
	return false;
    }

    /* Sets the second beam energy. Returns false at failure (negative argument).
     * */

    bool set_second_beam_energy(double e)
    {
	return MC_config::set_second_beam_energy(e);
    }
    bool MC_config::set_second_beam_energy(double e)
    {
	if(e>(double)0)
	{
	    MC_config::E2=e;
	    return true;
	}
	return false;
    }

    /* Sets the number of channel iterations and batch size at initialisation.
     * */

    void set_channel_init(std::size_t n1,std::size_t n2)
    {
	MC_config::set_channel_init(n1,n2);
    }
    void MC_config::set_channel_init(std::size_t n1,std::size_t n2)
    {
	MC_config::ch_iters=n1;
	MC_config::ch_batch=n2;
    }

    /* Sets the number of grid iterations and batch size at initialisation. */

    void set_grid_init(std::size_t n1,std::size_t n2)
    {
	MC_config::set_grid_init(n1,n2);
    }
    void MC_config::set_grid_init(std::size_t n1,std::size_t n2)
    {
	MC_config::grd_iters=n1;
	MC_config::grd_batch=n2;
    }

    /* Sets the number of events for subprocess cross section estimation. */

    void set_subprocess_events(std::size_t n)
    {
	MC_config::set_subprocess_events(n);
    }
    void MC_config::set_subprocess_events(std::size_t n)
    {
	MC_config::proc_evts=n;
    }

    /* Sets the automatic channel adaptation batch size. If set zero, automatic
     * channel adaptation is switched off. */

    void set_auto_channel_adapt(std::size_t n)
    {
	MC_config::set_auto_channel_adapt(n);
    }
    void MC_config::set_auto_channel_adapt(std::size_t n)
    {
	MC_config::auto_ch_batch=n;
    }

    /* Sets the automatic grid adaptation batch size. If set zero, automatic
     * grid adaptation is switched off. */

    void set_auto_grid_adapt(std::size_t n)
    {
	MC_config::set_auto_grid_adapt(n);
    }
    void MC_config::set_auto_grid_adapt(std::size_t n)
    {
	MC_config::auto_grd_batch=n;
    }

    /* Sets the automatic subprocess weight adaptation batch size. If set zero,
     * automatic subprocess weights adaption is switched off. */

    void set_auto_subprocess_adapt(std::size_t n)
    {
	MC_config::set_auto_subprocess_adapt(n);
    }
    void MC_config::set_auto_subprocess_adapt(std::size_t n)
    {
	MC_config::auto_proc_batch=n;
    }

    /* Sets the maximal number of events thrown not passing cuts before moving
     * on at initialisation of channels and grids. */

    bool set_max_init_rejects(std::size_t n)
    {
	return MC_config::set_max_init_rejects(n);
    }
    bool MC_config::set_max_init_rejects(std::size_t n)
    {
	if(n==0)
	{
	    return false;
	}
	MC_config::mx_evts=n;
	return true;
    }

    /* Controls whether adaptive grids should optimise the s-sampling
     * algorithms. If set true, the one-shot method is adopted for the
     * generation of s-pairs. */

    void set_adaptive_s_sampling(bool q)
    {
	MC_config::set_adaptive_s_sampling(q);
    }
    void MC_config::set_adaptive_s_sampling(bool q)
    {
	MC_config::s_grids=q;
	if(q)
	{
	    MC_config::mx_throws=1;
	}
    }

    /* Controls whether adaptive grids should optimise the t-sampling
     * algorithms. */

    void set_adaptive_t_sampling(bool q)
    {
	MC_config::set_adaptive_t_sampling(q);
    }
    void MC_config::set_adaptive_t_sampling(bool q)
    {
	MC_config::t_grids=q;
    }

    /* Controls whether adaptive grids should optimise the s-branching
     * angle-sampling algorithms. */

    void set_adaptive_angles(bool q)
    {
	MC_config::set_adaptive_angles(q);
    }
    void MC_config::set_adaptive_angles(bool q)
    {
	MC_config::theta_grids=q;
    }

    /* Sets the maximal number of bins used by parni instances. Returns false
    if the argument is zero. */

    bool set_grid_bins(std::size_t n)
    {
	return MC_config::set_grid_bins(n);
    }
    bool MC_config::set_grid_bins(std::size_t n)
    {
	return MC_config::grd_bins=n;
    }
    
    /* Sets the grid bin splitting criterium. */

    void set_grid_mode(grid_modes::type m)
    {
	MC_config::set_grid_mode(m);
    }
    void MC_config::set_grid_mode(grid_modes::type m)
    {
	MC_config::grd_mode=m;
    }

    /* Sets the pdf set name and number. Returns false if the argument is
    * invalid or the LHAPDF directory is not defined. */

    bool set_pdfs(const char* name,int n)
    {
	return MC_config::set_pdfs(name,n);
    }
    bool MC_config::set_pdfs(const char* name,int n)
    {
#if HAVE_LHAPDF_H_
	if(name!=NULL and n>0)
	{
	    if(MC_config::pdf_str!=name or MC_config::pdf_nr!=n)
	    {
		MC_config::init_pdf=false;
	    }
	    MC_config::pdf_str=name;
	    MC_config::pdf_nr=n;
	    return true;
	}
#endif
	return false;
    }

    /* Sets the multichannel threshold. Returns false if the argument is invalid
     * (not within (0,1]). */

    bool set_multichannel_threshold(double e)
    {
	return MC_config::set_multichannel_threshold(e);
    }
    bool MC_config::set_multichannel_threshold(double e)
    {
	if(e>=(double)0 and e<(double)1)
	{
	    MC_config::multi_thresh=e;
	    return true;
	}
	return false;
    }

    /* Sets the subprocess threshold. Returns false if the argument is invalid
     * (not within (0,1]). */

    bool set_subprocess_threshold(double e)
    {
	return MC_config::set_subprocess_threshold(e);
    }
    bool MC_config::set_subprocess_threshold(double e)
    {
	if(e>=(double)0 and e<(double)1)
	{
	    MC_config::proc_thresh=e;
	    return true;
	}
	return false;
    }

    /* Sets the multichannel adaptivity. */

    void set_multichannel_adaptivity(double a)
    {
	MC_config::set_multichannel_adaptivity(a);
    }
    void MC_config::set_multichannel_adaptivity(double a)
    {
	MC_config::multi_adapt=a;
    }

    /* Sets the subprocess weight adaptivity. */

    void set_subprocess_adaptivity(double a)
    {
	MC_config::set_subprocess_adaptivity(a);
    }
    void MC_config::set_subprocess_adaptivity(double a)
    {
	MC_config::proc_adapt=a;
    }

    /* Sets the partonic invariant mass sampling exponent. */

    void set_shat_exponent(double nu)
    {
	MC_config::set_shat_exponent(nu);
    }
    void MC_config::set_shat_exponent(double nu)
    {
	MC_config::nu_tau=nu;
    }

    /* Sets the massless s-channel propagator sampling exponent. */

    void set_timelike_exponent(double nu)
    {
	MC_config::set_timelike_exponent(nu);
    }
    void MC_config::set_timelike_exponent(double nu)
    {
	MC_config::nu_s=nu;
    }

    /* Sets the s-channel propagator sampling exponent for the given particle. */

    void set_timelike_exponent(const std::string& phi,double nu)
    {
	MC_config::set_timelike_exponent(phi,nu);
    }
    void MC_config::set_timelike_exponent(const std::string& phi,double nu)
    {
	MC_config::nu_s_phi[phi]=nu;
    }

    /* Sets the massless t-channel propagator sampling exponent. */

    void set_spacelike_exponent(double nu)
    {
	MC_config::set_spacelike_exponent(nu);
    }
    void MC_config::set_spacelike_exponent(double nu)
    {
	MC_config::nu_t=nu;
    }

    /* Sets the t-channel propagator sampling exponent for the given particle. */

    void set_spacelike_exponent(const std::string& phi,double nu)
    {
	MC_config::set_spacelike_exponent(phi,nu);
    }
    void MC_config::set_spacelike_exponent(const std::string& phi,double nu)
    {
	MC_config::nu_t_phi[phi]=nu;
    }

    /* Sets the auxilary multi-particle invariant mass sampling exponent. */

    void set_auxiliary_exponent(double nu)
    {
	MC_config::set_auxiliary_exponent(nu);
    }
    void MC_config::set_auxiliary_exponent(double nu)
    {
	MC_config::nu_u=nu;
    }

    /* Sets the fraction of discarded highest-weight events. */

    bool discard_weight_fraction(double eps)
    {
	return MC_config::discard_weight_fraction(eps);
    }
    bool MC_config::discard_weight_fraction(double eps)
    {
	if(eps>=(double)0 and eps<(double)1)
	{
	    MC_config::epsw=eps;
	    return true;
	}
	return false;
    }

    /* Sets the number of bins used to histogram the weights. */

    void set_weight_bins(std::size_t n)
    {
	MC_config::set_weight_bins(n);
    }
    void MC_config::set_weight_bins(std::size_t n)
    {
	MC_config::wbins=n;
    }

    /* Sets the number of Newton-Raphson iterations in massive RAMBO. Returns
     * false if the argument is invalid (zero). */

    bool set_NR_iterations(std::size_t n)
    {
	return MC_config::set_NR_iterations(n);
    }
    bool MC_config::set_NR_iterations(std::size_t n)
    {
	if(n==0)
	{
	    return false;
	}
	MC_config::NR_iters=n;
	return true;
    }

    /* Sets whether to adopt the alpha-s defined by the PDF's. */

    void set_pdf_alpha_s(bool q)
    {
	MC_config::set_pdf_alpha_s(q);
    }
    void MC_config::set_pdf_alpha_s(bool q)
    {
	MC_config::pdf_alpha_s=q;
    }

    /* Sets the pre-initialisation precision. */

    void set_pre_init_events(std::size_t n)
    {
	MC_config::set_pre_init_events(n);
    }
    void MC_config::set_pre_init_events(std::size_t n)
    {
	MC_config::init_evts=n;
    }

    /* Returns the helicity generator type. */

    helicity_generators::type helicity_generator_type()
    {
	return MC_config::helicity_generator_type();
    }
    helicity_generators::type MC_config::helicity_generator_type()
    {
	return MC_config::helgen;
    }

    /* Returns the colour generator type. */

    colour_generators::type colour_generator_type()
    {
	return MC_config::colour_generator_type();
    }
    colour_generators::type MC_config::colour_generator_type()
    {
	return MC_config::colgen;
    }
    
    /* Returns the initial state type. */
    
    initial_states::type initial_state_type()
    {
	return MC_config::initial_state_type();
    }
    initial_states::type MC_config::initial_state_type()
    {
	return MC_config::isgen;
    }

    /* Returns the phase space generator type. */

    phase_space_generators::type phase_space_generator_type()
    {
	return MC_config::phase_space_generator_type();
    }
    phase_space_generators::type MC_config::phase_space_generator_type()
    {
	return MC_config::psgen;
    }

    /* Returns whether the invariant masses are sampled backward. */

    bool backward_s_sampling()
    {
	return MC_config::backward_s_sampling();
    }
    bool MC_config::backward_s_sampling()
    {
	return (MC_config::psgen==phase_space_generators::recursive_backward_s or MC_config::psgen==phase_space_generators::recursive_backward_shat);
    }

    /* Returns whether the partonic invariant mass is sampled backward. */
    
    bool backward_shat_sampling()
    {
	return MC_config::backward_shat_sampling();
    }
    bool MC_config::backward_shat_sampling()
    {
	return (MC_config::psgen==phase_space_generators::recursive_backward_shat);
    }
    
    /* Returns the beam energy */
    
    double beam_energy(std::size_t i)
    {
	return MC_config::beam_energy(i);
    }
    double MC_config::beam_energy(std::size_t i)
    {
	if(i==1)
	{
	    return MC_config::E1;
	}
	if(i==2)
	{
	    return MC_config::E2;
	}
	return (double)0;
    }

    /* Returns the first beam energy. */

    double first_beam_energy()
    {
	return MC_config::first_beam_energy();
    }
    double MC_config::first_beam_energy()
    {
	return MC_config::E1;
    }

    /* Returns the second beam energy. */

    double second_beam_energy()
    {
	return MC_config::second_beam_energy();
    }
    double MC_config::second_beam_energy()
    {
	return MC_config::E2;
    }

    /* Returns the number of channel initialisation iterations. */

    std::size_t init_channel_iterations()
    {
	return MC_config::init_channel_iterations();
    }
    std::size_t MC_config::init_channel_iterations()
    {
	return MC_config::ch_iters;
    }
    
    /* Returns the channel initialisation iteration batch size. */
    
    std::size_t init_channel_batch()
    {
	return MC_config::init_channel_batch();
    }
    std::size_t MC_config::init_channel_batch()
    {
	return MC_config::ch_batch;
    }

    /* Returns the number of grid initialisation iterations. */

    std::size_t init_grid_iterations()
    {
	return MC_config::init_grid_iterations();
    }
    std::size_t MC_config::init_grid_iterations()
    {
	return MC_config::grd_iters;
    }
    
    /* Returns the grid initialisation iteration batch size. */
    
    std::size_t init_grid_batch()
    {
	return MC_config::init_grid_batch();
    }
    std::size_t MC_config::init_grid_batch()
    {
	return MC_config::grd_batch;
    }

    /* Returns the number of subprocess cross section estimation events. */

    std::size_t subprocess_events()
    {
	return MC_config::subprocess_events();
    }
    std::size_t MC_config::subprocess_events()
    {
	return MC_config::proc_evts;
    }

    /* Returns whether the channels are automatically adapted. */

    bool auto_channel_adapt()
    {
	return MC_config::auto_channel_adapt();
    }
    bool MC_config::auto_channel_adapt()
    {
	return (auto_ch_batch!=0);
    }
    
    /* Returns the automatic channel adaptation batch size. */
    
    std::size_t auto_channel_batch()
    {
	return MC_config::auto_channel_batch();
    }
    std::size_t MC_config::auto_channel_batch()
    {
	return MC_config::auto_ch_batch;
    }

    /* Returns whether the grids are automatically adapted. */

    bool auto_grid_adapt()
    {
	return MC_config::auto_grid_adapt();
    }
    bool MC_config::auto_grid_adapt()
    {
	return (auto_grd_batch!=0);
    }
    
    /* Returns the automatic grids adaptation batch size. */
    
    std::size_t auto_grid_batch()
    {
	return MC_config::auto_grid_batch();
    }
    std::size_t MC_config::auto_grid_batch()
    {
	return MC_config::auto_grd_batch;
    }

    /* Returns whether the subprocess weights are automatically adapted. */

    bool auto_subprocess_adapt()
    {
	return MC_config::auto_subprocess_adapt();
    }
    bool MC_config::auto_subprocess_adapt()
    {
	return (auto_proc_batch!=0);
    }
    
    /* Returns the automatic subprocess weights adaptation batch size. */
    
    std::size_t auto_subprocess_batch()
    {
	return MC_config::auto_subprocess_batch();
    }
    std::size_t MC_config::auto_subprocess_batch()
    {
	return MC_config::auto_proc_batch;
    }
    
    /* Returns the maximal number of thrown points at grid/channel
     * initialisation not passing the cuts before moving on. */
    
    std::size_t max_init_rejects()
    {
	return MC_config::max_init_rejects();
    }
    std::size_t MC_config::max_init_rejects()
    {
	return MC_config::mx_evts;
    }

    /* Returns whether the s-type propagators should be optimised with an
    * adaptive grid. */
    
    bool adaptive_s_sampling()
    {
	return MC_config::adaptive_s_sampling();
    }
    bool MC_config::adaptive_s_sampling()
    {
	return MC_config::s_grids;
    }

    /* Returns whether the t-type propagators should be optimised with an
    * adaptive grid. */
    
    bool adaptive_t_sampling()
    {
	return MC_config::adaptive_t_sampling();
    }
    bool MC_config::adaptive_t_sampling()
    {
	return MC_config::t_grids;
    }

    /* Returns whether the s_branching angles should be optimised with an
    * adaptive grid. */
    
    bool adaptive_angles()
    {
	return MC_config::adaptive_angles();
    }
    bool MC_config::adaptive_angles()
    {
	return MC_config::theta_grids;
    }

    /* Returns the maximal number of bins per 1-d grid. */

    std::size_t grid_bins()
    {
	return MC_config::grid_bins();
    }
    std::size_t MC_config::grid_bins()
    {
	return MC_config::grd_bins;
    }
    
    /* Returns the grid bin-splitting mode. */
    
    grid_modes::type grid_mode()
    {
	return MC_config::grid_mode();
    }
    grid_modes::type MC_config::grid_mode()
    {
	return grd_mode;
    }

    /* Returns the current pdf name. */

    const char* pdf_name()
    {
	return MC_config::pdf_name();
    }
    const char* MC_config::pdf_name()
    {
	return MC_config::pdf_str;
    }
    
    /* Returns the current pdf number. */
    
    int pdf_number()
    {
	return MC_config::pdf_number();
    }
    int MC_config::pdf_number()
    {
	return MC_config::pdf_nr;
    }

    /* Returns the multichannel threshold. */

    double multichannel_threshold()
    {
	return MC_config::multichannel_threshold();
    }
    double MC_config::multichannel_threshold()
    {
	return MC_config::multi_thresh;
    }

    /* Returns the subprocess threshold. */

    double subprocess_threshold()
    {
	return MC_config::subprocess_threshold();
    }
    double MC_config::subprocess_threshold()
    {
	return MC_config::proc_thresh;
    }
    
    /* Returns the multichannel adaptivity. */

    double multichannel_adaptivity()
    {
	return MC_config::multichannel_adaptivity();
    }
    double MC_config::multichannel_adaptivity()
    {
	return MC_config::multi_adapt;
    }
    
    /* Returns the subprocess adaptivity. */

    double subprocess_adaptivity()
    {
	return MC_config::subprocess_adaptivity();
    }
    double MC_config::subprocess_adaptivity()
    {
	return MC_config::proc_adapt;
    }
    
    /* Returns the partonic invariant mass sampling exponent. */
    
    double shat_exponent()
    {
	return MC_config::shat_exponent();
    }
    double MC_config::shat_exponent()
    {
	return MC_config::nu_tau;
    }

    /* Returns the timelike propagator sampling exponent. */

    double timelike_exponent()
    {
	return MC_config::timelike_exponent();
    }
    double MC_config::timelike_exponent()
    {
	return MC_config::nu_s;
    }

    /* Returns the timelike propagator sampling exponent for the given particle. */

    double timelike_exponent(const std::string& phi)
    {
	return MC_config::timelike_exponent(phi);
    }
    double MC_config::timelike_exponent(const std::string& phi)
    {
	return (MC_config::nu_s_phi.find(phi)==MC_config::nu_s_phi.end())?(MC_config::nu_s):(MC_config::nu_s_phi[phi]);
    }

    /* Returns the spacelike propagator sampling exponent. */

    double spacelike_exponent()
    {
	return MC_config::spacelike_exponent();
    }
    double MC_config::spacelike_exponent()
    {
	return MC_config::nu_t;
    }

    /* Returns the spacelike propagator sampling exponent for the given particle. */

    double spacelike_exponent(const std::string& phi)
    {
	return MC_config::spacelike_exponent(phi);
    }
    double MC_config::spacelike_exponent(const std::string& phi)
    {
	return (MC_config::nu_t_phi.find(phi)==MC_config::nu_t_phi.end())?(MC_config::nu_t):(MC_config::nu_t_phi[phi]);
    }

    /* Returns the multi-particle invariant mass sampling exponent. */

    double auxiliary_exponent()
    {
	return MC_config::auxiliary_exponent();
    }
    double MC_config::auxiliary_exponent()
    {
	return MC_config::nu_u;
    }

    /* Returns the fraction of discarded highest-weight events. */

    double discarded_weight_fraction()
    {
	return MC_config::discarded_weight_fraction();
    }
    double MC_config::discarded_weight_fraction()
    {
	return MC_config::epsw;
    }

    /* Returns the number of bins used to histogram the weights. */

    std::size_t weight_bins()
    {
	return MC_config::weight_bins();
    }
    std::size_t MC_config::weight_bins()
    {
	return MC_config::wbins;
    }
    
    /* Returns the nr of Newton-Raphson iterations in massive RAMBO. */
    
    std::size_t NR_iterations()
    {
	return MC_config::NR_iterations();
    }
    std::size_t MC_config::NR_iterations()
    {
	return MC_config::NR_iters;
    }

    /* Sets whether to adopt the alpha-s defined by the PDF's. */

    bool use_pdf_alpha_s()
    {
	return MC_config::use_pdf_alpha_s();
    }
    bool MC_config::use_pdf_alpha_s()
    {
	return MC_config::pdf_alpha_s;
    }

    /* Returns the pre-initialisation precision: */

    std::size_t pre_init_events()
    {
	return MC_config::pre_init_events();
    }
    std::size_t MC_config::pre_init_events()
    {
	return MC_config::init_evts;
    }

    /* Returns the s-pair generation mode: */

    s_pair_generation_modes::type s_pair_generation_mode()
    {
	return MC_config::s_pair_generation_mode();
    }
    s_pair_generation_modes::type MC_config::s_pair_generation_mode()
    {
	return MC_config::s_pair_genmode;
    }

    /* Sets the minimal invariant mass for particles (i,j). Returns the
       effectively inserted dimass. */

    double basic_cuts::set_m_min(int i,int j,double m)
    {
	return MC_config::set_m_min(i,j,m);
    }
    double MC_config::set_m_min(int i,int j,double m)
    {
	if(i>0 and j>0 and i!=j)
	{
	    std::pair<int,int> p1(i,j),p2(j,i);
	    MC_config::mmin[p1]=std::abs(m);
	    MC_config::mmin[p2]=std::abs(m);
	    return m;
	}
	return 0;
    }

    /* Sets the s-pair generation mode: */

    void set_s_pair_generation_mode(s_pair_generation_modes::type mode)
    {
	MC_config::set_s_pair_generation_mode(mode);
    }
    void MC_config::set_s_pair_generation_mode(s_pair_generation_modes::type mode)
    {
	MC_config::s_pair_genmode=mode;
    }

    /* Sets the minimal invariant mass-squared for particles (i,j). Returns the
       effectively inserted dimass. */

    double basic_cuts::set_s_min(int i,int j,double s)
    {
	if(s>0)
	{
	    return basic_cuts::set_m_min(i,j,std::sqrt(s));
	}
	return 0;
    }

    /* Sets all basic cut variables to their default values. */

    void basic_cuts::clear()
    {
	MC_config::clear_cuts();
    }
    void MC_config::clear_cuts()
    {
	mmin.clear();
    }

    /* Returns the minimal (i,j)-dimass. */

    double basic_cuts::m_min(int i,int j)
    {
	return MC_config::m_min(i,j);
    }
    double MC_config::m_min(int i,int j)
    {
	std::pair<int,int>p(i,j);
	return (MC_config::mmin.find(p)!=MC_config::mmin.end())?(MC_config::mmin[p]):(double)0;
    }

    helicity_generators::type MC_config::helgen=helicity_generators::uniform;
    colour_generators::type MC_config::colgen=colour_generators::flow_sampling;
    initial_states::type MC_config::isgen=initial_states::eplus_eminus;
    phase_space_generators::type MC_config::psgen=phase_space_generators::recursive;
    double MC_config::E1=104.5;
    double MC_config::E2=104.5;
    std::size_t MC_config::ch_iters=10;
    std::size_t MC_config::ch_batch=10000;
    std::size_t MC_config::grd_iters=100;
    std::size_t MC_config::grd_batch=100;
    std::size_t MC_config::proc_evts=100000;
    std::size_t MC_config::auto_ch_batch=0;
    std::size_t MC_config::auto_grd_batch=1000;
    std::size_t MC_config::auto_proc_batch=1000;
    std::size_t MC_config::mx_throws=10000;
    std::size_t MC_config::mx_evts=10000;
    bool MC_config::s_grids=false;
    bool MC_config::t_grids=false;
    bool MC_config::theta_grids=false;
    std::size_t MC_config::grd_bins=200;
    grid_modes::type MC_config::grd_mode=grid_modes::cumulant_weights;
    const char* MC_config::pdf_str="cteq6l.LHpdf";
    int MC_config::pdf_nr=1;
    bool MC_config::init_pdf=false;
    double MC_config::multi_thresh=0.01;
    double MC_config::proc_thresh=0.0001;
    double MC_config::multi_adapt=0.5;
    double MC_config::proc_adapt=0.5;
    double MC_config::nu_tau=0.75;
    double MC_config::nu_s=0.5;
    double MC_config::nu_t=0.9;
    double MC_config::nu_u=0.3;
    double MC_config::epsw=0;
    std::size_t MC_config::wbins=1000;
    std::size_t MC_config::NR_iters=10;
    bool MC_config::pdf_alpha_s=true;
    std::size_t MC_config::init_evts=10000;
    s_pair_generation_modes::type MC_config::s_pair_genmode=s_pair_generation_modes::hit_and_miss;
    
    std::map<std::pair<int,int>,double> MC_config::mmin;
    
    void MC_config::print_exponents()
    {
	std::cout<<"default timelike:   "<<nu_s<<std::endl;
	std::cout<<"default spacelike:  "<<nu_t<<std::endl;
	std::cout<<"auxiliary channels: "<<nu_u<<std::endl;
	std::cout<<"full process s-hat: "<<nu_tau<<std::endl;
	std::cout<<"specialised timelike values:"<<std::endl;
	for(std::map<std::string,double>::const_iterator it=nu_s_phi.begin();it!=nu_s_phi.end();++it)
	{
	    std::cout<<it->first<<"\t\t"<<it->second<<std::endl;
	}
	std::cout<<"specialised spacelike values:"<<std::endl;
	for(std::map<std::string,double>::const_iterator it=nu_t_phi.begin();it!=nu_t_phi.end();++it)
	{
	    std::cout<<it->first<<"\t\t"<<it->second<<std::endl;
	}
    }

    std::map<std::string,double> MC_config::make_nu_s_phi()
    {
	std::map<std::string,double> result;
	result["gamma"]=0.9;
	result["g"]=0.9;
	return result;
    }
    std::map<std::string,double> MC_config::make_nu_t_phi()
    {
	std::map<std::string,double> result;
	result["gamma"]=0.9;
	result["g"]=0.9;
	return result;
    }
    std::map<std::string,double> MC_config::nu_s_phi=MC_config::make_nu_s_phi();
    std::map<std::string,double> MC_config::nu_t_phi=MC_config::make_nu_t_phi();
}

