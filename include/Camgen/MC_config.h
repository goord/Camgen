//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file MC_config.h
    \brief Static configuration data for Camgen Monte Carlo module.
 */

#ifndef CAMGEN_MC_CONFIG_H_
#define CAMGEN_MC_CONFIG_H_

#include <map>
#include <string>
#include <cstddef>

namespace Camgen
{
    /// Helicity generator types.

    struct helicity_generators
    {
	enum type
	{
	    unused,
	    uniform,
	    longitudinal,
	    summation
	};
    };

    /// Colour generator types.
	
    struct colour_generators
    {
	enum type
	{
	    unused,
	    uniform,
	    summation,
	    flow_sampling
	};
    };
    
    /// Initial state generator types.

    struct initial_states
    {
	enum type
	{
	    unused,
	    partonic,
	    eplus_eminus,
	    proton_proton,
	    proton_antiproton,
	    antiproton_proton,
	    antiproton_antiproton,
	    proton_proton_xx,
	    proton_antiproton_xx,
	    antiproton_proton_xx,
	    antiproton_antiproton_xx
	};
    };

    /// Momentum generator types.

    struct phase_space_generators
    {
	enum type
	{
	    unused,
	    uniform,
	    recursive,
	    recursive_backward_s,
	    recursive_backward_shat
	};
    };

    /// Adaptive grid modes.

    struct grid_modes
    {
	enum type
	{
	    cumulant_weights,
	    variance_weights,
	    maximum_weights
	};
    };

    /// S-pair generation modes.
    
    struct s_pair_generation_modes
    {
	enum type
	{
	    asymmetric,
	    symmetric,
	    hit_and_miss
	};
    };

    /// Sets the default helicity generator.

    void set_helicity_generator_type(helicity_generators::type);

    /// Sets the default colour generator.

    void set_colour_generator_type(colour_generators::type);

    /// Sets the default initial state.

    void set_initial_state_type(initial_states::type);

    /// Sets the default phase space generator.

    void set_phase_space_generator_type(phase_space_generators::type);

    /// Sets the i-th beam energy, where label starts from 1,...

    bool set_beam_energy(int,double);

    /// Sets the first beam energy. Returns false at failure (negative argument).

    bool set_first_beam_energy(double);

    /// Sets the second beam energy. Returns false at failure (negative
    /// argument).
    
    bool set_second_beam_energy(double);

    /// Sets the number of channel iterations and batch size at initialisation.

    void set_channel_init(std::size_t,std::size_t);

    /// Sets the number of grid iterations and batch size at initialisation.

    void set_grid_init(std::size_t,std::size_t);

    /// Sets the number of events for subprocess cross section estimation.

    void set_subprocess_events(std::size_t);

    /// Sets the automatic channel adaptation batch size. If set zero, automatic
    /// channel adaptation is switched off.

    void set_auto_channel_adapt(std::size_t);

    /// Sets the automatic grid adaptation batch size. If set zero, automatic grid
    /// adaptation is switched off.

    void set_auto_grid_adapt(std::size_t);

    /// Sets the automatic subprocess weight adaptation batch size. If set zero,
    /// automatic subprocess weights adaption is switched off.

    void set_auto_subprocess_adapt(std::size_t);

    /// Sets the maximal number of events thrown not passing cuts before moving
    /// on at initialisation of channels and grids.

    bool set_max_init_rejects(std::size_t);

    /// Controls whether adaptive grids should optimise the s-sampling
    /// algorithms. If set true, the one-shot method is used for the generation
    /// of s-pairs.

    void set_adaptive_s_sampling(bool);

    /// Controls whether adaptive grids should optimise the t-sampling
    /// algorithms.

    void set_adaptive_t_sampling(bool);

    /// Controls whether adaptive grids should optimise the angle sampling in
    /// s-branchings.

    void set_adaptive_angles(bool);

    /// Sets the maximal number of bins used by parni instances. Returns false
    /// if the argument is zero.

    bool set_grid_bins(std::size_t);
    
    /// Sets the grid bin splitting criterium.

    void set_grid_mode(grid_modes::type);

    /// Sets the pdf set name and number. Returns false if the argument is
    /// invalid or the LHAPDF directory is not defined.

    bool set_pdfs(const char*,int);

    /// Sets the multichannel threshold. Returns false if the argument is
    /// invalid (not within (0,1]).

    bool set_multichannel_threshold(double);

    /// Sets the subprocess cross section threshold. Returns false if the
    /// argument is invalid (not within (0,1]).

    bool set_subprocess_threshold(double);

    /// Sets the multichannel adaptivity.

    void set_multichannel_adaptivity(double);

    /// Sets the subprocess weight adaptivity.

    void set_subprocess_adaptivity(double);

    /// Sets the partonic invariant mass sampling exponent.

    void set_shat_exponent(double);

    /// Sets the massless s-channel propagator sampling exponent.

    void set_timelike_exponent(double);

    /// Sets the s-channel propagator sampling exponent for the
    /// specific particle.

    void set_timelike_exponent(const std::string&,double);

    /// Sets the massless t-channel propagator sampling exponent.

    void set_spacelike_exponent(double);

    /// Sets the t-channel propagator exponent for the specific particle.

    void set_spacelike_exponent(const std::string&,double);

    /// Sets the auxilary multi-particle invariant mass sampling exponent.

    void set_auxiliary_exponent(double);

    /// Sets the maximal weight stabiliser. Returns false if the argument is not
    /// within [0,1].

    bool discard_weight_fraction(double);

    /// Sets the number of bins used to histogram event weights.

    void set_weight_bins(std::size_t);

    /// Sets the number of Newton-Raphson iterations in massive RAMBO. Returns false if
    /// the argument is invalid (zero).

    bool set_NR_iterations(std::size_t);

    /// Sets whether to adopt the alpha-s defined by the PDF's.

    void set_pdf_alpha_s(bool);

    /// Sets the pre-initialisation precision.

    void set_pre_init_events(std::size_t);

    /// Returns the helicity generator type.

    helicity_generators::type helicity_generator_type();

    /// Returns the colour generator type.

    colour_generators::type colour_generator_type();
    
    /// Returns the initial state type.
    
    initial_states::type initial_state_type();

    /// Returns the phase space generator type.

    phase_space_generators::type phase_space_generator_type();

    /// Returns whether the invariant masses are sampled backward.

    bool backward_s_sampling();

    /// Returns whether the partonic invariant mass is sampled backward.
    
    bool backward_shat_sampling();
    
    /// Returns the beam energy.
    
    double beam_energy(std::size_t);

    /// Returns the first beam energy.

    double first_beam_energy();
    
    /// Returns the second beam energy.
    
    double second_beam_energy();

    /// Returns the number of channel initialisation iterations.

    std::size_t init_channel_iterations();
    
    /// Returns the channel initialisation iteration batch size.
    
    std::size_t init_channel_batch();

    /// Returns the number of grid initialisation iterations.

    std::size_t init_grid_iterations();
    
    /// Returns the grid initialisation batch size.
    
    std::size_t init_grid_batch();

    /// Returns the number of subprocess cross section estimation events.

    std::size_t subprocess_events();

    /// Returns whether the channels are automatically adapted.

    bool auto_channel_adapt();
    
    /// Returns the automatic channel adaptation batch size.
    
    std::size_t auto_channel_batch();

    /// Returns whether the grids are automatically adapted.

    bool auto_grid_adapt();
    
    /// Returns the automatic grid adaptation batch size.
    
    std::size_t auto_grid_batch();

    /// Returns whether the subprocess weights are automatically adapted.

    bool auto_subprocess_adapt();
    
    /// Returns the automatic subprocess weights adaptation batch size.
    
    std::size_t auto_subprocess_batch();
    
    /// Returns the maximal number of thrown points at grid/channel
    /// initialisation not passing the cuts before moving on.
    
    std::size_t max_init_rejects();

    /// Returns whether the s-type propagators should be optimised with an
    /// adaptive grid.
    
    bool adaptive_s_sampling();

    /// Returns whether the t-type propagators should be optimised with an
    /// adaptive grid.
    
    bool adaptive_t_sampling();

    /// Returns whether the s-braching angles should be optimised with an
    /// adaptive grid.
    
    bool adaptive_angles();

    /// Returns the maximal number of bins per 1-d grid.

    std::size_t grid_bins();
    
    /// Returns the grid bin-splitting mode.
    
    grid_modes::type grid_mode();

    /// Returns the current pdf name.

    const char* pdf_name();
    
    /// Returns the current pdf number.
    
    int pdf_number();

    /// Returns the multichannel threshold.

    double multichannel_threshold();
    
    /// Returns the subprocess threshold.
    
    double subprocess_threshold();

    /// Returns the multichannel adaptivity.

    double multichannel_adaptivity();

    /// Returns the subprocess weight adaptivity.

    double subprocess_adaptivity();
    
    /// Returns the partonic invariant mass sampling exponent.
    
    double shat_exponent();

    /// Returns the timelike propagator sampling exponent.

    double timelike_exponent();

    /// Returns the timelike propagator sampling exponent for the given
    /// particle.

    double timelike_exponent(const std::string&);

    /// Returns the spacelike propagator sampling exponent.

    double spacelike_exponent();

    /// Returns the spacelike propagator sampling exponent for the given
    /// particle.

    double spacelike_exponent(const std::string&);

    /// Returns the multi-particle invariant mass sampling exponent.

    double auxiliary_exponent();

    /// Returns the ignored fraction of highest-weight events.

    double discarded_weight_fraction();

    /// Returns the number of bins for weights-histogramming.

    std::size_t weight_bins();
    
    /// Returns the nr of Newton-Raphson iterations in massive RAMBO.
    
    std::size_t NR_iterations();

    /// Returns whether to use the pdf's alpha_s

    bool use_pdf_alpha_s();

    /// Returns the pre-initialisation precision.

    std::size_t pre_init_events();

    /// Returns the s-pair generation mode.

    s_pair_generation_modes::type s_pair_generation_mode();

    /// Sets the s-pair generation mode.

    void set_s_pair_generation_mode(s_pair_generation_modes::type);
    
    /// Namespace for basic cut manipulation.

    namespace basic_cuts
    {
	/// Sets the minimal invariant mass for particles (i,j). Returns the
	/// effectively inserted minimal dimass.

	double set_m_min(int i,int j,double);

	/// Sets the minimal invariant mass-squared for particles (i,j). Returns
	/// the effectively inserted minimal dimass.

	double set_s_min(int i,int j,double);

	/// Sets all basic cut variables to their default values.
	
	void clear();

	/// Returns the minimal (i,j)-dimass.

	double m_min(int i,int j);
    }

    /// Class containing static data used by Camgen routines.

    class MC_config
    {
	public:

	    /* Public static modifiers: */

	    static void set_helicity_generator_type(helicity_generators::type);
	    static void set_colour_generator_type(colour_generators::type);
	    static void set_initial_state_type(initial_states::type);
	    static void set_phase_space_generator_type(phase_space_generators::type);
	    static bool set_beam_energy(int,double);
	    static bool set_first_beam_energy(double);
	    static bool set_second_beam_energy(double);
	    static void set_channel_init(std::size_t,std::size_t);
	    static void set_grid_init(std::size_t,std::size_t);
	    static void set_subprocess_events(std::size_t);
	    static void set_auto_channel_adapt(std::size_t);
	    static void set_auto_grid_adapt(std::size_t);
	    static void set_auto_subprocess_adapt(std::size_t);
	    static bool set_max_init_rejects(std::size_t);
	    static void set_adaptive_s_sampling(bool);
	    static void set_adaptive_t_sampling(bool);
	    static void set_adaptive_angles(bool);
	    static bool set_grid_bins(std::size_t);
	    static void set_grid_mode(grid_modes::type);
	    static bool set_pdfs(const char*,int);
	    static bool set_multichannel_threshold(double);
	    static bool set_subprocess_threshold(double);
	    static void set_multichannel_adaptivity(double);
	    static void set_subprocess_adaptivity(double);
	    static void set_shat_exponent(double);
	    static void set_timelike_exponent(double);
	    static void set_timelike_exponent(const std::string&,double);
	    static void set_spacelike_exponent(double);
	    static void set_spacelike_exponent(const std::string&,double);
	    static void set_auxiliary_exponent(double);
	    static bool discard_weight_fraction(double);
	    static void set_weight_bins(std::size_t);
	    static bool set_NR_iterations(std::size_t);
	    static void set_pdf_alpha_s(bool);
	    static void set_pre_init_events(std::size_t);
	    static void set_s_pair_generation_mode(s_pair_generation_modes::type);

	    /* Public static readout functions: */

	    static helicity_generators::type helicity_generator_type();
	    static colour_generators::type colour_generator_type();
	    static initial_states::type initial_state_type();
	    static phase_space_generators::type phase_space_generator_type();
	    static bool backward_s_sampling();
	    static bool backward_shat_sampling();
	    static double beam_energy(std::size_t);
	    static double first_beam_energy();
	    static double second_beam_energy();
	    static std::size_t init_channel_iterations();
	    static std::size_t init_channel_batch();
	    static std::size_t init_grid_iterations();
	    static std::size_t init_grid_batch();
	    static std::size_t subprocess_events();
	    static bool auto_channel_adapt();
	    static std::size_t auto_channel_batch();
	    static bool auto_grid_adapt();
	    static std::size_t auto_grid_batch();
	    static bool auto_subprocess_adapt();
	    static std::size_t auto_subprocess_batch();
	    static std::size_t max_init_rejects();
	    static bool adaptive_s_sampling();
	    static bool adaptive_t_sampling();
	    static bool adaptive_angles();
	    static std::size_t grid_bins();
	    static grid_modes::type grid_mode();
	    static const char* pdf_name();
	    static int pdf_number();
	    static double multichannel_threshold();
	    static double subprocess_threshold();
	    static double multichannel_adaptivity();
	    static double subprocess_adaptivity();
	    static double shat_exponent();
	    static double timelike_exponent();
	    static double timelike_exponent(const std::string&);
	    static double spacelike_exponent();
	    static double spacelike_exponent(const std::string&);
	    static double auxiliary_exponent();
	    static double discarded_weight_fraction();
	    static std::size_t weight_bins();
	    static std::size_t NR_iterations();
	    static bool use_pdf_alpha_s();
	    static std::size_t pre_init_events();
	    static s_pair_generation_modes::type s_pair_generation_mode();

	    /* Basic cut inserters: */

	    static double set_m_min(int,int,double);

	    /* Basic cut readout: */
	    
	    static double m_min(int,int);

	    /* Basic cut resetter: */

	    static void clear_cuts();
	    
	    static void print_exponents();

	private:

	    /* Static data: */

	    static helicity_generators::type helgen;
	    static colour_generators::type colgen;
	    static initial_states::type isgen;
	    static phase_space_generators::type psgen;
	    static double E1;
	    static double E2;
	    static std::size_t ch_iters;
	    static std::size_t ch_batch;
	    static std::size_t grd_iters;
	    static std::size_t grd_batch;
	    static std::size_t proc_evts;
	    static std::size_t auto_ch_batch;
	    static std::size_t auto_grd_batch;
	    static std::size_t auto_proc_batch;
	    static std::size_t mx_throws;
	    static std::size_t mx_evts;
	    static bool s_grids;
	    static bool t_grids;
	    static bool theta_grids;
	    static std::size_t grd_bins;
	    static grid_modes::type grd_mode;
	    static const char* pdf_str;
	    static int pdf_nr;
	    static bool init_pdf;
	    static double multi_thresh;
	    static double proc_thresh;
	    static double multi_adapt;
	    static double proc_adapt;
	    static double nu_tau;
	    static double nu_s;
	    static double nu_t;
	    static double nu_u;
	    static double epsw;
	    static std::size_t wbins;
	    static std::size_t NR_iters;
	    static bool pdf_alpha_s;
	    static std::size_t init_evts;
	    static s_pair_generation_modes::type s_pair_genmode;
	    
	    static std::map<std::string,double> nu_s_phi;
	    static std::map<std::string,double> nu_t_phi;
	    static std::map<std::pair<int,int>,double> mmin;

	    static std::map<std::string,double> make_nu_s_phi();
	    static std::map<std::string,double> make_nu_t_phi();
    };
}

#endif /*CAMGEN_MC_CONFIG_H_*/

