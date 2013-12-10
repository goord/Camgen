//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file proc_gen.h
    \brief Single-process event generator class template definition.
 */

#ifndef CAMGEN_PROC_GEN_H_
#define CAMGEN_PROC_GEN_H_

#include <Camgen/CM_algo.h>
#include <Camgen/helgen_fac.h>
#include <Camgen/colgen_fac.h>
#include <Camgen/psgen_fac.h>
#include <Camgen/gen_conf.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class event_generator;

    /// Single-process matrix-element event generator class.

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class process_generator: public MC_generator<typename model_t::value_type>,public ps_generator_base<model_t>,public phase_space_cut,public scale_expression<typename model_t::value_type>
    {
	friend class event_generator<model_t,N_in,N_out,rng_t>;
	typedef MC_generator<typename model_t::value_type> base_type;
	
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_t::value_type value_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef vector<value_type,model_t::dimension> momentum_type;
	    typedef std::size_t size_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename CM_algorithm<model_t,N_in,N_out>::tree_iterator CM_tree_iterator;
	    typedef typename CM_algorithm<model_t,N_in,N_out>::phase_space_type phase_space_type;
	    typedef typename ps_generator_factory<model_t,N_in,N_out,rng_t>::generator_type ps_generator_type;
	    typedef typename helicity_generator_factory<model_t,N_in,N_out,rng_t>::generator_type helicity_generator_type;
	    typedef typename colour_generator_factory<model_t,N_in,N_out,rng_t>::generator_type colour_generator_type;

	    /* Utility function: */

	    static bool accept(const value_type& f)
	    {
		if(f==f and f!=std::numeric_limits<value_type>::infinity() and !(f<(value_type)0))
		{
		    return true;
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"integrand "<<f<<" was not accepted by process generator"<<endlog;
		return false;
	    }

	    /* Public data members: */
	    /*----------------------*/

	    /// Picobarn conversion factor.
	    
	    static const value_type pb_conversion;

	    /// Utility process id.

	    const size_type id;

	    /// Symmetry factor.

	    const value_type symmetry_factor;

	    /* Public static factory methods: */
	    /*--------------------------------*/

	    /// Factory method reading a generator instance from the input
	    /// stream.
	    
	    static process_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,std::istream& is)
	    {
		std::string line;
		while(line!="<procgen>" and !is.eof())
		{
		    std::getline(is,line);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<"end of file reached before process generator flag was read"<<endlog;
		    return NULL;
		}
		process_generator<model_t,N_in,N_out,rng_t>* result=create_instance_noflags(algo,is,false);
		while(line!="</procgen>" and !is.eof())
		{
		    std::getline(is,line);
		}
		return result;
	    }

	    /// Factory method reading a generator instance from the input
	    /// file.
	    
	    static process_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,const std::string& filename)
	    {
		std::ifstream ifs;
		ifs.open(filename.c_str());
		if(!ifs.is_open())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"error opening input file "<<filename<<"--returning NULL."<<endlog;
		    return NULL;
		}
		return read(algo,ifs);
	    }

	    /// Factory method reading a generator instance from the input
	    /// stream and copying the cuts and scale expression from the
	    /// settings.

	    static process_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,generator_configuration<model_t>& settings,std::istream& is)
	    {
		process_generator<model_t,N_in,N_out,rng_t>* result=read(algo,is);
		result->insert_cut(&settings);
		result->insert_scale(&settings);
		return result;
	    }

	    /// Factory method reading a generator instance from the input
	    /// stream and copying the cuts and scale expression from the
	    /// settings.

	    static process_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,generator_configuration<model_t>& settings,const std::string& filename)
	    {
		std::ifstream ifs;
		ifs.open(filename.c_str());
		if(!ifs.is_open())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"error opening input file "<<filename<<"--returning NULL."<<endlog;
		    return NULL;
		}
		return read(algo,settings,ifs);
	    }

	    /* Public constructor: */
	    /*---------------------*/

	    /// Constructor configuring with static configurations settings.

	    process_generator(CM_tree_iterator it):id(0),symmetry_factor(it->symmetry_factor()),evt_counter(0),pos_evt_counter(0),tot_weight(0),alpha(1),amplitude(it),zero_me(it->count_diagrams()==(long long unsigned)0),me(1),ps_gen(NULL),ps_weight(1),ps_factor(1),hel_gen(NULL),hel_weight(1),hel_factor(1),col_gen(NULL),col_weight(1),col_factor(1),update_counter(0),auto_update(false),grid_adaptations(0),auto_grid_adapt(0),channel_adaptations(0),auto_channel_adapt(0),max_rejects(std::numeric_limits<size_type>::max()),pre_init_evts(0),ps_cut(NULL),scale(NULL)
	    {
		summed_spins.reset();
		summed_colours.reset();
		configure();
	    }

	    /// Constructor configuring with static configurations settings.
	    /// Takes the current subprocess in algo.

	    process_generator(CM_algorithm<model_t,N_in,N_out>& algo):id(0),symmetry_factor(algo.get_tree_iterator()->symmetry_factor()),evt_counter(0),pos_evt_counter(0),tot_weight(0),alpha(1),amplitude(algo.get_tree_iterator()),zero_me(amplitude->count_diagrams()==(long long unsigned)0),me(1),ps_gen(NULL),ps_weight(1),ps_factor(1),hel_gen(NULL),hel_weight(1),hel_factor(1),col_gen(NULL),col_weight(1),col_factor(1),update_counter(0),auto_update(false),grid_adaptations(0),auto_grid_adapt(0),channel_adaptations(0),auto_channel_adapt(0),max_rejects(std::numeric_limits<size_type>::max()),pre_init_evts(0),ps_cut(NULL),scale(NULL)
	    {
		summed_spins.reset();
		summed_colours.reset();
		configure();
	    }

	    /// Constructor configuring with configuration class settings.
	    
	    process_generator(CM_tree_iterator it,generator_configuration<model_t>& settings):id(0),symmetry_factor(it->symmetry_factor()),evt_counter(0),pos_evt_counter(0),tot_weight(0),alpha(1),amplitude(it),zero_me(it->count_diagrams()==(long long unsigned)0),me(1),ps_gen(NULL),ps_weight(1),ps_factor(1),hel_gen(NULL),hel_weight(1),hel_factor(1),col_gen(NULL),col_weight(1),col_factor(1),update_counter(0),auto_update(false),grid_adaptations(0),auto_grid_adapt(0),channel_adaptations(0),auto_channel_adapt(0),max_rejects(std::numeric_limits<size_type>::max()),pre_init_evts(0),ps_cut(NULL),scale(NULL)
	    {
		summed_spins.reset();
		summed_colours.reset();
		configure(settings);
	    }

	    /// Constructor configuring with configuration class settings.
	    /// Takes the current subprocess in algo.
	    
	    process_generator(CM_algorithm<model_t,N_in,N_out>& algo,generator_configuration<model_t>& settings):id(0),symmetry_factor(algo.get_tree_iterator()->symmetry_factor()),evt_counter(0),pos_evt_counter(0),tot_weight(0),alpha(1),amplitude(algo.get_tree_iterator()),zero_me(amplitude->count_diagrams()==(long long unsigned)0),me(1),ps_gen(NULL),ps_weight(1),ps_factor(1),hel_gen(NULL),hel_weight(1),hel_factor(1),col_gen(NULL),col_weight(1),col_factor(1),update_counter(0),auto_update(false),grid_adaptations(0),auto_grid_adapt(0),channel_adaptations(0),auto_channel_adapt(0),max_rejects(std::numeric_limits<size_type>::max()),pre_init_evts(0),ps_cut(NULL),scale(NULL)
	    {
		summed_spins.reset();
		summed_colours.reset();
		configure(settings);
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /// Destructor.
	    
	    virtual ~process_generator()
	    {
		if(ps_gen!=NULL)
		{
		    delete ps_gen;
		}
		if(hel_gen!=NULL)
		{
		    delete hel_gen;
		}
		if(col_gen!=NULL)
		{
		    delete col_gen;
		}
	    }

	    /* Public modifying member functions */
	    /*-----------------------------------*/

	    /// Adopts parameters and generator instances from the static
	    /// configuration data.
	    
	    void configure()
	    {
		set_helicity_generator();
		set_colour_generator();
		set_phase_space_generator();
		if(Camgen::auto_channel_adapt())
		{
		    set_auto_channel_adapt(auto_channel_batch());
		}
		if(Camgen::auto_grid_adapt())
		{
		    set_auto_grid_adapt(auto_grid_batch());
		}
		max_rejects=max_init_rejects();
		if(discarded_weight_fraction()>(value_type)0 and discarded_weight_fraction()<(value_type)1)
		{
		    this->bin_weights(weight_bins());
		    this->reduce_max_weight(discarded_weight_fraction());
		}
		set_pre_init_events(Camgen::pre_init_events());
	    }

	    /// Adopts parameters and generator instances from the configuration
	    /// object.
	    
	    void configure(generator_configuration<model_t>& settings)
	    {
		settings.lock_generator(this);
		basic_cuts::clear();
		settings.configure();
		configure();
		insert_scale(&settings);
		insert_cut(&settings);
	    }

	    /// Pre-initialisation method. Used by event_generator class to get
	    /// a first estimation of the cross section.

	    void pre_initialise(bool verbose=false)
	    {
		if(verbose)
		{
		    std::stringstream ss;
		    print_process(ss);
		    std::cout<<"pre-init subprocess "<<ss.str();
		    std::cout.flush();
		}
		size_type batch=pre_init_evts/10;
		for(size_type i=0;i<pre_init_evts;++i)
		{
		    throw_event();
		    this->refresh_cross_section();
		    if((i%batch==0) and verbose)
		    {
			std::cout<<'.';
			std::cout.flush();
		    }
		}
		if(verbose)
		{
		    std::cout<<this->cross_section()<<std::endl;
		}
	    }

	    /// Initialises with number of channel and grid iterations/batch
	    /// sizes defined in the static configuration.
	    
	    void initialise(bool verbose=false)
	    {
		initialise(init_channel_iterations(),init_channel_batch(),init_grid_iterations(),init_grid_batch(),verbose);
	    }

	    /// Initialises with the number of channel and grid iterations/batch
	    /// sizes defined in the configuration object.

	    void initialise(generator_configuration<model_t>& settings,bool verbose=false)
	    {
		settings.configure();
		initialise(verbose);
	    }

	    /// Initialises with channel_iter multichannel iterations of batch
	    /// size channel_batch and subsequently with grid_iters grid
	    /// adaptations of batch size grid_batch.

	    void initialise(size_type channel_iters,size_type channel_batch,size_type grid_iters,size_type grid_batch,bool verbose=false)
	    {
		if(zero_me)
		{
		    size_type ntot=channel_iters*channel_batch+grid_iters*grid_batch;
		    for(size_type i=0;i<ntot;++i)
		    {
			generate();
		    }
		    if(verbose)
		    {
			std::cout<<"no adaptations performed for zero matrix element."<<std::endl;
		    }
		    return;
		}
		size_type batch=channel_iters*channel_batch/10;
		size_type count=0;
		if(verbose)
		{
		    std::cout<<"adapting ps channels";
		    std::cout.flush();
		}
		for(size_type i=0;i<channel_iters;++i)
		{
		    for(size_type j=0;j<channel_batch;++j)
		    {
			throw_pos_weight_event();
			if((++count)%batch==0 and verbose)
			{
			    std::cout<<'.';
			    std::cout.flush();
			}
		    }
		    adapt_channels();
		}
		if(verbose)
		{
		    std::cout<<"done"<<std::endl;
		}
		batch=grid_iters*grid_batch/10;
		count=0;
		if(verbose)
		{
		    std::cout<<"adapting grids";
		    std::cout.flush();
		}
		for(size_type i=0;i<grid_iters;++i)
		{
		    for(size_type j=0;j<grid_batch;++j)
		    {
			throw_pos_weight_event();
			if((++count)%batch==0 and verbose)
			{
			    std::cout<<'.';
			    std::cout.flush();
			}
		    }
		    adapt_grids();
		}
		if(verbose)
		{
		    std::cout<<"done"<<std::endl;
		}
	    }

	    /// Sets the helicity generator type.

	    helicity_generator_type* set_helicity_generator(helicity_generators::type t=Camgen::helicity_generator_type())
	    {
		if(hel_gen!=NULL)
		{
		    delete hel_gen;
		}
		hel_gen=helicity_generator_factory<model_t,N_in,N_out,rng_t>::create_instance(amplitude,t);
		hel_factor=(hel_gen==NULL)?(value_type)1:(hel_gen->averaging_factor());
		if(hel_gen!=NULL)
		{
		    for(size_type i=0;i<(N_in+N_out);++i)
		    {
			summed_spins[i]=hel_gen->sum_helicity(i);
		    }
		}
		return hel_gen;
	    }

	    /// Sets the colour generator type.

	    colour_generator_type* set_colour_generator(colour_generators::type t=Camgen::colour_generator_type())
	    {
		if(col_gen!=NULL)
		{
		    delete col_gen;
		}
		col_gen=colour_generator_factory<model_t,N_in,N_out,rng_t>::create_instance(amplitude,t);
		col_factor=(col_gen==NULL)?(value_type)1:(col_gen->averaging_factor());
		if(col_gen!=NULL)
		{
		    for(size_type i=0;i<(N_in+N_out);++i)
		    {
			summed_colours[i]=col_gen->sum_colour(i);
		    }
		}
		return col_gen;
	    }

	    /// Sets the phase space generator type.

	    ps_generator_type* set_phase_space_generator(initial_states::type t1=Camgen::initial_state_type(),phase_space_generators::type t2=Camgen::phase_space_generator_type())
	    {
		if(ps_gen!=NULL)
		{
		    delete ps_gen;
		}
		ps_gen=ps_generator_factory<model_t,N_in,N_out,rng_t>::create_instance(amplitude,t1,t2);
		return ps_gen;
	    }

	    /// Sets the i-th beam energy.

	    bool set_beam_energy(int i,const value_type& E)
	    {
		if(ps_gen!=NULL)
		{
		    return ps_gen->set_beam_energy(i,E);
		}
		return false;
	    }

	    /// Sets the minimal ij-dimass cut.

	    bool set_m_min(int i,int j,const value_type& sqrts)
	    {
		if(ps_gen!=NULL)
		{
		    return ps_gen->set_m_min(i,j,sqrts);
		}
		return false;
	    }

	    /// Inserts a phase space cut.

	    void insert_cut(phase_space_cut* cut)
	    {
		ps_cut=cut;
		if(ps_gen!=NULL)
		{
		    ps_gen->insert_cut(cut);
		}
	    }

	    /// Inserts a scale expression.

	    void insert_scale(scale_expression<value_type>* expr)
	    {
		scale=expr;
		if(ps_gen!=NULL)
		{
		    ps_gen->insert_scale(expr);
		}
	    }

	    /// Refreshes minimal invariant mass generation parameters.

	    bool refresh_m_min()
	    {
		if(ps_gen!=NULL)
		{
		    for(size_type i=1;i<=N_out;++i)
		    {
			for(size_type j=i+1;j<=N_out;++j)
			{
			    ps_gen->set_m_min(i,j,basic_cuts::m_min(i,j));
			}
		    }
		    return ps_gen->refresh_m_min();
		}
		return true;
	    }

	    /// Refreshes the collider energy.

	    bool refresh_Ecm()
	    {
		if(ps_gen!=NULL)
		{
		    return ps_gen->refresh_Ecm();
		}
		return true;
	    }

	    /// Generates momenta, where the cuts are determined by the argument class.

	    bool generate_momenta()
	    {
		if(ps_gen!=NULL)
		{
		    if(ps_gen->generate())
		    {
			ps_factor=ps_gen->integrand();
			ps_weight=ps_gen->weight();
			return true;
		    }
		    else
		    {
			ps_factor=ps_gen->integrand();
			ps_weight=(value_type)0;
			return false;
		    }
		}
		else
		{
		    ps_weight=(value_type)1;
		}
		return true;
	    }

	    /// Evaluates the phase space weight, where the cuts are determined by the
	    /// argument class.

	    bool evaluate_phase_space_weight()
	    {
		if(ps_gen!=NULL)
		{
		    if(ps_gen->evaluate_weight())
		    {
			ps_factor=ps_gen->integrand();
			ps_weight=ps_gen->weight();
			return true;
		    }
		    else
		    {
			ps_factor=ps_gen->integrand();
			ps_weight=(value_type)0;
			return false;
		    }
		}
		else
		{
		    ps_weight=(value_type)1;
		}
		return true;
	    }

	    /// Generates helicities.
	    
	    bool generate_helicities()
	    {
		if(hel_gen!=NULL)
		{
		    if(hel_gen->generate())
		    {
			hel_weight=hel_gen->weight();
			return true;
		    }
		    else
		    {
			hel_weight=(value_type)0;
			return false;
		    }
		}
		else
		{
		    hel_weight=(value_type)1;
		}
		return true;
	    }

	    /// Evaluates the helicity weight.
	    
	    bool evaluate_helicity_weight()
	    {
		if(hel_gen!=NULL)
		{
		    if(hel_gen->evaluate_weight())
		    {
			hel_weight=hel_gen->weight();
			return true;
		    }
		    else
		    {
			hel_weight=(value_type)0;
			return false;
		    }
		}
		else
		{
		    hel_weight=(value_type)1;
		}
		return true;
	    }

	    /// Generates colours

	    bool generate_colours()
	    {
		if(col_gen!=NULL)
		{
		    if(col_gen->generate())
		    {
			col_weight=col_gen->weight();
			return true;
		    }
		    else
		    {
			col_weight=(value_type)0;
			return false;
		    }
		}
		else
		{
		    col_weight=(value_type)1;
		}
		return true;
	    }

	    /// Evaluates the colour weight.
	    
	    bool evaluate_colour_weight()
	    {
		if(col_gen!=NULL)
		{
		    if(col_gen->evaluate_weight())
		    {
			col_weight=col_gen->weight();
			return true;
		    }
		    else
		    {
			col_weight=(value_type)0;
			return false;
		    }
		}
		else
		{
		    col_weight=(value_type)1;
		}
		return true;
	    }

	    /// Evaluates the integrand.

	    void evaluate_amplitude()
	    {
		amplitude->reset();
		if(alpha_pdf)
		{
		    value_type as=ps_gen->alpha_s();
		    if(as!=value_type(-1))
		    {
			model_type::set_alpha_s(as);
		    }
		    else if(scale!=NULL)
		    {
			model_type::set_QCD_scale(scale->R_scale());
		    }
		}
		else if(scale!=NULL)
		{
		    model_type::set_QCD_scale(scale->R_scale());
		}
		if(summed_spins.any() or summed_colours.any())
		{
		    me=amplitude->evaluate(summed_spins,summed_colours);
		    return;
		}
		me=std::norm(amplitude->evaluate());
	    }

	    /// Generation method implementation, where the cuts are imposed by the
	    /// argument.

	    bool generate()
	    {
		throw_event();
		if(auto_update)
		{
		    update();
		    this->refresh_cross_section();
		}
		if(update_counter==0)
		{
		    return true;
		}
		if(auto_grid_adapt!=0)
		{
		    if(update_counter%auto_grid_adapt==0)
		    {
			adapt_grids();
		    }
		}
		if(auto_channel_adapt!=0)
		{
		    if(update_counter%auto_channel_adapt==0)
		    {
			adapt_channels();
		    }
		}
		return true;
	    }

	    /// Unweighted generation with cut object.

	    void generate_unweighted(bool verbose=false)
	    {
		if(verbose)
		{
		    std::cout<<"generating event for ";
		    print_process(std::cout);
		    std::cout<<"......";
		    std::cout.flush();
		}
		if(zero_me)
		{
		    generate();
		    return;
		}
		value_type rho=std::numeric_limits<value_type>::infinity();
		do
		{
		    generate();
		    if(!this->valid(tot_weight))
		    {
			continue;
		    }
		    rho=rn_stream::throw_number((value_type)0,this->max_w_eps);
		}
		while(rho>tot_weight);
		tot_weight=this->cross_section().value;
		this->refresh_max_weight();
		this->weight()=tot_weight;
		this->integrand()=(value_type)1;
		if(verbose)
		{
		    std::cout<<".....done."<<std::endl;
		}
	    }

	    /* Generates new event according to the given strategy argument: */

	    bool next_event(int strategy)
	    {
		if(std::abs(strategy)!=3)
		{
		    return generate();
		}
		generate_unweighted();
		return true;
	    }

	    /// Weight evaluation method, where the cuts are imposed by the argument.

	    bool evaluate_weight()
	    {
		if(!(evaluate_phase_space_weight() and evaluate_helicity_weight() and evaluate_colour_weight()))
		{
		    this->weight()=(value_type)0;
		    set_integrands(0);
		    return false;
		}
		this->weight()=ps_weight*hel_weight*col_weight;
		value_type f=pb_conversion*symmetry_factor*ps_factor*hel_factor*col_factor;
		if(f!=(value_type)0)
		{
		    evaluate_amplitude();
		    if(!accept(me))
		    {
			this->weight()=(value_type)0;
			set_integrands(0);
			return false;
		    }
		    f*=me;
		}
		set_integrands(f);
		if(auto_update)
		{
		    update();
		    this->refresh_cross_section();
		}
		if(update_counter==0)
		{
		    return true;
		}
		if(auto_grid_adapt!=0)
		{
		    if(update_counter%auto_grid_adapt==0)
		    {
			adapt_grids();
		    }
		}
		if(auto_channel_adapt!=0)
		{
		    if(update_counter%auto_channel_adapt==0)
		    {
			adapt_channels();
		    }
		}
		return true;
	    }

	    /// Implementation of the pass_cuts method.

	    bool pass()
	    {
		return (ps_gen==NULL)?false:(ps_gen->pass());
	    }

	    /// Updates the generators.

	    void update()
	    {
		bool q=false;
		if(ps_gen!=NULL)
		{
		    if(accept(ps_gen->integrand()))
		    {
			ps_gen->update();
			q=true;
		    }
		}
		if(hel_gen!=NULL)
		{
		    if(accept(hel_gen->integrand()))
		    {
			hel_gen->update();
			q=true;
		    }
		}
		if(col_gen!=NULL)
		{
		    if(accept(col_gen->integrand()))
		    {
			col_gen->update();
			q=true;
		    }
		}
		if(q)
		{
		    ++update_counter;
		}
	    }

	    /// Adapts the momentum generator grids.

	    void adapt_grids()
	    {
		if(ps_gen!=NULL)
		{
		    ps_gen->adapt_grids();
		    ++grid_adaptations;
		}
	    }

	    /// Adapts the momentum generator multichannels.

	    void adapt_channels()
	    {
		if(ps_gen!=NULL)
		{
		    ps_gen->adapt_channels();
		    ++channel_adaptations;
		}
	    }

	    /// Adapts the generators.

	    void adapt()
	    {
		bool q;
		if(ps_gen!=NULL)
		{
		    ps_gen->adapt();
		    q=true;
		}
		if(hel_gen!=NULL)
		{
		    hel_gen->adapt();
		    q=true;
		}
		if(col_gen!=NULL)
		{
		    col_gen->adapt();
		    q=true;
		}
		if(q)
		{
		    ++grid_adaptations;
		    ++channel_adaptations;
		}
	    }

	    /// Resets cross section, multichannel weights and adaptive grids of
	    /// phase space, helicity and colour generators.

	    void reset()
	    {
		this->base_type::reset();
		if(ps_gen!=NULL)
		{
		    ps_gen->reset();
		}
		if(hel_gen!=NULL)
		{
		    hel_gen->reset();
		}
		if(col_gen!=NULL)
		{
		    col_gen->reset();
		}
		evt_counter=0;
		pos_evt_counter=0;
		update_counter=0;
		grid_adaptations=0;
		channel_adaptations=0;
	    }

	    /// Resets cross sections of phase space, helicity and colour generators.

	    void reset_cross_section()
	    {
		this->base_type::reset_cross_section();
		if(ps_gen!=NULL)
		{
		    ps_gen->reset_cross_section();
		}
		if(hel_gen!=NULL)
		{
		    hel_gen->reset_cross_section();
		}
		if(col_gen!=NULL)
		{
		    col_gen->reset_cross_section();
		}
	    }

	    /// Refreshes phase space generator's internal parameters.

	    bool refresh_params()
	    {
		zero_me=(amplitude->count_diagrams()==(long long unsigned)0);
		if(ps_gen!=NULL)
		{
		    return ps_gen->refresh_params();
		}
		return true;
	    }

	    /// Sets the auto-update flag, so that grids, channel variances and
	    /// cross section are automatically updated every generation call.

	    void set_auto_update(bool q)
	    {
		auto_update=q;
	    }

	    /// Sets the automatic channel adaptation batch size. If nonzero, the
	    /// auto-update flag will also be set.

	    void set_auto_channel_adapt(size_type n)
	    {
		auto_channel_adapt=n;
		if(n!=0)
		{
		    set_auto_update(true);
		}
	    }

	    /// Unsets the automatic grid adaptation batch size.

	    void unset_auto_channel_adapt()
	    {
		auto_channel_adapt=0;
	    }

	    /// Sets the automatic grid adaptation batch size. If nonzero, the
	    /// auto-update flag will also be set.

	    void set_auto_grid_adapt(size_type n)
	    {
		auto_grid_adapt=n;
		if(n!=0)
		{
		    set_auto_update(true);
		}
	    }

	    /// Unsets the automatic grid adaptation batch size.

	    void unset_auto_grid_adapt()
	    {
		auto_grid_adapt=0;
	    }

	    /// Sets the optimal multichannel for the phase space generator.

	    void set_optimal_multichannel()
	    {
		if(ps_gen!=NULL)
		{
		    ps_gen->set_optimal_multichannel();
		}
	    }

	    /// Sets the pre-initialisation precision.

	    void set_pre_init_events(size_type n)
	    {
		pre_init_evts=n;
	    }

	    /// Sets whther to use the pdf alpha_s.

	    void set_pdf_alpha_s(bool q)
	    {
		alpha_pdf=q;
	    }

	    /// Collects weight and integrands

	    void collect()
	    {
		this->weight()=ps_weight*hel_weight*col_weight;
		value_type f=pb_conversion*symmetry_factor*ps_factor*hel_factor*col_factor;
		if(accept(me))
		{
		    set_integrands(me*f);
		}
		else
		{
		    set_integrands(0);
		}
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /// Returns the number of incoming particles.

	    size_type n_in() const
	    {
		return N_in;
	    }

	    /// Returns the number of outgoing particles.

	    size_type n_out() const
	    {
		return N_out;
	    }

	    /// Returns the cross section.

	    MC_integral<value_type> xsec() const
	    {
		return this->cross_section();
	    }

	    /// Returns the process id.

	    size_type process_id() const
	    {
		return id;
	    }

	    /// Returns the const reference to the i-th incoming momentum (no
	    /// range checking on i).

	    const momentum_type& p_in(size_type i) const
	    {
		return ps_gen->p_in(i);
	    }

	    /// Returns the const reference to the i-th incoming momentum (no
	    /// range checking on i).

	    const momentum_type& p_out(size_type i) const
	    {
		return ps_gen->p_out(i);
	    }

	    /// Returns the const reference to the i-th incoming mass (no
	    /// range checking on i).

	    value_type M_in(size_type i) const
	    {
		return ps_gen->M_in(i);
	    }

	    /// Returns the const reference to the i-th incoming mass (no
	    /// range checking on i).

	    value_type M_out(size_type i) const
	    {
		return ps_gen->M_out(i);
	    }

	    /// Returns the i-th incoming mass-squared (no range checking on i).

	    value_type s_in(size_type i) const
	    {
		return static_cast<ps_generator_base<model_type>*>(ps_gen)->s_in(i);
	    }

	    /// Returns the i-th outgoing mass-squared (no range checking on i).

	    value_type s_out(size_type i) const
	    {
		return ps_gen->s_out(i);
	    }
	    
	    /// Returns the i-th incoming mass (no range checking on i).

	    value_type m_in(size_type i) const
	    {
		return ps_gen->m_in(i);
	    }

	    /// Returns the i-th outgoing mass (no range checking on i).

	    value_type m_out(size_type i) const
	    {
		return ps_gen->m_out(i);
	    }

	    /// Returns the i-th beam id.

	    int beam_id(int i) const
	    {
		return ps_gen->beam_id(i);
	    }

	    /// Returns the cernlib pdf group number for the i-th incoming beam.

	    int pdfg(int i) const
	    {
		return ps_gen->pdfg(i);
	    }

	    /// Returns the cernlib pdf set number for the i-th incoming beam.

	    int pdfs(int i) const
	    {
		return ps_gen->pdfs(i);
	    }

	    /// Returns the i-th incoming particle id.

	    int id_in(size_type i) const
	    {
		return particle_in(i)->particle_type->get_pdg_id();
	    }

	    /// Returns the i-th outgoing particle id.

	    int id_out(size_type i) const
	    {
		return particle_out(i)->particle_type->get_pdg_id();
	    }

	    /// Method determining the colour connection for the event.

	    void fill_colours(std::vector<int>& c,std::vector<int>& cbar) const
	    {
		vector<int,N_in+N_out> a,abar;
		col_gen->LH_output(a,abar);
		c.assign(a.begin(),a.end());
		cbar.assign(abar.begin(),abar.end());
	    }

	    /// Returns the total hadronic CM-frame energy.

	    value_type Ecm() const
	    {
		return ps_gen->Ecm();
	    }

	    /// Returns the total hadronic invariant mass-squared.

	    value_type s_tot() const
	    {
		return ps_gen->s_tot();
	    }

	    /// Returns the total partonic CM-frame energy.

	    value_type Ecm_hat() const
	    {
		return ps_gen->Ecm_hat();
	    }

	    /// Returns the total partonic invariant mass-squared.

	    value_type s_in() const
	    {
		return ps_gen->s_in();
	    }

	    /// Returns the total event weight.

	    value_type w() const
	    {
		return tot_weight;
	    }

	    /// Returns the maximal event weight.

	    value_type max_w() const
	    {
		return this->max_weight();
	    }

	    /// Returns the event weight without flux and symmetry factors etc.
	    /// The first argument denotes whether helicities are included, the
	    /// second whether colours are.

	    value_type ps_weights(bool with_hels,bool with_cols) const
	    {
		value_type result=ps_weight;
		if(with_hels)
		{
		    result*=hel_weight;
		}
		if(with_cols)
		{
		    result*=col_weight;
		}
		return result;
	    }

	    /// Returns the flux, symmetry and conversion factors. The first
	    /// argument denotes whther helicities are included, the second
	    /// whether colours are.

	    value_type ps_factors(bool with_hels,bool with_cols) const
	    {
		value_type result=pb_conversion*symmetry_factor*ps_factor;
		if(with_hels)
		{
		    result*=hel_factor;
		}
		if(with_cols)
		{
		    result*=col_factor;
		}
		return result;

	    }

	    /// Returns the i-th beam energy.

	    value_type beam_energy(int i) const
	    {
		return ps_gen->beam_energy(i);
	    }

	    /// Returns the current factorisation scale.

	    value_type mu_F() const
	    {
		return ps_gen->mu_F();
	    }

	    /// Returns the factorisation scale.

	    value_type F_scale()
	    {
		return ps_gen->F_scale();
	    }

	    /// Returns the renormalisation scale.

	    value_type R_scale()
	    {
		return ps_gen->R_scale();
	    }

	    /// Returns the QCD scale.

	    value_type QCD_scale()
	    {
		return ps_gen->QCD_scale();
	    }

	    /// Returns the phase space generator.

	    const ps_generator_type* get_momentum_generator() const
	    {
		return ps_gen;
	    }

	    /// Returns the helicity generator.

	    const helicity_generator_type* get_helicity_generator() const
	    {
		return hel_gen;
	    }

	    /// Returns the colour generator.

	    const colour_generator_type* get_colour_generator() const
	    {
		return col_gen;
	    }

	    /// Returns the phase space weight.

	    const value_type& phase_space_weight() const
	    {
		return ps_weight;
	    }

	    /// Returns the initial state flux and phase space prefactors.

	    const value_type& phase_space_factor() const
	    {
		return ps_factor;
	    }

	    /// Returns the initial-state phase space weight.

	    value_type is_weight() const
	    {
		return (ps_gen==NULL)?0:(ps_gen->is_weight());
	    }

	    /// Returns the final-state phase space weight.

	    value_type fs_weight() const
	    {
		return (ps_gen==NULL)?0:(ps_gen->fs_weight());
	    }

	    /// Returns the helicity weight.

	    const value_type& helicity_weight() const
	    {
		return hel_weight;
	    }

	    /// Returns the helicity prefactors.

	    const value_type& helicity_factor() const
	    {
		return hel_factor;
	    }

	    /// Returns the colour weight

	    const value_type& colour_weight() const
	    {
		return col_weight;
	    }

	    /// Returns the helicity prefactors.

	    const value_type& colour_factor() const
	    {
		return col_factor;
	    }

	    /// Returns the matrix element

	    const value_type& matrix_element() const
	    {
		return me;
	    }

	    /// Returns the number of parameter updates performed since the last
	    /// adapt() call.

	    size_type updates() const
	    {
		return update_counter;
	    }

	    /// Returns the const i-th incoming particle phase space (no bound
	    //checking)

	    const phase_space_type* particle_in(size_type i) const
	    {
		return amplitude->get_phase_space(i);
	    }

	    /// Returns the const i-th incoming particle phase space (no bound
	    /// checking)

	    const phase_space_type* particle_out(size_type i) const
	    {
		return amplitude->get_phase_space(N_in+i);
	    }

	    /// Returns the i-th const particle phase space, where negative values
	    /// enumerate the incoming particles.

	    const phase_space_type* particle(int i) const
	    {
		return (i<0)?(particle_in(-i-1)):(particle_out(i-1));
	    }

	    /// Returns whether to use the pdf alpha_s.

	    bool use_pdf_alpha_s() const
	    {
		return alpha_pdf;
	    }

	    /// Returns the pre-initialisation precision.

	    value_type pre_init_events() const
	    {
		return pre_init_evts;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Prints the process being integrated.

	    std::ostream& print_process(std::ostream& os=std::cout) const
	    {
		amplitude->print_process(os);
		return os;
	    }

	    /// Calls the phase space generator printer method.

	    std::ostream& print_ps(std::ostream& os=std::cout) const
	    {
		if(ps_gen!=NULL)
		{
		    ps_gen->print(os);
		}
		return os;
	    }

	    /// Prints the status of the generator.
	    
	    std::ostream& print_status(std::ostream& os=std::cout) const
	    {
		os<<"###############################################################################################"<<std::endl;
		os<<"Nr of events generated:                            "<<std::scientific<<evt_counter<<std::endl;
		os<<"Nr of positive weight events generated:            "<<std::scientific<<pos_evt_counter<<std::endl;
		os<<"Nr of events contributing to cross-section:        "<<std::scientific<<this->calls()<<std::endl;
		os<<"Nr of updates performed:                           "<<std::scientific<<update_counter<<std::endl;
		os<<"Nr of grid adaptations performed:                  "<<std::scientific<<grid_adaptations<<std::endl;
		os<<"Nr of channel adaptations performed:               "<<std::scientific<<channel_adaptations<<std::endl;
		os<<"Monte Carlo efficiency (%):                        "<<std::scientific<<this->efficiency()<<std::endl;
		os<<"Cross section (pb):                                "<<std::scientific<<this->cross_section()<<std::endl;
		os<<"###############################################################################################"<<std::endl;
		return os;
	    }

	    /// Prints the ingredients of the returned weight.

	    std::ostream& print_weights(std::ostream& os=std::cout) const
	    {
		os<<"###############################################################################################"<<std::endl;
		os<<"Matrix element:                                    "<<std::scientific<<me<<std::endl;
		os<<"Phase space weight:                                "<<std::scientific<<ps_weight<<std::endl;
		os<<" -> initial state:					"<<std::scientific<<ps_gen->is_weight()<<std::endl;
		os<<" -> final state:					"<<std::scientific<<ps_gen->fs_weight()<<std::endl;
		os<<" -> flux factor:					"<<std::scientific<<ps_gen->flux_factor()<<std::endl;
		os<<"Helicity weight:                                   "<<std::scientific<<hel_weight<<std::endl;
		os<<"Colour weight:                                     "<<std::scientific<<col_weight<<std::endl;
		os<<"Phase space factor:                                "<<std::scientific<<ps_factor<<std::endl;
		os<<"Helicity factor:                                   "<<std::scientific<<hel_factor<<std::endl;
		os<<"Colour factor:                                     "<<std::scientific<<col_factor<<std::endl;
		os<<"Symmetry factor:                                   "<<std::scientific<<symmetry_factor<<std::endl;
		os<<"Full integrand:                                    "<<std::scientific<<this->integrand()<<std::endl;
		os<<"Full weight:                                       "<<std::scientific<<this->weight()<<std::endl;
		os<<"###############################################################################################"<<std::endl;
		return os;
	    }

	    /// Prints the generator cuts.

	    std::ostream& print_cuts(std::ostream& os=std::cout) const
	    {
		return (ps_cut==NULL)?os:(ps_cut->print(os));
	    }

	    /// Prints the generator settings.

	    std::ostream& print_settings(std::ostream& os=std::cout) const
	    {
		os<<std::setw(30)<<std::left<<"spin summmation:";
		for(size_type i=0;i<(N_in+N_out);++i)
		{
		    os<<summed_spins[i];
		}
		os<<std::endl;
		os<<std::setw(30)<<std::left<<"helicity generator:";
		if(hel_gen==NULL)
		{
		    os<<"none"<<std::endl;
		}
		else
		{
		    os<<hel_gen->type()<<std::endl;
		}
		os<<std::setw(30)<<std::left<<"colour summmation:";
		for(size_type i=0;i<(N_in+N_out);++i)
		{
		    os<<summed_colours[i];
		}
		os<<std::endl;
		os<<std::setw(30)<<std::left<<"colour generator:";
		if(col_gen==NULL)
		{
		    os<<"none"<<std::endl;
		}
		else
		{
		    os<<col_gen->type()<<std::endl;
		}
		os<<std::setw(30)<<std::left<<"momentum generator:";
		if(ps_gen==NULL)
		{
		    os<<"none"<<std::endl;
		}
		else
		{
		    os<<ps_gen->type()<<std::endl;
		    ps_gen->print_settings(os);
		}
		return os;
	    }

	    /// Overridden loading method.

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		safe_read(is,alpha);
		is>>evt_counter>>pos_evt_counter>>update_counter;
		is>>grid_adaptations>>channel_adaptations>>summed_spins>>summed_colours;
		is>>auto_update>>auto_grid_adapt>>auto_channel_adapt>>max_rejects>>alpha_pdf;
		return is;
	    }

	    /// Overridden saving method.

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<procgen>"<<std::endl;
		os<<id<<std::endl;
		print_process(os);
		os<<std::endl;
		if(ps_gen==NULL)
		{
		    os<<'N'<<std::endl;
		}
		else
		{
		    ps_gen->save(os);
		}
		if(col_gen==NULL)
		{
		    os<<'N'<<std::endl;
		}
		else
		{
		    col_gen->save(os);
		}
		if(hel_gen==NULL)
		{
		    os<<'N'<<std::endl;
		}
		else
		{
		    hel_gen->save(os);
		}
		this->base_type::save(os);
		safe_write(os,alpha);
		os<<"\t"<<evt_counter<<"\t"<<pos_evt_counter<<"\t"<<update_counter<<"\t";
		os<<grid_adaptations<<"\t"<<channel_adaptations<<std::endl;
		os<<summed_spins<<"\t"<<summed_colours;
		os<<auto_update<<"\t"<<auto_grid_adapt<<"\t"<<auto_channel_adapt<<"\t"<<max_rejects<<"\t"<<alpha_pdf<<std::endl;
		os<<"</procgen>"<<std::endl;
		return os;
	    }

	    /// Writes the object to the argument filename.

	    bool write(const std::string& filename) const
	    {
		std::ofstream ofs;
		ofs.open(filename.c_str());
		if(!ofs.is_open())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"error opening output file "<<filename<<endlog;
		    return false;
		}
		save(ofs);
		ofs.close();
		return true;
	    }

	private:

	    /* Private static methods: */
	    /*-------------------------*/

	    /* Factory method from input stream, no request for input stream
	     * process flags: */
	    
	    static process_generator<model_t,N_in,N_out,rng_t>* create_instance_noflags(CM_algorithm<model_t,N_in,N_out>& algo,std::istream& is,bool verbose=false)
	    {
		size_type id_=0;
		is>>id_>>std::ws;
		std::string procline;
		std::getline(is,procline);
		algo.set_process(procline);
		if(!algo.valid_process())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid process detected--returning NULL"<<endlog;
		    return NULL;
		}
		process_generator<model_t,N_in,N_out,rng_t>* result=new process_generator<model_t,N_in,N_out,rng_t>(algo.get_tree_iterator(),id_);
		if(verbose)
		{
		    std::cout<<"Created process generator instance..."<<std::endl;
		    std::cout<<"------------------"<<std::endl;
		    std::cout<<"Process:"<<std::endl;
		    std::cout<<"------------------"<<std::endl;
		    algo.get_tree_iterator()->print_process(std::cout);
		    std::cout<<std::endl;
		}
		is>>std::ws;
		if(is.peek()=='N')
		{
		    is.ignore(1);
		}
		else
		{
		    result->ps_gen=ps_generator_factory<model_t,N_in,N_out,rng_t>::create_instance(algo.get_tree_iterator(),is);
		}
		is>>std::ws;
		if(is.peek()=='N')
		{
		    is.ignore(1);
		}
		else
		{
		    result->col_gen=colour_generator_factory<model_t,N_in,N_out,rng_t,model_t::coloured>::create_instance(algo.get_tree_iterator(),is);
		    result->col_factor=(result->col_gen==NULL)?(value_type)1:(result->col_gen->averaging_factor());
		    if(result->col_gen!=NULL)
		    {
			for(size_type i=0;i<(N_in+N_out);++i)
			{
			    result->summed_colours[i]=result->col_gen->sum_colour(i);
			}
		    }
		}
		is>>std::ws;
		if(is.peek()=='N')
		{
		    is.ignore(1);
		}
		else
		{
		    result->hel_gen=helicity_generator_factory<model_t,N_in,N_out,rng_t>::create_instance(algo.get_tree_iterator(),is);
		    result->hel_factor=(result->hel_gen==NULL)?(value_type)1:(result->hel_gen->averaging_factor());
		    if(result->hel_gen!=NULL)
		    {
			for(size_type i=0;i<(N_in+N_out);++i)
			{
			    result->summed_spins[i]=result->hel_gen->sum_helicity(i);
			}
		    }
		}
		result->load(is);
		if(verbose)
		{
		    std::cout<<"------------------"<<std::endl;
		    std::cout<<"Settings:"<<std::endl;
		    std::cout<<"------------------"<<std::endl;
		    result->print_settings(std::cout);
		}
		return result;
	    }

	    /* Private data members */
	    /*----------------------*/

	    /* Event counter: */

	    size_type evt_counter;

	    /* Nonzero-weight event counter: */

	    size_type pos_evt_counter;

	    /* Full weight variable. */

	    value_type tot_weight;

	    /* Multichannel weight variable. */
	    
	    value_type alpha;

	    /* Camgen tree iterator for matrix-element calculation: */

	    CM_tree_iterator amplitude;

	    /* Helicity summation flags: */

	    std::bitset<N_in+N_out>summed_spins;

	    /* Colour summation flags: */

	    std::bitset<N_in+N_out>summed_colours;

	    /// Flag denoting whether the matrix element is zero.

	    bool zero_me;

	    /* Matrix element, final state symmetry factor: */

	    value_type me;

	    /* Momentum generator instance: */
	    
	    ps_generator_type* ps_gen;
	    
	    /* Phase space weight and flux factor: */
	    
	    value_type ps_weight,ps_factor;

	    /* Helicity generator instance: */

	    helicity_generator_type* hel_gen;
	    
	    /* Helicity generator weight and normalisation factor: */
	    
	    value_type hel_weight,hel_factor;

	    /* Colour generator instance: */
	    
	    colour_generator_type* col_gen;

	    /* Colour generator weight and normalisation factor: */

	    value_type col_weight,col_factor;

	    /* Update counter: */

	    size_type update_counter;

	    /* Automatic update flag: */

	    bool auto_update;

	    /* Grid adaptation counter: */

	    size_type grid_adaptations;

	    /* Automatic grid adaptation batch size: */

	    size_type auto_grid_adapt;

	    /* Channel adaptation counter: */

	    size_type channel_adaptations;

	    /* Automatic channel adaptation batch size: */

	    size_type auto_channel_adapt;

	    /* Maximal nr of rejections for positive event generation: */

	    size_type max_rejects;

	    /* Pre-initialisation precision: */

	    size_type pre_init_evts;

	    /* Phase space cuts depending on this instance: */

	    phase_space_cut* ps_cut;

	    /* QCD scale expression instance: */

	    scale_expression<value_type>* scale;

	    /* Boolean denoting whether to adopt pdf's alpha: */

	    bool alpha_pdf;

	    /* Private constructor, no configuration performed: */

	    process_generator(CM_tree_iterator it,size_type id_):id(id_),symmetry_factor(it->symmetry_factor()),evt_counter(0),pos_evt_counter(0),tot_weight(0),alpha(1),amplitude(it),zero_me(it->count_diagrams()==(long long unsigned)0),me(1),ps_gen(NULL),ps_weight(1),ps_factor(1),hel_gen(NULL),hel_weight(1),hel_factor(1),col_gen(NULL),col_weight(1),col_factor(1),update_counter(0),auto_update(false),grid_adaptations(0),auto_grid_adapt(0),channel_adaptations(0),auto_channel_adapt(0),max_rejects(std::numeric_limits<size_type>::max()),pre_init_evts(0),ps_cut(NULL),scale(NULL),alpha_pdf(Camgen::use_pdf_alpha_s())
	    {
		summed_spins.reset();
		summed_colours.reset();
	    }

	    /* Integrand definition helper: */

	    void set_integrands(const value_type& f)
	    {
		this->integrand()=f;
		tot_weight=this->weight()*f;
		if(hel_gen!=NULL)
		{
		    hel_gen->integrand()=f;
		}
		if(col_gen!=NULL)
		{
		    col_gen->integrand()=f;
		}
		if(ps_gen!=NULL)
		{
		    ps_gen->integrand()=f*hel_weight*col_weight;
		}
	    }

	    /* Event generation helper: */

	    bool throw_event()
	    {
		++evt_counter;
		if(!(generate_momenta() and generate_helicities() and generate_colours()))
		{
		    this->weight()=(value_type)0;
		    set_integrands(0);
		    return false;
		}
		else
		{
		    this->weight()=ps_weight*hel_weight*col_weight;
		    value_type f=pb_conversion*symmetry_factor*ps_factor*hel_factor*col_factor;
		    if(f!=(value_type)0)
		    {
			evaluate_amplitude();
			if(accept(me))
			{
			    set_integrands(me*f);
			}
			else
			{
			    set_integrands(0);
			}
		    }
		    else
		    {
			set_integrands(0);
		    }
		}
		if(tot_weight!=(value_type)0)
		{
		    ++pos_evt_counter;
		}
		return true;
	    }

	    /* Initialisation helper: */

	    bool throw_pos_weight_event()
	    {
		size_type n(0);
		do
		{
		    throw_event();
		    update();
		    ++n;
		}
		while(tot_weight==(value_type)0 and (n<max_rejects));
		return (tot_weight!=(value_type)0);
	    }
    };
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>const typename process_generator<model_t,N_in,N_out,rng_t>::value_type process_generator<model_t,N_in,N_out,rng_t>::pb_conversion=0.389e+09;
}

#endif /*CAMGEN_PROC_GEN_H_*/

