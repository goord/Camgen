//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_gen.h
    \brief Multi-process event generator class.
 */

#ifndef CAMGEN_EVT_GEN_H_
#define CAMGEN_EVT_GEN_H_

#include <Camgen/MC_config.h>
#include <Camgen/proc_gen.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class event_generator: public MC_generator<typename model_t::value_type>,
    												 public ps_generator_base<model_t>,
												 public phase_space_cut,
												 public scale_expression<typename model_t::value_type>
    {
	typedef MC_generator<typename model_t::value_type> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_t::value_type value_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef vector<value_type,spacetime_type::dimension> momentum_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename CM_algorithm<model_t,N_in,N_out>::tree_iterator CM_tree_iterator;
	    typedef process_generator<model_t,N_in,N_out,rng_t> process_generator_type;
	    typedef std::vector<process_generator_type*> process_container;
	    typedef typename process_container::iterator process_iterator;
	    typedef typename process_container::const_iterator const_process_iterator;
	    typedef typename process_generator_type::ps_generator_type ps_generator_type;
	    typedef typename process_generator_type::helicity_generator_type helicity_generator_type;
	    typedef typename process_generator_type::colour_generator_type colour_generator_type;
	    typedef MC_integral<value_type> cross_section_type;
	    typedef typename CM_algorithm<model_t,N_in,N_out>::phase_space_type phase_space_type;

	    /* Public static methods: */
	    /*------------------------*/

	    /// Process comparison according to ascending cross section.

	    static bool xsec_less(const process_generator_type* proc1,const process_generator_type* proc2)
	    {
		return proc1->cross_section().value<proc2->cross_section().value;
	    }

	    /// Process comparison class according to descending cross section.

	    static bool xsec_more(const process_generator_type* proc1,const process_generator_type* proc2)
	    {
		return proc1->cross_section().value>=proc2->cross_section().value;
	    }

	    /// Process comparison according to ascending calling frequency.

	    static bool alpha_less(const process_generator_type* proc1,const process_generator_type* proc2)
	    {
		return proc1->alpha<proc2->alpha;
	    }

	    /// Process comparison class according to descending calling frequence.

	    static bool alpha_more(const process_generator_type* proc1,const process_generator_type* proc2)
	    {
		return proc1->alpha>proc2->alpha;
	    }

	    /// Creates an event generator from a CM algorithm and an input
	    /// stream.
	    
	    static event_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,std::istream& is)
	    {
		if(algo.n_trees()==0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"event generator requested for empty process list--NULL returned"<<endlog;
		    return NULL;
		}
		event_generator<model_t,N_in,N_out,rng_t>* result=new event_generator<model_t,N_in,N_out,rng_t>(algo,false);
		std::string line;
		while(line!="<evtgen>" and !is.eof())
		{
		    std::getline(is,line);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before event generator flag was read"<<endlog;
		    return NULL;
		}
		result->load(is);
		while(line!="<procgen>" and !is.eof())
		{
		    std::getline(is,line);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before first process flag was read"<<endlog;
		    return NULL;
		}
		process_generator_type* p;
		while(line=="<procgen>" and !is.eof())
		{
		    p=process_generator_type::create_instance_noflags(algo,is,false);
		    if(p!=NULL)
		    {
			result->procs.push_back(p);
		    }
		    while(line!="</procgen>" and !is.eof())
		    {
			std::getline(is,line);
		    }
		    std::getline(is,line);
		}
		std::sort(result->procs.begin(),result->procs.end(),alpha_more);
		if(result->procs.size()!=0)
		{
		    result->sub_proc=result->procs.begin();
		}
		else
		{
		    result->sub_proc=result->procs.end();
		}
		return result;
	    }

	    /// Factory method reading a generator instance from the input
	    /// file.
	    
	    static event_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,const std::string& filename)
	    {
		std::ifstream ifs;
		ifs.open(filename.c_str());
		if(!ifs.is_open())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"error opening input file "<<filename<<"--NULL returned"<<endlog;
		    return NULL;
		}
		return read(algo,ifs);
	    }

	    /// Creates an event generator from a CM algorithm and an input
	    /// stream.
	    
	    static event_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,generator_configuration<model_t>& settings,std::istream& is)
	    {
		if(algo.n_trees()==0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"event generator requested for empty process list--NULL returned"<<endlog;
		    return NULL;
		}
		event_generator<model_t,N_in,N_out,rng_t>* result=new event_generator<model_t,N_in,N_out,rng_t>(algo,false);
		std::string line;
		while(line!="<evtgen>" and !is.eof())
		{
		    std::getline(is,line);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before event generator flag was read"<<endlog;
		    return NULL;
		}
		result->load(is);
		while(line!="<procgen>" and !is.eof())
		{
		    std::getline(is,line);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before first process flag was read"<<endlog;
		    return NULL;
		}
		process_generator_type* p;
		while(line=="<procgen>" and !is.eof())
		{
		    p=process_generator_type::create_instance_noflags(algo,is,false);
		    if(p!=NULL)
		    {
			p->insert_cut(&settings);
			p->insert_scale(&settings);
			result->procs.push_back(p);
		    }
		    while(line!="</procgen>" and !is.eof())
		    {
			std::getline(is,line);
		    }
		    std::getline(is,line);
		}
		std::sort(result->procs.begin(),result->procs.end(),alpha_more);
		if(result->procs.size()!=0)
		{
		    result->sub_proc=result->procs.begin();
		}
		else
		{
		    result->sub_proc=result->procs.end();
		}
		settings.lock_generator(result);
		result->ps_cut=&settings;
		result->scale=&settings;
		return result;
	    }

	    /// Factory method reading a generator instance from the input
	    /// stream and copying the cuts and scale expression from the
	    /// settings.

	    static event_generator<model_t,N_in,N_out,rng_t>* read(CM_algorithm<model_t,N_in,N_out>& algo,generator_configuration<model_t>& settings,const std::string& filename)
	    {
		std::ifstream ifs;
		ifs.open(filename.c_str());
		if(!ifs.is_open())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"error opening input file "<<filename<<"--NULL returned"<<endlog;
		    return NULL;
		}
		return read(algo,settings,ifs);
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /// Constructor from static configuration data.

	    event_generator(CM_algorithm<model_t,N_in,N_out>& algo):algorithm(algo),update_counter(0),auto_update(false),auto_proc_adapt(0),ps_cut(NULL),scale(NULL),subproc_params(1,0)
	    {
		procs.reserve(algorithm.n_trees());
		init();
		for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
		{
		    (*sub_proc)->configure();
		}
		if(auto_subprocess_adapt())
		{
		    set_auto_proc_adapt(auto_subprocess_batch());
		}
		subproc_params.first=subprocess_adaptivity();
		subproc_params.second=subprocess_threshold();
	    }

	    /// Constructor from configuration data in the settings object.

	    event_generator(CM_algorithm<model_t,N_in,N_out>& algo,generator_configuration<model_t>& settings):algorithm(algo),update_counter(0),auto_update(false),auto_proc_adapt(0),ps_cut(NULL),scale(NULL),subproc_params(1,0)
	    {
		procs.reserve(algorithm.n_trees());
		init();
		for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
		{
		    settings.set_generator(*sub_proc);
		    basic_cuts::clear();
		    settings.configure();
		    (*sub_proc)->configure();
		    (*sub_proc)->insert_scale(&settings);
		    (*sub_proc)->insert_cut(&settings);
		}
		sub_proc=procs.begin();
		settings.lock_generator(this);
		settings.configure();
		if(auto_subprocess_adapt())
		{
		    set_auto_proc_adapt(auto_subprocess_batch());
		}
		subproc_params.first=subprocess_adaptivity();
		subproc_params.second=subprocess_threshold();
	    }

	    /// Pre-initialisation method to obtain first subprocess cross
	    /// section estimates, determining which subprocesses are relevant.

	    void pre_initialise(bool verbose=false)
	    {
		for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
		{
		    (*sub_proc)->pre_initialise(verbose);
		}
		refresh_xsec();
		adapt_processes();
	    }

	    /// Initialiser method using static configuration data.

	    void initialise(bool verbose=false)
	    {
		for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
		{
		    if(verbose)
		    {
			std::stringstream ss;
			(*sub_proc)->print_process(ss);
			std::cout<<std::endl<<"init subprocess "<<ss.str()<<std::endl;
		    }
		    (*sub_proc)->initialise(verbose);
		    loop_process(subprocess_events(),verbose);
		}
		refresh_xsec();
		adapt_processes();
	    }

	    /// Initialiser method using generator settings data.

	    void initialise(generator_configuration<model_t>& settings,bool verbose=false)
	    {
		for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
		{
		    if(verbose)
		    {
			std::stringstream ss;
			(*sub_proc)->print_process(ss);
			std::cout<<std::endl<<"init subprocess "<<ss.str()<<std::endl;
		    }
		    settings.configure();
		    (*sub_proc)->initialise(verbose);
		    loop_process(subprocess_events(),verbose);
		}
		refresh_xsec();
		adapt_processes();
	    }

	    /// Initialiser method.

	    void initialise(size_type channel_iters,size_type channel_batch,size_type grid_iters,size_type grid_batch,size_type sub_proc_evts,bool verbose=false)
	    {
		for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
		{
		    if(verbose)
		    {
			std::stringstream ss;
			(*sub_proc)->print_process(ss);
			std::cout<<std::endl<<ss.str()<<std::endl;
		    }
		    (*sub_proc)->initialise(channel_iters,channel_batch,grid_iters,grid_batch,verbose);
		    loop_process(sub_proc_evts,verbose);
		}
		refresh_xsec();
		adapt_processes();
	    }

	    /// Destructor.

	    virtual ~event_generator()
	    {
		for(size_type i=0;i<procs.size();++i)
		{
		    delete procs[i];
		}
	    }

	    /* Public modifying member functions */
	    /*-----------------------------------*/

	    /// Sets the i-th beam energy.

	    bool set_beam_energy(int n,const value_type& E)
	    {
		bool q=true;
		for(size_type i=0;i<procs.size();++i)
		{
		    q&=(procs[i]->set_beam_energy(n,E));
		}
		return q;
	    }

	    /// Sets a minimal invariant dimass cut.

	    bool set_m_min(int i,int j,const value_type& sqrts)
	    {
		bool q=true;
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    q&=(*it)->set_m_min(i,j,sqrts);
		}
		return q;
	    }

	    /// Sets the auto-update flag for the event generator and all
	    /// subprocesses.

	    void set_auto_update(bool q)
	    {
		auto_update=q;
		for(size_type i=0;i<procs.size();++i)
		{
		    procs[i]->set_auto_update(q);
		}
	    }

	    /// Sets the automatic multichannel adaptation batch size for all
	    /// subprocess generators.

	    void set_auto_channel_adapt(size_type n)
	    {
		for(size_type i=0;i<procs.size();++i)
		{
		    procs[i]->set_auto_channel_adapt(n);
		}
	    }

	    /// Resets the automatic multichannel adaptation batch size for all
	    /// subprocess generators.

	    void unset_auto_channel_adapt()
	    {
		set_auto_channel_adapt(0);
	    }

	    /// Sets the automatic grid adaptation batch size for all
	    /// subprocess generators.

	    void set_auto_grid_adapt(size_type n)
	    {
		for(size_type i=0;i<procs.size();++i)
		{
		    procs[i]->set_auto_grid_adapt(n);
		}
	    }

	    /// Resets the automatic grid adaptation batch size for all
	    /// subprocess generators.

	    void unset_auto_grid_adapt()
	    {
		set_auto_grid_adapt(0);
	    }

	    /// Sets the automatic subprocess channel adaptation batch size.

	    void set_auto_proc_adapt(size_type n)
	    {
		auto_proc_adapt=n;
		if(n!=0)
		{
		    set_auto_update(true);
		}
	    }

	    /// Resets the automatic subprocess channel adaptation batch size.

	    void unset_auto_proc_adapt()
	    {
		auto_proc_adapt=0;
	    }

	    /// Insert a phase space cut.

	    void insert_cut(phase_space_cut* cut)
	    {
		ps_cut=cut;
		for(size_type i=0;i<procs.size();++i)
		{
		    procs[i]->insert_cut(cut);
		}
	    }

	    /// Inserts a scale expression.

	    void insert_scale(scale_expression<value_type>* expr)
	    {
		scale=expr;
		for(size_type i=0;i<procs.size();++i)
		{
		    procs[i]->insert_scale(expr);
		}
	    }

	    /// Generation method with phase space cut argument.
	    
	    bool generate()
	    {
		if(procs.size()==0)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type rho=rn_stream::throw_number();
		sub_proc=procs.begin();
		value_type r(0);
		while(r<rho and sub_proc!=procs.end())
		{
		    r+=(*sub_proc)->alpha;
		    ++sub_proc;
		}
		--sub_proc;
		if(sub_proc==procs.end())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"subprocess iterator overflow detected--no generation performed"<<endlog;
		    return false;
		}
		bool q=(*sub_proc)->generate();
		this->integrand()=(*sub_proc)->integrand();
		this->weight()=(*sub_proc)->weight()/(*sub_proc)->alpha;
		if(update_counter!=0 and auto_proc_adapt!=0 and update_counter%auto_proc_adapt==0)
		{
		    adapt_processes();
		}
		this->up_to_date=false;
		return q;
	    }

	    /// Generates an unweighted event with argument cut object.

	    void generate_unweighted(bool verbose=false)
	    {
		value_type rho=rn_stream::throw_number(0,this->cross_section().value);
		sub_proc=procs.begin();
		value_type r(0);
		while(r<rho and sub_proc!=procs.end())
		{
		    r+=((*sub_proc)->cross_section().value);
		    ++sub_proc;
		}
		--sub_proc;
		if(sub_proc==procs.end())
		{
		    log(log_level::warning)<<"Subprocesses iterator overflow detected...no generation performed"<<endlog;
		    return;
		}
		(*sub_proc)->generate_unweighted(verbose);
		if(update_counter!=0 and auto_proc_adapt!=0 and update_counter%auto_proc_adapt==0)
		{
		    adapt_processes();
		}
		this->up_to_date=false;
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

	    /// Returns whether the current event passes the cuts.

	    bool pass()
	    {
		return (sub_proc==procs.end())?false:((*sub_proc)->pass()); 
	    }

	    /// Resets all subprocess generators.

	    void reset()
	    {
		this->base_type::reset();
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->reset();
		    (*it)->alpha=(value_type)1/procs.size();
		}
		update_counter=0;
	    }

	    /// Resets all subprocess cross sections.

	    void reset_cross_section()
	    {
		this->base_type::reset_cross_section();
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->reset_cross_section();
		}
	    }

	    /// Total cross section refresher method.

	    void refresh_cross_section(bool with_integrand=true)
	    {
		refresh_xsec();
	    }

	    /// Refreshes minimal invariant masses.

	    bool refresh_m_min()
	    {
		bool q=true;
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    q&=(*it)->refresh_m_min();
		}
		return q;
	    }

	    /// Inserts the minimal invariant masses from the configuration
	    /// instance.

	    bool refresh_m_min(const generator_configuration<model_t>& config)
	    {
		bool q=true;
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    q&=(*it)->refresh_m_min(config);
		}
		return q;
	    }

	    /// Refreshes the collider energy.

	    bool refresh_Ecm()
	    {
		bool q=true;
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    q&=(*it)->refresh_Ecm();
		}
		return q;
	    }

	    /// Refreshes internal parameters.

	    bool refresh_params()
	    {
		bool q=true;
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    q&=(*it)->refresh_params();
		}
		return q;
	    }

	    /// Update method.

	    void update()
	    {
		(*sub_proc)->update();
		++update_counter;
	    }

	    /// Grid adaptation method.

	    void adapt_grids()
	    {
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->adapt_grids();
		}
	    }

	    /// Multichannel adaptation method.

	    void adapt_channels()
	    {
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->adapt_channels();
		}
	    }

	    /// Adaptation method.

	    void adapt()
	    {
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->adapt();
		}
		adapt_processes();
	    }

	    /// Adapts the process multichannel parameters to the variances of
	    /// the subprocesses. Removes processes with cross section + error
	    /// smaller than the process threshold times the total cross section.

	    void adapt_processes()
	    {
		value_type p=process_adaptivity();
		value_type norm(0);
		MC_integral<value_type>proc_xsec;
		value_type sigma=cross_section().value;
		value_type q=process_threshold()*sigma/procs.size();
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    proc_xsec=(*it)->cross_section();
		    if((proc_xsec.value+proc_xsec.error)>q)
		    {
			(*it)->alpha=std::pow((*it)->cross_section().error,p);
			norm+=(*it)->alpha;
		    }
		    else
		    {
			std::cout<<"deleting process with xsec "<<proc_xsec<<" smaller than "<<q<<std::endl;
			(*it)->alpha=(value_type)0;
		    }
		}
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->alpha/=norm;
		}
		std::sort(procs.begin(),procs.end(),alpha_more);
		process_iterator it=procs.begin();
		while(it!=procs.end() and (*it)->alpha!=(value_type)0)
		{
		    ++it;
		}
		if(it==procs.end() or it==procs.begin())
		{
		    return;
		}
		std::size_t newsize=it-procs.begin();
		process_iterator it2;
		while(it!=procs.end())
		{
		    algorithm.set_process((*it)->amplitude);
		    algorithm.remove_process();
		    delete (*it);
		    ++it;
		}
		procs.resize(newsize);
		sub_proc=procs.begin();
	    }

	    /// Sets up the weight histogram (default 1000 bins).

	    const weight_histogrammer<value_type>* bin_weights(size_type bins)
	    {
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->bin_weights(bins);
		}
		return NULL;
	    }

	    /// Sets the epsilon-reduction to the maximum weight.

	    value_type reduce_max_weight(const value_type& eps_)
	    {
		for(process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->base_type::reduce_max_weight(eps_);
		}
		return 0;
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /// Returns the total cross section.

	    MC_integral<value_type> cross_section() const
	    {
		if(this->up_to_date)
		{
		    return this->integral;
		}
		refresh_xsec();
		return this->integral;
	    }

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

	    /// Returns the number of subprocesses.

	    size_type processes() const
	    {
		return procs.size();
	    }

	    /// Returns the current process id.

	    size_type process_id() const
	    {
		return ((*sub_proc)->id);
	    }

	    /// Returns the i-th subprocess id.

	    size_type process_id(size_type i) const
	    {
		return (procs[i]->id);
	    }

	    /// Returns the total cross section.

	    MC_integral<value_type> xsec() const
	    {
		return cross_section();
	    }

	    /// Returns the i-th subprocess cross section.

	    MC_integral<value_type> xsec(size_type i) const
	    {
		return procs[i]->cross_section();
	    }

	    /// Returns the current subprocess cross section.

	    MC_integral<value_type> sub_xsec() const
	    {
		return (*sub_proc)->cross_section();
	    }

	    /// Returns the const reference to the i-th incoming momentum (no
	    /// range checking on i).

	    const momentum_type& p_in(size_type i) const
	    {
		return (*sub_proc)->p_in(i);
	    }

	    /// Returns the const reference to the i-th incoming momentum (no
	    /// range checking on i).

	    const momentum_type& p_out(size_type i) const
	    {
		return (*sub_proc)->p_out(i);
	    }

	    /// Returns the const reference to the i-th incoming mass (no
	    /// range checking on i).

	    value_type M_in(size_type i) const
	    {
		return (*sub_proc)->M_in(i);
	    }

	    /// Returns the const reference to the i-th incoming mass (no
	    /// range checking on i).

	    value_type M_out(size_type i) const
	    {
		return (*sub_proc)->M_out(i);
	    }

	    /// Returns the i-th incoming mass-squared (no range checking on i).

	    value_type s_in(size_type i) const
	    {
		return (*sub_proc)->s_in(i);
	    }

	    /// Returns the i-th outgoing mass-squared (no range checking on i).

	    value_type s_out(size_type i) const
	    {
		return (*sub_proc)->s_out(i);
	    }
	    
	    /// Returns the i-th incoming mass (no range checking on i).

	    value_type m_in(size_type i) const
	    {
		return (*sub_proc)->m_in(i);
	    }

	    /// Returns the i-th outgoing mass (no range checking on i).

	    value_type m_out(size_type i) const
	    {
		return (*sub_proc)->m_out(i);
	    }

	    /// Returns the i-th beam id.

	    int beam_id(int i) const
	    {
		return (*sub_proc)->beam_id(i);
	    }

	    /// Returns the cernlib pdf group number for the i-th incoming beam.

	    int pdfg(int i) const
	    {
		return (*sub_proc)->pdfg(i);
	    }

	    /// Returns the cernlib pdf set number for the i-th incoming beam.

	    int pdfs(int i) const
	    {
		return (*sub_proc)->pdfs(i);
	    }

	    /// Returns the i-th incoming particle id.

	    int id_in(size_type i) const
	    {
		return (*sub_proc)->id_in(i);
	    }

	    /// Returns the i-th outgoing particle id.

	    int id_out(size_type i) const
	    {
		return (*sub_proc)->id_out(i);
	    }

	    /// Method determining the colour connection for the event.

	    void fill_colours(std::vector<int>& c,std::vector<int>& cbar) const
	    {
		(*sub_proc)->fill_colours(c,cbar);
	    }

	    /// Returns the total hadronic CM-frame energy.

	    value_type Ecm() const
	    {
		return (*sub_proc)->Ecm();
	    }

	    /// Returns the total hadronic invariant mass-squared.

	    value_type s_tot() const
	    {
		return (*sub_proc)->s_tot();
	    }

	    /// Returns the total partonic CM-frame energy.

	    value_type Ecm_hat() const
	    {
		return (*sub_proc)->Ecm_hat();
	    }

	    /// Returns the total partonic invariant mass-squared.

	    value_type s_hat() const
	    {
		return (*sub_proc)->s_hat();
	    }

	    /// Returns the total event weight.

	    value_type w() const
	    {
		return (*sub_proc)->w();
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
		return (*sub_proc)->ps_weights(with_hels,with_cols);
	    }

	    /// Returns the flux, symmetry and conversion factors. The first
	    /// argument denotes whther helicities are included, the second
	    /// whether colours are.

	    value_type ps_factors(bool with_hels,bool with_cols) const
	    {
		return (*sub_proc)->ps_factors(with_hels,with_cols);
	    }

	    /// Returns the i-th beam energy.

	    value_type beam_energy(int i) const
	    {
		return (*sub_proc)->beam_energy(i);
	    }

	    /// Returns the current factorisation scale.

	    value_type mu_F() const
	    {
		return (*sub_proc)->mu_F();
	    }

	    /// Returns the factorisation scale.

	    value_type F_scale()
	    {
		return (*sub_proc)->F_scale();
	    }

	    /// Returns the renormalisation scale.

	    value_type R_scale()
	    {
		return (*sub_proc)->R_scale();
	    }

	    /// Returns the QCD scale.

	    value_type QCD_scale()
	    {
		return (*sub_proc)->QCD_scale();
	    }

	    /// Returns the const i-th incoming particle phase space (no bound
	    //checking)

	    const phase_space_type* particle_in(size_type i) const
	    {
		return (*sub_proc)->particle_in(i);
	    }

	    /// Returns the const i-th incoming particle phase space (no bound
	    //checking)

	    const phase_space_type* particle_out(size_type i) const
	    {
		return (*sub_proc)->particle_out(i);
	    }

	    /// Returns the i-th const particle phase space, where negative values
	    /// enumerate the incoming particles.

	    const phase_space_type* particle(int i) const
	    {
		return (*sub_proc)->particle(i);
	    }

	    /// Returns the const pointer to the n-th subprocess generator.

	    const process_generator_type* process(size_type n) const
	    {
		return procs[n];
	    }

	    /// Returns the lastly generated const subprocess generator.

	    const process_generator_type* process() const
	    {
		return *sub_proc;
	    }

	    /// Returns the const iterator to the first process generator.
	    
	    const_process_iterator begin_processes() const
	    {
		return procs.begin();
	    }

	    /// Returns the const iterator past the last process generator.
	    
	    const_process_iterator end_processes() const
	    {
		return procs.end();
	    }

	    /// Returns the subprocess weight adaptivity.

	    value_type process_adaptivity() const
	    {
		return subproc_params.first;
	    }

	    /// Returns the subprocess weight threshold.

	    value_type process_threshold() const
	    {
		return subproc_params.second;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Prints the subprocess cross sections.

	    std::ostream& print_cross_sections(std::ostream& os=std::cout) const
	    {
		for(const_process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    std::stringstream ss;
		    (*it)->print_process(ss);
		    os<<std::setw(60)<<ss.str();
		    os<<std::setw(15)<<(*it)->cross_section().value;
		    os<<std::setw(15)<<(*it)->cross_section().error;
		    os<<std::setw(15)<<(*it)->cross_section().error_error;
		    os<<std::setw(15)<<(*it)->alpha;
		    os<<std::setw(10)<<(*it)->calls()<<std::endl;
		}
		os<<"--------------------------------------------------------------------------"<<std::endl;
		os<<"Total: "<<cross_section()<<std::endl;
		return os;
	    }

	    /// Prints the status of the generator.
	    
	    std::ostream& print_status(std::ostream& os=std::cout) const
	    {
		size_type evt_counter=0,pos_evt_counter=0,calls=0,grid_adaptations=0,channel_adaptations=0;
		value_type efficiency;
		for(const_process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    evt_counter+=((*it)->evt_counter);
		    pos_evt_counter+=((*it)->pos_evt_counter);
		    calls+=((*it)->calls());
		    grid_adaptations+=((*it)->grid_adaptations);
		    channel_adaptations+=((*it)->channel_adaptations);
		    efficiency+=((*it)->efficiency()/procs.size());
		}
		os<<"###############################################################################################"<<std::endl;
		os<<"Nr of subprocesses:                                "<<std::scientific<<procs.size()<<std::endl;
		os<<"Nr of events generated:                            "<<std::scientific<<evt_counter<<std::endl;
		os<<"Nr of positive weight events generated:            "<<std::scientific<<pos_evt_counter<<std::endl;
		os<<"Nr of events contributing to cross-section:        "<<std::scientific<<calls<<std::endl;
		os<<"Nr of updates performed:                           "<<std::scientific<<update_counter<<std::endl;
		os<<"Nr of grid adaptations performed:                  "<<std::scientific<<grid_adaptations<<std::endl;
		os<<"Nr of channel adaptations performed:               "<<std::scientific<<channel_adaptations<<std::endl;
		os<<"Mean Monte Carlo efficiency (%):                   "<<std::scientific<<efficiency<<std::endl;
		os<<"Cross section (pb):                                "<<std::scientific<<cross_section()<<std::endl;
		os<<"###############################################################################################"<<std::endl;
	    }

	    /// Prints the generator cuts.

	    std::ostream& print_cuts(std::ostream& os=std::cout) const
	    {
		return (ps_cut==NULL)?os:(ps_cut->print(os));
	    }

	    /// Prints the generator settings.

	    std::ostream& print_settings(std::ostream& os=std::cout) const
	    {
		os<<std::setw(30)<<std::left<<"Nr of subprocesses:"<<procs.size()<<std::endl;
		os<<std::setw(30)<<std::left<<"auto-updates:"<<(auto_update?"yes":"no")<<std::endl;
		if(auto_update)
		{
		    os<<std::setw(30)<<std::left<<"auto-adapt batch size:"<<auto_proc_adapt<<std::endl;
		}
		os<<std::setw(30)<<std::left<<"process adaptivity:"<<subprocess_adaptivity()<<std::endl;
		os<<std::setw(30)<<std::left<<"process xsec threshold:"<<subprocess_threshold()<<std::endl;
		os<<"..........................."<<std::endl;
		os<<"Subprocess settings:	"<<std::endl;
		os<<"..........................."<<std::endl;
		if(procs.size()==0)
		{
		    os<<"no subprocesses"<<std::endl;
		}
		else
		{
		    (*(procs.begin()))->print_settings(os);
		}
		return os;
	    }

	    /// Overridden loading method.

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		is>>update_counter>>auto_update>>auto_proc_adapt>>subproc_params.first>>subproc_params.second;
		return is;
	    }

	    /// Overridden saving method.

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<evtgen>"<<std::endl;
		this->base_type::save(os);
		os<<update_counter<<"\t"<<auto_update<<"\t"<<auto_proc_adapt<<"\t"<<subproc_params.first<<"\t"<<subproc_params.second<<std::endl;
		for(const_process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    (*it)->save(os);
		}
		os<<"</evtgen>"<<std::endl;
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

	protected:

	    /// Camgen algorithm instance.

	    CM_algorithm<model_t,N_in,N_out>& algorithm;

	private:

	    /* Subprocess generators: */

	    process_container procs;

	    /* Subprocess iterator: */

	    process_iterator sub_proc;

	    /* Counts the number of updates since an adapt() call: */

	    size_type update_counter;

	    /* Automatic update flag: */

	    bool auto_update;

	    /* Automatic process channel adaptation flag: */

	    size_type auto_proc_adapt;

	    /* Phase space cut: */

	    phase_space_cut* ps_cut;

	    /* QCD scale expression object: */

	    scale_expression<value_type>* scale;

	    /* Subprocess adaptation parameters: */

	    std::pair<value_type,value_type> subproc_params;

	    /* Private constructor: */
	    
	    event_generator(CM_algorithm<model_t,N_in,N_out>& algo,bool alloc):algorithm(algo),update_counter(0),auto_update(false),auto_proc_adapt(0),ps_cut(NULL),scale(NULL),subproc_params(1,0)
	    {
		procs.reserve(algorithm.n_trees());
		if(alloc)
		{
		    init();
		}
		else
		{
		    sub_proc=procs.end();
		}
	    }

	    /* Constructor helper function: */

	    void init()
	    {
		if(algorithm.reset_process())
		{
		    size_type id=0;
		    do
		    {
			if(!algorithm.get_tree_iterator()->is_empty())
			{
			    procs.push_back(new process_generator_type(algorithm.get_tree_iterator(),id));
			    ++id;
			}
		    }
		    while(algorithm.next_process());
		    if(procs.size()!=0)
		    {
			value_type alpha=(value_type)1/(value_type)procs.size();
			for(size_type i=0;i<procs.size();++i)
			{
			    procs[i]->alpha=alpha;
			}
			sub_proc=procs.begin();
		    }
		    else
		    {
			sub_proc=procs.end();
		    }
		}
		else
		{
		    procs.clear();
		    sub_proc=procs.end();
		}
	    }

	    /* Helper function for initialisation: */

	    void loop_process(size_type sub_proc_evts,bool verbose)
	    {
		if(verbose)
		{
		    std::cout<<"estimating xsec";
		    std::cout.flush();
		    size_type batch=sub_proc_evts/10;
		    for(size_type i=0;i<sub_proc_evts;++i)
		    {
			(*sub_proc)->throw_event();
			(*sub_proc)->update();
			(*sub_proc)->refresh_cross_section();
			this->integrand()=(*sub_proc)->integrand();
			this->weight()=(*sub_proc)->weight()/(*sub_proc)->alpha;
			if(i%batch==0)
			{
			    std::cout<<'.';
			    std::cout.flush();
			}
		    }
		    std::cout<<"done"<<std::endl;
		}
		else
		{
		    for(size_type i=0;i<sub_proc_evts;++i)
		    {
			(*sub_proc)->generate();
			(*sub_proc)->update();
			(*sub_proc)->refresh_cross_section();
			this->integrand()=(*sub_proc)->integrand();
			this->weight()=(*sub_proc)->weight()/(*sub_proc)->alpha;
		    }
		}
		if(verbose)
		{
		    std::cout<<std::setw(15)<<(*sub_proc)->cross_section().value;
		    std::cout<<std::setw(15)<<(*sub_proc)->cross_section().error;
		    std::cout<<std::setw(15)<<(*sub_proc)->cross_section().error_error;
		    std::cout<<std::setw(10)<<(*sub_proc)->calls()<<std::endl;
		}
	    }

	    /* Total cross section refresher method. */

	    void refresh_xsec() const
	    {
		this->integral.value=(value_type)0;
		this->integral.error=(value_type)0;
		this->integral.error_error=(value_type)0;

		for(const_process_iterator it=procs.begin();it!=procs.end();++it)
		{
		    MC_integral<value_type>xsect((*it)->cross_section());
		    (this->integral)+=xsect.value;
		    this->integral.error+=(xsect.error*xsect.error);
		    this->integral.error_error+=(xsect.error_error*xsect.error_error);
		}
		this->integral.error=std::sqrt(this->integral.error);
		this->integral.error_error=std::sqrt(this->integral.error_error);
		this->up_to_date=true;
	    }
    };
}

#endif /*CAMGEN_EVT_GEN_H_*/

