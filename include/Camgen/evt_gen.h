//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_gen.h
  \brief Multi-process Monte Carlo event generator class template.
  */

#ifndef CAMGEN_EVT_GEN_H_
#define CAMGEN_EVT_GEN_H_

#include <Camgen/proc_gen.h>

namespace Camgen
{
    /* Forward declaration of process generator factory base: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class event_generator_factory_base;

    /// Multi-process matrix-element Monte Carlo event generator class.

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class event_generator: public MC_generator<typename model_t::value_type>,
                                                                                                 public MC_integrator_base<typename model_t::value_type>,
                                                                                                 public ps_generator_base<model_t>,
                                                                                                 public phase_space_cut,
                                                                                                 public scale_expression<typename model_t::value_type>
    {
        /* Forward declaration of process generator factory base: */

        friend class event_generator_factory_base<model_t,N_in,N_out,rng_t>;

        typedef MC_integrator_base<typename model_t::value_type> base_type;

        public:

        struct subprocess_type
        {
            typedef process_generator<model_t,N_in,N_out,rng_t> generator_type;

            generator_type* generator;
            typename generator_type::value_type alpha;
        };

        /* Type definitions: */

        typedef model_t model_type;
        typedef typename model_t::value_type value_type;
        typedef typename get_spacetime_type<model_t>::type spacetime_type;
        typedef vector<value_type,spacetime_type::dimension> momentum_type;
        typedef typename momentum_type::size_type size_type;
        typedef rng_t random_number_generator_type;
        typedef random_number_stream<value_type,rng_t> random_number_generator;
        typedef typename CM_algorithm<model_t,N_in,N_out>::tree_iterator CM_tree_iterator;
        typedef process_generator<model_t,N_in,N_out,rng_t> process_generator_type;
        typedef std::vector<subprocess_type> process_container;
        typedef typename process_container::iterator process_iterator;
        typedef typename process_container::const_iterator const_process_iterator;
        typedef MC_integral<value_type> cross_section_type;
        typedef typename CM_algorithm<model_t,N_in,N_out>::phase_space_type phase_space_type;

        /* Public static methods: */
        /*------------------------*/

        /// Process comparison according to ascending cross section.

        static bool xsec_less(const subprocess_type& proc1,const subprocess_type& proc2)
        {
            return proc1.generator->cross_section().value<proc2.generator->cross_section().value;
        }

        /// Process comparison class according to descending cross section.

        static bool xsec_more(const subprocess_type& proc1,const subprocess_type& proc2)
        {
            return proc1.generator->cross_section().value>=proc2.generator->cross_section().value;
        }

        /// Process comparison according to ascending calling frequency.

        static bool alpha_less(const subprocess_type& proc1,const subprocess_type& proc2)
        {
            return proc1.alpha<proc2.alpha;
        }

        /// Process comparison class according to descending calling frequency.

        static bool alpha_more(const subprocess_type& proc1,const subprocess_type& proc2)
        {
            return proc1.alpha>proc2.alpha;
        }

        /// Throws a random value between min and max.

        static value_type throw_number(const value_type& min,const value_type& max)
        {
            return random_number_generator::throw_number(min,max);
        }

        /* Public data: */
        /*--------------*/

        /// Subprocess weight adaptivity.

        value_type process_adaptivity;

        /// Subprocess weight threshold.

        value_type process_threshold;

        /* Public modifying member functions */
        /*-----------------------------------*/

        /// Destructor.

        ~event_generator()
        {
            for(size_type i=0;i<procs.size();++i)
            {
                delete procs[i].generator;
            }
        }

        /// Pre-initialisation method to obtain first subprocess cross
        /// section estimates, determining which subprocesses are relevant.

        void pre_initialise(size_type n_evts,bool verbose=false)
        {
            for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
            {
                sub_proc->generator->pre_initialise(n_evts,verbose);
            }
            this->refresh_cross_section();
            adapt_processes();
            this->base_type::reset();
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->reset();
            }
            update_counter=0;
        }

        /// Initialiser method.

        void initialise(size_type channel_iters,size_type channel_batch,size_type grid_iters,size_type grid_batch,size_type sub_proc_evts,bool verbose=false)
        {
            for(sub_proc=procs.begin();sub_proc!=procs.end();++sub_proc)
            {
                if(verbose)
                {
                    std::stringstream ss;
                    sub_proc->generator->print_process(ss);
                    std::cout<<std::endl<<ss.str()<<std::endl;
                }
                sub_proc->generator->initialise(channel_iters,channel_batch,grid_iters,grid_batch,verbose);
                loop_process(sub_proc_evts,verbose);
            }
            this->refresh_cross_section();
            adapt_processes();
        }

        /// Sets the i-th beam energy.

        bool set_beam_energy(int n,const value_type& E)
        {
            bool q=true;
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                q&=(it->generator->set_beam_energy(n,E));
            }
            return q;
        }

        /// Sets a minimal invariant dimass cut.

        bool set_m_min(int i,int j,const value_type& sqrts)
        {
            bool q=true;
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                q&=it->generator->set_m_min(i,j,sqrts);
            }
            return q;
        }

        /// Sets the auto-update flag for the event generator and all
        /// subprocesses.

        void set_auto_update(bool q)
        {
            auto_update=q;
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->set_auto_update(q);
            }
        }

        /// Sets the automatic multichannel adaptation batch size for all
        /// subprocess generators.

        void set_auto_channel_adapt(size_type n)
        {
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->set_auto_channel_adapt(n);
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
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->set_auto_grid_adapt(n);
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
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->insert_cut(cut);
            }
        }

        /// Inserts a scale expression.

        void insert_scale(scale_expression<value_type>* expr)
        {
            scale=expr;
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->insert_scale(expr);
            }
        }

        /// Generation method with phase space cut argument.

        bool generate()
        {
            if(procs.size()==0)
            {
                this->weight()=(value_type)0;
                return true;
            }
            value_type rho=throw_number((value_type)0,(value_type)1);
            sub_proc=procs.begin();
            value_type r(0);
            while(r<=rho and sub_proc!=procs.end())
            {
                r+=sub_proc->alpha;
                ++sub_proc;
            }
            --sub_proc;
            if(sub_proc==procs.end())
            {
                log(log_level::warning)<<CAMGEN_STREAMLOC<<"subprocess iterator overflow detected--no generation performed"<<endlog;
                return false;
            }
            bool q=sub_proc->generator->generate();
//            sub_proc->generator->refresh_cross_section();
            this->integrand()=sub_proc->generator->integrand();
            this->weight()=sub_proc->generator->weight()/sub_proc->alpha;
            if(update_counter!=0 and auto_proc_adapt!=0 and update_counter%auto_proc_adapt==0)
            {
                adapt_processes();
            }
            return q;
        }

        /// Evaluates the weight of the current sub-process.

        bool evaluate_weight()
        {
            if(sub_proc!=procs.end())
            {
                return sub_proc->generator->evaluate_weight();
            }
            return false;
        }

        /// Generates an unweighted event with argument cut object.

        void generate_unweighted(bool verbose=false)
        {
            value_type rho=throw_number(0,this->cross_section().value);
            sub_proc=procs.begin();
            value_type r(0);
            while(r<rho and sub_proc!=procs.end())
            {
                r+=(sub_proc->generator->cross_section().value);
                ++sub_proc;
            }
            --sub_proc;
            if(sub_proc==procs.end())
            {
                log(log_level::warning)<<"Subprocesses iterator overflow detected...no generation performed"<<endlog;
                return;
            }
            sub_proc->generator->generate_unweighted(verbose);
            if(update_counter!=0 and auto_proc_adapt!=0 and update_counter%auto_proc_adapt==0)
            {
                adapt_processes();
            }
        }

        void refresh_cross_section(bool with_integrand=true)
        {
            if(sub_proc!=procs.end())
            {
                sub_proc->generator->refresh_cross_section(with_integrand);
            }
            n_calls++;
            up_to_date=false;
        }

        /// Resets all subprocess generators.

        void reset()
        {
            this->base_type::reset();
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->reset();
                it->alpha=(value_type)1/procs.size();
            }
            update_counter=0;
        }

        /// Resets all subprocess cross sections.

        void reset_cross_section()
        {
            this->base_type::reset_cross_section();
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->reset_cross_section();
            }
        }

        /// Generates new event according to the given strategy argument.

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
            return (sub_proc==procs.end())?false:(sub_proc->generator->pass()); 
        }

        /// Refreshes minimal invariant masses.

        bool refresh_m_min()
        {
            bool q=true;
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                q&=it->generator->refresh_m_min();
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
                q&=it->generator->refresh_m_min(config);
            }
            return q;
        }

        /// Refreshes the collider energy.

        bool refresh_Ecm()
        {
            bool q=true;
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                q&=it->generator->refresh_Ecm();
            }
            return q;
        }

        /// Refreshes internal parameters.

        bool refresh_params()
        {
            bool q=true;
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                q&=it->generator->refresh_params();
            }
            return q;
        }

        /// Update method.

        void update()
        {
            if(sub_proc!=procs.end())
            {
                sub_proc->generator->update();
            }
            ++update_counter;
        }

        /// Grid adaptation method.

        void adapt_grids()
        {
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->adapt_grids();
            }
        }

        /// Multichannel adaptation method.

        void adapt_channels()
        {
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->adapt_channels();
            }
        }

        /// Adaptation method.

        void adapt()
        {
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->adapt();
            }
            adapt_processes();
        }

        /// Adapts the process multichannel parameters to the variances of
        /// the subprocesses. Removes processes with cross section + error
        /// smaller than the process threshold times the total cross section.

        void adapt_processes()
        {
            value_type norm(0);
            MC_integral<value_type>proc_xsec;
            value_type sigma=this->cross_section().value;
            value_type q=process_threshold*sigma/procs.size();
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                proc_xsec=it->generator->cross_section();
                if((proc_xsec.value+proc_xsec.error)>q)
                {
                    it->alpha=std::pow(it->generator->cross_section().error,process_adaptivity);
                    norm+=it->alpha;
                }
                else
                {
                    it->alpha=(value_type)0;
                }
            }
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->alpha/=norm;
            }
            std::sort(procs.begin(),procs.end(),alpha_more);
            process_iterator it=procs.begin();
            while(it!=procs.end() and it->alpha!=(value_type)0)
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
                algorithm.set_process(it->generator->amplitude);
                algorithm.remove_process();
                delete it->generator;
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
                it->generator->bin_weights(bins);
            }
            return NULL;
        }

        /// Sets the epsilon-reduction to the maximum weight.

        value_type reduce_max_weight(const value_type& eps_)
        {
            for(process_iterator it=procs.begin();it!=procs.end();++it)
            {
                it->generator->base_type::reduce_max_weight(eps_);
            }
            return 0;
        }

        /* Public readout methods */
        /*------------------------*/

        /// Cross section implementation.

        MC_integral<value_type> cross_section() const
        {
            if(!up_to_date)
            {
                integral_sum=compute_cross_section();
                up_to_date=true;
            }
            return integral_sum;
        }

        /// Returns the sum of subprocess cross sections.

        cross_section_type compute_cross_section() const
        {
            MC_integral<value_type> summed_xsec;

            for(const_process_iterator it=procs.begin();it!=procs.end();++it)
            {
                MC_integral<value_type>xsect(it->generator->cross_section());
                summed_xsec+=xsect.value;
                summed_xsec.error+=(xsect.error*xsect.error);
                summed_xsec.error_error+=(xsect.error_error*xsect.error_error);
            }
            summed_xsec.error=std::sqrt(summed_xsec.error);
            summed_xsec.error_error=std::sqrt(summed_xsec.error_error);
            return summed_xsec;
        }

        /// Returns the number of cross section updates.
        
        size_type calls() const
        {
            return n_calls;
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
            return sub_proc->generator->id;
        }

        /// Returns the i-th subprocess id.

        size_type process_id(size_type i) const
        {
            return procs[i].generator->id;
        }

        /// Returns the total cross section.

        cross_section_type xsec() const
        {
            return this->cross_section();
        }

        /// Returns the i-th subprocess cross section.

        MC_integral<value_type> xsec(size_type i) const
        {
            return procs[i].generator->cross_section();
        }

        /// Returns the current subprocess cross section.

        MC_integral<value_type> sub_xsec() const
        {
            return sub_proc->generator->cross_section();
        }

        /// Returns the const reference to the i-th incoming momentum (no
        /// range checking on i).

        const momentum_type& p_in(size_type i) const
        {
            return sub_proc->generator->p_in(i);
        }

        /// Returns the const reference to the i-th incoming momentum (no
        /// range checking on i).

        const momentum_type& p_out(size_type i) const
        {
            return sub_proc->generator->p_out(i);
        }

        /// Returns the const reference to the i-th incoming mass (no
        /// range checking on i).

        value_type M_in(size_type i) const
        {
            return sub_proc->generator->M_in(i);
        }

        /// Returns the const reference to the i-th incoming mass (no
        /// range checking on i).

        value_type M_out(size_type i) const
        {
            return sub_proc->generator->M_out(i);
        }

        /// Returns the i-th incoming mass-squared (no range checking on i).

        value_type s_in(size_type i) const
        {
            return sub_proc->generator->s_in(i);
        }

        /// Returns the i-th outgoing mass-squared (no range checking on i).

        value_type s_out(size_type i) const
        {
            return sub_proc->generator->s_out(i);
        }

        /// Returns the i-th incoming mass (no range checking on i).

        value_type m_in(size_type i) const
        {
            return sub_proc->generator->m_in(i);
        }

        /// Returns the i-th outgoing mass (no range checking on i).

        value_type m_out(size_type i) const
        {
            return sub_proc->generator->m_out(i);
        }

        /// Returns the i-th beam id.

        int beam_id(int i) const
        {
            return sub_proc->generator->beam_id(i);
        }

        /// Returns the cernlib pdf group number for the i-th incoming beam.

        int pdfg(int i) const
        {
            return sub_proc->generator->pdfg(i);
        }

        /// Returns the cernlib pdf set number for the i-th incoming beam.

        int pdfs(int i) const
        {
            return sub_proc->generator->pdfs(i);
        }

        /// Returns the i-th incoming particle id.

        int id_in(size_type i) const
        {
            return sub_proc->generator->id_in(i);
        }

        /// Returns the i-th outgoing particle id.

        int id_out(size_type i) const
        {
            return sub_proc->generator->id_out(i);
        }

        /// Method determining the colour connection for the event.

        void fill_colours(std::vector<int>& c,std::vector<int>& cbar) const
        {
            sub_proc->generator->fill_colours(c,cbar);
        }

        /// Returns the total hadronic CM-frame energy.

        value_type Ecm() const
        {
            return sub_proc->generator->Ecm();
        }

        /// Returns the total hadronic invariant mass-squared.

        value_type s_tot() const
        {
            return sub_proc->generator->s_tot();
        }

        /// Returns the total partonic CM-frame energy.

        value_type Ecm_hat() const
        {
            return sub_proc->generator->Ecm_hat();
        }

        /// Returns the total partonic invariant mass-squared.

        value_type s_in() const
        {
            return sub_proc->generator->s_in();
        }

        /// Returns the total event weight.

        value_type w() const
        {
            return sub_proc->generator->w();
        }

        /// Returns the maximal event weight.

        value_type max_w() const
        {
            std::vector<value_type> maxwts(procs.size());
            typename std::vector<value_type>::iterator it1=maxwts.begin();
            typename std::vector<value_type>::iterator it2=std::transform(begin_processes(),end_processes(),it1,get_max_weight);
            return *(std::max_element(it1,it2));
        }

        static value_type get_max_weight(const subprocess_type& p)
        {
            return p.generator->max_weight();
        }

        /// Returns the event weight without flux and symmetry factors etc.
        /// The first argument denotes whether helicities are included, the
        /// second whether colours are.

        value_type ps_weights(bool with_hels,bool with_cols) const
        {
            return sub_proc->generator->ps_weights(with_hels,with_cols);
        }

        /// Returns the flux, symmetry and conversion factors. The first
        /// argument denotes whther helicities are included, the second
        /// whether colours are.

        value_type ps_factors(bool with_hels,bool with_cols) const
        {
            return sub_proc->generator->ps_factors(with_hels,with_cols);
        }

        /// Returns the i-th beam energy.

        value_type beam_energy(int i) const
        {
            return sub_proc->generator->beam_energy(i);
        }

        /// Returns the current factorisation scale.

        value_type mu_F() const
        {
            return sub_proc->generator->mu_F();
        }

        /// Returns the factorisation scale.

        value_type F_scale()
        {
            return sub_proc->generator->F_scale();
        }

        /// Returns the renormalisation scale.

        value_type R_scale()
        {
            return sub_proc->generator->R_scale();
        }

        /// Returns the QCD scale.

        value_type QCD_scale()
        {
            return sub_proc->generator->QCD_scale();
        }

        /// Returns the const i-th incoming particle phase space (no bound
        //checking)

        const phase_space_type* particle_in(size_type i) const
        {
            return sub_proc->generator->particle_in(i);
        }

        /// Returns the const i-th incoming particle phase space (no bound
        //checking)

        const phase_space_type* particle_out(size_type i) const
        {
            return sub_proc->generator->particle_out(i);
        }

        /// Returns the i-th const particle phase space, where negative values
        /// enumerate the incoming particles.

        const phase_space_type* particle(int i) const
        {
            return sub_proc->generator->particle(i);
        }

        /// Returns the const pointer to the n-th subprocess generator.

        const process_generator_type* process(size_type n) const
        {
            return procs[n].generator;
        }

        /// Returns the lastly generated const subprocess generator.

        const process_generator_type* process() const
        {
            return sub_proc->generator;
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

        /* Serialization: */
        /*----------------*/

        /// Prints the subprocess cross sections.

        std::ostream& print_cross_sections(std::ostream& os=std::cout) const
        {
            for(const_process_iterator it=procs.begin();it!=procs.end();++it)
            {
                std::stringstream ss;
                it->generator->print_process(ss);
                os<<std::setw(60)<<ss.str();
                os<<std::setw(15)<<it->generator->cross_section().value;
                os<<std::setw(15)<<it->generator->cross_section().error;
                os<<std::setw(15)<<it->generator->cross_section().error_error;
                os<<std::setw(15)<<it->alpha;
                os<<std::setw(10)<<it->generator->calls()<<std::endl;
            }
            os<<"--------------------------------------------------------------------------"<<std::endl;
            os<<"Total: "<<this->cross_section()<<std::endl;
            return os;
        }

        /// Prints the status of the generator.

        std::ostream& print_status(std::ostream& os=std::cout) const
        {
            size_type evt_counter=0,pos_evt_counter=0,calls=0,grid_adaptations=0,channel_adaptations=0;
            value_type efficiency=0;
            for(const_process_iterator it=procs.begin();it!=procs.end();++it)
            {
                evt_counter+=(it->generator->evt_counter);
                pos_evt_counter+=(it->generator->pos_evt_counter);
                calls+=(it->generator->calls());
                grid_adaptations+=(it->generator->grid_adaptations);
                channel_adaptations+=(it->generator->channel_adaptations);
                efficiency+=(it->generator->efficiency()/procs.size());
            }
            os<<"###############################################################################################"<<std::endl;
            os<<"Nr of subprocesses:                                "<<std::scientific<<procs.size()<<std::endl;
            os<<"Nr of events generated:                            "<<std::scientific<<evt_counter<<std::endl;
            os<<"Nr of positive weight events generated:            "<<std::scientific<<pos_evt_counter<<std::endl;
            os<<"Nr of events contributing to cross-section:        "<<std::scientific<<calls<<std::endl;
            os<<"Nr of updates performed:                           "<<std::scientific<<update_counter<<std::endl;
            os<<"Nr of grid adaptations performed:                  "<<std::scientific<<grid_adaptations<<std::endl;
            os<<"Nr of channel adaptations performed:               "<<std::scientific<<channel_adaptations<<std::endl;
            os<<"Monte Carlo efficiency (%):                        "<<std::scientific<<efficiency<<std::endl;
            os<<"Cross section (pb):                                "<<std::scientific<<this->cross_section()<<std::endl;
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
            os<<std::setw(30)<<std::left<<"Nr of subprocesses:"<<procs.size()<<std::endl;
            os<<std::setw(30)<<std::left<<"auto-updates:"<<(auto_update?"yes":"no")<<std::endl;
            if(auto_update)
            {
                os<<std::setw(30)<<std::left<<"auto-adapt batch size:"<<auto_proc_adapt<<std::endl;
            }
            os<<std::setw(30)<<std::left<<"process adaptivity:"<<process_adaptivity<<std::endl;
            os<<std::setw(30)<<std::left<<"process xsec threshold:"<<process_threshold<<std::endl;
            os<<"..........................."<<std::endl;
            os<<"Subprocess settings:	"<<std::endl;
            os<<"..........................."<<std::endl;
            if(procs.size()==0)
            {
                os<<"no subprocesses"<<std::endl;
            }
            else
            {
                procs.begin()->generator->print_settings(os);
            }
            return os;
        }

        protected:

        /// Camgen algorithm instance.

        CM_algorithm<model_t,N_in,N_out>& algorithm;

        /// Constructor from static configuration data.

        event_generator(CM_algorithm<model_t,N_in,N_out>& algo):process_adaptivity(1),process_threshold(0),algorithm(algo),up_to_date(false),n_calls(0),update_counter(0),auto_update(false),auto_proc_adapt(0),ps_cut(NULL),scale(NULL){}

        private:

        /* Subprocess generators: */

        process_container procs;

        /* Subprocess iterator: */

        process_iterator sub_proc;

        /* Total cross section: */

        mutable MC_integral<value_type> integral_sum;

        /* Lazy cross section update flag: */

        mutable bool up_to_date;

        /* Counts the number of cross section updates: */

        size_type n_calls;

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
                    sub_proc->generator->throw_event();
                    sub_proc->generator->update();
                    sub_proc->generator->refresh_cross_section();
                    this->integrand()=sub_proc->generator->integrand();
                    this->weight()=sub_proc->generator->weight()/sub_proc->generator->alpha;
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
                    sub_proc->generator->generate();
                    sub_proc->generator->update();
                    sub_proc->generator->refresh_cross_section();
                    this->integrand()=sub_proc->generator->integrand();
                    this->weight()=sub_proc->generator->weight()/sub_proc->alpha;
                }
            }
            if(verbose)
            {
                std::cout<<std::setw(15)<<sub_proc->generator->cross_section().value;
                std::cout<<std::setw(15)<<sub_proc->generator->cross_section().error;
                std::cout<<std::setw(15)<<sub_proc->generator->cross_section().error_error;
                std::cout<<std::setw(10)<<sub_proc->generator->calls()<<std::endl;
            }
        }
    };
}

#endif /*CAMGEN_EVT_GEN_BASE_H_*/

