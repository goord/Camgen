//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PROCESS_TREE_H_
#define CAMGEN_PROCESS_TREE_H_

#include <list>
#include <set>
#include <Camgen/current_tree.h>
#include <Camgen/interaction.h>
#include <Camgen/bspart.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the process tree class template in Camgen. This class is at  *
 * the core of Camgen, and wraps all the interactions between the currents in *
 * the static current tree for a given subprocess.                             *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out>class process_tree
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	    
	    /* Particle, vertex and phase space type definitions: */
	    
	    typedef typename get_basic_types<model_t>::particle_type particle_type;
	    typedef typename get_basic_types<model_t>::vertex_type vertex_type;
	    typedef typename get_basic_types<model_t>::phase_space_type phase_space_type;
	    
	    /* Fusion map and fusion iterator type definitions: */
	    
	    typedef typename model_wrapper<model_t>::fusion_map fusion_map;
	    typedef typename model_wrapper<model_t>::fusion_iterator fusion_iterator;

	    /* Number of particles that the tree is built upon: */

	    static const std::size_t N_bits=N_in+N_out-1;
	    
	    /* Boolean denoting if we are in colour decomposition mode: */
	    
	    static const bool decomposes=get_colour_treatment<model_t>::decomposes;

	    /* Current, interaction and current tree type definitions: */

	    typedef current<model_t,N_bits,decomposes> current_type;
	    typedef interaction<model_t,N_bits,decomposes> interaction_type;
	    typedef current_tree<model_t,N_in,N_out> current_tree_type;
	    
	    /* Iterator definitions: */
	    
	    typedef typename std::vector<current_type>::iterator current_iterator;
	    typedef typename std::vector<current_type>::const_iterator const_current_iterator;
	    typedef typename std::list<interaction_type>::iterator interaction_iterator;
	    typedef typename std::list<interaction_type>::const_iterator const_interaction_iterator;
	    typedef typename std::list<interaction_type>::reverse_iterator reverse_interaction_iterator;

	    /* Declaring the CM_algorithm user interface class friend: */

	    friend class CM_algorithm<model_t,N_in,N_out>;

	    /* Static initialiser: */

	    static void initialise(std::size_t N_f)
	    {
		current_tree_type::initialise(N_f);
		N_final=current_tree_type::final_current();

		if(!initialised)
		{
		    // TODO: refactor to deal with dynamically changing model
		    // structure...

		    max_rank=model_wrapper<model_t>::maximal_vertex_rank();
		    flavours=model_wrapper<model_t>::flavours();
		    bit_string_partition<N_bits>::initialise(max_rank-1);
		    bit_string_partition<N_bits>::sort();
		    initialised=true;
		}
	    }

	    /* Constructor from a vector of incoming particles and a vector of
	     * outgoing particles: */

	    process_tree(const vector<const particle_type*,N_in>& in,const vector<const particle_type*,N_out>& out):empty(true),counter(0)
	    {

		/* Copying initial current addresses to the process tree: */

		init_currents.reserve(N_bits);
		current_iterator it=current_tree_type::begin();
		for(size_type i=0;i<N_in;++i)
		{
		    if(in[i]==NULL)
		    {
			log(log_level::error)<<CAMGEN_STREAMLOC<<"incoming particle "<<i<<" was not instantiated"<<endlog;
		    }
		    if(i != current_tree_type::final_current())
		    {
			init_currents.push_back(it+in[i]->get_flavour());
			it+=flavours;
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    if(out[i]==NULL)
		    {
			log(log_level::error)<<CAMGEN_STREAMLOC<<"outgoing particle "<<i<<" was not instantiated"<<endlog;
		    }
		    if(i != (N_final-N_in))
		    {
			init_currents.push_back(it+out[i]->get_flavour());
			it+=flavours;
		    }
		}

		/* Copying final current addresses to the process tree: */

		it=current_tree_type::end()-flavours;
		if(N_final<N_in)
		{
		    final_current=it+in[N_final]->get_flavour();
		}
		else
		{
		    final_current=it+out[N_final-N_in]->get_flavour();
		}
		evaluate_symmetry_factor();
	    }

	    /* Tree building initialisation function, marking the initial and
	     * final currents: */

	    void init_build()
	    {
		for(size_type i=0;i<init_currents.size();++i)
		{
		    init_currents[i]->mark();
		}
		B.set(N_bits-1);
		level=1;
		final_current->mark();
	    }

	    /* Tree building function: */

	    void build()
	    {
		/* Initialise: */

		init_build();
		
		/* Recursive build: */
		
		for(size_type i=0;i<(1<<N_bits)-N_bits-2;++i)
		{
		    build_next();
		    if(B.count()/(max_rank-1)>level)
		    {
			empty=true;
			break;
		    }
		}

		/* Final build: */

		if(N_bits/(max_rank-1)<=level)
		{
		    build_final();
		}

		/* Sorting the interactions: */

		if(!empty)
		{
		    interactions.sort();
		}
	    }

	    /* Cleaning algorithm, erasing redundant interactions: */

	    void clean()
	    {
		if(!empty)
		{
		    /* Unmark all currents: */

		    current_tree_type::unmark();

		    /* Starting from the last interaction, iteratively mark all
		     * participating currents and their interactions: */

		    reverse_interaction_iterator it=interactions.rbegin();
		    it->mark();
		    for(;it != interactions.rend();++it)
		    {
			it->mark_if();
		    }

		    /* If the initial currents don't end up being marked, the
		     * process is empty: */

		    for(size_type i=0;i<N_bits;++i)
		    {
			if(!(init_currents[i]->is_marked()))
			{
			    empty=true;
			    interactions.clear();
			    return;
			}
		    }

		    /* Remove the unmarked interactions: */

		    redundant<interaction_type>pred;
		    interactions.remove_if(pred);
		}
		else
		{
		    interactions.clear();
		}

		/* Unmark all currents again: */

		current_tree_type::unmark();
	    }

	    /* Current initialisation phase: */

	    void initialise_currents()
	    {
		/* Initialise external currents: */

		for(size_type i=0;i<N_bits;++i)
		{
		    init_currents[i]->initialise();
		}

		/* Initialise produced currents: */

		for(interaction_iterator it=interactions.begin();it != interactions.end();++it)
		{
		    it->get_produced_current()->initialise();
		}

		/* Initialise final current: */

		final_current->initialise();

		/* Set the contraction argument of the final current: */

		if(!empty)
		{
		    final_current->set_argument(&(interactions.back().get_produced_current()->amplitude));
		}
	    }

	    /* Momentum policy determination function: */

	    void assign_momenta()
	    {
		if(!empty)
		{
		    /* The first interaction should compute its produced
		     * momentum: */

		    interaction_iterator it=interactions.begin();
		    it->set_computes_momentum();
		    
		    /* Copy the momentum channel bit string, particle type and
		     * produced momentum to resp. b, phi and p: */

		    bit_string<N_bits>b=it->get_produced_bit_string();
		    const particle<model_t>* phi=it->get_produced_particle();
		    const momentum_type* p=&(it->get_produced_current()->momentum);
		    ++it;
		    
		    /* Iteratively contstruct the momentum policy: */
		    
		    for(;it != interactions.end();++it)
		    {
			/* If the momentum channel changes, let the new
			 * interaction compute the produced momentum, and the
			 * previous interaction propagate (unless it produces
			 * the final leg): */

			if(it->get_produced_bit_string()!=b)
			{
			    it->set_computes_momentum();
			    --it;
			    if(it->get_produced_bit_string().count()<N_bits)
			    {
				it->set_propagates();
			    }
			    ++it;

			    /* Copy the momentum channel, particle type and
			     * momentum address to b, phi and p: */

			    b=it->get_produced_bit_string();
			    phi=it->get_produced_particle();
			    p=&(it->get_produced_current()->momentum);
			}

			/* Within a momentum channel... */
			else
			{
			    /* If the produced particle type changes, the new
			     * interaction should assign the momentum from *p,
			     * and the previous interaction should let its
			     * produced current propagate (unless this is the
			     * final leg): */

			    if(it->get_produced_particle()!=phi)
			    {
				it->set_assigns_momentum(p);
				--it;
				if(it->get_produced_bit_string().count()<N_bits)
				{
				    it->set_propagates();
				}
				++it;

				phi=it->get_produced_particle();
			    }
			}
		    }
		}
	    }

	    /* Computes all the coupling flags (this is done every evaluation
	     * round...)*/

	    void compute_coupling_flags()
	    {
		for(interaction_iterator it=interactions.begin();it!=interactions.end();++it)
		{
		    it->compute_coupling_flag();
		}
	    }


	    /* Function computing the Fermi signs of the interactions: */

	    void set_Fermi_signs()
	    {
		for(interaction_iterator it=interactions.begin();it!=interactions.end();++it)
		{
		    /* Decompose the produced momentum channel and if the
		     * participating external momentum belongs to a fermion,
		     * store the leg number in the vector f: */
		    
		    std::vector<size_type>f;
		    for(typename std::vector<current_iterator>::const_iterator it2=it->begin();it2!=it->end();++it2)
		    {
			if(*it2 != it->get_produced_current())
			{
			    bit_string<N_bits>b=(*it2)->get_bit_string();
			    for(size_type i=0;i<N_bits;++i)
			    {
				if(b[i] and init_currents[i]->is_fermion())
				{
				    f.push_back(i);
				}
			    }
			}
		    }

		    /* The Fermi sign corresponds to the relative ordering of f,
		     * i.e. the sign of the permutation sorting f: */

		    if(f.size()>1)
		    {
			int sign=1;
			for(size_type i=0;i<f.size();++i)
			{
			    for(size_type j=0;j<i;++j)
			    {
				if(f[j]>f[i])
				{
				    sign=-sign;
				}
			    }
			}
			if(sign==-1)
			{
			    it->set_Fermi_sign();
			}
		    }
		}
	    }

	    /* Feynman diagram counting facility: */

	    long long unsigned count_diagrams()
	    {
		if(empty)
		{
		    return 0;
		}
		
		/* Reset all multiplicities: */

		for(current_iterator it=current_tree_type::begin();it != current_tree_type::end();++it)
		{
		    it->multiplicity=0;
		}

		/* Set all external current multiplicities to 1: */

		for(size_type i=0;i<init_currents.size();++i)
		{
		    init_currents[i]->multiplicity=1;
		}

		/* Let all interactions count multiplicities: */

		for(interaction_iterator it=interactions.begin();it != interactions.end();++it)
		{
		    it->count_diagrams();
		}
		return interactions.back().get_produced_current()->multiplicity;
	    }

	    /* Counts the number of participating interactions: */

	    long long unsigned count_current_combinations()
	    {
		long long unsigned n=0;
		for(interaction_iterator it=interactions.begin();it!=interactions.end();++it)
		{
		    if(it->is_coupled())
		    {
			++n;
		    }
		}
		return n;
	    }

	    /* Process printing facility: */

	    std::ostream& print_process(std::ostream& os) const
	    {
		std::vector<std::string>ext_parts;
		if(N_in+N_out>0)
		{
		    ext_parts.resize(N_bits);
		    for(size_type i=0;i<N_bits;++i)
		    {
			ext_parts[i]=init_currents[i]->get_particle_type()->get_name();
		    }
		    ext_parts.insert(ext_parts.begin()+N_final,final_current->get_particle_type()->get_name());
		}
		if(N_out>0)
		{
		    if(N_in>0)
		    {
			for(size_type i=0;i<N_in-1;++i)
			{
			    os<<ext_parts[i]<<",";
			}
			os<<ext_parts[N_in-1]<<" > ";
		    }
		    else
		    {
			os<<"none > ";
		    }
		    for(size_type i=N_in;i<N_bits;++i)
		    {
			os<<ext_parts[i]<<",";
		    }
		    os<<ext_parts[N_bits];
		}
		else
		{
		    if(N_in>0)
		    {
			for(size_type i=0;i<N_in-1;++i)
			{
			    os<<ext_parts[i]<<",";
			}
			os<<ext_parts[N_in-1]<<" > none";
		    }
		    else
		    {
			os<<"empty process";
			return os;
		    }
		}
		return os;
	    }

	    /* Process tree printing function (with internal subamplitudes): */

	    std::ostream& print(std::ostream& os) const
	    {
		print_process(os);
		os<<std::endl;
		
		for(size_type i=0;i<N_bits;++i)
		{
		    init_currents[i]->print(os);
		    os<<std::endl;
		    init_currents[i]->print_amp(os);
		    os<<std::endl<<std::endl;
		}
		if(!empty)
		{
		    const_interaction_iterator it=interactions.begin();
		    for(;it!=--interactions.end();++it)
		    {
			if(it->get_vertex()->is_coupled())
			{
			    it->print(os);
			    os<<std::endl;
			    if(it->propagates())
			    {
				it->get_produced_current()->print_amp(os);
				os<<std::endl<<std::endl;
			    }
			}
		    }
		    it->print(os);
		    os<<std::endl;
		    it->get_produced_current()->print_amp(os);
		    os<<std::endl<<std::endl;
		}
		final_current->print(os);
		os<<std::endl;
		final_current->print_amp(os);
		os<<std::endl;
		return os;
	    }

	    template<class rep_t>std::ostream& print_cf_transform(std::ostream& os) const
	    {
		print_process(os);
		os<<std::endl;
		for(size_type i=0;i<N_bits;++i)
		{
		    init_currents[i]->print(os);
		    os<<std::endl;
		    init_currents[i]->template print_cf_amp<rep_t>(os);
		    os<<std::endl<<std::endl;
		}
		if(!empty)
		{
		    const_interaction_iterator it=interactions.begin();
		    for(;it!=--interactions.end();++it)
		    {
			it->print(os);
			os<<std::endl;
			if(it->propagates())
			{
			    it->get_produced_current()->template print_cf_amp<rep_t>(os);
			    os<<std::endl<<std::endl;
			}
		    }
		    it->print(os);
		    os<<std::endl;
		    it->get_produced_current()->template print_cf_amp<rep_t>(os);
		    os<<std::endl<<std::endl;
		}
		final_current->print(os);
		os<<std::endl;
		final_current->template print_cf_amp<rep_t>(os);
		os<<std::endl;
		return os;
	    }

	    /* Process tree printing function (without subamplitudes): */

	    std::ostream& print_tree(std::ostream& os) const
	    {
		print_process(os);
		os<<std::endl;
		for(size_type i=0;i<N_bits;++i)
		{
		    init_currents[i]->print(os);
		    os<<std::endl;
		}
		for(const_interaction_iterator it=interactions.begin();it!=interactions.end();++it)
		{
		    it->print(os);
		    os<<std::endl;
		}
		final_current->print(os);
		os<<std::endl;
		return os;
	    }

	    /* Function printing all vertices participating in the amplitude: */

	    std::ostream& print_vertices(std::ostream& os) const
	    {
		if(!empty)
		{
		    std::set<const vertex_type*>vertices;
		    for(const_interaction_iterator it1=interactions.begin();it1!=interactions.end();++it1)
		    {
			if(it1->get_vertex()->is_coupled())
			{
			    vertices.insert(it1->get_vertex());
			}
		    }
		    for(typename std::set<const vertex_type*>::const_iterator it2=vertices.begin();it2!=vertices.end();++it2)
		    {
			os<<(*it2)->get_name()<<std::endl;
		    }
		}
		return os;
	    }

	    /* Request whether the process is empty: */

	    bool is_empty() const
	    {
		return empty;
	    }

	    std::string helicity_print(int i)
	    {
		if(i==-1)
		{
		    return "-";
		}
		if(i==0)
		{
		    return "0";
		}
		if(i==1)
		{
		    return "+";
		}
		return "nan";
	    }

	    /* Spin summed amplitude evaluation: */

	    r_value_type evaluate(std::bitset<N_in+N_out> summed_spins,std::bitset<N_in+N_out> summed_cols)
	    {
		if(empty){return 0;}
		r_value_type colsum_factor(1);
		for(size_type i=0;i<N_final;++i)
		{
		    init_currents[i]->init_dof_sum(summed_spins[i],summed_cols[i]);
		    if(summed_cols[i])
		    {
			colsum_factor*=(init_currents[i]->colour_sum_factor());
		    }
		}
		final_current->init_dof_sum(summed_spins[N_final],summed_cols[N_final]);
		if(summed_cols[N_final])
		{
		    colsum_factor*=(final_current->colour_sum_factor());
		}
		for(size_type i=N_final;i<init_currents.size();++i)
		{
		    init_currents[i]->init_dof_sum(summed_spins[i+1],summed_cols[i+1]);
		    if(summed_cols[i+1])
		    {
			colsum_factor*=(init_currents[i]->colour_sum_factor());
		    }
		}
		r_value_type summed_amplitude(0);
		size_type counter=0;
		do
		{
		    for(interaction_iterator it=interactions.begin();it!=interactions.end();++it)
		    {
			it->evaluate();
		    }
		    r_value_type M=std::norm(final_current->contract_wave_function());
		    summed_amplitude+=M;
		    ++counter;
		    for(interaction_iterator it=interactions.begin();it!=interactions.end();++it)
		    {
			it->reset();
		    }
		    (--interactions.end())->get_produced_current()->reset();
		}
		while(next_dof_configuration(summed_spins,summed_cols));
		return colsum_factor*summed_amplitude;
	    }

	    /* Tree evaluation function: */

	    value_type evaluate()
	    {
		if(empty){return value_type(0,0);}

		/* Evaluate initial external wave functions: */

		for(size_type i=0;i<init_currents.size();++i)
		{
		    init_currents[i]->evaluate();
		}

		/* Evaluate all recursive relations: */

		for(interaction_iterator it=interactions.begin();it != interactions.end();++it)
		{
		    it->evaluate();
		}

		/* Evaluate the final-particle wave function: */

		final_current->evaluate();

		/* Return the contraction of the final particle wave function
		 * with the final internal subamplitude: */

		return final_current->contract_wave_function();
	    }

	    /* Function evaluating the next interaction/wave function w.r.t. the
	     * counter: */

	    size_type evaluate_next()
	    {
		/* Compute external initial wave functions: */

		if(counter<N_final)
		{
		    init_currents[counter]->evaluate();
		    ++counter;
		    return counter;
		}
		else if(counter<N_bits)
		{
		    init_currents[counter]->evaluate();
		    ++counter;
		    return counter;
		}

		/* Compute internal subamplitudes: */

		else if(counter<(interactions.size()+N_bits))
		{
		    interaction_iterator it=interactions.begin();
		    for(size_type i=0;i<(counter-N_bits);++i)
		    {
			++it;
		    }
		    it->evaluate();
		    ++counter;
		    return counter;
		}

		/* Compute final wave function: */

		else
		{
		    final_current->evaluate();
		    counter=0;
		    return counter;
		}
	    }

	    /* Resetting algorithm: */

	    void reset()
	    {
		/* Reset external currents: */

		for(size_type i=0;i<N_final;++i)
		{
		    init_currents[i]->reset();
		}
		for(size_type i=N_final+1;i<N_in+N_out;++i)
		{
		    init_currents[i-1]->reset();
		}

		/* Reset interactions, which will effectively reset the internal
		 * produced currents: */

		for(interaction_iterator it=interactions.begin();it != interactions.end();++it)
		{
		    it->reset();
		}
		if(!empty)
		{
		    (--interactions.end())->get_produced_current()->reset();
		}
		
		/* Reset the final current: */
		
		final_current->reset();
	    }

	    /* Function returning the i-th external particle phase space: */

	    phase_space_type* get_phase_space(size_type i)
	    {
		if(i<N_final)
		{
		    return init_currents[i]->get_phase_space();
		}
		if(i==N_final)
		{
		    return final_current->get_phase_space();
		}
		if(i<N_in+N_out)
		{
		    return init_currents[i-1]->get_phase_space();
		}
		log(log_level::warning)<<"requested leg "<<i<<" out of range--NULL returned"<<endlog;
		return NULL;
	    }

	    /* Function returning the i-th external const particle phase space: */

	    const phase_space_type* get_phase_space(size_type i) const
	    {
		if(i<N_final)
		{
		    return init_currents[i]->get_phase_space();
		}
		if(i==N_final)
		{
		    return final_current->get_phase_space();
		}
		if(i<N_in+N_out)
		{
		    return init_currents[i-1]->get_phase_space();
		}
		log(log_level::warning)<<"requested leg "<<i<<" out of range--NULL returned"<<endlog;
		return NULL;
	    }

	    /// Perform handlebar operation on the i-th external particle.

	    void set_handlebar(size_type i)
	    {
		if(get_phase_space(i)!=NULL)
		{
		    get_phase_space(i)->set_handlebar();
		}
	    }

	    /// Discard handlebar operation on the i-th external particle.

	    void unset_handlebar(size_type i)
	    {
		if(get_phase_space(i)!=NULL)
		{
		    get_phase_space(i)->unset_handlebar();
		}
	    }

	    /// Returns whether the i-th external particle is a handlebar.

	    bool is_handlebar(size_type i)
	    {
		if(get_phase_space(i)!=NULL)
		{
		    return get_phase_space(i)->is_handlebar();
		}
		return false;
	    }

	    /* Function returning the vector of incoming particles: */

	    vector<const particle_type*,N_in> get_phi_in()
	    {
		vector<const particle_type*,N_in>p;
		for(size_type i=0;i<N_in;++i)
		{
		    p[i]=get_phase_space(i)->particle_type;
		}
		return p;
	    }

	    /* Function returning the vector of outgoing particles: */

	    vector<const particle_type*,N_out> get_phi_out()
	    {
		vector<const particle_type*,N_out>p;
		for(size_type i=N_in;i<N_in+N_out;++i)
		{
		    p[i-N_in]=get_phase_space(i)->particle_type;
		}
		return p;
	    }

	    /* Function returning the vector of incoming mass pointers: */

	    vector<const r_value_type*,N_in> get_m_in()
	    {
		vector<const r_value_type*,N_in>m;
		for(size_type i=0;i<N_in;++i)
		{
		    m[i]=get_phase_space(i)->particle_type->get_mass_address();
		}
		return m;
	    }

	    /* Function returning the vector of outgoing mass pointers: */

	    vector<const r_value_type*,N_out> get_m_out()
	    {
		vector<const r_value_type*,N_out>m;
		for(size_type i=N_in;i<N_in+N_out;++i)
		{
		    m[i-N_in]=get_phase_space(i)->particle_type->get_mass_address();
		}
		return m;
	    }

	    /* Function returning the vector of incoming momentum pointers: */

	    vector<momentum_type*,N_in> get_p_in()
	    {
		vector<momentum_type*,N_in>p;
		for(size_type i=0;i<N_in;++i)
		{
		    p[i]=&(get_phase_space(i)->momentum());
		}
		return p;
	    }

	    /* Function returning the vector of outgoing momentum pointers: */

	    vector<momentum_type*,N_out> get_p_out()
	    {
		vector<momentum_type*,N_out>p;
		for(size_type i=N_in;i<N_in+N_out;++i)
		{
		    p[i-N_in]=&(get_phase_space(i)->momentum());
		}
		return p;
	    }

	    /* Clearance function: */

	    void clear()
	    {
		level=0;
		B.reset();
		interactions.clear();
		empty=true;
	    }

	    /* Returns the symmetry factor: */

	    r_value_type symmetry_factor() const
	    {
		return (r_value_type)1/((r_value_type)symm_factor);
	    }

	    /* Computes the symmetry factor for cross section calculations: */

	    void evaluate_symmetry_factor()
	    {
		std::multiset<int>flavs;
		if(N_final<N_in)
		{
		    for(size_type i=N_in-1;i<init_currents.size();++i)
		    {
			flavs.insert(init_currents[i]->get_particle_type()->get_flavour());
		    }
		}
		else
		{
		    for(size_type i=N_in;i<init_currents.size();++i)
		    {
			flavs.insert(init_currents[i]->get_particle_type()->get_flavour());
		    }
		    flavs.insert(final_current->get_particle_type()->get_flavour());
		}
		symm_factor=1;
		std::multiset<int>::iterator it=flavs.begin();
		int n=*it;
		size_type counter=1;
		++it;
		for(;it!=flavs.end();++it)
		{
		    if(n==*it)
		    {
			++counter;
			symm_factor*=counter;
		    }
		    else
		    {
			counter=1;
		    }
		    n=*it;
		}
	    }

	    /* Returns a const iterator to the beginning of the list of interactions
	     * defining the tree: */

	    const_interaction_iterator interactions_begin() const
	    {
		return interactions.begin();
	    }

	    /* Returns a const iterator to the end of the list of interactions defining
	     * the tree: */

	    const_interaction_iterator interactions_end() const
	    {
		return interactions.end();
	    }


	private:
	    
	    /* Integer denoting the level (nr of bits in the produced momentum
	     * channel) of the interaction under consideration: */
	    
	    size_type level;

	    /* Initial external current iterators. These are the external
	     * currents corresponding to the external legs of the process,
	     * except the final current: */

	    std::vector<current_iterator>init_currents;
	    
	    /* Final current iterator: */
	    
	    current_iterator final_current;

	    /* List of interaction objects participating in the tree: */

	    std::list<interaction_type>interactions;
	    
	    /* Utility bit string: */
	    
	    bit_string<N_bits> B;

	    /* Boolean denoting if the process is nonexistent: */

	    bool empty;
	    
	    /* Utility counter: */
	    
	    size_type counter;

	    /* Process symmetry factor: */

	    size_type symm_factor;

	    /* Maximum vertex rank in the model under consideration: */

	    static size_type max_rank;
	    
	    /* Number of flavours in the model under consideration: */
	    
	    static size_type flavours;

	    /* Final current: */

	    static std::size_t N_final;

	    /* Initialisation tag: */

	    static bool initialised;
	    
	    /* Bit string partition table: */
	    
	    std::vector< std::vector< bit_string<N_bits> > >partition;

	    /* function moving to the next helicity configuration of external
	     * particles in a spin sum: */

	    bool next_dof_configuration(std::bitset<N_in+N_out>summed_spins,std::bitset<N_in+N_out>summed_cols)
	    {
		for(size_type n=0;n<N_final;++n)
		{
		    if(!init_currents[n]->evaluate_dof_sum(summed_spins[n],summed_cols[n]))
		    {
			return true;
		    }
		}
		if(!final_current->evaluate_dof_sum(summed_spins[N_final],summed_cols[N_final]))
		{
		    return true;
		}
		for(size_type n=N_final+1;n<N_in+N_out;++n)
		{
		    if(!init_currents[n-1]->evaluate_dof_sum(summed_spins[n],summed_cols[n]))
		    {
			return true;
		    }
		}
		return false;
	    }

	    /* Iterative tree building function: builds all interactions
	     * yielding produced currents in the next momentum channel: */

	    void build_next()
	    {
		/* Go to the next momentum channel: */

		B.next();
		
		/* Maximum number of incoming currents of the to-be-constructed
		 * interactions: */
		
		size_type legs=std::min(B.count(),max_rank-1);
		    
		/* Declare an abortion tag: */
		    
		bool abort;

		for(size_type n=2;n<=legs;++n)
		{
		    /* Make partitions of the momentum channel into n pieces: */
		    
		    partition=bit_string_partition<N_bits>::make_partition(B,n);

		    /* Construct the vector of current iterators: */

		    std::vector<current_iterator>iters(n);

		    /* Construct a vector of participating particle types: */

		    std::vector<const particle_type*>partvec(n);
		    
		    /* For every partition, search all possible combinations: */

		    for(size_type j=0;j<partition.size();++j)
		    {
			abort=false;
			
			/* If a momentum bit string posseses no marked currents,
			 * break the loop and go to the next partition: */
			
			for(size_type k=0;k<n;++k)
			{
			    iters[k]=current_tree_type::first_marked_current(partition[j][k]);
			    if(iters[k]==current_tree_type::end())
			    {
				abort=true;
				break;
			    }
			}

			/* Else, attempt combinations: */

			if(!abort)
			{
			    /* Loop over all combinations of marked currents
			     * with momentum channels determined by the
			     * partition: */

			    do
			    {
				/* Fill the particle vector: */

				for(size_type k=0;k<n;++k)
				{
				    partvec[k]=iters[k]->get_produced_particle();
				}

				/* Find a fusion corresponding to the particle
				 * vector; the pair consists of the first and
				 * last iterator pointing to the fusion class
				 * matching with partvec...*/

				std::pair<fusion_iterator,fusion_iterator>iterpair=model_wrapper<model_t>::find_fusion(partvec);

				/* Construct interaction objects from fusions:
				 * */

				for(fusion_iterator f_iter=iterpair.first; f_iter!=iterpair.second;++f_iter)
				{
				    /* Update the current level: */

				    level=B.count();

				    /* Reconstruct the produced current: */

				    current_iterator prod_curr=current_tree_type::find_current(B,f_iter->second.get_produced_particle());

				    /* Mark it for future iterations: */

				    prod_curr->mark();

				    /* Construct the interaction from the fusion
				     * class data: */

				    interaction_type node(f_iter->second.get_vertex(),f_iter->second.get_leg());

				    /* Copy the incoming current iterators: */

				    std::vector<current_iterator>in_currents(iters);

				    /* Reshuffle the incoming current iterators
				     * to match the vertex legs: */

				    f_iter->second.apply_ordering(in_currents);

				    /* Insert the produced current iterator: */

				    in_currents.insert(in_currents.begin()+f_iter->second.get_leg(),prod_curr);

				    /* Insert them into the interaction object:
				     * */

				    node.insert_currents(in_currents);

				    /* Let the interaction fetch its Feynman
				     * rule: */

				    node.fetch_Feynman_rule();

				    /* Add the interaction to the list: */

				    interactions.push_back(node);
				}
			    }
			    while(current_tree_type::next_marked_currents(iters));
			}
		    }
		}
	    }

	    /* Tree building function for the final interactions (with maximal
	     * bit string): */

	    void build_final()
	    {
		/* Go to the next momentum channel: */

		B.next();
		
		/* Maximum number of incoming currents of the to-be-constructed
		 * interactions: */
		
		size_type legs=std::min(B.count(),max_rank-1);

		/* Declare an abortion tag: */

		bool abort;

		for(size_type n=2;n<=legs;++n)
		{
		    /* Make partitions of the momentum channel into n pieces: */
		    
		    partition=bit_string_partition<N_bits>::make_partition(B,n);

		    /* Construct the vector of current iterators: */

		    std::vector<current_iterator>iters(n);

		    /* Construct a vector of participating particle types: */

		    std::vector<const particle_type*>partvec(n);

		    for(size_type j=0;j<partition.size();++j)
		    {
			abort=false;
			
			/* If a momentum bit string posseses no marked currents,
			 * break the loop and go to the next partition: */
			
			for(size_type k=0;k<n;++k)
			{
			    iters[k]=current_tree_type::first_marked_current(partition[j][k]);
			    if(iters[k]==current_tree_type::end())
			    {
				abort=true;
				break;
			    }
			}

			/* Else, attempt combinations: */

			if(!abort)
			{
			    do
			    {
				/* Fill the particle vector: */

				for(size_type k=0;k<n;++k)
				{
				    partvec[k]=iters[k]->get_produced_particle();
				}

				/* Find a fusion corresponding to the particle
				 * vector; the pair consists of the first and
				 * last iterator pointing to the fusion class
				 * matching with partvec...*/

				std::pair<fusion_iterator,fusion_iterator>iterpair= model_wrapper<model_t>::find_fusion(partvec);

				/* Construct interaction objects from fusions:
				 * */

				for(fusion_iterator f_iter=iterpair.first; f_iter!=iterpair.second;++f_iter)
				{
				    /* Proceed only if the produced particle can
				     * annihilate with the final particle: */

				    if(f_iter->second.get_produced_particle() == final_current->get_produced_particle()->get_anti_particle())
				    {
					/* The tree contains at least one
					 * diagram. */

					empty=false;

					/* Update the current level: */

					level=B.count();
				    
					/* Reconstruct the produced current: */

					current_iterator prod_curr = current_tree_type::find_current(B,f_iter->second.get_produced_particle());

					/* Mark it: */

					prod_curr->mark();

					/* Construct the interaction from the fusion
					 * class data: */

					interaction_type node(f_iter->second.get_vertex(),f_iter->second.get_leg());

					/* Copy the incoming current iterators: */

					std::vector<current_iterator>in_currents(iters);

					/* Reshuffle the incoming current iterators
					 * to match the vertex legs: */

					f_iter->second.apply_ordering(in_currents);

					/* Insert the produced current iterator: */

					in_currents.insert(in_currents.begin()+f_iter->second.get_leg(),prod_curr);

					/* Insert them into the interaction object:
					 * */

					node.insert_currents(in_currents);

					/* Let the interaction fetch its Feynman
					 * rule: */

					node.fetch_Feynman_rule();

					/* Add the interaction to the list: */

					interactions.push_back(node);
				    }
				}
			    }
			    while(current_tree_type::next_marked_currents(iters));
			}
		    }
		}
	    }

	    /* Function marking all currents participating in the tree: */

	    void mark_currents()
	    {
		for(size_type i=0;i<N_bits;++i)
		{
		    init_currents[i]->mark();
		}
		for(interaction_iterator it=interactions.begin();it != interactions.end();++it)
		{
		    it->produced_current->mark();
		}
		final_current->mark();
	    }

	    /* Function unmarking all currents participating in the tree: */

	    void unmark_currents()
	    {
		for(size_type i=0;i<N_bits;++i)
		{
		    init_currents[i]->unmark();
		}
		for(interaction_iterator it=interactions.begin();it != interactions.end();++it)
		{
		    it->produced_current->unmark();
		}
		final_current->unmark();
	    }
    };
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t process_tree<model_t,N_in,N_out>::N_bits;
    template<class model_t,std::size_t N_in,std::size_t N_out>const bool process_tree<model_t,N_in,N_out>::decomposes;
    template<class model_t,std::size_t N_in,std::size_t N_out>typename process_tree<model_t,N_in,N_out>::size_type process_tree<model_t,N_in,N_out>::max_rank=0;
    template<class model_t,std::size_t N_in,std::size_t N_out>typename process_tree<model_t,N_in,N_out>::size_type process_tree<model_t,N_in,N_out>::flavours=0;
    template<class model_t,std::size_t N_in,std::size_t N_out>std::size_t process_tree<model_t,N_in,N_out>::N_final=0;
    template<class model_t,std::size_t N_in,std::size_t N_out>bool process_tree<model_t,N_in,N_out>::initialised=false;

    /* Function requesting whether the tree is empty: */

    template<class model_t,std::size_t N_in,std::size_t N_out>bool is_empty(const process_tree<model_t,N_in,N_out>& t)
    {
	return t.is_empty();
    }
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_PROCESS_TREE_H_*/

