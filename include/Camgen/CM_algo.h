//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file CM_algo.h
    \brief CM_algorithm class template interface and implementation header.
 */

#ifndef CAMGEN_CM_ALGO_H_
#define CAMGEN_CM_ALGO_H_

#include <Camgen/process.h>
#include <Camgen/license_print.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the CM_algorithm class template. The CM_algorithm constitutes   *
 * Camgen's user interface for computing process amplitudes. It creates a list  *
 * of process trees, which can be selected via the set_process member function,  *
 * and then evaluated by calling evaluate(). It allows access to the external    *
 * particle's phase space by calling get_phase_space(), swapping degrees of      *
 * freedom, counting diagrams and memory usage etc.                              *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Amplitude-computing algorithm class template.
    /** Core class of the Camgen library for evaluating scattering amplitudes of
     * arbitrary processes in any QFT. The first template parameter denotes the
     * model class, the second and third the respective number of incoming and
     * outgoing particles. */

    template<class model_t,std::size_t N_in,std::size_t N_out>class CM_algorithm
    {
	public:
	    
	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	    
	    /* Process tree type definition: */
	    
	    typedef process_tree<model_t,N_in,N_out> tree_type;

	    /* Process type definition: */

	    typedef process<model_t,N_in,N_out> process_type;

	    /* Particle phase space type definition: */

	    typedef typename tree_type::phase_space_type phase_space_type;
	    
	    /* Current tree type definition: */

	    typedef typename tree_type::current_tree_type current_tree_type;

	    /* Tree iterator type definitions: */

	    typedef typename std::list<tree_type>::iterator tree_iterator;
	    typedef typename std::list<tree_type>::const_iterator const_tree_iterator;

	    /* Flavour vector type definition: */

	    typedef vector<size_type,N_in+N_out> flavour_vector;
	    typedef vector<int,N_in+N_out> pdg_id_vector;

	    /* Process iterator type definitions: */

	    typedef typename std::list<process_type>::iterator process_iterator;
	    typedef typename std::list<process_type>::const_iterator const_process_iterator;

	    /* Number of incoming and outgoing particles: */

	    static const std::size_t N_incoming=N_in;
	    static const std::size_t N_outgoing=N_out;
	    static const std::size_t N_external=N_in+N_out;

	    /// Trivial constructor.
	    /// The optional integer argument denotes the external particle
	    /// chosen by the algorithm as the final current in the recursion.
	    /// The amplitudes should not depend on this choice, unless there is
	    /// a bug in Camgen, or the model violates CPT invariance.  

	    CM_algorithm(std::size_t n=0):sorted_by_flavour(false),sorted_by_pdg_id(false)
	    {
		license_print::initialise();
		if(n>=N_external)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid choice "<<n<<" for final particle in "<<N_in<<" -> "<<N_out<<" process--proceeding with "<<N_external-1<<endlog;
		}
		N_final=std::min(n,N_external-1);
		tree_type::initialise(N_final);
		reset_ordering();
	    }

	    /// Constructor specifying the process.
	    /// The first argument should be of the form "phi1,...,phiN_in >
	    /// psi1,...,psiN_out", and the second argument is the optional final
	    /// current (see trivial constructor).

	    CM_algorithm(const std::string& str,std::size_t n=0):sorted_by_flavour(false),sorted_by_pdg_id(false)
	    {
		license_print::initialise();
		if(n>=N_external)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid choice "<<n<<" for final particle in "<<N_in<<" -> "<<N_out<<" process--proceeding with "<<N_external-1<<endlog;
		}
		N_final=std::min(n,N_external-1);
		tree_type::initialise(N_final);
		reset_ordering();
		process_type::add_process(processes,str);
	    }

	    /// Subprocess insertion method.
	    /// The argument should be of the form "phi1,...,phiN_in >
	    /// psi1,...,psiN_out". If the insertion was succesful, the function
	    /// returns true and moves the process iterator members to this subprocess entry.

	    bool add_process(const std::string& str)
	    {
		return process_type::add_process(processes,str);
	    }

	    /// Function loading the external particles.
	    /// This routine sorts the subprocess list by lexicographical
	    /// comparison of external particle flavours and removes identical
	    /// initial + final states. Then, it initializes a vertex tree for
	    /// each subprocess.

	    void load()
	    {
		/* Remove identical processes: */

		processes.sort();
		processes.unique();

		/* Construct trees and unordered processes for the remaining
		 * processes: */

		for(process_iterator it=processes.begin();it!=processes.end();++it)
		{
		    it->add_tree(trees);
		}
	
		/* Assign the process and tree iterators: */

		process_it=processes.begin();
		tree_it=process_it->get_tree();
	    }

	    /// Rearranges the process list by flavours.
	    /// This may speed up multi-process evaluations if the processes are
	    /// searched by internal Camgen flavours (=order of particle
	    /// insertion in the model class).

	    void sort_by_flavours()
	    {
		if(!sorted_by_flavour)
		{
		    processes.sort();
		    sorted_by_flavour=true;
		    sorted_by_pdg_id=false;
		}
	    }

	    /// Rearranges the process list by PDG id.
	    /// This may speed up multi-process evaluations if the processes are
	    /// searched by PDG id.

	    void sort_by_pdg_id()
	    {
		if(!sorted_by_pdg_id)
		{
		    processes.sort(pdg_id_comp<model_t,N_in,N_out>);
		    sorted_by_flavour=false;
		    sorted_by_pdg_id=true;
		}
	    }

	    /// Selects a subprocess by string argument.
	    /// The argument should be of the same form as in the constructor.
	    /// Because this string has to be decomposed and matched with
	    /// particles, this selection doesn't have optimal performance. For
	    /// faster selection, one should choose either sort_by_flavour or
	    /// sort_by_PDG_id and then use respectively set_process_flavours or
	    /// set_process_pdg_ids. Note that if the process is not found in
	    /// the list, the process and corresponding vertex tree will be
	    /// added automatically.

	    process_iterator set_process(const std::string& proc)
	    {
		vector<std::string,N_in+N_out>namevec;
		if(!decompose_process<N_in,N_out>(proc,namevec))
		{
		    process_it=processes.end();
		    tree_it=trees.end();
		    return process_it;
		}
		flavour_vector fv;
		const particle<model_t>* phi;
		for(size_type i=0;i<N_in+N_out;++i)
		{
		    phi=model_wrapper<model_t>::get_particle(namevec[i]);
		    if(phi==NULL)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed to set process "<<proc<<": particle "<<namevec[i]<<" not found in model"<<endlog;
			process_it=processes.end();
			tree_it=trees.end();
			return process_it;
		    }
		    else
		    {
			fv[i]=phi->get_flavour();
		    }
		}
		return set_process_flavours(fv);
	    }

	    /// Selects a process by flavour vector. If the list is sorted by
	    /// flavour, the search is optimal. If the process is not found in
	    /// the list, a subprocess and vertex tree will be added
	    /// automatically.

	    process_iterator set_process_flavours(const flavour_vector& fv)
	    {
		/* Construct a sorted copy of the argument: */

		flavour_vector sorted_fv(fv);
		std::sort(sorted_fv.begin(),sorted_fv.begin()+N_in);
		std::sort(sorted_fv.begin()+N_in,sorted_fv.end());
		
		/* Look for the sorted flavours in the process list: */

		if(sorted_by_flavour)
		{
		    process_it=std::lower_bound(processes.begin(),processes.end(),sorted_fv,flav_vec_comp<model_t,N_in,N_out>);
		}
		else
		{
		    process_it=processes.begin();
		    for(;process_it!=processes.end();++process_it)
		    {
			if(!(process_it->get_sorted_flavours()<sorted_fv)){break;}
		    }
		}

		/* Reset the ordering vector to 'no ordering' mode: */

		reset_ordering();
		
		/* If the lower bound search ends within the tree, there are two
		scenario's: */
		
		if(process_it!=processes.end())
		{
		    /* If the process was found, set the ordering to the correct
		     * values and move the tree iterator to the corresponding
		     * tree: */

		    tree_it=process_it->get_tree();
		    if(process_it->get_sorted_flavours()==sorted_fv)
		    {
			permutation<flavour_vector>perm(fv,process_it->get_flavours());
			perm(ordering);
			return process_it;
		    }

		    /* Else, attempt to insert the process in the list: */

		    process_it=process<model_t,N_in,N_out>::insert_process_by_flavour(processes,process_it,fv);
		    
		    /* If the insertion was succesful, insert a corresponding
		     * tree: */
		    
		    if(process_it!=processes.end())
		    {
			tree_it=process_it->insert_tree(trees,tree_it);
			tree_it->build();
			tree_it->clean();
			tree_it->set_Fermi_signs();
			tree_it->initialise_currents();
			tree_it->assign_momenta();
			return process_it;
		    }

		    /* If the insertion was not succesful, set the iterators to
		     * the end of the lists: */

		    tree_it=trees.end();
		    return process_it;
		}

		/* If the search result outside the process list, attempt to
		 * insert the process at the back: */

		process_it=process<model_t,N_in,N_out>::insert_process_by_flavour(processes,process_it,fv);
		
		/* If the insertion was succesful, insert the corresponding tree
		 * at the back of the tree list: */
		
		if(process_it != processes.end())
		{
		    tree_it=process_it->insert_tree(trees,trees.end());
		    tree_it->build();
		    tree_it->clean();
		    tree_it->set_Fermi_signs();
		    tree_it->initialise_currents();
		    tree_it->assign_momenta();
		}

		/* Else, set the tree iterator outside the list: */

		else
		{
		    tree_it=trees.end();
		}
		return process_it;
	    }

	    /// Selects a process by flavour PDG id. If the list is sorted by
	    /// PDG id's, the search is optimal. If the process is not found in
	    /// the list, a subprocess and vertex tree will be added
	    /// automatically.

	    process_iterator set_process_pdg_ids(const vector<int,N_external>& pv)
	    {
		/* Construct a sorted copy of the argument: */

		vector<int,N_external> sorted_pv(pv);
		std::sort(sorted_pv.begin(),sorted_pv.begin()+N_in);
		std::sort(sorted_pv.begin()+N_in,sorted_pv.end());
		
		/* Look for the sorted pdg id's in the process list: */
		
		if(sorted_by_pdg_id)
		{
		    process_it=std::lower_bound(processes.begin(),processes.end(),sorted_pv,pdg_id_vec_comp<model_t,N_in,N_out>);
		}
		else
		{
		    process_it=processes.begin();
		    for(;process_it!=processes.end();++process_it)
		    {
			if(!(process_it->get_sorted_pdg_ids()<sorted_pv)){break;}
		    }
		}

		/* Reset the ordering vector to 'no ordering' mode: */

		reset_ordering();
		
		/* If the lower bound search ends within the tree, there are two
		scenario's: */
		
		if(process_it!=processes.end())
		{
		    /* If the process was found, set the ordering to the correct
		     * values and move the tree iterator to the corresponding
		     * tree: */

		    tree_it=process_it->get_tree();
		    if(process_it->get_sorted_pdg_ids()==sorted_pv)
		    {
			permutation< vector<int,N_external> >perm(pv,process_it->get_pdg_ids());
			perm(ordering);
			tree_it=process_it->get_tree();
			return process_it;
		    }

		    /* Else, attempt to insert the process in the list: */

		    process_it=process<model_t,N_in,N_out>::insert_process_by_pdg_id(processes,process_it,pv);
		    
		    /* If the insertion was succesful, insert a corresponding
		     * tree: */
		    
		    if(process_it!=processes.end())
		    {
			tree_it=process_it->insert_tree(trees,tree_it);
			tree_it->build();
			tree_it->clean();
			tree_it->set_Fermi_signs();
			tree_it->initialise_currents();
			tree_it->assign_momenta();
			return process_it;
		    }

		    /* If the insertion was not succesful, set the iterators to
		     * the end of the lists: */

		    tree_it=trees.end();
		    return process_it;
		}

		/* If the search result outside the process list, attempt to
		 * insert the process at the back: */

		process_it=process<model_t,N_in,N_out>::insert_process_by_pdg_id(processes,process_it,pv);
		if(process_it != processes.end())
		{
		    tree_it=process_it->insert_tree(trees,trees.end());
		    tree_it->build();
		    tree_it->clean();
		    tree_it->set_Fermi_signs();
		    tree_it->initialise_currents();
		    tree_it->assign_momenta();
		}

		/* Else, set the tree iterator outside the list: */

		else
		{
		    tree_it=trees.end();
		}
		return process_it;
	    }

	    /// If the process iterator argument points to a process within the list,
	    /// sets the tree and process iterator to the corresponding tree,
	    /// otherwise no action is taken.

	    process_iterator set_process(process_iterator it)
	    {
		tree_iterator it2=trees.begin();
		process_iterator it3=processes.begin();
		while(it3!=it and it3!=processes.end())
		{
		    ++it2;
		    ++it3;
		}
		if(it3!=processes.end())
		{
		    tree_it=it2;
		    process_it=it3;
		}
		return process_it;
	    }

	    /// If the tree iterator argument points to a tree within the list,
	    /// sets the tree and process iterator to the corresponding tree,
	    /// otherwise no action is taken.

	    process_iterator set_process(tree_iterator it)
	    {
		tree_iterator it2=trees.begin();
		process_iterator it3=processes.begin();
		while(it2!=it and it2!=trees.end())
		{
		    ++it2;
		    ++it3;
		}
		if(it2!=trees.end())
		{
		    tree_it=it2;
		    process_it=it3;
		}
		return process_it;
	    }

	    /// Resets the process to the first in the list.

	    bool reset_process()
	    {
		tree_it=trees.begin();
		process_it=processes.begin();
		return (tree_it!=trees.end());
	    }

	    /// Moves to the next process in the list. Returns a boolean
	    /// denoting whether this is the final process.

	    bool next_process()
	    {
		++tree_it;
		++process_it;
		return (tree_it!=trees.end());
	    }

	    /// Returns whether the process is valid.

	    bool valid_process() const
	    {
		return (process_it!=processes.end());
	    }

	    /// Returns the iterator to the current vertex tree.

	    tree_iterator get_tree_iterator()
	    {
		return tree_it;
	    }

	    /// Returns the (const) iterator to the current subprocess.

	    process_iterator get_process_iterator()
	    {
		return process_it;
	    }

	    /// Returns the (const) iterator to the current vertex tree.

	    const_tree_iterator get_tree_iterator() const
	    {
		return tree_it;
	    }

	    /// Returns the (const) iterator to the current subprocess.

	    const_process_iterator get_process_iterator() const
	    {
		return process_it;
	    }

	    /// Returns the flavour of the i-th external particle in the current
	    /// subprocess.

	    int get_flavour(std::size_t i)
	    {
		if(process_it!=processes.end())
		{
		    return process_it->get_flavour(ordering[i]);
		}
		return -1;
	    }

	    /// Returns the PDG id of the i-th external particle in the current
	    /// subprocess.

	    int get_pdg_id(std::size_t i)
	    {
		if(process_it!=processes.end())
		{
		    return process_it->get_pdg_id(ordering[i]);
		}
		return -1;
	    }

	    /* Tree-building initialiser function: */

	    void init_construct()
	    {
		tree_it->init_build();
	    }

	    /* Tree-building iteration step: */

	    void construct_next_level()
	    {
		tree_it->build_next();
	    }

	    /// Constructs the vertex tree of the current subprocess.

	    void construct()
	    {
		if(tree_it!=trees.end())
		{
		    tree_it->build();
		    tree_it->clean();
		    tree_it->set_Fermi_signs();
		    tree_it->initialise_currents();
		    tree_it->assign_momenta();
		    tree_it->compute_coupling_flags();
		}
	    }

	    /// Constructs the vertex trees of all subprocesses.

	    void construct_trees()
	    {
		for(tree_iterator it=trees.begin();it!=trees.end();++it)
		{
		    it->build();
		    it->clean();
		}
		for(tree_iterator it=trees.begin();it != trees.end();++it)
		{
		    it->set_Fermi_signs();
		    it->initialise_currents();
		    it->assign_momenta();
		    it->compute_coupling_flags();
		}
		tree_it=trees.begin();
	    }

	    /// Removes current process from tree.

	    tree_iterator remove_process()
	    {
		if(tree_it!=trees.end())
		{
		    processes.erase(process_it);
		    trees.erase(tree_it);
		}
		return tree_it;
	    }

	    /// Removes subprocesses with empty vertex trees.

	    void remove_empty_processes()
	    {
		processes.remove_if(has_empty_tree<model_t,N_in,N_out>);
		trees.remove_if(is_empty<model_t,N_in,N_out>);
	    }

	    /// Counts indistinguishable processes.

	    size_type n_processes() const
	    {
		return processes.size();
	    }

	    /// Counts current trees.

	    size_type n_trees() const
	    {
		return trees.size();
	    }

	    /// Counts indistinguishable processes with nonempty vertex trees.

	    size_type count_nonempty_processes() const
	    {
		size_type n=0;
		for(const_tree_iterator it=trees.begin();it!=trees.end();++it)
		{
		    if(!(it->is_empty()))
		    {
			++n;
		    }
		}
		return n;
	    }

	    /// Counts Feynman diagrams of current subprocess.

	    long long unsigned count_diagrams()
	    {
		if(tree_it!=trees.end())
		{
		    return tree_it->count_diagrams();
		}
		return 0;
	    }

	    /// Counts the recursive complexity, i.e. the number of current
	    /// combinations in the process tree.

	    long long unsigned count_current_combinations()
	    {
		if(tree_it!=trees.end())
		{
		    return tree_it->count_current_combinations();
		}
		return 0;
	    }

	    /// Counts Feynman diagrams of all subproccesses.

	    long long unsigned count_all_diagrams()
	    {
		long long unsigned n=0;
		for(tree_iterator it=trees.begin();it != trees.end();++it)
		{
		    n+=(it->count_diagrams());
		}
		return n;
	    }

	    /// Counts the recursive complexity, i.e. the number of current
	    /// combinations in all process trees.

	    long long unsigned count_all_current_combinations()
	    {
		long long unsigned n=0;
		for(tree_iterator it=trees.begin();it != trees.end();++it)
		{
		    n+=(it->count_current_combinations());
		}
		return n;
	    }
	    
	    /// Evaluates current subprocess amplitude

	    value_type evaluate()
	    {
		if(tree_it != trees.end())
		{
		    tree_it->reset();
	            return tree_it->evaluate();
		}
		return value_type(0,0);
	    }
	    
	    /// Evaluates current subprocess amplitude squared

	    r_value_type evaluate2()
	    {
		if(tree_it != trees.end())
		{
		    tree_it->reset();
	            return std::norm(tree_it->evaluate());
		}
		return (r_value_type)0;
	    }

	    /// Evaluates spin-summed subprocess amplitude and returns the spin-summed squared amplitude.

	    r_value_type evaluate_spin_sum()
	    {
		std::bitset<N_in+N_out>s;
		s.set();
		std::bitset<N_in+N_out>c;
		c.reset();
		if(tree_it != trees.end())
		{
		    tree_it->reset();
		    return tree_it->evaluate(s,c);
		}
		return (r_value_type)0;
	    }

	    /// Evaluates colour-summed subprocess amplitude and returns the spin-summed squared amplitude.

	    r_value_type evaluate_colour_sum()
	    {
		std::bitset<N_in+N_out>s;
		s.reset();
		std::bitset<N_in+N_out>c;
		c.set();
		if(tree_it != trees.end())
		{
		    tree_it->reset();
		    return tree_it->evaluate(s,c);
		}
		return (r_value_type)0;
	    }

	    /// Evaluates the summed amplitude where inidividual particle
	    /// helicity and colour summation flags control the summing.

	    r_value_type evaluate_sum()
	    {
		if(tree_it != trees.end())
		{
		    tree_it->reset();
		    return tree_it->evaluate(summed_spins,summed_cols);
		}
		return (r_value_type)0;
	    }

	    /// Returns phase space object address of the i-th external particle
	    /// in the current subprocess.

	    phase_space_type* get_phase_space(size_type i)
	    {
		if(tree_it != trees.end())
		{
		    return tree_it->get_phase_space(ordering[i]);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requesting phase space of empty process--NULL returned."<<endlog;
		    return NULL;
		}
	    }

	    /// Returns const phase space object address of the i-th external particle
	    /// in the current subprocess.

	    const phase_space_type* get_phase_space(size_type i) const
	    {
		if(tree_it != trees.end())
		{
		    return tree_it->get_phase_space(ordering[i]);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"requesting phase space of empty process--NULL returned"<<endlog;
		    return NULL;
		}
	    }

	    /* Evaluates the next current combination in the subprocess: */

	    void evaluate_next()
	    {
		if(tree_it != trees.end())
		{
		    tree_it->evaluate_next();
		}
	    }

	    /// Prints the current vertex tree and interacting off-shell
	    /// currents tensors to the argument streaming object.

	    std::ostream& print(std::ostream& os) const
	    {
		tree_it->print(os);
		return os;
	    }

	    /// Prints the current vertex tree without current tensors to the
	    /// argument streaming object.

	    std::ostream& print_tree(std::ostream& os) const
	    {
		tree_it->print_tree(os);
		return os;
	    }

	    /// Prints the process list to the argument streaming object.

	    std::ostream& print_processes(std::ostream& os) const
	    {
		for(const_process_iterator it=processes.begin();it!=processes.end();++it)
		{
		    it->print(os);
		    os<<std::endl;
		}
		return os;
	    }

	    /// Prints the current process to the argument streaming object.

	    std::ostream& print_process(std::ostream& os) const
	    {
		if(process_it != processes.end())
		{
		    return tree_it->print_process(os);
		}
		return os;
	    }

	    /// Prints the participating vertices to the argument streaming object.

	    std::ostream& print_vertices(std::ostream& os) const
	    {
		if(tree_it != trees.end())
		{
		    return tree_it->print_vertices(os);
		}
		return os;
	    }

	    /* Output method printing the currents transformed to the
	     * colour-flow basis for gluons: */

	    template<class rep_t>std::ostream& print_cf_transform(std::ostream& os)
	    {
		process_it->print(os);
		tree_it->template print_cf_transform<rep_t>(os);
		return os;
	    }

	    /// Switches on spin summation for particle i.

	    void sum_spin(size_type i)
	    {
		summed_spins.set(i);
	    }

	    /// Sums all spins when calling evaluate_spin_sum.

	    void sum_spins()
	    {
		summed_spins.set();
	    }

	    /// Switches off spin summation for particle i.

	    void nsum_spin(size_type i)
	    {
		summed_spins.reset(i);
	    }

	    /// Switches off spin summation for all particles.

	    void nsum_spins()
	    {
		summed_spins.reset();
	    }

	    /// Switches on colour summation for particle i.

	    void sum_colour(size_type i)
	    {
		summed_cols.set(i);
	    }

	    /// Sums all colours when calling evaluate_spin_sum.

	    void sum_colours()
	    {
		summed_cols.set();
	    }

	    /// Switches off colour summation for particle i.

	    void nsum_colour(size_type i)
	    {
		summed_cols.reset(i);
	    }

	    /// Switches off spin summation for all particles.

	    void nsum_colours()
	    {
		summed_cols.reset();
	    }

	    /// Refreshes the current tree and reconstructs the subprocess trees.

	    void refresh()
	    {
		current_tree_type::refresh();
		for(tree_iterator it=trees.begin();it!=trees.end();++it)
		{
		    it->clear();
		}
		construct_trees();
	    }

	    /// Clears all data.

	    void clear()
	    {
		trees.clear();
		processes.clear();
	    }

	protected:
	    
	    /* List of subprocess trees: */

	    std::list<tree_type>trees;

	    /* List of corresponding processes: */
	    
	    std::list<process_type> processes;

	    /* Boolean denoting whether the processes are ordered by
	     * flavour configuration: */

	    bool sorted_by_flavour;

	    /* Boolean denoting whether the processes are ordered by
	     * pdg id: */

	    bool sorted_by_pdg_id;
	    
	    /* Iterator across the tree list: */
	    
	    tree_iterator tree_it;

	    /* Iterator across the process list: */

	    process_iterator process_it;

	    /* Ordering of the external particles w.r.t. the underlying ordered
	     * process: */

	    vector<std::size_t,N_external>ordering;

	    /* Helicity summation information: */

	    std::bitset<N_external>summed_spins;

	    /* Colour summation information: */

	    std::bitset<N_external>summed_cols;

	    /* Function setting the ordering to its default value, denoting an
	     * ordered process: */

	    void reset_ordering()
	    {
		for(size_type n=0;n<N_external;++n)
		{
		    ordering[n]=n;
		}
	    }

	    /* Final particle in current tree: */

	    static std::size_t N_final;
    };
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t CM_algorithm<model_t,N_in,N_out>::N_incoming;
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t CM_algorithm<model_t,N_in,N_out>::N_outgoing;
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t CM_algorithm<model_t,N_in,N_out>::N_external;
    template<class model_t,std::size_t N_in,std::size_t N_out>std::size_t CM_algorithm<model_t,N_in,N_out>::N_final=0;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_CM_ALGO_H_*/

