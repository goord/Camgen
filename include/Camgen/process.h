//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PROCESS_H_
#define CAMGEN_PROCESS_H_

#include <algorithm>
#include <Camgen/decomp_proc.h>
#include <Camgen/particle.h>
#include <Camgen/process_tree.h>
#include <Camgen/model.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Process class definition and declaration. The process class only contains a   *
 * vector of initial state particle addresses and final state particle           *
 * addresses. For lookup speed it also contains multisets of flavour and pdg-id  *
 * integers. Adding a particle name to a process may return an entire list of    *
 * processes as it also accepts particle family names.                           *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out>class process
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef particle<model_t> particle_type;
	    typedef typename get_basic_types<model_t>::size_type size_type;
	    typedef process_tree<model_t,N_in,N_out> tree_type;
	    typedef typename std::list<tree_type>::iterator tree_iterator;
	    typedef typename model_t::value_type r_value_type;

	    /* Number of in- and outgoing particles: */

	    static const std::size_t N_incoming=N_in;
	    static const std::size_t N_outgoing=N_out;
	    static const std::size_t N_external=N_in+N_out;

	    /* Static function adding processes to a list by a string input: */

	    static bool add_process(std::list< process<model_t,N_in,N_out> >& procs,const std::string& str)
	    {
		vector<std::string,N_in+N_out>namevec;
		if(!decompose_process<N_in,N_out>(str,namevec))
		{
		    return false;
		}
		
		/* Constructing vectors of corresponding incoming particles: */

		std::vector<vector<const particle_type*,N_in> >Istates(1);
		const particle_type* phi;
		for(typename std::vector<std::string>::size_type i=0;i<N_in;++i)
		{
		    if(model_wrapper<model_t>::contains_particle(namevec[i]))
		    {
			phi=model_wrapper<model_t>::get_particle(namevec[i]);
			for(typename std::vector<vector<const particle_type*,N_in> >::size_type j=0;j<Istates.size();++j)
			{
			    Istates[j][i]=phi;
			}
		    }
		    else if(model_wrapper<model_t>::contains_family(namevec[i]))
		    {
			const particle_family<model_t>* Ijet=model_wrapper<model_t>::get_family(namevec[i]);
			typename std::vector<vector<const particle_type*,N_in> >::size_type In=Istates.size();
			Istates.resize(In*Ijet->size());
			typename std::vector<vector<const particle_type*,N_in> >::size_type Il=0;
			for(typename particle_family<model_t>::size_type j=0;j<Ijet->size();++j)
			{
			    for(typename std::vector<vector<const particle_type*,N_in> >::size_type k=0;k<In;++k)
			    {
				for(typename vector<const particle_type*,N_in>::size_type l=0;l<i;++l)
				{
				    Istates[Il][l]=Istates[k][l];
				}
				Istates[Il][i]=(*Ijet)[j];
				++Il;
			    }
			}
		    }
		    else
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle or family "<<namevec[i]<<" not found, ignoring process "<<str<<endlog;
			return false;
		    }
		}
		
		/* Constructing vectors of corresponding outgoing particles: */

		std::vector<vector<const particle_type*,N_out> >Fstates(1);
		for(typename std::vector<std::string>::size_type i=0;i<N_out;++i)
		{
		    if(model_wrapper<model_t>::contains_particle(namevec[N_in+i]))
		    {
			const particle_type* phi=model_wrapper<model_t>::get_particle(namevec[N_in+i]);
			for(typename std::vector<vector<const particle_type*,N_out> >::size_type j=0;j<Fstates.size();++j)
			{
			    Fstates[j][i]=phi;
			}
		    }
		    else if(model_wrapper<model_t>::contains_family(namevec[N_in+i]))
		    {
			const particle_family<model_t>* Fjet=model_wrapper<model_t>::get_family(namevec[N_in+i]);
			typename std::vector<vector<const particle_type*,N_out> >::size_type Fn=Fstates.size();
			Fstates.resize(Fn*Fjet->size());
			typename std::vector<vector<const particle_type*,N_out> >::size_type Fl=0;
			
			for(typename particle_family<model_t>::size_type j=0;j<Fjet->size();++j)
			{
			    for(typename std::vector<vector<const particle_type*,N_out> >::size_type k=0;k<Fn;++k)
			    {
				for(typename vector<const particle_type*,N_out>::size_type l=0;l<i;++l)
				{
				    Fstates[Fl][l]=Fstates[k][l];
				}
				Fstates[Fl][i]=(*Fjet)[j];
				++Fl;
			    }
			}
		    }
		    else
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle or family "<<namevec[N_in+i]<<" not found, ignoring process "<<str<<endlog;
			return false;
		    }
		}

		/* Adding the constructed processes to the given list: */

		for(typename std::vector< vector<const particle_type*,N_in> >::size_type i=0;i<Istates.size();++i)
		{
		    for(typename std::vector< vector<const particle_type*,N_out> >::size_type j=0;j<Fstates.size();++j)
		    {
			procs.push_back(process<model_t,N_in,N_out>(Istates[i],Fstates[j]));
		    }
		}
		return true;
	    }

	    /* Static function inserting a process to a list by flavour-vector
	     * input: */

	    static typename std::list< process<model_t,N_in,N_out> >::iterator insert_process_by_flavour(std::list< process<model_t,N_in,N_out> >& processes,typename std::list< process<model_t,N_in,N_out> >::iterator proc_iter,const vector<size_type,N_external>& fv)
	    {
		process<model_t,N_in,N_out>p;
		const particle_type* phi;
		for(std::size_t i=0;i<N_in;++i)
		{
		    phi=model_wrapper<model_t>::get_particle(fv[i]);
		    if(phi==NULL)
		    {
			return processes.end();
		    }
		    p.IS_particles[i]=phi;
		    p.unsorted_pdg_ids[i]=phi->get_pdg_id();
		}
		for(std::size_t i=0;i<N_out;++i)
		{
		    phi=model_wrapper<model_t>::get_particle(fv[N_in+i]);
		    if(phi==NULL)
		    {
			return processes.end();
		    }
		    p.FS_particles[i]=phi;
		    p.unsorted_pdg_ids[N_in+i]=phi->get_pdg_id();
		}
		p.sorted_pdg_ids=p.unsorted_pdg_ids;
		std::sort(p.sorted_pdg_ids.begin(),p.sorted_pdg_ids.begin()+N_in);
		std::sort(p.sorted_pdg_ids.begin()+N_in,p.sorted_pdg_ids.end());
		p.unsorted_flavours=fv;
		p.sorted_flavours=fv;
		std::sort(p.sorted_flavours.begin(),p.sorted_flavours.begin()+N_in);
		std::sort(p.sorted_flavours.begin()+N_in,p.sorted_flavours.end());
		return processes.insert(proc_iter,p);
	    }

	    /* Static function inserting a process in a list by pdg-id vector
	     * input: */

	    static typename std::list< process<model_t,N_in,N_out> >::iterator insert_process_by_pdg_id(std::list< process<model_t,N_in,N_out> >& processes,typename std::list< process<model_t,N_in,N_out> >::iterator proc_iter,const vector<int,N_external>& pdgv)
	    {
		process<model_t,N_in,N_out>p;
		const particle_type* phi;
		for(std::size_t i=0;i<N_in;++i)
		{
		    phi=model_wrapper<model_t>::get_particle_by_pdg_id(pdgv[i]);
		    if(phi==NULL)
		    {
			return processes.end();
		    }
		    p.IS_particles[i]=phi;
		    p.unsorted_flavours[i]=phi->get_flavour();
		}
		for(std::size_t i=0;i<N_out;++i)
		{
		    phi=model_wrapper<model_t>::get_particle_by_pdg_id(pdgv[N_in+i]);
		    if(phi==NULL)
		    {
			return processes.end();
		    }
		    p.FS_particles[i]=phi;
		    p.unsorted_flavours[N_in+i]=phi->get_flavour();
		}
		p.sorted_flavours=p.unsorted_flavours;
		std::sort(p.sorted_flavours.begin(),p.sorted_flavours.begin()+N_in);
		std::sort(p.sorted_flavours.begin()+N_in,p.sorted_flavours.end());
		p.unsorted_pdg_ids=pdgv;
		p.sorted_pdg_ids=pdgv;
		std::sort(p.sorted_pdg_ids.begin(),p.sorted_pdg_ids.begin()+N_in);
		std::sort(p.sorted_pdg_ids.begin()+N_in,p.sorted_pdg_ids.end());
		return processes.insert(proc_iter,p);
	    }

	    /* Process comparison operators: */

	    bool operator == (const process<model_t,N_in,N_out>& p) const
	    {
		return (sorted_flavours==p.sorted_flavours);
	    }
	    bool operator != (const process<model_t,N_in,N_out>& p) const
	    {
		return (sorted_flavours!=p.sorted_flavours);
	    }
	    bool operator < (const process<model_t,N_in,N_out>& p) const
	    {
		return (sorted_flavours<p.sorted_flavours);
	    }
	    bool operator <= (const process<model_t,N_in,N_out>& p) const
	    {
		return !(p<(*this));
	    }
	    bool operator > (const process<model_t,N_in,N_out>& p) const
	    {
		return (p<(*this));
	    }
	    bool operator >= (const process<model_t,N_in,N_out>& p) const
	    {
		return !(this->operator<(p));
	    }

	    /* Function adding a process tree to a list of trees: */ 

	    std::list<tree_type>& add_tree(std::list<tree_type>& trees)
	    {
		trees.push_back(tree_type(IS_particles,FS_particles));
		tree_it=trees.end();
		--tree_it;
		return trees;
	    }

	    /* Function inserting a process tree to a list of trees: */ 

	    tree_iterator insert_tree(std::list<tree_type>& trees,tree_iterator it)
	    {
		tree_it=trees.insert(it,tree_type(IS_particles,FS_particles));
		return tree_it;
	    }

	    /* Tree iterator access: */

	    tree_iterator get_tree() const
	    {
		return tree_it;
	    }

	    /* Function returning the vector of flavour integers corresponding
	     * to the process: */

	    const vector<size_type,N_in+N_out>& get_flavours() const
	    {
		return unsorted_flavours;
	    }

	    /* Function returning the vector of sorted flavour integers
	     * corresponding to the process: */

	    const vector<size_type,N_in+N_out>& get_sorted_flavours() const
	    {
		return sorted_flavours;
	    }

	    /* Function returning the i-th flavour in the vector returned by the
	     * method above: */

	    size_type get_flavour(std::size_t n) const
	    {
		return unsorted_flavours[n];
	    }

	    /* Function returning a vector of pdg code integers corresponding
	     * to the process: */

	    const vector<int,N_external>& get_pdg_ids() const
	    {
		return unsorted_pdg_ids;
	    }

	    /* Function returning a vector of sorted pdg code integers corresponding
	     * to the process: */

	    const vector<int,N_external>& get_sorted_pdg_ids() const
	    {
		return sorted_pdg_ids;
	    }

	    /* Function returning the i-th pdg id in the vector returned by the
	     * method above: */

	    int get_pdg_id(std::size_t n) const
	    {
		return unsorted_pdg_ids[n];
	    }

	    /* Displaying functions: */
	    
	    std::ostream& print_particles(std::ostream& os) const
	    {
		if(N_in>0)
		{
		    for(std::size_t n=0;n<N_in-1;++n)
		    {
			os<<IS_particles[n]->get_name()<<",";
		    }
		    os<<IS_particles[N_in-1]->get_name()<<" --> ";
		}
		if(N_out>0)
		{
		    for(std::size_t n=0;n<N_out-1;++n)
		    {
			os<<FS_particles[n]->get_name()<<",";
		    }
		    os<<FS_particles[N_out-1]->get_name();
		}
		return os;
	    }
	    std::ostream& print(std::ostream& os) const
	    {
		std::stringstream s1,s2,s3;
		print_particles(s1);
		s2<<unsorted_flavours;
		s3<<unsorted_pdg_ids;
		os<<std::setw(7*(N_in+N_out)+4)<<std::left<<s1.str()<<std::setw(4*(N_in+N_out)+5)<<std::left<<s2.str()<<std::setw(7*(N_in+N_out)+5)<<std::left<<s3.str();
		return os;
	    }
	
	private:

	    /* Process constructors: */

	    process(){}
	    process(const vector<const particle_type*,N_in>& in,const vector<const particle_type*,N_out>& out):IS_particles(in),FS_particles(out)
	    {
		for(std::size_t n=0;n<N_in;++n)
		{
		    unsorted_flavours[n]=in[n]->get_flavour();
		    unsorted_pdg_ids[n]=in[n]->get_pdg_id();
		}
		for(std::size_t n=0;n<N_out;++n)
		{
		    unsorted_flavours[N_in+n]=out[n]->get_flavour();
		    unsorted_pdg_ids[N_in+n]=out[n]->get_pdg_id();
		}
		sorted_flavours=unsorted_flavours;
		std::sort(sorted_flavours.begin(),sorted_flavours.begin()+N_in);
		std::sort(sorted_flavours.begin()+N_in,sorted_flavours.end());
		sorted_pdg_ids=unsorted_pdg_ids;
		std::sort(sorted_pdg_ids.begin(),sorted_pdg_ids.begin()+N_in);
		std::sort(sorted_pdg_ids.begin()+N_in,sorted_pdg_ids.end());
	    }

	    /* Initial-state particle pointers: */

	    vector<const particle_type*,N_in> IS_particles;

	    /* Final-state particle pointers: */

	    vector<const particle_type*,N_out> FS_particles;
	    
	    /* Sorted and unsorted flavour vectors: */
	    
	    vector<size_type,N_external>unsorted_flavours,sorted_flavours;
	    
	    /* Sorted and unsorted pdg-id vectors: */
	    
	    vector<int,N_external>unsorted_pdg_ids,sorted_pdg_ids;

	    /* Process tree iterator: */

	    tree_iterator tree_it;
    };
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t process<model_t,N_in,N_out>::N_incoming;
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t process<model_t,N_in,N_out>::N_outgoing;
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t process<model_t,N_in,N_out>::N_external;

    /* Process comparison by flavour configuration: */

    template<class model_t,std::size_t N_in,std::size_t N_out>bool flav_comp(const process<model_t,N_in,N_out>& pr1,const process<model_t,N_in,N_out>& pr2)
    {
	return pr1.get_sorted_flavours()<pr2.get_sorted_flavours();
    }
    
    template<class model_t,std::size_t N_in,std::size_t N_out>bool flav_vec_comp(const process<model_t,N_in,N_out>& pr1,const vector<typename process<model_t,N_in,N_out>::size_type,process<model_t,N_in,N_out>::N_external>& fv)
    {
	return pr1.get_sorted_flavours()<fv;
    }

    /* Process comparison by pdg id configuration: */

    template<class model_t,std::size_t N_in,std::size_t N_out>bool pdg_id_comp(const process<model_t,N_in,N_out>& pr1,const process<model_t,N_in,N_out>& pr2)
    {
	if(pr1.get_sorted_pdg_ids()==pr2.get_sorted_pdg_ids())
	{
	    return pr1.get_sorted_flavours()<pr2.get_sorted_flavours();
	}
	return pr1.get_sorted_pdg_ids()<pr2.get_sorted_pdg_ids();
    }
    
    template<class model_t,std::size_t N_in,std::size_t N_out>bool pdg_id_vec_comp(const process<model_t,N_in,N_out>& pr1,const vector<int,process<model_t,N_in,N_out>::N_external>& pv)
    {
	return pr1.get_sorted_pdg_ids()<pv;
    }

    /* Predicate function determined by whether the correspoding tree is emtpy:
     * */

    template<class model_t,std::size_t N_in,std::size_t N_out>bool has_empty_tree(const process<model_t,N_in,N_out>& pr)
    {
	return pr.get_tree()->is_empty();
    }
}

#endif /*CAMGEN_PROCESS_H_*/

