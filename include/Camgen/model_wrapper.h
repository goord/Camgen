//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_MOD_WRAPPER_H_
#define CAMGEN_MOD_WRAPPER_H_

#include <map>
#include <Camgen/tensor.h>
#include <Camgen/forward_decs.h>
#include <Camgen/name_comp.h>
#include <Camgen/flav_comp.h>
#include <Camgen/vertex_comp.h>
#include <Camgen/part_family.h>
#include <Camgen/has_leg.h>
#include <Camgen/fusion_cont.h>
#include <Camgen/type_holders.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Model wrapper class definition. The model wrapper contains static lists of    *
 * all particles, vertices, fusion rules and family definitions corresponding to *
 * the model class template parameter. It contains various searching, insertion  *
 * and ersing algorithms for particles, vertices and jets.                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class model_wrapper
    {
	public:

	    /* Value and size type definitions: */

	    typedef typename get_basic_types<model_t>::value_type value_type;
	    typedef typename get_basic_types<model_t>::r_value_type r_value_type;
	    typedef typename get_basic_types<model_t>::size_type size_type;
	    
	    /* Particle type definitions: */
	    
	    typedef typename get_basic_types<model_t>::particle_type particle_type;

	    /* Vertex type definitions: */

	    typedef typename get_basic_types<model_t>::vertex_type vertex_type;

	    /* Fusion map type definition: */

	    typedef std::multimap<std::vector<const particle_type*>,fusion_class<model_t>,flavourvec_comp<particle_type> > fusion_map;
	    
	    /* Fusion map iterator definitions: */
	    
	    typedef typename fusion_map::iterator fusion_iterator;
	    typedef typename fusion_map::const_iterator const_fusion_iterator;

	    friend class model<model_t>;

	    /* Query for particle named 'name': */

	    static bool contains_particle(const std::string& name)
	    {
		name_select<particle_type>pred(name);
		return (std::find_if(particle_content.begin(),particle_content.end(),pred) != particle_content.end());
	    }

	    /* Particle with flavour 'flav' finder: */

	    static const particle_type* get_particle(const size_type flav)
	    {
		if(flav < particle_content.size())
		{
		    return particle_content[flav];
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle with flavour "<<flav<<" not found in model--returning NULL"<<endlog;
		    return NULL;
		}
	    }

	    /* Particle with name 'name' finder: */

	    static const particle_type* get_particle(const std::string& name)
	    {
		name_select<particle_type>pred(name);
		typename std::vector<particle_type*>::const_iterator it=std::find_if(particle_content.begin(),particle_content.end(),pred);
		if(it != particle_content.end())
		{
		    return *it;
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle with name "<<name<<" not found in model--returning NULL"<<endlog;
		    return NULL;
		}
	    }

	    /* Particle with pdg id finder: */

	    static const particle_type* get_particle_by_pdg_id(const int id)
	    {
		typename std::vector<particle_type*>::const_iterator it=particle_content.begin();
		for(;it!=particle_content.end();++it)
		{
		    if((*it)->get_pdg_id()==id)
		    {
			break;
		    }
		}
		if(it != particle_content.end())
		{
		    return *it;
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle with id "<<id<<" not found in model--returning NULL"<<endlog;
		    return NULL;
		}
	    }

	    /* Function adding a family definition with name 'str', containing
	     * particles in string 'phi' (particle name separated by commas): */

	    static void construct_family(const std::string& str,const std::string& phi)
	    {
		/* Construct family definition: */

		particle_family<model_t>* fam=new particle_family<model_t>(str);
		
		/* Identify particle names in second argument, and add them to
		 * the family: */
		
		std::string psi;
		const char* comma=",";
		const particle_type* part=NULL;
		for(size_type i=0;i<phi.length();++i)
		{
		    if(phi[i]==*comma)
		    {
			part=get_particle(psi);
			if(part != NULL)
			{
			    fam->add_particle(part);
			    psi.clear();
			}
		    }
		    else
		    {
			psi.push_back(phi[i]);
		    }
		}
		part=get_particle(psi);
		if(part != NULL)
		{
		    fam->add_particle(part);
		    psi.clear();
		}

		/* If the family is nonempty and not contained in the family
		 * definitions, add it: */

		if(fam->size() != 0)
		{
		    name_select< particle_family<model_t> >pred(str);
		    typename std::vector<particle_family<model_t>* >::iterator it=std::find_if(part_families.begin(),part_families.end(),pred);
		    if(it==part_families.end())
		    {
			part_families.push_back(fam);
		    }
		    else
		    {
			*it=fam;
		    }
		    model<model_t>::logfile()<<"particle family definition "<<*fam<<" added to to model..."<<std::endl;
		}
	    }

	    /* Jet definition with name 'name' finder: */

	    static const particle_family<model_t>* get_family(const std::string& name)
	    {
		name_select< particle_family<model_t> >pred(name);
		typename std::vector<particle_family<model_t>* >::const_iterator it=std::find_if(part_families.begin(),part_families.end(),pred);
		if(it != part_families.end())
		{
		    return *it;
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle family with name "<<name<<" not found in model--returning NULL"<<endlog;
		    return NULL;
		}
	    }

	    /* Jet query with name 'name': */

	    static bool contains_family(const std::string& name)
	    {
		name_select< particle_family<model_t> >pred(name);
		return (std::find_if(part_families.begin(),part_families.end(),pred) != part_families.end());
	    }

	    /* Maximal vertex rank: */

	    static size_type maximal_vertex_rank()
	    {
		return max_vertex_rank;
	    }

	    /* Number of particle flavours (=maximal particle flavour): */

	    static size_type flavours()
	    {
		return particle_content.size();
	    }

	    /* Particle fusion finder function: */

	    static std::pair<fusion_iterator,fusion_iterator>find_fusion(const std::vector<const particle_type*>& parts)
	    {
		return fusion_rules.equal_range(parts);
	    }

	    /* Fusion iterators: */

	    static const_fusion_iterator fusions_begin()
	    {
		return fusion_rules.begin();
	    }
	    static const_fusion_iterator fusions_end()
	    {
		return fusion_rules.end();
	    }

	    /* Particle content printing function: */

	    static std::ostream& print_particles(std::ostream& os)
	    {
		for(size_type i=0;i<particle_content.size();++i)
		{
		    os<<*particle_content[i]<<std::endl;
		}
		return os;
	    }

	    /* Jet content printing function: */

	    static std::ostream& print_families(std::ostream& os)
	    {
		for(size_type i=0;i<part_families.size();++i)
		{
		    os<<*part_families[i]<<std::endl;
		}
		return os;
	    }

	    /* Vertex content printing function: */

	    static std::ostream& print_vertices(std::ostream& os)
	    {
		for(size_type i=0;i<vertex_content.size();++i)
		{
		    os<<*vertex_content[i]<<std::endl;
		}
		return os;
	    }

	    /* Omitting a particle from the model: */

	    static void erase_particle(const std::string& str)
	    {
		const particle_type* phi=get_particle(str);

		/* str is a valid particle name: */

		if(phi != NULL)
		{
		    /* Erase the particle from the particle content: */

		    particle_content.erase(std::find(particle_content.begin(),particle_content.end(),phi));
		    model<model_t>::logfile()<<"particle "<<phi->get_name()<<" erased from model..."<<std::endl;
		    
		    /* Erase every vertex that contains the particle: */
		    
		    has_leg<model_t>pred(phi);

		    for(typename std::vector<vertex_type*>::iterator it=vertex_content.begin();it!=vertex_content.end();++it)
		    {
			if(pred(*it))
			{
			    fusion_iterator f_it=fusion_rules.begin();
			    while(f_it!=fusion_rules.end())
			    {
				if(f_it->second.get_vertex()==*it)
				{
				    fusion_rules.erase(f_it++);
				}
				else
				{
				    ++f_it;
				}
			    }
			    model<model_t>::logfile()<<"vertex "<<(*it)->get_name()<<" erased from model..."<<std::endl;
			}
		    }

		    typename std::vector<vertex_type*>::iterator it=std::remove_if(vertex_content.begin(),vertex_content.end(),pred);
		    vertex_content.resize(it-vertex_content.begin());

		    /* Removing the particle from the family definitions: */

		    for(size_type i=0;i<part_families.size();++i)
		    {
			part_families[i]->remove_particle(phi);
		    }

		    const particle_type* phibar=phi->get_anti_particle();
		    
		    if(phibar != phi)
		    {
			/* Removing the anti-particle from the particle list: */
			
			particle_content.erase(std::find(particle_content.begin(),particle_content.end(),phibar));
			model<model_t>::logfile()<<"particle "<<phibar->get_name()<<" erased from model..."<<std::endl;
			
			has_leg<model_t>pred(phibar);

			for(typename std::vector<vertex_type*>::iterator it=vertex_content.begin();it!=vertex_content.end();++it)
			{
			    if(pred(*it))
			    {
				fusion_iterator f_it=fusion_rules.begin();
				while(f_it!=fusion_rules.end())
				{
				    if(f_it->second.get_vertex()==*it)
				    {
					fusion_rules.erase(f_it++);
				    }
				    else
				    {
					++f_it;
				    }
				}
				model<model_t>::logfile()<<"vertex "<<(*it)->get_name()<<" erased from model..."<<std::endl;
			    }
			}

			/* Removing all vertices containing the anti-particle:
			 * */

			typename std::vector<vertex_type*>::iterator it=std::remove_if(vertex_content.begin(),vertex_content.end(),pred);
			vertex_content.resize(it-vertex_content.begin());

			/* Removing the particle from the family definitions: */

			for(size_type i=0;i<part_families.size();++i)
			{
			    part_families[i]->remove_particle(phibar);
			}
		    }
		}
	    }

	    /* Decoupling a particle from the model: */

	    static void decouple_particle(const std::string& str)
	    {
		particle_type* phi=get_private_particle(str);

		/* str is a valid particle name: */

		if(phi != NULL)
		{
		    if(phi->is_coupled())
		    {
			/* Decouple the particle: */

			phi->decouple();
			model<model_t>::logfile()<<"particle "<<phi->get_name()<<" decoupled from model..."<<std::endl;
			if(phi->get_anti_particle()!=phi)
			{
			    model<model_t>::logfile()<<"particle "<<phi->get_anti_particle()->get_name()<<" decoupled from model..."<<std::endl;
			}

			/* Decouple every vertex that contains the particle or
			 * anti-particle: */

			has_leg<model_t>pred(phi);
			for(typename std::vector<vertex_type*>::iterator it=vertex_content.begin();it!=vertex_content.end();++it)
			{
			    if(pred(*it) and (*it)->is_coupled())
			    {
				(*it)->decouple();
				model<model_t>::logfile()<<"vertex "<<(*it)->get_name()<<" decoupled from model..."<<std::endl;
			    }
			}
			if(phi->get_anti_particle()!=phi)
			{
			    has_leg<model_t>pred2(phi->get_anti_particle());
			    for(typename std::vector<vertex_type*>::iterator it=vertex_content.begin();it!=vertex_content.end();++it)
			    {
				if(pred2(*it) and (*it)->is_coupled())
				{
				    (*it)->decouple();
				    model<model_t>::logfile()<<"vertex "<<(*it)->get_name()<<" decoupled from model..."<<std::endl;
				}
			    }
			}
		    }
		}
	    }

	    /* Coupling a particle to the model: */

	    static void couple_particle(const std::string& str)
	    {
		particle_type* phi=get_private_particle(str);

		/* str is a valid particle name: */

		if(phi != NULL)
		{
		    if(!phi->is_coupled())
		    {
			/* Couple the particle: */

			phi->couple();
			model<model_t>::logfile()<<"particle "<<phi->get_name()<<" coupled to model..."<<std::endl;
			if(phi->get_anti_particle()!=phi)
			{
			    model<model_t>::logfile()<<"particle "<<phi->get_anti_particle()->get_name()<<" coupled to model..."<<std::endl;
			}
		    }
		}
	    }

	    /* Omitting all 3-vertices from the list which have as legs
	     * str1,str2 and str3: */

	    static void erase_vertex(const std::string& str1,const std::string& str2,const std::string& str3)
	    {
		const particle_type* phi1=get_particle(str1);
		const particle_type* phi2=get_particle(str2);
		const particle_type* phi3=get_particle(str3);

		if(phi1 != NULL and phi2 != NULL and phi3 != NULL)
		{
		    /* Construct leg vector: */

		    std::vector<const particle_type* >legs(3);
		    legs[0]=phi1;
		    legs[1]=phi2;
		    legs[2]=phi3;
		    has_legs<model_t>pred(legs);
		    typename std::vector<vertex_type*>::iterator v_it=vertex_content.begin();
		    fusion_iterator f_it;
		    
		    /* Erase all vertices and fusion rules having those legs: */
		    
		    while(v_it != vertex_content.end())
		    {
			if(pred(*v_it))
			{
			    f_it=fusion_rules.begin();
			    fusion_contains_vertex<model_t>fpred(*v_it);
			    while(f_it != fusion_rules.end())
			    {
				if(fpred(*f_it))
				{
				    fusion_rules.erase(f_it++);
				}
				else
				{
				    ++f_it;
				}
			    }
			    std::string vertex_name = (*v_it)->get_name();
			    vertex_content.erase(v_it++);
			    model<model_t>::logfile()<<"vertex "<<vertex_name<<" erased from model..."<<std::endl;
			}
			else
			{
			    ++v_it;
			}
		    }
		}
	    }

	    /* Decoupling all 3-vertices from the list which have as legs
	     * str1,str2 and str3: */

	    static void decouple_vertex(const std::string& str1,const std::string& str2,const std::string& str3)
	    {
		const particle_type* phi1=get_particle(str1);
		const particle_type* phi2=get_particle(str2);
		const particle_type* phi3=get_particle(str3);

		if(phi1 != NULL and phi2 != NULL and phi3 != NULL)
		{
		    /* Construct leg vector: */

		    std::vector<const particle_type* >legs(3);
		    legs[0]=phi1;
		    legs[1]=phi2;
		    legs[2]=phi3;
		    has_legs<model_t>pred(legs);
		    
		    /* Decouple all vertices having those legs: */
		    
		    for(typename std::vector<vertex_type*>::iterator v_it=vertex_content.begin();v_it!=vertex_content.end();++v_it)
		    {
			if(pred(*v_it) and (*v_it)->is_coupled())
			{
			    (*v_it)->decouple();
			    model<model_t>::logfile()<<"vertex "<<(*v_it)->get_name()<<" decoupled from model..."<<std::endl;
			}
		    }
		}
	    }

	    /* Coupling all 3-vertices from the list which have as legs
	     * str1,str2 and str3: */

	    static void couple_vertex(const std::string& str1,const std::string& str2,const std::string& str3)
	    {
		const particle_type* phi1=get_particle(str1);
		const particle_type* phi2=get_particle(str2);
		const particle_type* phi3=get_particle(str3);

		if(phi1 != NULL and phi2 != NULL and phi3 != NULL)
		{
		    /* Construct leg vector: */

		    std::vector<const particle_type* >legs(3);
		    legs[0]=phi1;
		    legs[1]=phi2;
		    legs[2]=phi3;
		    has_legs<model_t>pred(legs);
		    
		    /* Decouple all vertices having those legs: */
		    
		    for(typename std::vector<vertex_type*>::iterator v_it=vertex_content.begin();v_it!=vertex_content.end();++v_it)
		    {
			if(pred(*v_it) and !(*v_it)->is_coupled())
			{
			    (*v_it)->couple();
			    model<model_t>::logfile()<<"vertex "<<(*v_it)->get_name()<<" coupled to model..."<<std::endl;
			}
		    }
		}
	    }

	    /* Omitting all 4-vertices from the list which have as legs str1,
	     * str2, str3 and str4: */

	    static void erase_vertex(const std::string& str1,const std::string& str2,const std::string& str3,const std::string& str4)
	    {
		const particle_type* phi1=get_particle(str1);
		const particle_type* phi2=get_particle(str2);
		const particle_type* phi3=get_particle(str3);
		const particle_type* phi4=get_particle(str4);

		if(phi1 != NULL and phi2 != NULL and phi3 != NULL and phi4 != NULL)
		{
		    /* Construct leg vector: */

		    std::vector<const particle_type* >legs(4);
		    legs[0]=phi1;
		    legs[1]=phi2;
		    legs[2]=phi3;
		    legs[3]=phi4;
		    has_legs<model_t>pred(legs);
		    typename std::vector<vertex_type*>::iterator v_it=vertex_content.begin();
		    fusion_iterator f_it;
		    
		    /* Erase all vertices and fusion rules having those legs: */
		    
		    while(v_it != vertex_content.end())
		    {
			if(pred(*v_it))
			{
			    f_it=fusion_rules.begin();
			    fusion_contains_vertex<model_t>fpred(*v_it);
			    while(f_it != fusion_rules.end())
			    {
				if(fpred(*f_it))
				{
				    fusion_rules.erase(f_it++);
				}
				else
				{
				    ++f_it;
				}
			    }
			    std::string vertex_name = (*v_it)->get_name();
			    vertex_content.erase(v_it++);
			    model<model_t>::logfile()<<"vertex "<<vertex_name<<" erased from model..."<<std::endl;
			}
			else
			{
			    ++v_it;
			}
		    }
		}
	    }

	    /* Decoupling all 4-vertices from the list which have as legs
	     * str1,str2,str3 and str4: */

	    static void decouple_vertex(const std::string& str1,const std::string& str2,const std::string& str3,const std::string& str4)
	    {
		const particle_type* phi1=get_particle(str1);
		const particle_type* phi2=get_particle(str2);
		const particle_type* phi3=get_particle(str3);
		const particle_type* phi4=get_particle(str4);

		if(phi1 != NULL and phi2 != NULL and phi3 != NULL and phi4 != NULL)
		{
		    /* Construct leg vector: */

		    std::vector<const particle_type* >legs(4);
		    legs[0]=phi1;
		    legs[1]=phi2;
		    legs[2]=phi3;
		    legs[3]=phi4;
		    has_legs<model_t>pred(legs);
		    
		    /* Decouple all vertices having those legs: */
		    
		    for(typename std::vector<vertex_type*>::iterator v_it=vertex_content.begin();v_it!=vertex_content.end();++v_it)
		    {
			if(pred(*v_it) and (*v_it)->is_coupled())
			{
			    (*v_it)->decouple();
			    model<model_t>::logfile()<<"vertex "<<(*v_it)->get_name()<<" decoupled from model..."<<std::endl;
			}
		    }
		}
	    }

	    /* Function decoupling all vertices with zero couplings and coupling
	     * vertices with nonzero couplings: */

	    static void decouple_zero_vertices()
	    {
		for(typename std::vector<vertex_type*>::iterator it=vertex_content.begin();it!=vertex_content.end();++it)
		{
		    (*it)->decouple_if_zero();
		}
	    }

	    /* Coupling all 4-vertices from the list which have as legs
	     * str1,str2,str3 and str4: */

	    static void couple_vertex(const std::string& str1,const std::string& str2,const std::string& str3,const std::string& str4)
	    {
		const particle_type* phi1=get_particle(str1);
		const particle_type* phi2=get_particle(str2);
		const particle_type* phi3=get_particle(str3);
		const particle_type* phi4=get_particle(str4);

		if(phi1 != NULL and phi2 != NULL and phi3 != NULL and phi4 != NULL)
		{
		    /* Construct leg vector: */

		    std::vector<const particle_type* >legs(4);
		    legs[0]=phi1;
		    legs[1]=phi2;
		    legs[2]=phi3;
		    legs[3]=phi4;
		    has_legs<model_t>pred(legs);
		    
		    /* Decouple all vertices having those legs: */
		    
		    for(typename std::vector<vertex_type*>::iterator v_it=vertex_content.begin();v_it!=vertex_content.end();++v_it)
		    {
			if(pred(*v_it) and !(*v_it)->is_coupled())
			{
			    (*v_it)->couple();
			    model<model_t>::logfile()<<"vertex "<<(*v_it)->get_name()<<" coupled to model..."<<std::endl;
			}
		    }
		}
	    }

	    /* Modify the propagator of a particle names 'str' to the type prop_t: */

	    template<class prop_t>static void set_propagator(const std::string& str)
	    {
		particle_type* phi=get_private_particle(str);
		if(phi!=NULL)
		{
		    phi->template set_propagator<prop_t>();
		}
	    }

	    /* Functions inserting a masses/widths after model construction: */

	    static void set_mass(const std::string& name,const r_value_type* m)
	    {
		particle_type* phi=get_private_particle(name);
		if(phi!=NULL)
		{
		    phi->set_mass(m);
		}
	    }
	    static void set_massless(const std::string& name)
	    {
		particle_type* phi=get_private_particle(name);
		if(phi!=NULL)
		{
		    phi->set_mass(NULL);
		}
	    }
	    static void set_width(const std::string& name,const r_value_type* w)
	    {
		particle_type* phi=get_private_particle(name);
		if(phi!=NULL)
		{
		    phi->set_width(w);
		}
	    }
	    static void set_widthless(const std::string& name)
	    {
		particle_type* phi=get_private_particle(name);
		if(phi!=NULL)
		{
		    phi->set_width(NULL);
		}
	    }

	    /* Fusion rule printing function: */

	    static std::ostream& print_fusion_rules(std::ostream& os)
	    {
		for(fusion_iterator it=fusion_rules.begin();it != fusion_rules.end();++it)
		{
		    for(size_type i=0;i<it->first.size()-1;++i)
		    {
			os<<it->first[i]->get_name()<<",";
		    }
		    os<<it->first.back()->get_name()<<"--->"<<it->second.get_produced_particle()->get_name()<<std::endl;
		}
		return os;
	    }

	protected:

	    /* Particle inserter function: */

	    static void insert_particle(particle_type* phi)
	    {
		if(phi != NULL)
		{
		    name_select<particle_type>pred(phi);
		    typename std::vector<particle_type*>::iterator it=std::find_if(particle_content.begin(),particle_content.end(),pred);
		    if(it==particle_content.end())
		    {
			particle_content.push_back(phi);
		    }
		    else
		    {
			phi->flavour=(*it)->flavour;
			particle_type::flavour_nr--;
			*it=phi;
		    }
		    model<model_t>::logfile()<<"particle "<<phi->get_name()<<" added to model..."<<std::endl;
		}
	    }

	    /* Particle/antiparticle pair inserter function: */

	    static void insert_particle_pair(const std::pair<particle_type*,particle_type*>& phi)
	    {
		if(phi.first != NULL and phi.second != NULL)
		{
		    name_select<particle_type>pred1(phi.first);
		    typename std::vector<particle_type*>::iterator it=std::find_if(particle_content.begin(),particle_content.end(),pred1);
		    if(it==particle_content.end())
		    {
			particle_content.push_back(phi.first);
		    }
		    else
		    {
			phi.first->flavour=(*it)->flavour;
			particle_type::flavour_nr--;
			*it=phi.first;
		    }
		    name_select<particle_type>pred2(phi.second);
		    it=std::find_if(particle_content.begin(),particle_content.end(),pred2);
		    if(it==particle_content.end())
		    {
			particle_content.push_back(phi.second);
		    }
		    else
		    {
			phi.second->flavour=(*it)->flavour;
			particle_type::flavour_nr--;
			*it=phi.second;
		    }
		    model<model_t>::logfile()<<"particle/anti-particle pair "<<phi.first->get_name()<<','<<phi.second->get_name()<<" added to model..."<<std::endl;
		}
	    }

	    /* Vertex inserter function: */

	    static void insert_vertex(vertex_type* V)
	    {
		if(V != NULL)
		{
		    vertex_content.push_back(V);
		    vertex_comp<model_t>comp;
		    std::sort(vertex_content.begin(),vertex_content.end(),comp);
		    if(V->get_rank()>max_vertex_rank)
		    {
			max_vertex_rank=V->get_rank();
		    }
		    fusion_class<model_t>::insert_fusion_rules(V,fusion_rules);
		    model<model_t>::logfile()<<"vertex "<<V->get_name()<<" added to model..."<<std::endl;
		}
	    }

	private:
	    
	    /* Particle content of the model: */

	    static std::vector<particle_type*>particle_content;
	    
	    /* Vertex content of the model: */
	    
	    static std::vector<vertex_type*>vertex_content;
	    
	    /* Maximal vertex rank in the model: */
	    
	    static size_type max_vertex_rank;

	    /* List of particle family definitions in the model: */

	    static std::vector< particle_family<model_t>* > part_families;
	    
	    /* Map of fusion rules: */
	    
	    static fusion_map fusion_rules;

	    /* Non-const particle pointer finding function: */

	    static particle_type* get_private_particle(const std::string& name)
	    {
		name_select<particle_type>pred(name);
		typename std::vector<particle_type*>::const_iterator it=std::find_if(particle_content.begin(),particle_content.end(),pred);
		if(it != particle_content.end())
		{
		    return *it;
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle with name "<<name<<" not found in model--returning NULL"<<endlog;
		    return NULL;
		}
	    }
    };

    template<class model_t>std::vector<typename model_wrapper<model_t>::particle_type*> model_wrapper<model_t>::particle_content;
    template<class model_t>std::vector<typename model_wrapper<model_t>::vertex_type*> model_wrapper<model_t>::vertex_content;
    template<class model_t>typename model_wrapper<model_t>::size_type model_wrapper<model_t>::max_vertex_rank=0;
    template<class model_t>std::vector<particle_family<model_t>*> model_wrapper<model_t>::part_families;
    template<class model_t>typename model_wrapper<model_t>::fusion_map model_wrapper<model_t>::fusion_rules;
}

#endif /*MOD_WRAPPER_H_*/

