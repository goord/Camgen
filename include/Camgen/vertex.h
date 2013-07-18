//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VERTEX_H_
#define CAMGEN_VERTEX_H_

#include <bitset>
#include <Camgen/flav_comp.h>
#include <Camgen/charge_conj.h>
#include <Camgen/particle.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Vertex class definition. The vertex class contains information such as the    *
 * legs (vector of particle pointers), couplings (vector of complex number       *
 * pointers) and the Feynman rules (4-array of function pointers). For the       *
 * latter there are 3 more arrays for the case, containing the charge conjugated *
 * and reversed Feynman rules, which may have to be called in precesses          *
 * involving Majorana fermions. A vertex object dispatches the correct Feynman   *
 * rule function pointer to an interaction, which performs the actual callback   *
 * of the recursive relation.                                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition of the vertex class: */

    template<class model_t>class vertex
    {
	public:

	    /* The usual type defintions: */

	    DEFINE_BASIC_TYPES(model_t);
	    
	    /* Type definitions involving particle pointer vectors: */
	    
	    typedef particle<model_t> particle_type;
	    typedef typename std::vector<const particle_type*>::iterator particle_iterator;
	    typedef typename std::vector<const particle_type*>::const_iterator const_particle_iterator;

	    /* Type definition of a function pointer returning a string: */

	    typedef std::string (*formula_getter)();

	    /* Recursive relation function pointer type definition: */
	    
	    typedef typename get_basic_types<model_t>::vert_func vert_func;

	    /* Friend declarations: */

	    friend class has_leg<model_t>;
	    friend class has_legs<model_t>;

	    /* Static 3-vertex pointer creation method: */

	    template<class Feynrule_t>static vertex<model_t>* create(const particle_type* phi1,const particle_type* phi2,const particle_type* phi3,const value_type* c0,const value_type* c1=NULL,const value_type* c2=NULL)
	    {
		/* Check if all particles are valid: */
		
		if(phi1==NULL or phi2==NULL or phi3==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"NULL particle pointer argument in vertex detected--no vertex constructed"<<endlog;
		    return NULL;
		}

		/* Check if the rank of the Feynman rule matches: */

		if(Feynrule_t::rank!=3)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3);
		    log<<std::endl<<"vertex rank doesn't match number of legs--no vertex constructed"<<endlog;
		    return NULL;
		}

		/* Check if the number of given couplings matches the Feynman
		 * rule: */

		switch(Feynrule_t::params)
		{
		    case 1:
			if(c0==NULL)
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC;
			    alert(log,phi1,phi2,phi3);
			    log<<std::endl<<"incorrect number of couplings given--no vertex constructed"<<endlog;
			    return NULL;
			}
			break;
		    case 2:
			if(c0==NULL or c1==NULL)
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC;
			    alert(log,phi1,phi2,phi3);
			    log<<std::endl<<"incorrect number of couplings given--no vertex constructed"<<endlog;
			    return NULL;
			}
			break;
		    case 3:
			if(c0==NULL or c1==NULL or c2==NULL)
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC;
			    alert(log,phi1,phi2,phi3);
			    log<<std::endl<<"incorrect number of couplings given--no vertex constructed"<<endlog;
			    return NULL;
			}
			break;
		}

		/* Check if the index ranges of the participating particle
		 * subamplitudes match the feynman rule: */

		std::vector<size_type>r;
		if(evaluate<Feynrule_t>::get_index_ranges(0,r) != phi1->index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3);
		    log<<std::endl<<"particle "<<phi1->name<<" has incorrect tensor rank--no vertex constructed"<<endlog;
		    return NULL;
		}
		if(evaluate<Feynrule_t>::get_index_ranges(1,r) != phi2->index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3);
		    log<<std::endl<<"particle "<<phi2->name<<" has incorrect tensor rank--no vertex constructed"<<endlog;
		    return NULL;
		}
		if(evaluate<Feynrule_t>::get_index_ranges(2,r) != phi3->index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3);
		    log<<std::endl<<"particle "<<phi3->name<<" has incorrect tensor rank--no vertex constructed"<<endlog;
		    return NULL;
		}

		/* Proceed with vertex construction: */

		vertex<model_t>* vert=new vertex<model_t>(phi1,phi2,phi3);
		
		/* Insert non-nullary couplings: */
		
		vert->add_coupling(c0);
		if(Feynrule_t::params>1)
		{
		    vert->add_coupling(c1);
		}
		if(Feynrule_t::params>2)
		{
		    vert->add_coupling(c2);
		}

		/* Insert Feynman rule: */

		vert->insert_Feynman_rule<Feynrule_t>();
		
		/* Return the address: */

		return vert;
	    }

	    /* Static 4-vertex pointer creation method: */

	    template<class Feynrule_t>static vertex<model_t>* create(const particle_type* phi1,const particle_type* phi2,const particle_type* phi3,const particle_type* phi4,const value_type* c0,const value_type* c1=NULL,const value_type* c2=NULL)
	    {
		/* Check if all particles are valid: */
		
		if(phi1==NULL or phi2==NULL or phi3==NULL or phi4==NULL)
		{
		    log(log_level::warning)<<"NULL particle pointer argument in vertex detected--no vertex constructed"<<endlog;
		    return NULL;
		}

		/* Check if the rank of the Feynman rule matches: */

		if(Feynrule_t::rank!=4)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3,phi4);
		    log<<std::endl<<"vertex rank doesn't match number of legs--no vertex constructed"<<endlog;
		    return NULL;
		}

		/* Check if the number of given couplings matches the Feynman
		 * rule: */

		switch(Feynrule_t::params)
		{
		    case 1:
			if(c0==NULL)
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC;
			    alert(log,phi1,phi2,phi3,phi4);
			    log<<std::endl<<"incorrect number of couplings given--no vertex constructed"<<endlog;
			    return NULL;
			}
			break;
		    case 2:
			if(c0==NULL or c1==NULL)
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC;
			    alert(log,phi1,phi2,phi3,phi4);
			    log<<std::endl<<"incorrect number of couplings given--no vertex constructed"<<endlog;
			    return NULL;
			}
			break;
		    case 3:
			if(c0==NULL or c1==NULL or c2==NULL)
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC;
			    alert(log,phi1,phi2,phi3,phi4);
			    log<<std::endl<<"incorrect number of couplings given--no vertex constructed"<<endlog;
			    return NULL;
			}
			break;
		}

		/* Check if the index ranges of the participating particle
		 * subamplitudes match the feynman rule: */

		std::vector<size_type>r;
		if(evaluate<Feynrule_t>::get_index_ranges(0,r) != phi1->index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3,phi4);
		    log<<std::endl<<"particle "<<phi1->name<<" has incorrect tensor rank--no vertex constructed"<<endlog;
		    return NULL;
		}
		if(evaluate<Feynrule_t>::get_index_ranges(1,r) != phi2->index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3,phi4);
		    log<<std::endl<<"particle "<<phi2->name<<" has incorrect tensor rank--no vertex constructed"<<endlog;
		    return NULL;
		}
		if(evaluate<Feynrule_t>::get_index_ranges(2,r) != phi3->index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3,phi4);
		    log<<std::endl<<"particle "<<phi3->name<<" has incorrect tensor rank--no vertex constructed"<<endlog;
		    return NULL;
		}
		if(evaluate<Feynrule_t>::get_index_ranges(3,r) != phi4->index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC;
		    alert(log,phi1,phi2,phi3,phi4);
		    log<<std::endl<<"particle "<<phi4->name<<" has incorrect tensor rank--no vertex constructed"<<endlog;
		    return NULL;
		}

		/* Proceed with vertex construction: */

		vertex<model_t>* vert=new vertex<model_t>(phi1,phi2,phi3,phi4);
		
		/* Insert non-nullary couplings: */
		
		vert->add_coupling(c0);
		if(Feynrule_t::params>1)
		{
		    vert->add_coupling(c1);
		}
		if(Feynrule_t::params>2)
		{
		    vert->add_coupling(c2);
		}

		/* Insert Feynman rule: */

		vert->insert_Feynman_rule<Feynrule_t>();

		/* Return the address: */

		return vert;
	    }

	    /* Vertex rank readout: */

	    size_type get_rank() const
	    {
		return legs.size();
	    }

	    /* n-th leg readout: */

	    const particle_type* get_leg(size_type n) const
	    {
		CAMGEN_ERROR_IF((n>=legs.size()),"leg request out of range");
		return legs[n];
	    }

	    /* n-th sorted leg readout: */

	    const particle_type* get_sorted_leg(size_type n) const
	    {
		CAMGEN_ERROR_IF((n>=sorted_legs.size()),"sorted leg request out of range");
		return sorted_legs[n];
	    }

	    /* Leg vector readout: */

	    const std::vector<const particle_type*>& get_legs() const
	    {
		return legs;
	    }

	    /* Sorted leg vector readout: */

	    const std::vector<const particle_type*>& get_sorted_legs() const
	    {
		return sorted_legs;
	    }

	    /* Vertex name construction: creates a string with all the
	     * interaction particle names appended and appends 'vertex': */

	    std::string get_name() const
	    {
		std::string str;
		for(size_type i=0;i<legs.size();++i)
		{
		    str.append(legs[i]->name);
		}
		str.append("vertex");
		return str;
	    }

	    /* Feynman rule formula readout: */

	    std::string get_Feynman_rule() const
	    {
		return formula;
	    }

	    /* Number of coupling constants used by the vertex: */

	    size_type nr_of_couplings() const
	    {
		return couplings.size();
	    }
	    /* n-th coupling constant readout: */

	    const value_type& get_coupling(size_type n) const
	    {
		CAMGEN_ERROR_IF((n>=couplings.size()),"coupling constant request out of range");
		return *(couplings[n]);
	    }

	    const std::vector<const value_type*>& get_couplings() const
	    {
		return couplings;
	    }

	    /* Majorana-type readout: */
	    
	    size_type get_Majorana_type() const
	    {
		return Majorana_type;
	    }

	    /* Fermionic property readout: */

	    bool is_fermionic() const
	    {
		return fermionic;
	    }

	    /* Binary search whether a particle matches an incoming leg: */

	    bool is_incoming(const particle_type* phi) const
	    {
		return std::binary_search(sorted_legs.begin(),sorted_legs.end(),phi,comp);
	    }

	    /* Binary search whether a particle matches an outgoing leg: */

	    bool is_outgoing(const particle_type* phi) const
	    {
		return is_incoming(phi->get_anti_particle());
	    }

	    /* Comparison operators: determined by the lexicographical
	     * comparison of the legs and the Feynman rules: */

	    bool operator == (const vertex<model_t>& other) const
	    {
		return (legs==other.legs and Feynman_rule[0]==other.Feynman_rule[0]);
	    }
	    bool operator != (const vertex<model_t>& other) const
	    {
		return (legs!=other.legs or Feynman_rule[0]!=other.Feynman_rule[0]);
	    }
	    bool operator < (const vertex<model_t>& other) const
	    {
		if(legs==other.legs)
		{
		    return Feynman_rule[0]<other.Feynman_rule[0];
		}
		return std::lexicographical_compare(legs.begin(),legs.end(),other.legs.begin(),other.legs.end(),comp);
	    }
	    bool operator <= (const vertex<model_t>& other) const
	    {
		return !(other<(*this));
	    }
	    bool operator >= (const vertex<model_t>& other) const
	    {
		return !(*this<other);
	    }
	    bool operator > (const vertex<model_t>& other) const
	    {
		return (other<(*this));
	    }

	    /* Const iterators to the leg vector: */

	    const_particle_iterator begin_legs() const
	    {
		return legs.begin();
	    }
	    const_particle_iterator end_legs() const
	    {
		return legs.end();
	    }

	    /* Const iterator to the sorted leg vector: */

	    const_particle_iterator begin_sorted_legs() const
	    {
		return sorted_legs.begin();
	    }
	    const_particle_iterator end_sorted_legs() const
	    {
		return sorted_legs.end();
	    }

	    /* Recursive relation dispatching function: */
	    
	    vert_func dispatch_Feynman_rule(std::bitset<4>flow,size_type prod_part,bool& swap_fermions) const
	    {
	 	CAMGEN_ERROR_IF((prod_part>=legs.size()),"produced particle index out of range");

		switch(this->Majorana_type)
		{
		    case 0:

			/* Ordinary vertices: no current reversal */
			
			return Feynman_rule[prod_part];

		    case 1:

			/* anti-Dirac-Majorana vertex type: */
			
			switch(prod_part)
			{
			    /* Boson production: */

			    case 0:
				
				/* If the Majorana particle is outgoing, flip
				 * its fermion flow: */

				if(flow[2])
				{
				    return Feynman_rule_C[0];
				}

				/* Else, no current reversal: */

				else
				{
				    return Feynman_rule[0];
				}
				
			    /* Dirac fermion production: */

			    case 1:

				/* If the Majorana particle is outgoing, flip
				 * its fermion flow: */

				if(flow[2])
				{
				    return Feynman_rule_C[1];
				}

				/* Else, no current reversal: */

				else
				{
				    return Feynman_rule[1];
				}

			    /* Majorana fermion production: */
				
			    case 2:

				/* The produced Majorana fermion should carry
				 * the fermion flow along its momentum.
				 * Consequently the incoming anti-fermion or
				 * outgoing fermion needs to be flipped. Because
				 * a Dirac current is reversed, the reversed
				 * Feynman rule is dispatched */

				return reversed_Feynman_rule_C[2];
			}

		    case 2:

			/* Majorana-Dirac vertex type: */
			
			switch(prod_part)
			{
			    /* Boson production: */

			    case 0:

				/* If the Majorana particle is outgoing, the
				 * vertex behaves as an antifermion-fermion
				 * vertex: */

				if(flow[1])
				{
				    return Feynman_rule[0];
				}

				/* Else, the Majorana fermion is reversed: */

				else
				{
				    return Cc_Feynman_rule[0];
				}

			    /* Majorana fermion production: */

			    case 1:

				return Feynman_rule[1];

			    /* Dirac antifermion production: */

			    case 2:

				/* If the Majorana fermion is outgoing, it
				 * behaves as an antifermion and no reversal has
				 * to be performed: */

				if(flow[1])
				{
				    return Feynman_rule[2];
				}

				/* Otherwise, the Majorana fermion spinor has to
				 * be reversed: */

				else
				{
				    return Cc_Feynman_rule[2];
				}
			}

		    case 3:

			/* Majorana-Majorana vertex type: */

			switch(prod_part)
			{
			    /* Boson production: */

			    case 0:

				/* If both Majorana particles are incoming,
				 * reverse the first fermion flow: */

				if(!flow[1] and !flow[2])
				{
				    return Cc_Feynman_rule[0];
				}

				/* If the first particle is outgoing, no
				 * reversal has to be performed: */

				else if(flow[1] and !flow[2])
				{
				    return Feynman_rule[0];
				}

				/* If the second particle is outgoing, we
				 * swap the currents */

				else if(!flow[1] and flow[2])
				{
				    swap_fermions=true;
				    return Feynman_rule[0];
				}

				/* If the particles are outgoing, reverse the
				 * second fermion: */

				else
				{
				    return Feynman_rule_C[0];
				}

			    /* Majorana particle production: */

			    case 1:

				/* If the second fermion is outgoing, reverse
				 * it: */

				if(flow[2])
				{
				    return Feynman_rule_C[1];
				}

				/* Else, no reversal is performed: */

				else
				{
				    return Feynman_rule[1];
				}

			    case 2:

				/* The final particle is constructed by swapping
				 * the fermions and applying the second step: */

				swap_fermions=true;
				
				if(flow[1])
				{
				    return Feynman_rule_C[1];
				}
				else
				{
				    return Feynman_rule[1];
				}
			}
		}
		return NULL;
	    }

	    /* Function coupling the vertex to the model it resides in: */

	    void couple()
	    {
		coupled=true;
	    }

	    /* Function decoupling the vertex from the model it resides in: */

	    void decouple()
	    {
		coupled=false;
	    }

	    /* Function decoupling the vertex from the model if the couplings
	     * constants are zero: */

	    void decouple_if_zero()
	    {
		value_type zero(0,0);
		coupled=false;
		for(size_type i=0;i<couplings.size();++i)
		{
		    coupled|=(*(couplings[i])!=zero);
		}
	    }

	    /* Request whether the vertex participates in tree evaluations: */

	    bool is_coupled() const
	    {
		return coupled;
	    }
	
	private:

	    /* Leg vector. Sequence of particles (regarded as incoming)
	     * participating in the vertex. Maximally of length 4: */

	    std::vector<const particle_type*> legs;

	    /* Sorted leg vector. Leg vector ordered with respect to the
	     * particle flavour comparison: */

	    std::vector<const particle_type*> sorted_legs;
	    
	    /* Vector of couplings, implemented as pointers to complex numbers:
	     * */
	    
	    std::vector<const value_type*> couplings;

	    /* Static particle ponter comparison object: */

	    static flavour_comp< particle_type > comp;
	    
	    /* Majorana-type tag. 0 corresponds to a vertex without Majorana
	     * particles, 1 to a vertex with the second particle Majorana, 2 if
	     * the third particle is Majorana-like, and 3 for
	     * boson-Majorana-Majorana fermion type vertices: */
	    
	    size_type Majorana_type;

	    /* Fermion property tag. Is true if the vertex is
	     * boson-fermion-fermion type: */

	    bool fermionic;
	    
	    /* Function pointer returning the Feynman rule formula: */
	    
	    std::string formula;
	    
	    /* Recursive relation function pointers: */
	    
	    vert_func Feynman_rule[4];

	    /* Right charge-conjugate recursive relation function pointers: */

	    vert_func Feynman_rule_C[4];

	    /* Left charge-conjugate recursive relation function pointers: */

	    vert_func Cc_Feynman_rule[4];
	    
	    /* Right cgarge-conjugate reversed recursive relation function
	     * pointers: */
	    
	    vert_func reversed_Feynman_rule_C[4];

	    /* Flag denoting whether the vertex participates in tree
	     * evaluations: */

	    bool coupled;

	    /* Protected trivial constructor: */

	    vertex():Majorana_type(0),coupled(true)
	    {
		/* Initialise recursive relations to NULL: */

		for(size_type i=0;i<4;++i)
		{
		    Feynman_rule[i]=NULL;
		    Feynman_rule_C[i]=NULL;
		    Cc_Feynman_rule[i]=NULL;
		    reversed_Feynman_rule_C[i]=NULL;
		}
	    }

	    /* Constructor from 3 legs: */

	    vertex(const particle_type* phi1,const particle_type* phi2,const particle_type* phi3):Majorana_type(0),coupled(true)
	    {
		legs.push_back(phi1);
		legs.push_back(phi2);
		if(phi2->is_Majorana())
		{
		    Majorana_type=2;
		}
		legs.push_back(phi3);
		if(phi3->is_Majorana())
		{
		    Majorana_type+=1;
		}
		sorted_legs=legs;
		std::sort(sorted_legs.begin(),sorted_legs.end(),comp);
		
		/* Initialise recursive relations to NULL: */

		for(size_type i=0;i<4;++i)
		{
		    Feynman_rule[i]=NULL;
		    Feynman_rule_C[i]=NULL;
		    Cc_Feynman_rule[i]=NULL;
		    reversed_Feynman_rule_C[i]=NULL;
		}
	    }

	    /* Constructor from 4 legs: */

	    vertex(const particle_type* phi1,const particle_type* phi2,const particle_type* phi3,const particle_type* phi4):Majorana_type(0),coupled(true)
	    {
		legs.push_back(phi1);
		legs.push_back(phi2);
		legs.push_back(phi3);
		legs.push_back(phi4);

		sorted_legs=legs;
		std::sort(sorted_legs.begin(),sorted_legs.end(),comp);
		
		/* Initialise recursive relations to NULL: */

		for(size_type i=0;i<4;++i)
		{
		    Feynman_rule[i]=NULL;
		    Feynman_rule_C[i]=NULL;
		    Cc_Feynman_rule[i]=NULL;
		    reversed_Feynman_rule_C[i]=NULL;
		}
	    }

	    /* Trivial copy constructor: */

	    vertex(const vertex<model_t>& other){}
	   
	    /* Feynman rule insertion method: */

	    template<class Feynrule_t>void insert_Feynman_rule()
	    {
		fermionic=Feynrule_t::fermionic;
		formula=Feynrule_t::formula;
		typedef get_recursive_relation<model_t,get_colour_treatment<model_t,model_t::coloured>::decomposes> rr_getter;
		typedef typename rr_getter::template apply<Feynrule_t>::type eval_type;
		typedef typename rr_getter::template apply< charge_conj_right<Feynrule_t,Feynrule_t::fermionic> >::type eval_C_type;
		typedef typename rr_getter::template apply< charge_conj_left<Feynrule_t,Feynrule_t::fermionic> >::type Cc_eval_type;
		typedef typename rr_getter::template apply< reverse_conj_right<Feynrule_t,Feynrule_t::fermionic> >::type reverse_eval_C_type;
		
		eval_type::initialise();

		Feynman_rule[0]=&(eval_type::first);
		Feynman_rule[1]=&(eval_type::second);
		Feynman_rule[2]=&(eval_type::third);
		Feynman_rule[3]=&(eval_type::fourth);

		if(Majorana_type != 0)
		{
		    evaluate< charge_conj_right<Feynrule_t,Feynrule_t::fermionic> >::initialise();
		    Feynman_rule_C[0]=&(eval_C_type::first);
		    Feynman_rule_C[1]=&(eval_C_type::second);
		    Feynman_rule_C[2]=&(eval_C_type::third);
		    Feynman_rule_C[3]=&(eval_C_type::fourth);

		    evaluate< charge_conj_left<Feynrule_t,Feynrule_t::fermionic> >::initialise();
		    Cc_Feynman_rule[0]=&(Cc_eval_type::first);
		    Cc_Feynman_rule[1]=&(Cc_eval_type::second);
		    Cc_Feynman_rule[2]=&(Cc_eval_type::third);
		    Cc_Feynman_rule[3]=&(Cc_eval_type::fourth);
		}
		if(Majorana_type == 1)
		{
		    evaluate< reverse_conj_right<Feynrule_t,Feynrule_t::fermionic> >::initialise();
		    reversed_Feynman_rule_C[0]=&(reverse_eval_C_type::first);
		    reversed_Feynman_rule_C[1]=&(reverse_eval_C_type::second);
		    reversed_Feynman_rule_C[2]=&(reverse_eval_C_type::third);
		    reversed_Feynman_rule_C[3]=&(reverse_eval_C_type::fourth);
		}
	    }
	    /* Coupling constant insertion method: */
	    
	    void add_coupling(const value_type* c)
	    {
		couplings.push_back(c);
	    }

	    /* Logging functions: */

	    static logstream& alert(logstream& os,const particle_type* phi1,const particle_type* phi2,const particle_type* phi3)
	    {
		return os<<"failed to make "<<phi1->name<<phi2->name<<phi3->name<<"-vertex: ";
	    }
	    static logstream& alert(logstream& os,const particle_type* phi1,const particle_type* phi2,const particle_type* phi3,const particle_type* phi4)
	    {
		return os<<"failed to make "<<phi1->name<<phi2->name<<phi3->name<<phi4->name<<"-vertex: ";
	    }
    };

    /* Output operator: */

    template<class model_t>std::ostream& operator << (std::ostream& os,const vertex<model_t>& V)
    {
	os<<std::setw(40)<<std::left<<V.get_name();
	std::stringstream s;
	std::size_t n=V.nr_of_couplings()-1;
	for(std::size_t i=0;i<n;++i)
	{
	    s<<V.get_coupling(i)<<",";
	}
	s<<V.get_coupling(n);
	os<<std::setw(60)<<std::left<<s.str();
	if(V.is_fermionic())
	{
	    switch(V.get_Majorana_type())
	    {
		case 0:
		    os<<std::setw(8)<<std::left<<"(aDD)";
		    break;
		case 1:
		    os<<std::setw(8)<<std::left<<"(aDM)";
		    break;
		case 2:
		    os<<std::setw(8)<<std::left<<"(MD)";
		    break;
		case 3:
		    os<<std::setw(8)<<std::left<<"(MM)";
		    break;
	    }
	}
	else
	{
	    os<<std::setw(8)<<std::left<<"(B)";
	}
	os<<std::setw(8)<<std::left<<V.is_coupled();
	os<<std::endl<<"Feynman rule: "<<V.get_Feynman_rule();
	return os;
    }
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VERTEX_H_*/

