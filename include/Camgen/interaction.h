//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_INTERACTION_H_
#define CAMGEN_INTERACTION_H_

#include <Camgen/utils.h>
#include <Camgen/curr_iter_comp.h>
#include <Camgen/vertex.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the interaction classes in Camgen. There are two               *
 * implementations, depending on whether we are in colour decomposition mode or  *
 * not, and both classes are derived from a common base class. This              *
 * interaction_base contains the interaction vertex address, the vector of       *
 * interacting current addresses, an integer that denotes which of the currents  *
 * is produced, and some utility information on the momentum flows and           *
 * computation policy. The derived classes only differ in the callback of the    *
 * Feynman rule: the colour-decomposed interaction type contains an additional   *
 * vector which runs over all combinations of propagating colour singlets in the *
 * incoming current subamplitude tensors.                                        *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,std::size_t N>class interaction_base
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	    /* Particle type definition: */

	    typedef typename get_basic_types<model_t>::particle_type particle_type;
	    
	    /* Vertex type definition: */
	    
	    typedef typename get_basic_types<model_t>::vertex_type vertex_type;

	    /* Feynman rule definition: */

	    typedef typename get_basic_types<model_t>::vert_func vert_func;

	private:

	    /* Boolean denoting whether we are in colour decomposition mode: */

	    static const bool decomposes=get_colour_treatment<model_t,model_t::coloured>::decomposes;

	public:

	    /* Current type definition: */

	    typedef current<model_t,N,decomposes> current_type;

	    /* Current vector iterator definitions: */

	    typedef typename std::vector<current_type>::iterator current_iterator;

	    typedef typename std::vector<current_type>::const_iterator const_current_iterator;

	    /* Friend declaration of the process_tree class template: */

	    template<class mod,std::size_t N_in,std::size_t N_out>friend class process_tree;

	    /* Trivial constructor: */

	    interaction_base():vertex_t(NULL),produced_current(0),Feynman_rule(NULL),swap_fermions(false),Fermi_sign(1),flow(0),CM_tag(false),coupled(vertex_t->is_coupled()),produced_momentum(NULL)
	    {
		++object_counter;
	    }

	    /* Regular constructor: */

	    interaction_base(const vertex_type* v,size_type n=0):vertex_t(v),currents(v->get_rank()),amp_iters(v->get_rank()),momenta(v->get_rank(),NULL),produced_current(n),Feynman_rule(NULL),swap_fermions(false),Fermi_sign(1),flow(0),CM_tag(false),coupled(vertex_t->is_coupled()),produced_momentum(NULL)
	    {
		CAMGEN_ERROR_IF((n>=v->get_rank()),"produced particle number out of range");
		++object_counter;
	    }

	    /* Copy constructor: */

	    interaction_base(const interaction_base<model_t,N>& other):vertex_t(other.vertex_t),currents(other.currents),amp_iters(other.amp_iters),momenta(other.momenta),produced_current(other.produced_current),Feynman_rule(other.Feynman_rule),swap_fermions(false),Fermi_sign(other.Fermi_sign),flow(other.flow),CM_tag(other.CM_tag),coupled(other.coupled),prop_policy(other.prop_policy),produced_momentum(other.produced_momentum)
	    {
		++object_counter;
	    }

	    /* Clearance member function, resetting all data members to their
	     * default values: */

	    void clear()
	    {
		momenta.clear();
		Fermi_sign=1;
		vertex_t=NULL;
		produced_current=0;
		flow.reset();
		CM_tag=false;
		prop_policy.reset();
		produced_momentum=NULL;
		currents.clear();
		amp_iters.clear();
		Feynman_rule=NULL;
		swap_fermions=false;
	    }

	    /* Vertex rank readout: */

	    size_type get_rank() const
	    {
		return vertex_t==NULL?0:vertex_t->get_rank();
	    }

	    /* Vertex address readout: */

	    const vertex_type* get_vertex() const
	    {
		return vertex_t;
	    }

	    /* CM_tag readout: */

	    bool is_marked() const
	    {
		return CM_tag;
	    }

	    /* Request whether this interaction computes the momentum
	     * flowing out of the vertex: */

	    bool computes_momentum() const
	    {
		return prop_policy[0];
	    }

	    /* Function computing the momentum of the produced current: */

	    void compute_momentum()
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		
		if(prop_policy[0])
		{
		    currents[produced_current]->momentum.assign((r_value_type)0);
		    for(size_type i=0;i<produced_current;++i)
		    {
			currents[produced_current]->momentum+=currents[i]->momentum;
		    }
		    for(size_type i=produced_current+1;i<currents.size();++i)
		    {
			currents[produced_current]->momentum+=currents[i]->momentum;
		    }
		    return;
		}
		if(prop_policy[1])
		{
		    currents[produced_current]->momentum=*(produced_momentum);
		    return;
		}
	    }

	    /* Function computing the decoupling tag: */

	    void compute_coupling_flag()
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
	
		coupled=vertex_t->is_coupled();	
		if(coupled)
		{
		    for(size_type i=0;i<produced_current;++i)
		    {
			coupled&=currents[i]->coupled;
		    }
		    for(size_type i=produced_current+1;i<currents.size();++i)
		    {
			coupled&=currents[i]->coupled;
		    }
		}
		if(prop_policy[0] or prop_policy[1])
		{
		    currents[produced_current]->coupled=coupled;
		}
		else
		{
		    currents[produced_current]->coupled|=coupled;
		}
	    }

	    /* Request whether this interaction copies the momentum flowing out
	     * of the vertex from another current: */

	    bool assigns_momentum() const
	    {
		return prop_policy[1];
	    }

	    /* Request whether the produced current propagates after the
	     * recursive relation is applied: */

	    bool propagates() const
	    {
		return prop_policy[2];
	    }

	    /* Function propagating the produced current: */

	    void propagate()
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		
		if(prop_policy[2])
		{
		    currents[produced_current]->propagate();
		}
	    }

	    /* Let the interaction object compute the produced momentum from the
	     * incoming ones: */

	    void set_computes_momentum()
	    {
		prop_policy.set(0);
		prop_policy.reset(1);
	    }

	    /* Let the interaction object copy the produced momentum from the
	     * dereferenced argument: */

	    void set_assigns_momentum(const momentum_type* q)
	    {
		prop_policy.reset(0);
		prop_policy.set(1);
		produced_momentum=q;
	    }

	    /* Let the interaction object propagate the produced current after
	     * the recursive relation is performed: */

	    void set_propagates()
	    {
		prop_policy.set(2);
	    }

	    /* Set the Fermi sign to -1: */

	    void set_Fermi_sign()
	    {
		Fermi_sign=-1;
	    }

	    /* Set the Fermi sign to +1: */

	    void reset_Fermi_sign()
	    {
		Fermi_sign=1;
	    }

	    /* Function resetting the produced current: */

	    void reset()
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		
		if(this->prop_policy[2])
		{
		    currents[produced_current]->reset();
		}
	    }

	    /* Function checking if the currents match the vertex legs: */

	    bool match_currents(const std::vector<current_iterator>& p) const
	    {
		if(vertex_t==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"no matching procedure could be started for NULL vertex type--returning false"<<endlog;
		    return false;
		}
		if(p.size()!=vertex_t->get_rank())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"number of currents to be matched "<<p.size()<<"does not match vertex rank"<<vertex_t->get_rank()<<"--returning false"<<endlog;
		    return false;
		}
		const particle_type* phi;
		for(size_type i=0;i<p.size();++i)
		{
		    phi=p[i]->get_produced_particle();
		    if(phi==NULL)
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle produced by current iterator "<<i<<" yielded NULL--returning false"<<endlog;
			return false;
		    }
		    if(i==produced_current)
		    {
			phi=phi->get_anti_particle();
			if(phi==NULL)
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle produced by current iterator "<<i<<" yielded NULL--returning false"<<endlog;
			    return false;
			}
		    }
		    if(phi!=vertex_t->get_leg(i))
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle "<<phi->get_name()<<" does not match leg "<<i<<" of "<<vertex_t->get_name()<<"--returning false"<<endlog;
			return false;
		    }
		}
		bit_string<N>b;
		for(size_type i=0;i<p.size();++i)
		{
		    if(i!=produced_current)
		    {
			b|=(p[i]->get_bit_string());
		    }
		}
		if(b != p[produced_current]->get_bit_string())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"bitstrings ";
		    for(size_type i=0;i<p.size();++i)
		    {
			if(i!=produced_current)
			{
			    log<<p[i]->get_bit_string()<<" ";
			}
		    }
		    log<<" do not add up to produced string "<<p[produced_current]->get_bit_string()<<"--returning false"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Current iterator insertion method. Returns a boolean that
	     * denoting if the insertion was successful: */
	    
	    void insert_currents(const std::vector<current_iterator>& p)
	    {
		if(match_currents(p))
		{
		    for(size_type i=0;i<p.size();++i)
		    {
			currents[i]=p[i];
			momenta[i]=&(p[i]->momentum);
			amp_iters[i]=currents[i]->amplitude.begin();
			if(p[i]->outgoing)
			{
			    flow.set(i);
			}
		    }
		}
	    }

	    void fetch_Feynman_rule()
	    {
		Feynman_rule=vertex_t->dispatch_Feynman_rule(flow,produced_current,swap_fermions);
		if(swap_fermions)
		{
		    std::swap(amp_iters[1],amp_iters[2]);
		}
	    }

	    /* Function computing the memory usage of all interaction objects
	     * together in the program: */

	    static std::size_t memory_usage()
	    {
		return sizeof(interaction_base<model_t,N>)*object_counter;
	    }

	    bool is_coupled() const
	    {
		return coupled;
	    }

	    /* Iterators to the produced current: */

	    current_iterator get_produced_current()
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		return currents[produced_current];
	    }
	    const_current_iterator get_produced_current() const
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		return currents[produced_current];
	    }

	    /* Produced particle: */

	    const particle_type* get_produced_particle() const 
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		return currents[produced_current]->get_produced_particle();
	    }

	    /* Produced bit string: */

	    bit_string<N> get_produced_bit_string() const
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		return currents[produced_current]->get_bit_string();
	    }

	    /* Incoming particles: */

	    std::vector<const particle_type*> get_incoming_particles() const
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		std::vector<const particle_type*>result;
		result.reserve(currents.size()-1);
		for(size_type i=0;i<currents.size();++i)
		{
		    if(i!=produced_current)
		    {
			result.push_back(currents[i]->get_produced_particle());
		    }
		}
		return result;
	    }

	    /* Incoming bit strings: */

	    std::vector< bit_string<N> > get_incoming_bit_strings() const
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		std::vector< bit_string<N> >result;
		result.reserve(currents.size()-1);
		for(size_type i=0;i<currents.size();++i)
		{
		    if(i!=produced_current)
		    {
			result.push_back(currents[i]->get_bit_string());
		    }
		}
		return result;
	    }

	    /* Diagram counting function: */

	    void count_diagrams()
	    {
		if(vertex_t->is_coupled())
		{
		    long long unsigned n=1;
		    for(size_type i=0;i<produced_current;++i)
		    {
			n*=(currents[i]->multiplicity);
		    }
		    for(size_type i=produced_current+1;i<currents.size();++i)
		    {
			n*=(currents[i]->multiplicity);
		    }
		    if(produced_current<currents.size())
		    {
			currents[produced_current]->multiplicity+=n;
		    }
		    else
		    {
			log(log_level::warning)<<"produced current "<<produced_current<<" out of range for vertex "<<*vertex_t<<endlog;
		    }
		}
	    }

	    /* Constant iterators to the vector of current iterators: */

	    typename std::vector<current_iterator>::const_iterator begin() const
	    {
		return currents.begin();
	    }
	    typename std::vector<current_iterator>::const_iterator end() const
	    {
		return currents.end();
	    }

	    /* CM marking function: */

	    void mark()
	    {
		CM_tag=true;
		for(size_type i=0;i<currents.size();++i)
		{
		    currents[i]->mark();
		}
	    }

	    /* Conditional reverse marking function: if the produced current is
	     * marked, the incoming ones are: */

	    void mark_if()
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		if(currents[produced_current]->is_marked())
		{
		    CM_tag=true;
		    for(size_type i=0;i<currents.size();++i)
		    {
			currents[i]->mark();
		    }
		}
	    }

	    /* Comparison operators: */

	    bool operator == (const interaction_base<model_t,N>& other) const
	    {
		return (currents==other.currents and produced_current==other.produced_current and vertex_t==other.vertex_t); 
	    }
	    bool operator != (const interaction_base<model_t,N>& other) const
	    {
		return !(this->operator==(other));
	    }
	    bool operator < (const interaction_base<model_t,N>& other) const
	    {
		if(this->get_produced_current()==other.get_produced_current())
		{
		    return std::lexicographical_compare(currents.begin(),currents.end(),other.currents.begin(),other.currents.end(),current_comparison);
		}
		return this->get_produced_current()<other.get_produced_current();
	    }
	    bool operator > (const interaction_base<model_t,N>& other) const
	    {
		return other<*this;
	    }
	    bool operator <= (const interaction_base<model_t,N>& other) const
	    {
		return !(other<*this);
	    }
	    bool operator >= (const interaction_base<model_t,N>& other) const
	    {
		return !(this->operator<(other));
	    }

	    /* Output method: */

	    std::ostream& print(std::ostream& os) const
	    {
		CAMGEN_ERROR_IF((produced_current>=currents.size()),"produced current out of range");
		if(this->is_coupled())
		{
		    current_iterator out_curr=currents[produced_current];
		    std::vector<current_iterator> in_currs(currents);
		    in_currs.erase(in_currs.begin()+produced_current);

		    os<<*out_curr<<"\t<--\t";
		    std::stringstream s;
		    s<<"(";
		    for(size_type i=0;i<in_currs.size()-1;++i)
		    {
			s<<in_currs[i]->get_bit_string()<<",";
		    }
		    s<<in_currs.back()->get_bit_string()<<")";
		    os<<std::setw(3*N+5)<<s.str()<<"\t";
		    if(Fermi_sign==1)
		    {
			os<<"+"<<"\t";
		    }
		    else
		    {
			os<<"-"<<"\t";
		    }
		    os<<flow[0]<<flow[1]<<flow[2]<<flow[3]<<"\t"; os<<prop_policy[0]<<prop_policy[1]<<prop_policy[2];
		}
		return os;
	    }

	protected:

	    /* Vertex address of the interaction: */

	    const vertex_type* vertex_t;
	    
	    /* Vector of iterators pointing to the interacting currents: */
	    
	    std::vector<current_iterator>currents;
	    
	    /* Vector of tensor iterators in the interacting subamplitudes: */
	    
	    std::vector<iterator>amp_iters;

	    /* Momentum addresses of the interacting currents: */
	    
	    std::vector<const momentum_type*>momenta;

	    /* Index of the current produced by the interaction: */ 

	    size_type produced_current;

	    /* Recursive relation: */

	    vert_func Feynman_rule;

	    /* Boolean denoting whether the 2 fermionic tensor iterators should
	     * be swapped (only applicable for boson-Majorana-Majorana
	     * vertices): */

	    bool swap_fermions;

	    /* Fermi sign of the vertex: */

	    int Fermi_sign;
	    
	    /* Momentum flows of the legs: if a bit is set, the current is
	     * external and outgoing in the process: */
	    
	    std::bitset<4>flow;

	    /* Utility boolean; the bit is used by the tree cleaning algorithm,
	     * where it will be set if the vertex can be connected to the final
	     * current in the tree: */

	    bool CM_tag;

	    /* Utility boolean keeping track of decoupled vertices in the tree:
	     * */

	    bool coupled;

	    /* Momentum propagation policy: the first bit denotes whether the
	     * produced momentum should be computed from the incoming currents,
	     * the second bit whether the this momentum is assigned by another
	     * current's momentum, and the third bit denotes if the produced
	     * current will propagate after the recursion: */

	    std::bitset<3>prop_policy;

	    /* Pointer to a momentum vector identical to the produced current's
	     * momentum. If the prop_policy has first bit zero and second bit
	     * set, the outgoing current's momentum will be copied from the
	     * dereferenced vector: */

	    const momentum_type* produced_momentum;

	    /* Static integer counting the interaction objects created by
	     * Camgen: */

	    static std::size_t object_counter;

	    /* Static current iterator comparison object: */

	    static current_iter_comp<model_t,N>current_comparison;

	    /* Iterators to the vector of current iterators: */

	    typename std::vector<current_iterator>::iterator begin()
	    {
		return currents.begin();
	    }
	    typename std::vector<current_iterator>::iterator end()
	    {
		return currents.end();
	    }

    };
    template<class model_t,std::size_t N>const bool interaction_base<model_t,N>::decomposes;
    template<class model_t,std::size_t N>std::size_t interaction_base<model_t,N>::object_counter=0;
    template<class model_t,std::size_t N>current_iter_comp<model_t,N> interaction_base<model_t,N>::current_comparison;

    /* Interaction subclass for models without colour flow decomposition: */

    template<class model_t,std::size_t N>class interaction<model_t,N,false>: public interaction_base<model_t,N>
    {
	public:

	    /* Base type definition: */

	    typedef interaction_base<model_t,N> base_type;
	    
	    /* Useful type definitions derived from the base type: */
	    
	    typedef typename base_type::r_value_type r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::particle_type particle_type;
	    typedef typename base_type::vertex_type vertex_type;
	    typedef typename base_type::current_iterator current_iterator;
	    typedef typename base_type::const_current_iterator const_current_iterator;

	    /* Friend class declaration: */

	    template<class mod,std::size_t N_in,std::size_t N_out>friend class process_tree;

	    /* Trivial constructor: */

	    interaction():base_type(){}

	    /* Regular constructor: */

	    interaction(const vertex_type* v,unsigned n=0):base_type(v,n){}

	    /* Copy constructor: */

	    interaction(const interaction<model_t,N,false>& other):base_type(other){}

	    /* Function that computes the produced current momentum, applies the
	     * recursive relation and propagates the produced current: */

	    void evaluate()
	    {
		this->compute_momentum();
		this->compute_coupling_flag();
		if(this->coupled)
		{
		    this->Feynman_rule(value_type(this->Fermi_sign,0),this->vertex_t->get_couplings(),this->amp_iters,this->momenta);
		}
		this->propagate();
	    }
    };

    /* Interaction subclass for models with colour flow decomposition: */

    template<class model_t,std::size_t N>class interaction<model_t,N,true>: public interaction_base<model_t,N>
    {
	public:

	    /* Base type definition: */

	    typedef interaction_base<model_t,N> base_type;
	    
	    /* Useful type definitions derived from the base type: */
	    
	    typedef typename base_type::r_value_type r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::particle_type particle_type;
	    typedef typename base_type::vertex_type vertex_type;
	    typedef typename base_type::current_iterator current_iterator;
	    typedef typename base_type::const_current_iterator const_current_iterator;

	    /* Iterator set type definitions: */
	    
	    typedef typename current<model_t,N,true>::iterset iterset;
	    typedef typename current<model_t,N,true>::iterset_iterator iterset_iterator;
	    typedef typename current<model_t,N,true>::const_iterset_iterator const_iterset_iterator;

	    /* Friend class declaration: */

	    template<class mod,std::size_t N_in,std::size_t N_out>friend class process_tree;

	    /* Trivial constructor: */

	    interaction():base_type(){}

	    /* Regular constructor: */

	    interaction(const vertex_type* v,unsigned n=0):base_type(v,n),set_iters(v->get_rank()-1){}

	    interaction(const interaction<model_t,N,true>& other):base_type(other),set_iters(other.set_iters){}

	    /* Function that computes the produced current momentum, applies the
	     * recursive relation and propagates the produced current: */

	    void evaluate()
	    {
		this->compute_momentum();
		this->compute_coupling_flag();
		if(this->coupled)
		{
		    if(first_iters())
		    {
			do
			{
			    import_iters();
			    this->Feynman_rule(value_type(this->Fermi_sign,0),this->vertex_t->get_couplings(),this->amp_iters,this->momenta,this->currents[this->produced_current]->amp_iters);
			}
			while(next_iters());
		    }
		}
		this->propagate();
	    }
	private:

	    /* Vector of propagating colours in the produced current subamplitude
	     * tensor: */
	    
	    std::vector<iterset_iterator>set_iters;

	    /* Function shifting all the incoming tensor iterators to the first
	     * propagating colour: */

	    bool first_iters()
	    {
		for(size_type i=0;i<this->produced_current;++i)
		{
		    if(this->currents[i]->amp_iters.size()==0)
		    {
			return false;
		    }
		    set_iters[i]=this->currents[i]->cfd_begin();
		}
		for(size_type i=this->produced_current+1;i<this->currents.size();++i)
		{
		    if(this->currents[i]->amp_iters.size()==0)
		    {
			return false;
		    }
		    set_iters[i-1]=this->currents[i]->cfd_begin();
		}
		return true;
	    }

	    /* Function shifting the incoming set of propagating iterators to
	     * the next colour combination: */

	    bool next_iters()
	    {
		for(size_type n=0;n<this->produced_current;++n)
		{
		    ++set_iters[n];
		    if(set_iters[n] == this->currents[n]->cfd_end())
		    {
			set_iters[n]=this->currents[n]->cfd_begin();
		    }
		    else
		    {
			return true;
		    }
		}
		for(size_type n=this->produced_current;n<set_iters.size();++n)
		{
		    ++set_iters[n];
		    if(set_iters[n] == this->currents[n+1]->cfd_end())
		    {
			set_iters[n]=this->currents[n+1]->cfd_begin();
		    }
		    else
		    {
			return true;
		    }
		}
		return false;
	    }

	    /* Function copying the set of propagating iterators to the argument
	     * iterator vector of the Feynman rule: */

	    void import_iters()
	    {
		for(size_type i=0;i<this->produced_current;++i)
		{
		    this->amp_iters[i]=*(set_iters[i]);
		}
		this->amp_iters[this->produced_current]=this->currents[this->produced_current]->amplitude.begin();
		for(size_type i=this->produced_current;i<set_iters.size();++i)
		{
		    this->amp_iters[i+1]=*(set_iters[i]);
		}
	    }
    };

    /* Outstream operator overload: */

    template<class model_t,std::size_t N>std::ostream& operator << (std::ostream& os,const interaction_base<model_t,N>& I)
    {
	return I.print(os);
    }

    /* Unary function object returning whether the argument is marked: */

    template<class T>class redundant
    {
	public:
	    bool operator () (const T& v) const
	    {
		return !(v.is_marked());
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_INTERACTION_H_*/
