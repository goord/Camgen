//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_CURRENT_H_
#define CAMGEN_CURRENT_H_

#include <set>
#include <Camgen/forward_decs.h>
#include <Camgen/particle.h>
#include <Camgen/bit_string.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Off-shell current class template definition. The off-shell currents are the   *
 * basic building blocks of the Camgen process tree, together with the          *
 * interaction class. There are two current classes, and which one is invoked    *
 * depends on the compile-time boolen 'decomposes', which represents whether we  *
 * are running in colour-flow decomposition mode or not. Both types are derived  *
 * from a common base class, incorporating the members that do not depend on     *
 * this mode.                                                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition and declaration of the current base type, incorporating members which do
     * not depend on the colour-flow decomposition mode. The class template takes 2
     * parameters: the model type and the number of particles involved in the process: */
    
    template<class model_t,std::size_t N>class current_base
    {
	public:
	    
	    /* Usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	    
	    /* Definitions of the particle and phase space types: */

	    typedef typename get_basic_types<model_t>::particle_type particle_type;
	    typedef typename get_basic_types<model_t>::phase_space_type phase_space_type;

	    /* Friend declaration of the interaction base class: */

	    friend class interaction_base<model_t,N>;
	    
	    /* Friend declarations of current_tree, process_tree and CM_algorithm classes: */
	    
	    template<class mod_t,std::size_t N_in,std::size_t N_out>friend class process_tree;
	    template<class mod_t,std::size_t N_in,std::size_t N_out>friend class current_tree;
	    template<class mod_t,std::size_t N_in,std::size_t N_out>friend class CM_algorithm;

	    /* Default constructor: */

	    current_base():particle_t(NULL),CM_tag(false),coupled(true),outgoing(false),initialised(false),final_amplitude(NULL),phase_space(NULL){}

	    /* Regular constructor, specifying the type of particle propagated by the
	     * current, the bitstring-coded momentum channel and an outgoing-particle
	     * boolean tag: */

	    current_base(const particle_type* phi,const bit_string<N>& b,bool out):particle_t(phi),bitstring(b),CM_tag(false),coupled(true),initialised(false),final_amplitude(NULL),phase_space(NULL)
	    {
		if(b.count()<=1)
		{
		    outgoing=out;
		}
		else
		{
		    outgoing=false;
		}
	    }

	    /* Copy constructor: */

	    current_base(const current_base<model_t,N>& other):momentum(other.momentum),particle_t(other.particle_t),bitstring(other.bitstring),CM_tag(other.CM_tag),coupled(other.coupled),outgoing(other.outgoing),initialised(other.initialised),amplitude(other.amplitude),multiplicity(other.multiplicity),final_amplitude(other.final_amplitude),phase_space(NULL)
	    {
		if(other.phase_space!=NULL)
		{
		    allocate_phase_space();
		}
	    }
	    /* Destructor: */

	    ~current_base()
	    {
		if(phase_space!=NULL)
		{
		    delete phase_space;
		}
	    }
	    
	    /* Particle pointer readout: */
	    
	    const particle_type* get_particle_type() const
	    {
		return particle_t;
	    }

	    /* Equivalent incoming particle pointer readout: */

	    const particle_type* get_produced_particle() const
	    {
		if(!outgoing)
		{
		    return particle_t;
		}
		else
		{
		    return particle_t->get_anti_particle();
		}
	    }

	    /* Marker utility function: */

	    void mark()
	    {
		CM_tag=true;
	    }
	    
	    /* Unmarking utility function: */
	    
	    void unmark()
	    {
		CM_tag=false;
	    }

	    /* Marker utility tag readout: */

	    bool is_marked() const
	    {
		return CM_tag;
	    }

	    /* Momentum-channel bitstring readout: */

	    bit_string<N> get_bit_string() const
	    {
		return bitstring;
	    }

	    /* Multiplicity readout (used for counting diagrams): */

	    long long unsigned get_multiplicity() const
	    {
		return multiplicity;
	    }

	    /* Amplitude tensor readout: */

	    const tensor_type& get_amplitude() const
	    {
		return amplitude;
	    }

	    /* Fermionic property readout: */

	    bool is_fermion() const
	    {
		return particle_t->is_fermion();
	    }

	    /* Bosonic property readout: */

	    bool is_boson() const
	    {
		return !(this->is_fermion());
	    }

	    /* Direction of momentum readout: */

	    bool is_outgoing() const
	    {
		return outgoing;
	    }
	    
	    /* External property readout: */

	    bool is_external() const
	    {
		return (bitstring.count()<=1);
	    }
	    
	    /* Initialisation readout: */
	    
	    bool is_initialised() const
	    {
		return initialised;
	    }

	    /* Argument setting for the final contraction (only makes sense if the current
	     * is the final internal current matching the final external leg: */

	    void set_argument(tensor_type* t)
	    {
		if(bitstring.count() != 0)
		{
		    return;
		}
		final_amplitude=t;
	    }

	    /* Floating-point data output: */

	    std::ostream& print_amp(std::ostream& os) const
	    {
		if(initialised)
		{
		    os<<std::endl<<amplitude<<std::endl;
		    os<<"momentum: "<<momentum<<std::endl;
		}
		return os;
	    }

	    template<class rep_t>std::ostream& print_cf_amp(std::ostream& os) const
	    {
		if(initialised)
		{
		    if(amplitude.size()%rep_t::group_type::rank==0)
		    {
			std::size_t d=amplitude.size()/rep_t::group_type::rank;
			tensor_type cf_amp(3,rep_t::dimension,rep_t::dimension,d);
			iterator it0=cf_amp.begin();
			const_iterator it1=amplitude.begin();
			for(unsigned A=0;A<rep_t::dimension;++A)
			{
			    for(unsigned B=0;B<rep_t::dimension;++B)
			    {
				for(unsigned a=0;a<rep_t::group_type::rank;++a)
				{
				    for(unsigned mu=0;mu<d;++mu)
				    {
					it0[mu]+=rep_t::template implementation<r_value_type>::generator(a,A,B)*it1[mu];
				    }
				    it1+=d;
				}
				it1-=(d*(rep_t::group_type::rank));
				it0+=d;
			    }
			}
			os<<std::endl<<cf_amp<<std::endl;
			os<<"momentum: "<<momentum<<std::endl;
		    }
		    else
		    {
			os<<std::endl<<amplitude<<std::endl;
			os<<"momentum: "<<momentum<<std::endl;
		    }
		}
		return os;
	    }
	    
	    /* Iterator to the first entry of the subamplitude: */
	    
	    iterator get_amplitude_iterator()
	    {
		return amplitude.begin();
	    }

	    /* Comparison operators. Policy is that the momentum channel bitstring
	     * determines ordering before the particle type: */

	    bool operator == (const current_base<model_t,N>& other) const
	    {
		return (bitstring==other.bitstring and particle_t==other.particle_t);
	    }
	    bool operator != (const current_base<model_t,N>& other) const
	    {
		return !(this->operator==(other));
	    }
	    bool operator < (const current_base<model_t,N>& other) const
	    {
		if(bitstring==other.bitstring)
		{
		    CAMGEN_ERROR_IF((particle_t==NULL),"attempt to dereference a NULL particle pointer...");
		    return (particle_t->get_flavour() < other.particle_t->get_flavour());
		}
		else
		{
		    return bitstring<other.bitstring;
		}
	    }
	    bool operator > (const current_base<model_t,N>& other) const
	    {
		return other<(*this);
	    }
	    bool operator <= (const current_base<model_t,N>& other) const
	    {
		return !(other<(*this));
	    }
	    bool operator >= (const current_base<model_t,N>& other) const
	    {
		return !(this->operator<(other));
	    }

	protected:

	    /* Momentum vector flowing through the current (reversed for outgoing
	     * particles): */

	    momentum_type momentum;
	    
	    /* Particle corresponding to the current: */

	    const particle_type* particle_t;
	    
	    /* Bitstring denoting the momentum channel in the process: */
	    
	    bit_string<N> bitstring;

	    /* Utility boolean: */

	    bool CM_tag;

	    /* Decoupling utlilty boolean: */

	    bool coupled;
	    
	    /* Boolean denoting the momentum flow direction: */
	    
	    bool outgoing;

	    /* Initialisation tag: */

	    bool initialised;
	    
	    /* Tensor denoting the subamplitude: */
	    
	    tensor_type amplitude;

	    /* Integer denoting the multiplicity of the current, i.e. the number of
	     * Feynman diagrams it participates in: */

	    long long unsigned multiplicity;

	    /* Address of the final internal subamplitude. This member is NULL for all
	     * internal currents and non-final external currents...*/
	    
	    tensor_type* final_amplitude;
	    
	    /* Phase space object, only allocatable and callable for external currents: */

	    phase_space_type* phase_space;
	    
	    /* Allocate phase space object: */

	    void allocate_phase_space()
	    {
		if(phase_space==NULL)
		{
		    phase_space=new phase_space_type(particle_t,outgoing);
		}
	    }
	    
	    /* Allocate phase space object again: */

	    void reallocate_phase_space()
	    {
		if(phase_space!=NULL)
		{
		    delete phase_space;
		}
		phase_space=new phase_space_type(particle_t,outgoing);
	    }

	    /* Return phase space address: */

	    phase_space_type* get_phase_space()
	    {
		return phase_space;
	    }

	    /* Return phase space const address: */

	    const phase_space_type* get_phase_space() const
	    {
		return phase_space;
	    }

	    /* Momentum direction change: */

	    void set_incoming()
	    {
		if(outgoing)
		{
		    outgoing=false;
		    if(phase_space!=NULL)
		    {
			phase_space->set_incoming(); 
		    }
		}
	    }
	    
	    void set_outgoing()
	    {
		if(!outgoing)
		{
		    outgoing=true;
		    if(phase_space!=NULL)
		    {
			phase_space->set_outgoing(); 
		    }
		}
	    }
    };

    /* Derived current class for non-colour-decomposed calculations: */

    template<class model_t,std::size_t N>class current<model_t,N,false>: public current_base<model_t,N>
    {
	public:

	    /* Base type definition: */

	    typedef current_base<model_t,N> base_type;
	    
	    /* Useful type definitions all inherited: */
	    
	    typedef typename base_type::r_value_type r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::particle_type particle_type;

	    /* Non-colour-decomposed interaction class declared friend: */

	    friend class interaction<model_t,N,false>;

	    /* Default constructor: */

	    current():base_type(){}
	    
	    /* Regular constructor: */
	    
	    current(const particle_type* phi,const bit_string<N>& b,bool ext):base_type(phi,b,ext){}

	    /* Copy-constructor: */

	    current(const current<model_t,N,false>& other):base_type(other){}

	    /* Reset subamplitude to start new event: */

	    void reset()
	    {
		this->amplitude.reset();
	    }

	    /* Initialisation function: */

	    void initialise()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"initialisation called without particle type instance...");
		if(!this->initialised)
		{
		    /* Allocate subamplitude: */

		    this->particle_t->make_amplitude(this->amplitude);

		    /* Allocate phase space: */

		    if(this->is_external())
		    {
			this->allocate_phase_space();
		    }
		    this->initialised=true;
		}
	    }

	    /* External wave function computation: */

	    void evaluate()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function evaluation called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function evaluation called without phase-space instance...");
		
		this->phase_space->copy_momentum(this->momentum);
		if(this->particle_t->is_coupled())
		{
		    this->phase_space->fill_wave_function(this->amplitude.begin());
		}
	    }

	    /* Initialise wave function computation in summation mode: */

	    void init_dof_sum(bool spin_summing,bool col_summing)
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function spin sum initialisation called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function spin sum initialisation called without phase-space instance...");

		this->phase_space->copy_momentum(this->momentum);
		if(this->particle_t->is_coupled())
		{
		    this->phase_space->init_dof_sum(this->amplitude.begin(),spin_summing,col_summing);
		}
	    }
	    
	    /* External wave function computation in spin-summation mode: */

	    bool evaluate_dof_sum(bool spin_summing,bool col_summing)
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function spin sum evaluation called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function spin sum evaluation called without phase-space instance...");
		
		if(this->particle_t->is_coupled())
		{
		    if(spin_summing or col_summing)
		    {
			this->phase_space->sum_wave_function(this->amplitude.begin(),spin_summing,col_summing);
		    }
		    return this->phase_space->minimal_dof(spin_summing,col_summing);
		}
		return true;
	    }

	    /* Function copying the momentum in the phase space instance to the
	     * current, adding the correct sign when the particle is outgoing:
	     * */

	    void copy_momentum()
	    {
		CAMGEN_ERROR_IF((this->phase_space==NULL),"momentum copy called without phase-space instance...");

		this->phase_space->copy_momentum(this->momentum);
	    }

	    /* Wave-function-final subamplitude contraction: */

	    value_type contract_wave_function()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function contraction called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function contraction called without phase-space instance...");
		
		if(this->particle_t->is_coupled())
		{
		    return this->phase_space->contract(this->final_amplitude->begin(),this->amplitude.begin());
		}
		else
		{
		    return 0;
		}
	    }

	    /* Additional factor to obtain the correct colour sum: */

	    r_value_type colour_sum_factor() const
	    {
		return (r_value_type)1;
	    }

	    /* Propagation of the subamplitude: */

	    void propagate()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"propagator called for undefined particle type...");
		
		if(this->is_external())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"propagator called for internal current "<<*this<<endlog;
		}
		
		if(this->particle_t->is_coupled())
		{
		    this->particle_t->refresh_propagator(&(this->momentum));
		    this->particle_t->propagate(this->amplitude.begin(),this->amplitude.end());
		}
	    }

	    /* Output method: */

	    std::ostream& print(std::ostream& os) const
	    {
		if(this->particle_t!=NULL and this->particle_t->is_coupled())
		{
		    os<<this->bitstring<<"\t"<<std::setw(6)<<((this->outgoing)?"-":"+")<<std::setw(10)<<std::left<<this->particle_t->get_name();
		}
		return os;
	    }

	    /* Comparisaon operator overloads: */

	    bool operator == (const current<model_t,N,false>& other) const
	    {
		return this->base_type::operator==(other);
	    }
	    bool operator != (const current<model_t,N,false>& other) const
	    {
		return this->base_type::operator!=(other);
	    }
	    bool operator < (const current<model_t,N,false>& other) const
	    {
		return this->base_type::operator<(other);
	    }
	    bool operator > (const current<model_t,N,false>& other) const
	    {
		return this->base_type::operator>(other);
	    }
	    bool operator <= (const current<model_t,N,false>& other) const
	    {
		return this->base_type::operator<=(other);
	    }
	    bool operator >= (const current<model_t,N,false>& other) const
	    {
		return this->base_type::operator>=(other);
	    }
    };

    /* Derived current class for colour-decomposed calculations: */

    template<class model_t,std::size_t N>class current<model_t,N,true>: public current_base<model_t,N>
    {
	public:
	    
	    /* Base class type definition: */

	    typedef current_base<model_t,N> base_type;
	    
	    /* Useful type definitions inherited from base type: */
	    
	    typedef typename base_type::r_value_type r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::particle_type particle_type;

	    /* Additional type definitions involving stl sets of tensor
	     * iterators: */

	    typedef std::set<iterator> iterset;
	    typedef typename iterset::iterator iterset_iterator;
	    typedef typename iterset::const_iterator const_iterset_iterator;

	    /* Declaring the colour-decomposed interactioon class a friend: */
	    
	    friend class interaction<model_t,N,true>;

	    /* Default constructor: */

	    current():base_type(){}
	    
	    /* Regular constructor: */
	    
	    current(const particle_type* phi,const bit_string<N>& b,bool ext):base_type(phi,b,ext){}
	    
	    /* Copy constructor: */
	    
	    current(const current<model_t,N,true>& other):base_type(other),amp_iters(other.amp_iters){}
	    
	    /* Subamplitude resetting function, only resetting the propagating
	     * colour modes: */
	    
	    void reset()
	    {
		for(iterset_iterator it=amp_iters.begin();it != amp_iters.end();++it)
		{
		    this->particle_t->reset(*it);
		}
		amp_iters.clear();
	    }

	    /* Initialisation function: */

	    void initialise()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"initialisation called without particle type instance...");
		
		if(!this->initialised)
		{
		    /* Allocate subamplitude: */

		    this->particle_t->make_amplitude(this->amplitude);
		    
		    /* Allocate phase space: */
		    
		    if(this->is_external())
		    {
			this->allocate_phase_space();
		    }
		    this->initialised=true;
		}
	    }

	    /* External wave function computation: */

	    void evaluate()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function evaluation called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function evaluation called without phase-space instance...");
		
		this->phase_space->copy_momentum(this->momentum);
		if(this->particle_t->is_coupled())
		{
		    amp_iters.insert(this->phase_space->fill_wave_function(this->amplitude.begin()));
		}
	    }

	    /* Initialise wave function computation in summation mode: */

	    void init_dof_sum(bool spin_summing,bool col_summing)
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function spin sum initialisation called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function spin sum initialisation called without phase-space instance...");

		this->phase_space->copy_momentum(this->momentum);
		if(this->particle_t->is_coupled())
		{
		    amp_iters.insert(this->phase_space->init_dof_sum(this->amplitude.begin(),spin_summing,col_summing));
		}
	    }
	    
	    /* External wave function computation in spin-summation mode: */

	    bool evaluate_dof_sum(bool spin_summing,bool col_summing)
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function spin sum evaluation called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function spin sum evaluation called without phase-space instance...");
		
		if(this->particle_t->is_coupled())
		{
		    if(col_summing or spin_summing)
		    {
			amp_iters.clear();
			amp_iters.insert(this->phase_space->sum_wave_function(this->amplitude.begin(),spin_summing,col_summing));
		    }
		    return this->phase_space->minimal_dof(spin_summing,col_summing);
		}
		return true;
	    }

	    /* Function copying the momentum in the phase space instance to the
	     * current, adding the correct sign when the particle is outgoing:
	     * */

	    void copy_momentum()
	    {
		CAMGEN_ERROR_IF((this->phase_space==NULL),"momentum copy called without phase-space instance...");

		this->phase_space->copy_momentum(this->momentum);
	    }

	    /* Wave-function-final subamplitude contraction: */

	    value_type contract_wave_function()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"wave function contraction called without particle type instance...");
		CAMGEN_ERROR_IF((this->phase_space==NULL),"wave function contraction called without phase-space instance...");
		
		if(this->particle_t->is_coupled() and amp_iters.size()>0)
		{
		    return this->phase_space->contract(this->final_amplitude->begin(),this->amplitude.begin());
		}
		return value_type(0,0);
	    }

	    /* Additional factor to obtain the correct colour sum: */

	    r_value_type colour_sum_factor() const
	    {
		if(this->particle_t==NULL)
		{
		    return (r_value_type)1;
		}
		return (this->particle_t->get_decomposed_colours()==2)?((r_value_type)0.5):((r_value_type)1);
	    }

	    /* Propagation of the subamplitude: */

	    void propagate()
	    {
		CAMGEN_ERROR_IF((this->particle_t==NULL),"propagator called for undefined particle type...");
		
		if(this->is_external())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"propagator called for internal current "<<*this<<endlog;
		}
		
		this->particle_t->refresh_propagator(&(this->momentum));
		if(this->particle_t->is_coupled())
		{
		    for(iterset_iterator it=amp_iters.begin();it != amp_iters.end();++it)
		    {
			this->particle_t->propagate(*it);
		    }
		}
	    }

	    /* Output method: */

	    std::ostream& print(std::ostream& os) const
	    {
		if(this->particle_t!=NULL and this->particle_t->is_coupled())
		{
		    os<<this->bitstring<<"\t"<<std::setw(6)<<((this->outgoing)?"-":"+")<<std::setw(10)<<std::left<<this->particle_t->get_name()<<"\t";
		    if(this->particle_t->get_decomposed_colours()>0)
		    {
			for(const_iterset_iterator it=amp_iters.begin();it != amp_iters.end();++it)
			{
			    os<<"(";
			    for(size_type i=0;i<this->particle_t->get_decomposed_colours()-1;++i)
			    {
				os<<it->index(i)<<",";
			    }
			    os<<it->index(this->particle_t->get_decomposed_colours()-1)<<")";
			}
		    }
		    else
		    {
			os<<"(0)";
		    }
		}
		return os;
	    }

	    /* Overloaded comparison operators: */

	    bool operator == (const current<model_t,N,true>& other) const
	    {
		return this->base_type::operator==(other);
	    }
	    bool operator != (const current<model_t,N,true>& other) const
	    {
		return this->base_type::operator!=(other);
	    }
	    bool operator < (const current<model_t,N,true>& other) const
	    {
		return this->base_type::operator<(other);
	    }
	    bool operator > (const current<model_t,N,true>& other) const
	    {
		return this->base_type::operator>(other);
	    }
	    bool operator <= (const current<model_t,N,true>& other) const
	    {
		return this->base_type::operator<=(other);
	    }
	    bool operator >= (const current<model_t,N,true>& other) const
	    {
		return this->base_type::operator>=(other);
	    }

	    /* Constant iterators in the set of propagating colour modes: */

	    const_iterset_iterator cfd_begin() const
	    {
		return amp_iters.begin();
	    }
	    const_iterset_iterator cfd_end() const
	    {
		return amp_iters.end();
	    }

	private:

	    /* Set of propagating colour modes (i.e. iterators at nonzero
	     * subtensors of the subamplitude: */

	    iterset amp_iters;
	    
	    /* Non-constant iterators in the set of propagating colour modes: */
	    
	    iterset_iterator cfd_begin()
	    {
		return amp_iters.begin();
	    }
	    iterset_iterator cfd_end()
	    {
		return amp_iters.end();
	    }
    };

    /* Overloaded I/O-stream operator: */

    template<class model_t,std::size_t N,bool decomp>std::ostream& operator << (std::ostream& os,const current<model_t,N,decomp>& c)
    {
	return c.print(os);
    }
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_CURRENT_H_*/

