//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file phase_space.h
    \brief contains the phase space classes
 */

#ifndef CAMGEN_PHASE_SPACE_H_
#define CAMGEN_PHASE_SPACE_H_

#include <algorithm>
#include <Camgen/unused.h>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/forward_decs.h>
#include <Camgen/vector.h>
#include <Camgen/tensor.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definitions of the particle phase space types. All phase space classes are    *
 * derived from a common base class holding the (possibly zero-dimensional)      *
 * momentum degrees of freedom. Then there are the derived class template        *
 * containing the helicity and colour degrees of freedom. The eventual particle  *
 * phase space classes are derived from both these types, unless the number of   *
 * colours or the spacetime dimension is zero.                                   *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Phase space base class template. The template
    /// parameter denotes the model type. The data stored in this class are a
    /// momentum and a particle pointer.

    template<class model_t>class base_ps
    {
	public:

	    template<class M,std::size_t N>friend class current_base;
	    
	    /* Some useful type definitions: */

	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;

	    /* Particle type: */

	    const particle<model_t>* const particle_type; 
	    
	    /* Maximum helicity integer: */
	    
	    const int max_helicity;

	    /* Spacetime tensor size: */

	    const size_type spacetime_tensor_size;

	    /* Colour tensor size: */

	    const size_type colour_tensor_size;

	    /* amplitude size: */

	    const size_type tensor_size;

	    /* Constructor taking the particle type and momentum direction s
	     * arguments: */

	    base_ps(const particle<model_t>* phi,bool out=false):particle_type(phi),max_helicity(phi->is_boson()?((phi->get_spin().twice())/2):(phi->get_spin().twice())),spacetime_tensor_size(phi->get_spacetime_tensor_size()),colour_tensor_size(phi->get_colour_tensor_size()),tensor_size(spacetime_tensor_size*colour_tensor_size),outgoing(out),handlebar(false){}

	    /* Destructor: */

	    virtual ~base_ps(){}
	    
	    /// Returns a reference to the particle momentum.
	    
	    momentum_type& momentum()
	    {
		return p;
	    }
	    
	    /// Returns a constant reference to the particle momentum.
	    
	    const momentum_type& momentum() const
	    {
		return p;
	    }
	    
	    /// Returns a pointer to the particle momentum.
	    
	    momentum_type* get_momentum()
	    {
		return &p;
	    }
	    
	    /// Returns a const pointer to the particle momentum.
	    
	    const momentum_type* get_momentum() const
	    {
		return &p;
	    }

	    /// Returns a reference to the i-th component of the particle momentum.

	    r_value_type& momentum(size_type i)
	    {
		return p[i];
	    }

	    /// Returns a constant reference to the i-th component of the particle momentum.

	    const r_value_type& momentum(size_type i) const
	    {
		return p[i];
	    }

	    /// Returns the colour rank of the particle.

	    size_type colour_rank() const
	    {
		return particle_type->colour_index_ranges.size();
	    }

	    /// Returns the dimension of the i-th colour index of the particle.

	    size_type colour_range(size_type i) const
	    {
		return particle_type->colour_index_ranges[i];
	    }

	    /// Returns the flavour of the particle.

	    size_type flavour() const
	    {
		return particle_type->flavour;
	    }

	    /* Writing the momentum to the argument vector, flipping if the
	     * particle is outgoing: */

	    void copy_momentum(momentum_type& q)
	    {
		if(outgoing)
		{
		    q=-p;
		}
		else
		{
		    q=p;
		}
	    }

	    /// Returns the particle mass.

	    const r_value_type& mass() const
	    {
		return particle_type->get_mass();
	    }

	    /// Returns particle mass const pointer.

	    const r_value_type* get_mass() const
	    {
		return particle_type->get_mass_address();
	    }

	    /// Returns the number of helicities of the particle.

	    size_type helicities() const
	    {
		return particle_type->spin_dof();
	    }

	    /// Returns the maximal helicity of the particle.

	    int maximal_helicity() const
	    {
		return max_helicity;
	    }

	    bool massive() const
	    {
		return (particle_type->mass!=NULL);
	    }

	    bool massless() const
	    {
		return (particle_type->mass==NULL);
	    }

	    /// Returns true if the particle allows a zero-helicity state.

	    bool has_zero_helicity() const
	    {
		if(particle_type->mass==NULL)
		{
		    return (particle_type->zero_hel_outgoing_massless_state!=NULL);
		}
		return (particle_type->zero_hel_outgoing_massive_state!=NULL);
	    }

	    /// Returns whether the state is incoming.

	    bool is_incoming() const
	    {
		return !outgoing;
	    }

	    /// Returns whether the state is outgoing.

	    bool is_outgoing() const
	    {
		return outgoing;
	    }

	    /// Returns the handlebar flag.

	    bool is_handlebar() const
	    {
		return handlebar;
	    }

	    /// Replaces the polarisation vector by the momentum for spin-1
	    /// particles.

	    void set_handlebar()
	    {
		if(particle_type->get_spin()==1)
		{
		    handlebar=true;
		}
	    }

	    /// Resets the handlebar flag.

	    void unset_handlebar()
	    {
		handlebar=false;
	    }

	protected:

	    /* Particle momentum: */

	    momentum_type p;

	    /* Direction of the momentum flow: */

	    bool outgoing;

	    /* Gauge switch: */

	    bool handlebar;

	    /* Sets the particle incoming: */

	    void set_incoming()
	    {
		outgoing=false;
	    }

	    /* Sets the particle outgoing: */

	    void set_outgoing()
	    {
		outgoing=true;
	    }
    };

    /// Helicity phase space class template. Takes as template parameters the model type
    /// and the compile-time continuous_helicities boolean that should be
    ///  defined in model_t.

    template<class model_t,bool cont_helicities=model_t::continuous_helicities>class helicity_ps;

    /// Helicity phase space template specialisation for discrete-helicity
    /// model types.

    template<class model_t>class helicity_ps<model_t,false>
    {
	public:

	    /* Some useful type definitions: */

	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;

	    /* Wave function pointer type definition: */

	    typedef void(*wave_func)(iterator,const momentum_type*,const r_value_type*);

	    /* Phase space instance: */

	    const base_ps<model_t>* const ps;

	    /* Constructor taking a ps base instance as argument: */
	    
	    helicity_ps(const base_ps<model_t>* ps_):ps(ps_),hel(-(ps->max_helicity)){}
	    
	    /// Returns a reference to the helicity integer member.
	    
	    int& helicity()
	    {
		return hel;
	    }

	    /// Returns a constant reference to the helicity integer member.

	    const int& helicity() const
	    {
		return hel;
	    }

	    /* Wave function caller: */

	    iterator fill_wave_function(iterator it)
	    {
		if(ps->is_handlebar())
		{
		    for(size_type mu=0;mu<model_t::dimension;++mu)
		    {
			it[mu]=ps->momentum(mu);
		    }
		    return it;
		}
		const particle<model_t>* phi(ps->particle_type);
		if(ps->massless())
		{
		    if(ps->is_outgoing())
		    {
			if(hel>0)
			{
			    CAMGEN_ERROR_IF((hel>(int)(phi->pos_hel_outgoing_massless_states).size()),"positive helicity out of range");
			    (phi->pos_hel_outgoing_massless_states)[hel-1](it,&(ps->momentum()),NULL);
			    return it;
			}
			if(hel==0)
			{
			    CAMGEN_ERROR_IF((phi->zero_hel_outgoing_massless_state==NULL),"undefined massless spin-0 state encountered");
			    (phi->zero_hel_outgoing_massless_state)(it,&(ps->momentum()),NULL);
			    return it;
			}
			if(hel<0)
			{
			    CAMGEN_ERROR_IF((hel<-(int)(phi->neg_hel_outgoing_massless_states).size()),"negative helicity out of range");
			    (phi->neg_hel_outgoing_massless_states)[-hel-1](it,&(ps->momentum()),NULL);
			    return it;
			}
		    }
		    else
		    {
			if(hel>0)
			{
			    CAMGEN_ERROR_IF((hel>(int)(phi->pos_hel_incoming_massless_states).size()),"positive helicity out of range");
			    (phi->pos_hel_incoming_massless_states)[hel-1](it,&(ps->momentum()),NULL);
			    return it;
			}
			if(hel==0)
			{
			    CAMGEN_ERROR_IF((phi->zero_hel_incoming_massless_state==NULL),"undefined massless spin-0 state encountered");
			    (phi->zero_hel_incoming_massless_state)(it,&(ps->momentum()),NULL);
			    return it;
			}
			if(hel<0)
			{
			    CAMGEN_ERROR_IF((hel<-(int)(phi->neg_hel_incoming_massless_states).size()),"negative helicity out of range");
			    (phi->neg_hel_incoming_massless_states)[-hel-1](it,&(ps->momentum()),NULL);
			    return it;
			}
		    }
		}
		else
		{
		    if(ps->is_outgoing())
		    {
			if(hel>0)
			{
			    CAMGEN_ERROR_IF((hel>(int)(phi->pos_hel_outgoing_massive_states).size()),"positive helicity out of range");
			    (phi->pos_hel_outgoing_massive_states)[hel-1](it,&(ps->momentum()),phi->mass);
			    return it;
			}
			if(hel==0)
			{
			    CAMGEN_ERROR_IF((phi->zero_hel_outgoing_massive_state==NULL),"undefined massive spin-0 state encountered");
			    (phi->zero_hel_outgoing_massive_state)(it,&(ps->momentum()),phi->mass);
			    return it;
			}
			if(hel<0)
			{
			    CAMGEN_ERROR_IF((hel<-(int)(phi->neg_hel_outgoing_massive_states).size()),"negative helicity out of range");
			    (phi->neg_hel_outgoing_massive_states)[-hel-1](it,&(ps->momentum()),phi->mass);
			    return it;
			}
		    }
		    else
		    {
			if(hel>0)
			{
			    CAMGEN_ERROR_IF((hel>(int)(phi->pos_hel_incoming_massive_states).size()),"positive helicity out of range");
			    (phi->pos_hel_incoming_massive_states)[hel-1](it,&(ps->momentum()),phi->mass);
			    return it;
			}
			if(hel==0)
			{
			    CAMGEN_ERROR_IF((phi->zero_hel_incoming_massive_state==NULL),"undefined massive zero-helicity state encountered");
			    (phi->zero_hel_incoming_massive_state)(it,&(ps->momentum()),phi->mass);
			    return it;
			}
			if(hel<0)
			{
			    CAMGEN_ERROR_IF((hel<-(int)(phi->neg_hel_incoming_massive_states).size()),"negative helicity out of range");
			    (phi->neg_hel_incoming_massive_states)[-hel-1](it,&(ps->momentum()),phi->mass);
			    return it;
			}
		    }
		}
		return it;
	    }

	    /* Initialisation of helicity summation mode: */

	    iterator init_sum_wave_function(iterator it)
	    {
		hel=-ps->max_helicity;
		iterator itend=it+ps->spacetime_tensor_size;
		for(iterator it2=it;it2!=itend;++it2)
		{
		    *it2=value_type(0,0);
		}
		return fill_wave_function(it);
	    }

	    /* Wave function spin sum function: */
	    
	    iterator sum_wave_function(iterator it)
	    {
		++hel;
		if(hel==0 and !ps->particle_type->has_zero_helicity_state())
		{
		    hel=1;
		}
		if(hel>ps->max_helicity)
		{
		    hel=-hel+1;
		}
		for(size_type i=0;i<ps->spacetime_tensor_size;++i)
		{
		    it[i]=value_type(0,0);
		}
		return fill_wave_function(it);
	    }

	    /* Minimal helicity output boolean: */

	    bool minimal_helicity() const
	    {
		return (hel==-(ps->max_helicity));
	    }

	    /* Subamplitude contraction caller: */

	    value_type contract(const_iterator internal,const_iterator external)
	    {
		value_type amp(0,0);
		if(ps->is_outgoing())
		{
		    ps->particle_type->contraction(amp,internal,external);
		    return amp;
		}
		ps->particle_type->conj_contraction(amp,internal,external);
		return amp;
	    }

	    /* Printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<hel;
		return os;
	    }

	protected:

	    /* Helicity: */

	    int hel;
    };

    /// Continuous helicity data holder type.

    template<class value_t>class helicity_phases
    {
	public:

	    /* Useful type definition: */

	    typedef value_t r_value_type;
	    typedef std::complex<value_t> value_type;
	    typedef typename std::vector<value_type>::size_type size_type;

	    /// Returns a scalar helicity holder.

	    static helicity_phases<value_t> create_scalar()
	    {
		helicity_phases result;
		return result;
	    }

	    /// Returns a massless spin-s helicity holder.

	    static helicity_phases<value_t> create_massless(int s)
	    {
		helicity_phases result(s,false);
		return result;
	    }

	    /// Returns a massive spin-s helicity holder.

	    static helicity_phases<value_t> create_massive(int s)
	    {
		helicity_phases result(s,true);
		return result;
	    }

	    /// Flag denoting whether we have a zero helicity phase.

	    const bool zero_helicity_state;

	    /// Default constructor (constructs scalar state).
	    
	    helicity_phases():zero_hel(1,0),zero_helicity_state(true){}

	    /// Constructor with maximal helicity max_hel.

	    helicity_phases(int max_hel,bool massive=false):zero_helicity_state(massive)
	    {
		if(max_hel==0)
		{
		    zero_hel=value_type(1,0);
		}
		else
		{
		    pos_hels.assign(std::abs(max_hel),value_type(1,0));
		    neg_hels.assign(std::abs(max_hel),value_type(1,0));
		    zero_hel=massive?value_type(1,0):value_type(0,0);
		}
	    }
	    
	    /// Helicity phase access.
	    
	    value_type& helicity_phase(int i)
	    {
		if(i>0)
		{
		    return pos_hels[i-1];
		}
		if(i<0)
		{
		    return neg_hels[-i-1];
		}
		return zero_hel;
	    }
	    
	    /// Positive helicity phase access.
	    
	    value_type& positive_helicity_phase(int i)
	    {
		return pos_hels[i];
	    }
	    
	    /// Negative helicity phase access.
	    
	    value_type& negative_helicity_phase(int i)
	    {
		return neg_hels[i];
	    }
	    
	    /// Zero helicity phase access.
	    
	    value_type& zero_helicity_phase()
	    {
		return zero_hel;
	    }

	    /// Const helicity phase access.
	    
	    const value_type& helicity_phase(int i) const
	    {
		if(i>0)
		{
		    return pos_hels[i-1];
		}
		if(i<0)
		{
		    return neg_hels[-i-1];
		}
		return zero_hel;
	    }
	    
	    /// Const positive helicity phase access.
	    
	    const value_type& positive_helicity_phase(int i) const
	    {
		return pos_hels[i];
	    }
	    
	    /// Const negative helicity phase access.
	    
	    const value_type& negative_helicity_phase(int i) const
	    {
		return neg_hels[i];
	    }
	    
	    /// Const zero helicity phase access.
	    
	    const value_type& zero_helicity_phase() const
	    {
		return zero_hel;
	    }

	    /// Maximal helicity readout.

	    int max_helicity() const
	    {
		return pos_hels.size();
	    }

	    /// Nr. of helicity states.
	    
	    size_type n_helicities() const
	    {
		return zero_helicity_state?(2*pos_hels.size()+1):(2*pos_hels.size());
	    }

	    /// Sets the superposition to a pure helicity-n state.

	    helicity_phases<value_t>& set_pure_state(int n)
	    {
		pos_hels.assign(pos_hels.size(),value_type(0,0));
		neg_hels.assign(neg_hels.size(),value_type(0,0));
		if(n>0)
		{
		    zero_hel=value_type(0,0);
		    pos_hels[n-1].real()=(r_value_type)1;
		}
		else if(n==0)
		{
		    zero_hel=value_type(1,0);
		}
		else
		{
		    zero_hel=value_type(0,0);
		    neg_hels[-n-1].real()=(r_value_type)1;
		}
		return *this;
	    }

	    /* Printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		for(typename std::vector<value_type>::size_type i=0;i<neg_hels.size();++i)
		{
		    os<<neg_hels[i];
		}
		if(zero_helicity_state)
		{
		    os<<zero_hel;
		}
		for(typename std::vector<value_type>::size_type i=0;i<pos_hels.size();++i)
		{
		    os<<pos_hels[i];
		}
		return os;
	    }

	private:

	    std::vector<value_type>pos_hels;
	    value_type zero_hel;
	    std::vector<value_type>neg_hels;
    };
    
    /// Helicity phase space template specialisation for continuous-helicity
    ///  model types.
    
    template<class model_t>class helicity_ps<model_t,true>
    {
	public:

	    /* Some useful type definitions: */

	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;
	    
	    /* Wave function pointer type definition: */

	    typedef void(*wave_func)(const value_type&,iterator,const momentum_type*,const r_value_type*);

	    /* Phase space instance: */

	    const base_ps<model_t>* const ps;

	    /* Constructor taking a base ps instance as argument: */
	    
	    helicity_ps(const base_ps<model_t>* ps_):ps(ps_),hel(ps->max_helicity,ps->massive()),n(-(ps->max_helicity)){}

	    /// Returns a reference to the (complex) coefficient that multiplies the i-th helicity wave function.

	    value_type& helicity_phase(int i)
	    {
		return hel.helicity_phase(i);
	    }

	    /// Returns a constant reference to the (complex) coefficient that multiplies the i-th helicity wave function.

	    const value_type& helicity_phase(int i) const
	    {
		return hel.helicity_phase(i);
	    }

	    /* Wave function caller: */

	    iterator fill_wave_function(iterator it)
	    {
		if(ps->is_handlebar())
		{
		    for(size_type mu=0;mu<model_t::dimension;++mu)
		    {
			it[mu]=ps->momentum(mu);
		    }
		    return it;
		}
		const particle<model_t>* phi(ps->particle_type);
		if(ps->massless())
		{
		    if(ps->is_outgoing())
		    {
			for(int i=0;i<ps->max_helicity;++i)
			{
			    (phi->pos_hel_outgoing_massless_states)[i](hel.positive_helicity_phase(i),it,&(ps->momentum()),NULL);
			    (phi->neg_hel_outgoing_massless_states)[i](hel.negative_helicity_phase(i),it,&(ps->momentum()),NULL);
			}
			if(phi->zero_hel_outgoing_massless_state!=NULL)
			{
			    phi->zero_hel_outgoing_massless_state(hel.zero_helicity_phase(),it,&(ps->momentum()),NULL);
			}
			return it;
		    }
		    else
		    {
			for(int i=0;i<ps->max_helicity;++i)
			{
			    (phi->pos_hel_incoming_massless_states)[i](hel.positive_helicity_phase(i),it,&(ps->momentum()),NULL);
			    (phi->neg_hel_incoming_massless_states)[i](hel.negative_helicity_phase(i),it,&(ps->momentum()),NULL);
			}
			if(phi->zero_hel_incoming_massless_state!=NULL)
			{
			    phi->zero_hel_incoming_massless_state(hel.zero_helicity_phase(),it,&(ps->momentum()),NULL);
			}
			return it;
		    }
		}
		else
		{
		    if(ps->is_outgoing())
		    {
			for(int i=0;i<ps->max_helicity;++i)
			{
			    (phi->pos_hel_outgoing_massive_states)[i](hel.positive_helicity_phase(i),it,&(ps->momentum()),phi->mass);
			    (phi->neg_hel_outgoing_massive_states)[i](hel.negative_helicity_phase(i),it,&(ps->momentum()),phi->mass);
			}
			if(phi->zero_hel_outgoing_massive_state!=NULL)
			{
			    phi->zero_hel_outgoing_massive_state(hel.zero_helicity_phase(),it,&(ps->momentum()),phi->mass);
			}
			return it;
		    }
		    else
		    {
			for(int i=0;i<ps->max_helicity;++i)
			{
			    (phi->pos_hel_incoming_massive_states)[i](hel.positive_helicity_phase(i),it,&(ps->momentum()),phi->mass);
			    (phi->neg_hel_incoming_massive_states)[i](hel.negative_helicity_phase(i),it,&(ps->momentum()),phi->mass);
			}
			if(phi->zero_hel_incoming_massive_state!=NULL)
			{
			    phi->zero_hel_incoming_massive_state(hel.zero_helicity_phase(),it,&(ps->momentum()),phi->mass);
			}
			return it;
		    }
		}
		return it;
	    }
	    
	    /* Initialisation of helicity summation mode: */

	    iterator init_sum_wave_function(iterator it)
	    {
		n=-ps->max_helicity;
		hel.set_pure_state(n);
		if(ps->is_handlebar())
		{
		    for(size_type mu=0;mu<model_t::dimension;++mu)
		    {
			it[mu]=ps->momentum(mu);
		    }
		    return it;
		}
		for(size_type i=0;i<ps->spacetime_tensor_size;++i)
		{
		    it[i]=value_type(0,0);
		}
		const particle<model_t>* phi(ps->particle_type);
		if(ps->massless())
		{
		    if(ps->is_outgoing())
		    {
			if(n==0)
			{
			    (phi->zero_hel_outgoing_massless_state)(value_type(1,0),it,&(ps->momentum()),NULL);
			}
			if(n<0)
			{
			    (phi->neg_hel_outgoing_massless_states)[-n-1](value_type(1,0),it,&(ps->momentum()),NULL);
			}
		    }
		    else
		    {
			if(n==0)
			{
			    (phi->zero_hel_incoming_massless_state)(value_type(1,0),it,&(ps->momentum()),NULL);
			}
			if(n<0)
			{
			    (phi->neg_hel_incoming_massless_states)[-n-1](value_type(1,0),it,&(ps->momentum()),NULL);
			}
		    }
		}
		else
		{
		    if(ps->is_outgoing())
		    {
			if(n==0)
			{
			    (phi->zero_hel_outgoing_massive_state)(value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
			if(n<0)
			{
			    (phi->neg_hel_outgoing_massive_states)[-n-1](value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
		    }
		    else
		    {
			if(n==0)
			{
			    (phi->zero_hel_incoming_massive_state)(value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
			if(n<0)
			{
			    (phi->neg_hel_incoming_massive_states)[-n-1](value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
		    }
		}
		return it;
	    }
	    
	    /* Wave function spin sum function: */

	    iterator sum_wave_function(iterator it)
	    {
		++n;
		if(n==0 and !ps->particle_type->has_zero_helicity_state())
		{
		    n=1;
		}
		if(n>ps->max_helicity)
		{
		    n=-n+1;
		}
		hel.set_pure_state(n);
		if(ps->is_handlebar())
		{
		    for(size_type mu=0;mu<model_t::dimension;++mu)
		    {
			it[mu]=ps->momentum(mu);
		    }
		    return it;
		}
		for(size_type i=0;i<ps->spacetime_tensor_size;++i)
		{
		    it[i]=value_type(0,0);
		}
		const particle<model_t>* phi(ps->particle_type);
		if(ps->massless())
		{
		    if(ps->is_outgoing())
		    {
			if(n>0)
			{
			    (phi->pos_hel_outgoing_massless_states)[n-1](value_type(1,0),it,&(ps->momentum()),NULL);
			}
			if(n==0)
			{
			    (phi->zero_hel_outgoing_massless_state)(value_type(1,0),it,&(ps->momentum()),NULL);
			}
			if(n<0)
			{
			    (phi->neg_hel_outgoing_massless_states)[-n-1](value_type(1,0),it,&(ps->momentum()),NULL);
			}
		    }
		    else
		    {
			if(n>0)
			{
			    (phi->pos_hel_incoming_massless_states)[n-1](value_type(1,0),it,&(ps->momentum()),NULL);
			}
			if(n==0)
			{
			    (phi->zero_hel_incoming_massless_state)(value_type(1,0),it,&(ps->momentum()),NULL);
			}
			if(n<0)
			{
			    (phi->neg_hel_incoming_massless_states)[-n-1](value_type(1,0),it,&(ps->momentum()),NULL);
			}
		    }
		}
		else
		{
		    if(ps->is_outgoing())
		    {
			if(n>0)
			{
			    (phi->pos_hel_outgoing_massive_states)[n-1](value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
			if(n==0)
			{
			    (phi->zero_hel_outgoing_massive_state)(value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
			if(n<0)
			{
			    (phi->neg_hel_outgoing_massive_states)[-n-1](value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
		    }
		    else
		    {
			if(n>0)
			{
			    (phi->pos_hel_incoming_massive_states)[n-1](value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
			if(n==0)
			{
			    (phi->zero_hel_incoming_massive_state)(value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
			if(n<0)
			{
			    (phi->neg_hel_incoming_massive_states)[-n-1](value_type(1,0),it,&(ps->momentum()),phi->mass);
			}
		    }
		}
		return it;
	    }

	    /* Minimal helicity output boolean: */

	    bool minimal_helicity() const
	    {
		return (n==-(ps->max_helicity));
	    }

	    /* Sub-amplitude contraction caller: */

	    value_type contract(const_iterator internal,const_iterator external)
	    {
		value_type amp(0,0);
		if(ps->is_outgoing())
		{
		    ps->particle_type->contraction(amp,internal,external);
		    return amp;
		}
		ps->particle_type->conj_contraction(amp,internal,external);
		return amp;
	    }

	    /* Helicity holder address readout: */

	    helicity_phases<r_value_type>& get_helicity_phases()
	    {
		return hel;
	    }

	    /* Printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		hel.print(os);
		return os;
	    }

	protected:

	    /* Helicity holder instance: */

	    helicity_phases<r_value_type> hel;
	    
	    /* Position of the helicity summer in the helicity vector: */
	    
	    int n;
    };

    /// Colour phase space class template. Takes as template parameters the model type
    /// and the compile-time continuous_colours boolean that should be
    /// defined in model_t.

    template<class model_t,bool cont_colours>class colour_ps;
    
    /// Colour phase space template specialisation for discrete-colour model
    /// types.

    template<class model_t>class colour_ps<model_t,false>
    {
	public:

	    /* Some useful type definitions: */

	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;

	    /* Contraction function pointer type definition: */

	    typedef value_type(*contr_func)(const_iterator,const_iterator,const std::vector<size_type>&,size_type);

	    /* Phase space instance: */

	    const base_ps<model_t>* const ps;

	    /* Constructor taking the phase space base object as argument: */
	    
	    colour_ps(const base_ps<model_t>* ps_):ps(ps_),n(0)
	    {
		colours.resize(ps->particle_type->get_colour_rank(),0);
	    }

	    /// Returns a reference to the i-th colour of the particle.

	    size_type& colour(size_type i)
	    {
		CAMGEN_ERROR_IF((i>=colours.size()),"requested colour entry out of range");
		return colours[i];
	    }

	    /// Returns a constant reference to the i-th colour of the particle.

	    const size_type& colour(size_type i) const
	    {
		CAMGEN_ERROR_IF((i>=colours.size()),"requested colour entry out of range");
		return colours[i];
	    }

	    /// Returns a pointer to the full colour vector

	    std::vector<size_type>* get_colours()
	    {
		return &colours;
	    }

	    /* Wave function call preparation (moves the iterator to the position
	     * determined by the colours): */

	    iterator pre_fill_wave_function(iterator it)
	    {
		std::vector<size_type>swapped_cols(colours);
		const particle<model_t>* phi(ps->particle_type);
		phi->swap_colours(ps->is_outgoing(),swapped_cols);
		for(size_type i=0;i<phi->decomposed_colours;++i)
		{
		    it.forward(i,swapped_cols[i]);
		}
		iterator clone=it;
		for(size_type i=phi->decomposed_colours;i<swapped_cols.size();++i)
		{
		    clone.forward(i,swapped_cols[i]);
		}
		return clone;
	    }

	    /* Wave function call evaluation (trivial): */

	    iterator post_fill_wave_function(iterator it)
	    {
		return it;
	    }

	    /* Initialisation of colour summation mode: */

	    iterator init_sum_wave_function(iterator it)
	    {
		n=0;
		colours.assign(colours.size(),0);
		return it;
	    }

	    /* Colour wave function sum utility: */

	    iterator shift_wave_function(iterator it)
	    {
		it+=n*(ps->spacetime_tensor_size);
		++n;
		if(n==ps->colour_tensor_size)
		{
		    iterator it2=it;
		    it2.reset();
		    for(size_type i=0;i<ps->spacetime_tensor_size;++i)
		    {
			*it2=*it;
			*it=value_type(0,0);
			++it;
			++it2;
		    }
		    colours.assign(colours.size(),0);
		    n=0;
		    return it.reset();
		}
		else
		{
		    iterator it2=it+ps->spacetime_tensor_size;
		    for(size_type i=0;i<ps->spacetime_tensor_size;++i)
		    {
			*it2=*it;
			*it=value_type(0,0);
			++it;
			++it2;
		    }
		    typename std::vector<size_type>::const_reverse_iterator it4=ps->particle_type->get_colour_index_ranges().rbegin();
		    for(typename std::vector<size_type>::reverse_iterator it3=colours.rbegin();it3!=colours.rend();++it3)
		    {
			++(*it3);
			(*it3)%=(*it4);
			if(*it3!=0)
			{
			    break;
			}
			++it4;
		    }
		    return it;
		}
	    }

	    /* Contraction function caller: */

	    value_type contract(const_iterator internal,const_iterator external)
	    {
		std::vector<size_type>swapped_cols(colours);
		const particle<model_t>* phi(ps->particle_type);
		phi->swap_colours(ps->is_outgoing(),swapped_cols);
		if(ps->is_outgoing())
		{
		    return phi->contraction(internal,external,swapped_cols,0);
		}
		return phi->conj_contraction(internal,external,swapped_cols,0);
	    }

	    /* Minimal colour output boolean: */

	    bool minimal_colours() const
	    {
		return (n==0);
	    }

	    /* Printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<"( ";
		for(typename std::vector<size_type>::size_type i=0;i<colours.size();++i)
		{
		    os<<colours[i]<<" ";
		}
		os<<")";
		return os;
	    }

	protected:

	    /* Colours of the particle: */

	    std::vector<size_type> colours;

	    /* Colour sum variable: */

	    size_type n;
    };
    
    /// Colour phase space template specialisation for continuous-colour model
    /// types.

    template<class model_t>class colour_ps<model_t,true>
    {
	public:

	    /* Some useful type definitions: */

	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;

	    /* Contraction function pointer type definition: */

	    typedef void(*contr_func)(value_type&,const_iterator,const_iterator);

	    /* Phase space instance: */

	    const base_ps<model_t>* const ps;

	    /* Constructor taking the phase space base pointer as argument: */

	    colour_ps(const base_ps<model_t>* ps_):ps(ps_),n(0)
	    {
		colours.resize(ps->particle_type->get_colour_index_ranges());
		spacetime_part.resize(ps->particle_type->get_spacetime_index_ranges());
	    }

	    /// Returns a reference to the i-th entry in the colour tensor of
	    /// the particle.

	    value_type& colour_entry(size_type i)
	    {
		return colours[i];
	    }

	    /// Returns a constant reference to the i-th entry in the colour tensor of
	    /// the particle.

	    const value_type& colour_entry(size_type i) const
	    {
		return colours[i];
	    }

	    /// Returns a pointer to the colour tensor.

	    tensor_type* get_colours()
	    {
		return &colours;
	    }

	    /// Returns a reference to the (i0,...,in)-component in the colour tensor of
	    /// the particle.

	    value_type& colour_coeff(size_type i0,...)
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<colours.rank();++j)
		{
		    n+=(i*colours.block_size(j));
		    i=va_arg(I,size_type);
		}
		va_end(I);
		return colours[n];
	    }

	    /// Returns a constant reference to the (i0,...,in)-component in the colour tensor of
	    /// the particle.

	    const value_type& colour_coeff(size_type i0,...) const
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<colours.rank();++j)
		{
		    n+=(i*colours.block_size(j));
		    i=va_arg(I,size_type);
		}
		va_end(I);
		return colours[n];
	    }

	    /* Wave function call preparation (replace iterator by spacetime
	     * tensor's begin() iterator: */

	    iterator pre_fill_wave_function(iterator it)
	    {
		return it;
	    }

	    /* Wave function call evaluation (multiply all spacetime and colour
	     * components): */

	    iterator post_fill_wave_function(iterator it)
	    {
		for(size_type i=0;i<spacetime_part.size();++i)
		{
		    spacetime_part[i]=it[i];
		}
		for(size_type i=0;i<colours.size();++i)
		{
		    for(size_type j=0;j<spacetime_part.size();++j)
		    {
			(*it)=colours[i]*spacetime_part[j];
			++it;
		    }
		}
		return it.reset();
	    }

	    /* Initialisation of colour summation mode: */

	    iterator init_sum_wave_function(iterator it)
	    {
		n=0;
		colours[0]=value_type(1,0);
		colours.reset();
		return it.reset();
	    }

	    /* Colour wave function sum utility: */

	    iterator shift_wave_function(iterator it)
	    {
		it+=n*spacetime_part.size();
		colours[n]=value_type(0,0);
		++n;
		if(n==colours.size())
		{
		    iterator it2=it;
		    it2.reset();
		    for(size_type i=0;i<spacetime_part.size();++i)
		    {
			*it2=*it;
			*it=value_type(0,0);
			++it;
			++it2;
		    }
		    n=0;
		    colours[n]=value_type(1,0);
		    return it.reset();
		}
		else
		{
		    colours[n]=value_type(1,0);
		    iterator it2=it+spacetime_part.size();
		    for(size_type i=0;i<spacetime_part.size();++i)
		    {
			*it2=*it;
			*it=value_type(0,0);
			++it;
			++it2;
		    }
		    return it;
		}
	    }

	    /* Wave function contraction function call: */

	    value_type contract(const_iterator internal,const_iterator external)
	    {
		value_type amp(0,0);
		if(ps->is_outgoing())
		{
		    ps->particle_type->contraction(amp,internal,external);
		    return amp;
		}
		ps->particle_type->conj_contraction(amp,internal,external);
		return amp;
	    }

	    /* Minimal colour output boolean: */

	    bool minimal_colours() const
	    {
		return (n==0);
	    }

	    /* Printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<colours;
		return os;
	    }

	protected:

	    /* Colour tensor: */

	    tensor_type colours;

	    /* Spacetime factor of the external wave function: */

	    tensor_type spacetime_part;

	    /* Summation mode counter: */

	    size_type n;
    };

    /// Particle phase space template specialisation for colourful model types.

    template<class model_t,std::size_t dim>class particle_ps<model_t,dim,true>: public base_ps<model_t>,public helicity_ps<model_t,model_t::continuous_helicities>,public colour_ps<model_t,model_t::continuous_colours>
    {
	public:

	    /* Base class type definitions: */

	    typedef base_ps<model_t> momentum_base_type;
	    typedef helicity_ps<model_t,model_t::continuous_helicities> helicity_base_type;
	    typedef colour_ps<model_t,model_t::continuous_colours> colour_base_type;
	
	    /* Some useful type definitions: */

	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;
	    
	    /* Wave function and contraction function pointer type definitions:
	     * */
	    
	    typedef typename helicity_base_type::wave_func wave_func;
	    typedef typename colour_base_type::contr_func contr_func;

	    /* Constructor taking the particle type as argument: */
	    
	    particle_ps(const particle<model_t>* phi):momentum_base_type(phi),helicity_base_type(static_cast<const momentum_base_type*>(this)),colour_base_type(static_cast<const momentum_base_type*>(this)){}

	    /* Constructor taking the particle type and momentum direction as
	     * arguments: */

	    particle_ps(const particle<model_t>* phi,bool out):momentum_base_type(phi,out),helicity_base_type(static_cast<const momentum_base_type*>(this)),colour_base_type(static_cast<const momentum_base_type*>(this)){}
	    
	    /* Wave function call. The colour base class performs the
	     * preparation and evaluation, whereas the helicity base class the
	     * actual filling of the wave function: */

	    iterator fill_wave_function(iterator it)
	    {
		return this->post_fill_wave_function(this->helicity_base_type::fill_wave_function(this->pre_fill_wave_function(it)));
	    }

	    /* Initialisation of the wave function summation mode: */

	    iterator init_dof_sum(iterator it,bool hsum,bool csum)
	    {
		if((!hsum) and (!csum))
		{
		    return this->post_fill_wave_function(this->helicity_base_type::fill_wave_function(this->pre_fill_wave_function(it)));
		}
		if(hsum and !csum)
		{
		    return this->post_fill_wave_function(this->helicity_base_type::init_sum_wave_function(this->pre_fill_wave_function(it)));
		}
		if(csum and !hsum)
		{
		    return this->helicity_base_type::fill_wave_function(this->colour_base_type::init_sum_wave_function(it));
		}
		return this->helicity_base_type::init_sum_wave_function(this->colour_base_type::init_sum_wave_function(it));
	    }
	    
	    /* Wave function call in spin/colour summation mode. The second
	     * argument denotes whether we are summing over helicities, the
	     * third argument whether we are summing over colours. */

	    iterator sum_wave_function(iterator it,bool hsum,bool csum)
	    {
		if((!hsum) and (!csum))
		{
		    return this->post_fill_wave_function(this->helicity_base_type::fill_wave_function(this->pre_fill_wave_function(it)));
		}
		if(hsum and !csum)
		{
		    return this->post_fill_wave_function(this->helicity_base_type::sum_wave_function(this->pre_fill_wave_function(it)));
		}
		if(csum and !hsum)
		{
		    return this->shift_wave_function(it);
		}
		iterator it2=this->shift_wave_function(it);
		if(this->minimal_colours())
		{
		    return this->helicity_base_type::sum_wave_function(it2);
		}
		return it2;
	    }

	    /* Boolean function denoting whether the summation loop has ended.
	     * */

	    bool minimal_dof(bool hsum,bool csum) const
	    {
		if(hsum and !csum)
		{
		    return this->minimal_helicity();
		}
		if(csum and !hsum)
		{
		    return this->minimal_colours();
		}
		if(hsum and csum)
		{
		    return (this->minimal_helicity() and this->minimal_colours());
		}
		return true;
	    }

	    /* Contraction function call, inherited from the colour base class:
	     * */

	    value_type contract(const_iterator internal,const_iterator external)
	    {
		return colour_base_type::contract(internal,external);
	    }

	    /* Swapping method, interchanging all phase space variables with
	     * another particle phase space: */

	    void swap(particle_ps<model_t,dim,true>& other)
	    {
		std::swap(this->p,other.p);
		std::swap(this->hel,other.hel);
		std::swap(this->colours,other.colours);
	    }

	    /* Printing method: */

	    std::ostream& print_dofs(std::ostream& os) const
	    {
		os<<"[";
		this->helicity_base_type::print(os);
		os<<";";
		this->colour_base_type::print(os);
		os<<"]";
		return os;
	    }
    };

    /// Particle phase space template specialisation for colourless model types.

    template<class model_t,std::size_t dim>class particle_ps<model_t,dim,false>: public base_ps<model_t>,public helicity_ps<model_t,model_t::continuous_helicities>
    {
	public:

	    /* Base class type definitions: */

	    typedef base_ps<model_t> momentum_base_type;
	    typedef helicity_ps<model_t,model_t::continuous_helicities> helicity_base_type;
	    typedef unused colour_base_type;
	
	    /* Wave function and contraction function pointer type definitions:
	     * */
	    
	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;
	    
	    /* Wave function and contraction function pointer type definitions:
	     * */
	    
	    typedef typename helicity_base_type::wave_func wave_func;
	    typedef typename colour_ps<model_t,true>::contr_func contr_func;

	    /* Constructor taking the particle type as argument: */
	    
	    particle_ps(const particle<model_t>* phi):momentum_base_type(phi),helicity_base_type(static_cast<const momentum_base_type*>(this)){}

	    /* Constructor taking the particle type and momentum direction as
	     * arguments: */

	    particle_ps(const particle<model_t>* phi,bool out):momentum_base_type(phi,out),helicity_base_type(static_cast<const momentum_base_type*>(this)){}
	    
	    /* Wave function calling method. Invokes helicity base class wave
	     * function fill: */

	    iterator fill_wave_function(iterator it)
	    {
		return this->helicity_base_type::fill_wave_function(it);
	    }

	    /* Initialisation of the wave function summation mode: */

	    iterator init_dof_sum(iterator it,bool hsum,bool csum)
	    {
		if(hsum)
		{
		    return this->helicity_base_type::init_sum_wave_function(it);
		}
		return this->helicity_base_type::fill_wave_function(it);
	    }
	    
	    /* Wave function call in spin/colour summation mode. The second
	     * argument denotes whether we are summing over helicities, the
	     * third argument whether we are summing over colours. */

	    iterator sum_wave_function(iterator it,bool hsum,bool csum)
	    {
		if(hsum)
		{
		    return this->helicity_base_type::sum_wave_function(it);
		}
		return this->helicity_base_type::fill_wave_function(it);
	    }

	    /* Boolean function denoting whether the summation loop has ended.
	     * */

	    bool minimal_dof(bool hsum,bool csum) const
	    {
		if(hsum)
		{
		    return this->minimal_helicity();
		}
		return true;
	    }

	    /* Swapping method, interchanging all phase space variables with 
	     * another particle phase space: */

	    void swap(particle_ps<model_t,dim,false>& other)
	    {
		std::swap(this->p,other.p);
		std::swap(this->hel,other.hel);
	    }

	    /* Printing method: */

	    std::ostream& print_dofs(std::ostream& os) const
	    {
		os<<"[";
		this->helicity_base_type::print(os);
		os<<";0]";
		return os;
	    }
    };

    /// Particle phase space template specialisation for colourful zero-dimensional model types.

    template<class model_t>class particle_ps<model_t,0,true>: public base_ps<model_t>,public colour_ps<model_t,model_t::continuous_colours>
    {
	public:

	    /* Base class type definitions: */

	    typedef base_ps<model_t> momentum_base_type;
	    typedef colour_ps<model_t,model_t::continuous_colours> colour_base_type;
	    typedef unused helicity_base_type;
	
	    /* Some useful type definitions: */
	    
	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;
	    
	    /* Wave function and contraction function pointer type definitions:
	     * */

	    typedef typename helicity_ps<model_t,true>::wave_func wave_func;
	    typedef typename colour_base_type::contr_func contr_func;

	    /* Constructor taking the particle type as argument */

	    particle_ps(const particle<model_t>* phi):momentum_base_type(phi),colour_base_type(static_cast<const momentum_base_type*>(this)){}

	    /* Constructor taking the particle type and momentum direction as
	     * arguments: */

	    particle_ps(const particle<model_t>* phi,bool out):momentum_base_type(phi,out),colour_base_type(static_cast<const momentum_base_type*>(this)){}
	    
	    /* Wave function calling method. Invokes the wave function
	     * preparation in the base class, sets the iterator value to one,
	     * and calls the post_fill_wave_function method: */

	    iterator fill_wave_function(iterator it)
	    {
		*(this->pre_fill_wave_function(it))=value_type(1,0);
		return this->post_fill_wave_function(it);
	    }

	    /* Initialisation of the wave function summation mode: */

	    iterator init_dof_sum(iterator it,bool hsum,bool csum)
	    {
		if(csum)
		{
		    *it=value_type(1,0);
		    return it;
		}
		*(this->pre_fill_wave_function(it))=value_type(1,0);
		return this->post_fill_wave_function(it);
	    }
	    
	    /* Wave function call in spin/colour summation mode. The second
	     * argument denotes whether we are summing over helicities, the
	     * third argument whether we are summing over colours. */

	    iterator sum_wave_function(iterator it,bool hsum,bool csum)
	    {
		if(csum)
		{
		    if(this->minimal_colours())
		    {
			*it=value_type(1,0);
			return it;
		    }
		    else
		    {
			return this->shift_wave_function(it);
		    }
		}
		*(this->pre_fill_wave_function(it))=value_type(1,0);
		return this->post_fill_wave_function(it);
	    }

	    /* Boolean function denoting whether the summation loop has ended.
	     * */

	    bool minimal_dof(bool hsum,bool csum) const
	    {
		if(csum)
		{
		    return this->minimal_colours();
		}
		return true;
	    }

	    /* Contraction method is derived from the base class...*/

	    /* Swapping method, interchanging all phase space variables with 
	     * another particle phase space: */

	    void swap(particle_ps<model_t,0,true>& other)
	    {
		std::swap(this->p,other.p);
		std::swap(this->colours,other.colours);
	    }

	    /* Printing method: */

	    std::ostream& print_dofs(std::ostream& os) const
	    {
		os<<"[0;";
		this->colour_base_type::print(os);
		os<<"]";
		return os;
	    }
    };

    /// Particle phase space template specialisation for colourless zero-dimensional model types.
    
    template<class model_t>class particle_ps<model_t,0,false>: public base_ps<model_t>
    {
	public:

	    /* Base class type definitions: */

	    typedef base_ps<model_t> momentum_base_type;
	    typedef unused colour_base_type;
	    typedef unused helicity_base_type;
	    
	    /* Some useful type definitions: */
	    
	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef typename tensor_type::size_type size_type;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;
	    
	    /* Wave function and contraction function pointer type definitions:
	     * */

	    typedef typename helicity_ps<model_t,true>::wave_func wave_func;
	    typedef typename colour_ps<model_t,true>::contr_func contr_func;

	    /* Constructor taking the particle type as argument: */

	    particle_ps(const particle<model_t>* phi):momentum_base_type(phi){}

	    /* Constructor taking the particle type and momentum direction as
	     * arguments: */

	    particle_ps(const particle<model_t>* phi,bool out):momentum_base_type(phi,out){}

	    /* Number of scalar particle helicities: */

	    std::size_t helicities() const
	    {
		return 1;
	    }

	    /* Wave function for scalar singlets: */
	    
	    iterator fill_wave_function(iterator it)
	    {
		*it=value_type(1,0);
		return it;
	    }

	    /* Initialisation of the wave function summation mode: */

	    iterator init_dof_sum(iterator it,bool hsum,bool csum)
	    {
		*it=value_type(1,0);
		return it;
	    }
	    
	    /* Wave function call in spin/colour summation mode. The second
	     * argument denotes whether we are summing over helicities, the
	     * third argument whether we are summing over colours. */

	    iterator sum_wave_function(iterator it,bool hsum,bool csum)
	    {
		*it=value_type(1,0);
		return it;
	    }

	    /* Boolean function denoting whether the summation loop has ended.
	     * */

	    bool minimal_dof(bool hsum,bool csum) const
	    {
		return true;
	    }

	    /* Scalar singlet contraction method: */

	    value_type contract(const_iterator internal,const_iterator external)
	    {
		return (*internal)*(*external);
	    }

	    /* Swapping method is trivial... */

	    void swap(particle_ps<model_t,0,false>& other){}

	    /* Printing method: */

	    std::ostream& print_dofs(std::ostream& os) const
	    {
		os<<"[0]";
		return os;
	    }
    };
}

#endif /*CAMGEN_PHASE_SPACE_H_*/

