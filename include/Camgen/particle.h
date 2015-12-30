//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PARTICLE_H_
#define CAMGEN_PARTICLE_H_

#include <Camgen/spin.h>
#include <Camgen/Dirac_dim.h>
#include <Camgen/comp_contr.h>
#include <Camgen/charge_conj.h>
#include <Camgen/def_args.h>
#include <Camgen/part_painter.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Particle class definition. The particle class has only private constructors   *
 * and public static functions returning pointers to particle objects. In        *
 * Camgen, particles are only allocated in the model construction, from that    *
 * point on only addresses are copied as also the copy constructor and           *
 * assignment operator are private. The class holds the name, spin, mass, width  *
 * etc. and function pointers to the wave function filling algorithms,           *
 * propagators, contraction function and possibly their charge conjugated        *
 * counterparts.                                                                 *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition of the particle class template, where model_t is denotes the
     * model type the particle resides in: */

    template<class model_t>class particle
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	    /* Definitions of the propagator, propagator refresher, wave
	     * function and contraction function pointer types: */

	    typedef typename get_basic_types<model_t>::prop_func prop_func;
	    typedef typename get_basic_types<model_t>::prop_refresher prop_refresher;
	    typedef typename get_basic_types<model_t>::wave_func wave_func;
	    typedef typename get_basic_types<model_t>::contr_func contr_func;

	    /* Definition of the phase space type: */

	    typedef typename get_basic_types<model_t>::phase_space_type phase_space_type;

	    /* Declaring the vertex, flavour comparison classes, model wrapper
	     * and phase space class as friends: */

	    friend class vertex<model_t>;
	    friend class flavour_comp< particle<model_t> >;
	    friend class flavour_eq< particle<model_t> >;
	    friend class flavour_neq< particle<model_t> >;
	    friend class flavour_select< particle<model_t> >;
	    friend class flavourvec_comp< particle<model_t> >;
	    friend class model_wrapper<model_t>;
	    friend class particle_ps<model_t,model_t::dimension,model_t::coloured>;
	    friend class base_ps<model_t>;
	    friend class helicity_ps<model_t,false>;
	    friend class helicity_ps<model_t,true>;
	    friend class colour_ps<model_t,false>;
	    friend class colour_ps<model_t,true>;
	    template<class rep_t,std::size_t N>friend class particle_painter;

	    /* Static checking utility whether the given propagator and
	     * wave function types correspond to the spin argument: */

	    template<class particle_t>static bool check_spin(const std::string& name)
	    {
		/* Creating a tensor rank vector from the given spin: */

		std::vector<size_type>I1;
		for(size_type i=0;i<particle_t::particle_spin.Lorentz_rank();++i)
		{
		    I1.push_back(model_t::dimension);
		}
		for(size_type i=0;i<particle_t::particle_spin.Dirac_rank();++i)
		{
		    I1.push_back(Dirac_dim<model_t::dimension>::value);
		}

		/* Checking whether the wave function correponds to the given
		 * spin: */

		std::vector<size_type>I2;
		if(particle_t::get_index_ranges(I2) != I1)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid wave function assignment for particle "<<name<<endlog;
		    return false;
		}

		/* Checking whether the propagator corresponds to the given
		 * spin: */

		if(particle_t::propagator_type::get_index_ranges(I2) != I1)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid propagator assignment for particle "<<name<<endlog;
		    return false;
		}

		/* Checking whether the contraction routine corresponds to the given
		 * spin: */

		if(evaluate<typename particle_t::contraction_type>::fill_rank_vector(I2) != I1)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid wave function contraction assignment for particle "<<name<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Static creator of a pointer to a colour singlet of type
	     * particle_t with mass *m and width *w: */

	    template<class particle_t>static particle<model_t>* create(const std::string& name,const r_value_type* m=NULL,const r_value_type* w=NULL,int n=0)
	    {
		if(check_spin<particle_t>(name))
		{
		    particle<model_t>* phi=new particle<model_t>(name,particle_t::particle_spin,m,w,n);
		    phi->set_particle_type<particle_t>();
		    return phi;
		}
		else
		{
		    return NULL;
		}
	    }

	    /* Static creator of a pointer to a colour singlet
	     * particle/antiparticle pair of type
	     * particle_t with mass *m and width *w: */

	    template<class particle_t>static std::pair<particle<model_t>*,particle<model_t>*> create(const std::string& name0,const std::string& name1,const r_value_type* m=NULL,const r_value_type* w=NULL,int n=0)
	    {
		particle<model_t>* phi0=create<particle_t>(name0,m,w,n);
		particle<model_t>* phi1=create<typename particle_t::anti_particle_type>(name1,m,w,-n);
		if(phi0!=NULL and phi1!=NULL)
		{
		    phi0->anti_particle=phi1;
		    phi1->anti_particle=phi0;
		    if(phi0!=phi1)
		    {
			phi0->Majorana_tag=false;
			phi0->conj_contraction=phi1->contraction;
			phi1->Majorana_tag=false;
			phi1->conj_contraction=phi0->contraction;
		    }
		}
		return std::pair<particle<model_t>*,particle<model_t>*>(phi0,phi1);
	    }

	    /* Static creator of a pointer to a colour multiplet in
	     * representation rep_t of type particle_t with mass *m and width
	     * *w: */

	    template<class particle_t,class rep_t>static particle<model_t>* create(const std::string& name,const r_value_type* m=NULL,const r_value_type* w=NULL,int n=0)
	    {
		typedef compose_contraction<colour_tensor::d<rep_t,0,1>,typename particle_t::contraction_type> adj_contr_type;
		typedef typename model_t::colour_treatment::template make_contraction<adj_contr_type>::type contr_type;
		
		if(check_spin<particle_t>(name))
		{
		    particle<model_t>* phi=new particle<model_t>(name,particle_t::particle_spin,m,w,n);
		    phi->set_particle_type<particle_t>();
		    particle_painter<rep_t,rep_t::multiplicity>::template apply< particle<model_t> >(phi);
		    phi->set_contraction<contr_type>();
		    return phi;
		}
		else
		{
		    return NULL;
		}
	    }

	    /* Static creator of a pointer to a colour multiplet
	     * particle/antiparticle pair of type
	     * particle_t with mass *m and width *w: */

	    template<class particle_t,class rep_t>static std::pair<particle<model_t>*,particle<model_t>*> create(const std::string& name0,const std::string& name1,const r_value_type* m=NULL,const r_value_type* w=NULL,int n=0)
	    {
		particle<model_t>* phi0=create<particle_t,rep_t>(name0,m,w,n);
		particle<model_t>* phi1=create<typename particle_t::anti_particle_type,rep_t>(name1,m,w,-n);
		if(phi0!=NULL and phi1!=NULL)
		{
		    phi0->anti_particle=phi1;
		    phi1->anti_particle=phi0;
		    if(phi0!=phi1)
		    {
			phi0->Majorana_tag=false;
			phi0->conj_contraction=phi1->contraction;
			phi1->Majorana_tag=false;
			phi1->conj_contraction=phi0->contraction;
			for(size_type i=0;i<phi1->colour_numbers.size();++i)
			{
			    (phi1->colour_numbers[i])*=(-1);
			}
		    }
		}
		return std::pair<particle<model_t>*,particle<model_t>*>(phi0,phi1);
	    }

	    /* Static creator of a pointer to an auxiliary colour singlet of
	     * spin s (trivial propagator and wave function types): */

	    static particle<model_t>* create_auxiliary(const std::string& name,const spin& s=spin::integer(0))
	    {
		particle<model_t>* phi=new particle<model_t>(name,s,NULL,NULL,0);
		phi->auxiliary=true;
		phi->propagator=&trivial_prop;
		phi->anti_propagator=&trivial_prop;
		phi->refresh_prop=&trivial_refresh;
		phi->refresh_anti_prop=&trivial_refresh;
		return phi;
	    }

	    /* Static creator of pair of particle/antiparticle pointers to
	     * auxiliary colour singlets of spin s (trivial propagator and wave
	     * function types): */

	    static std::pair<particle<model_t>*,particle<model_t>*> create_auxiliary(const std::string& name0,const std::string& name1,const spin& s=spin::integer(0))
	    {
		particle<model_t>* phi0=create_auxiliary(name0,s);
		particle<model_t>* phi1=create_auxiliary(name1,s);
		phi0->anti_particle=phi1;
		phi1->anti_particle=phi0;
		return std::pair<particle<model_t>*,particle<model_t>*>(phi0,phi1);
	    }

	    /* Static creator of a pointer to an auxiliary colour multiplet of
	     * spin s in representation rep_t (trivial propagator and wave
	     * function types): */

	    template<class rep_t>static particle<model_t>* create_auxiliary(const std::string& name,const spin& s=spin::integer(0))
	    {
		particle<model_t>* phi=new particle<model_t>(name,s,NULL,NULL,0);
		phi->auxiliary=true;
		phi->propagator=&trivial_prop;
		phi->anti_propagator=&trivial_prop;
		phi->refresh_prop=&trivial_refresh;
		phi->refresh_anti_prop=&trivial_refresh;
		particle_painter<rep_t>::template apply< particle<model_t> >(phi);
		return phi;
	    }

	    /* Static creator of a pair of particle/antiparticle pointers to an
	     * auxiliary colour multiplet of spin s in representation rep_t
	     * (trivial propagator and wave function types): */

	    template<class rep_t>static std::pair<particle<model_t>*,particle<model_t>*> create_auxiliary(const std::string& name0,const std::string& name1,const spin& s=spin::integer(0))
	    {
		particle<model_t>* phi0=create_auxiliary<rep_t>(name0,s);
		particle<model_t>* phi1=create_auxiliary<rep_t>(name1,s);
		phi0->anti_particle=phi1;
		phi1->anti_particle=phi0;
		if(phi1!=phi0)
		{
		    for(size_type i=0;i<phi1->colour_numbers.size();++i)
		    {
			(phi1->colour_numbers[i])*=(-1);
		    }
		}
		return std::pair<particle<model_t>*,particle<model_t>*>(phi0,phi1);
	    }

	    /* Particle name readout: */

	    const std::string& get_name() const
	    {
		return name;
	    }

	    /* Camgen's internal flavour int readout: */

	    int get_flavour() const
	    {
		return flavour;
	    }

	    /* Mass value readout: */

	    r_value_type get_mass() const
	    {
		if(mass != NULL)
		{
		    return *mass;
		}
		else
		{
		    return 0;
		}
	    }

	    /* Function setting the mass pointer (default=no entry=massless) */

	    void set_mass(const r_value_type* m=NULL)
	    {
		mass=m;
		anti_particle->mass=m;
	    }

	    /* Mass address readout: */

	    const r_value_type* get_mass_address() const
	    {
		return mass;
	    }

	    /* Decay width value readout: */

	    r_value_type get_width() const
	    {
		if(width != NULL)
		{
		    return *width;
		}
		else
		{
		    return 0;
		}
	    }

	    /* Function setting the decay width pointer (default=no entry=zero
	     * width) */

	    void set_width(const r_value_type* w=NULL)
	    {
		width=w;
		anti_particle->width=w;
	    }

	    /* Width address readout: */

	    const r_value_type* get_width_address() const
	    {
		return width;
	    }

	    /* Particle data group id readout: */

	    int get_pdg_id() const
	    {
		return pdg_id;
	    }

	    /* Tensor rank readout: */

	    size_type get_rank() const
	    {
		return index_ranges.size();
	    }

	    /* Index range readout: */

	    size_type get_index_range(const int i) const
	    {
		return index_ranges[i];
	    }

	    /* Index ranges readout: */

	    const std::vector<size_type>& get_index_ranges() const
	    {
		return index_ranges;
	    }

	    /* Spacetime tensor rank readout: */

	    size_type get_spacetime_rank() const
	    {
		return spacetime_index_ranges.size();
	    }

	    /* Spacetime index ramge readout: */

	    size_type get_spacetime_index_range(const int i) const
	    {
		return spacetime_index_ranges[i];
	    }

	    /* Spacetime index ranges readout: */

	    const std::vector<size_type>& get_spacetime_index_ranges() const
	    {
		return spacetime_index_ranges;
	    }

	    /* Spacetime tensor size: */

	    size_type get_spacetime_tensor_size() const
	    {
		size_type n(1);
		for(size_type i=0;i<spacetime_index_ranges.size();++i)
		{
		    n*=(spacetime_index_ranges[i]);
		}
		return n;
	    }

	    /* Colour rank readout: */

	    size_type get_colour_rank() const
	    {
		return colour_index_ranges.size();
	    }

	    /* Colour index range readout: */

	    size_type get_colour_index_range(const int i) const
	    {
		return colour_index_ranges[i];
	    }

	    /* Colour index ranges readout: */

	    const std::vector<size_type>& get_colour_index_ranges() const
	    {
		return colour_index_ranges;
	    }

	    /* Colour tensor size: */

	    size_type get_colour_tensor_size() const
	    {
		size_type n(1);
		for(size_type i=0;i<colour_index_ranges.size();++i)
		{
		    n*=(colour_index_ranges[i]);
		}
		return n;
	    }

	    /* Function allocating a tensor to the appropriate shape to be a
	     * subamplitude: */

	    tensor_type& make_amplitude(tensor_type& t) const
	    {
		return t.resize(index_ranges);
	    }

	    /* Particle spin readout: */

	    spin get_spin() const
	    {
		return S;
	    }

	    /* Function computing the degrees of freedom (helicities) of the
	     * particle: */

	    size_type spin_dof() const
	    {
		if(fermion_nr==0)
		{
		    if(mass==NULL and S!=0)
		    {
			return S.twice();
		    }
		    return S.twice()+1;
		}
		else
		{
		    return S.twice()+1;
		}
	    }

	    /* Returns the physical colour degrees of freedom: */

	    size_type colour_dof() const
	    {
		return col_dim;
	    }

	    /* Fermionic property readout: */

	    bool is_fermion() const
	    {
		return (fermion_nr != 0);
	    }

	    /* Majorana property readout: */

	    bool is_Majorana() const
	    {
		return Majorana_tag;
	    }

	    /* Dirac property readout: */

	    bool is_Dirac() const
	    {
		return (fermion_nr == 1);
	    }

	    /* anti-Dirac property readout: */

	    bool is_antiDirac() const
	    {
		return (fermion_nr == -1);
	    }

	    /* Bosonic property readout: */

	    bool is_boson() const
	    {
		return (fermion_nr == 0);
	    }

	    /* Auxiliary-field property readout: */

	    bool is_auxiliary() const
	    {
		return auxiliary;
	    }

	    /* Massless property readout: */

	    bool is_massless() const
	    {
		return (mass==NULL);
	    }

	    /* Massive property readout: */

	    bool is_massive() const
	    {
		return (mass!=NULL);
	    }

	    /* Request whether the particle has a zero-helicity state: */

	    bool has_zero_helicity_state() const
	    {
		return (mass==NULL)?(zero_hel_incoming_massless_state!=NULL):(zero_hel_incoming_massive_state!=NULL);
	    }

	    /* Comparison operators (determined by the flavour): */

	    bool operator == (const particle<model_t>& other) const
	    {
		return (other.flavour==flavour);
	    }
	    bool operator != (const particle<model_t>& other) const
	    {
		return (other.flavour!=flavour);
	    }
	    bool operator < (const particle<model_t>& other) const
	    {
		return (flavour < other.flavour);
	    }

	    /* Function returning the address of the anti-particle: */

	    const particle<model_t>* get_anti_particle() const
	    {
		return anti_particle;
	    }

	    /* Fermion number readout: */

	    int get_fermion_number() const
	    {
		return fermion_nr;
	    }

	    /* Computing the denominator of a propagator: */

	    void refresh_propagator(const momentum_type* p) const
	    {
		if(coupled)
		{
		    refresh_prop(p,mass,width);
		}
	    }

	    /* Propagator caller with limiting iterator: */

	    void propagate(iterator first,iterator last) const
	    {
		propagator(first,last);
	    }

	    /* Propagator caller for a spacetime block: */

	    void propagate(iterator first) const
	    {
		propagator(first,first+block_size);
	    }

	    /* Colour type function (in the case of colour-flow it yields -1 for
	     * anti-quarks and the second index of gluons): */

	    int colour_type(size_type i) const
	    {
		if(i<colour_numbers.size())
		{
		    return colour_numbers[i];
		}
		else
		{
		    return 0;
		}
	    }

	    /* Decomposed colour number readout: */

	    size_type get_decomposed_colours() const
	    {
		return decomposed_colours;
	    }

	    /* Colour-anticolour swapping function: */

	    void swap_colours(bool q,std::vector<size_type>& cols) const
	    {
		if(q)
		{
		    for(size_type i=0;i<swapped_colours.size();++i)
		    {
			std::swap(cols[swapped_colours[i]],cols[swapped_anti_colours[i]]);
		    }
		}
	    }

	    /* Spacetime block iterator resetter: */

	    void reset(iterator it) const
	    {
		for(size_type i=0;i<block_size;++i)
		{
		    *it=0;
		    ++it;
		}
	    }

	    /* Function coupling the particle to the model it resides in: */

	    void couple()
	    {
		coupled=true;
		anti_particle->coupled=true;
	    }

	    /* Function decoupling the particle from the model it resides in: */

	    void decouple()
	    {
		coupled=false;
		anti_particle->coupled=false;
	    }

	    /* Request whether the particle participates in tree evaluations: */

	    bool is_coupled() const
	    {
		return coupled;
	    }
	
	protected:

	    /* Name of the particle: */

	    std::string name;
	    
	    /* Camgen's internal flavour number: */
	    
	    int flavour;

	    /* Particle spin: */

	    spin S;

	    /* Particle colour degrees of freedom: */

	    size_type col_dim;
	    
	    /* Particle mass (=NULL when massless): */
	    
	    const r_value_type* mass;

	    /* Particle decay width (=NULL when widthless): */

	    const r_value_type* width;

	    /* Particle fermion number (0=boson,1=fermion,-1=anti-fermion): */

	    int fermion_nr;

	    /* Auxiliary-field property: */

	    bool auxiliary;
	    
	    /* Address of the corresponding anti-particle: */
	    
	    particle<model_t>* anti_particle;

	    /* Majorana property tag: */

	    bool Majorana_tag;

	    /* Undecomposed block size: */

	    size_type block_size;

	    /* Number of decomposed colour indices: */

	    size_type decomposed_colours;

	    /* Colour indices to be swapped when the momentum is reversed: */

	    std::vector<size_type>swapped_colours;
	    std::vector<size_type>swapped_anti_colours;

	    /* Lorentz and Dirac index ranges: */

	    std::vector<size_type>spacetime_index_ranges;

	    /* Colour index ranges: */

	    std::vector<size_type>colour_index_ranges;
	    
	    /* Colour index properties (1=colour, -1=anticolour): */
	    
	    std::vector<int>colour_numbers;

	    /* All index ranges: */

	    std::vector<size_type>index_ranges;

	    /* Propagator function address: */

	    prop_func propagator;
	    
	    /* Propagator denominator function address: */
	    
	    prop_refresher refresh_prop;

	    /* Incoming massless wave function addresses: */

	    std::vector<wave_func>pos_hel_incoming_massless_states;
	    wave_func zero_hel_incoming_massless_state;
	    std::vector<wave_func>neg_hel_incoming_massless_states;

	    /* Outgoing massless wave function addresses: */

	    std::vector<wave_func>pos_hel_outgoing_massless_states;
	    wave_func zero_hel_outgoing_massless_state;
	    std::vector<wave_func>neg_hel_outgoing_massless_states;

	    /* Incoming massive wave function addresses: */

	    std::vector<wave_func>pos_hel_incoming_massive_states;
	    wave_func zero_hel_incoming_massive_state;
	    std::vector<wave_func>neg_hel_incoming_massive_states;

	    /* Outgoing massive wave function addresses: */

	    std::vector<wave_func>pos_hel_outgoing_massive_states;
	    wave_func zero_hel_outgoing_massive_state;
	    std::vector<wave_func>neg_hel_outgoing_massive_states;

	    /* Anti-particle propagator address: */

	    prop_func anti_propagator;
	    
	    /* Anti-particle propagator denominator address: */
	    
	    prop_refresher refresh_anti_prop;

	    /* Contraction function address: */

	    contr_func contraction;
	    
	    /* Charge conjugate contraction function address: */
	    
	    contr_func conj_contraction;

	    /* Static flavour number, increases every time a particle in the
	     * model is created. Hence the order of particle creation determines
	     * the internal flavour number: */

	    static int flavour_nr;

	    /* Particle data group id of the particle: */

	    int pdg_id;

	    /* Flag representing whether the particle is coupled (participates
	     * in tree evaluations) or not: */

	    bool coupled;

	    /* Constructor. By default, it creates non-auxiliary self-conjugate
	     * particles (in the case of fermions, Majorana particles). */

	    particle(const std::string& str,const spin& s,const r_value_type* m,const r_value_type* w,int n=0):name(str),flavour(flavour_nr),S(s),col_dim(1),mass(m),width(w),auxiliary(false),anti_particle(this),block_size(1),decomposed_colours(0),contraction(NULL),conj_contraction(NULL),pdg_id(n),coupled(true)
	    {
		/* Filling index ranges: */

		if(model_t::dimension > 0)
		{
		    for(size_type i=0;i<s.Lorentz_rank();++i)
		    {
			spacetime_index_ranges.push_back(model_t::dimension);
			block_size*=(model_t::dimension);
		    }
		}
		for(size_type i=0;i<s.Dirac_rank();++i)
		{
		    spacetime_index_ranges.push_back(Dirac_dim<model_t::dimension>::value);
		    block_size*=(Dirac_dim<model_t::dimension>::value);
		}
		index_ranges=spacetime_index_ranges;

		/* Determining fermionic properties: */

		if(s.fermionic())
		{
		    fermion_nr=1;
		    Majorana_tag=true;
		}
		else
		{
		    fermion_nr=0;
		    Majorana_tag=false;
		}

		/* Increasing the static flavour number: */

		++flavour_nr;
	    }

	    /* Function assigning a set of wave functions of type to the
	     * particle: */

	    template<class particle_t>void set_particle_type()
	    {
		/* Initialise the wave function class: */

		particle_t::initialise();

		/* Copy all helicity states: */
		
		pos_hel_incoming_massless_states.resize(particle_t::pos_hel_incoming_massless_states.size());
		std::copy(particle_t::pos_hel_incoming_massless_states.begin(),particle_t::pos_hel_incoming_massless_states.end(),pos_hel_incoming_massless_states.begin());
		
		zero_hel_incoming_massless_state=particle_t::zero_hel_incoming_massless_state;
		
		neg_hel_incoming_massless_states.resize(particle_t::neg_hel_incoming_massless_states.size());
		std::copy(particle_t::neg_hel_incoming_massless_states.begin(),particle_t::neg_hel_incoming_massless_states.end(),neg_hel_incoming_massless_states.begin());
		
		pos_hel_outgoing_massless_states.resize(particle_t::pos_hel_outgoing_massless_states.size());
		std::copy(particle_t::pos_hel_outgoing_massless_states.begin(),particle_t::pos_hel_outgoing_massless_states.end(),pos_hel_outgoing_massless_states.begin());
		
		zero_hel_outgoing_massless_state=particle_t::zero_hel_outgoing_massless_state;
		
		neg_hel_outgoing_massless_states.resize(particle_t::neg_hel_outgoing_massless_states.size());
		std::copy(particle_t::neg_hel_outgoing_massless_states.begin(),particle_t::neg_hel_outgoing_massless_states.end(),neg_hel_outgoing_massless_states.begin());

		pos_hel_incoming_massive_states.resize(particle_t::pos_hel_incoming_massive_states.size());
		std::copy(particle_t::pos_hel_incoming_massive_states.begin(),particle_t::pos_hel_incoming_massive_states.end(),pos_hel_incoming_massive_states.begin());
		
		zero_hel_incoming_massive_state=particle_t::zero_hel_incoming_massive_state;
		
		neg_hel_incoming_massive_states.resize(particle_t::neg_hel_incoming_massive_states.size());
		std::copy(particle_t::neg_hel_incoming_massive_states.begin(),particle_t::neg_hel_incoming_massive_states.end(),neg_hel_incoming_massive_states.begin());
		
		pos_hel_outgoing_massive_states.resize(particle_t::pos_hel_outgoing_massive_states.size());
		std::copy(particle_t::pos_hel_outgoing_massive_states.begin(),particle_t::pos_hel_outgoing_massive_states.end(),pos_hel_outgoing_massive_states.begin());
		
		zero_hel_outgoing_massive_state=particle_t::zero_hel_outgoing_massive_state;
		
		neg_hel_outgoing_massive_states.resize(particle_t::neg_hel_outgoing_massive_states.size());
		std::copy(particle_t::neg_hel_outgoing_massive_states.begin(),particle_t::neg_hel_outgoing_massive_states.end(),neg_hel_outgoing_massive_states.begin());
		
		/* Copy the static propagating member function addresses: */
		
		propagator=&(particle_t::propagator_type::evaluate_range);
		refresh_prop=&(particle_t::propagator_type::refresh);
		
		/* Initialise contraction and charge conjugate contraction
		 * classes: */

		evaluate<typename particle_t::contraction_type>::initialise();
		evaluate< charge_conj_contr<typename particle_t::contraction_type,particle_t::contraction_type::fermionic> >::initialise();
		
		/* Copy the static apply function address of both contraction
		 * and charge conjugate contraction class: */
		
		contraction=&evaluate<typename particle_t::contraction_type>::apply;
		if(Majorana_tag)
		{
		    conj_contraction=&evaluate<charge_conj_contr<typename particle_t::contraction_type,particle_t::contraction_type::fermionic> >::apply;
		}
		else
		{
		    conj_contraction=contraction;
		}
	    }

	    /* Function assigning a propagator of type prop_t to the particle:
	     * */

	    template<class prop_t>bool set_propagator()
	    {
		/* Initialise the propagator class: */

		prop_t::initialise();

		std::vector<size_type>I;

		prop_t::get_index_ranges(I);

		if(I!=spacetime_index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid propagator assignment for particle "<<name<<endlog;
		    return false;
		}
		
		/* Copy the static propagating member function addresses: */
		
		propagator=&prop_t::evaluate_range;
		refresh_prop=&prop_t::refresh;
		return true;
	    }

	    template<class contr_t>bool set_contraction()
	    {
		/* Initialise the contraction class: */

		evaluate<contr_t>::initialise();

		std::vector<size_type>I;

		evaluate<contr_t>::fill_rank_vector(I);

		if(I!=index_ranges)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid wave function contraction assignment for particle "<<name<<endlog;
		    return false;
		}
		
		/* Copy the static propagating member function addresses: */
		
		contraction=&evaluate<contr_t>::apply;
		if(Majorana_tag)
		{
		    conj_contraction=&evaluate<charge_conj_contr<contr_t,contr_t::fermionic> >::apply;
		}
		else
		{
		    conj_contraction=contraction;
		}
		return true;
	    }

	private:

	    /* Trivial copy constructor: */

	    particle(particle<model_t>& other){}

	    /* Trivial assignment operator: */

	    particle<model_t>& operator = (const particle<model_t>& other){}

	    /* Propagator-denominator function for auxiliary particles: */

	    static void trivial_refresh(const momentum_type* p,const r_value_type* m,const r_value_type* W){}
	    
	    /* Propagator function for auxiliary particles: */
	    
	    static void trivial_prop(iterator first,const iterator last)
	    {
		do
		{
		    (*first)*=value_type(0,1);
		    ++first;
		}
		while(first != last);
	    }
    };

    template<class model_t>int particle<model_t>::flavour_nr=0;

    template<class model_t>std::ostream& operator << (std::ostream& os,const Camgen::particle<model_t>& phi)
    {
	os<<std::setw(3)<<phi.get_flavour();
	os<<std::setw(10)<<phi.get_pdg_id();
	os<<std::setw(12)<<phi.get_name();
	os<<std::setw(6)<<phi.get_spin();
	os<<std::setw(12)<<phi.get_mass();
	os<<std::setw(12)<<phi.get_width();
	os<<std::setw(12)<<phi.get_anti_particle()->get_name();
	os<<std::setw(5)<<phi.get_fermion_number();
	if(phi.is_fermion())
	{
	    if(phi.is_Majorana())
	    {
		os<<std::setw(5)<<"M";
	    }
	    else
	    {
		if(phi.is_Dirac())
		{
		    os<<std::setw(5)<<"D";
		}
		else
		{
		    os<<std::setw(5)<<"AD";
		}
	    }
	}
	else
	{
	    os<<std::setw(5)<<"B";
	}
	if(phi.get_colour_rank()==0)
	{
	    os<<std::setw(12)<<"[1]";
	}
	else
	{
	    std::stringstream s;
	    unsigned n=phi.get_colour_rank()-1;
	    s<<"[";
	    for(unsigned i=0;i<n;++i)
	    {
		s<<phi.get_colour_index_range(i)<<",";
	    }
	    s<<phi.get_colour_index_range(n)<<"]";
	    os<<std::setw(12)<<s.str();
	}
	os<<std::setw(5)<<phi.is_coupled();
	return os;
    }
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_PARTICLE_H_*/

