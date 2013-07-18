//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file s_gen.h
    \brief Propagator sampler base class declaration and implementation.
 */

#ifndef CAMGEN_S_GEN_H_
#define CAMGEN_S_GEN_H_

#include <limits>
#include <Camgen/utils.h>
#include <Camgen/MC_gen.h>
#include <Camgen/rn_strm.h>
#include <Camgen/histogram.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and implementation of the invariant mass generator class. This  *
 * is an abstract base class objects sampling particle propagators.            *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t,class rng_t>class s_generator;
    template<class value_t,class rng_t>class adaptive_s_generator;
    template<class value_t,class rng_t,class sub_type>class inversion_s_generator;
    template<class value_t,class rng_t>class BW_s_generator;
    template<class value_t,class rng_t>class pl_s_generator;
    template<class value_t,class rng_t>class uni_s_generator;
    template<class value_t,class rng_t>class Dd_s_generator;
    template<class value_t,class rng_t>class inv_cosh_y_generator;
    template<class value_t,class rng_t>value_t integrate(const s_generator<value_t,rng_t>*,const s_generator<value_t,rng_t>*,const value_t&);
    
    /// Momentum channel base class for phase space generators.

    template<class value_t,class rng_t>class s_generator: public MC_generator<value_t>
    {
	friend class adaptive_s_generator<value_t,rng_t>;
	friend value_t integrate<value_t,rng_t>(const s_generator<value_t,rng_t>*,const s_generator<value_t,rng_t>*,const value_t&);
	
	typedef MC_generator<value_t> base_type;

	public:

	    /* Type definitions: */
	    
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;

	    /* Public static methods: */
	    /*------------------------*/

	    /// Static serialized factory method.
	    
	    static s_generator<value_t,rng_t>* create_instance(std::istream& is,const value_t* mass,const value_t* width,const value_t* nu)
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<s_gen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read--returning NULL"<<endlog;
		    return NULL;
		}
		std::string id;
		is>>id;
		bool adaptive=false;
		if(id[id.size()-1]=='*')
		{
		    adaptive=true;
		    id.resize(id.size()-1);
		}
		s_generator<value_t,rng_t>* sgen;
		if(id=="Dd")
		{
		    sgen=new Dd_s_generator<value_t,rng_t>(mass);
		}
		else if(id=="uni")
		{
		    sgen=new uni_s_generator<value_t,rng_t>();
		}
		else if(id=="pl")
		{
		    sgen=new pl_s_generator<value_t,rng_t>(mass,nu);
		}
		else if(id=="BW")
		{
		    sgen=new BW_s_generator<value_t,rng_t>(mass,width);
		}
		else if(id=="invcosh")
		{
		    sgen=new inv_cosh_y_generator<value_t,rng_t>(mass);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-generator id "<<id<<" not recognised--returning NULL."<<endlog;
		    std::string endflag;
		    do
		    {
			std::getline(is,endflag);
		    }
		    while(endflag!="</s_gen>" and !is.eof());
		    return NULL;
		}
		s_generator<value_t,rng_t>* result;
		if(adaptive)
		{
		    result=new adaptive_s_generator<value_t,rng_t>(sgen,500);
		}
		else
		{
		    result=sgen;
		}
		result->load(is);
		std::string endflag;
		do
		{
		    std::getline(is,endflag);
		}
		while(endflag!="</s_gen>" and !is.eof());
		return result;
	    }

	    /// Static serialized factory method with external invariant mass.
	    
	    static s_generator<value_t,rng_t>* create_instance(value_type* s,std::istream& is,const value_type* mass,const value_type* width,const value_type* nu)    
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<s_gen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read--returning NULL"<<endlog;
		    return NULL;
		}
		std::string id;
		is>>id;
		bool adaptive=false;
		if(id[id.size()-1]=='*')
		{
		    adaptive=true;
		    id.resize(id.size()-1);
		}
		s_generator<value_t,rng_t>* sgen;
		if(id=="Dd")
		{
		    sgen=new Dd_s_generator<value_t,rng_t>(s,mass);
		}
		else if(id=="uni")
		{
		    sgen=new uni_s_generator<value_t,rng_t>(s);
		}
		else if(id=="pl")
		{
		    sgen=new pl_s_generator<value_t,rng_t>(s,mass,nu);
		}
		else if(id=="BW")
		{
		    sgen=new BW_s_generator<value_t,rng_t>(s,mass,width);
		}
		else if(id=="invcosh")
		{
		    sgen=new inv_cosh_y_generator<value_t,rng_t>(s,mass);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-generator id "<<id<<" not recognised--returning NULL."<<endlog;
		    std::string endflag;
		    do
		    {
			std::getline(is,endflag);
		    }
		    while(endflag!="</s_gen>" and !is.eof());
		    return NULL;
		}
		s_generator<value_t,rng_t>* result;
		if(adaptive)
		{
		    result=new adaptive_s_generator<value_t,rng_t>(s,sgen);
		}
		else
		{
		    result=sgen;
		}
		result->load(is);
		std::string endflag;
		do
		{
		    std::getline(is,endflag);
		}
		while(endflag!="</s_gen>" and !is.eof());
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /// Default constructor.
	    
	    s_generator():smin(-std::numeric_limits<value_t>::infinity()),smax(std::numeric_limits<value_t>::infinity()),mmin(smin),mmax(smax),norm(0),invmass(new value_type),alloc_s(true){}
	    
	    /// Constructor with invariant mass pointer argument.
	    
	    s_generator(value_type* invmass_):smin(-std::numeric_limits<value_t>::infinity()),smax(std::numeric_limits<value_t>::infinity()),mmin(smin),mmax(smax),norm(0),invmass(invmass_),alloc_s(false){}

	    /// Copy constructor.

	    s_generator(const s_generator<value_t,rng_t>& other):MC_generator<value_t>(other),smin(other.smin),smax(other.smax),mmin(other.mmin),mmax(other.mmax),norm(0),alloc_s(other.alloc_s)
	    {
		invmass=alloc_s?(new value_type(*(other.invmass))):(other.invmass);
	    }

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    virtual ~s_generator()
	    {
		if(alloc_s)
		{
		    delete invmass;
		}
	    }

	    /* Public modifiers: */
	    /*--------------------*/

	    /// Sets the invariant mass address to the argument pointer.
	    
	    void set_s(value_type* s_)
	    {
		if(alloc_s)
		{
		    delete invmass;
		    alloc_s=false;
		}
		invmass=s_;
	    }

	    /// Sets the minimal generated invariant mass-squared.

	    bool set_s_min(const value_type& smin_)
	    {
		smin=smin_;
		mmin=sgn_sqrt(smin_);
		return refresh_s_min();
	    }

	    /// Sets the minimal minimal invariant mass-squared.

	    virtual bool set_s_min_min(const value_type& sminmin)
	    {
		return set_s_min(sminmin);
	    }

	    /// Sets the maximal generated invariant mass-squared.

	    bool set_s_max(const value_type& smax_)
	    {
		smax=smax_;
		mmax=sgn_sqrt(smax_);
		return refresh_s_max();
	    }

	    /// Sets the maximal maximal invariant mass-squared.

	    virtual bool set_s_max_max(const value_type& smaxmax)
	    {
		return set_s_max(smaxmax);
	    }

	    /// Sets the invariant mass-squared range.

	    bool set_s_range(const value_type& smin_,const value_type& smax_)
	    {
		smin=smin_;
		mmin=sgn_sqrt(smin_);
		smax=smax_;
		mmax=sgn_sqrt(smax_);
		return refresh_s_bounds();
	    }

	    /// Sets the minimal generated invariant mass.

	    bool set_m_min(const value_type& mmin_)
	    {
		mmin=mmin_;
		smin=sgn_sq(mmin_);
		return refresh_s_min();
	    }

	    /// Sets the maximal generated invariant mass.

	    bool set_m_max(const value_type& mmax_)
	    {
		mmax=mmax_;
		smax=sgn_sq(mmax_);
		return refresh_s_max();
	    }

	    /// Sets the invariant mass range.

	    bool set_m_range(const value_type& mmin_,const value_type& mmax_)
	    {
		mmin=mmin_;
		smin=sgn_sq(mmin_);
		mmax=mmax_;
		smax=sgn_sq(mmax_);
		return refresh_s_bounds();
	    }

	    /// Generation method.

	    virtual bool generate()
	    {
		if(!normalisable())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		*invmass=map(rn_stream::throw_number());
		this->weight()=get_fast_weight();
		return true;
	    }

	    /// Weight evaluation method.

	    virtual bool evaluate_weight()
	    {
		if(!normalisable() or !in_range())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->weight()=get_fast_weight();
		return true;
	    }

	    /// Mapping weight factor.

	    virtual value_type mapping_weight() const
	    {
		return this->weight();
	    }

	    /// Adaptive grid weight factor.

	    virtual value_type grid_weight() const
	    {
		return 1;
	    }

	    /// Re-evaluates internal parameters from external pointer values. Returns
	    /// false if the input parameters are invalid.

	    virtual bool refresh_params()
	    {
		refresh_norm();
		return normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on smin. Returns false if
	    /// the input parameters are invalid.

	    virtual bool refresh_s_min()
	    {
		refresh_norm();
		return normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on smax. Returns false if
	    /// the input parameters are invalid.

	    virtual bool refresh_s_max()
	    {
		refresh_norm();
		return normalisable();
	    }

	    /// Re-evaluates all internal parameters depending on smin and smax.
	    /// Returns false if the input parameters are invalid.

	    virtual bool refresh_s_bounds()
	    {
		refresh_s_min();
		refresh_s_max();
		return normalisable();
	    }

	    /// Updates adaptive grids given the total event weight.

	    virtual void update_weight(){}

	    /* Public constant methods: */
	    /*--------------------------*/

	    /// Instance cloning method. Returns the NULL pointer if not
	    /// overridden.

	    virtual s_generator<value_t,rng_t>* clone() const
	    {
		return NULL;
	    }

	    /// Maps the first argument, a random number in [0,1] to an
	    /// invariant mass.

	    virtual value_type map(const value_type&) const=0;

	    /// Maps the argument, an invariant mass-squared between [s-,s+], to
	    /// a random number in [0,1].

	    virtual value_type inverse_map(const value_type&) const=0;

	    /// Returns whether the requested density is nomalisable within
	    /// current interval.

	    bool normalisable() const
	    {
		return (norm>0 and norm!=std::numeric_limits<value_t>::infinity());
	    }

	    /// Returns whether the s-range is bounded.

	    bool bounded_s_range() const
	    {
		return !((smin==-std::numeric_limits<value_t>::infinity()) or (smax==std::numeric_limits<value_t>::infinity()));
	    }

	    /// Returns the probability density at the given invariant mass.

	    value_type density() const
	    {
		if(this->in_range())
		{
		    return this->norm/this->weight();
		}
		return 0;
	    }

	    /// Returns a reference to the invariant mass.
	    
	    value_type& s()
	    {
		return *invmass;
	    }

	    /// Returns a constant reference to the invariant mass.
	    
	    const value_type& s() const
	    {
		return *invmass;
	    }

	    /// Returns the invariant mass address.

	    value_type* get_s()
	    {
		return invmass;
	    }

	    /// Returns the constant invariant mass address.

	    const value_type* get_s() const
	    {
		return invmass;
	    }

	    /// Returns the minimal invariant mass-squared value.

	    const value_type& s_min() const
	    {
		return smin;
	    }

	    /// Returns the maximal invariant mass-squared value.

	    const value_type& s_max() const
	    {
		return smax;
	    }

	    /// Returns the minimal invariant mass value.

	    const value_type& m_min() const
	    {
		return mmin;
	    }

	    /// Returns the maximal invariant mass value.

	    const value_type& m_max() const
	    {
		return mmax;
	    }

	    /// Returns whether the invariant mass is within limits.

	    bool in_range() const
	    {
		return (*invmass>=smin and *invmass<=smax);
	    }

	    /// Double-dispatch integrator functions.
	    
	    virtual value_type integrate_with(const s_generator<value_t,rng_t>*,const value_type&) const=0;
	    virtual value_type integrate_with(const BW_s_generator<value_t,rng_t>*,const value_type&) const=0;
	    virtual value_type integrate_with(const pl_s_generator<value_t,rng_t>*,const value_type&) const=0;
	    virtual value_type integrate_with(const uni_s_generator<value_t,rng_t>*,const value_type&) const=0;
	    virtual value_type integrate_with(const Dd_s_generator<value_t,rng_t>*,const value_type&) const=0;

	    /* Testing methods: */
	    /*------------------*/

	    /// Checking function.

	    bool check()
	    {
		if(!in_range())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed check: invariant mass "<<s()<<" not within range ["<<smin<<','<<smax<<"]"<<endlog;
		    return false;
		}
		value_type w=this->weight();
		this->evaluate_weight();
		if(w!=this->weight())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"failed check: recomputed weight "<<this->weight()<<" does not equal generation weight "<<w<<endlog;
		    return false;
		}
		return true;
	    }

	    /// Testrun of N_events generations returning a histogram with
	    /// N_bins bins.

	    histogram<value_type> testrun(size_type N_events,size_type N_bins)
	    {
		histogram<value_type>hist(invmass,this->get_weight(),N_events);
		for(size_type n=0;n<N_events;++n)
		{
		    if(this->generate())
		    {
			check();
		    }
		    else
		    {
			this->weight()=(value_type)0;
		    }
		    hist.store();
		}
		hist.make(N_bins);
		return hist;
	    }

	    /* Serialization: */
	    /*----------------*/
	    
	    /// Returns a string containing the type.

	    virtual std::string type() const=0;

	    /// Returns the detailed type of the generator.

	    virtual std::string detailed_type() const
	    {
		return type();
	    }

	    /// Printing method:

	    virtual std::ostream& print(std::ostream& os) const
	    {
		os<<type();
		return os;
	    }

	    /// Loads derived class data (empty by default).

	    virtual std::istream& load_data(std::istream& is)
	    {
		return is;
	    }

	    /// Overridden loading method.

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		safe_read(is,smin);
		safe_read(is,smax);
		load_data(is);
		refresh_s_bounds();
		return is;
	    }

	    /// Saves derived class data (empty by default).

	    virtual std::ostream& save_data(std::ostream& os) const
	    {
		return os;
	    }

	    /// Overridden saving method.

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<s_gen>"<<std::endl;
		os<<type()<<std::endl;
		this->base_type::save(os);
		safe_write(os,smin);
		os<<"\t";
		safe_write(os,smax);
		os<<std::endl;
		save_data(os);
		os<<"</s_gen>"<<std::endl;
		return os;
	    }

	protected:
	    
	    /// Minimal and maximal invariant mass-squared.

	    value_type smin,smax;
	    
	    /// Minimal and maximal invariant mass.

	    value_type mmin,mmax;

	    /// Normalisation constant.

	    value_type norm;

	    /// Evaluates the (private) normalisation factor.

	    virtual void refresh_norm()=0;

	    /// Returns the weight of the current configuration without bound-
	    /// or norm-checking.

	    virtual value_type get_fast_weight() const=0;

	private:

	    /* Pointer to the generated invariant mass: */

	    value_type* invmass;

	    /* Boolean denoting whether the invariant mass address was allocated
	     * dynamically: */

	    bool alloc_s;
    };

    template<class value_t,class rng_t>value_t integrate(const s_generator<value_t,rng_t>* first,const s_generator<value_t,rng_t>* second,const value_t& sqrts)
    {
	return first->integrate_with(second,sqrts)/((first->norm)*(second->norm));
    }

    template<class value_t,class rng_t>value_t integrate_unnormalised(const s_generator<value_t,rng_t>* first,const s_generator<value_t,rng_t>* second,const value_t& sqrts)
    {
	return first->integrate_with(second,sqrts); 
    }
}

#include <Camgen/s_int.h>
#include <Camgen/Breit_Wigner.h>
#include <Camgen/power_law.h>
#include <Camgen/Dirac_delta.h>
#include <Camgen/uni_s.h>
#include <Camgen/inv_cosh.h>
#include <Camgen/sgen_grid.h>

#endif /*CAMGEN_S_GEN_H_*/

