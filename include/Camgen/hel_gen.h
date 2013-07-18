//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file hel_gen.h
    \brief helicity generator abstract base class template interface and implementation.
 */

#ifndef CAMGEN_HEL_GEN_H_
#define CAMGEN_HEL_GEN_H_

#include <bitset>
#include <vector>
#include <complex>
#include <Camgen/unused.h>
#include <Camgen/phase_space.h>
#include <Camgen/MC_obj_gen.h>
#include <Camgen/summation.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstact base class templates for helicity generators. To construct a helicity *
 * generator for Camgen, derive from this class and use the protected functions *
 * to access particle helicity variables. The class can be initialized in        *
 * 'standalone mode', where helicity variables are allocated, or in association  *
 * with a process tree in Camgen, where the helicity variables point to the     *
 * external particles' degrees of freedom.                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t,std::size_t N_in,std::size_t N_out,bool cont_hels>class helicity_generator;

    /// Helicity generator base class specialization for discrete helicities. A
    /// derived class should provide an implementation of the generate() and
    /// (optionally) evaluate_weight() member functions and have constructor
    /// which takes two arrays of numbers or a process_tree list iterator as
    /// arguments.

    template<class value_t,std::size_t N_in,std::size_t N_out>class helicity_generator<value_t,N_in,N_out,false>: public MC_object_generator<value_t,int,N_in+N_out,true>
    {
	typedef MC_object_generator<value_t,int,N_in+N_out,true> base_type;

	public:

	    /* Useful type definitions: */

	    typedef value_t value_type;
	    typedef std::size_t size_type;
	    typedef int object_type;
	    typedef typename base_type::integral_type integral_type;

	    /* Public static data: */
	    /*---------------------*/

	    /* Static integer constants: */

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static member functions: */
	    /*---------------------------------*/

	    /* Utility functions for named constructors */
	    
	    template<class model_t>static vector<int*,N_tot> make_helicity_vector(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		vector<int*,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i]=&(it->get_phase_space(i)->helicity());
		}
		return result;
	    }
	    template<class model_t>static vector<int,N_tot> make_max_helicity_vector(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		vector<int,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i]=it->get_phase_space(i)->maximal_helicity();
		}
		return result;
	    }
	    template<class model_t>static std::bitset<N_tot> make_zero_helicity_bitset(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		std::bitset<N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i]=it->get_phase_space(i)->has_zero_helicity();
		}
		return result;
	    }

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    virtual ~helicity_generator(){}

	    /* Public readout functions: */
	    /*---------------------------*/

	    /// Returns the reference to the helicity of the i-th particle.

	    int& helicity(size_type i)
	    {
		return this->object(i);
	    }

	    /// Returns the const reference to the helicity of the i-th particle.

	    const int& helicity(size_type i) const
	    {
		return this->object(i);
	    }

	    /// Returns the maximal helicity of particle i.

	    int maximal_helicity(size_type i) const
	    {
		return max_hels[i];
	    }

	    /// Returns whether particle i has a zero helicity state.

	    bool has_zero_helicity(size_type i) const
	    {
		return zero_hels[i];
	    }

	    /// Returns the initial-state helicity averaging factor for cross sections.

	    value_type averaging_factor() const
	    {
		return prefactor;
	    }

	    /// Returns a spin summation flag for the i-th particle. If true,
	    /// the process/event generator classes will sum over the
	    /// corresponding helicity when computing the matrix element.
	    /// Returns false by default.

	    virtual bool sum_helicity(int i) const
	    {
		return false;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Polymorphic type identifier.

	    virtual std::string type() const=0;

	    /// Loads derived class data from input stream (empty by default).

	    virtual std::istream& load_data(std::istream& is)
	    {
		return is;
	    }

	    /// Loads the helicity generator from input stream.

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		return load_data(is);
	    }

	    /// Writes derived class data to output stream (empty by default).

	    std::ostream& save_data(std::ostream& os) const
	    {
		return os;
	    }

	    /// Writes the helicity generator to the output stream.

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<helgen>"<<std::endl;
		os<<type()<<std::endl;
		this->base_type::save(os);
		save_data(os);
		os<<"</helgen>"<<std::endl;
		return os;
	    }

	    /// Prints the generated spin configuration to the argument
	    /// streaming object.

	    std::ostream& print(std::ostream& os=std::cout) const
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    if(this->object(i)==1)
		    {
			os<<"+  ";
		    }
		    else if(this->object(i)==-1)
		    {
			os<<"-  ";
		    }
		    else
		    {
			os<<this->object(i)<<"  ";
		    }
		}
		os<<">";
		for(size_type i=N_in;i<N_tot;++i)
		{
		    if(this->object(i)==1)
		    {
			os<<"  +";
		    }
		    else if(this->object(i)==-1)
		    {
			os<<"  -";
		    }
		    else
		    {
			os<<"  "<<this->object(i);
		    }
		}
		os<<"\t\t\t"<<this->weight()<<std::endl;
		return os;
	    }

	protected:

	    /* Protected members: */
	    /*--------------------*/

	    /// Spin-averaging prefactor.

	    value_type prefactor;

	    /* Protected constructors: */
	    /*-------------------------*/
	    
	    /// Non-allocating constructor. The first arguments is a vector of
	    /// integer addresses to store the helicity configuration in. The
	    /// second argument is a vector of maximal helicity values of all
	    /// external particles. The third argument is a bitset determining
	    /// whether the particles contain a zero-helicity state.

	    helicity_generator(const vector<int*,N_tot>& hels_,const vector<int,N_tot>& max_hels_,std::bitset<N_tot> zero_hels_):base_type(hels_),prefactor(1),zero_hels(zero_hels_)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    max_hels[i]=std::abs(max_hels_[i]);
		    if(max_hels[i]==0)
		    {
			zero_hels[i]=true;
		    }
		    if(i<N_in)
		    {
			prefactor/=(zero_hels[i]?(2*max_hels[i]+1):(2*max_hels[i]));
		    }
		}
	    }

	    /// Allocating constructor. The first argument is a vector of
	    /// maximal helicity values of all external particles. The second
	    /// argument is a bitset determining whether the particles contain
	    /// a zero-helicity state.

	    helicity_generator(const vector<int,N_tot>& max_hels_,std::bitset<N_tot> zero_hels_):prefactor(1),zero_hels(zero_hels_)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    max_hels[i]=std::abs(max_hels_[i]);
		    if(max_hels[i]==0)
		    {
			zero_hels[i]=true;
		    }
		    if(i<N_in)
		    {
			prefactor/=(zero_hels[i]?(2*max_hels[i]+1):(2*max_hels[i]));
		    }
		}
	    }

	private:

	    /* Private data: */
	    /*---------------*/

	    /* Maximal helicities. */

	    vector<int,N_tot>max_hels;

	    /* Booleans denoting external particles having a helicity 0 state.
	     * */

	    std::bitset<N_tot>zero_hels;
    };
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t helicity_generator<value_t,N_in,N_out,false>::N_tot;

    /// Helicity generator base class specialization for continuous helicities. A
    /// derived class should provide an implementation of the generate() and
    /// (optionally) evaluate_weight() member functions and have constructor
    /// which takes two arrays of numbers or a process_tree list iterator as
    /// arguments.

    template<class value_t,std::size_t N_in,std::size_t N_out>class helicity_generator<value_t,N_in,N_out,true>: public MC_object_generator<value_t,helicity_phases<value_t>,N_in+N_out,true>
    {
	typedef MC_object_generator<value_t,helicity_phases<value_t>,N_in+N_out,true> base_type;

	public:

	    /* Useful type definitions: */

	    typedef value_t value_type;
	    typedef std::complex<value_t> c_value_type;
	    typedef std::size_t size_type;
	    typedef typename base_type::object_type object_type;
	    typedef typename base_type::integral_type integral_type;

	    /* Public static data: */
	    /*---------------------*/
	    
	    /* Static integer constants: */

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static member functions: */
	    /*---------------------------------*/

	    /* Utility functions for named constructors */

	    template<class model_t>static vector<object_type*,N_tot> make_helicity_vector(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		vector<object_type*,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i]=&(it->get_phase_space(i)->get_helicity_phases());
		}
		return result;
	    }

	    /* Destructors: */
	    /*--------------*/

	    /// Destructor.
	    
	    virtual ~helicity_generator(){}

	    /* Public readout functions: */
	    /*---------------------------*/

	    /// Returns a reference to helicity phase l of particle i.
	    
	    c_value_type& helicity_phase(size_type i,int l)
	    {
		return this->get_object(i)->helicity_phase(l);
	    }

	    /// Returns a const reference to helicity phase l of particle i.
	    
	    const c_value_type& helicity_phase(size_type i,int l) const
	    {
		return this->get_object(i)->helicity_phase(l);
	    }

	    /// Returns maximal helicity of particle i.

	    int maximal_helicity(size_type i) const
	    {
		return this->get_object(i)->max_helicity();
	    }

	    /// Returns a boolean denoting whether particle i contains a zero
	    /// helicity state.

	    bool has_zero_helicity(size_type i) const
	    {
		return this->get_object(i)->zero_helicity_state;
	    }

	    /// Returns the initial-state helicity averaging factor for cross sections.

	    value_type averaging_factor() const
	    {
		return prefactor;
	    }

	    /// Returns a spin summation flag for the i-th particle. If true,
	    /// the process/event generator classes will sum over the
	    /// corresponding helicity when computing the matrix element.
	    /// Returns false by default.

	    virtual bool sum_helicity(int i) const
	    {
		return false;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Polymorphic type identifier.

	    virtual std::string type() const=0;

	    /// Loads derived class data from input stream (empty by default).

	    virtual std::istream& load_data(std::istream& is)
	    {
		return is;
	    }

	    /// Loads the helicity generator from input stream.

	    std::istream& load(std::istream& is)
	    {
		this->base_type::load(is);
		return load_data(is);
	    }

	    /// Writes derived class data to output stream (empty by default).

	    std::ostream& save_data(std::ostream& os) const
	    {
		return os;
	    }

	    /// Writes the helicity generator to the output stream.

	    std::ostream& save(std::ostream& os) const
	    {
		os<<"<helgen>"<<std::endl;
		os<<type()<<std::endl;
		this->base_type::save(os);
		save_data(os);
		os<<"</helgen>"<<std::endl;
		return os;
	    }

	    /// Prints the generated spin configuration to the argument
	    /// streaming object.

	    std::ostream& print(std::ostream& os=std::cout) const
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    os<<"[";
		    for(size_type l=0;l<this->get_object(i)->max_helicity();++l)
		    {
			os<<this->get_object(i)->negative_helicity_phase(l)<<",";
		    }
		    if(this->get_object(i)->zero_helicity_state)
		    {
			os<<this->get_object(i)->zero_helicity_phase();
		    }
		    else
		    {
			os<<"0";
		    }
		    for(size_type l=0;l<this->get_object(i)->max_helicity();++l)
		    {
			os<<","<<this->get_object(i)->positive_helicity_phase(l);
		    }
		    os<<"] ";
		}
		os<<">";
		for(size_type i=N_in;i<N_tot;++i)
		{
		    os<<"[";
		    for(size_type l=0;l<this->get_object(i)->max_helicity();++l)
		    {
			os<<this->get_object(i)->negative_helicity_phase(l)<<",";
		    }
		    if(this->get_object(i)->zero_helicity_state)
		    {
			os<<this->get_object(i)->zero_helicity_phase();
		    }
		    else
		    {
			os<<"0";
		    }
		    for(size_type l=0;l<this->get_object(i)->max_helicity();++l)
		    {
			os<<","<<this->get_object(i)->positive_helicity_phase(l);
		    }
		    os<<"] ";
		}
		os<<"\t\t\t"<<this->weight()<<std::endl;
		return os;
	    }

	protected:

	    /* Protected data: */
	    /*-----------------*/

	    value_type prefactor;

	    /* Protected constructors: */
	    /*-------------------------*/

	    /// Non-allocating constructor. The first argument is the vector of
	    /// positive helicity phases, the second argument the zero-helicity
	    /// phases and the third argument the vector of negative-helicity
	    /// phases.

	    helicity_generator(const vector<helicity_phases<value_type>*,N_tot>& hels_):base_type(hels_),prefactor(1)
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    prefactor/=(this->get_object(i)->n_helicities());
		}
	    }
    };
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t helicity_generator<value_t,N_in,N_out,true>::N_tot;

    /// Helicity-summing dummy class. Signals a full helicity sum of the
    /// amplitude in process/event generator instances.

    template<class value_t,std::size_t N_in,std::size_t N_out,bool cont_hels>class helicity_summer;
    
    template<class value_t,std::size_t N_in,std::size_t N_out>class helicity_summer<value_t,N_in,N_out,false>: public helicity_generator<value_t,N_in,N_out,false>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef helicity_generator<value_t,N_in,N_out,false> base_type;

	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/
	    
	    /* Factory method: */

	    template<class model_t>static helicity_summer<value_t,N_in,N_out,false>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		return new helicity_summer<value_t,N_in,N_out,false>(base_type::template make_helicity_vector<model_t>(it),base_type::template make_max_helicity_vector<model_t>(it),base_type::template make_zero_helicity_bitset<model_t>(it));
	    }

	    /* Public constructors: */
	    /*----------------------*/
	    
	    /* Constructor: */

	    helicity_summer(const vector<int*,N_tot>& hels_,const vector<int,N_tot>& max_hels_,std::bitset<N_tot> zero_hels_):base_type(hels_,max_hels_,zero_hels_){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Dummy generation method: */

	    bool generate()
	    {
		this->weight()=(value_t)1;
		return true;
	    }

	    /* Weight evaluation method: */
	    
	    bool evaluate_weight()
	    {
		this->weight()=(value_t)1;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Cloning method: */

	    helicity_summer<value_t,N_in,N_out,false>* clone() const
	    {
		return new helicity_summer<value_t,N_in,N_out,false>(*this);
	    }

	    /* Helicity summation flags set: */

	    bool sum_helicity(int i) const
	    {
		return true;
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "sum";
	    }
    };
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t helicity_summer<value_t,N_in,N_out,false>::N_tot;

    template<class value_t,std::size_t N_in,std::size_t N_out>class helicity_summer<value_t,N_in,N_out,true>: public helicity_generator<value_t,N_in,N_out,true>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/
	    
	    typedef helicity_generator<value_t,N_in,N_out,true> base_type;
	    typedef value_t value_type;
	    
	    /* Public static data: */
	    /*---------------------*/
	    
	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/
	    
	    /* Factory method: */

	    template<class model_t>static helicity_summer<value_t,N_in,N_out,true>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		return new helicity_summer<value_t,N_in,N_out,true>(base_type::template make_helicity_vector<model_t>(it));
	    }

	    /* Public constructors: */
	    /*----------------------*/
	    
	    /* Constructor: */
	    
	    helicity_summer(const vector<helicity_phases<value_type>*,N_tot>& hels_):base_type(hels_){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Dummy generation method: */

	    bool generate()
	    {
		this->weight()=(value_t)1;
		return true;
	    }

	    /* Weight evaluation method: */
	    
	    bool evaluate_weight()
	    {
		this->weight()=(value_t)1;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Cloning method: */

	    helicity_summer<value_t,N_in,N_out,true>* clone() const
	    {
		return new helicity_summer<value_t,N_in,N_out,true>(*this);
	    }

	    /* Helicity summation flags set: */

	    bool sum_helicity(int i) const
	    {
		return true;
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "sum";
	    }
    };
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t helicity_summer<value_t,N_in,N_out,true>::N_tot;
}

#endif /*CAMGEN_HEL_GEN_H_*/

