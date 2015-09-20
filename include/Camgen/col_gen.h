//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file col_gen.h
    \brief colour generator abstract base class template interface and implementation.
 */

#ifndef CAMGEN_COL_GEN_H_
#define CAMGEN_COL_GEN_H_

#include <vector>
#include <complex>
#include <Camgen/vector.h>
#include <Camgen/tensor.h>
#include <Camgen/forward_decs.h>
#include <Camgen/obj_alloc.h>
#include <Camgen/MC_gen.h>
#include <Camgen/summation.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class template definition and declaration for colour generators *
 * in Camgen.                                                                   *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t,std::size_t N_in,std::size_t N_out,bool continuous_colours>class colour_generator;
    
    /// Colour generator base class specialization for discrete colours. A
    /// derived class should provide an implementation of the generate() and
    /// (optionally) evaluate_weight() member functions and have constructor
    /// which takes a number of colours and an array of colour ranks, or a tree
    /// iterator as arguments.

    template<class value_t,std::size_t N_in,std::size_t N_out>class colour_generator<value_t,N_in,N_out,false>: public MC_generator<value_t>, public object_allocator<std::vector<std::size_t> >
    {
	typedef object_allocator< std::vector<std::size_t> > base_type;
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef value_t value_type;
	    typedef typename base_type::object_type object_type;
	    typedef typename base_type::size_type size_type;

	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Utility functions for named constructors: */

	    template<class model_t>static vector<std::vector<size_type>*,N_tot> make_colour_vector(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		vector<std::vector<size_type>*,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i]=it->get_phase_space(i)->get_colours();
		}
		return result;
	    }
	    template<class model_t>static vector<std::vector<size_type>,N_tot> make_colour_range_vector(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		vector<std::vector<size_type>,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i]=it->get_phase_space(i)->particle_type->get_colour_index_ranges();
		}
		return result;
	    }
	    template<class model_t>static value_type make_prefactor(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		value_type result(1);
		for(size_type i=0;i<N_in;++i)
		{
		    result/=it->get_phase_space(i)->particle_type->colour_dof();
		}
		return result;
	    }

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    virtual ~colour_generator(){}

	    /* Public const methods: */
	    /*-----------------------*/

	    /// Virtual cloning method.

	    virtual colour_generator<value_t,N_in,N_out,false>* clone() const=0;

	    /// Returns the colour rank of particle i.

	    size_type colour_rank(size_type i) const
	    {
		return col_ranges[i].size();
	    }

	    /// Returns the range of the j-th colour degree of freedom of
	    /// particle i

	    size_type colour_range(size_type i,size_type j) const
	    {
		return col_ranges[i][j];
	    }

	    /// Returns the reference to the j-th colour of particle i.

	    size_type& colour(size_type i,size_type j)
	    {
		return (this->object(i))[j];
	    }

	    /// Returns the const reference to the j-th colour of particle i.

	    const size_type& colour(size_type i,size_type j) const
	    {
		return (this->object(i))[j];
	    }

	    /// Returns a colour summation flag for the i-th particle. If true,
	    /// the process/event generator classes will sum over the
	    /// corresponding colour when computing the matrix element.
	    /// Returns false by default.

	    virtual bool sum_colour(int i) const
	    {
		return false;
	    }

	    /// Returns the initial state colours averaging factor.

	    value_type averaging_factor() const
	    {
		return prefactor;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Polymorphic type identifier.

	    virtual std::string type() const=0;

	    /// Virtual Les Houches event file output method.

	    virtual void LH_output(vector<int,N_tot>& c,vector<int,N_tot>& cbar) const
	    {
		c.assign(0);
		cbar.assign(0);
	    }

	    /// Prints the generated colour configuration to the argument
	    /// streaming object.

	    std::ostream& print(std::ostream& os=std::cout) const
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    os<<"(";
		    if(colour_rank(i)>0)
		    {
			for(size_type j=0;j<colour_rank(i)-1;++j)
			{
			    os<<colour(i,j)<<",";
			}
			os<<colour(i,colour_rank(i)-1);
		    }
		    os<<")";
		}
		os<<" > ";
		for(size_type i=N_in;i<N_tot;++i)
		{
		    os<<"(";
		    if(colour_rank(i)>0)
		    {
			for(size_type j=0;j<colour_rank(i)-1;++j)
			{
			    os<<colour(i,j)<<",";
			}
			os<<colour(i,colour_rank(i)-1);
		    }
		    os<<")";
		}
		os<<"\t\t\t"<<this->weight()<<std::endl;
		return os;
	    }

	protected:

	    /// Initial state colours averaging factor.
	    
	    value_type prefactor;

	    /* Protected constructors: */
	    /*-------------------------*/

	    /// Non-allocating constructor. The first argument denotes the
	    /// vector of color vector pointers to be filled by the generator. The
	    /// second argument is in array of colour range vectors.

	    colour_generator(const vector<object_type*,N_tot>& cols_,const vector<std::vector<size_type>,N_tot>& col_ranges_):base_type(cols_.begin(),cols_.end()),prefactor(1),col_ranges(col_ranges_){}

	    /// Non-allocating constructor with a global range argument.

	    colour_generator(const vector<object_type*,N_tot>& cols_,size_type Nc):base_type(cols_.begin(),cols_.end()),prefactor(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    col_ranges[i]=std::vector<size_type>(cols_[i]->size(),Nc);
		}
	    }

	    /// Standalone-mode constructor. The argument is an array of colour
	    /// range vectors.

	    colour_generator(const vector<std::vector<size_type>,N_tot>& col_ranges_):base_type(N_tot),col_ranges(col_ranges_),prefactor(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    this->object(i).resize(col_ranges[i].size(),0);
		}
	    }

	    /// Standalone-mode constructor. The first argument denotes the
	    /// number of colours, the second vector the colour ranks of the
	    /// external particles.

	    colour_generator(const vector<size_type,N_tot>& col_ranks_,size_type Nc):base_type(N_tot),prefactor(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    col_ranges[i]=std::vector<size_type>(col_ranks_[i],Nc);
		    this->object(i).resize(col_ranks_[i],0);
		}
	    }
	
	private:

	    /* Private data: */
	    /*---------------*/

	    /* External state colour range vectors. */

	    vector<std::vector<size_type>,N_tot>col_ranges;
    };
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t colour_generator<value_t,N_in,N_out,false>::N_tot;
    
    /// Colour generator base class specialization for continuous colours. A
    /// derived class should provide an implementation of the generate() and
    /// (optionally) evaluate_weight() member functions and have constructor
    /// which takes a number of colours and an array of colour ranks, or a tree
    /// iterator as arguments.

    template<class value_t,std::size_t N_in,std::size_t N_out>class colour_generator<value_t,N_in,N_out,true>: public object_allocator< tensor< std::complex<value_t> > >, public MC_generator<value_t>
    {
	typedef object_allocator<tensor< std::complex<value_t> > > base_type;

	public:
	    
	    /* Type definitions: */
	    /*-------------------*/

	    typedef value_t value_type;
	    typedef typename base_type::object_type object_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename object_type::value_type c_value_type;

	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Utility functions for named constructors: */

	    template<class model_t>static vector<object_type*,N_tot> make_colour_tensor_vector(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		vector<object_type*,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i]=it->get_phase_space(i)->get_colours();
		}
		return result;
	    }
	    template<class model_t>static value_type make_prefactor(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		value_type result(1);
		for(size_type i=0;i<N_in;++i)
		{
		    result/=it->get_phase_space(i)->particle_type->colour_dof();
		}
		return result;
	    }

	    /* Destructors: */
	    /*--------------*/

	    /// Destructor.

	    virtual ~colour_generator(){}

	    /* Public const methods: */
	    /*-----------------------*/

	    /// Virtual cloning method.

	    virtual colour_generator<value_t,N_in,N_out,true>* clone() const=0;

	    /// Returns the colour rank of particle i.

	    size_type colour_rank(size_type i) const
	    {
		return this->object(i).rank();
	    }

	    /// Returns the range of the j-th colour of particle i.

	    size_type colour_range(size_type i,size_type j) const
	    {
		return this->object(i).index_range(j);
	    }

	    /// Returns the reference to the j-th colour entry of particle i.

	    c_value_type& colour_entry(size_type i,size_type j)
	    {
		return (this->object(i))[j];
	    }

	    /// Returns the const reference to the j-th colour entry of particle i.

	    const c_value_type& colour_entry(size_type i,size_type j) const
	    {
		return (this->object(i))[j];
	    }

	    /// Returns the size of the colour tensor of particle i.

	    size_type colour_size(size_type i) const
	    {
		return this->object(i).size();
	    }
	    
	    /// Returns a reference to a colour component of particle k.

	    c_value_type& colour_coeff(size_type k,size_type i0,...)
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<colour_rank(k);++j)
		{
		    n+=(i*this->object(k).block_size(j));
		    i=va_arg(I,size_type);
		}
		va_end(I);
		return colour_entry(k,n);
	    }
	    
	    /// Returns a const reference to a colour component of particle k.

	    const c_value_type& colour_coeff(size_type k,size_type i0,...) const
	    {
		size_type n=0;
		va_list I;
		size_type i=i0;
		va_start(I,i0);
		for(size_type j=0;j<colour_rank(k);++j)
		{
		    n+=(i*this->object(k).block_size(j));
		    i=va_arg(I,size_type);
		}
		va_end(I);
		return colour_entry(k,n);
	    }

	    /// Returns a colour summation flag for the i-th particle. If true,
	    /// the process/event generator classes will sum over the
	    /// corresponding colour when computing the matrix element.
	    /// Returns false by default.

	    virtual bool sum_colour(int i) const
	    {
		return false;
	    }

	    /// Returns the initial state colours averaging factor.

	    value_type averaging_factor() const
	    {
		return prefactor;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /// Polymorphic type identifier.

	    virtual std::string type() const=0;

	    /// Virtual Les Houches event file output.

	    virtual void LH_output(vector<int,N_tot>& c,vector<int,N_tot>& cbar) const
	    {
		c.assign(0);
		cbar.assign(0);
	    }

	    /// Prints the generated colour configuration to the argument
	    /// streaming object.

	    std::ostream& print(std::ostream& os=std::cout) const
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    os<<this->object(i)<<std::endl;
		}
		os<<" > "<<std::endl;
		for(size_type i=N_in;i<N_tot;++i)
		{
		    os<<this->object(i)<<std::endl;
		}
		os<<"\t\t\t"<<this->weight()<<std::endl;
		return os;
	    }

	protected:

	    /* Protected data members: */
	    /*-------------------------*/

	    value_type prefactor;

	    /* Protected constructors: */
	    /*-------------------------*/

	    /// Non-allocating constructor, taking a vector of pointers to
	    /// colour tensors as argument.

	    colour_generator(const vector<object_type*,N_tot>& tensors):base_type(tensors.begin(),tensors.end()),prefactor(1){}

	    /// Allocating constructor. The argument is an array of colour
	    /// range vectors.

	    colour_generator(const vector<std::vector<size_type>,N_tot>& ranges):base_type(N_tot),prefactor(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    this->object(i).resize(ranges[i]);
		}
	    }
	    
	    /// Allocating constructor. The first argument denotes the
	    /// number of colours, the second denotes the colour ranks of
	    /// the external particles.

	    colour_generator(vector<size_type,N_tot>ranks,size_type Nc):base_type(N_tot),prefactor(1)
	    {
		std::vector<size_type> r;
		for(size_type i=0;i<N_tot;++i)
		{
		    r.resize(ranks[i],Nc);
		    this->object(i).resize(r);
		}
	    }
    };
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t colour_generator<value_t,N_in,N_out,true>::N_tot;

    /// Colour-summing dummy class. Signals a full colour sum of the
    /// amplitude in process/event generator instances.

    template<class value_t,std::size_t N_in,std::size_t N_out,bool cont_cols>class colour_summer;

    template<class value_t,std::size_t N_in,std::size_t N_out>class colour_summer<value_t,N_in,N_out,false>: public colour_generator<value_t,N_in,N_out,false>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef colour_generator<value_t,N_in,N_out,false> base_type;
	    typedef typename base_type::object_type object_type;
	    typedef typename base_type::size_type size_type;
	    
	    /* Public static data: */
	    /*---------------------*/
	    
	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Factory method: */

	    template<class model_t>static colour_summer<value_t,N_in,N_out,false>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		colour_summer<value_t,N_in,N_out,false>* result=new colour_summer<value_t,N_in,N_out,false>(base_type::template make_colour_vector<model_t>(it),base_type::template make_colour_range_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Constructor: */

	    colour_summer(const vector<object_type*,N_tot>& cols_,const vector<std::vector<size_type>,N_tot>& col_ranges_):base_type(cols_,col_ranges_){}

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
	
	    colour_summer<value_t,N_in,N_out,false>* clone() const
	    {
		return new colour_summer<value_t,N_in,N_out,false>(*this);
	    }

	    /* Colour summation flags set: */

	    bool sum_colour(int i) const
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
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t colour_summer<value_t,N_in,N_out,false>::N_tot;
    
    /* Colour summation class specialisation for continuous colours: */
    
    template<class value_t,std::size_t N_in,std::size_t N_out>class colour_summer<value_t,N_in,N_out,true>: public colour_generator<value_t,N_in,N_out,true>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef colour_generator<value_t,N_in,N_out,true> base_type;
	    typedef typename base_type::object_type object_type;
	    typedef typename base_type::size_type size_type;

	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Factory method: */

	    template<class model_t>static colour_summer<value_t,N_in,N_out,true>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		colour_summer<value_t,N_in,N_out,true>* result=new colour_summer<value_t,N_in,N_out,true>(base_type::template make_colour_tensor_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Constructor: */

	    colour_summer(const vector<object_type*,N_tot>& tensors):base_type(tensors){}

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

	    colour_summer<value_t,N_in,N_out,true>* clone() const
	    {
		return new colour_summer<value_t,N_in,N_out,true>(*this);
	    }

	    /* Helicity summation flags set: */

	    bool sum_colour(int i)
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
    template<class value_t,std::size_t N_in,std::size_t N_out>const std::size_t colour_summer<value_t,N_in,N_out,true>::N_tot;
}

#endif /*CAMGEN_COL_GEN_H_*/

