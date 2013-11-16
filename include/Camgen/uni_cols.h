//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file uni_cols.h
    \brief Uniform colour generators class interface and implementation.
 */

#ifndef CAMGEN_UNI_COLS_H_
#define CAMGEN_UNI_COLS_H_

#include <Camgen/uniform.h>
#include <Camgen/col_gen.h>
#include <Camgen/rn_strm.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Uniform colour generator class for discrete colour model types. *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Uniform colour generator class declaration. The first template parameter
    /// denotes the floating-point numerical type, the second and third integers
    /// resp. the initial and final state multiplicities, the fourth parameter
    /// the random number generator and the final boolean the continuous colours
    /// tag.

    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t,bool continuous_cols>class uniform_colours;
    
    /* Uniform colour generator template specialisation for discrete colours. */

    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>class uniform_colours<value_t,N_in,N_out,rng_t,false>: public colour_generator<value_t,N_in,N_out,false>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef value_t value_type;
	    typedef colour_generator<value_t,N_in,N_out,false> base_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef typename base_type::size_type size_type;

	    /* Static data: */
	    /*--------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Named constructor creating an instance from a CM algorithm tree iterator. */

	    template<class model_t>static uniform_colours<value_t,N_in,N_out,rng_t,false>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		uniform_colours<value_t,N_in,N_out,rng_t,false>* result=new uniform_colours<value_t,N_in,N_out,rng_t,false>(base_type::template make_colour_vector<model_t>(it),base_type::template make_colour_range_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Non-allocating constructor. */

	    uniform_colours(const vector<std::vector<size_type>*,N_tot>& cols_,const vector<std::vector<size_type>,N_tot>& ranges):base_type(cols_,ranges),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<ranges[i].size();++j)
		    {
			const_weight*=ranges[i][j];
		    }
		    if(ranges[i].size()>1)
		    {
			const_weight/=(ranges[i].size());
		    }
		}
	    }

	    /* Standalone-mode constructor. */

	    uniform_colours(const vector<std::vector<size_type>,N_tot>& ranges):base_type(ranges),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<ranges[i].size();++j)
		    {
			const_weight*=ranges[i][j];
		    }
		    if(ranges[i].size()>1)
		    {
			const_weight/=(ranges[i].size());
		    }
		}
	    }

	    /* Standalone-mode constructor. */

	    uniform_colours(size_type Nc,vector<size_type,N_tot>ranks):base_type(Nc,ranks),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<ranks[i];++j)
		    {
			const_weight*=Nc;
		    }
		    if(ranks[i]>1)
		    {
			const_weight/=(ranks[i]);
		    }
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generation method. Throws dice the generate the colours within
	     * the allowed ranges. The weight remains constant. */

	    bool generate()
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<this->colour_rank(i);++j)
		    {
			this->colour(i,j)=rn_stream::throw_dice(this->colour_range(i,j));
		    }
		}
		this->weight()=const_weight;
		return true;
	    }

	    /* Weight evaluation implementation, assigns the constant nr of
	     * colour states to the weight. */

	    bool evaluate_weight()
	    {
		this->weight()=const_weight;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method implementation: */

	    uniform_colours<value_t,N_in,N_out,rng_t,false>* clone() const
	    {
		return new uniform_colours<value_t,N_in,N_out,rng_t,false>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "uniform";
	    }

	private:

	    value_type const_weight;

    };
    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t uniform_colours<value_t,N_in,N_out,rng_t,false>::N_tot;
    
    /* Uniform colour generator template specialisation for continuous colours. */

    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>class uniform_colours<value_t,N_in,N_out,rng_t,true>: public colour_generator<value_t,N_in,N_out,true>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef value_t value_type;
	    typedef colour_generator<value_t,N_in,N_out,true> base_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef typename base_type::object_type object_type;
	    typedef typename base_type::size_type size_type;

	    /* Static data: */
	    /*--------------*/

	    static const std::size_t N_tot=N_in+N_out;
	    static const value_t twopi;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Named constructor creating an instance from a CM algorithm tree iterator. */

	    template<class model_t>static uniform_colours<value_t,N_in,N_out,rng_t,true>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		uniform_colours<value_t,N_in,N_out,rng_t,true>* result=new uniform_colours<value_t,N_in,N_out,rng_t,true>(base_type::template make_colour_tensor_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Non-allocating constructor. */

	    uniform_colours(const vector<object_type*,N_tot>& cols_):base_type(cols_),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    if(cols_[i]->rank()>1)
		    {
			const_weight/=(cols_[i]->rank());
		    }
		}
	    }

	    /* Allocating constructor: */

	    uniform_colours(const vector<std::vector<size_type>,N_tot>& ranges):base_type(ranges),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    if(ranges[i].size()>1)
		    {
			const_weight/=(ranges[i].size());
		    }
		}
	    }

	    /* Allocating constructor with one representation type: */

	    uniform_colours(vector<size_type,N_tot>ranks,size_type Nc):base_type(ranks),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    if(ranks[i]>1)
		    {
			const_weight/=(ranks[i]);
		    }
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generation method. Throws dice the generate the colours within
	     * the allowed ranges. The weight remains constant. */

	    bool generate()
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<this->colour_size(i);++j)
		    {
			this->colour_entry(i,j)=std::polar((value_type)1,rn_stream::throw_number(twopi));
		    }
		}
		this->weight()=const_weight;
		return true;
	    }

	    /* Weight evaluation implementation, assigns the constant nr of
	     * colour states to the weight. */

	    bool evaluate_weight()
	    {
		this->weight()=const_weight;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method implementation: */

	    uniform_colours<value_t,N_in,N_out,rng_t,true>* clone() const
	    {
		return new uniform_colours<value_t,N_in,N_out,rng_t,true>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "uniform";
	    }

	private:

	    value_type const_weight;

    };
    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t uniform_colours<value_t,N_in,N_out,rng_t,true>::N_tot;
    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>const value_t uniform_colours<value_t,N_in,N_out,rng_t,true>::twopi=2*std::acos(-(value_t)1);
}

#endif /*CAMGEN_UNI_COLS_H_*/

