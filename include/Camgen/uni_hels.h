//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file uni_hels.h
    \brief uniform helicity generator class template.
 */

#ifndef CAMGEN_UNI_HELS_H_
#define CAMGEN_UNI_HELS_H_

#include <Camgen/hel_gen.h>
#include <Camgen/rn_strm.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Uniform helicity generator classes. In the discrete case, the assignment      *
 * proceeds by dice throwing, in the continuous case by assigning random phases  *
 * to positive helicity states and their complex conjugates to their negative    *
 * counterparts, and if a longitudinal state exists, it gets coefficient 1.      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Uniform helicity generator class template declaration. */

    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t,bool q>class uniform_helicities;

    /* Uniform helicity specialisation for discrete helicities. */

    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>class uniform_helicities<value_t,N_in,N_out,rng_t,false>: public helicity_generator<value_t,N_in,N_out,false>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/
	    
	    typedef value_t value_type;
	    typedef helicity_generator<value_t,N_in,N_out,false> base_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef typename base_type::size_type size_type;

	    /* Public static data: */
	    /*---------------------*/

	    /* Static integer constants: */

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Factory method: */
	    
	    template<class model_t>static uniform_helicities<value_t,N_in,N_out,rng_t,false>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		return new uniform_helicities<value_t,N_in,N_out,rng_t,false>(base_type::template make_helicity_vector<model_t>(it),base_type::template make_max_helicity_vector<model_t>(it),base_type::template make_zero_helicity_bitset<model_t>(it));
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Non-allocating constructor (see helicity generator base class).
	     * Assigns the number of configurations to the weight. */

	    uniform_helicities(const vector<int*,N_tot>& hels_,const vector<int,N_tot>& max_hels_,std::bitset<N_tot> zero_hels_):base_type(hels_,max_hels_,zero_hels_),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    int j=max_hels_[i];
		    const_weight*=((zero_hels_[i])?(2*j+1):(2*j));
		}
	    }

	    /* Allocating constructor (see helicity generator base class).
	     * Assigns the number of configurations to the weight. */

	    uniform_helicities(const vector<int,N_tot>& max_hels_,std::bitset<N_tot> zero_hels_):base_type(max_hels_,zero_hels_),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    int j=max_hels_[i];
		    const_weight*=((zero_hels_[i])?(2*j+1):(2*j));
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generation operator implementation. Throws dice in the correct
	    ranges. Weight remains constant. */

	    bool generate()
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    int j=this->maximal_helicity(i);
		    if(j!=0)
		    {
			if(this->has_zero_helicity(i))
			{
			    this->helicity(i)=rn_stream::throw_dice(-j,j+1);
			}
			else
			{
			    this->helicity(i)=rn_stream::throw_dice(1,j+1);
			    if(rn_stream::throw_coin())
			    {
				this->helicity(i)=-this->helicity(i);
			    }
			}
		    }
		}
		this->weight()=(value_type)const_weight;
		return true;
	    }

	    /* Weight evaluation method. Assigns the constant weight. */

	    bool evaluate_weight()
	    {
		this->weight()=(value_type)const_weight;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method implementation: */

	    uniform_helicities<value_t,N_in,N_out,rng_t,false>* clone() const
	    {
		return new uniform_helicities<value_t,N_in,N_out,rng_t,false>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "uniform";
	    }

	private:

	    size_type const_weight;
    };
    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t uniform_helicities<value_t,N_in,N_out,rng_t,false>::N_tot;

    /* Uniform helicity specialisation for discrete helicities. */

    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>class uniform_helicities<value_t,N_in,N_out,rng_t,true>: public helicity_generator<value_t,N_in,N_out,true>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/
	    
	    typedef value_t value_type;
	    typedef std::complex<value_t> c_value_type;
	    typedef helicity_generator<value_t,N_in,N_out,true> base_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef typename base_type::size_type size_type;

	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Factory method: */
	    
	    template<class model_t>static uniform_helicities<value_t,N_in,N_out,rng_t,true>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		return new uniform_helicities<value_t,N_in,N_out,rng_t,true>(base_type::template make_helicity_vector<model_t>(it));
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Non-allocating constructor (see helicity generator base class).
	     * Weight remains 1. */

	    uniform_helicities(const vector<helicity_phases<value_type>*,N_tot>& hels_):base_type(hels_){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generation operator implementation. Generates complex phases for
	    * the positive helicity coefficients, the corresponding negative
	    * helicity state receive the complex conjugate phase. If a
	    * zero-helicity state exists, it gets coefficient 1. */

	    bool generate()
	    {
		for(size_type i=0;i<N_in+N_out;++i)
		{
		    int j=this->maximal_helicity(i);
		    for(int k=1;k<j+1;++k)
		    {
			c_value_type c=std::polar((value_type)1,rn_stream::throw_number(twopi));
			this->helicity_phase(i,k)=c;
			this->helicity_phase(i,-k)=std::conj(c);
		    }
		    if(this->has_zero_helicity(i))
		    {
			this->helicity_phase(i,0)=c_value_type(1,0);
		    }
		}
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Weight evaluation implementation. assigns unity to the weight
	    * variable. */

	    bool evaluate_weight()
	    {
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method implementation: */

	    uniform_helicities<value_t,N_in,N_out,rng_t,true>* clone() const
	    {
		return new uniform_helicities<value_t,N_in,N_out,rng_t,true>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "uniform";
	    }

	private:

	    static const value_type twopi;
    };
    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>const std::size_t uniform_helicities<value_t,N_in,N_out,rng_t,true>::N_tot;
    template<class value_t,std::size_t N_in,std::size_t N_out,class rng_t>const typename uniform_helicities<value_t,N_in,N_out,rng_t,true>::value_type uniform_helicities<value_t,N_in,N_out,rng_t,true>::twopi=2*std::acos(-1);
}

#endif /*CAMGEN_UNI_HELS_H_*/

