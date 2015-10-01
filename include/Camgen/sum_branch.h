//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Momentum summation branching class. Will only be constructed in 2->n  *
 * processes that contain s-type Feynman diagrams.                       *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CAMGEN_SUM_BRANCH_H_
#define CAMGEN_SUM_BRANCH_H_

#include <Camgen/ps_branching.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class sum_branching: public ps_branching<model_t,N_in,N_out,rng_t>
    {
	typedef ps_branching<model_t,N_in,N_out,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_t::value_type value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::channel_type channel_type;
	    typedef typename base_type::ps_channel_type ps_channel_type;
	    typedef typename channel_type::momentum_channel_type momentum_channel_type;

	    /* Second incoming momentum: */

	    const channel_type* const in2;

	    /* Constructor. */

	    sum_branching(channel_type* in1_,const channel_type* in2_,channel_type* out):base_type(in1_,out),in2(in2_){}

	    /* Destructor. */

	    ~sum_branching(){}

	    /* Trivial t-generation method: */

	    bool generate_t()
	    {
		return true;
	    }

	    bool generate_s()
	    {
		return true;
	    }

	    /* Generation implementation. */

	    bool generate_p()
	    {
		if(this->channel(0)->get_status()!=ps_channel_type::p_set)
		{
		    this->channel(0)->p()=this->p_in()+in2->p();
		    this->channel(0)->set_status_p_generated();
		}

		this->branching_weight=(value_type)1;
		return true;
	    }

	    /* Weight evaluation method. */

	    bool evaluate_branching_weight()
	    {
		this->branching_weight=(value_type)1;
		return true;
	    }

	    /* Returns whether the branchings are equivalent: */

	    bool equiv(const ps_branching<model_t,N_in,N_out,rng_t>* other) const
	    {
		const sum_branching<model_t,N_in,N_out,rng_t>* othercast=dynamic_cast<const sum_branching<model_t,N_in,N_out,rng_t>*>(other);
		if(othercast!=NULL)
		{
		    return (this->same_channels(other,false) and in2==othercast->in2);
		}
		return false;
	    }

	    /* Type output method. */

	    std::string type() const
	    {
		return "+";
	    }

	    momentum_type evaluate_p_in() const
	    {
		return this->p_out(0)-in2->p();
	    }
    };
    template<class spacetime_t,std::size_t N_out,class rng_t>class sum_branching<spacetime_t,1,N_out,rng_t>{};
}

#endif /*CAMGEN_SUM_BRANCH_H_*/

