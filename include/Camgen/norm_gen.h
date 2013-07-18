//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file norm_gen.h
    \brief normally distributed random number generator.
 */

#ifndef CAMGEN_NORMAL_GEN_H_
#define CAMGEN_NORMAL_GEN_H_

#include <sstream>
#include <Camgen/rn_strm.h>
#include <Camgen/plt_strm.h>
#include <Camgen/MC_obj_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Normally distributed random number generator. The class template uses a polar *
 * form of the Box-Muller method, avoiding calls to trigonometric functions.     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Normally distributed random number generator class definition. The first
    /// template argument denotes the numerical type, the second parameter denotes
    /// the random number stream.
    
    template<class value_t,class rng_t>class normal_generator: public MC_object_generator<value_t,value_t,1,true>
    {
	typedef MC_object_generator<value_t,value_t,1,true> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::object_type object_type;
	    typedef typename base_type::integral_type integral_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;

	    /* Prefactor of the distribution: */

	    static const value_type factor;

	    /// Mean value address.

	    const value_type* mu;

	    /// Variance address.

	    const value_type* sigma;

	    /// Allocating constructor with mean and variance address arguments.
	    
	    normal_generator(const value_type* mu_=NULL,const value_type* sigma_=NULL):mu(mu_),sigma(sigma_),cached(false){}

	    /// Non-allocating constructor with mean and variance addresses.

	    normal_generator(object_type* instance_,const value_type* mu_=NULL,const value_type* sigma_=NULL):base_type(instance_),mu(mu_),sigma(sigma_),cached(false){}

	    /// Cloning method.

	    normal_generator<value_t,rng_t>* clone() const
	    {
		return new normal_generator<value_t,rng_t>(*this);
	    }

	    /// Throwing operator.

	    bool generate()
	    {
		if(cached)
		{
		    this->object()=cache;
		    this->weight()=cached_weight;
		    cached=false;
		}
		else
		{
		    value_type w,x1,x2;
		    do
		    {
			x1=rn_stream::throw_number(-(value_type)1,(value_type)1);
			x2=rn_stream::throw_number(-(value_type)1,(value_type)1);
			w=x1*x1+x2*x2;
		    }
		    while(w>=(value_type)1);
		    value_type y=(sigma==NULL)?(std::sqrt(-(value_type)2*std::log(w)/w)):(std::sqrt(-(value_type)2*std::log(w)/w)*(*sigma));
		    cache=(mu==NULL)?(x1*y):(*mu+x1*y);
		    cached_weight=(sigma==NULL)?(std::pow(w,-x1*x1/w)/factor):(std::abs(*sigma)*std::pow(w,-x1*x1/w)/factor);
		    this->object()=(mu==NULL)?(x2*y):(*mu+x2*y);
		    this->weight()=(sigma==NULL)?(std::pow(w,-x2*x2/w)/factor):(std::abs(*sigma)*std::pow(w,-x2*x2/w)/factor);
		    cached=true;
		}
		return true;
	    }

	    /// Evaluates the weight of the current object.

	    bool evaluate_weight()
	    {
		value_type z=(mu==NULL)?(this->object()):(this->object()-*mu);
		value_type c;
		if(sigma!=NULL)
		{
		    c=factor/std::abs(*sigma);
		    z/=(*sigma);
		}
		else
		{
		    c=factor;
		}
		z*=((value_type)0.5*z);
		this->weight()=std::exp(z)/c;
		return true;
	    }

	    value_type mean() const
	    {
		return (mu==NULL)?(value_type)0:(*mu);
	    }

	    value_type variance() const
	    {
		return (sigma==NULL)?(value_type)1:(*sigma);
	    }

	    function_stream* pdf_plot(const value_type& x) const
	    {
		value_type norm=x*factor/std::abs(variance());
		std::ostringstream sstrm;
		if(mu==NULL)
		{
		    sstrm<<norm<<"*exp(-x**2";
		}
		else
		{
		    sstrm<<norm<<"*exp(-(x-"<<*mu<<")**2";
		}
		value_type denom=(value_type)2*std::pow(variance(),(int)2);
		sstrm<<"/"<<denom<<")";
		function_stream* result=new function_stream(sstrm.str());
		result->title="density";
		result->style="lines";
		return result;
	    }

	private:

	    value_type cache,cached_weight;
	    
	    /* Bool denoting whether the cache is full: */
	    
	    bool cached;
    };

    template<class value_t,class rng_t>const typename normal_generator<value_t,rng_t>::value_type normal_generator<value_t,rng_t>::factor=(value_t)1/std::sqrt((value_t)2*std::acos(-(value_t)1));
}

#endif /*CAMGEN_NORMAL_GEN_H_*/

