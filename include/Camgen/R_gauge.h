//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_R_GAUGE_H_
#define CAMGEN_R_GAUGE_H_

#include <vector>
#include <Camgen/width_scheme.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the R_xi gauge vector and scalar propagators.   *
 * In the current version of Camgen, all R_xi propagators in a process use the  *
 * same value of xi, defined in the class R_gauge<model_type>, and by default    *
 * set to zero. The vector propagator invoke the formula                         *
 *                                                                               *
 *   P(p)^{mu,nu} = -i(g^{mu,nu} - (1-xi)*p^{mu}p^{nu}/(p^2-xi*M^2)/(p^2-M^2)    *
 *                                                                               *
 * where M^2=m^2-i*m*Gamma, Gamma being the width. For the scalar propagator,    *
 * which is used by the would-be Goldstones of a broken symmetry, the xi appears *
 * in front of the mass squared in the denominator, making them decouple as xi   *
 * approaches infinity.                                                          *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


namespace Camgen
{
    /* Holder class of the global xi-parameter: */

    template<class model_t>class R_gauge
    {
	public:
	    static typename model_t::value_type xi;
    };

    template<class model_t>typename model_t::value_type R_gauge<model_t>::xi=0;

    /* R_xi-gauge propagator class for vector particles: */

    template<class model_t>class R_vector_gauge: public R_gauge<model_t>, public width_scheme<model_t>
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	
	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);
	    
	public:	

	    /* Propagating subamplitude size: */

	    static const size_type size=model_t::dimension;

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
	    }

	    /* Function returning the vector index ranges: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& v)
	    {
		v.clear();
		v.push_back(model_t::dimension);
		return v;
	    }
	    
	    /* Refresher function; computes the small and big denominators and stores the
	     * momentum pointer: */
	    
	    static void refresh(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		width_scheme<model_t>::evaluate(p,m,w);
	    }

	    /* Evaluate function: propagates the vector that 'first' points to: */

	    static void evaluate(iterator first)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		CAMGEN_ERROR_IF((first.range()<model_t::dimension),"tensor iterator out of range");
		value_type p2(width_scheme<model_t>::s,0);
		value_type a=value_type(1-R_gauge<model_t>::xi,0)*spacetime_type::dot(*p,first)/(p2-R_gauge<model_t>::xi*width_scheme<model_t>::gauge_mass2);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    (*first)-=a*(*p)[mu];
		    (*first)*=(-width_scheme<model_t>::denominator);
		    ++first;
		}
	    }

	    /* Evaluate function: propagates all vectors between 'first' and 'last': */

	    static void evaluate_range(iterator first,iterator last)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		CAMGEN_ERROR_IF((first.range()<model_t::dimension),"tensor iterator out of range");
		CAMGEN_ERROR_IF(((last-first)<0),"bounding iterator smaller than running iterator");
		CAMGEN_ERROR_IF(((last-first)%model_t::dimension != 0),"first and last iterator differ by non-integer multiple of vector size");
		value_type p2(width_scheme<model_t>::s,0);
		value_type r=value_type(1-R_gauge<model_t>::xi,0)/(p2-R_gauge<model_t>::xi*width_scheme<model_t>::gauge_mass2);
		value_type a;
		do
		{
		    a=spacetime_type::dot(*p,first)*r;
		    for(size_type mu=0;mu<model_t::dimension;++mu)
		    {
			(*first)-=a*(*p)[mu];
		    	(*first)*=(-width_scheme<model_t>::denominator);
			++first;
		    }
		}
		while(first!=last);
	    }
    };
    template<class model_t>const typename R_vector_gauge<model_t>::size_type R_vector_gauge<model_t>::size;

    /* R_xi-gauge propagator class for scalar particles: */

    template<class model_t>class R_scalar_gauge: public R_gauge<model_t>,public width_scheme<model_t>
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);

	public:	

	    /* Propagating subamplitude size: */

	    static const size_type size=1;

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
	    }

	    /* Propagating tensor index vector creation: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& v)
	    {
		v.clear();
		return v;
	    }
	    
	    /* Refresher function; computes the denominator: */
	    
	    static void refresh(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		width_scheme<model_t>::evaluate(p,m,w);
	    }

	    /* Evaluate function: propagates the scalar that 'first' points to: */

	    static void evaluate(iterator first)
	    {
		CAMGEN_ERROR_IF((first.range()==0),"tensor iterator out of range");
		(*first)*=(value_type(0,1)/(width_scheme<model_t>::s-R_gauge<model_t>::xi*width_scheme<model_t>::gauge_mass2));
	    }

	    /* Evaluate function: propagates all scalars between 'first' and 'last': */

	    static void evaluate_range(iterator first,iterator last)
	    {
		CAMGEN_ERROR_IF((first.range()==0),"tensor iterator out of range");
		CAMGEN_ERROR_IF(((last-first)<0),"bounding iterator smaller than running iterator");
		value_type denom=value_type(0,1)/(width_scheme<model_t>::s-R_gauge<model_t>::xi*width_scheme<model_t>::gauge_mass2);
		do
		{
		    (*first)*=denom;
		    ++first;
		}
		while(first!=last);
	    }
    };
    template<class model_t>const typename R_scalar_gauge<model_t>::size_type R_scalar_gauge<model_t>::size;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_R_GAUGE_H_*/

