//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_EUCLIDEAN_H_
#define CAMGEN_EUCLIDEAN_H_

#include <Camgen/spacetime.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Euclidean spacetime definition. For models in dimensions >0, spacetime        *
 * typedef is required for the inner products of momenta and Lorentz vectors. A  *
 * spacetime type wraps an implementation metafunction, taking the value type    *
 * and dimension as parameters and defining the metric tensor and inner product  *
 * between real- or complex-valued vectors. Note that in the present version of  *
 * Camgen, only scalar wave function are defined in Euclidean signature         *
 * spacetimes. Calling the functions for constructing base vectors for massless  *
 * basic spinors and spin vectors for massive spinors will result in an error    *
 * message and a program abortion.                                               *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen 
{
    /* Definition of the Euclidean spacetime wrapper class: */
    
    class Euclidean_type 
    { 
	public:
	    
	    /* Implementation metafunction class of Euclidean spacetime. This is
	     * a class derived from the spacetime base class template: */
	    
	    template<class value_t,std::size_t dim>class implementation: public spacetime<value_t,implementation<value_t,dim>,dim>
	    {
		public:

		    /* Type definition referring to the base class: */
	    	
		    typedef spacetime<value_t,implementation<value_t,dim>,dim> base_type;
	    	
		    /* Static constant boolean denoting whether the metric is
		     * diagonal (for optimisation purposes): */
	    	
		    static const bool diagonal=true;
	    	
		    /* Implementation of the metric tensor (entries are by
		     * default zero): */
		    
		    static void fill_metric() 
		    {
			
			for(std::size_t i=0;i<dim;++i) 
			{
			    base_type::g[i][i]=(value_t)(1); 
			}

		    }

		    /* Separate implementation of the inner product between
		     * vectors. The type member of product_type<T1,T2> assigns
		     * the correct value type (real if T1 and T2 are both
		     * real,complex if one of them is): */

		    template<class vector_t1,class vector_t2>static typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type dot(const vector_t1& t1,const vector_t2& t2) 
		    {
			typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type result=0; 
			for(std::size_t i=0;i<dim;++i) 
			{
			    result+=t1[i]*t2[i]; 
			}
			return result; 
		    }
	    }; 
    };

    template<class value_t,std::size_t dim>const bool Euclidean_type::implementation<value_t,dim>::diagonal;

    /* For optimisation purposes, a separate specialisation of the
     * four-dimensional case, unrolling the loops: */

    template<class value_t>class Euclidean_type::implementation<value_t,4>: public spacetime<value_t,implementation<value_t,4>,4>
    {
	public:
	    typedef spacetime<value_t,implementation<value_t,4>,4> base_type;		    
	    static const bool diagonal=true;
	    static void fill_metric()
	    {
		
		base_type::g[0][0]=(value_t)(1);
		base_type::g[1][1]=(value_t)(1);
		base_type::g[2][2]=(value_t)(1);
		base_type::g[3][3]=(value_t)(1);
		
	    }
	    template<class vector_t1,class vector_t2>static typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type dot(const vector_t1& t1,const vector_t2& t2)
	    {
		return (t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2]+t1[3]*t2[3]);
	    }
    };
    
    template<class value_t>const bool Euclidean_type::implementation<value_t,4>::diagonal;
}

#endif /*CAMGEN_EUCLIDEAN_H_*/

