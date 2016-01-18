//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SPACETIME_H_
#define CAMGEN_SPACETIME_H_

#include <vector>
#include <Camgen/debug.h>
#include <Camgen/vector.h>
#include <Camgen/product_type.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Spacetime base class template declaration and definition. The class takes 3   *
 * template parameters: the numerical type, the derived type defining the        *
 * implementation and the dimension. The class contains one static data member:  *
 * the metric, implemented as a D x D array. Also the inner product of vectors   *
 * is implemented, but usually is overloaded in the derived type for             *
 * opimisation purposes.                                                         *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Declaration and definition of the spacetime base class template: */
    
    template<class value_t,class type,std::size_t dim>class spacetime
    {
	public:
	    typedef value_t value_type;
	    static const std::size_t dimension=dim;
	    static const std::size_t index_range=dim;

	    /* Initialisation function: */

	    static void initialise()
	    {
		if(!initialised)
		{

		    /* Let the implementation subclass fill the metric tensor:
		     * */

		    type::fill_metric();
		    
		    initialised=true;

		}
	    }

	    /* Output method: */

	    static std::ostream& print(std::ostream& os)
	    {
		for(std::size_t i=0;i<dimension;++i)
		{
		    for(std::size_t j=0;j<dimension;++j)
		    {
			os<<g[i][j]<<"   ";
		    }
		    os<<std::endl;
		}
		return os;
	    }

	    /* Entries of the mteric tensor: */

	    static const value_type& metric(std::size_t i,std::size_t j)
	    {
		CAMGEN_ERROR_IF((i>=dimension),"index 1 out of range");
		CAMGEN_ERROR_IF((j>=dimension),"index 1 out of range");
		
		return g[i][j];
	    }

	    /* Inner product template member function: */

	    template<class vector_t1,class vector_t2>static typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type dot(const vector_t1& t1,const vector_t2& t2)
	    {
		typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type result=0;
		if(type::diagonal)
		{
		    for(std::size_t mu=0;mu<dimension;++mu)
		    {
			result+=(t1[mu]*g[mu][mu]*t2[mu]);
		    }
		    return result;
		}
		else
		{
		    for(std::size_t mu=0;mu<dimension;++mu)
		    {
			for(std::size_t nu=0;nu<dimension;++nu)
			{
			    result+=(t1[mu]*g[mu][nu]*t2[nu]);
			}
		    }
		    return result;
		}
	    }

	protected:
	    
	    /* Static data: the metric tensor: */

	    static value_t g[dim][dim];
	    
	private:

	    /* Initialisation tag: */

	    static bool initialised;
    };
    template<class value_t,class type,std::size_t dim>const std::size_t spacetime<value_t,type,dim>::dimension;
    template<class value_t,class type,std::size_t dim>const std::size_t spacetime<value_t,type,dim>::index_range;
    template<class value_t,class type,std::size_t dim>value_t spacetime<value_t,type,dim>::g[dim][dim]={{0}};
    template<class value_t,class type,std::size_t dim>bool spacetime<value_t,type,dim>::initialised=false;
}

#endif /*CAMGEN_SPACETIME_H_*/

