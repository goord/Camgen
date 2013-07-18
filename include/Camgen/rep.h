//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_REP_H_
#define CAMGEN_REP_H_

#include <Camgen/group.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Representation base class template declaration and definition. The class      *
 * takes 4 template parameters: the numerical type, the derived class providing  *
 * a function fill_generators(), defining the generator matrices of the          *
 * representation, the dimension of the representation space and the rank of the *
 * group. The class contains a checking function, which makes sure the           *
 * generators fulfill the algebra relations defined by the groups structure      *
 * constants.                                                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t,class type,std::size_t dim,std::size_t g_rank>class representation
    {
	public:
	    
	    /* The numerical type of the matrices are generically complex numbers: */

	    typedef std::complex<value_t> value_type;
	    
	    /* Constant static integers: */

	    static const std::size_t dimension=dim;
	    static const std::size_t group_rank=g_rank;
	    static const std::size_t index_range=dim;
	    static const std::size_t multiplicity=1;

	    /* Initialisation function: */
	    
	    static void initialise()
	    {
		if(!initialised)
		{

		    /* Let the subclass type define the generators: */

		    type::fill_generators();
		    
		    /* If the representation is Hermitian, fill the appropriate
		     * offdiagonal entries: */
		    
		    if(type::hermitian)
		    {
			for(std::size_t a=0;a<group_rank;++a)
			{
			    for(std::size_t i=0;i<dimension;++i)
			    {
				for(size_t j=i+1;j<dimension;++j)
				{
				    T[a][j][i]=std::conj(T[a][i][j]);
				}
			    }
			}
		    }

		    /* Validate the algebra: */

		    check_commutators();
		    initialised=true;

		}
	    }

	    /* Output method: */

	    static std::ostream& print(std::ostream& os)
	    {
		for(std::size_t a=0;a<group_rank;++a)
		{
		    os<<"T"<<a<<"=[ ";
		    for(std::size_t j=0;j<dimension;++j)
		    {
			os<<std::left;
			os<<std::setw(2);
			shortprint(os,T[a][0][j]);
			os<<" ";
		    }
		    os<<"]"<<std::endl;
		    for(std::size_t i=1;i<dimension;++i)
		    {
			os<<std::left<<"   [ ";
			for(std::size_t j=0;j<dimension;++j)
			{
			    os<<std::left;
			    os<<std::setw(2);
			    shortprint(os,T[a][i][j]);
			    os<<" ";
			}
			os<<"]"<<std::endl;
		    }
		    os<<std::endl;
		}
		return os;
	    }

	    /* Generator entry access: */

	    static const value_type& generator(std::size_t a,std::size_t i,std::size_t j)
	    {
		CAMGEN_ERROR_IF((a>=group_rank),"group index out of range");
		CAMGEN_ERROR_IF((i>=dimension),"matrix index 1 index out of range");
		CAMGEN_ERROR_IF((j>=dimension),"matrix index 2 out of range");
		return T[a][i][j];
	    }

	    /* Test whether the class is initialised: */

	    static bool is_initialised()
	    {
		return initialised;
	    }

	protected:
	    
	    /* Static data: a length-g_rank array of generator matrices: */

	    static std::complex<value_t> T[g_rank][dim][dim];
	
	private:
	    
	    /* Initialisation tag: */
	    
	    static bool initialised;

	    /* Commutator checking function: */

	    static void check_commutators()
	    {

		/* Initialise the structure constants: */

		type::group_implementation::initialise();
		
		/* Test whether the generators fulfill the algebra: */
		
		value_type x;
		value_type y;
		for(std::size_t a=0;a<group_rank;++a)
		{
		    for(std::size_t b=0;b<a;++b)
		    {
			for(std::size_t i=0;i<dimension;++i)
			{
			    for(std::size_t j=0;j<dimension;++j)
			    {
				x=0;
				for(std::size_t k=0;k<dimension;++k)
				{
				    x+=(T[a][i][k]*T[b][k][j]);
				    x-=(T[b][i][k]*T[a][k][j]);
				}
				y=0;
				for(std::size_t c=0;c<group_rank;++c)
				{
				    y+=(value_type(0,1)*(type::group_implementation::structure_constant(a,b,c))*T[c][i][j]);

				}
				if(!equals(x,y))
				{
				    log(log_level::warning)<<CAMGEN_STREAMLOC<<"commutator algebra not fulfilled by [T"<<a<<",T"<<b<<"]("<<i<<","<<j<<")"<<endlog;
				}
			    }
			}
		    }
		}
	    }
    };
    template<class value_t,class type,std::size_t dim,std::size_t g_rank>const std::size_t representation<value_t,type,dim,g_rank>::dimension;
    template<class value_t,class type,std::size_t dim,std::size_t g_rank>const std::size_t representation<value_t,type,dim,g_rank>::group_rank;
    template<class value_t,class type,std::size_t dim,std::size_t g_rank>const std::size_t representation<value_t,type,dim,g_rank>::index_range;
    template<class value_t,class type,std::size_t dim,std::size_t g_rank>std::complex<value_t> representation<value_t,type,dim,g_rank>::T[g_rank][dim][dim]={{{0}}};
    template<class value_t,class type,std::size_t dim,std::size_t g_rank>bool representation<value_t,type,dim,g_rank>::initialised=false;	
}

#endif /*CAMGEN_REP_H_*/

