//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef GROUP_H_
#define GROUP_H_

#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/c_utils.h>
#include <vector>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition and declaration of the group base class template. The class has  *
 * only one static data member: the structure constants array f. It can be     *
 * filled by derived class provideing an implementation of a specific group.   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Declaration and definition of the group class template. It takes three
     * template parameters: the numerical type, the derived type providing the
     * implementation and the rank: */

    template<class value_t,class type,std::size_t N>class group
    {
	public:

	    /* Some useful type definitions and static constants: */

	    typedef std::complex<value_t> value_type;
	    static const std::size_t rank=N;
	    static const std::size_t index_range=N;

	    /* Initialisation function: */

	    static void initialise()
	    {
		if(!initialised)
		{

		    /* Let the derived class fill the necessary structure
		     * constants: */
		    
		    type::fill_structure_constants();
		    
		    /* Use the antisymmetry properties to give values to the
		     * rest of the constants: */
		    
		    if(type::antisymmetric)
		    {
			for(std::size_t a=0;a<rank;++a)
			{
			    for(std::size_t b=a+1;b<rank;++b)
			    {
				for(std::size_t c=b+1;c<rank;++c)
				{
				    if(f[a][b][c] != (value_type)0)
				    {
					f[b][a][c]=-f[a][b][c];
					f[c][b][a]=-f[a][b][c];
					f[a][c][b]=-f[a][b][c];
					f[c][a][b]=f[a][b][c];
					f[b][c][a]=f[a][b][c];
				    }
				}
			    }
			}
		    }
		    else
		    {
			for(std::size_t c=0;c<rank;++c)
			{
			    for(std::size_t a=0;a<rank;++a)
			    {
				for(std::size_t b=0;b<a;++b)
				{
				    if(f[a][b][c] != (value_type)0)
				    {
					f[b][a][c]=-f[a][b][c];
				    }
				}
			    }
			}
		    }
		    check_Jacobi();
		    initialised=true;

		}
	    }

	    /* Output method: */

	    static std::ostream& print(std::ostream& os)
	    {
		os<<"commutators:"<<std::endl<<std::endl;
		for(std::size_t i=0;i<rank;++i)
		{
		    for(std::size_t j=i+1;j<rank;++j)
		    {
			os<<"[ T_"<<i<<" , "<<"T_"<<j<<" ] = ";
			std::vector<value_type>v;
			std::vector<std::size_t>I;
			for(std::size_t k=0;k<rank;++k)
			{
			    if(f[i][j][k] != (value_type)0)
			    {
				I.push_back(k);
				v.push_back(value_type(0,1)*f[i][j][k]);
			    }
			}
			if(v.size()>0)
			{
			    shortprint_first_prefactor(os,v[0]);
			    os<<"T_"<<I[0];
			    for(std::size_t k=1;k<v.size();++k)
			    {
				shortprint_prefactor(os,v[k]);
				os<<"T_"<<I[k];
			    }
			}
			else
			{
			    os<<"0";
			}
			os<<std::endl;
		    }
		}
		return os;
	    }

	    /* Access to the structure constants: */

	    static const value_type& structure_constant(std::size_t a,std::size_t b,std::size_t c)
	    {
		CAMGEN_ERROR_IF((a>=N),"index 1 out of range");
		CAMGEN_ERROR_IF((b>=N),"index 1 out of range");
		CAMGEN_ERROR_IF((c>=N),"index 1 out of range");
		return f[a][b][c];
	    }

	    /* Test whether the class is initialised: */

	    static bool is_initialised()
	    {
		return initialised;
	    }

	protected:

	    /* Static structure constants data: */

	    static std::complex<value_t>f[N][N][N];
	private:

	    /* Initialisation tag: */

	    static bool initialised;
	    
	    /* Jacobi-identity checking module: */
	    
	    static void check_Jacobi()
	    {
		value_type x;
		for(std::size_t a=0;a<rank;++a)
		{
		    for(std::size_t b=0;b<rank;++b)
		    {
			for(std::size_t c=0;c<rank;++c)
			{
			    for(std::size_t e=0;e<rank;++e)
			    {
				x=0;
				for(std::size_t d=0;d<rank;++d)
				{
				    x+=(f[a][d][e]*f[b][c][d]);
				    x+=(f[b][d][e]*f[c][a][d]);
				    x+=(f[c][d][e]*f[a][b][d]);
				}
				if(!equals(x.real(),(value_t)0) or !equals(x.imag(),(value_t)0))
				{
				    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Jacobi identity not fulfilled"<<"[T"<<a<<",[T"<<b<<",T"<<c<<"]] + "<<"[T"<<b<<",[T"<<c<<",T"<<a<<"]] + "<<"[T"<<c<<",[T"<<a<<",T"<<b<<"]]="<<x<<"*T"<<e<<"+..."<<endlog;
				}
			    }
			}
		    }
		}
	    }
    };

    template<class value_t,class type,std::size_t N>const std::size_t group<value_t,type,N>::rank;
    template<class value_t,class type,std::size_t N>const std::size_t group<value_t,type,N>::index_range;
    template<class value_t,class type,std::size_t N>std::complex<value_t>group<value_t,type,N>::f[N][N][N]={{{0}}};
    template<class value_t,class type,std::size_t N>bool group<value_t,type,N>::initialised=false;
}

#endif /*CAMGEN_GROUP_H_*/

