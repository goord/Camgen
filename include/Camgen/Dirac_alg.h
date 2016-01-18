//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DIRAC_ALG_H_
#define CAMGEN_DIRAC_ALG_H_

#include <Camgen/c_utils.h>
#include <Camgen/Dirac_dim.h>
#include <Camgen/vector.h>
#include <Camgen/tensor.h>
#include <Camgen/spacetime.h>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Base class of Dirac algebras. The Dirac matrices are implemented using the    *
 * "curiously recurring template pattern", where the derived class is a          *
 * template parameter. Therefore, a derived Dirac basis class only needs an      *
 * implementation of the fill_gamma_matrices(), fill_gamma_5() and               *
 * fill_C_matrix() members to obtain all recursive relations and bilinears in    *
 * this basis. However, such methods can also be overloaded to gain performance. *
 * In that case the Dirac_algebra class contains member functions that check if  *
 * the overloading implementation is correct.                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Preprocessor macro definitions to make the code shorter...*/

#define DIAGONAL type::spacetime_type::diagonal

#define METRIC(mu,nu) type::spacetime_type::metric(mu,mu)

#define S_FIRST(matrix,C_0,A_0,A_1,A_2)												\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_0[0]+=C_0*matrix[i][j]*A_1[i]*A_2[j];											\
    }																\
}

#define S_FIRST2(matrix0,matrix1,C_0,C_1,A_0,A_1,A_2)										\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_0[0]+=(C_0*matrix0[i][j]+C_1*matrix1[i][j])*A_1[i]*A_2[j];								\
    }																\
}

#define S_SECOND(matrix,C_0,A_0,A_1,A_2)											\
value_type z=C_0*A_0[0];													\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_1[i]+=(z*matrix[i][j]*A_2[j]);											\
    }																\
}

#define S_SECOND2(matrix0,matrix1,C_0,C_1,A_0,A_1,A_2)										\
value_type z0=C_0*A_0[0];													\
value_type z1=C_1*A_0[0];													\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_1[i]+=(z0*matrix0[i][j]+z1*matrix1[i][j])*A_2[j];									\
    }																\
}

#define S_THIRD(matrix,C_0,A_0,A_1,A_2)												\
value_type z=C_0*A_0[0];													\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_2[i]+=(z*matrix[j][i]*A_1[j]);											\
    }																\
}

#define S_THIRD2(matrix0,matrix1,C_0,C_1,A_0,A_1,A_2)										\
value_type z0=C_0*A_0[0];													\
value_type z1=C_1*A_0[0];													\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_2[i]+=(z0*matrix0[j][i]+z1*matrix1[j][i])*A_1[j];									\
    }																\
}

#define V_FIRST(matrix,C_0,A_0,A_1,A_2)												\
for(size_type mu=0;mu<spacetime_dimension;++mu)											\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    A_0[mu]+=C_0*matrix[mu][i][j]*A_1[i]*A_2[j];									\
	}															\
    }																\
}

#define V_FIRST2(matrix0,matrix1,C_0,C_1,A_0,A_1,A_2)										\
for(size_type mu=0;mu<spacetime_dimension;++mu)											\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    A_0[mu]+=(C_0*matrix0[mu][i][j]+C_1*matrix1[mu][i][j])*A_1[i]*A_2[j];						\
	}															\
    }																\
}

#define V_SECOND(matrix,C_0,A_0,A_1,A_2)											\
if(DIAGONAL)															\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		Vslash[i][j]+=matrix[mu][i][j]*(METRIC(mu,mu))*A_0[mu];								\
	    }															\
	}															\
    }																\
}																\
else																\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<spacetime_dimension;++nu)									\
		{														\
		    Vslash[i][j]+=matrix[mu][i][j]*(METRIC(mu,nu))*A_0[nu];							\
		}														\
	    }															\
	}															\
    }																\
}																\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_1[i]+=C_0*Vslash[i][j]*A_2[j];											\
    }																\
}

#define V_SECOND2(M_0,M_1,C_0,C_1,A_0,A_1,A_2)											\
if(DIAGONAL)															\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		Vslash[i][j]+=(C_0*M_0[mu][i][j]+C_1*M_1[mu][i][j])*(METRIC(mu,mu))*A_0[mu];					\
	    }															\
	}															\
    }																\
}																\
else																\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<spacetime_dimension;++nu)									\
		{														\
		    Vslash[i][j]+=(C_0*M_0[mu][i][j]+C_1*M_1[mu][i][j])*(METRIC(mu,nu))*A_0[nu];				\
		}														\
	    }															\
	}															\
    }																\
}																\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_1[i]+=Vslash[i][j]*A_2[j];												\
    }																\
}

#define V_THIRD(matrix,C_0,A_0,A_1,A_2)												\
if(DIAGONAL)															\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		Vslash[i][j]+=matrix[mu][i][j]*(METRIC(mu,mu))*A_0[mu];								\
	    }															\
	}															\
    }																\
}																\
else																\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<spacetime_dimension;++nu)									\
		{														\
		    Vslash[i][j]+=matrix[mu][i][j]*(METRIC(mu,nu))*A_0[nu];							\
		}														\
	    }															\
	}															\
    }																\
}																\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_2[i]+=C_0*A_1[j]*Vslash[j][i];											\
    }																\
}

#define V_THIRD2(M_0,M_1,C_0,C_1,A_0,A_1,A_2)											\
if(DIAGONAL)															\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		Vslash[i][j]+=(C_0*M_0[mu][i][j]+C_1*M_1[mu][i][j])*(METRIC(mu,mu))*A_0[mu];					\
	    }															\
	}															\
    }																\
}																\
else																\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<spacetime_dimension;++nu)									\
		{														\
		    Vslash[i][j]+=(C_0*M_0[mu][i][j]+C_1*M_1[mu][i][j])*(METRIC(mu,nu))*A_0[nu];				\
		}														\
	    }															\
	}															\
    }																\
}																\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_2[i]+=Vslash[j][i]*A_1[j];												\
    }																\
}

#define T_FIRST(matrix,C_0,A_0,A_1,A_2)												\
size_type n=0;															\
for(size_type mu=0;mu<spacetime_dimension;++mu)											\
{																\
    for(size_type nu=0;nu<spacetime_dimension;++nu)										\
    {																\
	for(size_type i=0;i<index_range;++i)											\
	{															\
	    for(size_type j=0;j<index_range;++j)										\
	    {															\
		A_0[n]+=C_0*matrix[mu][nu][i][j]*A_1[i]*A_2[j];									\
	    }															\
	}															\
	++n;															\
    }																\
}

#define T_SECOND(matrix,C_0,A_0,A_1,A_2)											\
size_type n;															\
if(DIAGONAL)															\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    n=0;														\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<spacetime_dimension;++nu)									\
		{														\
		    Vslash[i][j]+=matrix[mu][nu][i][j]*(METRIC(mu,mu))*(METRIC(nu,nu))*A_0[n];					\
		    ++n;													\
		}														\
	    }															\
	}															\
    }																\
}																\
else																\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    n=0;														\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<nu;++nu)											\
		{														\
		    n=0;													\
		    for(size_type sig=0;sig<spacetime_dimension;++sig)								\
		    {														\
			for(size_type rho=0;rho<spacetime_dimension;++rho)							\
			{													\
			    Vslash[i][j]+=matrix[mu][nu][i][j]*(METRIC(mu,sig))*(METRIC(nu,rho))*A_0[n];			\
			    ++n;												\
			}													\
		    }														\
		}														\
	    }															\
	}															\
    }																\
}																\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_1[i]+=C_0*Vslash[i][j]*A_2[j];											\
    }																\
}															

#define T_THIRD(matrix,C_0,A_0,A_1,A_2)												\
size_type n;															\
if(DIAGONAL)															\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    n=0;														\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<spacetime_dimension;++nu)									\
		{														\
		    Vslash[i][j]+=matrix[mu][nu][i][j]*(METRIC(mu,mu))*(METRIC(nu,nu))*A_0[n];					\
		    ++n;													\
		}														\
	    }															\
	}															\
    }																\
}																\
else																\
{																\
    for(size_type i=0;i<index_range;++i)											\
    {																\
	for(size_type j=0;j<index_range;++j)											\
	{															\
	    n=0;														\
	    Vslash[i][j]=0;													\
	    for(size_type mu=0;mu<spacetime_dimension;++mu)									\
	    {															\
		for(size_type nu=0;nu<spacetime_dimension;++nu)									\
		{														\
		    n=0;													\
		    for(size_type sig=0;sig<spacetime_dimension;++sig)								\
		    {														\
			for(size_type rho=0;rho<spacetime_dimension;++rho)							\
			{													\
			    Vslash[i][j]+=matrix[mu][nu][i][j]*(METRIC(mu,sig))*(METRIC(nu,rho))*A_0[n];			\
			    ++n;												\
			}													\
		    }														\
		}														\
	    }															\
	}															\
    }																\
}																\
for(size_type i=0;i<index_range;++i)												\
{																\
    for(size_type j=0;j<index_range;++j)											\
    {																\
	A_2[i]+=C_0*Vslash[j][i]*A_1[j];											\
    }																\
}

#define CHECK_SFF1(mat,gen)													\
fill_scal(gen);															\
check=scal;															\
type::mat##_first(c1,scal_iter,s1_iter,s2_iter);										\
mat##_first(c1,check_iter,s1_iter,s2_iter);											\
if(!equal_sequences(scal,check))												\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s1;															\
type::mat##_second(c1,scal_iter,s1_iter,s2_iter);										\
mat##_second(c1,scal_iter,check_iter,s2_iter);											\
if(!equal_sequences(s1,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s2;															\
type::mat##_third(c1,scal_iter,s1_iter,s2_iter);										\
mat##_third(c1,scal_iter,s1_iter,check_iter);											\
if(!equal_sequences(s2,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"third recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
return true

#define CHECK_SFF2(mat,gen)													\
fill_scal(gen);															\
check=scal;															\
type::mat##_first(c1,c2,scal_iter,s1_iter,s2_iter);										\
mat##_first(c1,c2,check_iter,s1_iter,s2_iter);											\
if(!equal_sequences(scal,check))												\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s1;															\
type::mat##_second(c1,c2,scal_iter,s1_iter,s2_iter);										\
mat##_second(c1,c2,scal_iter,check_iter,s2_iter);										\
if(!equal_sequences(s1,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s2;															\
type::mat##_third(c1,c2,scal_iter,s1_iter,s2_iter);										\
mat##_third(c1,c2,scal_iter,s1_iter,check_iter);										\
if(!equal_sequences(s2,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"third recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
return true

#define CHECK_VFF1(mat,gen)													\
fill_vec(gen);															\
check=vec;															\
type::mat##_first(c1,vec_iter,s1_iter,s2_iter);											\
mat##_first(c1,check_iter,s1_iter,s2_iter);											\
if(!equal_sequences(vec,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s1;															\
type::mat##_second(c1,vec_iter,s1_iter,s2_iter);										\
mat##_second(c1,vec_iter,check_iter,s2_iter);											\
if(!equal_sequences(s1,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s2;															\
type::mat##_third(c1,vec_iter,s1_iter,s2_iter);											\
mat##_third(c1,vec_iter,s1_iter,check_iter);											\
if(!equal_sequences(s2,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"third recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
return true

#define CHECK_VFF2(mat,gen)													\
fill_vec(gen);															\
check=vec;															\
type::mat##_first(c1,c2,vec_iter,s1_iter,s2_iter);										\
mat##_first(c1,c2,check_iter,s1_iter,s2_iter);											\
if(!equal_sequences(vec,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s1;															\
type::mat##_second(c1,c2,vec_iter,s1_iter,s2_iter);										\
mat##_second(c1,c2,vec_iter,check_iter,s2_iter);										\
if(!equal_sequences(s1,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s2;															\
type::mat##_third(c1,c2,vec_iter,s1_iter,s2_iter);										\
mat##_third(c1,c2,vec_iter,s1_iter,check_iter);											\
if(!equal_sequences(s2,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"third recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
return true

#define CHECK_TFF1(mat,gen)													\
fill_tens(gen);															\
check=tens;															\
type::mat##_first(c1,tens_iter,s1_iter,s2_iter);										\
mat##_first(c1,check_iter,s1_iter,s2_iter);											\
if(!equal_sequences(tens,check))												\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s1;															\
type::mat##_second(c1,tens_iter,s1_iter,s2_iter);										\
mat##_second(c1,tens_iter,check_iter,s2_iter);											\
if(!equal_sequences(s1,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s2;															\
type::mat##_third(c1,tens_iter,s1_iter,s2_iter);										\
mat##_third(c1,tens_iter,s1_iter,check_iter);											\
if(!equal_sequences(s2,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"third recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
return true

#define CHECK_TFF2(mat,gen)													\
fill_tens(gen);															\
check=tens;															\
type::mat##_first(c1,c2,tens_iter,s1_iter,s2_iter);										\
mat##_first(c1,c2,check_iter,s1_iter,s2_iter);											\
if(!equal_sequences(tens,check))												\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s1;															\
type::mat##_second(c1,c2,tens_iter,s1_iter,s2_iter);										\
mat##_second(c1,c2,tens_iter,check_iter,s2_iter);										\
if(!equal_sequences(s1,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
check=s2;															\
type::mat##_third(c1,c2,tens_iter,s1_iter,s2_iter);										\
mat##_third(c1,c2,tens_iter,s1_iter,check_iter);										\
if(!equal_sequences(s2,check))													\
{																\
    log(log_level::warning)<<CAMGEN_STREAMLOC<<"third recursive relation not correctly overloaded for "<<#mat<<endlog;		\
    return false;														\
}																\
return true

namespace Camgen
{
    /* Definition of the Dirac algebra base class template. The first template
     * parameter is the numerical type, the second parameter the subtype that
     * defines the implementation of the matrices and the third parameter is a
     * boolean denoting whether the spacetime dimension is even or odd: */

    template<class value_t,class type,std::size_t dim,bool q=(dim%2==0)>class Dirac_algebra;

    /* Definition of the Dirac algebra base class template in the case of even
     * dimensions, where there exists a chiral matrix: */

    template<class value_t,class type,std::size_t dim>class Dirac_algebra<value_t,type,dim,true>
    {
	public:

	    /* Some convenient type definitions: */

	    typedef value_t r_value_type;
	    typedef std::complex<value_t> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef vector<value_t,dim> momentum_type;
	    typedef typename tensor_type::size_type size_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    
	    /* Definitions of the matrix size (index_range) and the spacetime
	     * dimension: */
	    
	    static const std::size_t index_range=Dirac_dim<dim>::value;
	    static const std::size_t spacetime_dimension=dim;

	    /* Initialisation function: */

	    static void initialise()
	    {
		if(!initialised)
		{

		    /* Let the derived type's spacetime type fill the metric
		     * tensor:*/
		    
		    type::spacetime_type::initialise();
		    
		    /* Let the derived type fill the gamma matrices, the chiral
		     * matrix and charge conjugation matrix: */
		    
		    type::fill_gamma_matrices();
		    type::fill_gamma_5();
		    type::fill_C_matrix();
		    type::fill_D_matrix();

		    /* Construction of the left and right projectors and the
		     * required products with the charge conjugation matrix: */

		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    for(size_type k=0;k<index_range;++k)
			    {
				g_5_C[i][j]+=g_5[i][k]*C_matrix[k][j];
				Cc_g_5[i][j]+=std::conj(C_matrix[i][k])*g_5[k][j];
			    }
			    if(i==j)
			    {
				g_L[i][j]=(value_type(1,0)+g_5[i][j])/value_type(2,0);
				g_R[i][j]=(value_type(1,0)-g_5[i][j])/value_type(2,0);
			    }
			    else
			    {
				g_L[i][j]=(g_5[i][j])/value_type(2,0);
				g_R[i][j]=-(g_5[i][j])/value_type(2,0);
			    }
			    g_L_C[i][j]=(C_matrix[i][j]+g_5_C[i][j])/value_type(2,0);
			    g_R_C[i][j]=(C_matrix[i][j]-g_5_C[i][j])/value_type(2,0);
			    Cc_g_L[i][j]=(std::conj(C_matrix[i][j])+Cc_g_5[i][j])/value_type(2,0);
			    Cc_g_R[i][j]=(std::conj(C_matrix[i][j])-Cc_g_5[i][j])/value_type(2,0);
			}
		    }

		    /* Construction of the vector-current matrices and their
		     * products with the charge conjugation matrix: */

		    for(size_type mu=0;mu<spacetime_dimension;++mu)
		    {
			for(size_type i=0;i<index_range;++i)
			{
			    for(size_type j=0;j<index_range;++j)
			    {
				for(size_type k=0;k<index_range;++k)
				{
				    gg_5[mu][i][j]+=(g[mu][i][k]*g_5[k][j]);
				    g_C[mu][i][j]+=(g[mu][i][k]*C_matrix[k][j]);
				    Cc_g[mu][i][j]+=(std::conj(C_matrix[i][k])*g[mu][k][j]);
				}
				gg_L[mu][i][j]=(g[mu][i][j]+gg_5[mu][i][j])/value_type(2,0);
				gg_R[mu][i][j]=(g[mu][i][j]-gg_5[mu][i][j])/value_type(2,0);
			    }
			}
			for(size_type i=0;i<index_range;++i)
			{
			    for(size_type j=0;j<index_range;++j)
			    {
				for(size_type k=0;k<index_range;++k)
				{
				    gg_5_C[mu][i][j]+=(gg_5[mu][i][k]*C_matrix[k][j]);
				    Cc_gg_5[mu][i][j]+=(std::conj(C_matrix[i][k])*gg_5[mu][k][j]);
				}
				gg_L_C[mu][i][j]=(g_C[mu][i][j]+gg_5_C[mu][i][j])/value_type(2,0);
				gg_R_C[mu][i][j]=(g_C[mu][i][j]-gg_5_C[mu][i][j])/value_type(2,0);
				Cc_gg_L[mu][i][j]=(Cc_g[mu][i][j]+Cc_gg_5[mu][i][j])/value_type(2,0);
				Cc_gg_R[mu][i][j]=(Cc_g[mu][i][j]-Cc_gg_5[mu][i][j])/value_type(2,0);
			    }
			}
		    }

		    /* Construction of the tensor current matrices, i.e. the
		     * commutators of Dirac matrices, and their products with
		     * the charge conjugation matrix: */

		    for(size_type mu=0;mu<spacetime_dimension;++mu)
		    {
			for(size_type nu=0;nu<spacetime_dimension;++nu)
			{
			    for(size_type i=0;i<index_range;++i)
			    {
				for(size_type j=0;j<index_range;++j)
				{
				    for(size_type k=0;k<index_range;++k)
				    {
					g_comms[mu][nu][i][j]+=g[mu][i][k]*g[nu][k][j];
					g_comms[mu][nu][i][j]-=g[nu][i][k]*g[mu][k][j];
				    }
				}
			    }
			    for(size_type i=0;i<index_range;++i)
			    {
				for(size_type j=0;j<index_range;++j)
				{
				    for(size_type k=0;k<index_range;++k)
				    {
					g_comms_C[mu][nu][i][j]+=g_comms[mu][nu][i][k]*C_matrix[k][j];
					Cc_g_comms[mu][nu][i][j]+=std::conj(C_matrix[i][k])*g_comms[mu][nu][k][j];
				    }
				}
			    }
			}
		    }
		    initialised=true;
		}
	    }

	    /* Default implementation of the Dirac conjugation matrix, equal to
	     * the zeroth gamma matrix: */

	    static void fill_D_matrix()
	    {
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			D_matrix[i][j]=g[0][i][j];
		    }
		}
	    }

	    /* Check whether the gamma matrices fulfill the Clifford algebra: */
	    
	    static bool check_anticommutators()
	    {
		value_type c;
		for(size_type mu=0;mu<spacetime_dimension;++mu)
		{
		    for(size_type nu=0;nu<spacetime_dimension;++nu)
		    {
			value_t gmn=type::spacetime_type::metric(mu,nu);
			for(size_type i=0;i<index_range;++i)
			{
			    for(size_type j=0;j<index_range;++j)
			    {
				c=0;
				for(size_type k=0;k<index_range;++k)
				{
				    c+=g[mu][i][k]*g[nu][k][j];
				    c+=g[nu][i][k]*g[mu][k][j];
				}
				if(i==j)
				{
				    if(!equals(c,value_type(2*gmn,0)))
				    {
					log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac algebra not fulfilled: {g"<<mu<<",g"<<nu<<"}("<<i<<","<<j<<")="<<c<<endlog;
					return false;
				    }
				}
				else
				{
				    if(!equals(c,value_type(0,0)))
				    {
					log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac algebra not fulfilled: {g"<<mu<<",g"<<nu<<"}("<<i<<","<<j<<")="<<c<<endlog;
					return false;
				    }
				}
			    }
			}
		    }
		}
		return true;
	    }

	    /* Hermiticity and algebra properties of the chiral Dirac matrix
	     * checking function: */

	    static bool check_gamma_5()
	    {
		value_type c;

		/* Checking whether g5 squares to one: */

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			c=0;
			for(size_type k=0;k<index_range;++k)
			{
			    c+=g_5[i][k]*g_5[k][j];
			}
			if(i==j)
			{
			    if(!equals(c,value_type(1,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac algebra not fulfilled: g5^2("<<i<<","<<j<<")="<<c<<endlog;
			    	return false;
			    }
			}
			else
			{
			    if(!equals(c,value_type(0,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac algebra not fulfilled: g5^2("<<i<<","<<j<<")="<<c<<endlog;
			    	return false;
			    }
			}
		    }
		}

		/* Checking whether g5 is Hermitian: */

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<=i;++j)
		    {
			if(!equals(g_5[i][j],std::conj(g_5[j][i])))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"given g5 is not hermitian: g5("<<i<<","<<j<<") = "<<g_5[i][j]<<endlog;
			    return false;
			}
		    }
		}

		/* Checking whether g5 anticommutes with all gamma matrices: */

		for(size_type mu=0;mu<spacetime_dimension;++mu)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    c=0;
			    for(size_type k=0;k<index_range;++k)
			    {
				c+=g[mu][i][k]*g_5[k][j];
				c+=g_5[i][k]*g[mu][k][j];
			    }
			    if(!equals(c,value_type(0,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac algebra not fulfilled: {g5,g"<<mu<<"}("<<i<<","<<j<<")="<<c<<endlog;
				return false;
			    }
			}
		    }
		}
		return true;
	    }

	    /* Checking unitarity of the charge conjugation matrix: */

	    static bool check_C_matrix()
	    {
		value_type c=0;
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			c=0;
			for(size_type k=0;k<index_range;++k)
			{
			    c+=(std::conj(C_matrix[k][i])*C_matrix[k][j]);

			}
			if(i==j)
			{
			    if(!equals(c,value_type(1,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"conjugation matrix is not unitary: C^+C("<<i<<","<<j<<")="<<c<<endlog;
			    	return false;
			    }
			}
			else
			{
			    if(!equals(c,value_type(0,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"conjugation matrix is not unitary: C^+C("<<i<<","<<j<<")="<<c<<endlog;
			    	return false;
			    }
			}
		    }
		}
		return true;
	    }
	    static bool check_D_matrix()
	    {

	  	/* Checking Hermiticity of the Dirac conjugation matrix: */

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<=i;++j)
		    {
			if(!equals(D_matrix[j][i],std::conj(D_matrix[i][j])))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac conjugation matrix is not hermitian: D("<<i<<","<<j<<") = "<<D_matrix[i][j]<<endlog;
			    return false;
			}
		    }
		}

		/* Checking if the gamma matrices are self-adjoint under Dirac
		 * conjugation: */

		value_type c;
		for(size_type mu=0;mu<spacetime_dimension;++mu)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    c=0;
			    for(size_type k=0;k<index_range;++k)
			    {
				c+=g[mu][i][k]*D_matrix[k][j];
				c-=D_matrix[i][k]*std::conj(g[mu][j][k]);
			    }
			    if(!equals(c,value_type(0,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac matrix g"<<mu<<" detected that is not self-adjoint under Dirac conjugation..."<<endlog;
				return false;
			    }
			}
		    }
		}
		return true;
	    }
	    
	    /* Output method: */

	    static std::ostream& print(std::ostream& os)
	    {
		for(size_type mu=0;mu<spacetime_dimension;++mu)
		{
		    os<<"g"<<mu<<"=[ ";
		    for(size_type j=0;j<index_range;++j)
		    {
			os<<std::left;
			os<<std::setw(2);
			shortprint(os,g[mu][0][j]);
			os<<" ";
		    }
		    os<<"]"<<std::endl;
		    for(size_type i=1;i<index_range;++i)
		    {
			os<<std::left<<"   [ ";
			for(size_type j=0;j<index_range;++j)
			{
			    os<<std::left;
			    os<<std::setw(2);
			    shortprint(os,g[mu][i][j]);
			    os<<" ";
			}
			os<<"]"<<std::endl;
		    }
		    os<<std::endl;
		}
		os<<"g5=[ ";
		for(size_type j=0;j<index_range;++j)
		{
		    os<<std::left;
		    os<<std::setw(2);
		    shortprint(os,g_5[0][j]);
		    os<<" ";
		}
		os<<"]"<<std::endl;
		for(size_type i=1;i<index_range;++i)
		{
		    os<<std::left<<"   [ ";
		    for(size_type j=0;j<index_range;++j)
		    {
			os<<std::left;
			os<<std::setw(2);
			shortprint(os,g_5[i][j]);
			os<<" ";
		    }
		    os<<"]"<<std::endl;
		}
		os<<std::endl;
		os<<"C= [ ";
		for(size_type j=0;j<index_range;++j)
		{
		    os<<std::left;
		    os<<std::setw(2);
		    shortprint(os,C_matrix[0][j]);
		    os<<" ";
		}
		os<<"]"<<std::endl;
		for(size_type i=1;i<index_range;++i)
		{
		    os<<std::left<<"   [ ";
		    for(size_type j=0;j<index_range;++j)
		    {
			os<<std::left;
			os<<std::setw(2);
			shortprint(os,C_matrix[i][j]);
			os<<" ";
		    }
		    os<<"]"<<std::endl;
		}
		return os;
	    }

	    /* Entries of the unit matrix: */

	    static value_type Id(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		if(i!=j)
		{
		    return value_type(1,0);
		}
		return value_type(0,0);
	    }
	    
	    /* Basic spinor bilinear without prefactor x += bar{psi1}.psi2: */
	    
	    static void Id_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    x+=A_1[i]*A_2[i];
		}
	    }
	    
	    /* Basic spinor bilinear without prefactor x += bar{psi1}.psi2: */
	    
	    static void Id_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    (*A_0)+=C_0*A_1[i]*A_2[i];
		}
	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*psi2 */

	    static void Id_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    A_1[i]+=c*A_2[i];
		}
	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*phi*bar{psi1}: */

	    static void Id_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    A_2[i]+=c*A_1[i];
		}
	    }

	    /* Function checking whether the Id-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Id(generator& gen)
	    {
		fill_scal(gen);
		value_type x1(0,0);
		value_type x2(0,0);
		Id_first(x1,s1_iter,s2_iter);
		type::Id_first(x2,s1_iter,s2_iter);
		if(!equals(x1,x2))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation incorrectly overloaded for Id without coupling"<<endlog;
		    return false;
		}  
		CHECK_SFF1(Id,gen);
	    }

	    /* Entries of the Dirac conjugation matrix: */

	    static const value_type& M(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return D_matrix[i][j];
	    }

	    /* Entries of the charge conjugation matrix: */

	    static const value_type& C(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return C_matrix[i][j];
	    }
	    
	    /* Basic conjugate spinor bilinear without prefactor x +=
	     * bar{psi1}C.psi2: */

	    static void C_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			x+=A_1[i]*C_matrix[i][j]*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}C.psi2: */

	    static void C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(C_matrix,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*phi*C.psi2:
	     * */

	    static void C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(C_matrix,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi2*bar{psi1}.C: */

	    static void C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(C_matrix,C_0,A_0,A_1,A_2)
	    }

	    /* Left multiplication by C , psi1 += C.psi2: */

	    static void C_left(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(C_matrix[i][j]*A_1[j]);
		    }
		}
	    }

	    /* Right multiplication by C, psi1 += psi2.C: */

	    static void C_right(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(A_1[j]*C_matrix[j][i]);
		    }
		}
	    }

	    /* Function checking whether the C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_C(generator& gen)
	    {
		fill_scal(gen);
		value_type x1(0,0);
		value_type x2(0,0);
		C_first(x1,s1_iter,s2_iter);
		type::C_first(x2,s1_iter,s2_iter);
		if(!equals(x1,x2))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation incorrectly overloaded for C without coupling"<<endlog;
		    return false;
		}
		fill_scal(gen);
		s2.reset();
		check=s2;
		C_left(s2.begin(),s1.begin());
		type::C_left(check.begin(),s1.begin());
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"left charge conjugation by C incorrectly overloaded"<<endlog;
		    return false;
		}
		fill_scal(gen);
		s1.reset();
		check=s1;
		C_right(s1.begin(),s2.begin());
		type::C_right(check.begin(),s2.begin());
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"right charge conjugation by C incorrectly overloaded"<<endlog;
		    return false;
		}
		CHECK_SFF1(C,gen);
	    }

	    /* Entries of the complex conjugate charge conjugation matrix: */

	    static const value_type& Cc(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return std::conj(C_matrix)[i][j];
	    }

	    /* Basic conjugate spinor bilinear without prefactor x +=
	     * bar{psi1}.C^*.psi2: */
	    
	    static void Cc_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			x+=A_1[i]*std::conj(C_matrix[i][j])*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.psi2: */

	    static void Cc_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			(*A_0)+=C_0*std::conj(C_matrix[i][j])*A_1[i]*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.psi2: */

	    static void Cc_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_1[i]+=c*(std::conj(C_matrix[i][j]))*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*: */

	    static void Cc_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_2[i]+=c*(std::conj(C_matrix[j][i]))*A_1[j];
		    }
		}
	    }

	    /* Left multiplication by C^*, psi1 += C^*.psi2: */

	    static void Cc_left(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(std::conj(C_matrix[i][j])*A_1[j]);
		    }
		}
	    }

	    /* Right multiplication by C^*, psi1 += psi2.C^*: */

	    static void Cc_right(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(A_1[j]*std::conj(C_matrix[j][i]));
		    }
		}
	    }

	    /* Function checking whether the C*-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc(generator& gen)
	    {
		fill_scal(gen);
		value_type x1(0,0);
		value_type x2(0,0);
		Cc_first(x1,s1_iter,s2_iter);
		type::Cc_first(x2,s1_iter,s2_iter);
		if(!equals(x1,x2))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation incorrectly overloaded for C* without coupling"<<endlog;
		    return false;
		}  
		fill_scal(gen);
		s2.reset();
		check=s2;
		Cc_left(s2.begin(),s1.begin());
		type::Cc_left(check.begin(),s1.begin());
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"left charge conjugation by C* incorrectly overloaded"<<endlog;
		    return false;
		}
		fill_scal(gen);
		s1.reset();
		check=s1;
		Cc_right(s1.begin(),s2.begin());
		type::Cc_right(check.begin(),s2.begin());
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"right charge conjugation by C* incorrectly overloaded"<<endlog;
		    return false;
		}
		CHECK_SFF1(Cc,gen);
	    }

	    /* Entries of the gamma matrices: */

	    static const value_type& gamma(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return g[mu][i][j];
	    }

	    /* Spinor vertex recursive relation V(mu) += c*bar{psi1}.g(mu).psi2:
	     * */

	    static void g_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.psi2: */

	    static void g_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.psi2: */

	    static void g_second(const value_type& C_0,const momentum_type& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*bar{psi1}.Vslash:
	     * */
	    
	    static void g_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*bar{psi1}.Vslash:
	     * */
	    
	    static void g_third(const value_type& C_0,const momentum_type& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(g,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_g(generator& gen)
	    {
		CHECK_VFF1(g,gen);
	    }

	    /* Entries of g.C : */

	    static const value_type& gamma_C(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return g_C[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).C.psi2: */

	    static void g_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(g_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*Vslash.C.psi2: */

	    static void g_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(g_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.C: */

	    static void g_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(g_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_g_C(generator& gen)
	    {
		CHECK_VFF1(g_C,gen);
	    }

	    /* Entries of C*.g: */

	    static const value_type& Cc_gamma(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return Cc_g[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).psi2 */

	    static void Cc_g_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(Cc_g,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*C^*.Vslash.psi2: */

	    static void Cc_g_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(Cc_g,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*.Vslash: */

	    static void Cc_g_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(Cc_g,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.g-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_g(generator& gen)
	    {
		CHECK_VFF1(Cc_g,gen);
	    }

	    /* Entries of the chiral Dirac matrix: */

	    static const value_type& gamma_5(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return g_5[i][j];
	    }

	    /* Spinor vertex recursive relation phi += c*bar{psi1}.g5.psi2: */

	    static void g5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(g_5,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*g5.psi2: */

	    static void g5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(g_5,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*phi*bar{psi1}.g5:
	     * */

	    static void g5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(g_5,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g5-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_g5(generator& gen)
	    {
		CHECK_SFF1(g5,gen);
	    }

	    /* Entries of g5.C: */

	    static const value_type& gamma_5_C(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return g_5_C[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.g5.C.psi2: */

	    static void g5_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(g_5_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*g5.C.psi2: */

	    static void g5_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(g_5_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.g5.C: */

	    static void g5_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(g_5_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g5.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_g5_C(generator& gen)
	    {
		CHECK_SFF1(g5_C,gen);
	    }

	    /* Entries of C*.g5: */

	    static const value_type& Cc_gamma_5(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return Cc_g_5[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.g5.psi2: */

	    static void Cc_g5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(Cc_g_5,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.g5.psi2: */

	    static void Cc_g5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(Cc_g_5,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.g5: */

	    static void Cc_g5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(Cc_g_5,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.g5-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_g5(generator& gen)
	    {
		CHECK_SFF1(Cc_g5,gen);
	    }

	    /* Entries of (c0+c1*g5): */

	    static const value_type& gamma_VA(const value_type& C_0,const value_type& C_1,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		if(i==j)
		{
		    return C_0+C_1*g_5[i][j];
		}
		return C_1*g_5[i][j];
	    }

	    /* Spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(c0+c1*g5).psi2: */

	    static void gVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    (*A_0)+=C_0*A_1[i]*A_2[i];
		    for(size_type j=0;j<index_range;++j)
		    {
			(*A_0)+=(C_1*g_5[i][j])*A_1[i]*A_2[j];
		    }
		}
	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*(c0+c1*g5).psi2:
	     * */

	    static void gVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type z0=C_0*(*A_0);
		value_type z1=C_1*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    A_1[i]+=(z0*A_2[i]);
		    for(size_type j=0;j<index_range;++j)
		    {
			A_1[i]+=(z1*g_5[i][j]*A_2[j]);
		    }
		}
	    }

	    /* spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(c0+c1*g5): */

	    static void gVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type z0=C_0*(*A_0);
		value_type z1=C_1*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    A_2[i]+=(z0*A_1[i]);
		    for(size_type j=0;j<index_range;++j)
		    {
			A_2[i]+=(z1*g_5[j][i]*A_1[j]);
		    }
		}
	    }

	    /* Function checking whether the gVA-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gVA(generator& gen)
	    {
		CHECK_SFF2(gVA,gen);
	    }

	    /* Entries of (c0+c1*g5).C: */

	    static const value_type& gamma_VA_C(const value_type& C_0,const value_type& C_1,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return C_0*C_matrix[i][j]+C_1*g_5_C[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(c0+c1*g5).C.psi2: */
	    
	    static void gVA_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST2(C_matrix,g_5_C,C_0,C_1,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*(c0+c1*g5).C.psi2: */
	    
	    static void gVA_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND2(C_matrix,g_5_C,C_0,C_1,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(c0+c1*g5).C: */
	    
	    static void gVA_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD2(C_matrix,g_5_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the gVA.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gVA_C(generator& gen)
	    {
		CHECK_SFF2(gVA_C,gen);
	    }

	    /* Entries of C*.(c0+c1*g5): */

	    static const value_type& Cc_gamma_VA(const value_type& C_0,const value_type& C_1,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return C_0*std::conj(C_matrix[i][j])+C_1*Cc_g_5[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.(c0+c1*g5).psi2: */
	    
	    static void Cc_gVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			(*A_0)+=(C_0*std::conj(C_matrix[i][j])+C_1*Cc_g_5[i][j])*A_1[i]*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.(c0+c1*g5).psi2: */
	    
	    static void Cc_gVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type z0=C_0*(*A_0);
		value_type z1=C_1*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_1[i]+=(z0*std::conj(C_matrix[i][j])+z1*Cc_g_5[i][j])*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.(c0+c1*g5): */
	    
	    static void Cc_gVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type z0=C_0*(*A_0);
		value_type z1=C_1*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_2[i]+=(z0*std::conj(C_matrix[j][i])+z1*Cc_g_5[j][i])*A_1[j];
		    }
		}
	    }

	    /* Function checking whether the C*.gVA-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_gVA(generator& gen)
	    {
		CHECK_SFF2(Cc_gVA,gen);
	    }

	    /* Entries of the left projector: */

	    static const value_type& gamma_L(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return g_L[i][j];
	    }

	    /* Spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1+g5/2).psi2: */
	    
	    static void gL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(g_L,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*(1+g5/2).psi2: */
	    
	    static void gL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(g_L,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1+g5/2): */

	    static void gL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(g_L,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the gL-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gL(generator& gen)
	    {
		CHECK_SFF1(gL,gen);
	    }

	    /* Entries of PL.C: */

	    static const value_type& gamma_L_C(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return g_L_C[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1+g5/2).C.psi2: */

	    static void gL_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(g_L_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*(1+g5/2).C.psi2: */
	    
	    static void gL_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(g_L_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1+g5/2).C: */
	    
	    static void gL_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(g_L_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the gL.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gL_C(generator& gen)
	    {
		CHECK_SFF1(gL_C,gen);
	    }

	    /* Entries of C*.PL: */

	    static const value_type& Cc_gamma_L(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return Cc_g_L[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.(1+g5/2).psi2: */
	    
	    static void Cc_gL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(Cc_g_L,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.(1+g5/2).psi2: */
	    
	    static void Cc_gL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(Cc_g_L,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.(1+g5/2): */
	    
	    static void Cc_gL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(Cc_g_L,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.gL-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_gL(generator& gen)
	    {
		CHECK_SFF1(Cc_gL,gen);
	    }

	    /* Entries of the right projector: */

	    static const value_type& gamma_R(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return g_R[i][j];
	    }

	    /* Spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1-g5/2).psi2: */
	    
	    static void gR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(g_R,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*(1-g5/2).psi2: */
	    
	    static void gR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(g_R,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1-g5/2): */
	    
	    static void gR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(g_R,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the gR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gR(generator& gen)
	    {
		CHECK_SFF1(gR,gen);
	    }

	    /* Entries of PR.C: */

	    static const value_type& gamma_R_C(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return g_R_C[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1-g5/2).C.psi2: */
	    
	    static void gR_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(g_R_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*(1-g5/2).C.psi2: */
	    
	    static void gR_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(g_R_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1-g5/2).C: */
	    
	    static void gR_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		S_THIRD(g_R_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the gR.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gR_C(generator& gen)
	    {
		CHECK_SFF1(gR_C,gen);
	    }

	    /* Entries of C*.PR: */

	    static const value_type& Cc_gamma_R(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return Cc_g_R[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.(1-g5/2).psi2: */
	    
	    static void Cc_gR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(Cc_g_R,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.(1-g5/2).psi2: */
	    
	    static void Cc_gR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(Cc_g_R,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.(1-g5/2): */
	    
	    static void Cc_gR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD(Cc_g_R,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.gR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_gR(generator& gen)
	    {
		CHECK_SFF1(gR,gen);
	    }

	    /* Entries of the left-right combination: */

	    static const value_type& gamma_LR(const value_type& C_0,const value_type& C_1,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return C_0*g_L[i][j]+C_1*g_R[i][j];
	    }

	    /* Spinor vertex recursive relation phi +=
	     * bar{psi1}.(c0*PL+c1*PR).psi2: */
	    
	    static void gLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST2(g_L,g_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += phi*(c0*PL+c1*PR).psi2:
	     * */

	    static void gLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND2(g_L,g_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * phi*bar{psi1}.(c0*PL+c1*PR): */

	    static void gLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD2(g_L,g_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the gLR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gLR(generator& gen)
	    {
		CHECK_SFF2(gLR,gen);
	    }

	    /* Entries of (c0.PL+c1.PR).C: */

	    static const value_type& gamma_LR_C(const value_type& C_0,const value_type& C_1,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return C_0*g_L_C[i][j]+C_1*g_R_C[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * bar{psi1}.(c0*PL+c1*PR).C.psi2: */
	    
	    static void gLR_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST2(g_L_C,g_R_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * phi*(c0*PL+c1*PR).C.psi2: */

	    static void gLR_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND2(g_L_C,g_R_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * phi*bar{psi1}.(c0*PL+c1*PR).C: */

	    static void gLR_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD2(g_L_C,g_R_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the gLR.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gLR_C(generator& gen)
	    {
		CHECK_SFF2(gLR_C,gen);
	    }

	    /* Entries of C*.(c0*PL+c1* PR): */

	    static const value_type& Cc_gamma_LR(const value_type& C_0,const value_type& C_1,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return C_0*Cc_g_L[i][j]+C_1*Cc_g_R[i][j];
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * bar{psi1}.C^*.(c0*PL+c1*PR).psi2: */
	    
	    static void Cc_gLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST2(Cc_g_L,Cc_g_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * phi*C^*.(c0*PL+c1*PR).C.psi2: */

	    static void Cc_gLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND2(Cc_g_L,Cc_g_R,C_0,C_1,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * phi*bar{psi1}.C^*.(c0*PL+c1*PR): */
	    
	    static void Cc_gLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_THIRD2(Cc_g_L,Cc_g_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.gLR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_gLR(generator& gen)
	    {
		CHECK_SFF2(Cc_gLR,gen);
	    }

	    /* Entries of g.g5: */

	    static const value_type& ggamma_5(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return gg_5[mu][i][j];
	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).g5.psi2: */

	    static void gg5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(gg_5,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.g5.psi2: */

	    static void gg5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(gg_5,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.g5.psi2: */

	    static void gg5_second(const value_type& C_0,const momentum_type& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(gg_5,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.g5: */

	    static void gg5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(gg_5,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.g5-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gg5(generator& gen)
	    {
		CHECK_VFF1(gg5,gen);
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.g5: */

	    static void gg5_third(const value_type& C_0,const momentum_type& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(gg_5,C_0,A_0,A_1,A_2)
	    }

	    /* Entries of g.g5.C: */

	    static const value_type& ggamma_5_C(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return gg_5_C[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).g5.C.psi2: */
	    
	    static void gg5_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(gg_5_C,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.g5.C.psi2: */
	    
	    static void gg5_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(gg_5_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.g5.C: */

	    static void gg5_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(gg_5_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.g5.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gg5_C(generator& gen)
	    {
		CHECK_VFF1(gg5_C,gen);
	    }

	    /* Entries of C*.g5.g: */

	    static const value_type& Cc_ggamma_5(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return Cc_gg_5[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).g5.psi2: */

	    static void Cc_gg5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(Cc_gg_5,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*C^*.Vslash.g5.psi2: */
	    
	    static void Cc_gg5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(Cc_gg_5,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*.Vslash.g5: */
	    
	    static void Cc_gg5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(Cc_gg_5,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.g.g5-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_gg5(generator& gen)
	    {
		CHECK_VFF1(Cc_gg5,gen);
	    }

	    /* Entries of g.(c0+c1*g5): */

	    static const value_type& ggamma_VA(const value_type& C_0,const value_type& C_1,size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return C_0*g[mu][i][j]+C_1*gg_5[mu][i][j];
	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0+c1*g5).psi2: */

	    static void ggVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST2(g,gg_5,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += Vslash.(c0+c1*g5).psi2:
	     * */

	    static void ggVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND2(g,gg_5,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0+c1*g5): */

	    static void ggVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD2(g,gg_5,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.gVA-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggVA(generator& gen)
	    {
		CHECK_VFF2(ggVA,gen);
	    }

	    /* Entries of g.(c0+c1*g5).C: */

	    static const value_type& ggamma_VA_C(const value_type& C_0,const value_type& C_1,size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return C_0*g_C[mu][i][j]+C_1*gg_5_C[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0+c1*g5).C.psi2: */
	    
	    static void ggVA_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST2(g_C,gg_5_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * Vslash.(c0+c1*g5).C.psi2: */
	    
	    static void ggVA_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND2(g_C,gg_5_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0+c1*g5).C: */

	    static void ggVA_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD2(g_C,gg_5_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.gVA.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggVA_C(generator& gen)
	    {
		CHECK_VFF2(ggVA_C,gen);
	    }

	    /* Entries of C*.(c0+c1*g5): */

	    static const value_type& Cc_ggamma_VA(const value_type& C_0,const value_type& C_1,size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return C_0*Cc_g[mu][i][j]+C_1*Cc_gg_5[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.C^*.g(mu).(c0+c1*g5).psi2: */
	    
	    static void Cc_ggVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST2(Cc_g,Cc_gg_5,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * C^*.Vslash.(c0+c1*g5).psi2: */

	    static void Cc_ggVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND2(Cc_g,Cc_gg_5,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.C^*.Vslash.(c0+c1*g5): */

	    static void Cc_ggVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD2(Cc_g,Cc_gg_5,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.g.gVA-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_ggVA(generator& gen)
	    {
		CHECK_VFF2(Cc_ggVA,gen);
	    }

	    /* Entries of g.PL: */

	    static const value_type& ggamma_L(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return gg_L[mu][i][j];
	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1+g5/2).psi2: */

	    static void ggL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(gg_L,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1+g5/2).psi2: */

	    static void ggL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(gg_L,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1+g5/2): */

	    static void ggL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(gg_L,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.gL-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggL(generator& gen)
	    {
		CHECK_VFF1(ggL,gen);
	    }

	    /* Entries of g.PL.C: */

	    static const value_type& ggamma_L_C(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return gg_L_C[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1+g5/2).C.psi2: */

	    static void ggL_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(gg_L_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1+g5/2).C.psi2: */

	    static void ggL_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(gg_L_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1+g5/2).C: */

	    static void ggL_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(gg_L_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.gL.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggL_C(generator& gen)
	    {
		CHECK_VFF1(ggL_C,gen);
	    }

	    /* Entries of C*.g.PL: */

	    static const value_type& Cc_ggamma_L(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return Cc_gg_L[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).(1+g5/2).psi2: */
	    
	    static void Cc_ggL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(Cc_gg_L,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.C^*.(1+g5/2).psi2: */

	    static void Cc_ggL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(Cc_gg_L,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.C^*.(1+g5/2): */

	    static void Cc_ggL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(Cc_gg_L,C_0,A_0,A_1,A_2)
	    }
	    /* Function checking whether the C*.g.gL-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_ggL(generator& gen)
	    {
		CHECK_VFF1(Cc_ggL,gen);
	    }

	    /* Entries of g.PR: */

	    static const value_type& ggamma_R(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return gg_R[mu][i][j];
	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1-g5/2).psi2: */
	    
	    static void ggR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {	
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(gg_R,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1-g5/2).psi2: */
	    
	    static void ggR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(gg_R,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1-g5/2): */

	    static void ggR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(gg_R,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.gR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggR(generator& gen)
	    {
		CHECK_VFF1(ggR,gen);
	    }

	    /* Entries of PR.C: */

	    static const value_type& ggamma_R_C(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return gg_R_C[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1-g5/2).C.psi2: */

	    static void ggR_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(gg_R_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1-g5/2).C.psi2: */

	    static void ggR_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(gg_R_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1-g5/2).C: */

	    static void ggR_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(gg_R_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.gR.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggR_C(generator& gen)
	    {
		CHECK_VFF1(ggR_C,gen);
	    }

	    /* Entries of C*.g.PR: */

	    static const value_type& Cc_ggamma_R(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return Cc_gg_R[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).(1-g5/2).psi2: */
	    
	    static void Cc_ggR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST(Cc_gg_R,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.C^*.(1-g5/2).psi2: */

	    static void Cc_ggR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND(Cc_gg_R,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.C^*.(1-g5/2): */

	    static void Cc_ggR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD(Cc_gg_R,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Function checking whether the C*.g.gR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_ggR(generator& gen)
	    {
		CHECK_VFF1(Cc_ggR,gen);
	    }

	    /* Entries of g.(c0*PL+c1*PR): */

	    static const value_type& ggamma_LR(const value_type& C_0,const value_type& C_1,size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return C_0*gg_L[mu][i][j]+C_1*gg_R[mu][i][j];
	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0*PL+c1*PR).psi2: */
	    
	    static void ggLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST2(gg_L,gg_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 +=
	     * Vslash.(c0*PL+c1*PR).psi2: */

	    static void ggLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND2(gg_L,gg_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0*PL+c1*PR): */

	    static void ggLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD2(gg_L,gg_R,C_0,C_1,A_0,A_1,A_2)
	    }
	    
	    /* Function checking whether the g.gLR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggLR(generator& gen)
	    {
		CHECK_VFF2(ggLR,gen);
	    }


	    /* Entries of g.(c0*PL+c1*PR).C: */

	    static const value_type& ggamma_LR_C(const value_type& C_0,const value_type& C_1,size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return C_0*gg_L_C[mu][i][j]+C_1*gg_R_C[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0*PL+c1*PR).C.psi2: */
	    
	    static void ggLR_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST2(gg_L_C,gg_R_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * Vslash.(c0*PL+c1*PR).C.psi2: */

	    static void ggLR_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND2(gg_L_C,gg_R_C,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0*PL+c1*PR).C: */

	    static void ggLR_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD2(gg_L_C,gg_R_C,C_0,C_1,A_0,A_1,A_2)
	    }
	    
	    /* Function checking whether the g.gLR.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_ggLR_C(generator& gen)
	    {
		CHECK_VFF2(ggLR_C,gen);
	    }

	    /* Entries of C*.g.(c0*PL+c1*PR): */

	    static const value_type& Cc_ggamma_LR(const value_type& C_0,const value_type& C_1,size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");
		
		return C_0*Cc_gg_L[mu][i][j]+C_1*Cc_gg_R[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.C^*.g(mu).(c0*PL+c1*PR).psi2: */
	    
	    static void Cc_ggLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_FIRST2(Cc_gg_L,Cc_gg_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * C^*.Vslash.(c0*PL+c1*PR).psi2: */

	    static void Cc_ggLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_SECOND2(Cc_gg_L,Cc_gg_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.C^*.Vslash.(c0*PL+c1*PR): */

	    static void Cc_ggLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		V_THIRD2(Cc_gg_L,Cc_gg_R,C_0,C_1,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.g.gLR-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_ggLR(generator& gen)
	    {
		CHECK_VFF2(Cc_ggLR,gen);
	    }

	    /* Entries of [g,g]: */

	    static const value_type& gamma_commutator(size_type mu,size_type nu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((nu>=spacetime_dimension),"second index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"third index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"fourth index out of range");
		
		return g_comms[mu][nu][i][j];
	    }

	    /* Spinor vertex recursive relation value_t(mu,nu)
	     * +=c*bar{psi1}.[g(mu),g(nu)].psi2: */

	    static void gg_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_FIRST(g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1+=slash{value_t}.psi2: */

	    static void gg_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_SECOND(g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2}+=bar{psi1}.slash{value_t}: */

	    static void gg_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_THIRD(g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the [g,g]-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gg(generator& gen)
	    {
		CHECK_TFF1(gg,gen);
	    }

	    /* Entries of [g,g].C: */

	    static const value_type& gamma_commutator_C(size_type mu,size_type nu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((nu>=spacetime_dimension),"second index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"third index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"fourth index out of range");
		
		return g_comms_C[mu][nu][i][j];
	    }

	    /* Spinor vertex recursive relation value_t(mu,nu)
	     * +=c*bar{psi1}.[g(mu),g(nu)].C.psi2: */

	    static void gg_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_FIRST(g_comms_C,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1+=slash{value_t}.C.psi2: */

	    static void gg_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_SECOND(g_comms_C,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2}+=bar{psi1}.slash{value_t}.C: */

	    static void gg_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_THIRD(g_comms_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the [g,g].C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gg_C(generator& gen)
	    {
		CHECK_TFF1(gg_C,gen);
	    }

	    /* Entries of C*.[g,g]: */

	    static const value_type& Cc_gamma_commutator(size_type mu,size_type nu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((nu>=spacetime_dimension),"second index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"third index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"fourth index out of range");
		
		return Cc_g_comms[mu][nu][i][j];
	    }

	    /* Spinor vertex recursive relation value_t(mu,nu)
	     * +=c*bar{psi1}.C*.[g(mu),g(nu)].psi2: */

	    static void Cc_gg_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_FIRST(Cc_g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1+=C*.slash{value_t}.psi2: */

	    static void Cc_gg_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_SECOND(Cc_g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2}+=bar{psi1}.C*.slash{value_t}: */

	    static void Cc_gg_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		T_THIRD(Cc_g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the Cc_[g,g]-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_gg(generator& gen)
	    {
		CHECK_TFF1(Cc_gg,gen);
	    }

	    /* Massless Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion subamplitude to be propagated and p the momentum. */

	    static void massless_prop(const value_type& C_0,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		if(DIAGONAL)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[i][j]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }
	    
	    /* Correct massless propagator overloading checking function: */
	    
	    template<class generator>static bool check_massless_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s1;
		type::massless_prop(c1,s1_iter,momentum);
		massless_prop(c1,check_iter,momentum);
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massless propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Massless Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion-multiplet subamplitude to be propagated and p the momentum. */

	    static void massless_prop(const value_type& C_0,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		if(DIAGONAL)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[i][j]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }

	    /* Massive Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion subamplitude to be propagated, p the momentum and
	     * m the mass. */

	    static void massive_prop(const value_type& C_0,const value_type& m,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		
		if(DIAGONAL)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[i][j]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }
	    
	    /* Correct massive propagator overloading checking function: */
	    
	    template<class generator>static bool check_massive_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s1;
		type::massive_prop(c1,value_type(mass,0),s1_iter,momentum);
		massive_prop(c1,value_type(mass,0),check_iter,momentum);
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massive propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Massive Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion multiplet subamplitude to be propagated, p the
	     * momentum and m the mass. */

	    static void massive_prop(const value_type& C_0,const value_type& m,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[i][j]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }

	    /* Massless Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion subamplitude to be propagated
	     * and p the momentum. */

	    static void massless_anti_prop(const value_type& C_0,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		
		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[j][i]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }

	    /* Massless Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated and p the momentum. */

	    static void massless_anti_prop(const value_type& C_0,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[j][i]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }
	    
	    /* Correct massless anti-propagator overloading checking function: */
	    
	    template<class generator>static bool check_massless_anti_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s2;
		type::massless_anti_prop(c1,s2_iter,momentum);
		massless_anti_prop(c1,check_iter,momentum);
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massless anti-propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Massive Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated, p the momentum and m the mass. */

	    static void massive_anti_prop(const value_type& C_0,const value_type& m,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		
		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[j][i]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }
	    
	    /* Correct massive anti-propagator overloading checking function: */

	    template<class generator>static bool check_massive_anti_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s2;
		type::massive_anti_prop(c1,value_type(mass,0),s2_iter,momentum);
		massive_anti_prop(c1,value_type(mass,0),check_iter,momentum);
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massive anti-propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Massive Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated, p the momentum and m the mass. */

	    static void massive_anti_prop(const value_type& C_0,const value_type& m,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[j][i]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }

	    /* Dirac conjugation of spinor: */

	    static void make_bar(iterator it1,iterator it2)
	    {
		CAMGEN_ERROR_IF((it1.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((it2.range()<index_range),"tensor iterator 1 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    it1[i]=value_type(0,0);
		    for(size_type j=0;j<index_range;++j)
		    {
			it1[i]+=std::conj(it2[j])*D_matrix[j][i];
		    }
		}
	    }

	protected:

	    /* Static data contained in the Dirac algebra base class: */

	    /* The Dirac matrices: */

	    static value_type g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrices multiplied with the charge conjugation
	     * matrix: */
	    
	    static value_type g_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrices multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The chiral gamma matrix: */

	    static value_type g_5[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The chiral Dirac matrix multiplied with the charge conjugation
	     * matrix: */
	    
	    static value_type g_5_C[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The chiral Dirac matrix multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_g_5[Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The left projection matrix: */

	    static value_type g_L[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The left projection matrix multiplied with the charge conjugation
	     * matrix: */
	    
	    static value_type g_L_C[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The left projection matrix multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_g_L[Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The right projection matrix: */

	    static value_type g_R[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The right projection matrix multiplied with the charge conjugation
	     * matrix: */
	    
	    static value_type g_R_C[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The right projection multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_g_R[Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The chiral matrix vector: */

	    static value_type gg_5[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The chiral matrix vector multiplied with the charge conjugation
	     * matrix*/
	    
	    static value_type gg_5_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The chiral matrix vector multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_gg_5[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The left matrix vector: */

	    static value_type gg_L[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The left projection matrix vector multiplied with the charge conjugation
	     * matrix*/
	    
	    static value_type gg_L_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The left projection matrix vector multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_gg_L[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The right matrix vector: */

	    static value_type gg_R[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The right projection matrix vector multiplied with the charge conjugation
	     * matrix*/
	    
	    static value_type gg_R_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The right projection matrix vector multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_gg_R[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The Dirac matrix commutators: */
	    
	    static value_type g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrix commutators multiplied with the charge conjugation
	     * matrix*/
	    
	    static value_type g_comms_C[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrix commutators multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The charge conjugation matrix: */

	    static value_type C_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac conjugation matrix: */
	    
	    static value_type D_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* A Dirac matrix to store temporary slashed-vectors: */
	    
	    static value_type Vslash[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	
	private:

	    /* Initialisation tag: */

	    static bool initialised;
	    
	    /* A scalar, vector, tensor, 2 spinors and a real momentum
	     * to check overloaded recursive relations: */
	    
	    static tensor_type scal;
	    static tensor_type vec;
	    static tensor_type tens;
	    static tensor_type s1;
	    static tensor_type s2;
	    static momentum_type momentum;
	    static tensor_type check;

	    /* Corresponding iterators: */

	    static iterator scal_iter;
	    static iterator vec_iter;
	    static iterator tens_iter;
	    static iterator s1_iter;
	    static iterator s2_iter;
	    static iterator check_iter;

	    /* Two complex couplings and a real mass which will be assigned a
	     * random value each check: */
	    
	    static value_type c1;
	    static value_type c2;
	    static r_value_type mass;

	    /* Function randomly filling the tensors for a sff-vertex check: */

	    template<class generator>static void fill_scal(generator& gen)
	    {
		c1=value_type(gen(0,1),gen(0,1));
		c2=value_type(gen(0,1),gen(0,1));
		scal[0]=value_type(gen(0,1),gen(0,1));
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }

	    /* Function randomly filling the tensors for a vff-vertex check: */

	    template<class generator>static void fill_vec(generator& gen)
	    {
		c1=value_type(gen(0,1),gen(0,1));
		c2=value_type(gen(0,1),gen(0,1));
		for(size_type i=0;i<spacetime_dimension;++i)
		{
		    vec[i]=value_type(gen(0,1),gen(0,1));
		}
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }

	    /* Function randomly filling the tensors for a tff-vertex check: */

	    template<class generator>static void fill_tens(generator& gen)
	    {
		c1=value_type(gen(0,1),gen(0,1));
		c2=value_type(gen(0,1),gen(0,1));
		for(size_type i=0;i<spacetime_dimension*spacetime_dimension;++i)
		{
		    tens[i]=value_type(gen(0,1),gen(0,1));
		}
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }

	    /* Function randomly filling the momentum, mass and spinors */

	    template<class generator>static void fill_momentum(generator& gen)
	    {
		mass=gen(0,1);
		momentum.fill(gen,0,1);
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }
    };	

    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::D_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::C_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Vslash[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_5[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_5_C[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_g_5[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_L[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_L_C[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_g_L[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_R[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_R_C[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_g_R[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::gg_5[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::gg_5_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_gg_5[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::gg_L[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::gg_L_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_gg_L[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::gg_R[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::gg_R_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_gg_R[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{{0}}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::g_comms_C[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{{0}}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::Cc_g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{{0}}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::tensor_type Dirac_algebra<value_t,type,dim,true>::scal;
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::tensor_type Dirac_algebra<value_t,type,dim,true>::vec(1,dim);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::tensor_type Dirac_algebra<value_t,type,dim,true>::tens(2,dim,dim);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::tensor_type Dirac_algebra<value_t,type,dim,true>::s1(1,Dirac_dim<dim>::value);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::tensor_type Dirac_algebra<value_t,type,dim,true>::s2(1,Dirac_dim<dim>::value);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::momentum_type Dirac_algebra<value_t,type,dim,true>::momentum;
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::tensor_type Dirac_algebra<value_t,type,dim,true>::check;
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::iterator Dirac_algebra<value_t,type,dim,true>::scal_iter=Dirac_algebra<value_t,type,dim,true>::scal.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::iterator Dirac_algebra<value_t,type,dim,true>::vec_iter=Dirac_algebra<value_t,type,dim,true>::vec.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::iterator Dirac_algebra<value_t,type,dim,true>::tens_iter=Dirac_algebra<value_t,type,dim,true>::tens.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::iterator Dirac_algebra<value_t,type,dim,true>::s1_iter=Dirac_algebra<value_t,type,dim,true>::s1.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::iterator Dirac_algebra<value_t,type,dim,true>::s2_iter=Dirac_algebra<value_t,type,dim,true>::s2.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::iterator Dirac_algebra<value_t,type,dim,true>::check_iter=Dirac_algebra<value_t,type,dim,true>::check.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::c1(0,0);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,true>::value_type Dirac_algebra<value_t,type,dim,true>::c2(0,0);
    template<class value_t,class type,std::size_t dim>value_t Dirac_algebra<value_t,type,dim,true>::mass=(value_t)0;
    template<class value_t,class type,std::size_t dim>bool Dirac_algebra<value_t,type,dim,true>::initialised=false;
    template<class value_t,class type,std::size_t dim>const std::size_t Dirac_algebra<value_t,type,dim,true>::index_range;
    template<class value_t,class type,std::size_t dim>const std::size_t Dirac_algebra<value_t,type,dim,true>::spacetime_dimension;
    
    template<class value_t,class type,std::size_t dim>class Dirac_algebra<value_t,type,dim,false>
    {
	public:

	    /* Some convenient type definitions: */

	    typedef value_t r_value_type;
	    typedef std::complex<value_t> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef vector<value_t,dim> momentum_type;
	    typedef typename tensor_type::size_type size_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    
	    /* Definitions of the matrix size (index_range) and the spacetime
	     * dimension: */
	    
	    static const std::size_t index_range=Dirac_dim<dim>::value;
	    static const std::size_t spacetime_dimension=dim;

	    /* Initialisation function: */

	    static void initialise()
	    {
		if(!initialised)
		{

		    /* Let the derived type's spacetime type fill the metric
		     * tensor:*/
		    
		    type::spacetime_type::initialise();
		    
		    /* Let the derived type fill the gamma matrices, the chiral
		     * matrix and charge conjugation matrix: */
		    
		    type::fill_gamma_matrices();
		    type::fill_C_matrix();
		    type::fill_D_matrix();

		    /* Construction of the vector-current matrices and their
		     * products with the charge conjugation matrix: */

		    for(size_type mu=0;mu<spacetime_dimension;++mu)
		    {
			for(size_type i=0;i<index_range;++i)
			{
			    for(size_type j=0;j<index_range;++j)
			    {
				for(size_type k=0;k<index_range;++k)
				{
				    g_C[mu][i][j]+=(g[mu][i][k]*C_matrix[k][j]);
				    Cc_g[mu][i][j]+=(std::conj(C_matrix[i][k])*g[mu][k][j]);
				}
			    }
			}
		    }

		    /* Construction of the tensor current matrices, i.e. the
		     * commutators of Dirac matrices, and their products with
		     * the charge conjugation matrix: */

		    for(size_type mu=0;mu<spacetime_dimension;++mu)
		    {
			for(size_type nu=0;nu<spacetime_dimension;++nu)
			{
			    for(size_type i=0;i<index_range;++i)
			    {
				for(size_type j=0;j<index_range;++j)
				{
				    for(size_type k=0;k<index_range;++k)
				    {
					g_comms[mu][nu][i][j]+=g[mu][i][k]*g[nu][k][j];
					g_comms[mu][nu][i][j]-=g[nu][i][k]*g[mu][k][j];
				    }
				}
			    }
			    for(size_type i=0;i<index_range;++i)
			    {
				for(size_type j=0;j<index_range;++j)
				{
				    for(size_type k=0;k<index_range;++k)
				    {
					g_comms_C[mu][nu][i][j]+=g_comms[mu][nu][i][k]*C_matrix[k][j];
					Cc_g_comms[mu][nu][i][j]+=std::conj(C_matrix[i][k])*g_comms[mu][nu][k][j];
				    }
				}
			    }
			}
		    }
		    initialised=true;
		}
	    }

	    /* Check whether the gamma matrices fulfill the Clifford algebra: */
	    
	    static bool check_anticommutators()
	    {
		value_type c;
		for(size_type mu=0;mu<spacetime_dimension;++mu)
		{
		    for(size_type nu=0;nu<spacetime_dimension;++nu)
		    {
			value_t gmn=type::spacetime_type::metric(mu,nu);
			for(size_type i=0;i<index_range;++i)
			{
			    for(size_type j=0;j<index_range;++j)
			    {
				c=0;
				for(size_type k=0;k<index_range;++k)
				{
				    c+=g[mu][i][k]*g[nu][k][j];
				    c+=g[nu][i][k]*g[mu][k][j];
				}
				if(i==j)
				{
				    if(!equals(c,value_type(2*gmn,0)))
				    {
					log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac algebra not fulfilled: {g"<<mu<<",g"<<nu<<"}("<<i<<","<<j<<")="<<c<<endlog;
					return false;
				    }
				}
				else
				{
				    if(!equals(c,value_type(0,0)))
				    {
					log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac algebra not fulfilled: {g"<<mu<<",g"<<nu<<"}("<<i<<","<<j<<")="<<c<<endlog;
					return false;
				    }
				}
			    }
			}
		    }
		}
		return true;
	    }

	    /* Default implementation of the Dirac conjugation matrix, equal to
	     * the zeroth gamma matrix: */

	    static void fill_D_matrix()
	    {

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			D_matrix[i][j]=g[0][i][j];
		    }
		}

	    }

	    /* Checking unitarity of the charge conjugation matrix: */

	    static bool check_C_matrix()
	    {
		value_type c=0;
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			c=0;
			for(size_type k=0;k<index_range;++k)
			{
			    c+=(std::conj(C_matrix[k][i])*C_matrix[k][j]);

			}
			if(i==j)
			{
			    if(!equals(c,value_type(1,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"conjugation matrix is not unitary: C^+C("<<i<<","<<j<<")="<<c<<endlog;
			    	return false;
			    }
			}
			else
			{
			    if(!equals(c,value_type(0,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"conjugation matrix is not unitary: C^+C("<<i<<","<<j<<")="<<c<<endlog;
			    	return false;
			    }
			}
		    }
		}
		return true;
	    }

	    static bool check_D_matrix()
	    {
	  	/* Checking Hermiticity of the Dirac conjugation matrix: */

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<=i;++j)
		    {
			if(!equals(D_matrix[j][i],std::conj(D_matrix[i][j])))
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Dirac conjugation matrix is not hermitian: D("<<i<<","<<j<<") = "<<D_matrix[i][j]<<endlog;
			    return false;
			}
		    }
		}

		/* Checking if the gamma matrices are self-adjoint under Dirac
		 * conjugation: */

		value_type c;
		for(size_type mu=0;mu<spacetime_dimension;++mu)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    c=0;
			    for(size_type k=0;k<index_range;++k)
			    {
				c+=g[mu][i][k]*D_matrix[k][j];
				c-=D_matrix[i][k]*std::conj(g[mu][j][k]);
			    }
			    if(!equals(c,value_type(0,0)))
			    {
				log(log_level::warning)<<CAMGEN_STREAMLOC<<"Non-self-adjoint gamma matrix g"<<mu<<" detected"<<endlog;
				return false;
			    }
			}
		    }
		}
		return true;
	    }
	    
	    /* Output method: */

	    static std::ostream& print(std::ostream& os)
	    {
		for(size_type mu=0;mu<spacetime_dimension;++mu)
		{
		    os<<"g"<<mu<<"=[ ";
		    for(size_type j=0;j<index_range;++j)
		    {
			os<<std::left;
			os<<std::setw(2);
			shortprint(os,g[mu][0][j]);
			os<<" ";
		    }
		    os<<"]"<<std::endl;
		    for(size_type i=1;i<index_range;++i)
		    {
			os<<std::left<<"   [ ";
			for(size_type j=0;j<index_range;++j)
			{
			    os<<std::left;
			    os<<std::setw(2);
			    shortprint(os,g[mu][i][j]);
			    os<<" ";
			}
			os<<"]"<<std::endl;
		    }
		    os<<std::endl;
		}
		os<<"C= [ ";
		for(size_type j=0;j<index_range;++j)
		{
		    os<<std::left;
		    os<<std::setw(2);
		    shortprint(os,C_matrix[0][j]);
		    os<<" ";
		}
		os<<"]"<<std::endl;
		for(size_type i=1;i<index_range;++i)
		{
		    os<<std::left<<"   [ ";
		    for(size_type j=0;j<index_range;++j)
		    {
			os<<std::left;
			os<<std::setw(2);
			shortprint(os,C_matrix[i][j]);
			os<<" ";
		    }
		    os<<"]"<<std::endl;
		}
		return os;
	    }

	    /* Entries of the unit matrix: */

	    static value_type Id(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		if(i!=j)
		{
		    return value_type(1,0);
		}
		return value_type(0,0);
	    }
	    
	    /* Basic spinor bilinear without prefactor x += bar{psi1}.psi2: */
	    
	    static void Id_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    x+=A_1[i]*A_2[i];
		}
	    }
	    
	    /* Basic spinor bilinear without prefactor x += bar{psi1}.psi2: */
	    
	    static void Id_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    (*A_0)+=C_0*A_1[i]*A_2[i];
		}
	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*psi2 */

	    static void Id_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    A_1[i]+=c*A_2[i];
		}
	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*phi*bar{psi1}: */

	    static void Id_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    A_2[i]+=c*A_1[i];
		}
	    }

	    /* Function checking whether the Id-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Id(generator& gen)
	    {
		fill_scal(gen);
		value_type x1(0,0);
		value_type x2(0,0);
		Id_first(x1,s1_iter,s2_iter);
		type::Id_first(x2,s1_iter,s2_iter);
		if(x1 != x2)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation incorrectly overloaded for Id without coupling"<<endlog;
		    return false;
		}  
		CHECK_SFF1(Id,gen);
	    }

	    /* Entries of the Dirac conjugation matrix: */

	    static const value_type& M(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return D_matrix[i][j];
	    }

	    /* Entries of the charge conjugation matrix: */

	    static const value_type& C(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");
		
		return C_matrix[i][j];
	    }
	    
	    /* Basic conjugate spinor bilinear without prefactor x +=
	     * bar{psi1}C.psi2: */

	    static void C_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			x+=A_1[i]*C_matrix[i][j]*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}C.psi2: */

	    static void C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_FIRST(C_matrix,C_0,A_0,A_1,A_2)
	    }
	    static void C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)

	    /* Conjugate spinor vertex recursive relation psi1 += c*phi*C.psi2:
	     * */

	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");
		
		S_SECOND(C_matrix,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi2*bar{psi1}.C: */

	    static void C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		S_THIRD(C_matrix,C_0,A_0,A_1,A_2)
	    }

	    /* Left multiplication by C , psi1 += C.psi2: */

	    static void C_left(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(C_matrix[i][j]*A_1[j]);
		    }
		}
	    }

	    /* Right multiplication by C, psi1 += psi2.C: */

	    static void C_right(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(A_1[j]*C_matrix[j][i]);
		    }
		}
	    }

	    /* Function checking whether the C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_C(generator& gen)
	    {
		fill_scal(gen);
		value_type x1(0,0);
		value_type x2(0,0);
		C_first(x1,s1_iter,s2_iter);
		type::C_first(x2,s1_iter,s2_iter);
		if(x1 != x2)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation incorrectly overloaded for C without coupling"<<endlog;
		    return false;
		}
		fill_scal(gen);
		s2.reset();
		check=s2;
		C_left(s2.begin(),s1.begin());
		type::C_left(check.begin(),s1.begin());
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"left charge conjugation by C incorrectly overloaded"<<endlog;
		    return false;
		}
		fill_scal(gen);
		s1.reset();
		check=s1;
		C_right(s1.begin(),s2.begin());
		type::C_right(check.begin(),s2.begin());
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"right charge conjugation by C incorrectly overloaded"<<endlog;
		    return false;
		}
		CHECK_SFF1(C,gen);
	    }

	    /* Entries of the complex conjugate charge conjugation matrix: */

	    static const value_type& Cc(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=index_range),"first index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"second index out of range");

		return std::conj(C_matrix)[i][j];
	    }

	    /* Basic conjugate spinor bilinear without prefactor x +=
	     * bar{psi1}.C^*.psi2: */
	    
	    static void Cc_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			x+=A_1[i]*std::conj(C_matrix[i][j])*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.psi2: */

	    static void Cc_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			(*A_0)+=C_0*std::conj(C_matrix[i][j])*A_1[i]*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.psi2: */

	    static void Cc_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_1[i]+=c*(std::conj(C_matrix[i][j]))*A_2[j];
		    }
		}
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*: */

	    static void Cc_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		value_type c=C_0*(*A_0);
		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_2[i]+=c*(std::conj(C_matrix[j][i]))*A_1[j];
		    }
		}
	    }

	    /* Left multiplication by C^*, psi1 += C^*.psi2: */

	    static void Cc_left(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(std::conj(C_matrix[i][j])*A_1[j]);
		    }
		}
	    }

	    /* Right multiplication by C^*, psi1 += psi2.C: */

	    static void Cc_right(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		for(size_type i=0;i<index_range;++i)
		{
		    for(size_type j=0;j<index_range;++j)
		    {
			A_0[i]+=(A_1[j]*std::conj(C_matrix[j][i]));
		    }
		}
	    }

	    /* Function checking whether the C*-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc(generator& gen)
	    {
		fill_scal(gen);
		value_type x1(0,0);
		value_type x2(0,0);
		Cc_first(x1,s1_iter,s2_iter);
		type::Cc_first(x2,s1_iter,s2_iter);
		if(x1 != x2)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first recursive relation incorrectly overloaded for C* without coupling"<<endlog;
		    return false;
		}  
		fill_scal(gen);
		s2.reset();
		check=s2;
		Cc_left(s2.begin(),s1.begin());
		type::Cc_left(check.begin(),s1.begin());
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"left charge conjugation by C* incorrectly overloaded"<<endlog;
		    return false;
		}
		fill_scal(gen);
		s1.reset();
		check=s1;
		Cc_right(s1.begin(),s2.begin());
		type::Cc_right(check.begin(),s2.begin());
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"right charge conjugation by C* incorrectly overloaded"<<endlog;
		    return false;
		}
		CHECK_SFF1(Cc,gen);
	    }

	    /* Entries of the gamma matrices: */

	    static const value_type& gamma(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");

		return g[mu][i][j];
	    }

	    /* Spinor vertex recursive relation V(mu) += c*bar{psi1}.g(mu).psi2:
	     * */

	    static void g_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_FIRST(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.psi2: */

	    static void g_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_SECOND(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.psi2: */

	    static void g_second(const value_type& C_0,const momentum_type& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_SECOND(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*bar{psi1}.Vslash:
	     * */
	    
	    static void g_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_THIRD(g,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*bar{psi1}.Vslash:
	     * */
	    
	    static void g_third(const value_type& C_0,const momentum_type& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_THIRD(g,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_g(generator& gen)
	    {
		CHECK_VFF1(g,gen);
	    }

	    /* Entries of g.C : */

	    static const value_type& gamma_C(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");

		return g_C[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).C.psi2: */

	    static void g_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_FIRST(g_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*Vslash.C.psi2: */

	    static void g_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_SECOND(g_C,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.C: */

	    static void g_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_THIRD(g_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the g.C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_g_C(generator& gen)
	    {
		CHECK_VFF1(g_C,gen);
	    }

	    /* Entries of C*.g: */

	    static const value_type& Cc_gamma(size_type mu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"second index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"third index out of range");

		return Cc_g[mu][i][j];
	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).psi2 */

	    static void Cc_g_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_FIRST(Cc_g,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*C^*.Vslash.psi2: */

	    static void Cc_g_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_SECOND(Cc_g,C_0,A_0,A_1,A_2)
	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*.Vslash: */

	    static void Cc_g_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		V_THIRD(Cc_g,C_0,A_0,A_1,A_2)
	    }
	    
	    /* Function checking whether the C*.g-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_CC_g(generator& gen)
	    {
		CHECK_VFF1(Cc_g,gen);
	    }

	    /* Entries of [g,g]: */

	    static const value_type& gamma_commutator(size_type mu,size_type nu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((nu>=spacetime_dimension),"second index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"third index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"fourth index out of range");

		return g_comms[mu][nu][i][j];
	    }

	    /* Spinor vertex recursive relation value_t(mu,nu)
	     * +=c*bar{psi1}.[g(mu),g(nu)].psi2: */

	    static void gg_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_FIRST(g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1+=slash{value_t}.psi2: */

	    static void gg_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_SECOND(g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2}+=bar{psi1}.slash{value_t}: */

	    static void gg_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_THIRD(g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the [g,g]-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gg(generator& gen)
	    {
		CHECK_TFF1(gg,gen);
	    }

	    /* Entries of [g,g].C: */

	    static const value_type& gamma_commutator_C(size_type mu,size_type nu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((nu>=spacetime_dimension),"second index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"third index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"fourth index out of range");

		return g_comms_C[mu][nu][i][j];
	    }

	    /* Spinor vertex recursive relation value_t(mu,nu)
	     * +=c*bar{psi1}.[g(mu),g(nu)].C.psi2: */

	    static void gg_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_FIRST(g_comms_C,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1+=slash{value_t}.C.psi2: */

	    static void gg_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_SECOND(g_comms_C,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2}+=bar{psi1}.slash{value_t}.C: */

	    static void gg_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_THIRD(g_comms_C,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the [g,g].C-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_gg_C(generator& gen)
	    {
		CHECK_TFF1(gg_C,gen);
	    }

	    /* Entries of C*.[g,g]: */

	    static const value_type& Cc_gamma_commutator(size_type mu,size_type nu,size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((mu>=spacetime_dimension),"first index out of range");
		CAMGEN_ERROR_IF((nu>=spacetime_dimension),"second index out of range");
		CAMGEN_ERROR_IF((i>=index_range),"third index out of range");
		CAMGEN_ERROR_IF((j>=index_range),"fourth index out of range");

		return Cc_g_comms[mu][nu][i][j];
	    }

	    /* Spinor vertex recursive relation value_t(mu,nu)
	     * +=c*bar{psi1}.C*.[g(mu),g(nu)].psi2: */

	    static void Cc_gg_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_FIRST(Cc_g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation psi1+=C*.slash{value_t}.psi2: */

	    static void Cc_gg_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_SECOND(Cc_g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Spinor vertex recursive relation bar{psi2}+=bar{psi1}.C*.slash{value_t}: */

	    static void Cc_gg_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<spacetime_dimension*spacetime_dimension),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<index_range),"tensor iterator 2 out of range");

		T_THIRD(Cc_g_comms,C_0,A_0,A_1,A_2)
	    }

	    /* Function checking whether the C*.[g,g]-functions are correctly
	     * overloaded: */

	    template<class generator>static bool check_Cc_gg(generator& gen)
	    {
		CHECK_TFF1(Cc_gg,gen);
	    }

	    /* Massless Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion subamplitude to be propagated and p the momentum. */

	    static void massless_prop(const value_type& C_0,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");

		if(DIAGONAL)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[i][j]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }

	    /* Massless Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion-multiplet subamplitude to be propagated and p the momentum. */

	    static void massless_prop(const value_type& C_0,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		if(DIAGONAL)
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[i][j]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }
	    
	    /* Correct massless propagator overloading checking function: */
	    
	    template<class generator>static bool check_massless_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s1;
		type::massless_prop(c1,s1_iter,momentum);
		massless_prop(c1,check_iter,momentum);
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massless propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Massive Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion subamplitude to be propagated, p the momentum and
	     * m the mass. */

	    static void massive_prop(const value_type& C_0,const value_type& m,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");

		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[i][j]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }

	    /* Massive Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion multiplet subamplitude to be propagated, p the
	     * momentum and m the mass. */

	    static void massive_prop(const value_type& C_0,const value_type& m,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]+=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[i][j]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }
	    
	    /* Correct massive propagator overloading checking function: */
	    
	    template<class generator>static bool check_massive_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s1;
		type::massive_prop(c1,value_type(mass,0),s1_iter,momentum);
		massive_prop(c1,value_type(mass,0),check_iter,momentum);
		if(!equal_sequences(s1,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massive propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Massless Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion subamplitude to be propagated
	     * and p the momentum. */

	    static void massless_anti_prop(const value_type& C_0,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");

		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[j][i]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }

	    /* Massless Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated and p the momentum. */

	    static void massless_anti_prop(const value_type& C_0,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[j][i]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }
	    
	    /* Correct massless anti-propagator overloading checking function: */
	    
	    template<class generator>static bool check_massless_anti_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s2;
		type::massless_anti_prop(c1,s2_iter,momentum);
		massless_anti_prop(c1,check_iter,momentum);
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massless anti-propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Massive Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated, p the momentum and m the mass. */

	    static void massive_anti_prop(const value_type& C_0,const value_type& m,iterator A_0,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");

		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		for(size_type i=0;i<index_range;++i)
		{
		    temp[i]=0;
		    for(size_type j=0;j<index_range;++j)
		    {
			temp[i]+=Vslash[j][i]*(A_0[j]);
		    }
		}
		for(size_type i=0;i<index_range;++i)
		{
		    A_0[i]=C_0*temp[i];
		}
	    }

	    /* Massive Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated, p the momentum and m the mass. */

	    static void massive_anti_prop(const value_type& C_0,const value_type& m,iterator A_0,iterator A_1,const momentum_type& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<index_range),"tensor iterator 1 out of range");

		if(DIAGONAL)
		{	
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)	
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,mu))*p[mu];
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		else
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			for(size_type j=0;j<index_range;++j)
			{
			    Vslash[i][j]=0;	
			    for(size_type mu=0;mu<spacetime_dimension;++mu)
			    {	
				for(size_type nu=0;nu<spacetime_dimension;++nu)
				{
				    Vslash[i][j]-=g[mu][i][j]*(METRIC(mu,nu))*p[nu];
				}
			    }
			}
			Vslash[i][i]+=m;
		    }
		}
		value_type temp[index_range];
		do
		{
		    for(size_type i=0;i<index_range;++i)
		    {
			temp[i]=0;
			for(size_type j=0;j<index_range;++j)
			{
			    temp[i]+=Vslash[j][i]*(A_0[j]);
			}
		    }
		    for(size_type i=0;i<index_range;++i)
		    {
			(*A_0)=C_0*temp[i];
			++A_0;
		    }
		}
		while(A_0 != A_1);
	    }
	    
	    /* Correct massive anti-propagator overloading checking function: */
	    
	    template<class generator>static bool check_massive_anti_prop(generator& gen)
	    {
		fill_momentum(gen);
		check=s2;
		type::massive_anti_prop(c1,value_type(mass,0),s2_iter,momentum);
		massive_anti_prop(c1,value_type(mass,0),check_iter,momentum);
		if(!equal_sequences(s2,check))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"massive anti-propagator incorrectly overloaded"<<endlog;
		    return false;
		}
		return true;
	    }

	    /* Dirac conjugation of spinor: */

	    static void make_bar(iterator it1,iterator it2)
	    {
		CAMGEN_ERROR_IF((it1.range()<index_range),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((it2.range()<index_range),"tensor iterator 1 out of range");

		for(size_type i=0;i<index_range;++i)
		{
		    it1[i]=value_type(0,0);
		    for(size_type j=0;j<index_range;++j)
		    {
			it1[i]+=std::conj(it2[j])*D_matrix[j][i];
		    }
		}
	    }
	protected:
	    /* Static data contained in the Dirac algebra base class: */

	    /* The Dirac matrices: */

	    static value_type g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrices multiplied with the charge conjugation
	     * matrix: */
	    
	    static value_type g_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrices multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrix commutators: */
	    
	    static value_type g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrix commutators multiplied with the charge conjugation
	     * matrix*/
	    
	    static value_type g_comms_C[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac matrix commutators multiplied from the left with the complex
	     * conjugate of the charge conj. matrix: */
	    
	    static value_type Cc_g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value];

	    /* The charge conjugation matrix: */

	    static value_type C_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* The Dirac conjugation matrix: */
	    
	    static value_type D_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	    
	    /* A Dirac matrix to store temporary slashed-vectors: */
	    
	    static value_type Vslash[Dirac_dim<dim>::value][Dirac_dim<dim>::value];
	
	private:

	    /* Initialisation tag: */

	    static bool initialised;
	    
	    /* A scalar, vector, tensor, 2 spinors and a real momentum
	     * to check overloaded recursive relations: */
	    
	    static tensor_type scal;
	    static tensor_type vec;
	    static tensor_type tens;
	    static tensor_type s1;
	    static tensor_type s2;
	    static momentum_type momentum;
	    static tensor_type check;

	    /* Corresponding iterators: */

	    static iterator scal_iter;
	    static iterator vec_iter;
	    static iterator tens_iter;
	    static iterator s1_iter;
	    static iterator s2_iter;
	    static iterator check_iter;

	    /* Two complex couplings and a real mass which will be assigned a
	     * random value each check: */
	    
	    static value_type c1;
	    static value_type c2;
	    static r_value_type mass;

	    /* Function randomly filling the tensors for a sff-vertex check: */

	    template<class generator>static void fill_scal(generator& gen)
	    {
		c1=value_type(gen(0,1),gen(0,1));
		c2=value_type(gen(0,1),gen(0,1));
		scal[0]=value_type(gen(0,1),gen(0,1));
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }

	    /* Function randomly filling the tensors for a vff-vertex check: */

	    template<class generator>static void fill_vec(generator& gen)
	    {
		c1=value_type(gen(0,1),gen(0,1));
		c2=value_type(gen(0,1),gen(0,1));
		for(size_type i=0;i<spacetime_dimension;++i)
		{
		    vec[i]=value_type(gen(0,1),gen(0,1));
		}
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }

	    /* Function randomly filling the tensors for a tff-vertex check: */

	    template<class generator>static void fill_tens(generator& gen)
	    {
		c1=value_type(gen(0,1),gen(0,1));
		c2=value_type(gen(0,1),gen(0,1));
		for(size_type i=0;i<spacetime_dimension*spacetime_dimension;++i)
		{
		    tens[i]=value_type(gen(0,1),gen(0,1));
		}
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }

	    /* Function randomly filling the momentum, mass and spinors */

	    template<class generator>static void fill_momentum(generator& gen)
	    {
		mass=gen(0,1);
		momentum.fill(gen,0,1);
		for(size_type i=0;i<index_range;++i)
		{
		    s1[i]=value_type(gen(0,1),gen(0,1));
		    s2[i]=value_type(gen(0,1),gen(0,1));
		}
	    }
    };	

    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::g_C[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::Cc_g[dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{0}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::D_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::C_matrix[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::Vslash[Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{0}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{{0}}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::g_comms_C[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{{0}}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::Cc_g_comms[dim][dim][Dirac_dim<dim>::value][Dirac_dim<dim>::value]={{{{0}}}};
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::tensor_type Dirac_algebra<value_t,type,dim,false>::scal;
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::tensor_type Dirac_algebra<value_t,type,dim,false>::vec(1,dim);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::tensor_type Dirac_algebra<value_t,type,dim,false>::tens(2,dim,dim);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::tensor_type Dirac_algebra<value_t,type,dim,false>::s1(1,Dirac_dim<dim>::value);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::tensor_type Dirac_algebra<value_t,type,dim,false>::s2(1,Dirac_dim<dim>::value);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::momentum_type Dirac_algebra<value_t,type,dim,false>::momentum;
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::tensor_type Dirac_algebra<value_t,type,dim,false>::check;
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::iterator Dirac_algebra<value_t,type,dim,false>::scal_iter=Dirac_algebra<value_t,type,dim,false>::scal.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::iterator Dirac_algebra<value_t,type,dim,false>::vec_iter=Dirac_algebra<value_t,type,dim,false>::vec.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::iterator Dirac_algebra<value_t,type,dim,false>::tens_iter=Dirac_algebra<value_t,type,dim,false>::tens.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::iterator Dirac_algebra<value_t,type,dim,false>::s1_iter=Dirac_algebra<value_t,type,dim,false>::s1.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::iterator Dirac_algebra<value_t,type,dim,false>::s2_iter=Dirac_algebra<value_t,type,dim,false>::s2.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::iterator Dirac_algebra<value_t,type,dim,false>::check_iter=Dirac_algebra<value_t,type,dim,false>::check.begin();
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::c1(0,0);
    template<class value_t,class type,std::size_t dim>typename Dirac_algebra<value_t,type,dim,false>::value_type Dirac_algebra<value_t,type,dim,false>::c2(0,0);
    template<class value_t,class type,std::size_t dim>value_t Dirac_algebra<value_t,type,dim,false>::mass=(value_t)0;
    template<class value_t,class type,std::size_t dim>bool Dirac_algebra<value_t,type,dim,false>::initialised=false;
    template<class value_t,class type,std::size_t dim>const std::size_t Dirac_algebra<value_t,type,dim,false>::index_range;
    template<class value_t,class type,std::size_t dim>const std::size_t Dirac_algebra<value_t,type,dim,false>::spacetime_dimension;
}

#endif /*CAMGEN_DIRAC_ALG_H_*/

