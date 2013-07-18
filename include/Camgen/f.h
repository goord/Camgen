//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_F_H_
#define CAMGEN_F_H_

#include <Camgen/comp_vert.h>
#include <Camgen/eval.h>
#include <Camgen/f_helper.h>
#include <Camgen/adj_rep.h>
#include <Camgen/su(n).h>
#include <Camgen/ppranknum.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Structure constant colour structure definition and declaration. Moreover,     *
 * this file contains the specialisations of the evaluate and cfd_evaluate class *
 * templates for vertices which contain a composition with structure constants.  *
 * There is also a definition of the CF_f class, which denotes the structure     *
 * constant contraction in the colour flow treatment.                            *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	/* Structure constants class for generic group indices: */

	template<class group_t,std::size_t I=0,std::size_t J=1,std::size_t K=2>class f
	{
	    public:

		/* Self-referring type definition: */

		typedef f<group_t,I,J,K> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef f<group_t,I,J,K> CF_type;

		/* Number of indices: */

		static const std::size_t rank=3;
		
		/* Total colour tensor size: */
		
		static const std::size_t tensor_size = (group_t::rank)*(group_t::rank)*(group_t::rank);
		
		/* Contracted indices: */
		
		static const std::size_t contractions[3];

		/* Ranges of contracted indices: */

		static const std::size_t ranges[3];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<class group_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<group_t,I,J,K>::rank;
	template<class group_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<group_t,I,J,K>::tensor_size;
	template<class group_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<group_t,I,J,K>::multiplicity;
	template<class group_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<group_t,I,J,K>::contractions[3]={I,J,K};
	template<class group_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<group_t,I,J,K>::ranges[3]={group_t::rank,group_t::rank,group_t::rank};
	template<class group_t,std::size_t I,std::size_t J,std::size_t K>const std::string f<group_t,I,J,K>::formula="f(a1,a2,a3)";

	/* Declaration of the colour-flow structure constant class: */

	template<std::size_t N,std::size_t I=0,std::size_t J=1,std::size_t K=2>class CF_f;

	/* Structure constant specialisation for SU(N): */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>class f<SU<N>,I,J,K>
	{
	    public:

		/* Self-referring type definition: */

		typedef f<SU<N>,I,J,K> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef CF_f<N,I,J,K> CF_type;

		/* Number of indices: */

		static const std::size_t rank=3;
		
		/* Total colour tensor size: */
		
		static const std::size_t tensor_size = (N*N-1)*(N*N-1)*(N*N-1);
		
		/* Contracted indices: */
		
		static const std::size_t contractions[3];

		/* Ranges of contracted indices: */

		static const std::size_t ranges[3]; 

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<SU<N>,I,J,K>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<SU<N>,I,J,K>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<SU<N>,I,J,K>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<SU<N>,I,J,K>::contractions[3]={I,J,K};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t f<SU<N>,I,J,K>::ranges[3]={N*N-1,N*N-1,N*N-1};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::string f<SU<N>,I,J,K>::formula="f(a1,a2,a3)";
	
	/* Definition of the colour-flow structure constant class: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>class CF_f
	{
	    public:

		/* Self-referring type definition: */

		typedef CF_f<N,I,J,K> type;

		/* Number of indices: */

		static const std::size_t rank=3;

		/* Total colour tensor size: */

		static const std::size_t tensor_size=N*N*N*N*N*N;
		
		/* Contracted indices: */
		
		static const std::size_t contractions[3];

		/* Contracted index ranges: */

		static const std::size_t ranges[3];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_f<N,I,J,K>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_f<N,I,J,K>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_f<N,I,J,K>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_f<N,I,J,K>::contractions[3]={I,J,K};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_f<N,I,J,K>::ranges[3]={N*N,N*N,N*N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::string CF_f<N,I,J,K>::formula="2i(d(A1,B2)d(A2,B3)d(A3,B1)-d(A1,B3)d(A2,B1)d(A3,B2))";
    }

    /* Specialisation of the evaluate class template for vertices composed with
     * structure constants: */

    template<class group_t,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::f<group_t,I,J,K>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::f<group_t,I,J,K>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;
	    
	    /* Group implementation type definition: */
	    
	    DEFINE_GROUP_TYPE(model_type,group_t);

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
		colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::initialise();
	    }

	    /* Checking function, adding group index ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K)
		{
		    v.insert(v.begin(),group_t::rank);
		}
		return v;
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(i)															\
for(size_type a=0;a<colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::size();++a)							\
{																	\
iters[I]+=colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::jump(0,a);									\
iters[J]+=colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::jump(1,a);									\
iters[K]+=colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::jump(2,a);									\
evaluate<Feynrule_t>::RANK_NUMERAL(i)(colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::value(a)*factor,couplings,iters,momenta);	\
}																	\
iters[I]-=colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::reset(0);									\
iters[J]-=colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::reset(1);									\
iters[K]-=colour_tensor::f_helper<group_t,Feynrule_t,I,J,K>::reset(2);									\

	    /* Recursive relations: */

	    static void first(ARG_LIST)
	    {
		RECREL(0);
	    }
	    static void second(ARG_LIST)
	    {
		RECREL(1);
	    }
	    static void third(ARG_LIST)
	    {
		RECREL(2);
	    }
	    static void fourth(ARG_LIST)
	    {
		RECREL(3);
	    }
#undef RECREL
    };

#include <Camgen/col_macros.h>

    /* Specialisation of the evaluate class template for the colour-flow version
     * of Feynman rules composed with structure constants: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::CF_f<N,I,J,K>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_f<N,I,J,K>,Feynrule_t> vertex_type;
	    
	    /* The usual type definitions: */

	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Group type definition: */

	    DEFINE_GROUP_TYPE(model_type,SU<N>);

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
	    }

	    /* Checking function, adding group index ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(n)											\
value_type z(0,-1);											\
z*=factor;												\
if(I==n)												\
{													\
    GG_to_G_CHAIN(I,J,K,n,z);										\
    GG_to_G_CHAIN(I,K,J,n,-z);										\
    return;												\
}													\
if(J==n)												\
{													\
    GG_to_G_CHAIN(J,I,K,n,-z);										\
    GG_to_G_CHAIN(J,K,I,n,z);										\
    return;												\
}													\
if(K==n)												\
{													\
    GG_to_G_CHAIN(K,I,J,n,z);										\
    GG_to_G_CHAIN(K,J,I,n,-z);										\
    return;												\
}													\
z*=((r_value_type)2);											\
GGG_to_S_CHAIN(J,I,K,n,-z);										\
GGG_to_S_CHAIN(K,I,J,n,z);										\
return
	    
	    /* Recursive relations: */

	    static void first(ARG_LIST)
	    {
		RECREL(0);
	    }
	    static void second(ARG_LIST)
	    {
		RECREL(1);
	    }
	    static void third(ARG_LIST)
	    {
		RECREL(2);
	    }
	    static void fourth(ARG_LIST)
	    {
		RECREL(3);
	    }

#undef RECREL

	    static const size_type sizes[3][4];
    };

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>const typename evaluate< compose_vertices<colour_tensor::CF_f<N,I,J,K>,Feynrule_t> >::size_type evaluate< compose_vertices<colour_tensor::CF_f<N,I,J,K>,Feynrule_t> >::sizes[3][4]={{Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]},{N*Feynrule_t::sizes[0],N*Feynrule_t::sizes[1],N*Feynrule_t::sizes[2],N*Feynrule_t::sizes[3]},{N*N*Feynrule_t::sizes[0],N*N*Feynrule_t::sizes[1],N*N*Feynrule_t::sizes[2],N*N*Feynrule_t::sizes[3]}};

    /* Specialisation of the cfd_evaluate class template for Feynman rule
     * classes composed with structure constants: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_f<N,I,J,K>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_f<N,I,J,K>,Feynrule_t> vertex_type;

	    /* The usual type definition: */

	    typedef typename cfd_evaluate<Feynrule_t>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_t>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_t>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_t>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Group type definition: */

	    DEFINE_GROUP_TYPE(model_type,SU<N>);

	    /* Initialisation phase: */

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_t>::initialise();
	    }

	    /* Checking function, adding group index ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(n)												\
value_type z=value_type(0,1)*factor;										\
if(I==n)													\
{														\
    size_type b=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;							\
    size_type a=((iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N)-b)/N;					\
    size_type d=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;							\
    size_type c=((iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N)-d)/N;					\
    if(b == c and a != d)											\
    {														\
	iters[I]+=(a*N+d)*Feynrule_t::sizes[I];									\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);			\
	iters[I]-=(a*N+d)*Feynrule_t::sizes[I];									\
	return;													\
    }														\
    if(a == d and b != c)											\
    {														\
	iters[I]+=((c*N+b)*Feynrule_t::sizes[I]);								\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
	iters[I]-=((c*N+b)*Feynrule_t::sizes[I]);								\
	return;													\
    }														\
    if(a == d and b == c and a != b)										\
    {														\
	iterator temp_it=iters[I];										\
	iters[I]+=(a*(N+1)*Feynrule_t::sizes[I]);								\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);			\
	iters[I]=temp_it+(b*(N+1)*Feynrule_t::sizes[I]);							\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
	iters[I]=temp_it;											\
	return;													\
    }														\
    return;													\
}														\
if(J==n)													\
{														\
    size_type b=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;							\
    size_type a=((iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N)-b)/N;					\
    size_type d=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;							\
    size_type c=((iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N)-d)/N;					\
    if(b == c and a != d)											\
    {														\
	iters[J]+=(a*N+d)*Feynrule_t::sizes[J];									\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
	iters[J]-=(a*N+d)*Feynrule_t::sizes[J];									\
	return;													\
    }														\
    if(a == d and b != c)											\
    {														\
	iters[J]+=((c*N+b)*Feynrule_t::sizes[J]);								\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);			\
	iters[J]-=((c*N+b)*Feynrule_t::sizes[J]);								\
	return;													\
    }														\
    if(a == d and b == c and a != b)										\
    {														\
	iterator temp_it=iters[J];										\
	iters[J]+=(a*(N+1)*Feynrule_t::sizes[J]);								\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
	iters[J]=temp_it+(b*(N+1)*Feynrule_t::sizes[J]);							\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);			\
	iters[J]=temp_it;											\
	return;													\
    }														\
    return;													\
}														\
if(K==n)													\
{														\
    size_type b=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;							\
    size_type a=((iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N)-b)/N;					\
    size_type d=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;							\
    size_type c=((iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N)-d)/N;					\
    if(b == c and a != d)											\
    {														\
	iters[K]+=(a*N+d)*Feynrule_t::sizes[K];									\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);			\
	iters[K]-=(a*N+d)*Feynrule_t::sizes[K];									\
	return;													\
    }														\
    if(a == d and b != c)											\
    {														\
	iters[K]+=((c*N+b)*Feynrule_t::sizes[K]);								\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
	iters[K]-=((c*N+b)*Feynrule_t::sizes[K]);								\
	return;													\
    }														\
    if(a == d and b == c and a != b)										\
    {														\
	iterator temp_it=iters[K];										\
	iters[K]+=(a*(N+1)*Feynrule_t::sizes[K]);								\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);			\
	iters[K]=temp_it+(b*(N+1)*Feynrule_t::sizes[K]);							\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
	iters[K]=temp_it;											\
	return;													\
    }														\
    return;													\
}														\
else														\
{    														\
    z*=(r_value_type)2;												\
    size_type b=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;							\
    size_type a=((iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N)-b)/N;					\
    size_type d=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;							\
    size_type c=((iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N)-d)/N;					\
    size_type f=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;							\
    size_type e=((iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N)-f)/N;					\
    if(a==b and b==c and c==d and d==e and e==f)								\
    {														\
	return;													\
    }														\
    if(a==f and b==c and d==e)											\
    {														\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);			\
    }														\
    if(a==d and b==e and c==f)											\
    {														\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
    }														\
}														\

	    static void first(CFD_ARG_LIST)
	    {
		RECREL(0);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		RECREL(1);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		RECREL(2);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		RECREL(3);
	    }

#undef RECREL
    };
}

/* If the Yang-Mills three-vertex is included, a template specialisation of the
 * composition of the vvv-vertex with the structure constants: */

#ifdef CAMGEN_VVV_H_
#include <Camgen/ggg.h>
#endif /*CAMGEN_VVV_H_*/

#include <Camgen/undef_args.h>

#endif /*CAMGEN_F_H_*/

