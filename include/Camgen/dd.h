//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DD_H_
#define CAMGEN_DD_H_

#include <vector>
#include <Camgen/adj_rep.h>
#include <Camgen/su(n).h>
#include <Camgen/comp_vert.h>
#include <Camgen/eval.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Double-Kronecker delta colour structure declaration and definition. Moreover, *
 * this file contains the specialisations of the evaluate and cfd_evaluate class *
 * templates for vertices which contain a composition with a colour-index double *
 * Kronecker delta. There is also a declaration and definition of the CF_dd      *
 * class, which denotes the double-Kronecker delta structure for adjoint         *
 * representation indices in the colour-flow representation.                     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {

	/* Double-Kronecker delta class for generic representation indices: */
	
	template<class rep_t,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class dd
	{
	    public:

		/* Self-referring type definition: */

		typedef dd<rep_t,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef dd<rep_t,I,J,K,L> CF_type;

		/* Vertex rank: */

		static const std::size_t rank=4;
		
		/* Total tensor size: */
		
		static const std::size_t tensor_size = (rep_t::dimension)*(rep_t::dimension)*(rep_t::dimension)*(rep_t::dimension);
		
		/* Contracted indices: */
		
		static const std::size_t contractions[4];
		
		/* Ranges of contracted indices: */
		
		static const std::size_t ranges[4]; 

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<rep_t,I,J,K,L>::rank;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<rep_t,I,J,K,L>::tensor_size;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<rep_t,I,J,K,L>::multiplicity;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<rep_t,I,J,K,L>::contractions[4]={I,J,K,L};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<rep_t,I,J,K,L>::ranges[4]={rep_t::dimension,rep_t::dimension,rep_t::dimension,rep_t::dimension};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string dd<rep_t,I,J,K,L>::formula="d(a1,a2)d(a3,a4)";

	/* Declaration of the colour-flow double-Kronecker delta: */

	template<std::size_t N,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class CF_dd;

	/* Specialisation in the case of the SU(N) adjoint representation: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class dd<adjoint_rep< SU<N> >,I,J,K,L>
	{
	    public:

		/* Self-referring type definition: */

		typedef dd<adjoint_rep< SU<N> >,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef CF_dd<N,I,J,K,L> CF_type;

		/* Vertex rank: */

		static const std::size_t rank=4;
		
		/* Total tensor size: */
		
		static const std::size_t tensor_size=(N*N-1)*(N*N-1)*(N*N-1)*(N*N-1);
		
		/* Contracted indices: */
		
		static const std::size_t contractions[4];
		
		/* Ranges of contracted indices: */
		
		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<adjoint_rep<SU<N> >,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<adjoint_rep<SU<N> >,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<adjoint_rep<SU<N> >,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<adjoint_rep<SU<N> >,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd<adjoint_rep<SU<N> >,I,J,K,L>::ranges[4]={N*N-1,N*N-1,N*N-1,N*N-1};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string dd<adjoint_rep<SU<N> >,I,J,K,L>::formula="d(a1,a2)d(a3,a4)";

	/* Definition of the colour-flow double-Kronecker delta: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class CF_dd
	{
	    public:

		/* Self-referring type definition: */

		typedef CF_dd<N,I,J,K,L> type;

		/* Vertex rank: */

		static const std::size_t rank=4;
		
		/* Total tensor size: */
		
		static const std::size_t tensor_size=N*N*N*N*N*N*N*N;
		
		/* Contracted indices: */
		
		static const std::size_t contractions[4];
		
		/* Ranges of contracted indices: */
		
		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd<N,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd<N,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd<N,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd<N,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd<N,I,J,K,L>::ranges[4]={N*N,N*N,N*N,N*N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string CF_dd<N,I,J,K,L>::formula="4(d(A1,B2)d(A2,B1)d(A3,B4)d(A4,B3)-d(A1,B1)d(A2,B2)d(A3,B4)d(A4,B3)/N-d(A1,B2)d(A2,B1)d(A3,B3)d(A4,B4)/N+d(A1,B1)d(A2,B2)d(A3,B3)d(A4,B4)/N^2)";
    }
    
    /* Specialisation for the evaluate class template for generic Kronecker delta
     * compositions: */

    template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Initialisation phase, to compute the necessary subtensor-sizes: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
	    }

	    /* Checking function, adding representation ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),rep_t::dimension);
		}
		return v;
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(i)									\
for(size_type a=0;a<rep_t::dimension;++a)						\
{											\
    for(size_type b=0;b<rep_t::dimension;++b)						\
    {											\
	evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta);		\
	iters[K]+=sizes[0][K];								\
	iters[L]+=sizes[0][L];								\
    }											\
    iters[I]+=sizes[0][I];								\
    iters[J]+=sizes[0][J];								\
    iters[K]-=sizes[1][K];								\
    iters[L]-=sizes[1][L];								\
}											\
iters[I]-=sizes[1][I];									\
iters[J]-=sizes[1][J];									\

	    /* Recursion relations: */

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
	
	    static const size_type sizes[2][4];
    };

    template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >::sizes[2][4]={{Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]},{rep_t::dimension*Feynrule_t::sizes[0],rep_t::dimension*Feynrule_t::sizes[1],rep_t::dimension*Feynrule_t::sizes[2],rep_t::dimension*Feynrule_t::sizes[3]}};

    /* Evaluate class template specialisation for the colour-flow double
     * Kronecker delta colour structures: */
    
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Initialisation phase, to compute the necessary subtensor-sizes: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
	    }

	    /* Checking function, adding representation ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(i)											\
if(I==i or J==i)											\
{													\
for(size_type a=0;a<N*N;++a)										\
{													\
    for(size_type b=0;b<N;++b)										\
    {													\
	for(size_type c=0;c<N;++c)									\
	{												\
	    evaluate<Feynrule_t>::RANK_NUMERAL(i)((r_value_type)2*factor,couplings,iters,momenta);	\
	    iters[K]+=sizes[0][K];									\
	    iters[L]+=sizes[1][L];									\
	}												\
	iters[L]+=(sizes[0][L]-sizes[2][L]);								\
    }													\
    iters[K]-=sizes[2][K];										\
    iters[L]-=sizes[1][L];										\
    iters[I]+=sizes[0][I];										\
    iters[J]+=sizes[0][J];										\
}													\
iters[I]-=sizes[2][I];											\
iters[J]-=sizes[2][J];											\
return;													\
}													\
for(size_type a=0;a<N*N;++a)										\
{													\
for(size_type b=0;b<N;++b)		    								\
{													\
    for(size_type c=0;c<N;++c)										\
    {													\
	evaluate<Feynrule_t>::RANK_NUMERAL(i)((r_value_type)2*factor,couplings,iters,momenta);		\
	iters[I]+=sizes[0][I];										\
	iters[J]+=sizes[1][J];										\
    }													\
    iters[J]+=(sizes[0][J]-sizes[2][J]);								\
}													\
iters[I]-=sizes[2][I];											\
iters[J]-=sizes[1][J];											\
iters[K]+=sizes[0][K];											\
iters[L]+=sizes[0][L];											\
}													\
iters[K]-=sizes[2][K];											\
iters[L]-=sizes[2][L];	    										\

	    /* Recursion relations: */

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
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::sizes[3][4]={{Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]},{N*Feynrule_t::sizes[0],N*Feynrule_t::sizes[1],N*Feynrule_t::sizes[2],N*Feynrule_t::sizes[3]},{N*N*Feynrule_t::sizes[0],N*N*Feynrule_t::sizes[1],N*N*Feynrule_t::sizes[2],N*N*Feynrule_t::sizes[3]}};

    /* Specialisation of the cfd_evaluate class template for the double Kronecker delta
     * colour structures over SU(N)-fundamental representation indices: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename cfd_evaluate<Feynrule_t>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_t>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_t>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_t>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Initialisation phase, to compute the necessary subtensor-sizes: */

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_t>::initialise();
	    }

	    /* Checking function, adding representation ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro definition: */

#define	RECREL(i)													\
if(I==i)														\
{															\
    if((iters[K].get_offset()/Feynrule_t::sizes[K])%N==(iters[L].get_offset()/Feynrule_t::sizes[L])%N)			\
    {															\
	iters[I]+=((iters[J].get_offset()/Feynrule_t::sizes[J])%N)*Feynrule_t::sizes[I];				\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta,produced_iters);			\
	iters[I]-=((iters[J].get_offset()/Feynrule_t::sizes[J])%N)*Feynrule_t::sizes[I];				\
    }															\
    return;														\
}															\
if(J==i)														\
{															\
    if((iters[K].get_offset()/Feynrule_t::sizes[K])%N==(iters[L].get_offset()/Feynrule_t::sizes[L])%N)			\
    {															\
	iters[J]+=((iters[I].get_offset()/Feynrule_t::sizes[I])%N)*Feynrule_t::sizes[J];				\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta,produced_iters);			\
	iters[J]-=((iters[I].get_offset()/Feynrule_t::sizes[I])%N)*Feynrule_t::sizes[J];				\
    }															\
    return;														\
}															\
if(K==i)														\
{															\
    if((iters[I].get_offset()/Feynrule_t::sizes[I])%N==(iters[J].get_offset()/Feynrule_t::sizes[J])%N)			\
    {															\
	iters[K]+=((iters[L].get_offset()/Feynrule_t::sizes[L])%N)*Feynrule_t::sizes[K];				\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta,produced_iters);			\
	iters[K]-=((iters[L].get_offset()/Feynrule_t::sizes[L])%N)*Feynrule_t::sizes[K];				\
    }															\
    return;														\
}															\
if(L==i)														\
{															\
    if((iters[I].get_offset()/Feynrule_t::sizes[I])%N==(iters[J].get_offset()/Feynrule_t::sizes[J])%N)			\
    {															\
	iters[L]+=((iters[K].get_offset()/Feynrule_t::sizes[K])%N)*Feynrule_t::sizes[L];				\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta,produced_iters);			\
	iters[L]-=((iters[K].get_offset()/Feynrule_t::sizes[K])%N)*Feynrule_t::sizes[L];				\
    }															\
    return;														\
}															\

	    /*Recursion relations: */

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

    /* Specialisation of the cfd_evaluate class template for double Kronecker delta
     * colour structures over SU(N)-adjoint representation indices: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename cfd_evaluate<Feynrule_t>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_t>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_t>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_t>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Initialisation phase, to compute the necessary subtensor-sizes: */

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_t>::initialise();
	    }

	    /* Checking function, adding representation ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(n)														\
value_type z=(r_value_type)2*factor;												\
if(I==n)															\
{																\
    size_type l=(iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N);								\
    size_type AK=l/N;														\
    size_type BK=l%N;														\
    l=(iters[L].get_offset()/Feynrule_t::sizes[L])%(N*N);									\
    bool p=(AK==l%N and BK==l/N);												\
    bool q=(AK==BK and l%N==l/N);												\
    if(p or q)															\
    {																\
	if(p and q)														\
	{															\
	    z*=((r_value_type)1-N_inv);												\
	}															\
	if(!p and q)														\
	{															\
	    z*=(-N_inv);													\
	}															\
    }																\
    else															\
    {																\
	return;															\
    }																\
    l=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);									\
    size_type AJ=l/N;														\
    size_type BJ=l%N;														\
    if(AJ==BJ)															\
    {																\
	z*=(-N_inv);														\
	size_type step=(N+1)*Feynrule_t::sizes[I];										\
	for(size_type a=0;a<AJ;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[I]+=step;													\
	}															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-r_value_type(N-1)*z,couplings,iters,momenta,produced_iters);			\
	iters[I]+=step;														\
	for(size_type a=AJ+1;a<N;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[I]+=step;													\
	}															\
	iters[I]-=(N*step);													\
	return;															\
    }																\
    iters[I]+=l*Feynrule_t::sizes[I];												\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
    iters[I]-=l*Feynrule_t::sizes[I];												\
    return;															\
}																\
if(J==n)															\
{																\
    size_type l=(iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N);								\
    size_type AK=l/N;														\
    size_type BK=l%N;														\
    l=(iters[L].get_offset()/Feynrule_t::sizes[L])%(N*N);									\
    bool p=(AK==l%N and BK==l/N);												\
    bool q=(AK==BK and l%N==l/N);												\
    if(p or q)															\
    {																\
	if(p and q)														\
	{															\
	    z*=((r_value_type)1-N_inv);												\
	}															\
	if(!p and q)														\
	{															\
	    z*=(-N_inv);													\
	}															\
    }																\
    else															\
    {																\
	return;															\
    }																\
    l=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);									\
    size_type AI=l/N;														\
    size_type BI=l%N;														\
    if(AI==BI)															\
    {																\
	z*=(-N_inv);														\
	size_type step=(N+1)*Feynrule_t::sizes[J];										\
	for(size_type a=0;a<AI;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[J]+=step;													\
	}															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-r_value_type(N-1)*z,couplings,iters,momenta,produced_iters);			\
	iters[J]+=step;														\
	for(size_type a=AI+1;a<N;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[J]+=step;													\
	}															\
	iters[J]-=(N*step);													\
	return;															\
    }																\
    iters[J]+=l*Feynrule_t::sizes[J];												\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
    iters[J]-=l*Feynrule_t::sizes[J];												\
    return;															\
}																\
if(K==n)															\
{																\
    size_type l=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);								\
    size_type AI=l/N;														\
    size_type BI=l%N;														\
    l=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);									\
    bool p=(AI==l%N and BI==l/N);												\
    bool q=(AI==BI and l%N==l/N);												\
    if(p or q)															\
    {																\
	if(p and q)														\
	{															\
	    z*=((r_value_type)1-N_inv);												\
	}															\
	if(!p and q)														\
	{															\
	    z*=(-N_inv);													\
	}															\
    }																\
    else															\
    {																\
	return;															\
    }																\
    l=(iters[L].get_offset()/Feynrule_t::sizes[L])%(N*N);									\
    size_type AL=l/N;														\
    size_type BL=l%N;														\
    if(AL==BL)															\
    {																\
	z*=(-N_inv);														\
	size_type step=(N+1)*Feynrule_t::sizes[K];										\
	for(size_type a=0;a<AL;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[K]+=step;													\
	}															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-r_value_type(N-1)*z,couplings,iters,momenta,produced_iters);			\
	iters[K]+=step;														\
	for(size_type a=AL+1;a<N;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[K]+=step;													\
	}															\
	iters[K]-=(N*step);													\
	return;															\
    }																\
    iters[K]+=l*Feynrule_t::sizes[K];												\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
    iters[K]-=l*Feynrule_t::sizes[K];												\
    return;															\
}																\
if(L==n)															\
{																\
    size_type l=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);								\
    size_type AI=l/N;														\
    size_type BI=l%N;														\
    l=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);									\
    bool p=(AI==l%N and BI==l/N);												\
    bool q=(AI==BI and l%N==l/N);												\
    if(p or q)															\
    {																\
	if(p and q)														\
	{															\
	    z*=((r_value_type)1-N_inv);												\
	}															\
	if(!p and q)														\
	{															\
	    z*=(-N_inv);													\
	}															\
    }																\
    else															\
    {																\
	return;															\
    }																\
    l=(iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N);									\
    size_type AK=l/N;														\
    size_type BK=l%N;														\
    if(AK==BK)															\
    {																\
	z*=(-N_inv);														\
	size_type step=(N+1)*Feynrule_t::sizes[L];										\
	for(size_type a=0;a<AK;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[L]+=step;													\
	}															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-r_value_type(N-1)*z,couplings,iters,momenta,produced_iters);			\
	iters[L]+=step;														\
	for(size_type a=AK+1;a<N;++a)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[L]+=step;													\
	}															\
	iters[L]-=(N*step);													\
	return;															\
    }																\
    iters[L]+=l*Feynrule_t::sizes[L];												\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
    iters[L]-=l*Feynrule_t::sizes[L];												\
    return;															\
}																\

	    /* Recursion relations: */

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
	
	    /* Utility constant: */

	    static const r_value_type N_inv;
    };

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::r_value_type cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::N_inv=(typename cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::r_value_type)1/(typename cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::r_value_type)N;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_DD_H_*/

