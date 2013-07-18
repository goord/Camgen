//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FF_CONTR_H_
#define CAMGEN_FF_CONTR_H_

#include <Camgen/comp_vert.h>
#include <Camgen/eval.h>
#include <Camgen/ff_contr_helper.h>
#include <Camgen/adj_rep.h>
#include <Camgen/su(n).h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition and declaration of the contracted structure constants colour       *
 * tensor                                                                        *
 *                                                                               *
 * 			    f^{a,b}_{e}f^{e,c,d}                                 *
 *                                                                               *
 * and a specialisation of the evaluate and cfd_evaluate class templates for     *
 * vertices composed with the above colour structure. The is also a separate     *
 * spacialisation for a composition with the Yang-Mills 4-vertex, which cannot   *
 * be written as a simple tensor product of colour and a spacetime Feynman rule. *
 * This specialisation is implemented in the file gggg.h, included at the        *
 * bottom.                                                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	/* Contracted structure constants class for generic group indices: */

	template<class group_t,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class ff_contr
	{
	    public:

		/* Self-referring type definition: */

		typedef ff_contr<group_t,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef ff_contr<group_t,I,J,K,L> CF_type;

		/* Number of indices: */

		static const std::size_t rank=4;
		
		/* Total colour tensor size: */
		
		static const std::size_t tensor_size=(group_t::rank)*(group_t::rank)*(group_t::rank)*(group_t::rank);
		
		/* Contracted indices: */
		
		static const std::size_t contractions[4];

		/* Ranges of contracted indices: */

		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<class group_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<group_t,I,J,K,L>::rank;
	template<class group_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<group_t,I,J,K,L>::tensor_size;
	template<class group_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<group_t,I,J,K,L>::multiplicity;
	template<class group_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<group_t,I,J,K,L>::contractions[4]={I,J,K,L};
	template<class group_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<group_t,I,J,K,L>::ranges[4]={group_t::rank,group_t::rank,group_t::rank,group_t::rank};
	template<class group_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string ff_contr<group_t,I,J,K,L>::formula="f(a1,a2,e)f(a3,a4,e)";

	/* Declaration of the colour-flow structure constant class: */

	template<std::size_t N,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class CF_ff_contr;

	/* Specialisation of the contracted structure constants class for SU(N): */
	
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class ff_contr<SU<N>,I,J,K,L>
	{
	    public:

		/* Self-referring type definition: */

		typedef ff_contr<SU<N>,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef CF_ff_contr<N,I,J,K,L> CF_type;

		/* Number of indices: */

		static const std::size_t rank=4;
		
		/* Total colour tensor size: */
		
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
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<SU<N>,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<SU<N>,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<SU<N>,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<SU<N>,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t ff_contr<SU<N>,I,J,K,L>::ranges[4]={N*N-1,N*N-1,N*N-1,N*N-1};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string ff_contr<SU<N>,I,J,K,L>::formula="f(a1,a2,e)f(a3,a4,e)";

	/* Definition of the colour-flow counterpart: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class CF_ff_contr
	{
	    public:

		/* Self-referring type definition: */

		typedef CF_ff_contr<N,I,J,K,L> type;

		/* Number of indices: */

		static const std::size_t rank=4;
		
		/* Total colour tensor size: */
		
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
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_ff_contr<N,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_ff_contr<N,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_ff_contr<N,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_ff_contr<N,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_ff_contr<N,I,J,K,L>::ranges[4]={N*N,N*N,N*N,N*N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string CF_ff_contr<N,I,J,K,L>::formula="2(d(A1,B2)d(A2,B4)d(A3,B1)d(A4,B3)-d(A1,B2)d(A2,B3)d(A3,B4)d(A4,B1)+d(A1,B3)d(A2,B1)d(A3,B4)d(A4,B2)-d(A1,B4)d(A2,B1)d(A3,B2)d(A4,B3))";
    }

    /* Specialisation of the evaluate class template for vertices composed with
     * contracted double structure constants: */

    template<class group_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::ff_contr<group_t,I,J,K,L>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::ff_contr<group_t,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Group type definition: */

	    DEFINE_GROUP_TYPE(model_type,group_t);

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
		colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::initialise();
	    }

	    /* Checking function, adding group index ranges to the
	     * contracted indices in the integer vector: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),group_t::rank);
		}
		return v;
	    }

	    /* Convenient preprocessor macro definition, evaluating the recursive
	     * relations using the helper class: */

#define RECREL(n)																	\
for(size_type i=0;i<colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::size();++i)								\
{																			\
    iters[I]+=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::jump(0,i);									\
    iters[J]+=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::jump(1,i);									\
    iters[K]+=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::jump(2,i);									\
    iters[L]+=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::jump(3,i);									\
    evaluate<Feynrule_t>::RANK_NUMERAL(n)(colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::value(i)*factor,couplings,iters,momenta);		\
}																			\
iters[I]-=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::reset(0);										\
iters[J]-=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::reset(1);										\
iters[K]-=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::reset(2);										\
iters[L]-=colour_tensor::ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::reset(3);										\

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

    /* Specialisation of the evaluate class template for vertices composed with the
     * colour-flow counterpart of the ff_contr colour tensor: */
    
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_ff_contr<N,I,J,K,L>,Feynrule_t> vertex_type;
	    
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
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(n)				\
if(I==n)					\
{						\
    GGG_to_G_CHAIN(I,J,K,L,n,-factor);		\
    GGG_to_G_CHAIN(I,J,L,K,n,factor);		\
    GGG_to_G_CHAIN(I,K,L,J,n,factor);		\
    GGG_to_G_CHAIN(I,L,K,J,n,-factor);		\
    return;					\
}						\
if(J==n)					\
{						\
    GGG_to_G_CHAIN(J,I,K,L,n,factor);		\
    GGG_to_G_CHAIN(J,I,L,K,n,-factor);		\
    GGG_to_G_CHAIN(J,K,L,I,n,-factor);		\
    GGG_to_G_CHAIN(J,L,K,I,n,factor);		\
    return;					\
}						\
if(K==n)					\
{						\
    GGG_to_G_CHAIN(K,L,I,J,n,-factor);		\
    GGG_to_G_CHAIN(K,L,J,I,n,factor);		\
    GGG_to_G_CHAIN(K,I,J,L,n,factor);		\
    GGG_to_G_CHAIN(K,J,I,L,n,-factor);		\
    return;					\
}						\
if(L==n)					\
{						\
    GGG_to_G_CHAIN(L,K,I,J,n,factor);		\
    GGG_to_G_CHAIN(L,K,J,I,n,-factor);		\
    GGG_to_G_CHAIN(L,I,J,K,n,-factor);		\
    GGG_to_G_CHAIN(L,J,I,K,n,factor);		\
    return;					\
}						\
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

	    /* Subtensor sizes: */

	    static const size_type sizes[3][4];
    };

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename evaluate< compose_vertices<colour_tensor::CF_ff_contr<N,I,J,K,L>,Feynrule_t> >::size_type evaluate< compose_vertices<colour_tensor::CF_ff_contr<N,I,J,K,L>,Feynrule_t> >::sizes[3][4]={{Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]},{N*Feynrule_t::sizes[0],N*Feynrule_t::sizes[1],N*Feynrule_t::sizes[2],N*Feynrule_t::sizes[3]},{N*N*Feynrule_t::sizes[0],N*N*Feynrule_t::sizes[1],N*N*Feynrule_t::sizes[2],N*N*Feynrule_t::sizes[3]}};

    /* Specialisation of the cfd_evaluate class template for Feynman rule
     * classes composed with structure constants: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_ff_contr<N,I,J,K,L>,Feynrule_t> vertex_type;

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
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(n)												\
if(I==n)													\
{														\
    size_type SJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);						\
    size_type AJ=SJ/N;												\
    size_type BJ=SJ%N;												\
    size_type SK=(iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N);						\
    size_type AK=SK/N;												\
    size_type BK=SK%N;												\
    size_type SL=(iters[L].get_offset()/Feynrule_t::sizes[L])%(N*N);						\
    size_type AL=SL/N;												\
    size_type BL=SL%N;												\
    if(AK==BL and BK==AL and AK==BK)										\
    {														\
	return;													\
    }														\
    if(AK==BL)													\
    {														\
	if(AJ==BK)												\
	{													\
	    iters[I]+=(AL*N+BJ)*Feynrule_t::sizes[I];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[I]-=(AL*N+BJ)*Feynrule_t::sizes[I];								\
	}													\
	if(BJ==AL)												\
	{													\
	    iters[I]+=(AJ*N+BK)*Feynrule_t::sizes[I];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[I]-=(AJ*N+BK)*Feynrule_t::sizes[I];								\
	}													\
    }														\
    if(BK==AL)													\
    {														\
	if(AJ==BL)												\
	{													\
	    iters[I]+=(AK*N+BJ)*Feynrule_t::sizes[I];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[I]-=(AK*N+BJ)*Feynrule_t::sizes[I];								\
	}													\
	if(BJ==AK)												\
	{													\
	    iters[I]+=(AJ*N+BL)*Feynrule_t::sizes[I];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[I]-=(AJ*N+BL)*Feynrule_t::sizes[I];								\
	}													\
    }														\
    return;													\
}														\
if(J==n)													\
{														\
    size_type SI=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);						\
    size_type AI=SI/N;												\
    size_type BI=SI%N;												\
    size_type SK=(iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N);						\
    size_type AK=SK/N;												\
    size_type BK=SK%N;												\
    size_type SL=(iters[L].get_offset()/Feynrule_t::sizes[L])%(N*N);						\
    size_type AL=SL/N;												\
    size_type BL=SL%N;												\
    if(AK==BL and BK==AL and AK==BK)										\
    {														\
	return;													\
    }														\
    if(AK==BL)													\
    {														\
	if(AI==BK)												\
	{													\
	    iters[J]+=(AL*N+BI)*Feynrule_t::sizes[J];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[J]-=(AL*N+BI)*Feynrule_t::sizes[J];								\
	}													\
	if(BI==AL)												\
	{													\
	    iters[J]+=(AI*N+BK)*Feynrule_t::sizes[J];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[J]-=(AI*N+BK)*Feynrule_t::sizes[J];								\
	}													\
    }														\
    if(BK==AL)													\
    {														\
	if(AI==BL)												\
	{													\
	    iters[J]+=(AK*N+BI)*Feynrule_t::sizes[J];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[J]-=(AK*N+BI)*Feynrule_t::sizes[J];								\
	}													\
	if(BI==AK)												\
	{													\
	    iters[J]+=(AI*N+BL)*Feynrule_t::sizes[J];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[J]-=(AI*N+BL)*Feynrule_t::sizes[J];								\
	}													\
    }														\
    return;													\
}														\
if(K==n)													\
{														\
    size_type SI=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);						\
    size_type AI=SI/N;												\
    size_type BI=SI%N;												\
    size_type SJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);						\
    size_type AJ=SJ/N;												\
    size_type BJ=SJ%N;												\
    size_type SL=(iters[L].get_offset()/Feynrule_t::sizes[L])%(N*N);						\
    size_type AL=SL/N;												\
    size_type BL=SL%N;												\
    if(AI==BJ and BI==AJ and AI==BI)										\
    {														\
	return;													\
    }														\
    if(AI==BJ)													\
    {														\
	if(AL==BI)												\
	{													\
	    iters[K]+=(AJ*N+BL)*Feynrule_t::sizes[K];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[K]-=(AJ*N+BL)*Feynrule_t::sizes[K];								\
	}													\
	if(BL==AJ)												\
	{													\
	    iters[K]+=(AL*N+BI)*Feynrule_t::sizes[K];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[K]-=(AL*N+BI)*Feynrule_t::sizes[K];								\
	}													\
    }														\
    if(BI==AJ)													\
    {														\
	if(AL==BJ)												\
	{													\
	    iters[K]+=(AI*N+BL)*Feynrule_t::sizes[K];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[K]-=(AI*N+BL)*Feynrule_t::sizes[K];								\
	}													\
	if(BL==AI)												\
	{													\
	    iters[K]+=(AL*N+BJ)*Feynrule_t::sizes[K];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[K]-=(AL*N+BJ)*Feynrule_t::sizes[K];								\
	}													\
    }														\
    return;													\
}														\
if(L==n)													\
{														\
    size_type SI=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);						\
    size_type AI=SI/N;												\
    size_type BI=SI%N;												\
    size_type SJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);						\
    size_type AJ=SJ/N;												\
    size_type BJ=SJ%N;												\
    size_type SK=(iters[K].get_offset()/Feynrule_t::sizes[K])%(N*N);						\
    size_type AK=SK/N;												\
    size_type BK=SK%N;												\
    if(AI==BJ and BI==AJ and AI==BI)										\
    {														\
	return;													\
    }														\
    if(AI==BJ)													\
    {														\
	if(AK==BI)												\
	{													\
	    iters[L]+=(AJ*N+BK)*Feynrule_t::sizes[L];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[L]-=(AJ*N+BK)*Feynrule_t::sizes[L];								\
	}													\
	if(BK==AJ)												\
	{													\
	    iters[L]+=(AK*N+BI)*Feynrule_t::sizes[L];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[L]-=(AK*N+BI)*Feynrule_t::sizes[L];								\
	}													\
    }														\
    if(BI==AJ)													\
    {														\
	if(AK==BJ)												\
	{													\
	    iters[L]+=(AI*N+BK)*Feynrule_t::sizes[L];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-factor,couplings,iters,momenta,produced_iters);		\
	    iters[L]-=(AI*N+BK)*Feynrule_t::sizes[L];								\
	}													\
	if(BK==AI)												\
	{													\
	    iters[L]+=(AK*N+BJ)*Feynrule_t::sizes[L];								\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);		\
	    iters[L]-=(AK*N+BJ)*Feynrule_t::sizes[L];								\
	}													\
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

/* Specialisation for the composition with the Yang-Mills four-vertex: */

#ifdef CAMGEN_VVVV_H_
#include <Camgen/gggg.h>
#endif /*CAMGEN_VVVV_H_*/

#include <Camgen/undef_args.h>

#endif /*CAMGEN_FF_CONTR_H_*/
