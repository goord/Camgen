//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_TT_H_
#define CAMGEN_TT_H_

#include <Camgen/eval.h>
#include <Camgen/comp_vert.h>
#include <Camgen/su(n).h>
#include <Camgen/TT_helper.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the group generator products colour structures. *
 * Included is also the specialisation of the evaluate and cfd_evaluate class    *
 * templates for vertices composed with this colour structure.                   *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	/* Definition of the generator product colour structure in generic
	 * representations: */

	template<class rep_t,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class TT
	{
	    public:

		/* Self-referring type definition: */

		typedef TT<rep_t,I,J,K,L> type;
		
		/* Colour-flow representation counterpart: */
		
		typedef TT<rep_t,I,J,K,L> CF_type;

		/* Rank of the colour vertex: */

		static const std::size_t rank=4;
		
		/* Total size of the colour tensor: */
		
		static const std::size_t tensor_size=(rep_t::group_type::rank)*(rep_t::group_type::rank)*(rep_t::dimension)*(rep_t::dimension);

		/* Contracted fields: */

		static const std::size_t contractions[4];

		/* Colour index ranges: */

		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<rep_t,I,J,K,L>::rank;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<rep_t,I,J,K,L>::tensor_size;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<rep_t,I,J,K,L>::multiplicity;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<rep_t,I,J,K,L>::contractions[4]={I,J,K,L};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<rep_t,I,J,K,L>::ranges[4]={rep_t::group_type::rank,rep_t::group_type::rank,rep_t::dimension,rep_t::dimension};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string TT<rep_t,I,J,K,L>::formula="T(a,A,C)T(b,C,B)";

	/* Declaration of the colour-flow mode counterpart: */

	template<std::size_t N,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class CF_TT;

	/* Specialisation of the colour structure for the SU(N) fundamental
	 * representation: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class TT<fundamental_rep< SU<N> >,I,J,K,L>
	{
	    public:

		/* Self-referring type definition: */

		typedef TT<fundamental_rep< SU<N> >,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef CF_TT<N,I,J,K,L> CF_type;

		/* Rank of the colour vertex: */

		static const std::size_t rank=4;
		
		/* Size of the colour tensor: */
		
		static const std::size_t tensor_size=(N*N-1)*(N*N-1)*N*N;

		/* Contracted indices: */

		static const std::size_t contractions[4];
		
		/* Contracted index ranges: */
		
		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<fundamental_rep< SU<N> >,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<fundamental_rep< SU<N> >,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<fundamental_rep< SU<N> >,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<fundamental_rep< SU<N> >,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT<fundamental_rep< SU<N> >,I,J,K,L>::ranges[4]={N*N-1,N*N-1,N,N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string TT<fundamental_rep< SU<N> >,I,J,K,L>::formula="T(a,A,C)T(b,C,B)";

	/* Definition of th colour-flow mode TT-vertex colour structure: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class CF_TT
	{
	    public:

		/* Self-referring type definition: */

		typedef CF_TT<N,I,J,K,L> type;

		/* Rank of the colour vertex: */

		static const std::size_t rank=4;
		
		/* Size of the colour tensor: */
		
		static const std::size_t tensor_size=N*N*N*N*N*N;

		/* Contracted fields: */

		static const std::size_t contractions[4];
		
		/* Contracted index ranges: */
		
		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_TT<N,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_TT<N,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_TT<N,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_TT<N,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_TT<N,I,J,K,L>::ranges[4]={N*N,N*N,N,N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string CF_TT<N,I,J,K,L>::formula="(d(A1,B2)d(A2,B4)d(A3,B1)-d(A1,B1)d(A2,B4)d(A3,B2)/N-d(A1,B4)d(A2,B2)d(A3,B1)/N+d(A1,B1)d(A2,B2)d(A3,B4)/N^2)";
    }

    /* Specialisation of the evaluate class template for vertices composed with
     * the TT-colour structure in generic representations: */

    template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::TT<rep_t,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::TT<rep_t,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* The usual type definitions: */
	    
	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	private:

	    /* Representation type definition: */

	    DEFINE_REP_TYPE(model_type,rep_t);

	public:

	    /* Checking function returning the index ranges of interacting
	     * subamplitude tensors: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J)
		{
		    v.insert(v.begin(),rep_t::group_type::rank);
		}
		if(n==K or n==L)
		{
		    v.insert(v.begin(),rep_t::dimension);
		}
		return v;
	    }

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
		colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::initialise();
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(n)																\
for(size_type i=0;i<colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::size();++i)								\
{																		\
    iters[I]+=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::jump(0,i);									\
    iters[J]+=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::jump(1,i);									\
    iters[K]+=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::jump(2,i);									\
    iters[L]+=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::jump(3,i);									\
    evaluate<Feynrule_t>::RANK_NUMERAL(n)(colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::value(i)*factor,couplings,iters,momenta);		\
}																		\
iters[I]-=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::reset(0);										\
iters[J]-=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::reset(1);										\
iters[K]-=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::reset(2);										\
iters[L]-=colour_tensor::TT_helper<rep_t,Feynrule_t,I,J,K,L>::reset(3)

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
    };

#undef RECREL

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> vertex_type;

	    /* The usual type definitions: */

	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Checking function returning the index ranges of interacting
	     * subamplitude tensors: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J)
		{
		    v.insert(v.begin(),2,N);
		}
		if(n==K or n==L)
		{
		    v.insert(v.begin(),N);
		}
		return v;
	    }

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(n)										\
if(I==n)											\
{												\
    value_type z=(r_value_type)0.5*factor/(r_value_type)N;					\
    white_partI.reset();									\
    iterator temp=iters[I];									\
    iters[I]=white_partI.begin();								\
    for(size_type c=0;c<N;++c)									\
    {												\
	for(size_type d=0;d<N;++d)								\
	{											\
	    evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);			\
	    iters[J]+=sizes[0][J];								\
	    iters[L]+=sizes[0][L];								\
	}											\
	iters[L]-=sizes[1][L];									\
	iters[K]+=sizes[0][K];									\
    }												\
    iters[J]-=sizes[2][J];									\
    iters[K]-=sizes[1][K];									\
    iters[I]=temp;										\
    z*=(r_value_type)N;										\
    for(size_type b=0;b<N;++b)									\
    {												\
	for(size_type a=0;a<b;++a)								\
	{											\
	    for(size_type c=0;c<N;++c)								\
	    {											\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);		\
		iters[J]+=sizes[0][J];								\
		iters[L]+=sizes[0][L];								\
	    }											\
	    iters[L]-=sizes[1][L];								\
	    iters[I]+=sizes[1][I];								\
	}											\
	for(size_type mu=0;mu<sizes[0][I];++mu)							\
	{											\
	    iters[I][mu]-=white_partI[mu];							\
	}											\
	for(size_type a=b;a<N;++a)								\
	{											\
	    for(size_type c=0;c<N;++c)								\
	    {											\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);		\
		iters[J]+=sizes[0][J];								\
		iters[L]+=sizes[0][L];								\
	    }											\
	    iters[L]-=sizes[1][L];								\
	    iters[I]+=sizes[1][I];								\
	}											\
	iters[I]-=sizes[2][I];									\
	iters[J]-=sizes[2][J];									\
	iters[I]+=sizes[0][I];									\
	iters[K]+=sizes[0][K];									\
    }												\
    iters[I]-=sizes[1][I];									\
    iters[K]-=sizes[1][K];									\
    return;											\
}												\
if(J==n)											\
{												\
    value_type z=(r_value_type)0.5*factor/(r_value_type)N;					\
    white_partJ.reset();									\
    iterator temp=iters[J];									\
    iters[J]=white_partJ.begin();								\
    for(size_type c=0;c<N;++c)									\
    {												\
	for(size_type d=0;d<N;++d)								\
	{											\
	    evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);			\
	    iters[I]+=sizes[0][I];								\
	    iters[L]+=sizes[0][L];								\
	}											\
	iters[K]+=sizes[0][K];									\
	iters[L]-=sizes[1][L];									\
    }												\
    iters[I]-=sizes[2][I];									\
    iters[K]-=sizes[1][K];									\
    iters[J]=temp;										\
    z*=(r_value_type)N;										\
    for(size_type a=0;a<N;++a)									\
    {												\
	for(size_type b=0;b<a;++b)								\
	{											\
	    for(size_type c=0;c<N;++c)								\
	    {											\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);		\
		iters[I]+=sizes[1][I];								\
		iters[K]+=sizes[0][K];								\
	    }											\
	    iters[I]-=sizes[2][I];								\
	    iters[K]-=sizes[1][K];								\
	    iters[I]+=sizes[0][I];								\
	    iters[J]+=sizes[0][J];								\
	}											\
	for(size_type mu=0;mu<sizes[0][J];++mu)							\
	{											\
	    iters[J][mu]-=white_partJ[mu];							\
	}											\
	for(size_type b=a;b<N;++b)								\
	{											\
	    for(size_type c=0;c<N;++c)								\
	    {											\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);		\
		iters[I]+=sizes[1][I];								\
		iters[K]+=sizes[0][K];								\
	    }											\
	    iters[I]-=sizes[2][I];								\
	    iters[K]-=sizes[1][K];								\
	    iters[I]+=sizes[0][I];								\
	    iters[J]+=sizes[0][J];								\
	}											\
	iters[I]-=sizes[1][I];									\
	iters[L]+=sizes[0][L];									\
    }												\
    iters[J]-=sizes[2][J];									\
    iters[L]-=sizes[1][L];									\
    return;											\
}												\
if(K==n)											\
{												\
    for(size_type a=0;a<N;++a)									\
    {												\
	for(size_type b=0;b<N;++b)								\
	{											\
	    for(size_type c=0;c<N;++c)								\
	    {											\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta);		\
		iters[J]+=sizes[0][J];								\
		iters[L]+=sizes[0][L];								\
	    }											\
	    iters[I]+=sizes[0][I];								\
	    iters[L]-=sizes[1][L];								\
	}											\
	iters[J]-=sizes[2][J];									\
	iters[K]+=sizes[0][K];									\
    }												\
    iters[I]-=sizes[2][I];									\
    iters[K]-=sizes[1][K];									\
    return;											\
}												\
if(L==n)											\
{												\
    for(size_type a=0;a<N;++a)									\
    {												\
	for(size_type c=0;c<N;++c)								\
	{											\
	    for(size_type b=0;b<N;++b)								\
	    {											\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta);		\
		iters[J]+=sizes[0][J];								\
		iters[L]+=sizes[0][L];								\
	    }											\
	    iters[I]+=sizes[0][I];								\
	    iters[L]-=sizes[1][L];								\
	}											\
	iters[J]-=sizes[2][J];									\
	iters[K]+=sizes[0][K];									\
    }												\
    iters[I]-=sizes[2][I];									\
    iters[K]-=sizes[1][K];									\
}												\
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

	private:

	    /* Interacting (sub)tensor size matrix: */

	    static const size_type sizes[3][4];

	    /* Temporary data storage tensors: */

	    static tensor_type white_partI,white_partJ;

	    /* Utility index range holder vector: */

	    static std::vector<size_type> utilvec;
    };
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>std::vector<typename evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::size_type> evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::utilvec;
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>typename evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::tensor_type evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::white_partI(evaluate<Feynrule_t>::get_index_ranges(I,evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::utilvec));
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>typename evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::tensor_type evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::white_partJ(evaluate<Feynrule_t>::get_index_ranges(J,evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::utilvec));
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::sizes[3][4]={{Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]},{N*Feynrule_t::sizes[0],N*Feynrule_t::sizes[1],N*Feynrule_t::sizes[2],N*Feynrule_t::sizes[3]},{N*N*Feynrule_t::sizes[0],N*N*Feynrule_t::sizes[1],N*N*Feynrule_t::sizes[2],N*N*Feynrule_t::sizes[3]}};
    
    /* Specialisation of the cfd_evaluate class template for the colour-flow
     * Feynman rule version: */

template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> vertex_type;

	    /* The usual type definitions: */

	    typedef typename cfd_evaluate<Feynrule_t>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_t>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_t>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_t>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Checking function returning the index ranges of interacting
	     * subamplitude tensors: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I or n==J)
		{
		    v.insert(v.begin(),N*N);
		}
		if(n==K or N==L)
		{
		    v.insert(v.begin(),N);
		}
		return v;
	    }

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(n)														\
if(I==n)															\
{																\
    size_type SJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);								\
    size_type AJ=SJ/N;														\
    size_type BJ=SJ%N;														\
    size_type AL=(iters[L].get_offset()/Feynrule_t::sizes[L])%N;								\
    value_type z=(r_value_type)0.5*factor;											\
    size_type Aeff;														\
    if(BJ==AL and BJ!=AJ)													\
    {																\
	Aeff=AJ;														\
    }																\
    else if(BJ==AJ and BJ!=AL)													\
    {																\
	Aeff=AL;														\
	z*=(-N_inv);														\
    }																\
    else if(BJ==AJ and BJ==AL)													\
    {																\
	Aeff=AJ;														\
	z-=(N_inv*z);														\
    }																\
    else															\
    {																\
	return;															\
    }																\
    size_type BK=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;								\
    if(BK==Aeff)														\
    {																\
	size_type step=(N+1)*Feynrule_t::sizes[I];										\
	z*=(-N_inv);														\
	for(size_type i=0;i<Aeff;++i)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[I]+=step;													\
	}															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z*value_type(N-1,0),couplings,iters,momenta,produced_iters);			\
	iters[I]+=step;														\
	for(size_type i=Aeff+1;i<N;++i)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[I]+=step;													\
	}															\
	iters[I]-=(N*step);													\
	return;															\
    }																\
    else															\
    {																\
	iters[I]+=(Aeff*N+BK)*Feynrule_t::sizes[I];										\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
	iters[I]-=(Aeff*N+BK)*Feynrule_t::sizes[I];										\
    }																\
    return;															\
}																\
if(J==n)															\
{																\
    size_type SI=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);								\
    size_type AI=SI/N;														\
    size_type BI=SI%N;														\
    size_type AK=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;								\
    value_type z=(r_value_type)0.5*factor;											\
    size_type Beff;														\
    if(AI==AK and AI!=BI)													\
    {																\
	Beff=BI;														\
    }																\
    else if(AI==BI and AI!=AK)													\
    {																\
	Beff=AK;														\
	z*=(-N_inv);														\
    }																\
    else if(AI==BI and AI==AK)													\
    {																\
	Beff=BI;														\
	z-=(N_inv*z);														\
    }																\
    else															\
    {																\
	return;															\
    }																\
    size_type AL=(iters[L].get_offset()/Feynrule_t::sizes[L])%N;								\
    if(AL==Beff)														\
    {																\
	size_type step=(N+1)*Feynrule_t::sizes[J];										\
	z*=(-N_inv);														\
	for(size_type i=0;i<Beff;++i)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[J]+=step;													\
	}															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z*value_type(N-1,0),couplings,iters,momenta,produced_iters);			\
	iters[J]+=step;														\
	for(size_type i=Beff+1;i<N;++i)												\
	{															\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);				\
	    iters[J]+=step;													\
	}															\
	iters[J]-=(N*step);													\
	return;															\
    }																\
    else															\
    {																\
	iters[J]+=(AL*N+Beff)*Feynrule_t::sizes[J];										\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
	iters[J]-=(AL*N+Beff)*Feynrule_t::sizes[J];										\
    }																\
    return;															\
}																\
if(K==n)															\
{																\
    size_type SJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);								\
    size_type AJ=SJ/N;														\
    size_type BJ=SJ%N;														\
    size_type AL=(iters[L].get_offset()/Feynrule_t::sizes[L])%N;								\
    size_type Aeff;														\
    value_type z=factor;													\
    if(BJ==AL and BJ!=AJ)													\
    {																\
	Aeff=AJ;														\
    }																\
    else if(BJ==AJ and BJ!=AL)													\
    {																\
	Aeff=AL;														\
	z*=(-N_inv);														\
    }																\
    else if(BJ==AL and BJ==AJ)													\
    {																\
	Aeff=AJ;														\
	z-=(z*N_inv);														\
    }																\
    else															\
    {																\
	return;															\
    }																\
    size_type SI=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);								\
    size_type AI=SI/N;														\
    size_type BI=SI%N;														\
    if(BI==Aeff and BI!=AI)													\
    {																\
	iters[K]+=(AI*Feynrule_t::sizes[K]);											\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
	iters[K]-=(AI*Feynrule_t::sizes[K]);											\
	return;															\
    }																\
    if(BI==AI and BI!=Aeff)													\
    {																\
	iters[K]+=(Aeff*Feynrule_t::sizes[K]);											\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-N_inv*z,couplings,iters,momenta,produced_iters);				\
	iters[K]-=(Aeff*Feynrule_t::sizes[K]);											\
	return;															\
    }																\
    if(BI==AI and BI==Aeff)													\
    {																\
	iters[K]+=(AI*Feynrule_t::sizes[K]);											\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z-N_inv*z,couplings,iters,momenta,produced_iters);				\
	iters[K]-=(AI*Feynrule_t::sizes[K]);											\
	return;															\
    }																\
}																\
if(L==n)															\
{																\
    size_type SI=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);								\
    size_type AI=SI/N;														\
    size_type BI=SI%N;														\
    size_type BK=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;								\
    size_type Beff;														\
    value_type z=factor;													\
    if(AI==BK and AI!=BI)													\
    {																\
	Beff=BI;														\
    }																\
    else if(AI==BI and AI!=BK)													\
    {																\
	Beff=BK;														\
	z*=(-N_inv);														\
    }																\
    else if(AI==BK and AI==BI)													\
    {																\
	Beff=BK;														\
	z-=(z*N_inv);														\
    }																\
    else															\
    {																\
	return;															\
    }																\
    size_type SJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);								\
    size_type AJ=SJ/N;														\
    size_type BJ=SJ%N;														\
    if(AJ==Beff and AJ!=BJ)													\
    {																\
	iters[L]+=(BJ*Feynrule_t::sizes[L]);											\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
	iters[L]-=(BJ*Feynrule_t::sizes[L]);											\
	return;															\
    }																\
    if(AJ==BJ and AJ!=Beff)													\
    {																\
	iters[L]+=(Beff*Feynrule_t::sizes[L]);											\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-N_inv*z,couplings,iters,momenta,produced_iters);				\
	iters[L]-=(Beff*Feynrule_t::sizes[L]);											\
	return;															\
    }																\
    if(AJ==BJ and AJ==Beff)													\
    {																\
	iters[L]+=(BJ*Feynrule_t::sizes[L]);											\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z-N_inv*z,couplings,iters,momenta,produced_iters);				\
	iters[L]-=(BJ*Feynrule_t::sizes[L]);											\
	return;															\
    }																\
}																\
return

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

	    /* Inverse number of colours: */

	    static const r_value_type N_inv;
    };
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename cfd_evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::r_value_type cfd_evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::N_inv=(typename cfd_evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::r_value_type)1/(typename cfd_evaluate<compose_vertices<colour_tensor::CF_TT<N,I,J,K,L>,Feynrule_t> >::r_value_type)N;

}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_TT_H_*/
