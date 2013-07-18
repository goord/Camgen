//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_TT_CONTR_H_
#define CAMGEN_TT_CONTR_H_

#include <Camgen/eval.h>
#include <Camgen/comp_vert.h>
#include <Camgen/su(n).h>
#include <Camgen/TT_contr_helper.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the contracted group generators colour  *
 * structures.  Included is also the specialisation of the evaluate and  *
 * cfd_evaluate class templates for vertices composed with this colour   *
 * structure.                                                            *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	/* Definition of the generator product colour structure in generic
	 * representations: */

	template<class rep_t,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class TT_contr
	{
	    public:

		/* Self-referring type definition: */

		typedef TT_contr<rep_t,I,J,K,L> type;
		
		/* Colour-flow representation counterpart: */
		
		typedef TT_contr<rep_t,I,J,K,L> CF_type;

		/* Rank of the colour vertex: */

		static const std::size_t rank=4;
		
		/* Total size of the colour tensor: */
		
		static const std::size_t tensor_size=(rep_t::dimension)*(rep_t::dimension)*(rep_t::dimension)*(rep_t::dimension);

		/* Contracted fields: */

		static const std::size_t contractions[4];

		/* Colour index ranges: */

		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<rep_t,I,J,K,L>::rank;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<rep_t,I,J,K,L>::tensor_size;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<rep_t,I,J,K,L>::multiplicity;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<rep_t,I,J,K,L>::contractions[4]={I,J,K,L};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<rep_t,I,J,K,L>::ranges[4]={rep_t::dimension,rep_t::dimension,rep_t::dimension,rep_t::dimension};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string TT_contr<rep_t,I,J,K,L>::formula="T(a,A,B)T(a,C,D)";

	/* Specialisation of the colour structure for the SU(N) fundamental
	 * representation: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class TT_contr<fundamental_rep< SU<N> >,I,J,K,L>
	{
	    public:

		/* Self-referring type definition: */

		typedef TT_contr<fundamental_rep< SU<N> >,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef TT_contr<fundamental_rep< SU<N> >,I,J,K,L> CF_type;

		/* Rank of the colour vertex: */

		static const std::size_t rank=4;
		
		/* Size of the colour tensor: */
		
		static const std::size_t tensor_size=N*N*N*N;

		/* Contracted indices: */

		static const std::size_t contractions[4];
		
		/* Contracted index ranges: */
		
		static const std::size_t ranges[4];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<fundamental_rep< SU<N> >,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<fundamental_rep< SU<N> >,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<fundamental_rep< SU<N> >,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<fundamental_rep< SU<N> >,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t TT_contr<fundamental_rep< SU<N> >,I,J,K,L>::ranges[4]={N,N,N,N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string TT_contr<fundamental_rep< SU<N> >,I,J,K,L>::formula="1/2*(d(A1,B4)d(A3,B2)-d(A1,B4)d(A3,B2)/N)";
    }

    /* Specialisation of the evaluate class template for vertices composed with
     * the {T,T}-colour structure in generic representations: */

    template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::TT_contr<rep_t,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::TT_contr<rep_t,I,J,K,L>,Feynrule_t> vertex_type;
	    
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
		if(n==I or n==J or n==K or n==L)
		{
		    v.insert(v.begin(),rep_t::dimensions);
		}
		return v;
	    }

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
		colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::initialise();
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(n)																	\
for(size_type i=0;i<colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::size();++i)								\
{																			\
    iters[I]+=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::jump(0,i);									\
    iters[J]+=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::jump(1,i);									\
    iters[K]+=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::jump(2,i);									\
    iters[L]+=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::jump(3,i);									\
    evaluate<Feynrule_t>::RANK_NUMERAL(n)(colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::value(i)*factor,couplings,iters,momenta);		\
}																			\
iters[I]-=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::reset(0);										\
iters[J]-=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::reset(1);										\
iters[K]-=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::reset(2);										\
iters[L]-=colour_tensor::TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::reset(3)

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
    
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> vertex_type;

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
		if(n==I or n==J or n==K or n==L)
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

	    /* Convenient preprocessor macro definitions: */

#define CHAIN(n,z,I,J,K,L)									\
for(size_type a=0;a<N;++a)									\
{												\
    for(size_type b=0;b<N;++b)									\
    {												\
	evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);			\
	iters[K]+=sizes[0][K];									\
	iters[L]+=sizes[0][L];									\
    }												\
    iters[I]+=sizes[0][I];									\
    iters[J]+=sizes[0][J];									\
    iters[K]-=sizes[1][K];									\
    iters[L]-=sizes[1][L];									\
}												\
iters[I]-=sizes[1][I];										\
iters[J]-=sizes[1][J]

#define RECREL(n)										\
if(I==n)											\
{    												\
    value_type z=(r_value_type)0.5*factor;							\
    CHAIN(n,z,I,L,J,K);										\
    z=-z/(r_value_type)N;									\
    CHAIN(n,z,I,J,K,L);										\
    return;											\
}												\
if(J==n)											\
{												\
    value_type z=(r_value_type)0.5*factor;							\
    CHAIN(n,z,J,K,I,L);										\
    z=-z/(r_value_type)N;									\
    CHAIN(n,z,J,I,K,L);										\
    return;											\
}												\
if(K==n)											\
{												\
    value_type z=(r_value_type)0.5*factor;							\
    CHAIN(n,z,K,J,I,L);										\
    z=-z/(r_value_type)N;									\
    CHAIN(n,z,K,L,I,J);										\
    return;											\
}												\
if(L==n)											\
{												\
    value_type z=(r_value_type)0.5*factor;							\
    CHAIN(n,z,L,I,J,K);										\
    z=-z/(r_value_type)N;									\
    CHAIN(n,z,L,K,I,J);										\
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
#undef CHAIN
	
	private:

	    /* Interacting (sub)tensor size matrix: */

	    static const size_type sizes[2][4];
    };
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::sizes[2][4]={{Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]},{N*Feynrule_t::sizes[0],N*Feynrule_t::sizes[1],N*Feynrule_t::sizes[2],N*Feynrule_t::sizes[3]}};
    
    /* Specialisation of the cfd_evaluate class template for the colour-flow
     * Feynman rule version: */

template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> vertex_type;

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
		if(n==I or n==J or n==K or n==L)
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

#define RECREL(n)																		\
if(I==n)																			\
{																				\
    size_type AJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;												\
    size_type AK=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;												\
    size_type AL=(iters[L].get_offset()/Feynrule_t::sizes[L])%N;												\
    if(AJ==AK and AK!=AL)																	\
    {																				\
	iters[I]+=(AL*Feynrule_t::sizes[I]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor,couplings,iters,momenta,produced_iters);						\
	iters[I]-=(AL*Feynrule_t::sizes[I]);															\
	return;																			\
    }																				\
    if(AK==AL and AK!=AJ)																	\
    {																				\
	iters[I]+=(AJ*Feynrule_t::sizes[I]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);					\
	iters[I]-=(AJ*Feynrule_t::sizes[I]);															\
	return;																			\
    }																				\
    if(AJ==AK and AK==AL)																	\
    {																				\
	iters[I]+=(AJ*Feynrule_t::sizes[I]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);		\
	iters[I]-=(AJ*Feynrule_t::sizes[I]);															\
	return;																			\
    }																				\
}																				\
if(J==n)																			\
{																				\
    size_type AI=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;												\
    size_type AK=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;												\
    size_type AL=(iters[L].get_offset()/Feynrule_t::sizes[L])%N;												\
    if(AI==AL and AL!=AK)																	\
    {																				\
	iters[J]+=(AK*Feynrule_t::sizes[J]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor,couplings,iters,momenta,produced_iters);						\
	iters[J]-=(AK*Feynrule_t::sizes[J]);															\
	return;																			\
    }																				\
    if(AL==AK and AL!=AI)																	\
    {																				\
	iters[J]+=(AI*Feynrule_t::sizes[J]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);					\
	iters[J]-=(AI*Feynrule_t::sizes[J]);															\
	return;																			\
    }																				\
    if(AI==AL and AL==AK)																	\
    {																				\
	iters[J]+=(AI*Feynrule_t::sizes[J]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);		\
	iters[J]-=(AI*Feynrule_t::sizes[J]);															\
	return;																			\
    }																				\
}																				\
if(K==n)																			\
{																				\
    size_type AI=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;												\
    size_type AJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;												\
    size_type AL=(iters[L].get_offset()/Feynrule_t::sizes[L])%N;												\
    if(AI==AL and AL!=AJ)																	\
    {																				\
	iters[K]+=(AJ*Feynrule_t::sizes[K]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor,couplings,iters,momenta,produced_iters);						\
	iters[K]-=(AJ*Feynrule_t::sizes[K]);															\
	return;																			\
    }																				\
    if(AI==AJ and AI!=AL)																	\
    {																				\
	iters[K]+=(AL*Feynrule_t::sizes[K]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);					\
	iters[K]-=(AL*Feynrule_t::sizes[K]);															\
	return;																			\
    }																				\
    if(AI==AJ and AJ==AL)																	\
    {																				\
	iters[K]+=(AI*Feynrule_t::sizes[K]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);		\
	iters[K]-=(AI*Feynrule_t::sizes[K]);															\
	return;																			\
    }																				\
}																				\
if(L==n)																			\
{																				\
    size_type AI=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;												\
    size_type AJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;												\
    size_type AK=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;												\
    if(AJ==AK and AJ!=AI)																	\
    {																				\
	iters[L]+=(AI*Feynrule_t::sizes[L]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor,couplings,iters,momenta,produced_iters);						\
	iters[L]-=(AI*Feynrule_t::sizes[L]);															\
	return;																			\
    }																				\
    if(AI==AJ and AI!=AK)																	\
    {																				\
	iters[L]+=(AK*Feynrule_t::sizes[L]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);					\
	iters[L]-=(AK*Feynrule_t::sizes[L]);															\
	return;																			\
    }																				\
    if(AI==AJ and AJ==AK)																	\
    {																				\
	iters[L]+=(AI*Feynrule_t::sizes[L]);															\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)((r_value_type)0.5*factor-(r_value_type)0.5*N_inv*factor,couplings,iters,momenta,produced_iters);		\
	iters[L]-=(AI*Feynrule_t::sizes[L]);															\
	return;																			\
    }																				\
}																				\

	    /* Recursive relations: */

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
	    
	    /* Inverse of the number of colours: */

	    static const r_value_type N_inv;
    };
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>const typename cfd_evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::r_value_type cfd_evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::N_inv=(typename cfd_evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::r_value_type)1/(typename cfd_evaluate<compose_vertices<colour_tensor::TT_contr<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::r_value_type)N;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_TT_CONTR_H_*/
