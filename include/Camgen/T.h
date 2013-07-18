//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_T_H_
#define CAMGEN_T_H_

#include <Camgen/eval.h>
#include <Camgen/comp_vert.h>
#include <Camgen/su(n).h>
#include <Camgen/T_helper.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Group generator colour structure declaration and definition. Included are     *
 * also the specialisations of the evaluate and cfd_evaluate class templates for *
 * vertices composed with the T-matrix colour structure.                         *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {

	/* Declaration and definition of the T-matrix colour structure for generic
	 * representations: */

	template<class rep_t,std::size_t I=0,std::size_t J=1,std::size_t K=2>class T
	{
	    public:

		/* Self-referring type definition: */

		typedef T<rep_t,I,J,K> type;

		/* Colour-flow mode counterpart: */

		typedef T<rep_t,I,J,K> CF_type;

		/* Number of indices: */

		static const std::size_t rank=3;
		
		/* Colour tensor size: */
		
		static const std::size_t tensor_size=(rep_t::group_type::rank)*(rep_t::dimension)*(rep_t::dimension);
		
		/* Contracted fields: */
		
		static const std::size_t contractions[3];

		/* Index ranges: */

		static const std::size_t ranges[3];

		/* Output to model logfile: */		

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<rep_t,I,J,K>::rank;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<rep_t,I,J,K>::tensor_size;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<rep_t,I,J,K>::multiplicity;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<rep_t,I,J,K>::contractions[3]={I,J,K};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<rep_t,I,J,K>::ranges[3]={rep_t::group_type::rank,rep_t::dimension,rep_t::dimension};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K>const std::string T<rep_t,I,J,K>::formula="T(a,A,B)";

	template<std::size_t N,std::size_t I=0,std::size_t J=1,std::size_t K=2>class CF_T;

	/* Specialisation of the T-matrix colour structure for the SU(N) fundamental
	 * representation: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>class T<fundamental_rep< SU<N> >,I,J,K>
	{
	    public:

		/* Self-referring type definition: */

		typedef T<fundamental_rep< SU<N> >,I,J,K> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef CF_T<N,I,J,K> CF_type;

		/* Number of indices: */

		static const std::size_t rank=3;
		
		/* Colour tensor size: */
		
		static const std::size_t tensor_size=(N*N-1)*N*N;
		
		/* Contracted fields: */
		
		static const std::size_t contractions[3];

		/* Index ranges: */

		static const std::size_t ranges[3];

		/* Output to model logfile: */		

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<fundamental_rep< SU<N> >,I,J,K>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<fundamental_rep< SU<N> >,I,J,K>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<fundamental_rep< SU<N> >,I,J,K>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<fundamental_rep< SU<N> >,I,J,K>::contractions[3]={I,J,K};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t T<fundamental_rep< SU<N> >,I,J,K>::ranges[3]={N*N-1,N,N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::string T<fundamental_rep< SU<N> >,I,J,K>::formula="T(a,A,B)";

	/* T-matrix colour structure in the colour-flow representation of the gluons: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>class CF_T
	{
	    public:

		/* Self-referring type definition: */

		typedef CF_T<N,I,J,K> type;
		
		/* Number of interacting fields: */
		
		static const std::size_t rank=3;

		/* Colour tensor size: */

		static const std::size_t tensor_size=N*N*N*N;

		/* Contracted fields: */

		static const std::size_t contractions[3];
		
		/* Index ranges: */
		
		static const std::size_t ranges[3];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_T<N,I,J,K>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_T<N,I,J,K>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_T<N,I,J,K>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_T<N,I,J,K>::contractions[3]={I,J,K};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::size_t CF_T<N,I,J,K>::ranges[3]={N*N,N,N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K>const std::string CF_T<N,I,J,K>::formula="(d(A1,A3)d(A2,B1)-d(A1,B1)d(A2,A3)/N)";
    }

    /* Specialisation of the evaluate class template for vertices composed with
     * the T-colour structure in generic representations: */

    template<class rep_t,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::T<rep_t,I,J,K>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::T<rep_t,I,J,K>,Feynrule_t> vertex_type;
	    
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
		if(n==I)
		{
		    v.insert(v.begin(),rep_t::group_type::rank);
		}
		if(n==J or n==K)
		{
		    v.insert(v.begin(),rep_t::dimension);
		}
		return v;
	    }

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<Feynrule_t>::initialise();
		colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::initialise();
	    }

	    /* Useful preprocessor macro definition: */

#define RECREL(n)															\
for(size_type i=0;i<colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::size();++i)							\
{																	\
    iters[I]+=colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::jump(0,i);								\
    iters[J]+=colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::jump(1,i);								\
    iters[K]+=colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::jump(2,i);								\
    evaluate<Feynrule_t>::RANK_NUMERAL(n)(colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::value(i)*factor,couplings,iters,momenta);	\
}																	\
iters[I]-=colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::reset(0);									\
iters[J]-=colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::reset(1);									\
iters[K]-=colour_tensor::T_helper<rep_t,Feynrule_t,I,J,K>::reset(2);									\

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

    /* Specialisation of the evaluate class template for vertices composed with
     * the colour-flow version of the T-colour structure: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> vertex_type;
	    
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

	    DEFINE_REP_TYPE(model_type,fundamental_rep< SU<N> >);

	public:

	    /* Function returning the index ranges of interacting subamplitude
	     * tensors: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I)
		{
		    v.insert(v.begin(),2,N);
		}
		if(n==J or n==K)
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

#define RECREL(n)									\
if(I==n)										\
{											\
    value_type f=(r_value_type)0.5*factor;						\
    white_part.reset();									\
    iterator temp=iters[I];								\
    iters[I]=white_part.begin();							\
    for(size_type c=0;c<N;++c)								\
    {											\
	evaluate<Feynrule_t>::RANK_NUMERAL(n)(f,couplings,iters,momenta);		\
	iters[J]+=sizes[0][J];								\
	iters[K]+=sizes[0][K];								\
    }											\
    iters[J]-=sizes[1][J];								\
    iters[K]-=sizes[1][K];								\
    iters[I]=temp;									\
    for(size_type a=0;a<N;++a)								\
    {											\
	for(size_type b=0;b<a;++b)							\
	{										\
	    evaluate<Feynrule_t>::RANK_NUMERAL(n)(f,couplings,iters,momenta);		\
	    iters[I]+=sizes[0][I];							\
	    iters[J]+=sizes[0][J];							\
	}										\
	for(size_type mu=0;mu<sizes[0][I];++mu)						\
	{										\
	    iters[I][mu]-=white_part[mu]/value_type(N,0);				\
	}										\
	for(size_type b=a;b<N;++b)							\
	{										\
	    evaluate<Feynrule_t>::RANK_NUMERAL(n)(f,couplings,iters,momenta);		\
	    iters[I]+=sizes[0][I];							\
	    iters[J]+=sizes[0][J];							\
	}										\
	iters[J]-=sizes[1][J];								\
	iters[K]+=sizes[0][K];								\
    }											\
    iters[I]-=sizes[2][I];								\
    iters[K]-=sizes[1][K];								\
    return;										\
}											\
else											\
{											\
    for(size_type a=0;a<N;++a)								\
    {											\
	for(size_type b=0;b<N;++b)							\
	{										\
	    evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta);	\
	    iters[I]+=sizes[0][I];							\
	    iters[K]+=sizes[0][K];							\
	}										\
	iters[K]-=sizes[1][K];								\
	iters[J]+=sizes[0][J];								\
    }											\
    iters[I]-=sizes[2][I];								\
    iters[J]-=sizes[1][J];								\
}											\

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

	    /* Memory storage tensor for the trace part of the vertex: */

	    static tensor_type white_part;

	    /* Storage array of the (sub-)tensor sizes: */

	    static const size_type sizes[3][4];

	    /* Utility index range holder vector: */

	    static std::vector<size_type> utilvec;
    };
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>std::vector<typename evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::size_type> evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::utilvec;
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>typename evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::tensor_type evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::white_part(evaluate<Feynrule_t>::get_index_ranges(I,evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::utilvec));
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::sizes[3][4]={{Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]},{N*Feynrule_t::sizes[0],N*Feynrule_t::sizes[1],N*Feynrule_t::sizes[2],N*Feynrule_t::sizes[3]},{N*N*Feynrule_t::sizes[0],N*N*Feynrule_t::sizes[1],N*N*Feynrule_t::sizes[2],N*N*Feynrule_t::sizes[3]}};

    /* Specialisation of the cfd_evaluate class template for vertices composed with
     * the colour-flow version of the T-colour structure: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> vertex_type;
	    
	    /* The usual type definition: */
	    
	    typedef typename cfd_evaluate<Feynrule_t>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_t>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_t>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_t>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_t>::momentum_type momentum_type;

	private:

	    /* Representation type definition: */

	    DEFINE_REP_TYPE(model_type,fundamental_rep< SU<N> >);

	public:

	    /* Function returning the index ranges of the interacting
	     * subamplitude tensors: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		evaluate<Feynrule_t>::get_index_ranges(n,v);
		if(n==I)
		{
		    v.insert(v.begin(),2,N);
		}
		if(n==J or n==K)
		{
		    v.insert(v.begin(),N);
		}
		return v;
	    }

	    /* Initialisation phase: */

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_t>::initialise();
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(n)															\
value_type z=(r_value_type)0.5*factor;													\
size_type a,b,c,d;															\
if(I==n)																\
{																	\
    c=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;											\
    d=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;											\
    if(c!=d)																\
    {																	\
	iters[I]+=(d*N+c)*Feynrule_t::sizes[I];												\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);						\
	iters[I]-=(d*N+c)*Feynrule_t::sizes[I];												\
	return;																\
    }																	\
    z*=(-N_inv);															\
    for(size_type i=0;i<N;++i)														\
    {																	\
	if(i==c)															\
	{																\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z+(r_value_type)0.5*factor,couplings,iters,momenta,produced_iters);		\
	}																\
	else																\
	{																\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);					\
	}																\
	iters[I]+=((N+1)*Feynrule_t::sizes[I]);												\
    }																	\
    iters[I]-=((N*N+N)*Feynrule_t::sizes[I]);												\
    return;																\
}																	\
if(J==n)																\
{																	\
    b=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;											\
    a=((iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N)-b)/N;										\
    c=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;											\
    if(b==c and a!=b)															\
    {																	\
	iters[J]+=a*Feynrule_t::sizes[J];												\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);					\
	iters[J]-=a*Feynrule_t::sizes[J];												\
	return;																\
    }																	\
    if(a==b and b!=c)															\
    {																	\
	iters[J]+=c*Feynrule_t::sizes[J];												\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-N_inv*factor,couplings,iters,momenta,produced_iters);				\
	iters[J]-=c*Feynrule_t::sizes[J];												\
	return;																\
    }																	\
    if(a==b and b==c)															\
    {																	\
	iters[J]+=a*Feynrule_t::sizes[J];												\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor-N_inv*factor,couplings,iters,momenta,produced_iters);				\
	iters[J]-=a*Feynrule_t::sizes[J];												\
	return;																\
    }																	\
    return;																\
}																	\
if(K==n)																\
{																	\
    b=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;											\
    a=((iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N)-b)/N;										\
    c=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;											\
    if(a==c and a!=b)															\
    {																	\
	iters[K]+=b*Feynrule_t::sizes[K];												\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);					\
	iters[K]-=b*Feynrule_t::sizes[K];												\
	return;																\
    }																	\
    if(a==b and b!=c)															\
    {																	\
	iters[K]+=c*Feynrule_t::sizes[K];												\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-N_inv*factor,couplings,iters,momenta,produced_iters);				\
	iters[K]-=c*Feynrule_t::sizes[K];												\
	return;																\
    }																	\
    if(a==b and b==c)															\
    {																	\
	iters[K]+=a*Feynrule_t::sizes[K];												\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor-N_inv*factor,couplings,iters,momenta,produced_iters);				\
	iters[K]-=a*Feynrule_t::sizes[K];												\
	return;																\
    }																	\
    return;																\
}																	\
else																	\
{																	\
    b=(iters[I].get_offset()/Feynrule_t::sizes[I])%N;											\
    a=((iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N)-b)/N;										\
    c=(iters[J].get_offset()/Feynrule_t::sizes[J])%N;											\
    d=(iters[K].get_offset()/Feynrule_t::sizes[K])%N;											\
    if(a==c and b==d and a!=b)														\
    {																	\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);					\
	return;																\
    }																	\
    if(a==b and c==d and (a!=c or b!=d))												\
    {																	\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-N_inv*factor,couplings,iters,momenta,produced_iters);				\
	return;																\
    }    																\
    if(a==c and b==d and a==b)														\
    {																	\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor-N_inv*factor,couplings,iters,momenta,produced_iters);				\
    }																	\
}																	\
	    
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

	    static const r_value_type N_inv;
    };
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,class Feynrule_t>const typename cfd_evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::r_value_type cfd_evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::N_inv=(typename cfd_evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::r_value_type)1/(typename cfd_evaluate<compose_vertices<colour_tensor::CF_T<N,I,J,K>,Feynrule_t> >::r_value_type)N;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_T_H_*/
