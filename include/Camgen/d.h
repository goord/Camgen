//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_D_H_
#define CAMGEN_D_H_

#include <vector>
#include <Camgen/ppranknum.h>
#include <Camgen/comp_vert.h>
#include <Camgen/adj_rep.h>
#include <Camgen/su(n).h>
#include <Camgen/eval.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Kronecker delta colour structure declaration and definition. Moreover, this *
 * file contains the specialisations of the evaluate and cfd_evaluate class    *
 * templates for vertices which contain a composition with a colour-index      *
 * Kronecker delta. There is also a declaration and definition of the CF_d     *
 * class, which denotes the Kronecker delta for adjoint-representation indices *
 * in the colour-flow representation.                                          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Contains Camgen's color structures.

    namespace colour_tensor
    {

	/* Kronecker delta class for generic representation indices: */
	
	template<class rep_t,std::size_t I=0,std::size_t J=1>class d
	{
	    public:

		/* Self-referring type definition: */

		typedef d<rep_t,I,J> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef d<rep_t,I,J> CF_type;

		/* Vertex rank: */

		static const std::size_t rank=2;
		
		/* Total tensor size: */
		
		static const std::size_t tensor_size=(rep_t::dimension)*(rep_t::dimension);
		
		/* Contracted indices: */
		
		static const std::size_t contractions[2];
		
		/* Ranges of contracted indices: */
		
		static const std::size_t ranges[2];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<class rep_t,std::size_t I,std::size_t J>const std::size_t d<rep_t,I,J>::rank;
	template<class rep_t,std::size_t I,std::size_t J>const std::size_t d<rep_t,I,J>::tensor_size;
	template<class rep_t,std::size_t I,std::size_t J>const std::size_t d<rep_t,I,J>::multiplicity;
	template<class rep_t,std::size_t I,std::size_t J>const std::size_t d<rep_t,I,J>::contractions[2]={I,J};
	template<class rep_t,std::size_t I,std::size_t J>const std::size_t d<rep_t,I,J>::ranges[2]={rep_t::dimension,rep_t::dimension};
	template<class rep_t,std::size_t I,std::size_t J>const std::string d<rep_t,I,J>::formula="d(a1,a2)";

	/* Declaration of the colour-flow Kronecker delta: */

	template<std::size_t N,std::size_t I=0,std::size_t J=1>class CF_d;

	/* Specialisation in the case of the SU(N) adjoint representation: */

	template<std::size_t N,std::size_t I,std::size_t J>class d<adjoint_rep< SU<N> >,I,J>
	{
	    public:

		/* Self-referring type definition: */

		typedef d<adjoint_rep< SU<N> >,I,J> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef CF_d<N,I,J> CF_type;

		/* Vertex rank: */

		static const std::size_t rank=2;
		
		/* Total tensor size: */
		
		static const std::size_t tensor_size=(N*N-1)*(N*N-1);
		
		/* Contracted indices: */
		
		static const std::size_t contractions[2];
		
		/* Ranges of contracted indices: */
		
		static const std::size_t ranges[2];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t d<adjoint_rep< SU<N> >,I,J>::rank;
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t d<adjoint_rep< SU<N> >,I,J>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t d<adjoint_rep< SU<N> >,I,J>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t d<adjoint_rep< SU<N> >,I,J>::contractions[2]={I,J};
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t d<adjoint_rep< SU<N> >,I,J>::ranges[2]={N*N-1,N*N-1};
	template<std::size_t N,std::size_t I,std::size_t J>const std::string d<adjoint_rep< SU<N> >,I,J>::formula="d(a1,a2)";

	/* Definition of the colour-flow Kronecker delta: */
	
	template<std::size_t N,std::size_t I,std::size_t J>class CF_d
	{
	    public:

		/* Self-referring type definition: */

		typedef CF_d<N,I,J> type;

		/* Vertex rank: */

		static const std::size_t rank=2;
		
		/* Total tensor size: */
		
		static const std::size_t tensor_size=N*N*N*N;
		
		/* Contracted indices: */
		
		static const std::size_t contractions[2];
		
		/* Ranges of contracted indices: */
		
		static const std::size_t ranges[2];

		/* Output to model logfile: */

		static const std::string formula;

		/* Utility integer for composition: */

		static const std::size_t multiplicity=1;
	};
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t CF_d<N,I,J>::rank;
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t CF_d<N,I,J>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t CF_d<N,I,J>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t CF_d<N,I,J>::contractions[2]={I,J};
	template<std::size_t N,std::size_t I,std::size_t J>const std::size_t CF_d<N,I,J>::ranges[2]={N*N,N*N};
	template<std::size_t N,std::size_t I,std::size_t J>const std::string CF_d<N,I,J>::formula="2(d(A1,B2)d(A2,B1)-d(A1,B1)d(A2,B2)/N)";
    }
    
    /* Specialisation for the evaluate class template for generic Kronecker delta
     * compositions: */

    template<class rep_t,std::size_t I,std::size_t J,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> vertex_type;
	    
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
		if(n==I or n==J)
		{
		    v.insert(v.begin(),rep_t::dimension);
		}
		return v;
	    }

	    /* Convenient preprocessor macro: */

#define RECREL(i)									\
	    for(size_type a=0;a<rep_t::dimension;++a)					\
	    {										\
		evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta);	\
		iters[I]+=Isize0;							\
		iters[J]+=Jsize0;							\
	    }										\
	    iters[I]-=Isize1;								\
	    iters[J]-=Jsize1;								\

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
	
	private:

	    /* Static subtensor sizes: */

	    static const size_type Isize0;
	    static const size_type Jsize0;
	    static const size_type Isize1;
	    static const size_type Jsize1;
    };

    template<class rep_t,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::Isize0=Feynrule_t::sizes[I];
    template<class rep_t,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::Jsize0=Feynrule_t::sizes[J];
    template<class rep_t,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::Isize1=rep_t::dimension*Feynrule_t::sizes[I];
    template<class rep_t,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::d<rep_t,I,J>,Feynrule_t> >::Jsize1=rep_t::dimension*Feynrule_t::sizes[J];

    /* Evaluate class template specialisation for colour-flow Kronecker delta
     * colour structures: */

    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> vertex_type;
	    
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
		if(n==I or n==J)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(i)											\
if(I==i or J==i)											\
{													\
    for(size_type a=0;a<N*N;++a)									\
    {													\
	evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta);				\
	iters[I]+=Isize0;										\
	iters[J]+=Jsize0;										\
    }													\
    iters[I]-=Isize2;											\
    iters[J]-=Jsize2;											\
    return;												\
}													\
for(size_type a=0;a<N;++a)										\
{													\
    for(size_type b=0;b<N;++b)										\
    {													\
	evaluate<Feynrule_t>::RANK_NUMERAL(i)((r_value_type)2*factor,couplings,iters,momenta);		\
	iters[I]+=Isize0;										\
	iters[J]+=Jsize1;										\
    }													\
    iters[J]+=(Jsize0-Jsize2);										\
}													\
iters[I]-=Isize2;											\
iters[J]-=Jsize1;											\

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
	
	private:

	    /* Static subtensor sizes: */

	    static const size_type Isize0;
	    static const size_type Isize1;
	    static const size_type Isize2;
	    static const size_type Jsize0;
	    static const size_type Jsize1;
	    static const size_type Jsize2;
    };

    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::Isize0=Feynrule_t::sizes[I];
    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::Jsize0=Feynrule_t::sizes[J];
    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::Isize1=N*Feynrule_t::sizes[I];
    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::Jsize1=N*Feynrule_t::sizes[J];
    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::Isize2=N*N*Feynrule_t::sizes[I];
    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>const typename evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::size_type evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >::Jsize2=N*N*Feynrule_t::sizes[J];

    /* Specialisation of the cfd_evaluate class template for Kronecker delta
     * colour structures over SU(N)-fundamental representation indices: */

    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::d<fundamental_rep< SU<N> >,I,J>,Feynrule_t> >
    {
	public:

	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::d<fundamental_rep< SU<N> >,I,J>,Feynrule_t> vertex_type;
	    
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
		if(n==I or n==J)
		{
		    v.insert(v.begin(),N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro definition: */

#define	RECREL(i)											\
if(I==i)												\
{													\
    iters[I]+=((iters[J].get_offset()/Feynrule_t::sizes[J])%N)*Feynrule_t::sizes[I];			\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta,produced_iters);		\
    iters[I]-=((iters[J].get_offset()/Feynrule_t::sizes[J])%N)*Feynrule_t::sizes[I];			\
    return;												\
}													\
if(J==i)												\
{													\
    iters[J]+=((iters[I].get_offset()/Feynrule_t::sizes[I])%N)*Feynrule_t::sizes[J];			\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta,produced_iters);		\
    iters[J]-=((iters[I].get_offset()/Feynrule_t::sizes[I])%N)*Feynrule_t::sizes[J];			\
    return;												\
}													\
if((iters[I].get_offset()/Feynrule_t::sizes[I])%N==(iters[J].get_offset()/Feynrule_t::sizes[J])%N)	\
{													\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(i)(factor,couplings,iters,momenta,produced_iters);		\
}													\
return;													\

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

    /* Specialisation of the cfd_evaluate class template for Kronecker delta
     * colour structures over SU(N)-adjoint representation indices: */

    template<std::size_t N,std::size_t I,std::size_t J,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> >
    {
	public:

	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::CF_d<N,I,J>,Feynrule_t> vertex_type;
	    
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
		if(n==I or n==J)
		{
		    v.insert(v.begin(),2,N);
		}
		return v;
	    }

	    /* Convenient preprocessor macro definition: */

#define RECREL(n)												\
if(I==n)													\
{														\
    size_type l=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);						\
    size_type a=l/N;												\
    size_type b=l%N;												\
    if(a==b)													\
    {														\
	value_type z=factor/((r_value_type)N);									\
	size_type step=(N+1)*Feynrule_t::sizes[I];								\
	for(size_type i=0;i<a;++i)										\
	{													\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);		\
	    iters[I]+=step;											\
	}													\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor-z,couplings,iters,momenta,produced_iters);		\
	iters[I]+=step;												\
	for(size_type i=a+1;i<N;++i)										\
	{													\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);		\
	    iters[I]+=step;											\
	}													\
	iters[I]-=(N*step);											\
	return;													\
    }														\
    iters[I]+=(l*Feynrule_t::sizes[I]);										\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);			\
    iters[I]-=(l*Feynrule_t::sizes[I]);										\
    return;													\
}														\
if(J==n)													\
{														\
    size_type l=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);						\
    size_type a=l/N;												\
    size_type b=l%N;												\
    if(a==b)													\
    {														\
	value_type z=factor/((r_value_type)N);									\
	size_type step=(N+1)*Feynrule_t::sizes[J];								\
	for(size_type i=0;i<a;++i)										\
	{													\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);		\
	    iters[J]+=step;											\
	}													\
	cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor-z,couplings,iters,momenta,produced_iters);		\
	iters[J]+=step;												\
	for(size_type i=a+1;i<N;++i)										\
	{													\
	    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z,couplings,iters,momenta,produced_iters);		\
	    iters[J]+=step;											\
	}													\
	iters[J]-=(N*step);											\
	return;													\
    }														\
    iters[J]+=(l*Feynrule_t::sizes[J]);										\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(factor,couplings,iters,momenta,produced_iters);			\
    iters[J]-=(l*Feynrule_t::sizes[J]);										\
    return;													\
}														\
size_type lI=(iters[I].get_offset()/Feynrule_t::sizes[I])%(N*N);						\
size_type aI=lI/N;												\
size_type bI=lI%N;												\
size_type lJ=(iters[J].get_offset()/Feynrule_t::sizes[J])%(N*N);						\
size_type aJ=lJ/N;												\
size_type bJ=lJ%N;												\
value_type z=(r_value_type)2*factor;										\
if(aI==bJ and bI==aJ and aI!=aJ)										\
{														\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta,produced_iters);			\
    return;													\
}														\
if(aI==bI and aJ==bJ and aI!=bJ)										\
{														\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(-z/(r_value_type)N,couplings,iters,momenta,produced_iters);	\
    return;													\
}														\
if(aI==bJ and bI==aJ and aI==bI)										\
{														\
    cfd_evaluate<Feynrule_t>::RANK_NUMERAL(n)(z-z/(r_value_type)N,couplings,iters,momenta,produced_iters);	\
    return;													\
}														\
return

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
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_D_H_*/

