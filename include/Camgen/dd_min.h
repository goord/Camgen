//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DD_MIN_H_
#define CAMGEN_DD_MIN_H_

#include <Camgen/dd.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Anti-symmetrised double-Kronecker delta colour structure declaration and        *
 * definition. Moreover, this file contains the specialisations of the evaluate    *
 * and cfd_evaluate class templates for vertices which contain a composition       *
 * with a colour-index anti-symmetrised double Kronecker delta. There is also a    *
 * declaration and definition of the CF_dd class, which denotes the antisymmetric  *
 * double-Kronecker delta structure for adjoint representation indices in the      *
 * colour-flow representation.                                                     *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {

	/* Anti-symmetrised double-Kronecker delta class for generic representation
	 * indices: */
	
	template<class rep_t,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class dd_min
	{
	    public:

		/* Self-referring type definition: */

		typedef dd_min<rep_t,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef dd_min<rep_t,I,J,K,L> CF_type;

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
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<rep_t,I,J,K,L>::rank;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<rep_t,I,J,K,L>::tensor_size;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<rep_t,I,J,K,L>::multiplicity;
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<rep_t,I,J,K,L>::contractions[4]={I,J,K,L};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<rep_t,I,J,K,L>::ranges[4]={rep_t::dimension,rep_t::dimension,rep_t::dimension,rep_t::dimension};
	template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string dd_min<rep_t,I,J,K,L>::formula="(d(a1,a2)d(a3,a4)-d(a1,a4)d(a3,a2))";

	/* Declaration of the colour-flow anti-symmetrised double-Kronecker
	 * delta: */

	template<std::size_t N,std::size_t I=0,std::size_t J=1,std::size_t K=2,std::size_t L=3>class CF_dd_min;

	/* Specialisation in the case of the SU(N) adjoint representation: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class dd_min<adjoint_rep< SU<N> >,I,J,K,L>
	{
	    public:

		/* Self-referring type definition: */

		typedef dd_min<adjoint_rep< SU<N> >,I,J,K,L> type;
		
		/* Colour-flow mode counterpart: */
		
		typedef CF_dd_min<N,I,J,K,L> CF_type;

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
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<adjoint_rep<SU<N> >,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<adjoint_rep<SU<N> >,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<adjoint_rep<SU<N> >,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<adjoint_rep<SU<N> >,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t dd_min<adjoint_rep<SU<N> >,I,J,K,L>::ranges[4]={N*N-1,N*N-1,N*N-1,N*N-1};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string dd_min<adjoint_rep<SU<N> >,I,J,K,L>::formula="(d(a1,a2)d(a3,a4)-d(a1,a4)d(a3,a2))";

	/* Definition of the colour-flow anti-symmetrised double Kronecker
	 * delta: */

	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class CF_dd_min
	{
	    public:

		/* Self-referring type definition: */

		typedef CF_dd_min<N,I,J,K,L> type;

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
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd_min<N,I,J,K,L>::rank;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd_min<N,I,J,K,L>::tensor_size;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd_min<N,I,J,K,L>::multiplicity;
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd_min<N,I,J,K,L>::contractions[4]={I,J,K,L};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::size_t CF_dd_min<N,I,J,K,L>::ranges[4]={N*N,N*N,N*N,N*N};
	template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L>const std::string CF_dd_min<N,I,J,K,L>::formula="4(d(A1,B2)d(A2,B1)d(A3,B4)d(A4,B3)-d(A1,B4)d(A2,B3)d(A3,B2)d(A4,B1)-d(A1,B1)d(A2,B2)d(A3,B4)d(A4,B3)/N+d(A1,B1)d(A2,B3)d(A3,B2)d(A4,B4)/N-d(A1,B2)d(A2,B1)d(A3,B3)d(A4,B4)/N+d(A1,B4)d(A2,B2)d(A3,B3)d(A4,B1)/N)";
    }
    
    /* Specialisation for the evaluate class template for generic
     * anti-symmetrised double Kronecker delta compositions: */

    template<class rep_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::dd_min<rep_t,I,J,K,L>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::dd_min<rep_t,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >::initialise();
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,L,K,J>,Feynrule_t> >::initialise();
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

	    /* Recursive relations, simply applying the two double-delta terms
	     * consecutively: */

	    static void first(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >::first(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,L,K,J>,Feynrule_t> >::first(-factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >::second(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,L,K,J>,Feynrule_t> >::second(-factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >::third(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,L,K,J>,Feynrule_t> >::third(-factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,J,K,L>,Feynrule_t> >::fourth(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::dd<rep_t,I,L,K,J>,Feynrule_t> >::fourth(-factor,couplings,iters,momenta);
	    }
    };

    /* Evaluate class template specialisation for the colour-flow
     * anti-symmetrised double Kronecker delta colour structures: */
    
    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class evaluate<compose_vertices<colour_tensor::CF_dd_min<N,I,J,K,L>,Feynrule_t> >
    {
	public:
	    
	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::CF_dd_min<N,I,J,K,L>,Feynrule_t> vertex_type;
	    
	    /* Usual type definitions: */
	    
	    typedef typename evaluate<Feynrule_t>::model_type model_type;
	    typedef typename evaluate<Feynrule_t>::size_type size_type;
	    typedef typename evaluate<Feynrule_t>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_t>::value_type value_type;
	    typedef typename evaluate<Feynrule_t>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_t>::iterator iterator;
	    typedef typename evaluate<Feynrule_t>::momentum_type momentum_type;

	    /* Initialisation phase: */

	    static void initialise()
	    {
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::initialise();
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::initialise();
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

	    /* Recursive relations, simply applying the two double-delta terms
	     * consecutively: */
	    
	    static void first(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::first(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::first(-factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::second(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::second(-factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::third(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::third(-factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::fourth(factor,couplings,iters,momenta);
		evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::fourth(-factor,couplings,iters,momenta);
	    }
    };

    /* Specialisation of the cfd_evaluate class template for the
     * anti-symmetrised double Kronecker delta colour structures over
     * SU(N)-fundamental representation indices: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::dd_min<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::dd_min<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> vertex_type;
	    
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
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::initialise();
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,L,K,J>,Feynrule_t> >::initialise();
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

	    /* Recursive relations, simply applying the two double-delta terms
	     * consecutively: */
	    
	    static void first(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::first(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,L,K,J>,Feynrule_t> >::first(-factor,couplings,iters,momenta,produced_iters);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::second(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,L,K,J>,Feynrule_t> >::second(-factor,couplings,iters,momenta,produced_iters);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::third(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,L,K,J>,Feynrule_t> >::third(-factor,couplings,iters,momenta,produced_iters);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,J,K,L>,Feynrule_t> >::fourth(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::dd<fundamental_rep< SU<N> >,I,L,K,J>,Feynrule_t> >::fourth(-factor,couplings,iters,momenta,produced_iters);
	    }
    };

    /* Specialisation of the cfd_evaluate class template for anti-symmetrised
     * double Kronecker delta colour structures over SU(N)-adjoint
     * representation indices: */

    template<std::size_t N,std::size_t I,std::size_t J,std::size_t K,std::size_t L,class Feynrule_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_dd_min<N,I,J,K,L>,Feynrule_t> >
    {
	public:

	    /* Vertex class type definition: */

	    typedef compose_vertices<colour_tensor::CF_dd_min<N,I,J,K,L>,Feynrule_t> vertex_type;
	    
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
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::initialise();
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::initialise();
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
	    
	    /* Recursive relations, simply applying the two double-delta terms
	     * consecutively: */
	    
	    static void first(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::first(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::first(-factor,couplings,iters,momenta,produced_iters);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::second(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::second(-factor,couplings,iters,momenta,produced_iters);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::third(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::third(-factor,couplings,iters,momenta,produced_iters);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,J,K,L>,Feynrule_t> >::fourth(factor,couplings,iters,momenta,produced_iters);
		cfd_evaluate<compose_vertices<colour_tensor::CF_dd<N,I,L,K,J>,Feynrule_t> >::fourth(-factor,couplings,iters,momenta,produced_iters);
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_DD_PLUS_H_*/

