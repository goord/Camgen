//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/license_print.h>
#include <Camgen/stdrand.h>
#include <Camgen/rn_strm.h>
#include <Camgen/su(n).h>
#include <Camgen/adj_rep.h>
#include <Camgen/Minkowski.h>
#include <Camgen/adjoint.h>
#include <Camgen/col_flow.h>
#include <Camgen/ssss.h>
#include <Camgen/vvv.h>
#include <Camgen/vvvv.h>
#include <Camgen/f.h>
#include <Camgen/ff_contr.h>
#include <Camgen/d.h>
#include <Camgen/dd.h>
#include <Camgen/dd_plus.h>
#include <Camgen/dd_min.h>
#include <Camgen/T.h>
#include <Camgen/TT.h>
#include <Camgen/TT_min.h>
#include <Camgen/TT_plus.h>
#include <Camgen/TT_contr.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for continuous-colours versions of Camgen's colour          *
 * structures, including the gluon vertex specialisations. The testing procedure *
 * compares the results of the adjoint and colour flow versions, applied on      *
 * randomly filled subamplitudes. If both no gluons are involved, e.g. the       *
 * TT_contr type, the SU<N> specialisation is tested against a clone of the      *
 * fundamental representation.                                                   *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define VALUE_T double

namespace Camgen
{
    namespace test_utils
    {
	/* Adjoint-representation-type holder model class: */

	class testadjQCD
	{
	    public:
		const static std::size_t dimension=4;
		typedef VALUE_T value_type;
		typedef Minkowski_type spacetime_type;
		typedef adjoint colour_treatment;
		const static bool continuous_colours=true;
		const static bool continuous_helicities=true;
		const static bool coloured=true;
	};

	/* Colour-flow-type holder model class: */

	class testcfQCD
	{
	    public:
		const static std::size_t dimension=4;
		typedef VALUE_T value_type;
		typedef Minkowski_type spacetime_type;
		typedef colour_flow colour_treatment;
		const static bool continuous_colours=true;
		const static bool continuous_helicities=true;
		const static bool coloured=true;
	};

	/* Function converting a adjoint-rep. tensor to a colour-flow tensor (with spacetime
	 * sizes D) : */

	template<class rep_type>void convert_to_cf(tensor< std::complex<VALUE_T> >::iterator it0,tensor< std::complex<VALUE_T> >::const_iterator it1,std::size_t D=1)
	{
	    for(unsigned A=0;A<rep_type::index_range;++A)
	    {
		for(unsigned B=0;B<rep_type::index_range;++B)
		{
		    for(unsigned a=0;a<rep_type::group_implementation::index_range;++a)
		    {
			for(unsigned mu=0;mu<D;++mu)
			{
			    it0[mu]+=rep_type::generator(a,A,B)*it1[mu];
			}
			it1+=D;
		    }
		    it1-=(D*(rep_type::group_implementation::index_range));
		    it0+=D;
		}
	    }
	    it0-=(D*(rep_type::index_range)*(rep_type::index_range));
	}

	/* Colour structure testing class: */

	template<class rng_t,class group_t>class colour_tester
	{
	    public:

		/* Useful type definitions: */

		typedef VALUE_T r_value_type;
		typedef std::complex<r_value_type> value_type;
		typedef tensor<value_type> tensor_type;
		typedef typename tensor_type::iterator iterator;
		typedef typename tensor_type::const_iterator const_iterator;
		typedef vector<r_value_type,4> momentum_type;
		typedef typename fundamental_rep<group_t>::template implementation<VALUE_T> rep_type;
		typedef random_number_stream<r_value_type,rng_t> generator_type;

	    private:

		/* Useful constants: */

		const static std::size_t N=fundamental_rep<group_t>::dimension;
		const static std::size_t K=group_t::rank;

	    public:

		/* Constructor, assigning the iterators correctly: */

		colour_tester():adj_iters(4),cf_iters(4),c(gen(0,1),gen(0,1)),couplings(1,(value_type*)NULL),momenta(4,(momentum_type*)NULL)
	    {
		rep_type::initialise();
		for(unsigned i=0;i<4;++i)
		{
		    adj_iters[i]=adj_tens[i].begin();
		    cf_iters[i]=cf_tens[i].begin();
		    momenta[i]=&p[i];
		    p[i].fill(gen,0,1);
		}
		couplings[0]=&c;
	    }

		/* Method filling the tensor members randomly: */

		bool fill()
		{
		    for(unsigned i=0;i<4;++i)
		    {
			if(!fill(i))
			{
			    return false;
			}
		    }
		    return true;
		}

		/* Method filling adjoint-type tensor member i randomly and computing the
		 * corresponding colour flow tensor: */

		bool fill(unsigned i)
		{
		    adj_tens[i].fill(gen,value_type(1,1));
		    if(adj_trans[i])
		    {
			if(cf_tens[i].size()/(N*N)!=adj_tens[i].size()/K)
			{
			    std::cerr<<"adjoint rep. not mapped to fundamental rep. squared under colour flow transformation..."<<std::endl;
			    return false;
			}
			convert_to_cf<rep_type>(cf_iters[i],adj_iters[i],adj_tens[i].size()/K);
		    }
		    else
		    {
			if(adj_tens[i].size()!=cf_tens[i].size())
			{
			    std::cerr<<"not-adjoint rep. not invariant under colour flow transformation..."<<std::endl;
			    return false;
			}
			for(unsigned j=0;j<adj_tens[i].size();++j)
			{
			    cf_tens[i][j]=adj_tens[i][j];
			}
		    }
		    return true;
		}

		/* Clearing functions: */

		void clear()
		{
		    for(unsigned i=0;i<4;++i)
		    {
			adj_tens[i].clear();
			cf_tens[i].clear();
		    }
		}
		void clear(unsigned i)
		{
		    adj_tens[i].clear();
		    cf_tens[i].clear();
		}

		/* Resetting functions: */

		void reset()
		{
		    for(unsigned i=0;i<4;++i)
		    {
			adj_tens[i].reset();
			cf_tens[i].reset();
		    }
		}
		void reset(unsigned i)
		{
		    adj_tens[i].reset();
		    cf_tens[i].reset();
		}

		/* Function testing whether the adjoint version of coltens combined
		 * with the Feynrule spacetime vertex rule yields the same result as
		 * the colour-flow version: */

		template<class coltens,template<class model_t>class Feynrule_t>bool test()
		{
		    /* Adjoint vertex definition: */

		    typedef evaluate<compose_vertices<coltens,Feynrule_t<testadjQCD> > > adj_eval_type;
		    typedef typename adj_eval_type::size_type adj_size_type;
		    typedef typename adj_eval_type::value_type adj_value_type;
		    typedef typename adj_eval_type::momentum_type adj_momentum_type;
		    typedef typename adj_eval_type::iterator adj_iterator;
		    adj_eval_type::initialise();

		    /* Colour-flow vertex definition: */

		    typedef evaluate<compose_vertices<typename coltens::CF_type,Feynrule_t<testcfQCD> > > cf_eval_type;
		    typedef typename cf_eval_type::size_type cf_size_type;
		    typedef typename cf_eval_type::value_type cf_value_type;
		    typedef typename cf_eval_type::momentum_type cf_momentum_type;
		    typedef typename cf_eval_type::iterator cf_iterator;
		    cf_eval_type::initialise();

		    for(unsigned i=0;i<4;++i)
		    {
			adj_trans[i]=false;
		    }

		    /* Indentifying adjoint-representation subamplitudes: */

		    for(unsigned i=0;i<coltens::rank;++i)
		    {
			if(coltens::ranges[i]==K)
			{
			    adj_trans[coltens::contractions[i]]=true;
			}
		    }

		    clear();

		    /* Resizing the testing tensors: */

		    std::vector<adj_size_type>adj_indexvec;
		    std::vector<cf_size_type>cf_indexvec;
		    for(unsigned i=0;i<4;++i)
		    {
			adj_tens[i].resize(adj_eval_type::get_index_ranges((adj_size_type)i,adj_indexvec));
			cf_tens[i].resize(cf_eval_type::get_index_ranges((cf_size_type)i,cf_indexvec));
		    }

		    /* Randomly filling the testing tensors: */

		    if(!fill()){return false;}
		    value_type factor(gen(0,1),gen(0,1));
		    tensor_type check;

		    /* First recursive relation check: */

		    reset(0);
		    check=cf_tens[0];
		    adj_eval_type::first((const adj_value_type&)factor,(const std::vector<const adj_value_type*>&)couplings,(std::vector<adj_iterator>&)adj_iters,(const std::vector<const adj_momentum_type*>&)momenta);
		    cf_eval_type::first((const cf_value_type&)factor,(const std::vector<const cf_value_type*>&)couplings,(std::vector<cf_iterator>&)cf_iters,(const std::vector<const cf_momentum_type*>&)momenta);
		    cf_convert(check.begin(),0);
		    if(!equal_sequences(check,cf_tens[0]))
		    {
			std::cerr<<"Mismatch adjoint<->colour flow vertices in 1st recursive relation."<<std::endl;
			return false;
		    }

		    /* Second recursive relation check: */

		    reset(1);
		    check=cf_tens[1];
		    adj_eval_type::second((const adj_value_type&)factor,(const std::vector<const adj_value_type*>&)couplings,(std::vector<adj_iterator>&)adj_iters,(const std::vector<const adj_momentum_type*>&)momenta);
		    cf_eval_type::second((const cf_value_type&)factor,(const std::vector<const cf_value_type*>&)couplings,(std::vector<cf_iterator>&)cf_iters,(const std::vector<const cf_momentum_type*>&)momenta);
		    cf_convert(check.begin(),1);
		    if(!equal_sequences(check,cf_tens[1]))
		    {
			std::cerr<<"Mismatch adjoint<->colour flow vertices in 2nd recursive relation."<<std::endl;
			return false;
		    }

		    /* Third recursive relation check: */

		    reset(2);
		    check=cf_tens[2];
		    adj_eval_type::third((const adj_value_type&)factor,(const std::vector<const adj_value_type*>&)couplings,(std::vector<adj_iterator>&)adj_iters,(const std::vector<const adj_momentum_type*>&)momenta);
		    cf_eval_type::third((const cf_value_type&)factor,(const std::vector<const cf_value_type*>&)couplings,(std::vector<cf_iterator>&)cf_iters,(const std::vector<const cf_momentum_type*>&)momenta);
		    cf_convert(check.begin(),2);
		    if(!equal_sequences(check,cf_tens[2]))
		    {
			std::cerr<<"Mismatch adjoint<->colour flow vertices in 3rd recursive relation."<<std::endl;
			return false;
		    }

		    /* Fourth recursive relation check: */

		    reset(3);
		    check=cf_tens[3];
		    adj_eval_type::fourth((const adj_value_type&)factor,(const std::vector<const adj_value_type*>&)couplings,(std::vector<adj_iterator>&)adj_iters,(const std::vector<const adj_momentum_type*>&)momenta);
		    cf_eval_type::fourth((const cf_value_type&)factor,(const std::vector<const cf_value_type*>&)couplings,(std::vector<cf_iterator>&)cf_iters,(const std::vector<const cf_momentum_type*>&)momenta);
		    cf_convert(check.begin(),3);
		    if(!equal_sequences(check,cf_tens[3]))
		    {
			std::cerr<<"Mismatch adjoint<->colour flow vertices in 4th recursive relation."<<std::endl;
			return false;
		    }
		    return true;
		}

		/* Function testing whether coltens1 and coltens2 combined with the
		 * spacetime vertex rule Feynrule_t yield the same result: */

		template<class coltens1,class coltens2,template<class model_t>class Feynrule_t>bool comp_test()
		{
		    typedef evaluate<compose_vertices<coltens1,Feynrule_t<testadjQCD> > > eval_type;
		    typedef typename eval_type::size_type size_type;
		    typedef typename eval_type::value_type value_type;
		    typedef typename eval_type::momentum_type adj_momentum_type;
		    typedef typename eval_type::iterator iterator;
		    eval_type::initialise();

		    typedef evaluate<compose_vertices<coltens2,Feynrule_t<testadjQCD> > > check_eval_type;
		    check_eval_type::initialise();

		    for(unsigned i=0;i<4;++i)
		    {
			adj_trans[i]=false;
		    }

		    clear();
		    std::vector<size_type>indexvec;
		    for(unsigned i=0;i<4;++i)
		    {
			adj_tens[i].resize(eval_type::get_index_ranges((size_type)i,indexvec));
			cf_tens[i]=adj_tens[i];
		    }
		    fill();
		    value_type factor(gen(0,1),gen(0,1));

		    reset(0);
		    eval_type::first(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)adj_iters,(const std::vector<const momentum_type*>&)momenta);
		    check_eval_type::first(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)cf_iters,(const std::vector<const momentum_type*>&)momenta);
		    if(!equal_sequences(adj_tens[0],cf_tens[0]))
		    {
			std::cerr<<"Incorrectly specialised 1st recursive relation."<<std::endl;
			return false;
		    }
		    reset(1);
		    eval_type::second(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)adj_iters,(const std::vector<const momentum_type*>&)momenta);
		    check_eval_type::second(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)cf_iters,(const std::vector<const momentum_type*>&)momenta);
		    if(!equal_sequences(adj_tens[1],cf_tens[1]))
		    {
			std::cerr<<"Incorrectly specialised 2nd recursive relation."<<std::endl;
			return false;
		    }
		    reset(2);
		    eval_type::third(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)adj_iters,(const std::vector<const momentum_type*>&)momenta);
		    check_eval_type::third(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)cf_iters,(const std::vector<const momentum_type*>&)momenta);
		    if(!equal_sequences(adj_tens[2],cf_tens[2]))
		    {
			std::cerr<<"Incorrectly specialised 3rd recursive relation."<<std::endl;
			return false;
		    }
		    reset(3);
		    eval_type::fourth(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)adj_iters,(const std::vector<const momentum_type*>&)momenta);
		    check_eval_type::fourth(factor,(const std::vector<const value_type*>&)couplings,(std::vector<iterator>&)cf_iters,(const std::vector<const momentum_type*>&)momenta);
		    if(!equal_sequences(adj_tens[3],cf_tens[3]))
		    {
			std::cerr<<"Incorrectly specialised 4th recursive relation."<<std::endl;
			return false;
		    }
		    return true;
		}

		/* Colour-flow converter: */

		void cf_convert(iterator it,unsigned i) const
		{
		    if(adj_trans[i])
		    {
			convert_to_cf<rep_type>(it,adj_iters[i],adj_tens[i].size()/K);
		    }
		    else
		    {
			for(unsigned j=0;j<adj_tens[i].size();++j)
			{
			    it[j]+=adj_tens[i][j];
			}
		    }
		}

	    private:
		/* Random number generator: */

		generator_type gen;

		/* Adjoint-representation tensors: */ 

		tensor_type adj_tens[4];

		/* Colour-flow representation tensors: */

		tensor_type cf_tens[4];

		/* Iterators to the adj. rep. tensors: */

		std::vector<iterator>adj_iters;

		/* Iterators to the col. flow rep. tensors: */

		std::vector<iterator>cf_iters;

		/* Coupling: */

		value_type c;
		std::vector<value_type*>couplings;

		/* Momenta: */

		momentum_type p[4];
		std::vector<momentum_type*>momenta;

		/* Bools denoting the adjoint-representation tensors: */

		bool adj_trans[4];
	};

	/* Copy of the SU(N)-fundamental representation: */

	template<class G>class fundam_rep_clone{};

	template<std::size_t N>class fundam_rep_clone< SU<N> >
	{
	    public:

		/* Group wrapper class type definion: */

		typedef SU<N> group_type;

		/* Dimension of the representation: */

		static const std::size_t dimension=N;

		/* Implementation class template: */

		template<class value_t>class implementation: public representation<value_t,implementation<value_t>,N,N*N-1>
	    {
		public:

		    /* Type reference to the base class template: */

		    typedef representation<value_t,implementation<value_t>,N,N*N-1> base_type;

		    /* Group implementation class type definition: */

		    typedef typename SU<N>::template implementation<value_t> group_implementation;

		    /* Hermiticity tag set true: */

		    static const bool hermitian=true;

		    /* Instead of providing a fill_generators() function, the
		     * initialisation function is overloaded (omitting the
		     * checks usually performed): */

		    static void initialise()
		    {
			if(!initialised)
			{
			    fundamental_rep< SU<N> >::template implementation<value_t>::initialise();
			    for(std::size_t a=0;a<N*N-1;++a)
			    {
				for(std::size_t i=0;i<N;++i)
				{
				    for(std::size_t j=0;j<N;++j)
				    {
					base_type::T[a][i][j]=fundamental_rep< SU<N> >::template implementation<value_t>::generator(a,i,j);
				    }
				}
			    }
			    initialised=true;
			}
		    }
		private:
		    static bool initialised;
	    };
	};
	template<std::size_t N>const std::size_t fundam_rep_clone< SU<N> >::dimension;
	template<std::size_t N>template<class value_t>const bool fundam_rep_clone< SU<N> >::implementation<value_t>::hermitian;
	template<std::size_t N>template<class value_t>bool fundam_rep_clone< SU<N> >::implementation<value_t>::initialised=false;
    }
}

/* Testing program: */

#define N_COLS 3

#define GAUGE_GROUP SU<N_COLS>

using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing continuous colour structures......................"<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    test_utils::colour_tester<std::random,GAUGE_GROUP>coltest;
    std::cerr<<"Checking continuous d<I,J> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::d<adjoint_rep<GAUGE_GROUP>,1,2>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;
    
    std::cerr<<"Checking continuous dd<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::dd<adjoint_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;
    
    std::cerr<<"Checking continuous dd_plus<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::dd_plus<adjoint_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;
    
    std::cerr<<"Checking continuous dd_min<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::dd_min<adjoint_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous T<I,J,K> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::T<fundamental_rep<GAUGE_GROUP>,0,1,2>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous TT<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::TT<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous TT_plus<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::TT_plus<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous TT_min<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::TT_min<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous TT_contr<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.comp_test< colour_tensor::TT_contr<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,colour_tensor::TT_contr<test_utils::fundam_rep_clone<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous f<I,J,K> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::f<GAUGE_GROUP,0,1,2>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous ff_contr<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::ff_contr<GAUGE_GROUP,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous ggg vertex..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::f<GAUGE_GROUP,0,1,2>,vvv >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking continuous gggg vertex..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::ff_contr<GAUGE_GROUP,0,1,2>,vvvv >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;
}

