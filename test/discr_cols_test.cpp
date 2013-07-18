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
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for discrete-colours versions of Camgen's colour        *
 * structures, including the gluon 4-vertex specialisations. The testing     *
 * procedure compares the results of the continuous and discrete colour flow *
 * versions, applied on randomly filled subamplitudes.                       *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define VALUE_T double

namespace Camgen
{
    namespace test_utils
    {
	/* Continuous colour-flow model class: */

	class testQCD
	{
	    public:
		const static std::size_t dimension=4;
		typedef VALUE_T value_type;
		typedef Minkowski_type spacetime_type;
		typedef colour_flow colour_treatment;
		const static bool coloured=true;
		const static bool continuous_helicities=true;
		const static bool continuous_colours=false;
	};

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

		colour_tester():cont_iters(4),discr_iters(4),c(gen(0,1),gen(0,1)),couplings(1,(value_type*)NULL),momenta(4,(momentum_type*)NULL)
	    {
		rep_type::initialise();
		for(unsigned i=0;i<4;++i)
		{
		    cont_iters[i]=cont_tens[i].begin();
		    discr_iters[i]=discr_tens[i].begin();
		    momenta[i]=&p[i];
		    p[i].fill(gen,0,1);
		    adj_reps[i]=false;
		    fund_reps[i]=false;
		}
		couplings[0]=&c;
	    }

		/* Method filling spacetime subtensor member i : */

		void fill(unsigned i)
		{
		    if(adj_reps[i])
		    {
			if((discr_iters[i].get_offset()/st_sizes[i])%N==(discr_iters[i].get_offset()/st_sizes[i])/N)
			{
			    for(unsigned mu=0;mu<st_sizes[i];++mu)
			    {
				discr_iters[i][mu]=value_type(gen(0,1),gen(0,1));
			    }
			    for(unsigned aa=0;aa<N;++aa)
			    {
				if(aa==(discr_iters[i].get_offset()/st_sizes[i])%N)
				{
				    for(unsigned mu=0;mu<st_sizes[i];++mu)
				    {
					cont_iters[i][mu]=r_value_type(N-1)*discr_iters[i][mu]/(r_value_type)N;
				    }
				}
				else
				{
				    for(unsigned mu=0;mu<st_sizes[i];++mu)
				    {
					cont_iters[i][mu]=-discr_iters[i][mu]/(r_value_type)N;
				    }
				}
				cont_iters[i]+=((N+1)*st_sizes[i]);
			    }
			    cont_iters[i]=cont_tens[i].begin();
			    return;
			}
			else
			{
			    cont_iters[i]=cont_tens[i].begin()+discr_iters[i].get_offset();
			    for(unsigned mu=0;mu<st_sizes[i];++mu)
			    {
				discr_iters[i][mu]=value_type(gen(0,1),gen(0,1));
				cont_iters[i][mu]=discr_iters[i][mu];
			    }
			    cont_iters[i]=cont_tens[i].begin();
			    return;
			}
		    }
		    else
		    {
			cont_iters[i]=cont_tens[i].begin()+discr_iters[i].get_offset();
			for(unsigned mu=0;mu<st_sizes[i];++mu)
			{
			    discr_iters[i][mu]=value_type(gen(0,1),gen(0,1));
			    cont_iters[i][mu]=discr_iters[i][mu];
			}
			cont_iters[i]=cont_tens[i].begin();
			return;
		    }
		}

		/* Clearing functions: */

		void clear()
		{
		    for(unsigned i=0;i<4;++i)
		    {
			clear(i);
		    }
		}
		void clear(unsigned i)
		{
		    cont_tens[i].clear();
		    discr_tens[i].clear();
		    discr_iters[i]=discr_tens[i].begin();
		}

		/* Resetting functions: */

		void reset()
		{
		    for(unsigned i=0;i<4;++i)
		    {
			reset(i);
		    }
		}
		void reset(unsigned i)
		{
		    cont_tens[i].reset();
		    for(unsigned mu=0;mu<st_sizes[i];++mu)
		    {
			discr_iters[i][mu]=value_type(0,0);
		    }
		}

		void print_cols(std::ostream& os) const
		{
		    for(unsigned i=0;i<4;++i)
		    {
			if(adj_reps[i])
			{
			    os<<"("<<discr_iters[i].get_offset()/st_sizes[i]/N<<","<<(discr_iters[i].get_offset()/st_sizes[i])%N<<")";
			}
			else if(fund_reps[i])
			{
			    os<<"("<<discr_iters[i].get_offset()/st_sizes[i]<<")";
			}
			else
			{
			    os<<"()";
			}
		    }
		    os<<std::endl;
		}

		/* Function testing whether the continuous version of coltens combined
		 * with the Feynrule spacetime vertex rule yields the same result as
		 * the discrete colour version: */

		template<class coltens,template<class model_t>class Feynrule_t>bool test()
		{
		    typedef compose_vertices<typename coltens::CF_type,Feynrule_t<testQCD> > vertex_type;
		    /* Continuous-colours vertex definition: */

		    typedef evaluate<vertex_type> cont_eval_type;

		    /* Discrete-colours vertex definition: */

		    typedef cfd_evaluate<vertex_type> discr_eval_type;

		    /* The usual type definitions: */

		    DEFINE_BASIC_TYPES(testQCD);

		    /* Initialisation: */

		    cont_eval_type::initialise();
		    discr_eval_type::initialise();
		    for(unsigned i=0;i<4;++i)
		    {
			st_sizes[i]=Feynrule_t<testQCD>::sizes[i];
		    }
		    clear();

		    /* Indentifying adjoint and fundamental representation subamplitudes: */

		    for(unsigned i=0;i<4;++i)
		    {
			adj_reps[i]=false;
			fund_reps[i]=false;

		    }
		    for(unsigned i=0;i<coltens::rank;++i)
		    {
			if(coltens::ranges[i]==K)
			{
			    adj_reps[coltens::contractions[i]]=true;
			}
			if(coltens::ranges[i]==N)
			{
			    fund_reps[coltens::contractions[i]]=true;
			}
		    }

		    if(vertex_type::rank==3)
		    {
			/* Resizing the testing tensors: */

			std::vector<size_type>indexvec;
			for(unsigned i=0;i<3;++i)
			{
			    cont_tens[i].resize(cont_eval_type::get_index_ranges(i,indexvec));
			    discr_tens[i].resize(indexvec);
			}

			value_type factor(generator_type::throw_number(),generator_type::throw_number());

			/* First recursive relation check: */

			std::set<typename discr_eval_type::iterator>iterset;
			for(;discr_iters[1]!=discr_tens[1].end();discr_iters[1]+=st_sizes[1])
			{
			    fill(1);
			    for(;discr_iters[2]!=discr_tens[2].end();discr_iters[2]+=st_sizes[2])
			    {
				fill(2);
				cont_eval_type::first((const typename cont_eval_type::value_type)factor,(const std::vector<const typename cont_eval_type::value_type*>&)couplings,(std::vector<typename cont_eval_type::iterator>&)cont_iters,(const std::vector<const typename cont_eval_type::momentum_type*>&)momenta);
				discr_eval_type::first((const typename discr_eval_type::value_type)factor,(const std::vector<const typename discr_eval_type::value_type*>&)couplings,(std::vector<typename discr_eval_type::iterator>&)discr_iters,(const std::vector<const typename discr_eval_type::momentum_type*>&)momenta,iterset);
				if(!equal_sequences(cont_tens[0],discr_tens[0]))
				{
				    std::cerr<<"Mismatch continuous<->discrete colour flow vertices in 1st recursive relation."<<std::endl;
				    std::cerr<<"Failed colour configuration: ";
				    print_cols(std::cerr);
				    return false;
				}
				cont_tens[0].reset();
				discr_tens[0].reset();
				reset(2);
			    }
			    discr_iters[2]=discr_tens[2].begin();
			    reset(1);
			}
			discr_iters[1]=discr_tens[1].begin();

			/* Second recursive relation check: */

			for(;discr_iters[0]!=discr_tens[0].end();discr_iters[0]+=st_sizes[0])
			{
			    fill(0);
			    for(;discr_iters[2]!=discr_tens[2].end();discr_iters[2]+=st_sizes[2])
			    {
				fill(2);
				cont_eval_type::second((const typename cont_eval_type::value_type)factor,(const std::vector<const typename cont_eval_type::value_type*>&)couplings,(std::vector<typename cont_eval_type::iterator>&)cont_iters,(const std::vector<const typename cont_eval_type::momentum_type*>&)momenta);
				discr_eval_type::second((const typename discr_eval_type::value_type)factor,(const std::vector<const typename discr_eval_type::value_type*>&)couplings,(std::vector<typename discr_eval_type::iterator>&)discr_iters,(const std::vector<const typename discr_eval_type::momentum_type*>&)momenta,iterset);
				if(!equal_sequences(cont_tens[1],discr_tens[1]))
				{
				    std::cerr<<"Mismatch continuous<->discrete colour flow vertices in 2nd recursive relation."<<std::endl;
				    std::cerr<<"Failed colour configuration: ";
				    print_cols(std::cerr);
				    return false;
				}
				cont_tens[1].reset();
				discr_tens[1].reset();
				reset(2);
			    }
			    discr_iters[2]=discr_tens[2].begin();
			    reset(0);
			}
			discr_iters[0]=discr_tens[0].begin();

			/* Third recursive relation check: */

			for(;discr_iters[0]!=discr_tens[0].end();discr_iters[0]+=st_sizes[0])
			{
			    fill(0);
			    for(;discr_iters[1]!=discr_tens[1].end();discr_iters[1]+=st_sizes[1])
			    {
				fill(1);
				cont_eval_type::third((const typename cont_eval_type::value_type)factor,(const std::vector<const typename cont_eval_type::value_type*>&)couplings,(std::vector<typename cont_eval_type::iterator>&)cont_iters,(const std::vector<const typename cont_eval_type::momentum_type*>&)momenta);
				discr_eval_type::third((const typename discr_eval_type::value_type)factor,(const std::vector<const typename discr_eval_type::value_type*>&)couplings,(std::vector<typename discr_eval_type::iterator>&)discr_iters,(const std::vector<const typename discr_eval_type::momentum_type*>&)momenta,iterset);
				if(!equal_sequences(cont_tens[2],discr_tens[2]))
				{
				    std::cerr<<"Mismatch continuous<->discrete colour flow vertices in 3rd recursive relation."<<std::endl;
				    std::cerr<<"Failed colour configuration: ";
				    print_cols(std::cerr);
				    std::cerr<<cont_tens[2]<<","<<discr_tens[2]<<std::endl;
				    return false;
				}
				cont_tens[2].reset();
				discr_tens[2].reset();
				reset(1);
			    }
			    discr_iters[1]=discr_tens[1].begin();
			    reset(0);
			}
			discr_iters[0]=discr_tens[0].begin();
		    }
		    if(vertex_type::rank==4)
		    {
			/* Resizing the testing tensors: */

			std::vector<size_type>indexvec;
			for(unsigned i=0;i<4;++i)
			{
			    cont_tens[i].resize(cont_eval_type::get_index_ranges(i,indexvec));
			    discr_tens[i].resize(indexvec);
			}

			value_type factor(generator_type::throw_number(),generator_type::throw_number());

			/* First recursive relation check: */

			std::set<typename discr_eval_type::iterator>iterset;
			for(;discr_iters[1]!=discr_tens[1].end();discr_iters[1]+=st_sizes[1])
			{
			    fill(1);
			    for(;discr_iters[2]!=discr_tens[2].end();discr_iters[2]+=st_sizes[2])
			    {
				fill(2);
				for(;discr_iters[3]!=discr_tens[3].end();discr_iters[3]+=st_sizes[3])
				{
				    fill(3);
				    cont_eval_type::first((const typename cont_eval_type::value_type)factor,(const std::vector<const typename cont_eval_type::value_type*>&)couplings,(std::vector<typename cont_eval_type::iterator>&)cont_iters,(const std::vector<const typename cont_eval_type::momentum_type*>&)momenta);
				    discr_eval_type::first((const typename discr_eval_type::value_type)factor,(const std::vector<const typename discr_eval_type::value_type*>&)couplings,(std::vector<typename discr_eval_type::iterator>&)discr_iters,(const std::vector<const typename discr_eval_type::momentum_type*>&)momenta,iterset);
				    if(!equal_sequences(cont_tens[0],discr_tens[0]))
				    {
					std::cerr<<"Mismatch continuous<->discrete colour flow vertices in 1st recursive relation."<<std::endl;
					std::cerr<<"Failed colour configuration: ";
					print_cols(std::cerr);
					std::cerr<<cont_tens[0]<<","<<discr_tens[0]<<std::endl;
					return false;
				    }
				    cont_tens[0].reset();
				    discr_tens[0].reset();
				    reset(3);
				}
				discr_iters[3]=discr_tens[3].begin();
				reset(2);
			    }
			    discr_iters[2]=discr_tens[2].begin();
			    reset(1);
			}
			discr_iters[1]=discr_tens[1].begin();

			/* Second recursive relation check: */

			for(;discr_iters[0]!=discr_tens[0].end();discr_iters[0]+=st_sizes[0])
			{
			    fill(0);
			    for(;discr_iters[2]!=discr_tens[2].end();discr_iters[2]+=st_sizes[2])
			    {
				fill(2);
				for(;discr_iters[3]!=discr_tens[3].end();discr_iters[3]+=st_sizes[3])
				{
				    fill(3);
				    cont_eval_type::second((const typename cont_eval_type::value_type)factor,(const std::vector<const typename cont_eval_type::value_type*>&)couplings,(std::vector<typename cont_eval_type::iterator>&)cont_iters,(const std::vector<const typename cont_eval_type::momentum_type*>&)momenta);
				    discr_eval_type::second((const typename discr_eval_type::value_type)factor,(const std::vector<const typename discr_eval_type::value_type*>&)couplings,(std::vector<typename discr_eval_type::iterator>&)discr_iters,(const std::vector<const typename discr_eval_type::momentum_type*>&)momenta,iterset);
				    if(!equal_sequences(cont_tens[1],discr_tens[1]))
				    {
					std::cerr<<"Mismatch continuous<->discrete colour flow vertices in 2nd recursive relation."<<std::endl;
					std::cerr<<"Failed colour configuration: ";
					print_cols(std::cerr);
					std::cerr<<cont_tens[1]<<","<<discr_tens[1]<<std::endl;
					return false;
				    }
				    cont_tens[1].reset();
				    discr_tens[1].reset();
				    reset(3);
				}
				discr_iters[3]=discr_tens[3].begin();
				reset(2);
			    }
			    discr_iters[2]=discr_tens[2].begin();
			    reset(0);
			}
			discr_iters[0]=discr_tens[0].begin();

			/* Third recursive relation check: */

			for(;discr_iters[0]!=discr_tens[0].end();discr_iters[0]+=st_sizes[0])
			{
			    fill(0);
			    for(;discr_iters[1]!=discr_tens[1].end();discr_iters[1]+=st_sizes[1])
			    {
				fill(1);
				for(;discr_iters[3]!=discr_tens[3].end();discr_iters[3]+=st_sizes[3])
				{
				    fill(3);
				    cont_eval_type::third((const typename cont_eval_type::value_type)factor,(const std::vector<const typename cont_eval_type::value_type*>&)couplings,(std::vector<typename cont_eval_type::iterator>&)cont_iters,(const std::vector<const typename cont_eval_type::momentum_type*>&)momenta);
				    discr_eval_type::third((const typename discr_eval_type::value_type)factor,(const std::vector<const typename discr_eval_type::value_type*>&)couplings,(std::vector<typename discr_eval_type::iterator>&)discr_iters,(const std::vector<const typename discr_eval_type::momentum_type*>&)momenta,iterset);
				    if(!equal_sequences(cont_tens[2],discr_tens[2]))
				    {
					std::cerr<<"Mismatch continuous<->discrete colour flow vertices in 3rd recursive relation."<<std::endl;
					std::cerr<<"Failed colour configuration: ";
					print_cols(std::cerr);
					std::cerr<<cont_tens[2]<<","<<discr_tens[2]<<std::endl;
					return false;
				    }
				    cont_tens[2].reset();
				    discr_tens[2].reset();
				    reset(3);
				}
				discr_iters[3]=discr_tens[3].begin();
				reset(1);
			    }
			    discr_iters[1]=discr_tens[1].begin();
			    reset(0);
			}
			discr_iters[0]=discr_tens[0].begin();

			/* Fourth recursive relation check: */

			for(;discr_iters[0]!=discr_tens[0].end();discr_iters[0]+=st_sizes[0])
			{
			    fill(0);
			    for(;discr_iters[1]!=discr_tens[1].end();discr_iters[1]+=st_sizes[1])
			    {
				fill(1);
				for(;discr_iters[2]!=discr_tens[2].end();discr_iters[2]+=st_sizes[2])
				{
				    fill(2);
				    cont_eval_type::fourth((const typename cont_eval_type::value_type)factor,(const std::vector<const typename cont_eval_type::value_type*>&)couplings,(std::vector<typename cont_eval_type::iterator>&)cont_iters,(const std::vector<const typename cont_eval_type::momentum_type*>&)momenta);
				    discr_eval_type::fourth((const typename discr_eval_type::value_type)factor,(const std::vector<const typename discr_eval_type::value_type*>&)couplings,(std::vector<typename discr_eval_type::iterator>&)discr_iters,(const std::vector<const typename discr_eval_type::momentum_type*>&)momenta,iterset);
				    if(!equal_sequences(cont_tens[3],discr_tens[3]))
				    {
					std::cerr<<"Mismatch continuous<->discrete colour flow vertices in 4th recursive relation."<<std::endl;
					std::cerr<<"Failed colour configuration: ";
					print_cols(std::cerr);
					std::cerr<<cont_tens[3]<<","<<discr_tens[3]<<std::endl;
					return false;
				    }
				    cont_tens[3].reset();
				    discr_tens[3].reset();
				    reset(2);
				}
				discr_iters[2]=discr_tens[2].begin();
				reset(1);
			    }
			    discr_iters[1]=discr_tens[1].begin();
			    reset(0);
			}
			discr_iters[0]=discr_tens[0].begin();
		    }
		    return true;
		}

	    private:

		/* Random number generator: */

		generator_type gen;

		/* Cont tensors: */ 

		tensor_type cont_tens[4];

		/* Discr tensors: */

		tensor_type discr_tens[4];

		/* Iterators to the cont tensors: */

		std::vector<iterator> cont_iters;

		/* Iterators to the discr tensors: */

		std::vector<iterator> discr_iters;

		/* Coupling: */

		value_type c;
		std::vector<value_type*>couplings;

		/* Momenta: */

		momentum_type p[4];
		std::vector<momentum_type*>momenta;

		/* Bools denoting the adjoint and fundamental representation tensors: */

		bool adj_reps[4];
		bool fund_reps[4];

		/* Spacetime subtensor sizes: */

		unsigned st_sizes[4];
	};
    }
}

#include <Camgen/undef_args.h>

/* Testing program: */

#define N_COLS 3

#define GAUGE_GROUP SU<N_COLS>

using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing discrete colour structures........................"<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    test_utils::colour_tester<std::random,GAUGE_GROUP>coltest;
    std::cerr<<"Checking discrete d<I,J> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::d<adjoint_rep<GAUGE_GROUP>,1,2>,ssss >())
    {
	return 1;
    }

    std::cerr<<"..........done"<<std::endl;
    std::cerr<<"Checking discrete dd<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::dd<adjoint_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }

    std::cerr<<"..........done"<<std::endl;
    std::cerr<<"Checking discrete dd_plus<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::dd_plus<adjoint_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete dd_min<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::dd_min<adjoint_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete T<I,J,K> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::T<fundamental_rep<GAUGE_GROUP>,0,1,2>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete TT<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::TT<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete TT_plus<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::TT_plus<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete TT_min<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::TT_min<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete TT_contr<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::TT_contr<fundamental_rep<GAUGE_GROUP>,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete f<I,J,K> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::f<GAUGE_GROUP,0,1,2>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete ff_contr<I,J,K,L> colour structure..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::ff_contr<GAUGE_GROUP,0,1,2,3>,ssss >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete-colours ggg vertex..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::f<GAUGE_GROUP,0,1,2>,vvv >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;

    std::cerr<<"Checking discrete-colours gggg vertex..........";
    std::cerr.flush();
    if(!coltest.test< colour_tensor::ff_contr<GAUGE_GROUP,0,1,2,3>,vvvv >())
    {
	return 1;
    }
    std::cerr<<"..........done"<<std::endl;
}

