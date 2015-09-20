//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_TEST_GEN_H_
#define CAMGEN_TEST_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility class making uniform momentum, helicity and colour generators to  *
 * fill up the degrees of freedom of test amplitudes.                        *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/uni_evtgen_fac.h>
#include <Camgen/stdrand.h>

namespace Camgen
{
    namespace test_utils
    {
	template<class model_type,int N_in,int N_out>class test_generator_builder;

	template<class model_type,int N_out>class test_generator_builder<model_type,1,N_out>
	{
	    public:

		typedef typename CM_algorithm<model_type,1,N_out>::tree_iterator tree_iterator;

		static process_generator<model_type,1,N_out,std::random>* create_generator(tree_iterator it)
		{
		    set_initial_state_type(initial_states::partonic);
		    set_phase_space_generator_type(phase_space_generators::uniform);
		    set_helicity_generator_type(helicity_generators::uniform);
		    set_colour_generator_type(colour_generators::flow_sampling);
		    uniform_process_generator_factory<model_type,1,N_out,std::random> factory;
		    return factory.create_generator(it);
		}

		static void fill(tree_iterator it)
		{
		    process_generator<model_type,1,N_out,std::random>* gen=create_generator(it);
		    gen->generate();
		    delete gen;
		}

		static void fill(CM_algorithm<model_type,1,N_out>& algo)
		{
		    fill(algo.get_tree_iterator());
		}
	};

	template<class model_type,int N_out>class test_generator_builder<model_type,2,N_out>
	{
	    public:

		typedef typename CM_algorithm<model_type,2,N_out>::tree_iterator tree_iterator;
		typedef typename model_type::value_type value_type;

		static process_generator<model_type,2,N_out,std::random>* create_generator(tree_iterator it,const value_type& E)
		{
		    set_initial_state_type(initial_states::partonic);
		    set_phase_space_generator_type(phase_space_generators::uniform);
		    set_helicity_generator_type(helicity_generators::uniform);
		    set_colour_generator_type(colour_generators::flow_sampling);
		    set_beam_energy(1,0.5*E);
		    set_beam_energy(2,0.5*E);
		    uniform_process_generator_factory<model_type,2,N_out,std::random> factory;
		    return factory.create_generator(it);
		}

		static void fill(tree_iterator it,const value_type& E)
		{
		    process_generator<model_type,2,N_out,std::random> gen=create_generator(it,E);
		    gen->generate();
		    delete gen;
		}

		static void fill(CM_algorithm<model_type,2,N_out>& algo,const value_type& E)
		{
		    fill(algo.get_tree_iterator(),E);
		}
	};
    }
}

#endif /*CAMGEN_TEST_GEN_H_*/

