//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <QEDPbdh.h>
#include <QEDPbKsdh.h>
#include <QEDPbpsdh.h>
#include <QEDWbdh.h>
#include <QEDWbKsdh.h>
#include <QEDWbpsdh.h>
#include <test_gen.h>
#include <Camgen/CM_algo.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Facility checking the spin-summed amplitudes for different spin vector  *
 * choices, using final states of tau leptons in QED processes...          *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef double value_type;

    license_print::disable();
    unsigned N_events=1000;
    value_type Ecm=30;
    std::vector<value_type>matrix_els(N_events);
    std::string process;

    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing massive spinor constructions in QED..............."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    process="tau+,tau- > gamma,gamma";

    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    matrix_els[i]=algo.evaluate_spin_sum();
	}
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbKsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbKsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbKsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Pauli basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbpsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbpsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbpsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Pauli basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis helicity spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: helicity spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbKsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbKsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbKsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbpsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbpsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbpsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="tau+,tau- > tau+,tau-";
    
    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    matrix_els[i]=algo.evaluate_spin_sum();
	}
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbKsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbKsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbKsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Pauli basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbpsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbpsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbpsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Pauli basis yield different result for event "<<i<<std::endl;
		std::cerr<<"polarised basis: "<<algo.evaluate_spin_sum()<<", helicity basis: "<<matrix_els[i]<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis helicity spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: helicity spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbKsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbKsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbKsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbpsdh,2,2>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbpsdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbpsdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="tau+,tau- > tau+,tau-,gamma";

    {
	CM_algorithm<QEDPbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    matrix_els[i]=algo.evaluate_spin_sum();
	}
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbKsdh,2,3>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbKsdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbKsdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Pauli basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbpsdh,2,3>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbpsdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbpsdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Pauli basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis helicity spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: helicity spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbKsdh,2,3>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbKsdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDWbKsdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbpsdh,2,3>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbpsdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDWbpsdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    process="tau+,tau- > tau+,tau-,tau+,tau-";

    {
	CM_algorithm<QEDPbdh,2,4>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    matrix_els[i]=algo.evaluate_spin_sum();
	}
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbKsdh,2,4>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbKsdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbKsdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Pauli basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbpsdh,2,4>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDPbpsdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbpsdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Pauli basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Pauli basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbdh,2,4>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis helicity spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: helicity spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbKsdh,2,4>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbKsdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDWbKsdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis KS-spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: KS-spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbpsdh,2,4>algo(process);
	algo.load();
	algo.construct();
	algo.sum_spins();
	process_generator<QEDWbpsdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDWbpsdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis polarised spinors in process "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    if(!equals(algo.evaluate_spin_sum(),matrix_els[i]))
	    {
		std::cerr<<"Error: polarised spinors in Weyl basis yield different result for event "<<i<<std::endl;
		return 1;
	    }
	}
	std::cerr<<"..........done."<<std::endl;
    }
}

