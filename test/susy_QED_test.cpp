//
// This file is part of the CAMORRA library.
// Copyright (C) 2010 Gijs van den Oord.
// CAMORRA is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/susy_QED.h>
#include <susy_QEDWb.h>
#include <test_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for the computations of supersymmetric QED-amplitudes with *
 * Camgen. The tests consist of comparing results in different gauges, with   *
 * different bases of gamma matrices and permutations of the final particle in *
 * the tree, and some (anti-)symmetry properties of the amplitude.             *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef susy_QED::value_type value_type;

    license_print::disable();
    Camgen::log.enable_level=log_level::error;

    susy_QED::M_seL=230.8;
    susy_QED::M_seR=157.5;
    susy_QED::M_smuL=230.8;
    susy_QED::M_smuR=157.5;
    susy_QED::M_stauL=152.2;
    susy_QED::M_stauR=232.4;
    susy_QED::M_chi=117.9;

    susy_QEDWb::M_seL=susy_QED::M_seL;
    susy_QEDWb::M_seR=susy_QED::M_seR;
    susy_QEDWb::M_smuL=susy_QED::M_smuL;
    susy_QEDWb::M_smuR=susy_QED::M_smuR;
    susy_QEDWb::M_stauL=susy_QED::M_stauL;
    susy_QEDWb::M_stauR=susy_QED::M_stauR;
    susy_QEDWb::M_chi=susy_QED::M_chi;
    
    unsigned N_events=1000;
    value_type Ecm=1000;
    std::vector< std::complex<value_type> >matrix_els(N_events);
    std::string process;

    process="e+,e- > ~mu_L+,~mu_L-";
    
    {
	CM_algorithm<susy_QED,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="e+,e- > ~e_L+,~e_L-";
    
    {
	CM_algorithm<susy_QED,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="e-,~chi_10 > e-,~chi_10";
    
    {
	CM_algorithm<susy_QED,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		algo.print(std::cout);
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > ~e_L-,~e_L-";
    
    {
	CM_algorithm<susy_QED,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="e+,e- > ~mu_L+,~mu_L-,gamma";
    
    {
	CM_algorithm<susy_QED,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="e+,e- > ~e_L+,~e_L-,gamma";
    
    {
	CM_algorithm<susy_QED,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > ~e_L-,~e_L-,gamma";
    
    {
	CM_algorithm<susy_QED,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="e+,e- > ~mu_L+,~mu_L-,~chi_10";
    
    {
	CM_algorithm<susy_QED,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="e+,e- > ~e_L+,~e_L-,~chi_10";
    
    {
	CM_algorithm<susy_QED,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > ~e_L-,~e_L-,~chi_10";
    
    {
	CM_algorithm<susy_QED,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="e+,e- > ~mu_L+,~mu_L-,mu+,mu-";
    
    {
	CM_algorithm<susy_QED,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 5 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="e+,e- > ~e_L+,~e_L-,~e_R+,~e_R-";
    
    {
	CM_algorithm<susy_QED,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 5 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > ~e_L-,~e_L-,gamma,gamma";
    
    {
	CM_algorithm<susy_QED,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 5 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > ~e_L-,~e_L-,gamma,~chi_10";
    
    {
	CM_algorithm<susy_QED,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 5 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > e-,e-,~chi_10,~chi_10";
    
    {
	CM_algorithm<susy_QED,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking gauge invariance for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    susy_QED::set_unitary_gauge();
	    std::complex<double>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_R_xi_gauge();
	    std::complex<double>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cout<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cout<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    susy_QED::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<susy_QED,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QED,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QED,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QED,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking particle 5 as final current for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": taking particle 5 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<susy_QEDWb,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QEDWb,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QEDWb,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cout.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<double>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cout<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cout<<".........done."<<std::endl;
    }
}
