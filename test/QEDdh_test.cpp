//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <QEDPbdh.h>
#include <QEDWbdh.h>
#include <QEDPbch.h>
#include <QEDWbch.h>
#include <test_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for the computations of QED-amplitudes with Camgen. The     *
 * tests consist of comparing results in different gauges, with different bases  *
 * of gamma matrices and permutations of the final particle in the tree, and     *
 * some (anti-)symmetry properties of the amplitude.                             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef QEDPbdh::value_type value_type;

    license_print::disable();
    Camgen::log.enable_level=log_level::error;

    /* Global definitions: */

    unsigned N_events=1000;
    value_type Ecm=100;
    std::vector< std::complex<value_type> >matrix_els(N_events);
    std::string process="e+,e- > e+,e-";

    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing QED amplitudes with discrete helicities..........."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;

    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QEDPbdh::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude: ";
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDWbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDWbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    process="e+,e- > gamma,gamma";

    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QEDPbdh::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QEDPbdh,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDWbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDWbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking final-state symmetry for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>amp=algo.evaluate();
	    std::swap(algo.get_phase_space(2)->helicity(),algo.get_phase_space(3)->helicity());
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::complex<value_type>ampswap=algo.evaluate();

	    if(!equals(amp.real(),ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state photons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state photons"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > e-,e-";

    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QEDPbdh::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDWbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking initial-state antisymmetry for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>amp=algo.evaluate();
	    std::swap(algo.get_phase_space(0)->helicity(),algo.get_phase_space(1)->helicity());
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::complex<value_type>ampswap=algo.evaluate();

	    if(!equals(amp.real(),-ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping initial-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),-ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping initial-state electrons"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,2,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking final-state antisymmetry for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>amp=algo.evaluate();
	    std::swap(algo.get_phase_space(2)->helicity(),algo.get_phase_space(3)->helicity());
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::complex<value_type>ampswap=algo.evaluate();

	    if(!equals(amp.real(),-ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),-ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="e+,e- > e+,e-,gamma";

    {
	CM_algorithm<QEDPbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QEDPbdh::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QEDPbdh,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDWbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDWbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="e+,e- > gamma,gamma,gamma";

    {
	CM_algorithm<QEDPbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QEDPbdh::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QEDPbdh,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDWbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDWbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking final-state symmetry for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>amp=algo.evaluate();

	    std::swap(algo.get_phase_space(2)->helicity(),algo.get_phase_space(3)->helicity());
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::complex<value_type>ampswap=algo.evaluate();

	    if(!equals(amp.real(),ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    std::swap(algo.get_phase_space(2)->helicity(),algo.get_phase_space(4)->helicity());
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    ampswap=algo.evaluate();

	    if(!equals(amp.real(),ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    std::swap(algo.get_phase_space(3)->helicity(),algo.get_phase_space(4)->helicity());
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(4)->momentum());
	    ampswap=algo.evaluate();

	    if(!equals(amp.real(),ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="e-,e- > e-,e-,gamma";

    {
	CM_algorithm<QEDPbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QEDPbdh::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDWbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDWbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking initial-state antisymmetry for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>amp=algo.evaluate();
	    std::swap(algo.get_phase_space(0)->helicity(),algo.get_phase_space(1)->helicity());
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::complex<value_type>ampswap=algo.evaluate();

	    if(!equals(amp.real(),-ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping initial-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),-ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping initial-state electrons"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,3,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking final-state antisymmetry for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>amp=algo.evaluate();
	    std::swap(algo.get_phase_space(2)->helicity(),algo.get_phase_space(3)->helicity());
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::complex<value_type>ampswap=algo.evaluate();

	    if(!equals(amp.real(),-ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),-ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="e+,e- > e+,e-,e+,e-";

    {
	CM_algorithm<QEDPbdh,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QEDPbdh::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and R-xi gauge..."<<std::endl;
		return 1;
	    }
	    QEDPbdh::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QEDPbdh,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 1 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 1 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 2 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 2 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 3 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 4 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 4 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDPbdh,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 5 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 5 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 5 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QEDWbdh,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDWbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDWbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl-basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": Weyl basis computation yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QEDPbdh,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QEDPbdh,2,4,std::random>* gen=test_utils::test_generator_builder<QEDPbdh,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking final-state antisymmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>amp=algo.evaluate();
	    std::swap(algo.get_phase_space(2)->helicity(),algo.get_phase_space(4)->helicity());
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    std::complex<value_type>ampswap=algo.evaluate();

	    if(!equals(amp.real(),-ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),-ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    std::swap(algo.get_phase_space(3)->helicity(),algo.get_phase_space(5)->helicity());
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(5)->momentum());
	    ampswap=algo.evaluate();
	    if(!equals(amp.real(),ampswap.real()))
	    {
		std::cerr<<"event "<<i<<": amplitude's real part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	    if(!equals(amp.imag(),ampswap.imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude's imaginary part not antisymmetric under swapping final-state electrons"<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::cerr<<".........done."<<std::endl;
    }
}


