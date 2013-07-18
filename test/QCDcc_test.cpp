//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <QCDPbchabcc.h>
#include <QCDPbchcfcc.h>
#include <QCDPbchabcc_clone.h>
#include <QCDPbchcfcc_clone.h>
#include <test_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for the computations of QCD-amplitudes with Camgen. The     *
 * tests consist of comparing results in different gauges, with different bases  *
 * of gamma matrices, permutations of the final particle in the tree,            *
 * (anti-)symmetry properties of the amplitude and using adjoint versus          *
 * colour-flow gluons.                                                           *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef QCDPbchabcc::value_type value_type;
    
    license_print::disable();
    Camgen::log.enable_level=log_level::error;

    /* Global definitions: */
    
    unsigned N_events=1000;
    value_type Ecm=100;
    std::vector< std::complex<value_type> >matrix_els(N_events);
    std::string process="g,g > g,g";

    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing QCD amplitudes with continuous colours............"<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;

    {
	CM_algorithm<QCDPbchabcc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchabcc::set_unitary_gauge();
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
	    QCDPbchabcc::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QCDPbchabcc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc_clone,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc_clone,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc_clone,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking 4-gluon vertex in "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": model with 4-gluon vertex yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": model with 4-gluon vertex yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::swap(algo.get_phase_space(0)->helicity_phase(-1),algo.get_phase_space(1)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(0)->helicity_phase(1),algo.get_phase_space(1)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(0)->colour_entry(a),algo.get_phase_space(1)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 0 and 1 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 0 and 1 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(2)->colour_entry(a),algo.get_phase_space(3)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchabcc,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc_clone,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc_clone,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc_clone,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow 4-gluon vertex in "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment with 4-gluon vertex yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment with 4-gluon vertex yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    process = "g,g > g,g,g";

    {
	CM_algorithm<QCDPbchabcc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchabcc::set_unitary_gauge();
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
	    QCDPbchabcc::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QCDPbchabcc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchabcc_clone,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc_clone,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc_clone,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking 4-gluon vertex in "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": model with 4-gluon vertex yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": model with 4-gluon vertex yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::swap(algo.get_phase_space(0)->helicity_phase(-1),algo.get_phase_space(1)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(0)->helicity_phase(1),algo.get_phase_space(1)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(0)->colour_entry(a),algo.get_phase_space(1)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 0 and 1 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 0 and 1 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(2)->colour_entry(a),algo.get_phase_space(3)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(2)->colour_entry(a),algo.get_phase_space(4)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 4 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 4 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(3)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(3)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(3)->colour_entry(a),algo.get_phase_space(4)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 3 and 4 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 3 and 4 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchabcc,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc_clone,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc_clone,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc_clone,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow 4-gluon vertex in "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment with 4-gluon vertex yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment with 4-gluon vertex yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    process = "g,g > g,g,g,g";

    {
	CM_algorithm<QCDPbchabcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchabcc::set_unitary_gauge();
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
	    QCDPbchabcc::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QCDPbchabcc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchabcc_clone,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc_clone,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc_clone,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking 4-gluon vertex in "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": model with 4-gluon vertex yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": model with 4-gluon vertex yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::swap(algo.get_phase_space(0)->helicity_phase(-1),algo.get_phase_space(1)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(0)->helicity_phase(1),algo.get_phase_space(1)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(0)->colour_entry(a),algo.get_phase_space(1)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 0 and 1 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 0 and 1 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(2)->colour_entry(a),algo.get_phase_space(3)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(2)->colour_entry(a),algo.get_phase_space(4)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 4 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 4 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(5)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(5)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(5)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(2)->colour_entry(a),algo.get_phase_space(5)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 5 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 2 and 5 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(3)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(3)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(3)->colour_entry(a),algo.get_phase_space(4)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 3 and 4 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 3 and 4 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(5)->momentum());
	    std::swap(algo.get_phase_space(3)->helicity_phase(-1),algo.get_phase_space(5)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(3)->helicity_phase(1),algo.get_phase_space(5)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(3)->colour_entry(a),algo.get_phase_space(5)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 3 and 5 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 3 and 5 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(4)->momentum(),algo.get_phase_space(5)->momentum());
	    std::swap(algo.get_phase_space(4)->helicity_phase(-1),algo.get_phase_space(5)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(4)->helicity_phase(1),algo.get_phase_space(5)->helicity_phase(1));
	    for(unsigned a=0;a<(QCDPbchabcc::N_c*QCDPbchabcc::N_c-1);++a)
	    {
		std::swap(algo.get_phase_space(4)->colour_entry(a),algo.get_phase_space(5)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 4 and 5 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping gluons 4 and 5 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchabcc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 5 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 5 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 5 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc_clone,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc_clone,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc_clone,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow 4-gluon vertex in "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment with 4-gluon vertex yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment with 4-gluon vertex yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    process = "u,ubar > u,ubar";

    {
	CM_algorithm<QCDPbchabcc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchabcc::set_unitary_gauge();
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
	    QCDPbchabcc::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QCDPbchabcc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QCDPbchabcc,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    process = "u,ubar > u,ubar,g";

    {
	CM_algorithm<QCDPbchabcc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchabcc::set_unitary_gauge();
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
	    QCDPbchabcc::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QCDPbchabcc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    {
	CM_algorithm<QCDPbchabcc,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    process = "u,ubar > u,ubar,u,ubar";

    {
	CM_algorithm<QCDPbchabcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchabcc::set_unitary_gauge();
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
	    QCDPbchabcc::set_R_xi_gauge();
	    std::complex<value_type>amprg=algo.evaluate();
	    if(!equals(ampfg.real(),amprg.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),amprg.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    QCDPbchabcc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    for(unsigned a=0;a<QCDPbchabcc::N_c;++a)
	    {
		std::swap(algo.get_phase_space(2)->colour_entry(a),algo.get_phase_space(4)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": amplitude real part not antisymmetric under swapping quarks 2 and 4"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude imaginary part not antisymmetric under swapping quarks 2 and 4"<<std::endl;
		return 1;
	    }
	}
	random_number_stream<value_type,std::random>::reset_engine();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(5)->momentum());
	    std::swap(algo.get_phase_space(3)->helicity_phase(-1),algo.get_phase_space(5)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(3)->helicity_phase(1),algo.get_phase_space(5)->helicity_phase(1));
	    for(unsigned a=0;a<QCDPbchabcc::N_c;++a)
	    {
		std::swap(algo.get_phase_space(3)->colour_entry(a),algo.get_phase_space(5)->colour_entry(a));
	    }
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": amplitude real part not antisymmetric under swapping quarks 3 and 5"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": amplitude imaginary part not antisymmetric under swapping quarks 3 and 5"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchabcc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchabcc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    { 
	CM_algorithm<QCDPbchcfcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking colour-flow amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
}

