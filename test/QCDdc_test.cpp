//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <QCDPbchabcc.h>
#include <QCDPbchcfcc.h>
#include <QCDPbchcfdc.h>
#include <QCDPbchcfdc_clone.h>
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
    std::cout<<"testing QCD amplitudes with discrete colours.............."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;

    {
	CM_algorithm<QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchcfdc::set_unitary_gauge();
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
	    QCDPbchcfdc::set_R_xi_gauge();
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
	    QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchcfdc_clone,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc_clone,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc_clone,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::swap(algo.get_phase_space(0)->helicity_phase(-1),algo.get_phase_space(1)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(0)->helicity_phase(1),algo.get_phase_space(1)->helicity_phase(1));
	    std::swap(algo.get_phase_space(0)->colour(0),algo.get_phase_space(1)->colour(0));
	    std::swap(algo.get_phase_space(0)->colour(1),algo.get_phase_space(1)->colour(1));
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
	    std::swap(algo.get_phase_space(2)->colour(0),algo.get_phase_space(3)->colour(0));
	    std::swap(algo.get_phase_space(2)->colour(1),algo.get_phase_space(3)->colour(1));
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<QCDPbchcfcc,2,2>algocheck(process);
	algocheck.load();
	algocheck.construct();

	unsigned N=QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<algo.N_external;++k)
	    {
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
		algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		for(unsigned a=0;a<N*N;++a)
		{
		    algocheck.get_phase_space(k)->colour_entry(a)=0;
		}
		unsigned I=algo.get_phase_space(k)->colour(0);
		unsigned J=algo.get_phase_space(k)->colour(1);
		if(I==J)
		{
		    for(unsigned a=0;a<N*N;a+=(N+1))
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=-(value_type)1/(value_type)N;
		    }
		    algocheck.get_phase_space(k)->colour_entry(I*(N+1))+=(value_type)1;
		}
		else
		{
		    if(k<2)
		    {
			algocheck.get_phase_space(k)->colour_entry(I*N+J)=(value_type)1;
		    }
		    else
		    {
			algocheck.get_phase_space(k)->colour_entry(J*N+I)=(value_type)1;
		    }
		}
	    }
            QCDPbchcfcc::set_QCD_scale(QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="g,g > g,g,g";

    {
	CM_algorithm<QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchcfdc::set_unitary_gauge();
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
	    QCDPbchcfdc::set_R_xi_gauge();
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
	CM_algorithm<QCDPbchcfdc_clone,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc_clone,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc_clone,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::swap(algo.get_phase_space(0)->helicity_phase(-1),algo.get_phase_space(1)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(0)->helicity_phase(1),algo.get_phase_space(1)->helicity_phase(1));
	    std::swap(algo.get_phase_space(0)->colour(0),algo.get_phase_space(1)->colour(0));
	    std::swap(algo.get_phase_space(0)->colour(1),algo.get_phase_space(1)->colour(1));
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
	    std::swap(algo.get_phase_space(2)->colour(0),algo.get_phase_space(3)->colour(0));
	    std::swap(algo.get_phase_space(2)->colour(1),algo.get_phase_space(3)->colour(1));
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
	    std::swap(algo.get_phase_space(2)->colour(0),algo.get_phase_space(4)->colour(0));
	    std::swap(algo.get_phase_space(2)->colour(1),algo.get_phase_space(4)->colour(1));
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
	    std::swap(algo.get_phase_space(3)->colour(0),algo.get_phase_space(4)->colour(0));
	    std::swap(algo.get_phase_space(3)->colour(1),algo.get_phase_space(4)->colour(1));
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<QCDPbchcfcc,2,3>algocheck(process);
	algocheck.load();
	algocheck.construct();

	unsigned N=QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<algo.N_external;++k)
	    {
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
		algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		for(unsigned a=0;a<N*N;++a)
		{
		    algocheck.get_phase_space(k)->colour_entry(a)=0;
		}
		unsigned I=algo.get_phase_space(k)->colour(0);
		unsigned J=algo.get_phase_space(k)->colour(1);
		if(I==J)
		{
		    for(unsigned a=0;a<N*N;a+=(N+1))
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=-(value_type)1/(value_type)N;
		    }
		    algocheck.get_phase_space(k)->colour_entry(I*(N+1))+=(value_type)1;
		}
		else
		{
		    if(k<2)
		    {
			algocheck.get_phase_space(k)->colour_entry(I*N+J)=(value_type)1;
		    }
		    else
		    {
			algocheck.get_phase_space(k)->colour_entry(J*N+I)=(value_type)1;
		    }
		}
	    }
            QCDPbchcfcc::set_QCD_scale(QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="g,g > g,g,g,g";

    {
	CM_algorithm<QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchcfdc::set_unitary_gauge();
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
	    QCDPbchcfdc::set_R_xi_gauge();
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
	CM_algorithm<QCDPbchcfdc_clone,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc_clone,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc_clone,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::swap(algo.get_phase_space(0)->helicity_phase(-1),algo.get_phase_space(1)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(0)->helicity_phase(1),algo.get_phase_space(1)->helicity_phase(1));
	    std::swap(algo.get_phase_space(0)->colour(0),algo.get_phase_space(1)->colour(0));
	    std::swap(algo.get_phase_space(0)->colour(1),algo.get_phase_space(1)->colour(1));
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
	    std::swap(algo.get_phase_space(2)->colour(0),algo.get_phase_space(3)->colour(0));
	    std::swap(algo.get_phase_space(2)->colour(1),algo.get_phase_space(3)->colour(1));
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
	    std::swap(algo.get_phase_space(2)->colour(0),algo.get_phase_space(4)->colour(0));
	    std::swap(algo.get_phase_space(2)->colour(1),algo.get_phase_space(4)->colour(1));
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
	    std::swap(algo.get_phase_space(2)->colour(0),algo.get_phase_space(5)->colour(0));
	    std::swap(algo.get_phase_space(2)->colour(1),algo.get_phase_space(5)->colour(1));
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
	    std::swap(algo.get_phase_space(3)->colour(0),algo.get_phase_space(4)->colour(0));
	    std::swap(algo.get_phase_space(3)->colour(1),algo.get_phase_space(4)->colour(1));
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
	    std::swap(algo.get_phase_space(3)->colour(0),algo.get_phase_space(5)->colour(0));
	    std::swap(algo.get_phase_space(3)->colour(1),algo.get_phase_space(5)->colour(1));
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
	    std::swap(algo.get_phase_space(4)->colour(0),algo.get_phase_space(5)->colour(0));
	    std::swap(algo.get_phase_space(4)->colour(1),algo.get_phase_space(5)->colour(1));
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<QCDPbchcfcc,2,4>algocheck(process);
	algocheck.load();
	algocheck.construct();

	unsigned N=QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<algo.N_external;++k)
	    {
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
		algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		for(unsigned a=0;a<N*N;++a)
		{
		    algocheck.get_phase_space(k)->colour_entry(a)=0;
		}
		unsigned I=algo.get_phase_space(k)->colour(0);
		unsigned J=algo.get_phase_space(k)->colour(1);
		if(I==J)
		{
		    for(unsigned a=0;a<N*N;a+=(N+1))
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=-(value_type)1/(value_type)N;
		    }
		    algocheck.get_phase_space(k)->colour_entry(I*(N+1))+=(value_type)1;
		}
		else
		{
		    if(k<2)
		    {
			algocheck.get_phase_space(k)->colour_entry(I*N+J)=(value_type)1;
		    }
		    else
		    {
			algocheck.get_phase_space(k)->colour_entry(J*N+I)=(value_type)1;
		    }
		}
	    }
            QCDPbchcfcc::set_QCD_scale(QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    process="u,ubar > u,ubar";

    {
	CM_algorithm<QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchcfdc::set_unitary_gauge();
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
	    QCDPbchcfdc::set_R_xi_gauge();
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
	    QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchcfdc,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<QCDPbchcfcc,2,2>algocheck(process);
	algocheck.load();
	algocheck.construct();

	unsigned N=QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<algo.N_external;++k)
	    {
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
		algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		unsigned I=algo.get_phase_space(k)->colour(0);
		for(unsigned a=0;a<N;++a)
		{
		    if(a==I)
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=(value_type)1;
		    }
		    else
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=0;
		    }
		}
	    }
            QCDPbchcfcc::set_QCD_scale(QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="u,ubar > u,ubar,g";

    {
	CM_algorithm<QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchcfdc::set_unitary_gauge();
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
	    QCDPbchcfdc::set_R_xi_gauge();
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
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
	CM_algorithm<QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<QCDPbchcfcc,2,3>algocheck(process);
	algocheck.load();
	algocheck.construct();

	unsigned N=QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<algo.N_external;++k)
	    {
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
		algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		if(k==algo.N_external-1)
		{
		    for(unsigned a=0;a<N*N;++a)
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=0;
		    }
		    unsigned I=algo.get_phase_space(k)->colour(0);
		    unsigned J=algo.get_phase_space(k)->colour(1);
		    if(I==J)
		    {
			for(unsigned a=0;a<N*N;a+=(N+1))
			{
			    algocheck.get_phase_space(k)->colour_entry(a)=-(value_type)1/(value_type)N;
			}
			algocheck.get_phase_space(k)->colour_entry(I*(N+1))+=(value_type)1;
		    }
		    else
		    {
			if(k<2)
			{
			    algocheck.get_phase_space(k)->colour_entry(I*N+J)=(value_type)1;
			}
			else
			{
			    algocheck.get_phase_space(k)->colour_entry(J*N+I)=(value_type)1;
			}
		    }
		}
		else
		{
		    unsigned I=algo.get_phase_space(k)->colour(0);
		    for(unsigned a=0;a<N;++a)
		    {
			if(a==I)
			{
			    algocheck.get_phase_space(k)->colour_entry(a)=(value_type)1;
			}
			else
			{
			    algocheck.get_phase_space(k)->colour_entry(a)=(value_type)0;
			}
		    }
		}
	    }
            QCDPbchcfcc::set_QCD_scale(QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="u,ubar > u,ubar,u,ubar";

    {
	CM_algorithm<QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    QCDPbchcfdc::set_unitary_gauge();
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
	    QCDPbchcfdc::set_R_xi_gauge();
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
	    QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking crossing symmetries for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.evaluate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    std::swap(algo.get_phase_space(2)->colour(0),algo.get_phase_space(4)->colour(0));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping quarks 2 and 4 yields not minus sign"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping quarks 2 and 4 yields not minus sign"<<std::endl;
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
	    std::swap(algo.get_phase_space(3)->colour(0),algo.get_phase_space(5)->colour(0));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping quarks 3 and 5 yields not minus sign"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping quarks 3 and 5 yields not minus sign"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 3 as final current for "<<process<<"..........";
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
	CM_algorithm<QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<QCDPbchcfcc,2,4>algocheck(process);
	algocheck.load();
	algocheck.construct();

	unsigned N=QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<algo.N_external;++k)
	    {
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
		algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		unsigned I=algo.get_phase_space(k)->colour(0);
		for(unsigned a=0;a<N;++a)
		{
		    if(a==I)
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=(value_type)1;
		    }
		    else
		    {
			algocheck.get_phase_space(k)->colour_entry(a)=(value_type)0;
		    }
		}
	    }
            QCDPbchcfcc::set_QCD_scale(QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": discrete and continuous colour calculations yield different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
}

