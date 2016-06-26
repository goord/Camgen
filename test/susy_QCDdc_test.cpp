//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/CM_algo.h>
#include <susy_QCDPbchabcc.h>
#include <susy_QCDPbchcfcc.h>
#include <susy_QCDPbchcfdc.h>
#include <test_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for the computations of susyQCD-amplitudes with Camgen. The *
 * tests consist of comparing results in different gauges, with different bases  *
 * of gamma matrices, permutations of the final particle in the tree,            *
 * (anti-)symmetry properties of the amplitude and using adjoint versus          *
 * colour-flow gluons.                                                           *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef susy_QCDPbchabcc::value_type value_type;

    license_print::disable();
    Camgen::log.enable_level=log_level::error;

    susy_QCDPbchabcc::M_sdL=666.2;
    susy_QCDPbchabcc::M_sdR=639.0;
    susy_QCDPbchabcc::M_suL=660.3;
    susy_QCDPbchabcc::M_suR=644.3;
    susy_QCDPbchabcc::M_ssL=666.2;
    susy_QCDPbchabcc::M_ssR=639.0;
    susy_QCDPbchabcc::M_scL=660.3;
    susy_QCDPbchabcc::M_scR=644.3;
    susy_QCDPbchabcc::M_sbL=599.0;
    susy_QCDPbchabcc::M_sbR=636.6;
    susy_QCDPbchabcc::M_stL=446.9;
    susy_QCDPbchabcc::M_stR=670.9;
    susy_QCDPbchabcc::M_sg =717.5;

    susy_QCDPbchcfcc::M_sdL=susy_QCDPbchabcc::M_sdL;
    susy_QCDPbchcfcc::M_sdR=susy_QCDPbchabcc::M_sdR;
    susy_QCDPbchcfcc::M_suL=susy_QCDPbchabcc::M_suL;
    susy_QCDPbchcfcc::M_suR=susy_QCDPbchabcc::M_suR;
    susy_QCDPbchcfcc::M_ssL=susy_QCDPbchabcc::M_ssL;
    susy_QCDPbchcfcc::M_ssR=susy_QCDPbchabcc::M_ssR;
    susy_QCDPbchcfcc::M_scL=susy_QCDPbchabcc::M_scL;
    susy_QCDPbchcfcc::M_scR=susy_QCDPbchabcc::M_scR;
    susy_QCDPbchcfcc::M_sbL=susy_QCDPbchabcc::M_sbL;
    susy_QCDPbchcfcc::M_sbR=susy_QCDPbchabcc::M_sbR;
    susy_QCDPbchcfcc::M_stL=susy_QCDPbchabcc::M_stL;
    susy_QCDPbchcfcc::M_stR=susy_QCDPbchabcc::M_stR;
    susy_QCDPbchcfcc::M_sg =susy_QCDPbchabcc::M_sg;

    susy_QCDPbchcfdc::M_sdL=susy_QCDPbchabcc::M_sdL;
    susy_QCDPbchcfdc::M_sdR=susy_QCDPbchabcc::M_sdR;
    susy_QCDPbchcfdc::M_suL=susy_QCDPbchabcc::M_suL;
    susy_QCDPbchcfdc::M_suR=susy_QCDPbchabcc::M_suR;
    susy_QCDPbchcfdc::M_ssL=susy_QCDPbchabcc::M_ssL;
    susy_QCDPbchcfdc::M_ssR=susy_QCDPbchabcc::M_ssR;
    susy_QCDPbchcfdc::M_scL=susy_QCDPbchabcc::M_scL;
    susy_QCDPbchcfdc::M_scR=susy_QCDPbchabcc::M_scR;
    susy_QCDPbchcfdc::M_sbL=susy_QCDPbchabcc::M_sbL;
    susy_QCDPbchcfdc::M_sbR=susy_QCDPbchabcc::M_sbR;
    susy_QCDPbchcfdc::M_stL=susy_QCDPbchabcc::M_stL;
    susy_QCDPbchcfdc::M_stR=susy_QCDPbchabcc::M_stR;
    susy_QCDPbchcfdc::M_sg =susy_QCDPbchabcc::M_sg;
    
    unsigned N_events=1000;
    value_type Ecm=3000;
    std::vector< std::complex<value_type> >matrix_els(N_events);
    std::string process;

    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing susyQCD amplitudes with discrete colours.........."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    process="g,g > ~g,~g";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,2>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    process="g,g > ~u_L+,~u_L-";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
		std::cerr<<ampfg<<"<------->"<<matrix_els[i]<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,2>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<2;++k)
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
	    for(unsigned k=2;k<4;++k)
	    {
		algocheck.get_phase_space(k)->helicity_phase(0)=1;
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="u,u > ~u_R+,~u_L+";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
		std::cerr<<ampfg<<"<------->"<<matrix_els[i]<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,2>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<4;++k)
	    {
		if(k<2)
		{
		    algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		    algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		}
		else
		{
		    algocheck.get_phase_space(k)->helicity_phase(0)=1;
		}
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="g,g > ~g,~g,g";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,3>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    process="g,g > g,~u_L+,~u_L-";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
		std::cerr<<ampfg<<"<------->"<<matrix_els[i]<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,3>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<3;++k)
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
	    for(unsigned k=3;k<5;++k)
	    {
		algocheck.get_phase_space(k)->helicity_phase(0)=1;
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="u,u > ~u_R+,~u_L+,g";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
		std::cerr<<ampfg<<"<------->"<<matrix_els[i]<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,3>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<4;++k)
	    {
		if(k<2)
		{
		    algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		    algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		}
		else
		{
		    algocheck.get_phase_space(k)->helicity_phase(0)=1;
		}
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
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
	    algocheck.get_phase_space(4)->momentum()=algo.get_phase_space(4)->momentum();
	    algocheck.get_phase_space(4)->helicity_phase(-1)=algo.get_phase_space(4)->helicity_phase(-1);
	    algocheck.get_phase_space(4)->helicity_phase(1)=algo.get_phase_space(4)->helicity_phase(1);
	    for(unsigned a=0;a<N*N;++a)
	    {
		algocheck.get_phase_space(4)->colour_entry(a)=0;
	    }
	    unsigned I=algo.get_phase_space(4)->colour(0);
	    unsigned J=algo.get_phase_space(4)->colour(1);
	    if(I==J)
	    {
		for(unsigned a=0;a<N*N;a+=(N+1))
		{
		    algocheck.get_phase_space(4)->colour_entry(a)=-(value_type)1/(value_type)N;
		}
		algocheck.get_phase_space(4)->colour_entry(I*(N+1))+=(value_type)1;
	    }
	    else
	    {
		algocheck.get_phase_space(4)->colour_entry(J*N+I)=(value_type)1;
	    }
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    process="g,g > ~g,~g,g,g";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,4>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    process="g,g > g,g,~u_L+,~u_L-";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
		std::cerr<<ampfg<<"<------->"<<matrix_els[i]<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,4>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<4;++k)
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
	    for(unsigned k=4;k<6;++k)
	    {
		algocheck.get_phase_space(k)->helicity_phase(0)=1;
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="u,u > ~u_L+,~u_L+,~u_L+,~u_L-";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,4>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<6;++k)
	    {
		if(k<2)
		{
		    algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		    algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		}
		else
		{
		    algocheck.get_phase_space(k)->helicity_phase(0)=1;
		}
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		algo.print(std::cerr);
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();

    process="u,u > ~u_R+,~u_L+,~u_L+,~u_R-";
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    susy_QCDPbchcfdc::set_unitary_gauge();
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
	    susy_QCDPbchcfdc::set_R_xi_gauge();
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
	    susy_QCDPbchcfdc::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<susy_QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<susy_QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<susy_QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	CM_algorithm<susy_QCDPbchcfcc,2,4>algocheck(process);
	algocheck.load();
	algocheck.construct();
	unsigned N=susy_QCDPbchcfdc::N_c;
	std::cerr<<"Checking discrete against continuous colours for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    for(unsigned k=0;k<6;++k)
	    {
		if(k<2)
		{
		    algocheck.get_phase_space(k)->helicity_phase(-1)=algo.get_phase_space(k)->helicity_phase(-1);
		    algocheck.get_phase_space(k)->helicity_phase(1)=algo.get_phase_space(k)->helicity_phase(1);
		}
		else
		{
		    algocheck.get_phase_space(k)->helicity_phase(0)=1;
		}
		algocheck.get_phase_space(k)->momentum()=algo.get_phase_space(k)->momentum();
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
            susy_QCDPbchcfcc::set_QCD_scale(susy_QCDPbchcfdc::QCD_scale);
	    std::complex<value_type>ampfg=algo.evaluate();
	    std::complex<value_type>ampcheck=algocheck.evaluate();
	    if(!equals(ampfg.real(),ampcheck.real()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different real part of amplitude"<<std::endl;
		algo.print(std::cerr);
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampcheck.imag()))
	    {
		std::cerr<<"event "<<i<<": colour flow treatment yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
}

