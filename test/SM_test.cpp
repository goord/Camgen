//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/CM_algo.h>
#include <Camgen/stdrand.h>
#include <SMWbhsch.h>
#include <test_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for the standard model class. The tests consist of rotating  *
 * the final current, crossing symmetries and (especially for testing the weak   *
 * sector) gauge invariance.                                                     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef model_type::value_type value_type;

    license_print::disable();
    Camgen::log.enable_level=log_level::error;

    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing standard-model amplitudes........................."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    value_type Ecm=1000;
    unsigned N_events=1000;
    std::vector< std::complex<value_type> >matrix_els(N_events);
    model_type::discard_widths();
    SMWbhsch::discard_widths();

    std::string process="u,dbar > mu+,nu_mu";

    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="d,dbar > mu+,mu-";

    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > W+,W-";

    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > gamma,gamma";

    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
	    width_scheme<model_type>::set_complex_mass_scheme();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > Z,Z";

    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
    
    process="W+,W- > h0,h0";

    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
    
    process="Z,Z > h0,h0";

    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 0 and 1 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(0)->momentum(),algo.get_phase_space(1)->momentum());
	    std::swap(algo.get_phase_space(0)->helicity_phase(-1),algo.get_phase_space(1)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(0)->helicity_phase(1),algo.get_phase_space(1)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 0 and 1 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 0 and 1 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,2>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,2>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,2,std::random>* gen=test_utils::test_generator_builder<model_type,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > Z,Z,Z";

    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
	    std::complex<value_type>ampug=algo.evaluate();
	    if(!equals(ampfg.real(),ampug.real()))
	    {
		std::cerr<<"event "<<i<<": real parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		std::cerr<<ampfg<<"\t"<<ampug<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),ampug.imag()))
	    {
		std::cerr<<"event "<<i<<": imaginary parts of amplitude differ in Feynman and unitary gauge..."<<std::endl;
		return 1;
	    }
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 4 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 4 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 4 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 3 and 4 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(4)->momentum());
	    std::swap(algo.get_phase_space(3)->helicity_phase(-1),algo.get_phase_space(4)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(3)->helicity_phase(1),algo.get_phase_space(4)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 3 and 4 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 3 and 4 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > Z,Z,h0";

    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();

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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > Z,h0,h0";

    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 3 and 4 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(4)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 3 and 4 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 3 and 4 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > h0,h0,h0";

    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 4 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(4)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 3 and 4 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(3)->momentum(),algo.get_phase_space(4)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > Z,Z,h0,h0";

    {
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 4 and 5 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(4)->momentum(),algo.get_phase_space(5)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 4 and 5 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 4 and 5 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="W+,W- > Z,Z,h0,h0,W+,W-";

    {
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 2 and 3 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(2)->momentum(),algo.get_phase_space(3)->momentum());
	    std::swap(algo.get_phase_space(2)->helicity_phase(-1),algo.get_phase_space(3)->helicity_phase(-1));
	    std::swap(algo.get_phase_space(2)->helicity_phase(1),algo.get_phase_space(3)->helicity_phase(1));
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 2 and 3 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Swapping particles 4 and 5 for "<<process<<"..........";
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::swap(algo.get_phase_space(4)->momentum(),algo.get_phase_space(5)->momentum());
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 4 and 5 yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": swapping particle 4 and 5 yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,6>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,6);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 6 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 6 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 6 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<model_type,2,6>algo(process,7);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 7 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 7 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 7 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="e+,e- > u,dbar,W-";

    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="u,dbar > W+,Z,gamma";

    {
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,3>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,3>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,3,std::random>* gen=test_utils::test_generator_builder<model_type,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
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

    process="e+,e- > u,dbar,mu-,nu_mubar";

    {
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	   gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,4>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,4>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<model_type,2,4,std::random>* gen=test_utils::test_generator_builder<model_type,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<SMWbhsch,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<SMWbhsch,2,4,std::random>* gen=test_utils::test_generator_builder<SMWbhsch,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": computation in Weyl basis yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": computation in Weyl basis yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="e+,e- > u,dbar,mu-,nu_mubar,e+,e-";

    {
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,6>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,6);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 6 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 6 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 6 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<model_type,2,6>algo(process,7);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 7 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 7 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 7 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<SMWbhsch,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<SMWbhsch,2,6,std::random>* gen=test_utils::test_generator_builder<SMWbhsch,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis computation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": computation in Weyl basis yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": computation in Weyl basis yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();

    process="g,g > b,bbar,mu-,nu_mubar,e+,nu_e";

    {
	CM_algorithm<model_type,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking gauge invariance for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    model_type::set_unitary_gauge();
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
	    model_type::set_R_xi_gauge();
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
	    model_type::set_Feynman_gauge();
	    matrix_els[i]=ampfg;
	}
	std::cerr<<".........done."<<std::endl;
    }

    random_number_stream<value_type,std::random>::reset_engine();
    
    {
	CM_algorithm<model_type,2,6>algo(process,1);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,2);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,3);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,4);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,5);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
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
	CM_algorithm<model_type,2,6>algo(process,6);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 6 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 6 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 6 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<model_type,2,6>algo(process,7);
	algo.load();
	algo.construct();
	process_generator<model_type,2,6,std::random>* gen=test_utils::test_generator_builder<model_type,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking particle 7 as final current for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),-matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": taking particle 7 as final current yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),-matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": taking particle 7 as final current yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    random_number_stream<value_type,std::random>::reset_engine();
    
    { 
	CM_algorithm<SMWbhsch,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<SMWbhsch,2,6,std::random>* gen=test_utils::test_generator_builder<SMWbhsch,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cerr<<"Checking Weyl basis comoputation for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    std::complex<value_type>ampfg=algo.evaluate();
	    if(!equals(ampfg.real(),matrix_els[i].real()))
	    {
		std::cerr<<"event "<<i<<": computation in Weyl basis yields different real part of amplitude"<<std::endl;
		return 1;
	    }
	    if(!equals(ampfg.imag(),matrix_els[i].imag()))
	    {
		std::cerr<<"event "<<i<<": computation in Weyl basis yields different imaginary part of amplitude"<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    return 0;
}

