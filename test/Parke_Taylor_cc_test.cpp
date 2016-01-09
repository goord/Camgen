//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <QCDPbdhcfcc.h>
#include <test_gen.h>
#include <Camgen/license_print.h>
#include <Camgen/CM_algo.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Facility testing the vanishing Parke-Taylor amplitudes in QCD: any gluon  *
 * amplitude with all helicities equal or all but one helicity equal must    *
 * vanish...                                                                 *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;


int main()
{
    typedef QCDPbdhcfcc::value_type value_type;
    license_print::disable();

    unsigned N_events=500;
    value_type Ecm=500;

    std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing Parke-Taylor gluon amplitudes in QCD with continuous colours.............."<<std::endl;
    std::cout<<"----------------------------------------------------------------------------------"<<std::endl;

    std::string process="g,g > g,g";
    {
	CM_algorithm<QCDPbdhcfcc,2,2>algo(process);
	algo.load();
	algo.construct();

	std::cerr<<"Checking (+ + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	process_generator<QCDPbdhcfcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbdhcfcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=1;
	    algo.get_phase_space(1)->helicity()=1;
	    algo.get_phase_space(2)->helicity()=1;
	    algo.get_phase_space(3)->helicity()=1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=1;
	    algo.get_phase_space(2)->helicity()=1;
	    algo.get_phase_space(3)->helicity()=1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    process="g,g > g,g,g";
    
    {
	CM_algorithm<QCDPbdhcfcc,2,3>algo(process);
	algo.load();
	algo.construct();

	std::cerr<<"Checking (+ + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	process_generator<QCDPbdhcfcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbdhcfcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + - + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + - + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - + - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - + - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    process="g,g > g,g,g,g";
    N_events/=2;
    
    {
	CM_algorithm<QCDPbdhcfcc,2,4>algo(process);
	algo.load();
	algo.construct();

	std::cerr<<"Checking (+ + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	process_generator<QCDPbdhcfcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbdhcfcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + - + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + - + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + - + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + - + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - + - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - + - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - + - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - + - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    process="g,g > g,g,g,g,g";
    
    {
	CM_algorithm<QCDPbdhcfcc,2,5>algo(process);
	algo.load();
	algo.construct();

	std::cerr<<"Checking (+ + + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	process_generator<QCDPbdhcfcc,2,5,std::random>* gen=test_utils::test_generator_builder<QCDPbdhcfcc,2,5>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + - + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + - + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + - + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + - + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + - + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + - + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + + - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + + + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - + - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - + - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - + - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - + - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - + - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - + - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
    
    process="g,g > g,g,g,g,g,g";
    N_events/=2;
    
    {
	CM_algorithm<QCDPbdhcfcc,2,6>algo(process);
	algo.load();
	algo.construct();

	std::cerr<<"Checking (+ + + + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	process_generator<QCDPbdhcfcc,2,6,std::random>* gen=test_utils::test_generator_builder<QCDPbdhcfcc,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + + + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + + + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - + + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - + + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + - + + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + - + + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + - + + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + - + + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + - + + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + - + + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + + - + +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + - + +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + + + - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + + - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ + + + + + + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()= 1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ + + + + + + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (+ - - - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()= 1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (+ - - - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- + - - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()= 1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- + - - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - + - - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - + - - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - + - - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()= 1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - + - - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - + - - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()= 1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - + - - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - + - -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()= 1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - + - -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - - + -) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()= 1;
	    algo.get_phase_space(7)->helicity()=-1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - - + -) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
	std::cerr<<"Checking (- - - - - - - +) amplitudes for "<<process<<"..........";
	std::cerr.flush();
	for(unsigned i=0;i<N_events;++i)
	{
	    gen->generate();
	    algo.get_phase_space(0)->helicity()=-1;
	    algo.get_phase_space(1)->helicity()=-1;
	    algo.get_phase_space(2)->helicity()=-1;
	    algo.get_phase_space(3)->helicity()=-1;
	    algo.get_phase_space(4)->helicity()=-1;
	    algo.get_phase_space(5)->helicity()=-1;
	    algo.get_phase_space(6)->helicity()=-1;
	    algo.get_phase_space(7)->helicity()= 1;
	    std::complex<value_type>amp=algo.evaluate();
	    if(!equals(std::abs(amp),(value_type)0))
	    {
		std::cerr<<"event "<<i<<": (- - - - - - - +) helicity configuration yields nonzero amplitude..."<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
}

