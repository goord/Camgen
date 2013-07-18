//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/CM_algo.h>
#include <QCDPbdhcfcc.h>
#include <susy_QCDPbdhcfcc.h>
#include <test_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Facility testing the relation between 6g and 4g2~g supersymmetric QCD       *
 * helicity amplitudes (see eq. (17) in Z. Kunszt, "Combined use of the CALKUL *
 * method and N=1 supersymmetry to calculate QCD six-parton processes", Nucl.  *
 * Phys. B271 pp.333.                                                          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef QCDPbdhcfcc::value_type value_type;
    typedef massless_spinor_factory<get_Dirac_algebra_type<QCDPbdhcfcc>::type,QCDPbdhcfcc::beam_direction,QCDPbdhcfcc::dimension> spinfac;

    license_print::disable();
    unsigned N_events=1000;
    value_type Ecm=500;

    std::cout<<"---------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing susy relation between 6g and 4g2~g helicity amplitudes......."<<std::endl;
    std::cout<<"---------------------------------------------------------------------"<<std::endl;

    std::cerr<<"Checking supersymmetric Ward identity for (+ + + - - +) helicity configurations..........";

    CM_algorithm<QCDPbdhcfcc,2,4>qcd_algo("g,g > g,g,g,g");
    qcd_algo.load();
    qcd_algo.construct();

    CM_algorithm<susy_QCDPbdhcfcc,2,4>susy_qcd_algo("g,g > g,g,~g,~g");
    susy_QCDPbdhcfcc::set_massless("~g");
    susy_qcd_algo.load();
    susy_qcd_algo.construct();
    process_generator<QCDPbdhcfcc,2,4,std::random>* gen;
    gen=test_utils::test_generator_builder<QCDPbdhcfcc,2,4>::create_generator(qcd_algo.get_tree_iterator(),Ecm);

    unsigned N_c=QCDPbdhcfcc::N_c;

    for(unsigned i=0;i<N_events;++i)
    {
	gen->generate();
	for(unsigned j=0;j<6;++j)
	{
	    qcd_algo.get_phase_space(j)->helicity()=(j==3 or j==4)?-1:1; // compare special helicity amplitudes
	    susy_qcd_algo.get_phase_space(j)->momentum()=qcd_algo.get_phase_space(j)->momentum();
	    susy_qcd_algo.get_phase_space(j)->helicity()=qcd_algo.get_phase_space(j)->helicity();
	    for(unsigned a=0;a<N_c;++a)
	    {
		for(unsigned b=0;b<N_c;++b)
		{
		    susy_qcd_algo.get_phase_space(j)->colour_entry(a*N_c+b)=qcd_algo.get_phase_space(j)->colour_entry(a*N_c+b);
		}
	    }
	}
	
	std::complex<value_type> factor=spinfac::s(qcd_algo.get_phase_space(3)->momentum(),qcd_algo.get_phase_space(4)->momentum())/spinfac::s(qcd_algo.get_phase_space(3)->momentum(),qcd_algo.get_phase_space(5)->momentum());
	
	if(!equals(std::abs(qcd_algo.evaluate()),std::abs(factor*susy_qcd_algo.evaluate())))
	{
	    std::cerr<<"Supersymmetric Ward identity violated by event "<<i<<std::endl;
	    return 1;
	}
    }
    std::cerr<<".........done."<<std::endl;
    return 0;
}

