//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <ctime>
#include <Camgen/plt_config.h>
#include <Camgen/plt_script.h>
#include <Camgen/file_utils.h>
#include <Camgen/CM_algo.h>
#include <test_gen.h>
#include <QCDPbchabcc.h>
#include <QCDPbchcfcc.h>
#include <QCDPbchcfdc.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * *
 * Performance tests using multi-gluon amplitudes  *
 *                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef double value_type;

    license_print::disable();

    std::cout<<"-----------------------------------------------------"<<std::endl;
    std::cout<<"testing multi-gluon amplitude performance............"<<std::endl;
    std::cout<<"-----------------------------------------------------"<<std::endl;
    
    bool have_gp=plot_config::gnuplot_path!=NULL;
    file_utils::create_directory("test_output/speed_test");
    std::string filename("test_output/speed_test/speed_test");
    std::string fext=have_gp?".eps":".dat/.gp";
    value_type nglu,qcdadj,qcdcf,qcddc;
    
    data_wrapper* data_adj=have_gp?new data_wrapper(&nglu,&qcdadj):new data_wrapper(filename+".dat",&nglu,&qcdadj);
    data_wrapper* data_ccf=have_gp?new data_wrapper(&nglu,&qcdcf):new data_wrapper(filename+".dat",&nglu,&qcdcf);
    data_wrapper* data_dcf=have_gp?new data_wrapper(&nglu,&qcddc):new data_wrapper(filename+".dat",&nglu,&qcddc);
    
    value_type Ecm=500;
    value_type waiting_t=(value_type)10*CLOCKS_PER_SEC;
    unsigned Nmax=100;
    std::string process;

    process="g,g > g,g";
    nglu=4;

    {
	CM_algorithm<QCDPbchabcc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,2>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,2,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,2>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    process="g,g > g,g,g";
    nglu=5;

    {
	CM_algorithm<QCDPbchabcc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,3>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,3,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,3>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    process="g,g > g,g,g,g";
    nglu=6;

    {
	CM_algorithm<QCDPbchabcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,4,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,4>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_ccf->fill();
    }

    process="g,g > g,g,g,g,g";
    nglu=7;

    {
	CM_algorithm<QCDPbchabcc,2,5>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,5,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,5>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,5>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,5,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,5>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,5>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,5,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,5>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    process="g,g > g,g,g,g,g,g";
    nglu=8;

    {
	CM_algorithm<QCDPbchabcc,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,6,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,6,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,6>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,6,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,6>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    process="g,g > g,g,g,g,g,g,g";
    nglu=9;

    {
	CM_algorithm<QCDPbchabcc,2,7>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,7,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,7>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,7>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,7,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,7>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,7>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,7,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,7>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    process="g,g > g,g,g,g,g,g,g,g";
    nglu=10;

    {
	CM_algorithm<QCDPbchabcc,2,8>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,8,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,8>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,8>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,8,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,8>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,8>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,8,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,8>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    process="g,g > g,g,g,g,g,g,g,g,g";
    nglu=11;

    {
	CM_algorithm<QCDPbchabcc,2,9>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,9,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,9>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,9>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,9,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,9>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,9>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,9,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,9>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    process="g,g > g,g,g,g,g,g,g,g,g,g";
    nglu=12;

    {
	CM_algorithm<QCDPbchabcc,2,10>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchabcc,2,10,std::random>* gen=test_utils::test_generator_builder<QCDPbchabcc,2,10>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing adjoint QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdadj=timsec/N;
	data_adj->fill();
    }

    {
	CM_algorithm<QCDPbchcfcc,2,10>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfcc,2,10,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfcc,2,10>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcdcf=timsec/N;
	data_ccf->fill();
    }

    {
	CM_algorithm<QCDPbchcfdc,2,10>algo(process);
	algo.load();
	algo.construct();
	process_generator<QCDPbchcfdc,2,10,std::random>* gen=test_utils::test_generator_builder<QCDPbchcfdc,2,10>::create_generator(algo.get_tree_iterator(),Ecm);
	std::cout<<"Timing discrete colour-flow QCD for "<<process;
	std::cout.flush();
	long unsigned n=0;
	long unsigned N=0;
	std::clock_t start=std::clock();
	while((std::clock()-start)<waiting_t and N<Nmax)
	{
	    if((std::clock()-start)>n)
	    {
		std::cout<<".";
		std::cout.flush();
		n+=CLOCKS_PER_SEC;
	    }
	    gen->generate();
	    algo.evaluate();
	    ++N;
	}
	double timsec=((double)(std::clock()-start))/CLOCKS_PER_SEC;
	std::cout<<"done. Nr of events:"<<N<<", wall time: "<<timsec<<std::endl;
	qcddc=timsec/N;
	data_dcf->fill();
    }

    data_adj->write();
    data_stream* datastr_adj=new data_stream(data_adj,"1","2");
    datastr_adj->style="lines";
    datastr_adj->title="adjoint";

    data_ccf->write();
    data_stream* datastr_ccf=new data_stream(data_ccf,"1","2");
    datastr_ccf->style="lines";
    datastr_ccf->title="cont. color flow";

    data_dcf->write();
    data_stream* datastr_dcf=new data_stream(data_dcf,"1","2");
    datastr_dcf->style="lines";
    datastr_dcf->title="color flow";

    plot_script* plot=new plot_script(filename,"postscript color");
    plot->add_plot(datastr_adj);
    plot->add_plot(datastr_ccf);
    plot->add_plot(datastr_dcf);
    plot->ylog=true;
    plot->xlabel="Nr. of gluons";
    plot->ylabel="time per ME (s)";
    plot->plot();
    delete plot;
    std::cerr<<".........done, file "<<filename+fext<<" written."<<std::endl;
}

