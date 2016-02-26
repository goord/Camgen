//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/plt_config.h>
#include <Camgen/file_utils.h>
#include <Camgen/stdrand.h>
#include <Camgen/rn_strm.h>
#include <Camgen/val_gen_fac.h>
#include <s_int_tester.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for s-branching phase space integrals. *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing constrained density integrals...................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    
    /* Some useful type definitions: */

    typedef double value_type;
    typedef std::random rng_type;

    typedef Breit_Wigner<value_type,rng_type> Breit_Wigner_type;
    typedef power_law<value_type,rng_type> power_law_type;
    typedef uniform_value_generator<value_type,rng_type> uniform_type;
    typedef Dirac_delta<value_type,rng_type> Dirac_delta_type;

    file_utils::create_directory("test_output/s_int_test");
    
    //////////////////////////////////////////////////////////////////////////////

    /* Number of evaluated integrals: */
    
    std::size_t N_samples = 100;

    /* Monte Carlo points per batch: */

    std::size_t N_MC = 10000;
    
    /* Some useful constants: */
    
    //////////////////////////////////////////////////////////////////////////////

    value_type M_max	 = 500;
    value_type M1	 = 81.5;
    value_type W1	 = 2.5;
    value_type nu1	 = 1.05;
    value_type M_min1    = 5;
    value_type M2	 = 120;
    value_type W2	 = 5;
    value_type nu2	 = 1.75;
    value_type M_min2    = 10;
    
    //////////////////////////////////////////////////////////////////////////////

    {
	std::cerr<<"Checking integral Breit-Wigner("<<M1<<","<<W1<<") * Breit-Wigner("<<M2<<","<<W2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/BWBW1_int");
	Breit_Wigner_type* channel1=new Breit_Wigner_type(&M1,&W1);
	channel1->set_lower_bound(M_min1*M_min1);
	Breit_Wigner_type* channel2=new Breit_Wigner_type(&M2,&W2);
	channel2->set_lower_bound(M_min2*M_min2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M_min2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////

    M1	 	= 15.;
    W1	 	= 200.;
    M2	 	= 100;
    W2	 	= 5;
    
    /////////////////////////////////////////////////////////////////////////////
    
    {
	std::cerr<<"Checking integral Breit-Wigner("<<M1<<","<<W1<<") * Breit-Wigner("<<M2<<","<<W2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/BWBW2_int");
	Breit_Wigner_type* channel1=new Breit_Wigner_type(&M1,&W1);
	channel1->set_lower_bound(M_min1*M_min1);
	Breit_Wigner_type* channel2=new Breit_Wigner_type(&M2,&W2);
	channel2->set_lower_bound(M_min2*M_min2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M_min2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }
    
    //////////////////////////////////////////////////////////////////////////////

    M1	 	= 91.;
    W1	 	= 1.5;
    M2	 	= 91.;
    W2	 	= 1.5;
    
    /////////////////////////////////////////////////////////////////////////////
    
    {
	std::cerr<<"Checking integral Breit-Wigner("<<M1<<","<<W1<<") * Breit-Wigner("<<M2<<","<<W2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/BWBW3_int");
	Breit_Wigner_type* channel1=new Breit_Wigner_type(&M1,&W1);
	channel1->set_lower_bound(0.25*M1*M1);
	Breit_Wigner_type* channel2=new Breit_Wigner_type(&M2,&W2);
	channel2->set_lower_bound(0.25*M1*M1);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M_min2,M1+M2+50.,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////

    M1=81.5;
    M2=3.;
    
    //////////////////////////////////////////////////////////////////////////////

    {
	std::cerr<<"Checking integral Breit-Wigner("<<M1<<","<<W1<<") * power-law("<<M2<<","<<nu2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/BWpl_int");
	Breit_Wigner_type* channel1=new Breit_Wigner_type(&M1,&W1);
	channel1->set_lower_bound(M_min1*M_min1);
	power_law_type* channel2=new power_law_type(&M2,&nu2);
	channel2->set_lower_bound(M_min2*M_min2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M_min2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }

    {
	std::cerr<<"Checking integral Breit-Wigner("<<M1<<","<<W1<<") * uniform()............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/BWuni_int");
	Breit_Wigner_type* channel1=new Breit_Wigner_type(&M1,&W1);
	channel1->set_lower_bound(M_min1*M_min1);
	uniform_type* channel2=new uniform_type();
	channel2->set_lower_bound(M_min2*M_min2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M_min2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }
    
    {
	std::cerr<<"Checking integral Breit-Wigner("<<M1<<","<<W1<<") * Dirac-delta("<<M2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/BWDd_int");
	Breit_Wigner_type* channel1=new Breit_Wigner_type(&M1,&W1);
	channel1->set_lower_bound(M_min1*M_min1);
	Dirac_delta_type* channel2=new Dirac_delta_type(&M2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////

    M1=5;
    M_min1=10;
    M2=2.5;
    M_min2=5;
    
    //////////////////////////////////////////////////////////////////////////////

    {
	std::cerr<<"Checking integral power-law("<<M1<<","<<nu1<<") * power-law("<<M2<<","<<nu2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/plpl_int");
	power_law_type* channel1=new power_law_type(&M1,&nu1);
	channel1->set_lower_bound(M_min1*M_min1);
	power_law_type* channel2=new power_law_type(&M2,&nu2);
	channel2->set_lower_bound(M_min2*M_min2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M_min2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }
    
    //////////////////////////////////////////////////////////////////////////////

    M2=120;
    
    //////////////////////////////////////////////////////////////////////////////

    {
	std::cerr<<"Checking integral power-law("<<M1<<","<<nu1<<") * uniform()............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/pluni_int");
	power_law_type* channel1=new power_law_type(&M1,&nu1);
	channel1->set_lower_bound(M_min1*M_min1);
	uniform_type* channel2=new uniform_type();
	channel2->set_lower_bound(M_min2*M_min2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M_min2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }

    {
	std::cerr<<"Checking integral power-law("<<M1<<","<<nu1<<") * Dirac_delta("<<M2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/plDd_int");
	power_law_type* channel1=new power_law_type(&M1,&nu1);
	channel1->set_lower_bound(M_min1*M_min1);
	Dirac_delta_type* channel2=new Dirac_delta_type(&M2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1+M2,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }

    {
	std::cerr<<"Checking integral uniform() * Dirac_delta("<<M2<<")............";
	std::cerr.flush();
	std::string file("test_output/s_int_test/uDd_int");
	uniform_type* channel1=new uniform_type();
	channel1->set_lower_bound(M_min1*M_min1);
	Dirac_delta_type* channel2=new Dirac_delta_type(&M2);
	s_int_tester<value_type,std::random>tester(channel1,channel2);
	tester.run(M_min1,M_max,N_samples,N_MC,file,"postscript color");
	delete channel1;
	delete channel2;
	std::string output=plot_config::gnuplot_path==NULL?file+".dat/.gp":file+".eps";
	std::cerr<<"...........done, file "<<output<<" written."<<std::endl;
    }
}

