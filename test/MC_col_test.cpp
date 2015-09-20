//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/license_print.h>
#include <Camgen/CM_algo.h>
#include <Camgen/ps_copy.h>
#include <Camgen/stdrand.h>
#include <Camgen/uni_psgen_fac.h>
#include <Camgen/helgen_fac.h>
#include <Camgen/uni_cols.h>
#include <Camgen/qcd_cols.h>
#include <Camgen/plt_config.h>
#include <Camgen/plt_script.h>
#include <Camgen/QCD.h>
#include <Camgen/SM.h>
#include <QCDPbchcfcc.h>
#include <QCDPbchabcc.h>
#include <QCDPbchabdc.h>

using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing colour summation algorithms......................................"<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    
    typedef double value_type;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    /* Number of events per testrun: */

    std::size_t N_events=10000;

    /* Number of sampled data points in cross section plots: */

    std::size_t N_points=100;

    /* Center-of-mass energy definition: */

    value_type Ecm=100;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////

    std::size_t N_batch=N_events/N_points;
    set_initial_state_type(initial_states::partonic);
    set_phase_space_generator_type(phase_space_generators::uniform);
    set_helicity_generator_type(helicity_generators::uniform);
    set_first_beam_energy(0.5*Ecm);
    set_second_beam_energy(0.5*Ecm);
   
    {
	typedef QCD model_type;
	typedef QCD::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,2,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,2,std::random> hel_factory;
	
	std::string process("u,ubar > g,g");
	std::cerr<<"Checking full colour summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,2,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	algo.sum_colours();
	value_type M=algo.evaluate_sum();
	value_type M2=algo.evaluate_colour_sum();
	if(!equals(M,M2))
	{
	    std::cerr<<"Colour sum does not equal dof sum..."<<std::endl;
	    return 1;
	}
	value_type M3(0);
	for(std::size_t a1=0;a1<model_type::N_c;++a1)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->colour(0)=a1;
	    for(std::size_t a2=0;a2<model_type::N_c;++a2)
	    {              
		algo.get_tree_iterator()->get_phase_space(1)->colour(0)=a2;
		for(std::size_t a3=0;a3<model_type::N_c;++a3)
		{              
		    algo.get_tree_iterator()->get_phase_space(2)->colour(0)=a3;
		    for(std::size_t a4=0;a4<model_type::N_c;++a4)
		    {              
			algo.get_tree_iterator()->get_phase_space(2)->colour(1)=a4;
			for(std::size_t a5=0;a5<model_type::N_c;++a5)
			{              
			    algo.get_tree_iterator()->get_phase_space(3)->colour(0)=a5;
			    for(std::size_t a6=0;a6<model_type::N_c;++a6)
			    {
				algo.get_tree_iterator()->get_phase_space(3)->colour(1)=a6;
				M3+=algo.evaluate2();
			    }
			}
		    }
		}
	    }
	}
	M3/=(value_type)4;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    {
	typedef QCD model_type;
	typedef model_type::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,3,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,3,std::random> hel_factory;
	
	std::string process("u,ubar > u,ubar,g");
	std::cerr<<"Checking full colour summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,3,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	algo.sum_colours();
	value_type M=algo.evaluate_sum();
	value_type M2=algo.evaluate_colour_sum();
	if(!equals(M,M2))
	{
	    std::cerr<<"Colour sum does not equal dof sum..."<<std::endl;
	    return 1;
	}
	value_type M3(0);
	for(std::size_t a1=0;a1<model_type::N_c;++a1)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->colour(0)=a1;
	    for(std::size_t a2=0;a2<model_type::N_c;++a2)
	    {              
		algo.get_tree_iterator()->get_phase_space(1)->colour(0)=a2;
		for(std::size_t a3=0;a3<model_type::N_c;++a3)
		{              
		    algo.get_tree_iterator()->get_phase_space(2)->colour(0)=a3;
		    for(std::size_t a4=0;a4<model_type::N_c;++a4)
		    {              
			algo.get_tree_iterator()->get_phase_space(3)->colour(0)=a4;
			for(std::size_t a5=0;a5<model_type::N_c;++a5)
			{              
			    algo.get_tree_iterator()->get_phase_space(4)->colour(0)=a5;
			    for(std::size_t a6=0;a6<model_type::N_c;++a6)
			    {
				algo.get_tree_iterator()->get_phase_space(4)->colour(1)=a6;
				M3+=algo.evaluate2();
			    }
			}
		    }
		}
	    }
	}
	M3/=(value_type)2;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    {
	typedef QCD model_type;
	typedef model_type::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,4,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,4,std::random> hel_factory;
	
	std::string process("u,ubar > u,ubar,g,g");
	std::cerr<<"Checking full colour summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,4,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	algo.sum_colours();
	value_type M=algo.evaluate_sum();
	value_type M2=algo.evaluate_colour_sum();
	if(!equals(M,M2))
	{
	    std::cerr<<"Colour sum does not equal dof sum..."<<std::endl;
	    return 1;
	}
	value_type M3(0);
	for(std::size_t a1=0;a1<model_type::N_c;++a1)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->colour(0)=a1;
	    for(std::size_t a2=0;a2<model_type::N_c;++a2)
	    {              
		algo.get_tree_iterator()->get_phase_space(1)->colour(0)=a2;
		for(std::size_t a3=0;a3<model_type::N_c;++a3)
		{              
		    algo.get_tree_iterator()->get_phase_space(2)->colour(0)=a3;
		    for(std::size_t a4=0;a4<model_type::N_c;++a4)
		    {              
			algo.get_tree_iterator()->get_phase_space(3)->colour(0)=a4;
			for(std::size_t a5=0;a5<model_type::N_c;++a5)
			{              
			    algo.get_tree_iterator()->get_phase_space(4)->colour(0)=a5;
			    for(std::size_t a6=0;a6<model_type::N_c;++a6)
			    {
				algo.get_tree_iterator()->get_phase_space(4)->colour(1)=a6;
				for(std::size_t a7=0;a7<model_type::N_c;++a7)
				{
				    algo.get_tree_iterator()->get_phase_space(5)->colour(0)=a7;
				    for(std::size_t a8=0;a8<model_type::N_c;++a8)
				    {
					algo.get_tree_iterator()->get_phase_space(5)->colour(1)=a8;
					M3+=algo.evaluate2();
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	M3/=(value_type)4;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    {
	typedef QCD model_type;
	typedef model_type::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,2,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,2,std::random> hel_factory;
	
	std::string process("u,ubar > g,g");
	std::cerr<<"Checking partial colour summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,2,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	algo.get_tree_iterator()->get_phase_space(0)->colour(0)=1;
	algo.get_tree_iterator()->get_phase_space(1)->colour(0)=2;
	algo.sum_colour(2);
	algo.sum_colour(3);
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(std::size_t a3=0;a3<model_type::N_c;++a3)
	{              
	    algo.get_tree_iterator()->get_phase_space(3)->colour(1)=a3;
	    for(std::size_t a4=0;a4<model_type::N_c;++a4)
	    {              
		algo.get_tree_iterator()->get_phase_space(3)->colour(0)=a4;
		for(std::size_t a5=0;a5<model_type::N_c;++a5)
		{              
		    algo.get_tree_iterator()->get_phase_space(2)->colour(1)=a5;
		    for(std::size_t a6=0;a6<model_type::N_c;++a6)
		    {
			algo.get_tree_iterator()->get_phase_space(2)->colour(0)=a6;
			M3+=algo.evaluate2();
		    }
		}
	    }
	}
	M3/=(value_type)4;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    {
	typedef QCD model_type;
	typedef model_type::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,3,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,3,std::random> hel_factory;
	
	std::string process("u,ubar > u,ubar,g");
	std::cerr<<"Checking partial colour summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,3,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	algo.get_tree_iterator()->get_phase_space(0)->colour(0)=1;
	algo.get_tree_iterator()->get_phase_space(1)->colour(0)=2;
	algo.sum_colour(2);
	algo.sum_colour(3);
	algo.sum_colour(4);
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(std::size_t a3=0;a3<model_type::N_c;++a3)
	{              
	    algo.get_tree_iterator()->get_phase_space(2)->colour(0)=a3;
	    for(std::size_t a4=0;a4<model_type::N_c;++a4)
	    {              
		algo.get_tree_iterator()->get_phase_space(3)->colour(0)=a4;
		for(std::size_t a5=0;a5<model_type::N_c;++a5)
		{              
		    algo.get_tree_iterator()->get_phase_space(4)->colour(0)=a5;
		    for(std::size_t a6=0;a6<model_type::N_c;++a6)
		    {
			algo.get_tree_iterator()->get_phase_space(4)->colour(1)=a6;
			M3+=algo.evaluate2();
		    }
		}
	    }
	}
	M3/=(value_type)2;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    {
	typedef QCD model_type;
	typedef model_type::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,4,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,4,std::random> hel_factory;
	
	std::string process("u,ubar > u,ubar,g,g");
	std::cerr<<"Checking partial colour summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,4,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	algo.get_tree_iterator()->get_phase_space(0)->colour(0)=1;
	algo.get_tree_iterator()->get_phase_space(1)->colour(0)=2;
	algo.sum_colour(2);
	algo.sum_colour(3);
	algo.sum_colour(4);
	algo.sum_colour(5);
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(std::size_t a3=0;a3<model_type::N_c;++a3)
	{              
	    algo.get_tree_iterator()->get_phase_space(2)->colour(0)=a3;
	    for(std::size_t a4=0;a4<model_type::N_c;++a4)
	    {              
		algo.get_tree_iterator()->get_phase_space(3)->colour(0)=a4;
		for(std::size_t a5=0;a5<model_type::N_c;++a5)
		{              
		    algo.get_tree_iterator()->get_phase_space(4)->colour(0)=a5;
		    for(std::size_t a6=0;a6<model_type::N_c;++a6)
		    {
			algo.get_tree_iterator()->get_phase_space(4)->colour(1)=a6;
			for(std::size_t a7=0;a7<model_type::N_c;++a7)
			{
			    algo.get_tree_iterator()->get_phase_space(5)->colour(0)=a7;
			    for(std::size_t a8=0;a8<model_type::N_c;++a8)
			    {
				algo.get_tree_iterator()->get_phase_space(5)->colour(1)=a8;
				M3+=algo.evaluate2();
			    }
			}
		    }
		}
	    }
	}
	M3/=(value_type)4;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    {
	typedef QCD model_type;
	typedef model_type::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,2,std::random> ps_factory;
	
	std::string process("u,ubar > g,g");
	std::cerr<<"Checking full colour+helicity summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	algo.sum_colours();
	algo.sum_spins();
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(std::size_t a1=0;a1<model_type::N_c;++a1)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->colour(0)=a1;
	    for(int h1=-1;h1<2;h1+=2)
	    {
		algo.get_tree_iterator()->get_phase_space(0)->helicity_phase(h1)=std::complex<value_type>(1,0);
		algo.get_tree_iterator()->get_phase_space(0)->helicity_phase(-h1)=std::complex<value_type>(0,0);
		for(std::size_t a2=0;a2<model_type::N_c;++a2)
		{		
		    algo.get_tree_iterator()->get_phase_space(1)->colour(0)=a2;
		    for(int h2=-1;h2<2;h2+=2)
		    {
			algo.get_tree_iterator()->get_phase_space(1)->helicity_phase(h2)=std::complex<value_type>(1,0);
			algo.get_tree_iterator()->get_phase_space(1)->helicity_phase(-h2)=std::complex<value_type>(0,0);
			for(std::size_t a3=0;a3<model_type::N_c;++a3)
			{		    
			    algo.get_tree_iterator()->get_phase_space(2)->colour(0)=a3;
			    for(int h3=-1;h3<2;h3+=2)
			    {
				algo.get_tree_iterator()->get_phase_space(2)->helicity_phase(h3)=std::complex<value_type>(1,0);
				algo.get_tree_iterator()->get_phase_space(2)->helicity_phase(-h3)=std::complex<value_type>(0,0);
				for(std::size_t a4=0;a4<model_type::N_c;++a4)
				{              
				    algo.get_tree_iterator()->get_phase_space(2)->colour(1)=a4;
				    for(int h4=-1;h4<2;h4+=2)
				    {
					algo.get_tree_iterator()->get_phase_space(3)->helicity_phase(h4)=std::complex<value_type>(1,0);
					algo.get_tree_iterator()->get_phase_space(3)->helicity_phase(-h4)=std::complex<value_type>(0,0);
					for(std::size_t a5=0;a5<model_type::N_c;++a5)
					{              
					    algo.get_tree_iterator()->get_phase_space(3)->colour(0)=a5;
					    for(std::size_t a6=0;a6<model_type::N_c;++a6)
					    {
						algo.get_tree_iterator()->get_phase_space(3)->colour(1)=a6;
						M3+=algo.evaluate2();
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	M3/=(value_type)4;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    {
	typedef QCD model_type;
	typedef model_type::value_type value_type;
	typedef uniform_ps_generator_factory<model_type,2,3,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,3,std::random> hel_factory;
	
	std::string process("u,ubar > u,ubar,g");
	std::cerr<<"Checking partial colour+helicity summation in QCD for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	algo.sum_colour(0);
	algo.sum_colour(1);
	algo.get_tree_iterator()->get_phase_space(2)->colour(0)=1;
	algo.get_tree_iterator()->get_phase_space(3)->colour(0)=2;
	algo.sum_colour(4);
	value_type h=std::sqrt((value_type)0.5);
	algo.get_tree_iterator()->get_phase_space(0)->helicity_phase(-1)=h;
	algo.get_tree_iterator()->get_phase_space(0)->helicity_phase(1)=-h;
	algo.get_tree_iterator()->get_phase_space(1)->helicity_phase(-1)=h;
	algo.get_tree_iterator()->get_phase_space(1)->helicity_phase(1)=h;
	algo.sum_spin(2);
	algo.sum_spin(3);
	algo.sum_spin(4);
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(std::size_t a1=0;a1<model_type::N_c;++a1)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->colour(0)=a1;
	    for(std::size_t a2=0;a2<model_type::N_c;++a2)
	    {		
		algo.get_tree_iterator()->get_phase_space(1)->colour(0)=a2;
		for(int h3=-1;h3<2;h3+=2)
		{
		    algo.get_tree_iterator()->get_phase_space(2)->helicity_phase(h3)=std::complex<value_type>(1,0);
		    algo.get_tree_iterator()->get_phase_space(2)->helicity_phase(-h3)=std::complex<value_type>(0,0);
		    for(int h4=-1;h4<2;h4+=2)
		    {
			algo.get_tree_iterator()->get_phase_space(3)->helicity_phase(h4)=std::complex<value_type>(1,0);
			algo.get_tree_iterator()->get_phase_space(3)->helicity_phase(-h4)=std::complex<value_type>(0,0);
			for(std::size_t a5=0;a5<model_type::N_c;++a5)
			{              
			    algo.get_tree_iterator()->get_phase_space(4)->colour(0)=a5;
			    for(std::size_t a6=0;a6<model_type::N_c;++a6)
			    {
				algo.get_tree_iterator()->get_phase_space(4)->colour(1)=a6;
				for(int h6=-1;h6<2;h6+=2)
				{
				    algo.get_tree_iterator()->get_phase_space(4)->helicity_phase(h6)=std::complex<value_type>(1,0);
				    algo.get_tree_iterator()->get_phase_space(4)->helicity_phase(-h6)=std::complex<value_type>(0,0);
				    M3+=algo.evaluate2();
				}
			    }
			}
		    }
		}
	    }
	}
	M3/=(value_type)2;
	if(!equals(M,M3))
	{
	    std::cerr<<"Colour sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing colour sampling algorithms......................................."<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    
    bool have_gp=plot_config::gnuplot_path!=NULL;
    
    {
	typedef QCDPbchabdc model_type;
	typedef QCD model_type2;
	typedef uniform_ps_generator_factory<model_type,2,2,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,2,std::random> hel_factory;
	
	std::string process("g,g > dbar,d");
	std::string filename("plots/colMC1");
	std::cerr<<"Checking adjoint vs. colour flow generation in QCD for "<<process<<"..........";
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,2,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	CM_algorithm<model_type2,2,2>algo2(process);
	algo2.load();
	algo2.construct();
	copy_momenta<model_type,model_type2,2,2>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	copy_helicities<model_type,model_type2,2,2>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	value_type summed_result=algo.evaluate_colour_sum();
	value_type summed_result2=algo2.evaluate_colour_sum();
	MC_integrator<value_type>* colgen1=new MC_generator_wrapper<value_type>(adjoint_QCD<value_type,2,2,3,std::random,false>::create_instance<model_type>(algo.get_tree_iterator()));
	MC_integrator<value_type>* colgen2=new MC_generator_wrapper<value_type>(colour_flow_QCD<value_type,2,2,3,std::random,false>::create_instance<model_type2>(algo2.get_tree_iterator()));
	value_type n,x1,y1,x2,y2;
	data_wrapper* data=have_gp?(new data_wrapper(&n,&x1,&y1)):(new data_wrapper(filename+".dat",&n,&x1,&y1));
	data->add_leaf(&x2);
	data->add_leaf(&y2);
	for(std::size_t i=0;i<N_events;++i)
	{
	    colgen1->generate();
	    colgen1->integrand()=algo.evaluate2();
	    colgen1->refresh_cross_section();
	    colgen2->generate();
	    colgen2->integrand()=algo2.evaluate2();
	    colgen2->refresh_cross_section();
	    if(i!=0 and i%N_batch==0)
	    {
		n=(value_type)i;
		x1=colgen1->cross_section().value;
		y1=colgen1->cross_section().error;
		x2=colgen2->cross_section().value;
		y2=colgen2->cross_section().error;
		data->fill();
	    }
	}
	data->write();
	data_stream* datastr1=new data_stream(data,"1","2","3");
	datastr1->title="adjoint MC sum";
	datastr1->style="yerrorbars";
	data_stream* datastr2=new data_stream(data,"1","4","5");
	datastr2->title="col flow MC sum";
	datastr2->style="yerrorbars";
	plot_script* graph=new plot_script("MC over discrete colours","postscript color");
	graph->file=filename;
	graph->add_y_tic(summed_result,"adjoint sum");
	graph->add_y_tic(summed_result2,"col flow sum");
	graph->grid=true;
	graph->add_plot(datastr1);
	graph->add_plot(datastr2);
	graph->plot();
	delete colgen1;
	delete colgen2;
	delete psgen;
	delete graph;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    {
	typedef QCDPbchabdc model_type;
	typedef QCDPbchabcc model_type2;
	typedef uniform_ps_generator_factory<model_type,2,3,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,3,std::random> hel_factory;
	
	std::string process("u,ubar > d,dbar,g");
	std::string filename("plots/colMC2");
	std::cerr<<"Checking adjoint colour generation in QCD for "<<process<<"..........";
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,3,model_type::continuous_helicities>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	CM_algorithm<model_type2,2,3>algo2(process);
	algo2.load();
	algo2.construct();
	copy_momenta<model_type,model_type2,2,3>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	copy_helicities<model_type,model_type2,2,3>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	value_type summed_result=algo.evaluate_colour_sum();
	MC_integrator<value_type>* colgen1=new MC_generator_wrapper<value_type>(adjoint_QCD<value_type,2,3,3,std::random,false>::create_instance<model_type>(algo.get_tree_iterator()));
	MC_integrator<value_type>* colgen2=new MC_generator_wrapper<value_type>(adjoint_QCD<value_type,2,3,3,std::random,true>::create_instance<model_type2>(algo2.get_tree_iterator()));
	value_type n,x1,y1,x2,y2;
	data_wrapper* data=have_gp?(new data_wrapper(&n,&x1,&y1)):(new data_wrapper(filename+".dat",&n,&x1,&y1));
	data->add_leaf(&x2);
	data->add_leaf(&y2);
	for(std::size_t i=0;i<N_events;++i)
	{
	    colgen1->generate();
	    colgen1->integrand()=algo.evaluate2();
	    colgen1->refresh_cross_section();
	    colgen2->generate();
	    colgen2->integrand()=algo2.evaluate2();
	    colgen2->refresh_cross_section();
	    if(i!=0 and i%N_batch==0)
	    {
		n=(value_type)i;
		x1=colgen1->cross_section().value;
		y1=colgen1->cross_section().error;
		x2=colgen2->cross_section().value;
		y2=colgen2->cross_section().error;
		data->fill();
	    }
	}
	data->write();
	data_stream* datastr1=new data_stream(data,"1","2","3");
	datastr1->title="discr MC sum";
	datastr1->style="yerrorbars";
	data_stream* datastr2=new data_stream(data,"1","4","5");
	datastr2->title="cont MC sum";
	datastr2->style="yerrorbars";
	plot_script* graph=new plot_script("Adjoint rep MC","postscript color");
	graph->file=filename;
	graph->add_y_tic(summed_result,"sum");
	graph->grid=true;
	graph->add_plot(datastr1);
	graph->add_plot(datastr2);
	graph->plot();
	delete colgen1;
	delete colgen2;
	delete psgen;
	delete graph;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    {
	typedef QCD model_type;
	typedef QCDPbchcfcc model_type2;
	typedef uniform_ps_generator_factory<model_type,2,4,std::random> ps_factory;
	typedef helicity_generator_factory<model_type,2,4,std::random> hel_factory;
	
	std::string process("g,g > d,dbar,u,ubar");
	std::string filename("plots/colMC3");
	std::cerr<<"Checking colour-flow colour generation in QCD for "<<process<<"...........";
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* psgen=ps_factory::create_generator(algo.get_tree_iterator());
	psgen->generate();
	helicity_generator<value_type,2,4,true>* helgen=hel_factory::create_generator(algo.get_tree_iterator());
	helgen->generate();
	algo.evaluate();
	CM_algorithm<model_type2,2,4>algo2(process);
	algo2.load();
	algo2.construct();
	copy_momenta<model_type,model_type2,2,4>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	copy_helicities<model_type,model_type2,2,4>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	model_type::value_type summed_result=algo.evaluate_colour_sum();
	MC_integrator<value_type>* colgen0=new MC_generator_wrapper<value_type>(uniform_colours<value_type,2,4,std::random,false>::create_instance<model_type>(algo.get_tree_iterator()));
	MC_integrator<value_type>* colgen1=new MC_generator_wrapper<value_type>(colour_flow_QCD<value_type,2,4,3,std::random,false>::create_instance<model_type>(algo.get_tree_iterator()));
	MC_integrator<value_type>* colgen2=new MC_generator_wrapper<value_type>(colour_flow_QCD<value_type,2,4,3,std::random,true>::create_instance<model_type2>(algo2.get_tree_iterator()));
	value_type n,x1,y1,x2,y2,x3,y3;
	data_wrapper* data=have_gp?(new data_wrapper(&n,&x1,&y1)):(new data_wrapper(filename+".dat",&n,&x1,&y1));
	data->add_leaf(&x2);
	data->add_leaf(&y2);
	data->add_leaf(&x3);
	data->add_leaf(&y3);
	for(std::size_t i=0;i<N_events;++i)
	{
	    colgen0->generate();
	    colgen0->integrand()=algo.evaluate2();
	    colgen0->refresh_cross_section();
	    colgen1->generate();
	    colgen1->integrand()=algo.evaluate2();
	    colgen1->refresh_cross_section();
	    colgen2->generate();
	    colgen2->integrand()=algo2.evaluate2();
	    colgen2->refresh_cross_section();
	    if(i!=0 and i%N_batch==0)
	    {
		n=(value_type)i;
		x1=(value_type)0.25*colgen0->cross_section().value;
		y1=(value_type)0.25*colgen0->cross_section().error;
		x2=colgen1->cross_section().value;
		y2=colgen1->cross_section().error;
		x3=colgen2->cross_section().value;
		y3=colgen2->cross_section().error;
		data->fill();
	    }
	}
	data->write();
	data_stream* datastr1=new data_stream(data,"1","2","3");
	datastr1->title="uni discr MC sum";
	datastr1->style="yerrorbars";
	data_stream* datastr2=new data_stream(data,"1","4","5");
	datastr2->title="col-flow discr MC sum";
	datastr2->style="yerrorbars";
	data_stream* datastr3=new data_stream(data,"1","6","7");
	datastr3->title="col-flow cont MC sum";
	datastr3->style="yerrorbars";
	plot_script* graph=new plot_script("Colour-flow rep MC","postscript color");
	graph->file=filename;
	graph->add_y_tic(summed_result,"sum");
	graph->grid=true;
	graph->add_plot(datastr1);
	graph->add_plot(datastr2);
	graph->add_plot(datastr3);
	graph->plot();
	delete colgen0;
	delete colgen1;
	delete colgen2;
	delete helgen;
	delete psgen;
	delete graph;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
}

