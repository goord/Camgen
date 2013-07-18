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
#include <Camgen/uni_hels.h>
#include <Camgen/part_is.h>
#include <Camgen/rambo.h>
#include <Camgen/plt_config.h>
#include <Camgen/plt_script.h>
#include <Camgen/SM.h>
#include <QEDPbdh.h>
#include <QEDPbch.h>

using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing helicity summation algorithms...................................."<<std::endl;
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
   
    {
	typedef QEDPbdh model_type;
	typedef rambo<model_type,2,2,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-");
	std::cerr<<"Checking full helicity summation in QED for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	algo.sum_spins();
	value_type M=algo.evaluate_sum();
	value_type M2=algo.evaluate_spin_sum();
	if(!equals(M,M2))
	{
	    std::cerr<<"Helicity sum does not equal dof sum..."<<std::endl;
	    return 1;
	}
	value_type M3(0);
	for(int h1=-1;h1<2;h1+=2)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->helicity()=h1;
	    for(int h2=-1;h2<2;h2+=2)
	    {
		algo.get_tree_iterator()->get_phase_space(1)->helicity()=h2;
		for(int h3=-1;h3<2;h3+=2)
		{
		    algo.get_tree_iterator()->get_phase_space(2)->helicity()=h3;
		    for(int h4=-1;h4<2;h4+=2)
		    {
			algo.get_tree_iterator()->get_phase_space(3)->helicity()=h4;
			M3+=algo.evaluate2();
		    }
		}
	    }
	}
	if(!equals(M,M3))
	{
	    std::cerr<<"Helicity sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }

    {
	typedef QEDPbdh model_type;
	typedef rambo<model_type,2,3,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-,gamma");
	std::cerr<<"Checking full helicity summation in QED for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	algo.sum_spins();
	value_type M=algo.evaluate_sum();
	value_type M2=algo.evaluate_spin_sum();
	if(!equals(M,M2))
	{
	    std::cerr<<"Helicity sum does not equal dof sum..."<<std::endl;
	    return 1;
	}
	value_type M3(0);
	for(int h1=-1;h1<2;h1+=2)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->helicity()=h1;
	    for(int h2=-1;h2<2;h2+=2)
	    {
		algo.get_tree_iterator()->get_phase_space(1)->helicity()=h2;
		for(int h3=-1;h3<2;h3+=2)
		{
		    algo.get_tree_iterator()->get_phase_space(2)->helicity()=h3;
		    for(int h4=-1;h4<2;h4+=2)
		    {
			algo.get_tree_iterator()->get_phase_space(3)->helicity()=h4;
			for(int h5=-1;h5<2;h5+=2)
			{
			    algo.get_tree_iterator()->get_phase_space(4)->helicity()=h5;
			    M3+=algo.evaluate2();
			}
		    }
		}
	    }
	}
	if(!equals(M,M3))
	{
	    std::cerr<<"Helicity sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }

    {
	typedef QEDPbdh model_type;
	typedef rambo<model_type,2,4,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,gamma,e-,gamma");
	std::cerr<<"Checking full helicity summation in QED for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	algo.sum_spins();
	value_type M=algo.evaluate_sum();
	value_type M2=algo.evaluate_spin_sum();
	if(!equals(M,M2))
	{
	    std::cerr<<"Helicity sum does not equal dof sum..."<<std::endl;
	    return 1;
	}
	value_type M3(0);
	for(int h1=-1;h1<2;h1+=2)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->helicity()=h1;
	    for(int h2=-1;h2<2;h2+=2)
	    {
		algo.get_tree_iterator()->get_phase_space(1)->helicity()=h2;
		for(int h3=-1;h3<2;h3+=2)
		{
		    algo.get_tree_iterator()->get_phase_space(2)->helicity()=h3;
		    for(int h4=-1;h4<2;h4+=2)
		    {
			algo.get_tree_iterator()->get_phase_space(3)->helicity()=h4;
			for(int h5=-1;h5<2;h5+=2)
			{
			    algo.get_tree_iterator()->get_phase_space(4)->helicity()=h5;
			    for(int h6=-1;h6<2;h6+=2)
			    {
				algo.get_tree_iterator()->get_phase_space(5)->helicity()=h6;
				M3+=algo.evaluate2();
			    }
			}
		    }
		}
	    }
	}
	if(!equals(M,M3))
	{
	    std::cerr<<"Helicity sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }

    {
	typedef QEDPbdh model_type;
	typedef rambo<model_type,2,2,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-");
	std::cerr<<"Checking partial helicity summation in QED for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	algo.sum_spin(2);
	algo.sum_spin(3);
	algo.get_tree_iterator()->get_phase_space(0)->helicity()=1;
	algo.get_tree_iterator()->get_phase_space(1)->helicity()=-1;
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(int h3=-1;h3<2;h3+=2)
	{
	    algo.get_tree_iterator()->get_phase_space(2)->helicity()=h3;
	    for(int h4=-1;h4<2;h4+=2)
	    {
		algo.get_tree_iterator()->get_phase_space(3)->helicity()=h4;
		M3+=algo.evaluate2();
	    }
	}
	if(!equals(M,M3))
	{
	    std::cerr<<"Helicity sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }

    {
	typedef QEDPbdh model_type;
	typedef rambo<model_type,2,3,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-,gamma");
	std::cerr<<"Checking partial helicity summation in QED for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	algo.sum_spin(2);
	algo.sum_spin(3);
	algo.get_tree_iterator()->get_phase_space(0)->helicity()=1;
	algo.get_tree_iterator()->get_phase_space(1)->helicity()=-1;
	algo.get_tree_iterator()->get_phase_space(4)->helicity()=-1;
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(int h3=-1;h3<2;h3+=2)
	{
	    algo.get_tree_iterator()->get_phase_space(2)->helicity()=h3;
	    for(int h4=-1;h4<2;h4+=2)
	    {
		algo.get_tree_iterator()->get_phase_space(3)->helicity()=h4;
		M3+=algo.evaluate2();
	    }
	}
	if(!equals(M,M3))
	{
	    std::cerr<<"Helicity sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }

    {
	typedef QEDPbdh model_type;
	typedef rambo<model_type,2,4,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-,gamma,gamma");
	std::cerr<<"Checking partial helicity summation in QED for "<<process<<"..........";
	std::cerr.flush();
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	algo.sum_spin(0);
	algo.sum_spin(1);
	algo.sum_spin(2);
	algo.sum_spin(4);
	algo.get_tree_iterator()->get_phase_space(3)->helicity()=1;
	algo.get_tree_iterator()->get_phase_space(5)->helicity()=-1;
	value_type M=algo.evaluate_sum();
	value_type M3(0);
	for(int h1=-1;h1<2;h1+=2)
	{
	    algo.get_tree_iterator()->get_phase_space(0)->helicity()=h1;
	    for(int h2=-1;h2<2;h2+=2)
	    {
		algo.get_tree_iterator()->get_phase_space(1)->helicity()=h2;
		for(int h3=-1;h3<2;h3+=2)
		{
		    algo.get_tree_iterator()->get_phase_space(2)->helicity()=h3;
		    for(int h5=-1;h5<2;h5+=2)
		    {
			algo.get_tree_iterator()->get_phase_space(4)->helicity()=h5;
			M3+=algo.evaluate2();
		    }
		}
	    }
	}
	if(!equals(M,M3))
	{
	    std::cerr<<"Helicity sum "<<M<<" does not equal explicit sum "<<M3<<std::endl;
	    return 1;
	}
	delete psgen;
	std::cerr<<"..........done."<<std::endl;
    }
    
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing helicity sampling algorithms......................................"<<std::endl;
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
    
    bool have_gp=plot_config::gnuplot_path!=NULL;
    
    {
	typedef QEDPbdh model_type;
	typedef QEDPbch model_type2;
	typedef rambo<model_type,2,2,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-");
	std::string filename("plots/helMC1");
	std::cerr<<"Checking uniform helicity generation in QED for "<<process<<"...........";
	CM_algorithm<model_type,2,2>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,2>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	CM_algorithm<model_type2,2,2>algo2(process);
	algo2.load();
	algo2.construct();
	copy_momenta<model_type,model_type2,2,2>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	value_type summed_result=algo.evaluate_spin_sum();
	uniform_helicities<value_type,2,2,std::random,false>* helgen1=uniform_helicities<value_type,2,2,std::random,false>::create_instance<model_type>(algo.get_tree_iterator());
	uniform_helicities<value_type,2,2,std::random,true>* helgen2=uniform_helicities<value_type,2,2,std::random,true>::create_instance<model_type2>(algo2.get_tree_iterator());
	value_type n,x1,y1,x2,y2;
	data_wrapper* data=have_gp?(new data_wrapper(&n,&x1,&y1)):(new data_wrapper(filename+".dat",&n,&x1,&y1));
	data->add_leaf(&x2);
	data->add_leaf(&y2);
	for(std::size_t i=0;i<N_events;++i)
	{
	    helgen1->generate();
	    helgen1->integrand()=algo.evaluate2();
	    helgen1->refresh_cross_section();
	    helgen2->generate();
	    helgen2->integrand()=algo2.evaluate2();
	    helgen2->refresh_cross_section();
	    if(i!=0 and i%N_batch==0)
	    {
		n=(value_type)i;
		x1=helgen1->cross_section().value;
		y1=helgen1->cross_section().error;
		x2=helgen2->cross_section().value;
		y2=helgen2->cross_section().error;
		data->fill();
	    }
	}
	data->write();
	data_stream* datastr1=new data_stream(data,"1","2","3");
	datastr1->title="discrete MC";
	datastr1->style="yerrorbars";
	data_stream* datastr2=new data_stream(data,"1","4","5");
	datastr2->title="continuous MC";
	datastr2->style="yerrorbars";
	plot_script* graph=new plot_script("Monte Carlo over helicities","postscript color");
	graph->file=filename;
	graph->add_y_tic(summed_result,"sum");
	graph->grid=true;
	graph->add_plot(datastr1);
	graph->add_plot(datastr2);
	graph->plot();
	delete helgen1;
	delete helgen2;
	delete psgen;
	delete graph;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    {
	typedef QEDPbdh model_type;
	typedef QEDPbch model_type2;
	typedef rambo<model_type,2,3,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-,gamma");
	std::string filename("plots/helMC2");
	std::cerr<<"Checking uniform helicity generation in QED for "<<process<<"...........";
	CM_algorithm<model_type,2,3>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,3>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	CM_algorithm<model_type2,2,3>algo2(process);
	algo2.load();
	algo2.construct();
	copy_momenta<model_type,model_type2,2,3>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	value_type summed_result=algo.evaluate_spin_sum();
	uniform_helicities<value_type,2,3,std::random,false>* helgen1=uniform_helicities<value_type,2,3,std::random,false>::create_instance<model_type>(algo.get_tree_iterator());
	uniform_helicities<value_type,2,3,std::random,true>* helgen2=uniform_helicities<value_type,2,3,std::random,true>::create_instance<model_type2>(algo2.get_tree_iterator());
	value_type n,x1,y1,x2,y2;
	data_wrapper* data=have_gp?(new data_wrapper(&n,&x1,&y1)):(new data_wrapper(filename+".dat",&n,&x1,&y1));
	data->add_leaf(&x2);
	data->add_leaf(&y2);
	for(std::size_t i=0;i<N_events;++i)
	{
	    helgen1->generate();
	    helgen1->integrand()=algo.evaluate2();
	    helgen1->refresh_cross_section();
	    helgen2->generate();
	    helgen2->integrand()=algo2.evaluate2();
	    helgen2->refresh_cross_section();
	    if(i!=0 and i%N_batch==0)
	    {
		n=(value_type)i;
		x1=helgen1->cross_section().value;
		y1=helgen1->cross_section().error;
		x2=helgen2->cross_section().value;
		y2=helgen2->cross_section().error;
		data->fill();
	    }
	}
	data->write();
	data_stream* datastr1=new data_stream(data,"1","2","3");
	datastr1->title="discrete MC";
	datastr1->style="yerrorbars";
	data_stream* datastr2=new data_stream(data,"1","4","5");
	datastr2->title="continuous MC";
	datastr2->style="yerrorbars";
	plot_script* graph=new plot_script("Monte Carlo over helicities","postscript color");
	graph->file=filename;
	graph->add_y_tic(summed_result,"sum");
	graph->grid=true;
	graph->add_plot(datastr1);
	graph->add_plot(datastr2);
	graph->plot();
	delete helgen1;
	delete helgen2;
	delete psgen;
	delete graph;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }

    {
	typedef QEDPbdh model_type;
	typedef QEDPbch model_type2;
	typedef rambo<model_type,2,4,std::random> ps_gen;
	typedef partonic_is<model_type,2> init_state;
	
	std::string process("e+,e- > e+,e-,gamma,gamma");
	std::string filename("plots/helMC3");
	std::cerr<<"Checking uniform helicity generation in QED for "<<process<<"...........";
	CM_algorithm<model_type,2,4>algo(process);
	algo.load();
	algo.construct();
	ps_generator<model_type,2,4>* psgen=ps_gen::create_instance(algo.get_tree_iterator(),new init_state(Ecm,Ecm));
	psgen->generate();
	CM_algorithm<model_type2,2,4>algo2(process);
	algo2.load();
	algo2.construct();
	copy_momenta<model_type,model_type2,2,4>(algo.get_tree_iterator(),algo2.get_tree_iterator());
	value_type summed_result=algo.evaluate_spin_sum();
	uniform_helicities<value_type,2,4,std::random,false>* helgen1=uniform_helicities<value_type,2,4,std::random,false>::create_instance<model_type>(algo.get_tree_iterator());
	uniform_helicities<value_type,2,4,std::random,true>* helgen2=uniform_helicities<value_type,2,4,std::random,true>::create_instance<model_type2>(algo2.get_tree_iterator());
	value_type n,x1,y1,x2,y2;
	data_wrapper* data=have_gp?(new data_wrapper(&n,&x1,&y1)):(new data_wrapper(filename+".dat",&n,&x1,&y1));
	data->add_leaf(&x2);
	data->add_leaf(&y2);
	for(std::size_t i=0;i<N_events;++i)
	{
	    helgen1->generate();
	    helgen1->integrand()=algo.evaluate2();
	    helgen1->refresh_cross_section();
	    helgen2->generate();
	    helgen2->integrand()=algo2.evaluate2();
	    helgen2->refresh_cross_section();
	    if(i!=0 and i%N_batch==0)
	    {
		n=(value_type)i;
		x1=helgen1->cross_section().value;
		y1=helgen1->cross_section().error;
		x2=helgen2->cross_section().value;
		y2=helgen2->cross_section().error;
		data->fill();
	    }
	}
	data->write();
	data_stream* datastr1=new data_stream(data,"1","2","3");
	datastr1->title="discrete MC";
	datastr1->style="yerrorbars";
	data_stream* datastr2=new data_stream(data,"1","4","5");
	datastr2->title="continuous MC";
	datastr2->style="yerrorbars";
	plot_script* graph=new plot_script("Monte Carlo over helicities","postscript color");
	graph->file=filename;
	graph->add_y_tic(summed_result,"sum");
	graph->grid=true;
	graph->add_plot(datastr1);
	graph->add_plot(datastr2);
	graph->plot();
	delete helgen1;
	delete helgen2;
	delete psgen;
	delete graph;
	std::string output=have_gp?(filename+".eps"):(filename+".dat/gp");
	std::cerr<<"............done, "<<output<<" written."<<std::endl;
    }
}

