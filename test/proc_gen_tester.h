//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef PROC_GEN_TESTER_H_
#define PROC_GEN_TESTER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing class template for process generators. Creates an equivalent  *
 * uniform Monte Carlo generator and compares the cross sections within  *
 * error bounds.                                                         *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/procgen_fac.h>
#include <ps_gen_tester.h>

namespace Camgen
{
    // Class that instantiates a process generator and a corresponding uniform generator. Upon executing the 'run'
    // function a comparison is carried out to check the cross section.

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class process_generator_tester
    {
	public:

	    typedef process_generator<model_t,N_in,N_out,rng_t> process_generator_type;
	    typedef process_generator_factory<model_t,N_in,N_out,rng_t> process_generator_factory_type;
	    typedef typename process_generator_factory_type::CM_tree_iterator CM_tree_iterator;
	    typedef typename model_t::value_type value_type;
	    typedef std::size_t size_type;

	    /* Output file name: */

	    const std::string filename;

	    // Factory creating the process generators.

	    process_generator_factory_type* proc_gen_fac;

	    // Process generator to be tested.

	    process_generator_type* proc_gen;
	    
	    // Uniform process generator te be tested against.
	    
	    process_generator_type* uni_gen;

	    // Constructor. Takes an amplitude subprocess iterator as argument.
	    
	    process_generator_tester(CM_tree_iterator amplitude, const std::string& fname):filename(fname),proc_gen_fac(new process_generator_factory_type())
	    {
		proc_gen=proc_gen_fac->create_generator(amplitude);
		proc_gen->set_auto_update(true);
		
		phase_space_generators::type psgen_type=phase_space_generator_type();
		set_phase_space_generator_type(phase_space_generators::uniform);
		
		uni_gen=proc_gen_fac->create_generator(amplitude);
		uni_gen->set_auto_update(true);
		
		set_phase_space_generator_type(psgen_type);
	    }

	    // Sets the i-th incoming particle energy

	    bool set_beam_energy(int i,const value_type& E)
	    {
		return (proc_gen->set_beam_energy(i,E) and uni_gen->set_beam_energy(i,E));
	    }

	    // Runs the test.

	    bool run(size_type N_events,size_type N_points,bool verbose=false,bool skip_plot=false,bool skip_check=false)
	    {
		size_type N_batch=N_events/N_points;
		value_type n,xsec1,xsec1err,xsec2,xsec2err;
		bool have_gp=plot_config::gnuplot_path!=NULL;
		data_wrapper* data=have_gp?new data_wrapper(&n,&xsec1,&xsec1err):new data_wrapper(filename+".dat",&n,&xsec1,&xsec1err);
		data->add_leaf(&xsec2);
		data->add_leaf(&xsec2err);

		proc_gen->refresh_Ecm();
		uni_gen->refresh_Ecm();
		proc_gen->refresh_m_min();
		uni_gen->refresh_m_min();

		for(size_type i=0;i<N_events;++i)
		{
		    bool ps_generated=proc_gen->generate();

		    if(ps_generated)
		    {
			if(proc_gen->integrand()==0)
			{
			    proc_gen->print_amplitude(std::cerr);
			}
			if(!skip_check and !proc_gen->get_momentum_generator()->check())
			{
			    proc_gen->get_momentum_generator()->print(std::cerr);
			    return false;
			}
			if(proc_gen->weight()!=proc_gen->weight())
			{
			    Camgen::log<<"Invalid phase space weight encountered: "<<proc_gen->weight()<<endlog;
			    return false;
			}
			if(proc_gen->integrand()!=proc_gen->integrand())
			{
			    Camgen::log<<"Invalid phase space integrand encountered: "<<proc_gen->integrand()<<endlog;
			    return false;
			}
		    }
		    else
		    {
//			proc_gen->get_momentum_generator()->print(std::cerr);
		    }
		    if(verbose)
		    {
			std::cerr<<i<<": w = "<<proc_gen->weight()<<", M ="<<proc_gen->integrand()<<std::endl;
		    }
		    proc_gen->refresh_cross_section();
		    if(proc_gen->weight()!=(value_type)0)
		    {
			value_type w=proc_gen->weight();
			proc_gen->evaluate_weight();
			value_type w_check=proc_gen->weight();
			if(!equals(w_check,w,(value_type)1,(value_type)0.1))
			{
			    Camgen::log<<"Weight re-evaluation lead to different result: "<<w<<" not equal to  "<<w_check<<endlog;
			}
		    }
		    
		    bool rambo_generated=uni_gen->generate();
		    if(rambo_generated)
		    {
			if(!skip_check and !uni_gen->get_momentum_generator()->check())
			{
			    return false;
			}
		    }
		    uni_gen->refresh_cross_section();

		    if(i!=0 and i%N_batch==0)
		    {
			n=(value_type)i;
			value_type n2=n*(n-(value_type)1);
			xsec1=(proc_gen->cross_section()).value;
			xsec1err=(proc_gen->cross_section()).error;
			xsec2=(uni_gen->cross_section()).value;
			xsec2err=(uni_gen->cross_section()).error;
			data->fill();
		    }
		}

		if(!skip_plot)
		{
		    data->write();
		    data_stream* datstr1=new data_stream(data,"1","2","3");
		    datstr1->title="MC";
		    datstr1->style="yerrorbars";
		    data_stream* datstr2=new data_stream(data,"1","4","5");
		    datstr2->title="RAMBO";
		    datstr2->style="yerrorbars";
		    plot_script* plot=new plot_script("Cross section","postscript enhanced color");
		    plot->file=filename;
		    plot->add_plot(datstr1);
		    plot->add_plot(datstr2);
		    double dx=N_events/5;
		    for(size_type i=0;i<5;++i)
		    {
			plot->add_x_tic(i*dx);
		    }
		    plot->plot();
		    delete plot;
		}

		return equals(proc_gen->cross_section(),uni_gen->cross_section());
	    }

	    // Destructor.

	    ~process_generator_tester()
	    {
		delete proc_gen_fac;
		delete proc_gen;
		delete uni_gen;
	    }
    };

}

#endif /*PROC_GEN_TESTER_H_*/

