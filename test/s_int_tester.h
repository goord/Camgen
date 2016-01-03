//
// This file is part of the CAMORRA library.
// Copyright (C) 2010 Gijs van den Oord.
// CAMORRA is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef S_INT_TESTER_H_
#define S_INT_TESTER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing class for phase space integrals of constrained products of pairs of   *
 * invariant mass densities. The facility tests the hard-coded integrals against *
 * a Monte-Carlo run with unconstrained generated s-pairs.                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/multi_plot.h>
#include <Camgen/ss_gen.h>
#include <Camgen/Dirac_delta.h>
#include <Camgen/uni_val_gen.h>
#include <Camgen/MC_int.h>

namespace Camgen
{
    template<class value_t,class rng_t>class s_int_tester
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef value_generator<value_t,rng_t> s_generator_type;
	    typedef Dirac_delta<value_t,rng_t> Dirac_delta_type;
	    typedef uniform_value_generator<value_t,rng_t> uniform_type;

	    /* Constructor: */

	    s_int_tester(s_generator_type* sgen1_,s_generator_type* sgen2_):sgen1(sgen1_),sgen2(sgen2_),s1_on_shell(dynamic_cast<Dirac_delta_type*>(sgen1_)!=NULL),s2_on_shell(dynamic_cast<Dirac_delta_type*>(sgen2_)!=NULL)
	    {
		if(s1_on_shell)
		{
		    const value_type* m=static_cast<Dirac_delta_type*>(sgen1_)->m;
		    value_type mass=(m==NULL)?0:(*m);
		    sgen1->set_lower_bound(mass*mass);
		    sgen1_check=new Dirac_delta_type(m);
		    sgen1_check->set_lower_bound(mass*mass);
		}
		else
		{
		    sgen1_check=new uniform_type();
		    sgen1_check->set_lower_bound(sgen1->lower_bound());
		}
		if(s2_on_shell)
		{
		    const value_type* m=static_cast<Dirac_delta_type*>(sgen2_)->m;
		    value_type mass=(m==NULL)?0:(*m);
		    sgen2->set_lower_bound(mass*mass);
		    sgen2_check=new Dirac_delta_type(m);
		    sgen2_check->set_lower_bound(mass*mass);
		}
		else
		{
		    sgen2_check=new uniform_type();
		    sgen2_check->set_lower_bound(sgen2->lower_bound());
		}

		MC_generator<value_t>* ssgen1=new symmetric_s_generator_composition<value_t,rng_t>(sgen1,sgen2,s);
		ssgen_multichannel=new MC_generator_wrapper<value_t>(ssgen1);
		ssgen_multichannel_wrapper=new MC_generator_wrapper<value_t>(ssgen_multichannel,false);
		
		MC_generator<value_t>* ssgen2=new symmetric_s_pair_generator<value_t,rng_t>(sgen1,sgen2,s);
		ssgen_composition=new MC_generator_wrapper<value_t>(ssgen2);
		ssgen_composition_wrapper=new MC_generator_wrapper<value_t>(ssgen_composition,false);

		MC_generator<value_t>* ssgen3=new symmetric_s_pair_generator<value_t,rng_t>(sgen1_check,sgen2_check,s);
		ssgen_check=new MC_generator_wrapper<value_t>(ssgen3);
		ssgen_check_wrapper=new MC_generator_wrapper<value_t>(ssgen_check,false);
	    }

	    /* Destructor: */

	    ~s_int_tester()
	    {
		delete ssgen_multichannel_wrapper;
		delete ssgen_composition_wrapper;
		delete ssgen_multichannel;
		delete ssgen_composition;
		delete ssgen_check_wrapper;
		delete ssgen_check;
		delete sgen1_check;
		delete sgen2_check;
	    }

	    /* Refreshes upper and lower limits: */

	    bool refresh_s_bounds()
	    {
		value_type smin=sgen1->lower_bound()+sgen2->lower_bound()+2*std::sqrt(sgen1->lower_bound()*sgen2->lower_bound());
		
		if(s<=smin)
		{
		    return false;
		}

		value_type smax1=s+sgen2->lower_bound()-2*std::sqrt(s*sgen2->lower_bound());
		value_type smax2=s+sgen1->lower_bound()-2*std::sqrt(s*sgen1->lower_bound());

		bool q=true;

		q&=sgen1_check->set_lower_bound(sgen1->lower_bound());
		q&=sgen2_check->set_lower_bound(sgen2->lower_bound());

		q&=sgen1->set_upper_bound(smax1);
		q&=sgen2->set_upper_bound(smax2);
		
		q&=sgen1_check->set_upper_bound(smax1);
		q&=sgen2_check->set_upper_bound(smax2);

		return q;
	    }

	    void run(std::size_t nevts)
	    {
		sgen1_check->reset();
		sgen2_check->reset();
		ssgen_multichannel->reset();
		ssgen_composition->reset();
		ssgen_check->reset();
		ssgen_multichannel_wrapper->reset();
		ssgen_composition_wrapper->reset();
		ssgen_check_wrapper->reset();

		if(refresh_s_bounds())
		{
		    for(std::size_t i=0;i<nevts;++i)
		    {
			if(ssgen_multichannel_wrapper->generate())
			{
			    ssgen_multichannel->integrand()=(value_t)1;
			    ssgen_multichannel_wrapper->integrand()=sgen1->density()*sgen2->density();
			}
			ssgen_multichannel->refresh_cross_section();
			ssgen_multichannel_wrapper->refresh_cross_section();
			if(ssgen_composition_wrapper->generate())
			{
			    ssgen_composition->integrand()=(value_t)1;
			    ssgen_composition_wrapper->integrand()=sgen1->density()*sgen2->density();
			}
			ssgen_composition->refresh_cross_section();
			ssgen_composition_wrapper->refresh_cross_section();
			if(ssgen_check_wrapper->generate())
			{
			    sgen1->value()=sgen1_check->value();
			    sgen1->evaluate_weight();
			    sgen2->value()=sgen2_check->value();
			    sgen2->evaluate_weight();
			    ssgen_check->integrand()=1;
			    ssgen_check_wrapper->integrand()=sgen1->density()*sgen2->density();
			}
			ssgen_check->refresh_cross_section();
			ssgen_check_wrapper->refresh_cross_section();
		    }
		}

		MC_norm1=ssgen_multichannel->cross_section().value;
		MC_err1=ssgen_multichannel->cross_section().error;

		MC_norm2=ssgen_composition->cross_section().value;
		MC_err2=ssgen_composition->cross_section().error;

		MC_norm3=ssgen_check->cross_section().value;
		MC_err3=ssgen_check->cross_section().error;
		
		MC_norm4=ssgen_multichannel_wrapper->cross_section().value;
		MC_err4=ssgen_multichannel_wrapper->cross_section().error;

		MC_norm5=ssgen_composition_wrapper->cross_section().value;
		MC_err5=ssgen_composition_wrapper->cross_section().error;

		MC_norm6=ssgen_check_wrapper->cross_section().value;
		MC_err6=ssgen_check_wrapper->cross_section().error;
	    }

	    /* Performs N1 runs of each N2 Monte Carlo points, with invariant
	     * mass ranging between sqrtsmin and sqrtsmax. If the data file is
	     * opened, the method will write every sample result. */

	    void run(const value_type& sqrtsmin,const value_type& sqrtsmax,std::size_t N1,std::size_t N2,const std::string& filename,const char* term=NULL)
	    {
#ifdef GNUPLOTPATH
		data_wrapper* data=new data_wrapper(&sqrts,&psvol);
#else
		data_wrapper* data=new data_wrapper(filename+".dat",&sqrts,&psvol);
#endif
		data->add_leaf(&MC_norm1);
		data->add_leaf(&MC_err1);
		data->add_leaf(&MC_norm2);
		data->add_leaf(&MC_err2);
		data->add_leaf(&MC_norm3);
		data->add_leaf(&MC_err3);
		data->add_leaf(&MC_norm4);
		data->add_leaf(&MC_err4);
		data->add_leaf(&MC_norm5);
		data->add_leaf(&MC_err5);
		data->add_leaf(&MC_norm6);
		data->add_leaf(&MC_err6);
		
		sqrts=std::min(sqrtsmin,sqrtsmax);
		value_type delta=std::abs(sqrtsmax-sqrtsmin)/N1;
		for(std::size_t i=0;i<N1;++i)
		{
		    s=sqrts*sqrts;
		    psvol=ps_volume();
		    run(N2);
		    data->fill();
		    sqrts+=delta;
		}

		data_stream* exstr1=new data_stream(data,"1","2");
		exstr1->title="exact PS vol";
		exstr1->style="lines";
		
		data_stream* mcstr1=new data_stream(data,"1","3","4");
		mcstr1->title="MC PS vol (bnd check)";
		mcstr1->style="yerrorbars";
		
		data_stream* mcstr2=new data_stream(data,"1","5","6");
		mcstr2->title="MC PS vol";
		mcstr2->style="yerrorbars";
		
		data_stream* mcstr3=new data_stream(data,"1","7","8");
		mcstr3->title="uni MC PS vol";
		mcstr3->style="yerrorbars";
		
		plot_script* psplot=new plot_script(filename);
		psplot->ylog=true;
		psplot->add_plot(exstr1);
		psplot->add_plot(mcstr1);
		psplot->add_plot(mcstr2);
		psplot->add_plot(mcstr3);

		data_stream* mcstr4=new data_stream(data,"1","9","10");
		mcstr4->title="MC integral (bnd check)";
		mcstr4->style="yerrorbars";

		data_stream* mcstr5=new data_stream(data,"1","11","12");
		mcstr5->title="MC integral";
		mcstr5->style="yerrorbars";

		data_stream* mcstr6=new data_stream(data,"1","13","14");
		mcstr6->title="uni MC integral";
		mcstr6->style="yerrorbars";
		
		plot_script* intplot=new plot_script(filename);
		intplot->ylog=true;
		intplot->add_plot(mcstr4);
		intplot->add_plot(mcstr5);
		intplot->add_plot(mcstr6);

		multi_plot plt(1,2,filename,term);
		plt.add_plot(psplot,0,0);
		plt.add_plot(intplot,0,1);
		plt.plot();

		delete psplot;
		delete intplot;
	    }
	
	private:

	    value_type s,sqrts;
	    
	    s_generator_type* sgen1;
	    s_generator_type* sgen2;

	    MC_integrator<value_t>* ssgen_multichannel;
	    MC_integrator<value_t>* ssgen_multichannel_wrapper;
	    MC_integrator<value_t>* ssgen_composition;
	    MC_integrator<value_t>* ssgen_composition_wrapper;
	    
	    s_generator_type* sgen1_check;
	    s_generator_type* sgen2_check;

	    MC_integrator<value_t>* ssgen_check;
	    MC_integrator<value_t>* ssgen_check_wrapper;

	    bool s1_on_shell,s2_on_shell;
	    std::string filename;

	    value_type psvol;
	    value_type MC_norm1,MC_err1;
	    value_type MC_norm2,MC_err2;
	    value_type MC_norm3,MC_err3;
	    value_type MC_norm4,MC_err4;
	    value_type MC_norm5,MC_err5;
	    value_type MC_norm6,MC_err6;

	    /* Returns the current phase space volume: */

	    value_type ps_volume() const
	    {
		if(s<sgen1->lower_bound()+sgen2->lower_bound()+2*std::sqrt(sgen1->lower_bound()*sgen2->lower_bound()))
		{
		    return 0;
		}
		if(!s1_on_shell && !s2_on_shell)
		{
		    value_type smax2=s-sgen1->lower_bound();
		    value_type smax1=s-sgen2->lower_bound();
		    return smax2*(0.5*smax2-4*sqrts*std::sqrt(smax2)/(value_type)3+smax1);
		}
		if(s1_on_shell && !s2_on_shell)
		{
		    return s+sgen1->lower_bound()-sgen2->lower_bound()-2*std::sqrt(s*sgen1->lower_bound());
		}
		if(s2_on_shell && !s1_on_shell)
		{
		    return s+sgen2->lower_bound()-sgen1->lower_bound()-2*std::sqrt(s*sgen2->lower_bound());
		}
		return 1;
	    }
    };
}

#endif /*S_INT_TESTER_H_*/

