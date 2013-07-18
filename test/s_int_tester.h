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

#include <Camgen/s_gen.h>
#include <Camgen/multiplot.h>

namespace Camgen
{
    template<class value_t,class rng_t>class s_int_tester: public MC_generator<value_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;

	    /* Total invariant mass: */

	    value_type sqrts;

	    /* Maximal number of invariant masses thrown per event: */

	    std::size_t max_throws;

	    /* Constructor: */

	    s_int_tester(s_generator<value_t,rng_t>* sgen1_,s_generator<value_t,rng_t>* sgen2_):max_throws(10000),sgen1(sgen1_),sgen2(sgen2_),needs_refreshing(true),check_bounds(true)
	    {
		if(dynamic_cast<Dd_s_generator<value_t,rng_t>*>(sgen1_)!=NULL)
		{
		    const value_type* m=static_cast<Dd_s_generator<value_t,rng_t>*>(sgen1_)->m;
		    value_type mass=(m==NULL)?0:(*m);
		    sgen1->set_m_min(mass);
		    check_sgen1=new Dd_s_generator<value_t,rng_t>(m);
		    check_sgen1->set_m_min(mass);
		}
		else
		{
		    check_sgen1=new uni_s_generator<value_t,rng_t>;
		}
		if(dynamic_cast<Dd_s_generator<value_t,rng_t>*>(sgen2_)!=NULL)
		{
		    const value_type* m=static_cast<Dd_s_generator<value_t,rng_t>*>(sgen2_)->m;
		    value_type mass=(m==NULL)?0:(*m);
		    sgen2->set_m_min(mass);
		    check_sgen2=new Dd_s_generator<value_t,rng_t>(m);
		    check_sgen2->set_m_min(mass);
		}
		else
		{
		    check_sgen2=new uni_s_generator<value_t,rng_t>;
		}
	    }

	    /* Destructor: */

	    ~s_int_tester()
	    {
		delete check_sgen1;
		delete check_sgen2;
	    }

	    bool refresh_s_bounds()
	    {
		if(needs_refreshing)
		{
		    if(sqrts<(value_type)0)
		    {
			return false;
		    }
		    if(!(sqrts>sgen1->m_min()+sgen2->m_min()))
		    {
			return false;
		    }
		    return (sgen1->set_m_max(sqrts-sgen2->m_min()) and sgen2->set_m_max(sqrts-sgen1->m_min()));
		}
		return true;
	    }

	    bool generate()
	    {
		if(!refresh_s_bounds())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type s=sqrts*sqrts;
		const value_type& s1=sgen1->s();
		const value_type& s2=sgen2->s();
		value_type stot;
		if(check_bounds)
		{
		    std::size_t n=0;
		    do
		    {
			sgen1->generate();
			sgen2->generate();
			stot=s1+s2+(value_type)2*std::sqrt(s1*s2);
			++n;
		    }
		    while((stot>s) and (n<max_throws));
		    if(stot>s)
		    {
			Camgen::log(log_level::warning)<<"No succesful s-pair generated after "<<n<<" throws, returning weight 0."<<endlog;
			this->weight()=(value_type)0;
			return false;
		    }
		    this->weight()=integrate(sgen1,sgen2,sqrts)*(sgen1->weight())*(sgen2->weight());
		}
		else
		{
		    sgen1->generate();
		    sgen2->generate();
		    value_type stot=s1+s2+(value_type)2*std::sqrt(s1*s2);
		    if(stot>s)
		    {
			this->weight()=(value_type)0;
			return false;
		    }
		    this->weight()=(sgen1->weight())*(sgen2->weight());
		}
		return true;
	    }

            bool evaluate_weight()
            {
		if(!refresh_s_bounds())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type s=sqrts*sqrts;
		value_type s1=sgen1->s();
		value_type s2=sgen2->s();
		if(s1+s2+std::sqrt(s1*s2)>s)
		{
		    Camgen::log(log_level::warning)<<"Energy-conservation violation encountered in s-branching...returning weight 0"<<endlog;
		    this->weight()=(value_type)0;
		    return false;
		}
		bool q=(sgen1->evaluate_weight() and sgen2->evaluate_weight());
		if(!q)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		if(check_bounds)
		{
		    this->weight()=integrate(sgen1,sgen2,sqrts)*(sgen1->weight())*(sgen2->weight());
		}
		else
		{
		    this->weight()=(sgen1->weight())*(sgen2->weight());
		}
		return true;
            }

	    /* Performs a MC run at the current invariant mass: */

	    bool run(std::size_t N)
	    {
		if(refresh_s_bounds())
		{
		    norm2=integrate_unnormalised(sgen1,sgen2,sqrts);
		    needs_refreshing=false;
		}
		else
		{
		    return false;
		}
		check_sgen1->set_m_range(sgen1->m_min(),sgen1->m_max());
		check_sgen2->set_m_range(sgen2->m_min(),sgen2->m_max());
		norm1=integrate_unnormalised(check_sgen1,check_sgen2,sqrts);
		
		/* PS volume calculation: */
		
		this->reset();
		this->integrand()=(value_type)1;
		check_bounds=true;
		for(std::size_t i=0;i<N;++i)
		{
		    generate();
		    this->refresh_cross_section();
		}
		MC_norm1=(this->cross_section()).value;
		MC_err1=(this->cross_section()).error;

		/* Density integral calculation: */

		this->reset();
		check_bounds=false;
		for(std::size_t i=0;i<N;++i)
		{
		    if(generate())
		    {
			this->integrand()=sgen1->density()*sgen2->density();
		    }
		    this->refresh_cross_section();
		}

		MC_norm2=(this->cross_section()).value;
		MC_err2=(this->cross_section()).error;

		return true;
	    }

	    /* Performs N1 runs of each N2 Monte Carlo points, with invariant
	     * mass ranging between sqrtsmin and sqrtsmax. If the data file is
	     * opened, the method will write every sample result. */

	    void run(const value_type& sqrtsmin,const value_type& sqrtsmax,std::size_t N1,std::size_t N2,const std::string& filename,const char* term=NULL)
	    {
#ifdef GNUPLOTPATH
		data_wrapper* data=new data_wrapper(&sqrts,&norm1);
#else
		data_wrapper* data=new data_wrapper(filename+".dat",&sqrts,&norm1);
#endif
		data->add_leaf(&MC_norm1);
		data->add_leaf(&MC_err1);
		data->add_leaf(&norm2);
		data->add_leaf(&MC_norm2);
		data->add_leaf(&MC_err2);
		
		value_type delta=std::abs(sqrtsmax-sqrtsmin)/N1;
		sqrts=std::min(sqrtsmin,sqrtsmax);
		for(std::size_t i=0;i<N1;++i)
		{
		    needs_refreshing=true;
		    if(run(N2))
		    {
			data->fill();
		    }
		    sqrts+=delta;
		}

		data_stream* exstr1=new data_stream(data,"1","2");
		exstr1->title="exact PS vol";
		exstr1->style="lines";
		
		data_stream* mcstr1=new data_stream(data,"1","3","4");
		mcstr1->title="MC PS vol";
		mcstr1->style="yerrorbars";
		
		plot_script* psplot=new plot_script(filename);
		psplot->ylog=true;
		psplot->add_plot(exstr1);
		psplot->add_plot(mcstr1);

		data_stream* exstr2=new data_stream(data,"1","5");
		exstr2->title="exact integral";
		exstr2->style="lines";
		
		data_stream* mcstr2=new data_stream(data,"1","6","7");
		mcstr2->title="MC integral";
		mcstr2->style="yerrorbars";
		
		plot_script* intplot=new plot_script(filename);
		intplot->add_plot(exstr2);
		intplot->add_plot(mcstr2);

		multi_plot<1,2>plt(filename,term);
		plt.add_plot(psplot,0,0);
		plt.add_plot(intplot,0,1);
		plt.plot();

		delete psplot;
		delete intplot;
	    }
	
	private:
	    
	    s_generator<value_t,rng_t>* sgen1;
	    s_generator<value_t,rng_t>* sgen2;
	    s_generator<value_t,rng_t>* check_sgen1;
	    s_generator<value_t,rng_t>* check_sgen2;
	    bool needs_refreshing,check_bounds;
	    std::string filename;
	    value_type norm1,MC_norm1,MC_err1,norm2,MC_norm2,MC_err2;
    };
}

#endif /*S_INT_TESTER_H_*/

