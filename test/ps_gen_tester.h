//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef PS_GEN_TESTER_H_
#define PS_GEN_TESTER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing class template for phase space Monte-Carlo generators. The 'run'  *
 * function performs a check on every event and creates a plot comparing     *
 * respectively the phase space volume and the cross section with uniformly  *
 * generated points.                                                         *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/plt_config.h>
#include <Camgen/type_holders.h>
#include <Camgen/isgen_fac.h>
#include <Camgen/psgen_fac.h>
#include <Camgen/multi_plot.h>

namespace Camgen
{
    template<class model_t,std::size_t N_out>class pT_min_cut: public phase_space_cut
    {
	public:

	    typedef typename model_t::value_type value_type;

	    ps_generator_base<model_t>* generator;

	    value_type value;

	    int particle;

	    pT_min_cut(const value_type& value_,int i=0):generator(NULL),value(value_),particle(i){}

	    bool pass()
	    {
		if(generator==NULL)
		{
		    return true;
		}
		if(particle==0)
		{
		    for(int i=1;i<=(int)N_out;++i)
		    {
			if(generator->pT(i)<value)
			{
			    return false;
			}
		    }
		    return true;
		}
		return !(generator->pT(particle)<value);
	    }

	    std::ostream& print(std::ostream& os) const
	    {
		if(particle==0)
		{
		    os<<"pT >= "<<value;
		}
		else
		{
		    os<<"pT("<<particle<<") >= "<<value;
		}
		return os;
	    }
    };

    /* Monte-Carlo phase space generator class template declaration: */

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_generator_tester
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_t::value_type value_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
	    typedef vector<value_type,spacetime_type::dimension> momentum_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef ps_generator_factory<model_t,N_in,N_out,rng_t> gen_factory;
	    typedef typename CM_algorithm<model_t,N_in,N_out>::tree_iterator CM_tree_iterator;

	    /* Output file name: */

	    const std::string filename;
	    
	    /* Process symmetry factor: */
	    
	    const value_type sym_factor;
	    
	    /* Amplitude process tree iterator: */
	    
	    CM_tree_iterator amplitude;

	    /* Phase space generator: */

	    ps_generator<model_t,N_in,N_out>* ps_gen;
	    
	    /* Uniform phase space generator for comparison purposes: */
	    
	    ps_generator<model_t,N_in,N_out>* uni_gen;

	    /* psgen adaptation batch size: */

	    size_type adapt_batch_size;

	    /* Constructor from a CM-tree iterator, an initial-state generator instance
	     * and a file name: */

	    ps_generator_tester(CM_tree_iterator amplitude_,const std::string& filename_,initial_states::type isgen_type=Camgen::initial_state_type(),phase_space_generators::type psgen_type=Camgen::phase_space_generator_type()):filename(filename_),sym_factor(amplitude_->symmetry_factor()),amplitude(amplitude_),adapt_batch_size(0),pTcut1(NULL),pTcut2(NULL)
	    {
		ps_gen=gen_factory::create_generator(amplitude,isgen_type,psgen_type);
		uni_gen=gen_factory::create_generator(amplitude,isgen_type,phase_space_generators::uniform);
	    }

	    /* Sets the i-th beam energy: */

	    bool set_beam_energy(int i,const value_type& E)
	    {
		return (ps_gen->set_beam_energy(i,E) and uni_gen->set_beam_energy(i,E));
	    }

	    /* Run method, performing N_events calls to the MC generator and RAMBO, and
	     * recording N_points cross section along the run: */

	    bool run(size_type N_events,size_type N_points,bool verbose=false,bool skip_plot=false,bool skip_check=false)
	    {
		size_type N_batch=N_events/N_points;
		value_type n,psvol1,psvol1err,psvol2,psvol2err,xsec1,xsec1err,xsec2,xsec2err;
		bool have_gp=plot_config::gnuplot_path!=NULL;
		data_wrapper* data=have_gp?new data_wrapper(&n,&psvol1,&psvol1err):new data_wrapper(filename+".dat",&n,&psvol1,&psvol1err);
		data->add_leaf(&psvol2);
		data->add_leaf(&psvol2err);
		data->add_leaf(&xsec1);
		data->add_leaf(&xsec1err);
		data->add_leaf(&xsec2);
		data->add_leaf(&xsec2err);

		value_type w1sum=0;
		value_type w1sqsum=0;
		value_type w2sum=0;
		value_type w2sqsum=0;

		ps_gen->refresh_Ecm();
		uni_gen->refresh_Ecm();
		ps_gen->refresh_m_min();
		uni_gen->refresh_m_min();

		for(size_type i=0;i<N_events;++i)
		{
		    bool ps_generated=ps_gen->generate();

		    if(ps_generated)
		    {
			if(!skip_check and !ps_gen->check())
			{
			    return false;
			}
		    }
		    value_type w1=(ps_generated and ps_gen->pass())?ps_gen->weight():(value_type)0;
		    w1sum+=w1;
		    w1sqsum+=(w1*w1);
		    amplitude->reset();
		    if(w1!=(value_type)0)
		    {
			ps_gen->integrand()*=std::norm(amplitude->evaluate())*sym_factor;
			if(ps_gen->weight()!=ps_gen->weight())
			{
			    ps_gen->print(std::cerr);
			    Camgen::log<<"Invalid phase space weight encountered: "<<ps_gen->weight()<<endlog;
			    return false;
			}
			if(ps_gen->integrand()!=ps_gen->integrand())
			{
			    amplitude->print(std::cerr);
			    ps_gen->print(std::cerr);
			    Camgen::log<<"Invalid phase space integrand encountered: "<<ps_gen->integrand()<<endlog;
			    return false;
			}
		    }
		    else
		    {
			ps_gen->integrand()=(value_type)0;
		    }
		    ps_gen->update();
		    if(verbose)
		    {
			std::cerr<<i<<": w = "<<ps_gen->weight()<<", M ="<<ps_gen->integrand()<<std::endl;
		    }
		    if(adapt_batch_size>0 and i%adapt_batch_size==0)
		    {
			ps_gen->adapt();
		    }
		    ps_gen->refresh_cross_section();
		    if(w1!=(value_type)0)
		    {
			ps_gen->evaluate_weight();
			value_type w_check=ps_gen->weight();
			if(!equals(w_check,w1))
			{
			    ps_gen->print(std::cerr);
			    Camgen::log<<"Weight re-evaluation lead to different result: "<<w1<<" not equal to  "<<w_check<<endlog;
			}
		    }
		    
		    bool rambo_generated=uni_gen->generate();
		    if(rambo_generated)
		    {
			if(!skip_check and !uni_gen->check())
			{
			    return false;
			}
		    }
		    value_type w2=(rambo_generated and uni_gen->pass())?uni_gen->weight():(value_type)0;
		    
		    w2sum+=w2;
		    w2sqsum+=(w2*w2);
		    amplitude->reset();
		    if(w2!=(value_type)0)
		    {
			uni_gen->integrand()*=std::norm(amplitude->evaluate())*sym_factor;
		    }
		    else
		    {
			ps_gen->integrand()=(value_type)0;
		    }
		    uni_gen->update();
		    uni_gen->refresh_cross_section();

		    if(i!=0 and i%N_batch==0)
		    {
			n=(value_type)i;
			value_type n2=n*(n-(value_type)1);
			psvol1=w1sum/n;
			psvol1err=std::sqrt(std::abs((w1sqsum-w1sum*w1sum/n)/n2));
			psvol2=w2sum/n;
			psvol2err=std::sqrt(std::abs((w2sqsum-w2sum*w2sum/n)/n2));
			xsec1=(ps_gen->cross_section()).value;
			xsec1err=(ps_gen->cross_section()).error;
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
		    plot_script* psvolplot=new plot_script("Phase space volume");
		    psvolplot->add_plot(datstr1);
		    psvolplot->add_plot(datstr2);
		    double dx=N_events/5;
		    for(size_type i=0;i<5;++i)
		    {
			psvolplot->add_x_tic(i*dx);
		    }

		    data_stream* datstr3=new data_stream(data,"1","6","7");
		    datstr3->title="MC";
		    datstr3->style="yerrorbars";
		    data_stream* datstr4=new data_stream(data,"1","8","9");
		    datstr4->title="RAMBO";
		    datstr4->style="yerrorbars";
		    plot_script* xsecplot=new plot_script("Cross section");
		    xsecplot->add_plot(datstr3);
		    xsecplot->add_plot(datstr4);
		    for(size_type i=0;i<5;++i)
		    {
			xsecplot->add_x_tic(i*dx);
		    }

		    multi_plot plotall(1,2,filename,"postscript enhanced color");
		    plotall.add_plot(psvolplot,0,0);
		    plotall.add_plot(xsecplot,0,1);
		    plotall.plot();

		    delete xsecplot;
		}

		delete data;

		MC_integral<value_type>psvol_rambo(psvol1,psvol1err);
		MC_integral<value_type>psvol_tree(psvol2,psvol2err);

		return equals(ps_gen->cross_section(),uni_gen->cross_section()) or equals(psvol_tree,psvol_rambo);
	    }

	    /* Destructor: */

	    ~ps_generator_tester()
	    {
		if(pTcut1!=NULL)
		{
		    delete pTcut1;
		}
		if(pTcut2!=NULL)
		{
		    delete pTcut2;
		}
		delete ps_gen;
		delete uni_gen;
	    }

	    /* Sets a minimal dimass: */

	    bool set_m_min(size_type i,size_type j,const value_type& x)
	    {
		bool q=true;
		if(i!=0 and j!=0 and i<=N_out and j<=N_out)
		{
		    q&=ps_gen->set_m_min(i,j,x);
		    q&=uni_gen->set_m_min(i,j,x);
		}
		return false;
	    }

	    /* Sets a minimal pT: */

	    bool set_pT_min(const value_type& value,int i=0)
	    {
		if(pTcut1!=NULL)
		{
		    delete pTcut1;
		}
		pTcut1 = new pT_min_cut<model_t,N_out>(value,i);
		pTcut1->generator=uni_gen;
		uni_gen->insert_cut(pTcut1);

		if(pTcut2!=NULL)
		{
		    delete pTcut2;
		}
		pTcut2 = new pT_min_cut<model_t,N_out>(value,i);
		pTcut2->generator=ps_gen;
		ps_gen->insert_cut(pTcut2);
		return true;
	    }

	    /* Prints the amplitude: */

	    std::ostream& print_amplitude(std::ostream& os) const
	    {
		amplitude->print(os);
		return os;
	    }

	    /* Prints the phase space generator: */

	    std::ostream& print_generator(std::ostream& os) const
	    {
		ps_gen->print(os);
		return os;
	    }

	    std::ostream& print_status(std::ostream& os) const
	    {
		os<<"RAMBO status:"<<std::endl;
		uni_gen->print_status(std::cout);
		os<<"generator status: "<<std::endl;
		ps_gen->print_status(std::cout);
		return os;
	    }

	private:

	    pT_min_cut<model_t,N_out>* pTcut1;
	    pT_min_cut<model_t,N_out>* pTcut2;
    };
}

#endif /*PS_GEN_TESTER_H_*/

