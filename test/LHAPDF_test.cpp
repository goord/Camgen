//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <config.h>
#include <Camgen/plt_strm.h>
#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/procgen_fac.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Program plotting the parton density functions wrapped into the pdf_wrapper  *
 * and comparing x-parameterised hadronic initial states with y-parameterised  *
 * parton momentum generation.                                                 *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    typedef SM model_type;
    typedef SM::value_type value_type;
    typedef std::random rn_engine;

    license_print::disable();

#if HAVE_LHAPDF_H_

    std::cout<<"-------------------------------------------"<<std::endl;
    std::cout<<"testing hadronic initial states............"<<std::endl;
    std::cout<<"-------------------------------------------"<<std::endl;
    
    set_pdfs("cteq6l.LHpdf",1);
    pdf_wrapper::initialise(pdf_name(),pdf_number());
    bool have_gp=plot_config::gnuplot_path!=NULL;

    {
	std::string filename("plots/LHAPDF");
	std::string fext=have_gp?".eps":".dat/.gp";
	value_type x,g,u,ubar,d,dbar,s,sbar,c,cbar,b,bbar;

	data_wrapper* data=have_gp?new data_wrapper(&x,&g):new data_wrapper(filename+".dat",&x,&g);
	data->add_leaf(&u);
	data->add_leaf(&ubar);
	data->add_leaf(&d);
	data->add_leaf(&dbar);
	data->add_leaf(&c);
	data->add_leaf(&cbar);
	data->add_leaf(&s);
	data->add_leaf(&sbar);
	data->add_leaf(&b);
	data->add_leaf(&bbar);

	double mu=model_type::M_Z;

	std::cerr<<"Evaluating pdf sets for all partons at MZ"<<'.';

	for(int i=0;i<1000;++i)
	{
	    if(i%100==0)
	    {
		std::cerr<<'.';
	    }
	    x=((double)i)/1000;
	    if(x<pdf_wrapper::xmin())
	    {
		continue;
	    }
	    if(x>pdf_wrapper::xmax())
	    {
		break;
	    }
	    g=pdf_wrapper::xf(x,0,mu);
	    u=pdf_wrapper::xf(x,2,mu);
	    ubar=pdf_wrapper::xf(x,-2,mu);
	    d=pdf_wrapper::xf(x,1,mu);
	    dbar=pdf_wrapper::xf(x,-1,mu);
	    c=pdf_wrapper::xf(x,4,mu);
	    cbar=pdf_wrapper::xf(x,-4,mu);
	    s=pdf_wrapper::xf(x,3,mu);
	    sbar=pdf_wrapper::xf(x,-3,mu);
	    b=pdf_wrapper::xf(x,5,mu);
	    bbar=pdf_wrapper::xf(x,-5,mu);

	    data->fill();
	}

	data->write();

	data_stream* datstr_g=new data_stream(data,"1","2");
	datstr_g->title="g";
	datstr_g->style="lines";
	data_stream* datstr_u=new data_stream(data,"1","3");
	datstr_u->title="u";
	datstr_u->style="lines";
	data_stream* datstr_ubar=new data_stream(data,"1","4");
	datstr_ubar->title="ubar";
	datstr_ubar->style="lines";
	data_stream* datstr_d=new data_stream(data,"1","5");
	datstr_d->title="d";
	datstr_d->style="lines";
	data_stream* datstr_dbar=new data_stream(data,"1","6");
	datstr_dbar->title="dbar";
	datstr_dbar->style="lines";
	data_stream* datstr_c=new data_stream(data,"1","7");
	datstr_c->title="c";
	datstr_c->style="lines";
	data_stream* datstr_cbar=new data_stream(data,"1","8");
	datstr_cbar->title="cbar";
	datstr_cbar->style="lines";
	data_stream* datstr_s=new data_stream(data,"1","9");
	datstr_s->title="s";
	datstr_s->style="lines";
	data_stream* datstr_sbar=new data_stream(data,"1","10");
	datstr_sbar->title="sbar";
	datstr_sbar->style="lines";
	data_stream* datstr_b=new data_stream(data,"1","11");
	datstr_b->title="b";
	datstr_b->style="lines";
	data_stream* datstr_bbar=new data_stream(data,"1","12");
	datstr_bbar->title="bbar";
	datstr_bbar->style="lines";

	plot_script* lhapdfplot=new plot_script(filename,"postscript enhanced color");
	lhapdfplot->title="PDFs at mu_F = MZ";
	lhapdfplot->xlabel="x";
	lhapdfplot->ylabel="xf(x)";
	lhapdfplot->ylog=true;
	lhapdfplot->add_plot(datstr_g);
	lhapdfplot->add_plot(datstr_u);
	lhapdfplot->add_plot(datstr_ubar);
	lhapdfplot->add_plot(datstr_d);
	lhapdfplot->add_plot(datstr_dbar);
	lhapdfplot->add_plot(datstr_c);
	lhapdfplot->add_plot(datstr_cbar);
	lhapdfplot->add_plot(datstr_s);
	lhapdfplot->add_plot(datstr_sbar);
	lhapdfplot->add_plot(datstr_b);
	lhapdfplot->add_plot(datstr_bbar);
	lhapdfplot->plot();
	std::cerr<<"done, created "<<filename<<fext<<'.'<<std::endl;
    }

    {
	std::string filename("plots/had_is");
	std::string fext=have_gp?".eps":".dat/.gp";
	value_type ecm,xsec_xx,xsec_xx_err,xsec_sy,xsec_sy_err;

	data_wrapper* data=have_gp?new data_wrapper(&ecm,&xsec_xx,&xsec_xx_err):new data_wrapper(filename+".dat",&ecm,&xsec_xx,&xsec_xx_err);
	data->add_leaf(&xsec_sy);
	data->add_leaf(&xsec_sy_err);

	set_helicity_generator_type(helicity_generators::uniform);
	set_colour_generator_type(colour_generators::flow_sampling);
	set_grid_init(10,100);
	set_auto_grid_adapt(100);

	std::string process("u,dbar > e+,nu_e");

	CM_algorithm<model_type,2,2>amplitude(process);
	amplitude.load();
	amplitude.construct_trees();

	set_phase_space_generator_type(phase_space_generators::uniform);

	scale_expression<value_type>* scale=new scale_wrapper<value_type>(&model_type::M_Z);

	set_initial_state_type(initial_states::proton_proton_xx);

	process_generator_factory<model_type,2,2,rn_engine>factory;

	process_generator<model_type,2,2,rn_engine>* gen_xx=factory.create_generator(amplitude.get_tree_iterator());
	gen_xx->insert_scale(scale);

	set_initial_state_type(initial_states::proton_proton);
	process_generator<model_type,2,2,rn_engine>* gen_sy=factory.create_generator(amplitude.get_tree_iterator());
	gen_sy->insert_scale(scale);

	int n=0;

	std::cerr<<"Checking xx-parameterised against sy-parameterised pdf-sampling for "<<process<<'.';

	for(double Ecm=10;Ecm<210;Ecm+=10)
	{
	    double Ebeam=0.5*Ecm;
	    gen_xx->set_beam_energy(-1,Ebeam);
	    gen_xx->set_beam_energy(-2,Ebeam);
	    if(n==0)
	    {
		initialise(gen_xx);
	    }
	    gen_sy->set_beam_energy(-1,Ebeam);
	    gen_sy->set_beam_energy(-2,Ebeam);
	    if(n==0)
	    {
		initialise(gen_sy);
	    }
	    for(int i=0;i<100000;++i)
	    {
		gen_xx->generate();
		gen_sy->generate();
	    }

	    ecm=Ecm;
	    xsec_xx=gen_xx->cross_section().value;
	    xsec_xx_err=gen_xx->cross_section().error;
	    xsec_sy=gen_sy->cross_section().value;
	    xsec_sy_err=gen_sy->cross_section().error;
	    data->fill();

	    gen_xx->reset_cross_section();
	    gen_sy->reset_cross_section();

	    ++n;

	    if(n%2==0)
	    {
		std::cerr<<'.';
	    }
	}
	data->write();
	data_stream* datstr_xx=new data_stream(data,"1","2","3");
	datstr_xx->title="x-param";
	datstr_xx->style="yerrorbars";
	data_stream* datstr_sy=new data_stream(data,"1","4","5");
	datstr_sy->title="y-param";
	datstr_sy->style="yerrorbars";
	plot_script* lhapdfplot=new plot_script(filename,"postscript enhanced color");
	lhapdfplot->title=process;
	lhapdfplot->xlabel="Ecm (GeV)";
	lhapdfplot->ylabel="sigma (pb)";
	lhapdfplot->ylog=true;
	lhapdfplot->add_plot(datstr_xx);
	lhapdfplot->add_plot(datstr_sy);
	lhapdfplot->plot();
	delete scale;
	delete gen_xx;
	delete gen_sy;
	std::cerr<<"done, created "<<filename<<fext<<'.'<<std::endl;
    }

#else

    std::cerr<<"Skipping LHAPDF checks: no pdf library located."<<std::endl;

#endif

}

