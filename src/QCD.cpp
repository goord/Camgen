//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/QCD.h>
#include <Camgen/su(n).h>
#include <Camgen/vvv.h>
#include <Camgen/vvvv.h>
#include <Camgen/vff.h>
#include <Camgen/f.h>
#include <Camgen/ff_contr.h>
#include <Camgen/T.h>
#include <Camgen/vector_particle.h>
#include <Camgen/fermion.h>

namespace Camgen
{
    /* Compile-time constants: */

    const std::size_t QCD::dimension;
    const bool QCD::coloured;
    const bool QCD::continuous_colours;
    const bool QCD::continuous_helicities;
    const int QCD::beam_direction;
    const std::size_t QCD::N_c;

    /* Input parameters: */

    bool QCD::four_gluon_vertex=false;
    int QCD::gauge=0;
    const QCD::value_type QCD::alpha=-(QCD::value_type)1;
    QCD::value_type QCD::alpha_s=(QCD::value_type)SM_params::alpha_s;
    QCD::value_type QCD::QCD_scale=(QCD::value_type)SM_params::QCD_scale;
    const QCD::value_type QCD::pi=std::acos(-(QCD::value_type)1);
    QCD::value_type QCD::g_s=std::sqrt((QCD::value_type)4*QCD::pi*QCD::alpha_s);
    QCD::value_type QCD::M_u=(QCD::value_type)0;
    QCD::value_type QCD::M_d=(QCD::value_type)0;
    QCD::value_type QCD::M_c=(QCD::value_type)SM_params::M_c;
    QCD::value_type QCD::M_s=(QCD::value_type)0;
    QCD::value_type QCD::M_t=(QCD::value_type)SM_params::M_t;
    QCD::value_type QCD::M_b=(QCD::value_type)SM_params::M_b;

    /* Values of the couplings: */

    std::complex<QCD::value_type>QCD::ggg(QCD::g_s,0);
    std::complex<QCD::value_type>QCD::gggg(0,-QCD::g_s*QCD::g_s);
    std::complex<QCD::value_type>QCD::Tgg(0,std::sqrt((QCD::value_type)0.5)*QCD::g_s);
    std::complex<QCD::value_type>QCD::guu(0,QCD::g_s);
    std::complex<QCD::value_type>QCD::gdd(0,QCD::g_s);
    std::complex<QCD::value_type>QCD::gcc(0,QCD::g_s);
    std::complex<QCD::value_type>QCD::gss(0,QCD::g_s);
    std::complex<QCD::value_type>QCD::gtt(0,QCD::g_s);
    std::complex<QCD::value_type>QCD::gbb(0,QCD::g_s);

    /* Constructor: */

    QCD::QCD()
    {
	/* Gluon definition: */

	switch(gauge)
	{
	    case 0:
		add_gluon<Feynman_gauge,SU<N_c> >();
		break;
	    case 1:
		add_gluon<unitary_gauge,SU<N_c> >();
		break;
	    case 2:
		add_gluon<R_vector_gauge,SU<N_c> >();
		break;
	    default:
		add_gluon<Feynman_gauge,SU<N_c> >();
	}
	
	/* Gluon self-interactions: */
	
	add_vertex<colour_tensor::f< SU<N_c> >,vvv>("g","g","g",&ggg);
	if(four_gluon_vertex)
	{
	    add_vertex<colour_tensor::ff_contr< SU<N_c> >,vvvv>("g","g","g","g",&gggg);
	}
	else
	{
	    add_fast_4g_vertex< SU<N_c> >("g",&Tgg);
	}
	
	/* Quark definitions: */
	
	(M_d>(value_type)0)?add_quarks< SU<N_c> >("d","dbar",&M_d,1):add_quarks< SU<N_c> >("d","dbar",1);
	(M_u>(value_type)0)?add_quarks< SU<N_c> >("u","ubar",&M_u,2):add_quarks< SU<N_c> >("u","ubar",2);
	(M_s>(value_type)0)?add_quarks< SU<N_c> >("s","sbar",&M_s,3):add_quarks< SU<N_c> >("s","sbar",3);
	(M_c>(value_type)0)?add_quarks< SU<N_c> >("c","cbar",&M_c,4):add_quarks< SU<N_c> >("c","cbar",4);
	(M_b>(value_type)0)?add_quarks< SU<N_c> >("b","bbar",&M_b,5):add_quarks< SU<N_c> >("b","bbar",5);
	(M_t>(value_type)0)?add_quarks< SU<N_c> >("t","tbar",&M_t,6):add_quarks< SU<N_c> >("t","tbar",6);
	
	/* Gluon-quark interactions: */

	add_vertex< colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","dbar","d",&gdd);
	add_vertex< colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","ubar","u",&guu);
	add_vertex< colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","sbar","s",&gss);
	add_vertex< colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","cbar","c",&gcc);
	add_vertex< colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","bbar","b",&gbb);
	add_vertex< colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","tbar","t",&gtt);

	/* Quark family definition: */

	construct_family("q","u,d,c,s");

	/* Quark family definition including b's: */

	construct_family("Q","u,d,c,s,b");

	/* Quark family definition including b's and tops: */

	construct_family("Qt","u,d,c,s,b,t");

	/* Antiquark family definition: */

	construct_family("qbar","dbar,ubar,sbar,cbar");

	/* Antiquark family definition, including b's: */

	construct_family("Qbar","dbar,ubar,sbar,cbar,bbar");

	/* Antiquark family definition, including b's and tops: */

	construct_family("Qtbar","dbar,ubar,sbar,cbar,bbar,tbar");

	/* Up-type quark family definition: */

	construct_family("q_up","u,c");

	/* Up-type quark family definition, including b's: */

	construct_family("Q_up","u,c");

	/* Up-type quark family definition, including b's and tops: */

	construct_family("Qt_up","u,c,t");

	/* Up-type antiquark family definition: */

	construct_family("qbar_up","ubar,cbar");

	/* Up-type antiquark family definition, including b's: */

	construct_family("Qbar_up","ubar,cbar");

	/* Up-type antiquark family definition, including b's and tops: */

	construct_family("Qtbar_up","ubar,cbar,tbar");

	/* Down-type quark family definition: */

	construct_family("q_down","d,s");

	/* Down-type quark family definition, including b's: */

	construct_family("Q_down","d,s,b");

	/* Down-type quark family definition, including b's and tops: */

	construct_family("Qt_down","d,s,b");

	/* Down-type antiquark family definition: */

	construct_family("qbar_down","dbar,sbar");

	/* Down-type antiquark family definition, including b's: */

	construct_family("Qbar_down","dbar,sbar,bbar");

	/* Down-type antiquark family definition, including b's and tops: */

	construct_family("Qtbar_down","dbar,sbar,bbar");

	/* Positively-charged quark family definition: */

	construct_family("q+","u,dbar,c,sbar");

	/* Positively-charged quark family definition including b's: */

	construct_family("Q+","u,dbar,c,sbar,bbar");

	/* Positively-charged quark family definition including b's and
	 * tops: */

	construct_family("Qt+","u,dbar,c,sbar,bbar,t");

	/* Negatively-charged quark family definition: */

	construct_family("q-","d,ubar,s,cbar");

	/* Negatively-charged quark family definition, including b's: */

	construct_family("Q-","d,ubar,s,cbar,b");

	/* Negatively-charged quark family definition, including b's and
	 * tops: */

	construct_family("Qt-","d,ubar,s,cbar,b,tbar");

	/* Proton constituent parton family definition: */

	construct_family("p","u,ubar,d,dbar,c,cbar,s,sbar,g");

	/* Proton constituent parton family definition, including b's: */

	construct_family("P","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

	/* QCD jet parton family definition: */

	construct_family("j","u,ubar,d,dbar,c,cbar,s,sbar,g");

	/* QCD jet parton family definition, including b's: */

	construct_family("J","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

	/* QCD jet parton family definition, including b's and tops: */

	construct_family("Jt","u,ubar,d,dbar,c,cbar,s,sbar,t,tbar,b,bbar,g");
    }
	    
    /* Output of the number of QCD colours: */

    std::size_t QCD::QCD_colours()
    {
	return N_c;
    }

    /* Dummy function setting the EM fine-structure constant: */

    void QCD::set_alpha(const value_type& a){}

    /* Function setting alpha strong and computing the couplings: */

    void QCD::set_alpha_s(const value_type& a)
    {
	alpha_s=a;
	refresh_couplings();
    }

    /* Function setting the scale and computing alpha_s (at LO) and the
     * couplings: */

    void QCD::set_QCD_scale(const value_type& mu)
    {
	QCD_scale=mu;
	value_type rho=mu/(value_type)SM_params::QCD_scale;
	alpha_s=(value_type)SM_params::alpha_s/((value_type)1+(value_type)SM_params::alpha_s*(11*N_c-12)*std::log(rho)/((value_type)6*pi));
	refresh_couplings();
    }
	    
    /* Function computing the couplings from the current alpha_s value:
     * */
	    
    void QCD::refresh_couplings()
    {
	g_s=std::sqrt((value_type)4*pi*alpha_s);
	ggg=std::complex<value_type>(g_s,0);
	gggg=std::complex<value_type>(0,-g_s*g_s);
	Tgg=std::complex<value_type>(0,std::sqrt((value_type)0.5)*g_s);
	guu=std::complex<value_type>(0,g_s);
	gdd=guu;
	gcc=guu;
	gss=guu;
	gtt=guu;
	gbb=guu;
    }

    /* Refresh quark masses (to use when along the run massless fermions are
     * to be made massive): */

    void QCD::refresh_fermion_masses()
    {
	if(initialised())
	{
	    (M_u>(value_type)0)?set_mass("u",&M_u):set_massless("u");
	    (M_d>(value_type)0)?set_mass("d",&M_d):set_massless("d");
	    (M_c>(value_type)0)?set_mass("c",&M_c):set_massless("c");
	    (M_s>(value_type)0)?set_mass("s",&M_s):set_massless("s");
	}
    }

    /* Function setting the gluon-propagator to the Feynman gauge: */

    void QCD::set_Feynman_gauge()
    {
	if(initialised() and gauge!=0)
	{
	    set_propagator<Feynman_gauge>("g");
	}
	gauge=0;
    }

    /* Function setting the gluon-propagator to the unitary gauge: */

    void QCD::set_unitary_gauge()
    {
	if(initialised() and gauge!=1)
	{
	    set_propagator<unitary_gauge>("g");
	}
	gauge=1;
    }

    /* Function setting the gluon-propagator to the R-xi gauge: */

    void QCD::set_R_xi_gauge()
    {
	if(initialised() and gauge!=2)
	{
	    set_propagator<R_vector_gauge>("g");
	}
	gauge=2;
    }
	    
    /* Function setting the xi-parameter: */
	    
    void QCD::set_xi(const value_type& x)
    {
	R_gauge<QCD>::xi=x;
    }

    /* Function switching to a 4-gluon vertex Lagrangian: */

    void QCD::set_4_gluon_vertex()
    {
	if(initialised() and !four_gluon_vertex)
	{
	    erase_vertex("H_qcd","g","g");
	    add_vertex<colour_tensor::ff_contr< SU<N_c> >,vvvv>("g","g","g","g",&gggg);
	}
	four_gluon_vertex=true;
    }

    /* Function switching to an auxiliary tensor field description of
     * the 4-gluon vertex: */

    void QCD::set_auxiliary_QCD_field()
    {
	if(initialised() and four_gluon_vertex)
	{
	    erase_vertex("g","g","g","g");
	    add_fast_4g_vertex< SU<N_c> >("g",&Tgg);
	}
	four_gluon_vertex=false;
    }

    /* Reset all parameters using the default input parameters in the
     * SM_params.h file: */

    void QCD::set_default_params()
    {
	M_u=(value_type)0;
	M_d=(value_type)0;
	M_c=(value_type)SM_params::M_c;
	M_s=(value_type)0;
	M_t=(value_type)SM_params::M_t;
	M_b=(value_type)SM_params::M_b;
	refresh_fermion_masses();
	QCD_scale=(value_type)SM_params::QCD_scale;
	alpha_s=(value_type)SM_params::alpha_s;
	refresh_couplings();
	set_Feynman_gauge();
	set_auxiliary_QCD_field();
    }
}

