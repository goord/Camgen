//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_QCD_BASE_H_
#define CAMGEN_QCD_BASE_H_

#include <Camgen/su(n).h>
#include <Camgen/vvv.h>
#include <Camgen/vvvv.h>
#include <Camgen/vff.h>
#include <Camgen/f.h>
#include <Camgen/ff_contr.h>
#include <Camgen/T.h>
#include <Camgen/vector_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/model.h>
#include <Camgen/SM_params.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the QCD_base class template in Camgen.       *
 * Publicly deriving a model class model_t from QCD_base<model_t> amounts to   *
 * adding QCD to the model. It should be noted that the value type of base_QCD *
 * has to be castable to that of the derived model.                            *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */                                                                            

namespace Camgen
{
    template<class model_t,class value_t>class QCD_base: public model<model_t>
    {
	public:
	    
	    /* Value type definition (is always casted to subtype's value type):
	     * */

	    typedef value_t value_type;
	    
	    /* EM fine structure constant (unused): */

	    static const value_type alpha;

	    /* alpha strong: */

	    static value_type alpha_s;
	    
	    /* QCD scale parameter: */
	    
	    static value_type QCD_scale;

	    /* QCD coupling: */

	    static value_type g_s;
	    
	    /* Pi: */
	    
	    static const value_type pi;

	    /* Quark masses: */

	    static value_type M_u,M_d,M_c,M_s,M_t,M_b;

	    /* Triple-gluon coupling: */

	    static std::complex<value_type>ggg;
	    
	    /* Quadruple-gluon coupling: */
	    
	    static std::complex<value_type>gggg;

	    /* Tensor-gluon coupling: */

	    static std::complex<value_type>Tgg;
	    
	    /* Gluon-quark couplings: */
	    
	    static std::complex<value_type>guu,gdd,gcc,gss,gtt,gbb;

	    /* Constructor: */

	    QCD_base()
	    {
		/* Gluon definition: */

		switch(gauge)
		{
		    case 0:
			model<model_t>::template add_gluon< Feynman_gauge,SU<model_t::N_c> >();
			break;
		    case 1:
			model<model_t>::template add_gluon< unitary_gauge,SU<model_t::N_c> >();
			break;
		    case 2:
			model<model_t>::template add_gluon< R_vector_gauge,SU<model_t::N_c> >();
			break;
		    default:
			model<model_t>::template add_gluon< Feynman_gauge,SU<model_t::N_c> >();
		}

		/* Gluon self-interactions: */

		model<model_t>::template add_vertex<colour_tensor::f< SU<model_t::N_c> >,vvv>("g","g","g",&ggg);
		if(four_gluon_vertex)
		{
		    model<model_t>::template add_vertex<colour_tensor::ff_contr< SU<model_t::N_c> >,vvvv>("g","g","g","g",&gggg);
		}
		else
		{
		    model<model_t>::template add_fast_4g_vertex< SU<model_t::N_c> >("g",&Tgg);
		}

		/* Quark definitions: */

		(M_d>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("d","dbar",&M_d,1):model<model_t>::template add_quarks< SU<model_t::N_c> >("d","dbar",1);
		(M_u>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("u","ubar",&M_u,2):model<model_t>::template add_quarks< SU<model_t::N_c> >("u","ubar",2);
		(M_s>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("s","sbar",&M_s,3):model<model_t>::template add_quarks< SU<model_t::N_c> >("s","sbar",3);
		(M_c>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("c","cbar",&M_c,4):model<model_t>::template add_quarks< SU<model_t::N_c> >("c","cbar",4);
		(M_b>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("b","bbar",&M_b,5):model<model_t>::template add_quarks< SU<model_t::N_c> >("b","bbar",5);
		(M_t>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("t","tbar",&M_t,6):model<model_t>::template add_quarks< SU<model_t::N_c> >("t","tbar",6);

		/* Gluon-quark interactions: */

		model<model_t>::template add_vertex< colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","dbar","d",&gdd);
		model<model_t>::template add_vertex< colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","ubar","u",&guu);
		model<model_t>::template add_vertex< colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","sbar","s",&gss);
		model<model_t>::template add_vertex< colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","cbar","c",&gcc);
		model<model_t>::template add_vertex< colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","bbar","b",&gbb);
		model<model_t>::template add_vertex< colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","tbar","t",&gtt);

		/* Quark family definition: */

		model<model_t>::construct_family("q","u,d,c,s");

		/* Quark family definition including b's: */

		model<model_t>::construct_family("Q","u,d,c,s,b");

		/* Quark family definition including b's and tops: */

		model<model_t>::construct_family("Qt","u,d,c,s,b,t");

		/* Antiquark family definition: */

		model<model_t>::construct_family("qbar","dbar,ubar,sbar,cbar");

		/* Antiquark family definition, including b's: */

		model<model_t>::construct_family("Qbar","dbar,ubar,sbar,cbar,bbar");

		/* Antiquark family definition, including b's and tops: */

		model<model_t>::construct_family("Qtbar","dbar,ubar,sbar,cbar,bbar,tbar");

		/* Up-type quark family definition: */

		model<model_t>::construct_family("q_up","u,c");

		/* Up-type quark family definition, including b's: */

		model<model_t>::construct_family("Q_up","u,c");

		/* Up-type quark family definition, including b's and tops: */

		model<model_t>::construct_family("Qt_up","u,c,t");

		/* Up-type antiquark family definition: */

		model<model_t>::construct_family("qbar_up","ubar,cbar");

		/* Up-type antiquark family definition, including b's: */

		model<model_t>::construct_family("Qbar_up","ubar,cbar");

		/* Up-type antiquark family definition, including b's and tops: */

		model<model_t>::construct_family("Qtbar_up","ubar,cbar,tbar");

		/* Down-type quark family definition: */

		model<model_t>::construct_family("q_down","d,s");

		/* Down-type quark family definition, including b's: */

		model<model_t>::construct_family("Q_down","d,s,b");

		/* Down-type quark family definition, including b's and tops: */

		model<model_t>::construct_family("Qt_down","d,s,b");

		/* Down-type antiquark family definition: */

		model<model_t>::construct_family("qbar_down","dbar,sbar");

		/* Down-type antiquark family definition, including b's: */

		model<model_t>::construct_family("Qbar_down","dbar,sbar,bbar");

		/* Down-type antiquark family definition, including b's and tops: */

		model<model_t>::construct_family("Qtbar_down","dbar,sbar,bbar");
		
		/* Positively-charged quark family definition: */

		model<model_t>::construct_family("q+","u,dbar,c,sbar");
		
		/* Positively-charged quark family definition including b's: */

		model<model_t>::construct_family("Q+","u,dbar,c,sbar,bbar");
		
		/* Positively-charged quark family definition including b's and
		 * tops: */

		model<model_t>::construct_family("Qt+","u,dbar,c,sbar,bbar,t");

		/* Negatively-charged quark family definition: */

		model<model_t>::construct_family("q-","d,ubar,s,cbar");

		/* Negatively-charged quark family definition, including b's: */

		model<model_t>::construct_family("Q-","d,ubar,s,cbar,b");

		/* Negatively-charged quark family definition, including b's and
		 * tops: */

		model<model_t>::construct_family("Qt-","d,ubar,s,cbar,b,tbar");

		/* Proton constituent parton family definition: */

		model<model_t>::construct_family("p","u,ubar,d,dbar,c,cbar,s,sbar,g");

		/* Proton constituent parton family definition, including b's: */

		model<model_t>::construct_family("P","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

		/* QCD jet parton family definition: */

		model<model_t>::construct_family("j","u,ubar,d,dbar,c,cbar,s,sbar,g");

		/* QCD jet parton family definition, including b's: */

		model<model_t>::construct_family("J","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

		/* QCD jet parton family definition, including b's and tops: */

		model<model_t>::construct_family("Jt","u,ubar,d,dbar,c,cbar,s,sbar,t,tbar,b,bbar,g");
	    }
	    
	    /* Output of the number of QCD colours: */

	    std::size_t QCD_colours()
	    {
		return model_t::N_c;
	    }

	    /* Refresh quark masses (to use when along the run massless fermions are
	     * to be made massive): */

	    void refresh_fermion_masses()
	    {
		if(model<model_t>::initialised())
		{
		    (M_u>(value_type)0)?model<model_t>::set_mass("u",&M_u):model<model_t>::set_massless("u");
		    (M_d>(value_type)0)?model<model_t>::set_mass("d",&M_d):model<model_t>::set_massless("d");
		    (M_c>(value_type)0)?model<model_t>::set_mass("c",&M_c):model<model_t>::set_massless("c");
		    (M_s>(value_type)0)?model<model_t>::set_mass("s",&M_s):model<model_t>::set_massless("s");
		}
	    }

	    /* Dummy function setting the EM fine-structure constant: */

	    static void set_alpha(const value_type& a){}
	    
	    /* Function setting alpha strong and computing the couplings: */
	    
	    static void set_alpha_s(const value_type& a)
	    {
		alpha_s=a;
		refresh_couplings();
	    }

	    /* Function setting the scale and computing alpha_s (at LO) and the
	     * couplings: */

	    static void set_QCD_scale(const value_type& mu)
	    {
		QCD_scale=mu;
		alpha_s=SM_params::alpha_s/((value_type)1+SM_params::alpha_s*(11*model_t::N_c-12)*std::log(mu/SM_params::QCD_scale)/((value_type)6*pi));
		refresh_couplings();
	    }
	    
	    /* Function computing the couplings from the current alpha_s value:
	     * */
	    
	    static void refresh_couplings()
	    {
		g_s=std::sqrt((value_type)4*pi*alpha_s);
		ggg=std::complex<value_type>(g_s,0);
		gggg=std::complex<value_type>(0,-g_s*g_s);
		Tgg=std::complex<value_type>(0,std::sqrt((value_type)0.5)*g_s);
		guu=std::complex<value_type>(0,g_s);
		gdd=std::complex<value_type>(0,g_s);
		gcc=std::complex<value_type>(0,g_s);
		gss=std::complex<value_type>(0,g_s);
		gtt=std::complex<value_type>(0,g_s);
		gbb=std::complex<value_type>(0,g_s);
	    }

	    /* Function setting the gluon-propagator to the Feynman gauge: */

	    static void set_Feynman_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=0)
		{
		    model<model_t>::template set_propagator<Feynman_gauge>("g");
		}
		gauge=0;
	    }
	    
	    /* Function setting the gluon-propagator to the unitary gauge: */
	    
	    static void set_unitary_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=1)
		{
		    model<model_t>::template set_propagator<unitary_gauge>("g");
		}
		gauge=1;
	    }

	    /* Function setting the gluon-propagator to the R-xi gauge: */

	    static void set_R_xi_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=2)
		{
		    model<model_t>::template set_propagator<R_vector_gauge>("g");
		}
		gauge=2;
	    }
	    
	    /* Function setting the xi-parameter: */
	    
	    static void set_xi(const value_type& x)
	    {
		R_gauge<model_t>::xi=x;
	    }

	    /* Function switching to a 4-gluon vertex Lagrangian: */

	    void set_4_gluon_vertex()
	    {
		if(model<model_t>::initialised() and !four_gluon_vertex)
		{
		    model<model_t>::erase_vertex("H_qcd","g","g");
		    model<model_t>::template add_vertex<colour_tensor::ff_contr< SU<model_t::N_c> >,vvvv>("g","g","g","g",&gggg);
		}
		four_gluon_vertex=true;
	    }

	    /* Function switching to an auxiliary tensor field description of
	     * the 4-gluon vertex: */

	    void set_auxiliary_QCD_field()
	    {
		if(model<model_t>::initialised() and four_gluon_vertex)
		{
		    model<model_t>::erase_vertex("g","g","g","g");
		    model<model_t>::template add_fast_4g_vertex< SU<model_t::N_c> >("g",&Tgg);
		}
		four_gluon_vertex=false;
	    }

	    /* Reset all parameters using the default input parameters in the
	     * SM_params.h file: */

	    void set_default_params()
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
	
	private:

	    /* Gauge tag: */

	    static int gauge;

	    /* 4-gluon vertex tag (true: model constructs a 4-gluon vertex,
	     * false: model uses the antisymmetric tensor field): */

	    static bool four_gluon_vertex;
    };

    /* Initialisation of input parameters: */

    template<class model_t,class value_t>bool QCD_base<model_t,value_t>::four_gluon_vertex=false;
    template<class model_t,class value_t>int QCD_base<model_t,value_t>::gauge=0;
    template<class model_t,class value_t>const value_t QCD_base<model_t,value_t>::alpha=-(value_t)1;
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::alpha_s=(value_t)SM_params::alpha_s;
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::QCD_scale=(value_t)SM_params::QCD_scale;
    template<class model_t,class value_t>const value_t QCD_base<model_t,value_t>::pi=std::acos(-(value_t)1);
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::g_s=std::sqrt((value_t)4*QCD_base<model_t,value_t>::pi*QCD_base<model_t,value_t>::alpha_s);
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::M_u=(value_t)0;
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::M_d=(value_t)0;
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::M_c=(value_t)SM_params::M_c;
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::M_s=(value_t)0;
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::M_t=(value_t)SM_params::M_t;
    template<class model_t,class value_t>value_t QCD_base<model_t,value_t>::M_b=(value_t)SM_params::M_b;

    /* Initialisation values of the couplings: */

    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::ggg(QCD_base<model_t,value_t>::g_s,0);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::gggg(0,-QCD_base<model_t,value_t>::g_s*QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::Tgg(0,std::sqrt((value_t)0.5)*QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::guu(0,QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::gdd(0,QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::gcc(0,QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::gss(0,QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::gtt(0,QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>QCD_base<model_t,value_t>::gbb(0,QCD_base<model_t,value_t>::g_s);
}

#endif /*CAMGEN_QCD_BASE_H_*/

