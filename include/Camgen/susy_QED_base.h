//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SUSY_QED_BASE_H_
#define CAMGEN_SUSY_QED_BASE_H_

#include <Camgen/model.h>
#include <Camgen/vector_particle.h>
#include <Camgen/scalar_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/vff.h>
#include <Camgen/vss.h>
#include <Camgen/ssvv.h>
#include <Camgen/sff.h>
#include <Camgen/sffL.h>
#include <Camgen/sffR.h>
#include <Camgen/ssss.h>
#include <Camgen/SM_params.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the susy_QED_base class template in Camgen.    *
 * Publicly deriving a model class model_t from susy_QED_base<model_t> amounts   *
 * to adding supersymmetric QED to the model. It should be noted that the value  *
 * type of base_QED has to be castable to that of the derived model.             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */                                                                            

namespace Camgen
{
    template<class model_t,class value_t>class susy_QED_base: public model<model_t>
    {
	public:

	    typedef value_t value_type;
	    
	    /* Electromagnetic fine structure constant: */

	    static value_type alpha;

	    /* Strong coupling constant (not used): */

	    static const value_type alpha_s;

	    /* QCD scale parameter (not used): */

	    static const value_type QCD_scale;

	    /* Unit charge: */

	    static value_type Q_e;

	    /* Pi: */

	    static const value_type pi;

	    /* Masses: */

	    static value_type M_e,M_mu,M_tau,M_seL,M_seR,M_smuL,M_smuR,M_stauL,M_stauR,M_chi;

	    /* Widths: */

	    static value_type W_seL,W_seR,W_smuL,W_smuR,W_stauL,W_stauR,W_chi;

	    /* Photon-fermion-fermion couplings: */

	    static std::complex<value_type>gammaee,gammamumu,gammatautau;

	    /* Photon-slepton-slepton couplings: */

	    static std::complex<value_type>gammaseLseL,gammaseRseR,gammasmuLsmuL,gammasmuRsmuR,gammastauLstauL,gammastauRstauR;

	    /* Slepton-photino-lepton couplings: */

	    static std::complex<value_type>seLchie,seRchie,smuLchimu,smuRchimu,stauLchitau,stauRchitau;

	    /* Antislepton-antilepton-photino couplings: */

	    static std::complex<value_type>seLechi,seRechi,smuLmuchi,smuRmuchi,stauLtauchi,stauRtauchi;
	    
	    /* Slepton-slepton-photon-photon coupling: */

	    static std::complex<value_type>seLseLgammagamma,seRseRgammagamma,smuLsmuLgammagamma,smuRsmuRgammagamma,stauLstauLgammagamma,stauRstauRgammagamma;
	    
	    /* Quartic slepton couplings: */
	    
	    static std::complex<value_type>seLseLseLseL,seRseRseRseR,smuLsmuLsmuLsmuL,smuRsmuRsmuRsmuR,stauLstauLstauLstauL,stauRstauRstauRstauR;
	    static std::complex<value_type>seLseLseRseR,smuLsmuLsmuRsmuR,stauLstauLstauRstauR;

	    /* Constructor: */

	    susy_QED_base()
	    {
		/* Photon definition: */

		switch(gauge)
		{
		    case 0:
			model<model_t>::template add_vector<Feynman_gauge>("gamma",22);
			break;
		    case 1:
			model<model_t>::template add_vector<unitary_gauge>("gamma",22);
			break;
		    case 2:
			model<model_t>::template add_vector<R_vector_gauge>("gamma",22);
			break;
		    default:
			model<model_t>::template add_vector<Feynman_gauge>("gamma",22);
		}

		/* Lepton definitions: */

		(M_e>(value_type)0)?model<model_t>::add_fermions("e-","e+",&M_e,11):model<model_t>::add_fermions("e-","e+",11);
		(M_mu>(value_type)0)?model<model_t>::add_fermions("mu-","mu+",&M_mu,13):model<model_t>::add_fermions("mu-","mu+",13);
		(M_tau>(value_type)0)?model<model_t>::add_fermions("tau-","tau+",&M_tau,15):model<model_t>::add_fermions("tau-","tau+",15);

		/* Photino definition: */

		(M_chi>(value_type)0)?model<model_t>::add_fermion("~chi_10",&M_chi,&W_chi,1000022):model<model_t>::add_fermion("~chi_10",1000022);

		/* Slepton definitions: */

		(M_seL>(value_type)0)?model<model_t>::add_scalars("~e_L-","~e_L+",&M_seL,&W_seL,1000011):model<model_t>::add_scalars("~e_L-","~e_L+",1000011);
		(M_seR>(value_type)0)?model<model_t>::add_scalars("~e_R-","~e_R+",&M_seR,&W_seR,2000011):model<model_t>::add_scalars("~e_R-","~e_R+",2000011);
		(M_smuL>(value_type)0)?model<model_t>::add_scalars("~mu_L-","~mu_L+",&M_smuL,&W_smuL,1000013):model<model_t>::add_scalars("~mu_L-","~mu_L+",1000013);
		(M_smuR>(value_type)0)?model<model_t>::add_scalars("~mu_R-","~mu_R+",&M_smuR,&W_smuR,2000013):model<model_t>::add_scalars("~mu_R-","~mu_R+",2000013);
		(M_stauL>(value_type)0)?model<model_t>::add_scalars("~tau_L-","~tau_L+",&M_stauL,&W_stauL,1000015):model<model_t>::add_scalars("~tau_L-","~tau_L+",1000015);
		(M_stauR>(value_type)0)?model<model_t>::add_scalars("~tau_R-","~tau_R+",&M_stauR,&W_stauR,2000015):model<model_t>::add_scalars("~tau_R-","~tau_R+",2000015);

		/* Vertices: */

		model<model_t>::template add_vertex<vff>("gamma","e+","e-",&gammaee);
		model<model_t>::template add_vertex<vss>("gamma","~e_L-","~e_L+",&gammaseLseL);
		model<model_t>::template add_vertex<vss>("gamma","~e_R-","~e_R+",&gammaseRseR);
		model<model_t>::template add_vertex<sffR>("~e_L+","~chi_10","e-",&seLchie);
		model<model_t>::template add_vertex<sffL>("~e_L-","e+","~chi_10",&seLechi);
		model<model_t>::template add_vertex<sffL>("~e_R+","~chi_10","e-",&seRchie);
		model<model_t>::template add_vertex<sffR>("~e_R-","e+","~chi_10",&seRechi);
		model<model_t>::template add_vertex<ssvv>("~e_L-","~e_L+","gamma","gamma",&seLseLgammagamma);
		model<model_t>::template add_vertex<ssvv>("~e_R+","~e_R-","gamma","gamma",&seRseRgammagamma);
		model<model_t>::template add_vertex<ssss>("~e_L-","~e_L+","~e_L-","~e_L+",&seLseLseLseL);
		model<model_t>::template add_vertex<ssss>("~e_R+","~e_R-","~e_R+","~e_R-",&seRseRseRseR);
		model<model_t>::template add_vertex<ssss>("~e_L-","~e_L+","~e_R+","~e_R-",&seLseLseRseR);

		model<model_t>::template add_vertex<vff>("gamma","mu+","mu-",&gammamumu);
		model<model_t>::template add_vertex<vss>("gamma","~mu_L-","~mu_L+",&gammasmuLsmuL);
		model<model_t>::template add_vertex<vss>("gamma","~mu_R-","~mu_R+",&gammasmuRsmuR);
		model<model_t>::template add_vertex<sffR>("~mu_L+","~chi_10","mu-",&smuLchimu);
		model<model_t>::template add_vertex<sffL>("~mu_L-","mu+","~chi_10",&smuLmuchi);
		model<model_t>::template add_vertex<sffL>("~mu_R+","~chi_10","mu-",&smuRchimu);
		model<model_t>::template add_vertex<sffR>("~mu_R-","mu+","~chi_10",&smuRmuchi);
		model<model_t>::template add_vertex<ssvv>("~mu_L-","~mu_L+","gamma","gamma",&smuLsmuLgammagamma);
		model<model_t>::template add_vertex<ssvv>("~mu_R+","~mu_R-","gamma","gamma",&smuRsmuRgammagamma);
		model<model_t>::template add_vertex<ssss>("~mu_L-","~mu_L+","~mu_L-","~mu_L+",&smuLsmuLsmuLsmuL);
		model<model_t>::template add_vertex<ssss>("~mu_R+","~mu_R-","~mu_R+","~mu_R-",&smuRsmuRsmuRsmuR);
		model<model_t>::template add_vertex<ssss>("~mu_L-","~mu_L+","~mu_R+","~mu_R-",&smuLsmuLsmuRsmuR);

		model<model_t>::template add_vertex<vff>("gamma","tau+","tau-",&gammatautau);
		model<model_t>::template add_vertex<vss>("gamma","~tau_L-","~tau_L+",&gammastauLstauL);
		model<model_t>::template add_vertex<vss>("gamma","~tau_R-","~tau_R+",&gammastauRstauR);
		model<model_t>::template add_vertex<sffR>("~tau_L+","~chi_10","tau-",&stauLchitau);
		model<model_t>::template add_vertex<sffL>("~tau_L-","tau+","~chi_10",&stauLtauchi);
		model<model_t>::template add_vertex<sffL>("~tau_R+","~chi_10","tau-",&stauRchitau);
		model<model_t>::template add_vertex<sffR>("~tau_R-","tau+","~chi_10",&stauRtauchi);
		model<model_t>::template add_vertex<ssvv>("~tau_L-","~tau_L+","gamma","gamma",&stauLstauLgammagamma);
		model<model_t>::template add_vertex<ssvv>("~tau_R+","~tau_R-","gamma","gamma",&stauRstauRgammagamma);
		model<model_t>::template add_vertex<ssss>("~tau_L-","~tau_L+","~tau_L-","~tau_L+",&stauLstauLstauLstauL);
		model<model_t>::template add_vertex<ssss>("~tau_R+","~tau_R-","~tau_R+","~tau_R-",&stauRstauRstauRstauR);
		model<model_t>::template add_vertex<ssss>("~tau_L-","~tau_L+","~tau_R+","~tau_R-",&stauLstauLstauRstauR);
		
		/* Positively-charged lepton family definition: */

		model<model_t>::construct_family("l+","e+,mu+");

		/* Positively-charged lepton family definition, including taus: */

		model<model_t>::construct_family("L+","e+,mu+,tau+");

		/* Negatively-charged lepton family definition: */

		model<model_t>::construct_family("l-","e-,mu-");

		/* Negatively-charged lepton family definition, including taus: */

		model<model_t>::construct_family("L-","e-,mu-,tau-");

		/* Left-handed positive slepton family definitions: */

		model<model_t>::construct_family("~L_L+","~e_L+,~mu_L+,~tau_L+");

		/* Left-handed negative slepton family definitions: */

		model<model_t>::construct_family("~L_L-","~e_L-,~mu_L-,~tau_L-");

		/* Right-handed positive slepton family definitions: */

		model<model_t>::construct_family("~L_R+","~e_R+,~mu_R+,~tau_R+");

		/* Right-handed negative slepton family definitions: */

		model<model_t>::construct_family("~L_R-","~e_R-,~mu_R-,~tau_R-");

		/* Positive slepton family definitions: */

		model<model_t>::construct_family("~L+","~e_L+,~e_R+,~mu_L+,~mu_R+,~tau_L+,~tau_R+");

		/* Negative slepton family definitions: */

		model<model_t>::construct_family("~L-","~e_L-,~e_R-,~mu_L-,~mu_R-,~tau_L-,~tau_R-");

		/* Left-handed slepton family definitions: */

		model<model_t>::construct_family("~L_L","~e_L-,~e_L+,~mu_L-,~mu_L+,~tau_L-,~tau_L+");

		/* Right-handed slepton family definitions: */

		model<model_t>::construct_family("~L_R","~e_R-,~e_R+,~mu_R-,~mu_R+,~tau_R-,~tau_R+");
	    }

	    /* Function computing the couplings from the input parameters: */

	    static void refresh_couplings()
	    {
		Q_e=std::sqrt((value_type)4*pi*alpha);

		gammaee=std::complex<value_type>(0,Q_e);
		gammamumu=gammaee;
		gammatautau=gammaee;

		gammaseLseL=std::complex<value_type>(0,Q_e);
		gammaseRseR=gammaseLseL;

		gammasmuLsmuL=gammaseLseL;
		gammasmuRsmuR=gammaseLseL;

		gammastauLstauL=gammaseLseL;
		gammastauRstauR=gammaseLseL;

		seLchie=std::complex<value_type>(0,std::sqrt((value_type)2)*Q_e);
		seRchie=seLchie;
		seLechi=-seLchie;
		seRechi=-seLchie;

		smuLchimu=seLchie;
		smuRchimu=seLchie;
		smuLmuchi=-seLchie;
		smuRmuchi=-seLchie;

		stauLchitau=seLchie;
		stauRchitau=seLchie;
		stauLtauchi=-seLchie;
		stauRtauchi=-seLchie;


		seLseLgammagamma=std::complex<value_type>(0,(value_type)2*Q_e*Q_e);
		seRseRgammagamma=seLseLgammagamma;
		smuLsmuLgammagamma=seLseLgammagamma;
		smuRsmuRgammagamma=seLseLgammagamma;
		stauLstauLgammagamma=seLseLgammagamma;
		stauRstauRgammagamma=seLseLgammagamma;

		seLseLseLseL=-seLseLgammagamma;
		seRseRseRseR=seLseLseLseL;
		smuLsmuLsmuLsmuL=seLseLseLseL;
		smuRsmuRsmuRsmuR=seLseLseLseL;
		stauLstauLstauLstauL=seLseLseLseL;
		stauRstauRstauRstauR=seLseLseLseL;

		seLseLseRseR=(value_type)0.5*seLseLgammagamma;
		smuLsmuLsmuRsmuR=seLseLseRseR;
		stauLstauLstauRstauR=seLseLseRseR;
	    }

	    /* Function refreshing the particle masses: */

	    static void refresh_masses()
	    {
		if(model<model_t>::initialised())
		{
		    (M_e>(value_type)0)?model<model_t>::set_mass("e-",&M_e):model<model_t>::set_massless("e-");
		    (M_mu>(value_type)0)?model<model_t>::set_mass("mu-",&M_mu):model<model_t>::set_massless("mu-");
		    (M_tau>(value_type)0)?model<model_t>::set_mass("tau-",&M_tau):model<model_t>::set_massless("tau-");
		    (M_chi>(value_type)0)?model<model_t>::set_mass("~chi_10",&M_chi):model<model_t>::set_massless("~chi_10");
		    (M_seL>(value_type)0)?model<model_t>::set_mass("~e_L-",&M_seL):model<model_t>::set_massless("~e_L-");
		    (M_seR>(value_type)0)?model<model_t>::set_mass("~e_R-",&M_seR):model<model_t>::set_massless("~e_R-");
		    (M_smuL>(value_type)0)?model<model_t>::set_mass("~mu_L-",&M_smuL):model<model_t>::set_massless("~mu_L-");
		    (M_smuR>(value_type)0)?model<model_t>::set_mass("~mu_R-",&M_smuR):model<model_t>::set_massless("~mu_R-");
		    (M_stauL>(value_type)0)?model<model_t>::set_mass("~tau_L-",&M_stauL):model<model_t>::set_massless("~tau_L-");
		    (M_stauR>(value_type)0)?model<model_t>::set_mass("~tau_R-",&M_stauR):model<model_t>::set_massless("~tau_R-");
		}
	    }

	    /* Function refreshing the particle decay widths: */

	    static void refresh_widths()
	    {
		if(model<model_t>::initialised())
		{
		    (W_chi>(value_type)0)?model<model_t>::set_width("~gamma",&M_chi):model<model_t>::set_widthless("~gamma");
		    (W_seL>(value_type)0)?model<model_t>::set_width("~e_L-",&M_seL):model<model_t>::set_widthless("~e_L-");
		    (W_seR>(value_type)0)?model<model_t>::set_width("~e_R-",&M_seR):model<model_t>::set_widthless("~e_R-");
		    (W_smuL>(value_type)0)?model<model_t>::set_width("~mu_L-",&M_smuL):model<model_t>::set_widthless("~mu_L-");
		    (W_smuR>(value_type)0)?model<model_t>::set_width("~mu_R-",&M_smuR):model<model_t>::set_widthless("~mu_R-");
		    (W_stauL>(value_type)0)?model<model_t>::set_width("~tau_L-",&M_stauL):model<model_t>::set_widthless("~tau_L-");
		    (W_stauR>(value_type)0)?model<model_t>::set_width("~tau_R-",&M_stauR):model<model_t>::set_widthless("~tau_R-");
		}
	    }

	    /* Function returning the QCD colours: */

	    static std::size_t QCD_colours()
	    {
		return 1;
	    }

	    /* EM fine structure input: */

	    static void set_alpha(const value_type& a)
	    {
		alpha=a;
		refresh_couplings();
	    }

	    /* Function setting the strong coupling: */

	    static void set_alpha_s(const value_type& a){}

	    /* Function setting the QCD scale: */

	    static void set_QCD_scale(const value_type& mu){}

	    /* Function setting the photon propagator to the unitary gauge: */

	    static void set_unitary_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=1)
		{
		    model<model_t>::template set_propagator<unitary_gauge>("gamma");
		}
		gauge=1;
	    }

	    /* Function setting the photon propagator to the Feynman gauge: */

	    static void set_Feynman_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=0)
		{
		    model<model_t>::template set_propagator<Feynman_gauge>("gamma");
		}
		gauge=0;
	    }

	    /* Function setting the photon propagator to the R-xi gauge: */

	    static void set_R_xi_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=2)
		{
		    model<model_t>::template set_propagator<R_vector_gauge>("gamma");
		}
		gauge=2;
	    }

	    /* Function setting the gauge parameter: */

	    static void set_xi(const value_type& x)
	    {
		R_gauge<model_t>::xi=x;
	    }

	private:

	    /* Photon gauge switch: */

	    static int gauge;
    };

    /* Const static unused constants: */

    template<class model_t,class value_t>const value_t susy_QED_base<model_t,value_t>::alpha_s=-(value_t)1;
    template<class model_t,class value_t>const value_t susy_QED_base<model_t,value_t>::QCD_scale=-(value_t)1;

    /* Real computed parameter initialisations: */

    template<class model_t,class value_t>const value_t susy_QED_base<model_t,value_t>::pi=std::acos(-(value_t)1);
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::alpha=(value_t)SM_params::alpha;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::Q_e=std::sqrt((value_t)4*susy_QED_base<model_t,value_t>::pi*susy_QED_base<model_t,value_t>::alpha);

    /* Particle masses and widths: */

    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_e=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_mu=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_tau=(value_t)SM_params::M_tau;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_seL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_seR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_smuL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_smuR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_stauL=(value_t)SM_params::M_tau;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_stauR=(value_t)SM_params::M_tau;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::M_chi=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::W_seL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::W_seR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::W_smuL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::W_smuR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::W_stauL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::W_stauR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QED_base<model_t,value_t>::W_chi=(value_t)0;
    
    /* Coupling constant initialisations: */

    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammaee(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammamumu(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammatautau(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammaseLseL(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammaseRseR(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammasmuLsmuL(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammasmuRsmuR(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammastauLstauL(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::gammastauRstauR(0,susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seLchie(0,std::sqrt((value_t)2)*susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seRchie=susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seLechi=-susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seRechi=-susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuLchimu=susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuRchimu=susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuLmuchi=-susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuRmuchi=-susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauLchitau=susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauRchitau=susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauLtauchi=-susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauRtauchi=-susy_QED_base<model_t,value_t>::seLchie;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seLseLgammagamma(0,(value_t)2*susy_QED_base<model_t,value_t>::Q_e*susy_QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seRseRgammagamma=susy_QED_base<model_t,value_t>::seLseLgammagamma;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuLsmuLgammagamma=susy_QED_base<model_t,value_t>::seLseLgammagamma;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuRsmuRgammagamma=susy_QED_base<model_t,value_t>::seLseLgammagamma;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauLstauLgammagamma=susy_QED_base<model_t,value_t>::seLseLgammagamma;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauRstauRgammagamma=susy_QED_base<model_t,value_t>::seLseLgammagamma;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seLseLseLseL=-susy_QED_base<model_t,value_t>::seLseLgammagamma;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seRseRseRseR=susy_QED_base<model_t,value_t>::seLseLseLseL;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuLsmuLsmuLsmuL=susy_QED_base<model_t,value_t>::seLseLseLseL;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuRsmuRsmuRsmuR=susy_QED_base<model_t,value_t>::seLseLseLseL;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauLstauLstauLstauL=susy_QED_base<model_t,value_t>::seLseLseLseL;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauRstauRstauRstauR=susy_QED_base<model_t,value_t>::seLseLseLseL;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::seLseLseRseR=(value_t)0.5*susy_QED_base<model_t,value_t>::seLseLgammagamma;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::smuLsmuLsmuRsmuR=susy_QED_base<model_t,value_t>::seLseLseRseR;
    template<class model_t,class value_t>std::complex<value_t>susy_QED_base<model_t,value_t>::stauLstauLstauRstauR=susy_QED_base<model_t,value_t>::seLseLseRseR;

    /* Gauge tag initialisation (to Feynman gauge): */

    template<class model_t,class value_t>int susy_QED_base<model_t,value_t>::gauge=0;
}

#endif /*CAMGEN_SUSY_QED_BASE_H_*/

