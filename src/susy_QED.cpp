//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/susy_QED.h>
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

namespace Camgen
{
    /* Compile-time data: */

    const std::size_t susy_QED::dimension;
    const bool susy_QED::coloured;
    const bool susy_QED::continuous_helicities;
    const int susy_QED::beam_direction;

    /* Static unused constants: */

    const susy_QED::value_type susy_QED::alpha_s=-(susy_QED::value_type)1;
    const susy_QED::value_type susy_QED::QCD_scale=-(susy_QED::value_type)1;

    /* Real computed parameter initialisations: */

    susy_QED::value_type susy_QED::alpha=(susy_QED::value_type)SM_params::alpha;
    const susy_QED::value_type susy_QED::pi=std::acos(-(susy_QED::value_type)1);
    susy_QED::value_type susy_QED::Q_e=std::sqrt((susy_QED::value_type)4*susy_QED::pi*susy_QED::alpha);

    /* Particle masses and widths: */

    susy_QED::value_type susy_QED::M_e=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::M_mu=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::M_tau=(susy_QED::value_type)SM_params::M_tau;
    susy_QED::value_type susy_QED::M_seL=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::M_seR=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::M_smuL=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::M_smuR=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::M_stauL=(susy_QED::value_type)SM_params::M_tau;
    susy_QED::value_type susy_QED::M_stauR=(susy_QED::value_type)SM_params::M_tau;
    susy_QED::value_type susy_QED::M_chi=(susy_QED::value_type)0;

    susy_QED::value_type susy_QED::W_seL=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::W_seR=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::W_smuL=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::W_smuR=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::W_stauL=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::W_stauR=(susy_QED::value_type)0;
    susy_QED::value_type susy_QED::W_chi=(susy_QED::value_type)0;
    
    /* Coupling constant initialisations: */

    std::complex<susy_QED::value_type>susy_QED::gammaee(0,susy_QED::Q_e);
    std::complex<susy_QED::value_type>susy_QED::gammamumu(0,susy_QED::Q_e);
    std::complex<susy_QED::value_type>susy_QED::gammatautau(0,susy_QED::Q_e);
    
    std::complex<susy_QED::value_type>susy_QED::gammaseLseL(0,susy_QED::Q_e);
    std::complex<susy_QED::value_type>susy_QED::gammaseRseR(0,susy_QED::Q_e);
    
    std::complex<susy_QED::value_type>susy_QED::gammasmuLsmuL(0,susy_QED::Q_e);
    std::complex<susy_QED::value_type>susy_QED::gammasmuRsmuR(0,susy_QED::Q_e);
    
    std::complex<susy_QED::value_type>susy_QED::gammastauLstauL(0,susy_QED::Q_e);
    std::complex<susy_QED::value_type>susy_QED::gammastauRstauR(0,susy_QED::Q_e);
    
    std::complex<susy_QED::value_type>susy_QED::seLchie(0,std::sqrt((susy_QED::value_type)2)*susy_QED::Q_e);
    std::complex<susy_QED::value_type>susy_QED::seRchie=susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::seLechi=-susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::seRechi=-susy_QED::seLchie;

    std::complex<susy_QED::value_type>susy_QED::smuLchimu=susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::smuRchimu=susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::smuLmuchi=-susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::smuRmuchi=-susy_QED::seLchie;
    
    std::complex<susy_QED::value_type>susy_QED::stauLchitau=susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::stauRchitau=susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::stauLtauchi=-susy_QED::seLchie;
    std::complex<susy_QED::value_type>susy_QED::stauRtauchi=-susy_QED::seLchie;

    
    std::complex<susy_QED::value_type>susy_QED::seLseLgammagamma(0,(susy_QED::value_type)2*susy_QED::Q_e*susy_QED::Q_e);
    std::complex<susy_QED::value_type>susy_QED::seRseRgammagamma=susy_QED::seLseLgammagamma;
    std::complex<susy_QED::value_type>susy_QED::smuLsmuLgammagamma=susy_QED::seLseLgammagamma;
    std::complex<susy_QED::value_type>susy_QED::smuRsmuRgammagamma=susy_QED::seLseLgammagamma;
    std::complex<susy_QED::value_type>susy_QED::stauLstauLgammagamma=susy_QED::seLseLgammagamma;
    std::complex<susy_QED::value_type>susy_QED::stauRstauRgammagamma=susy_QED::seLseLgammagamma;

    std::complex<susy_QED::value_type>susy_QED::seLseLseLseL=-susy_QED::seLseLgammagamma;
    std::complex<susy_QED::value_type>susy_QED::seRseRseRseR=susy_QED::seLseLseLseL;
    std::complex<susy_QED::value_type>susy_QED::smuLsmuLsmuLsmuL=susy_QED::seLseLseLseL;
    std::complex<susy_QED::value_type>susy_QED::smuRsmuRsmuRsmuR=susy_QED::seLseLseLseL;
    std::complex<susy_QED::value_type>susy_QED::stauLstauLstauLstauL=susy_QED::seLseLseLseL;
    std::complex<susy_QED::value_type>susy_QED::stauRstauRstauRstauR=susy_QED::seLseLseLseL;
    
    std::complex<susy_QED::value_type>susy_QED::seLseLseRseR=(susy_QED::value_type)0.5*susy_QED::seLseLgammagamma;
    std::complex<susy_QED::value_type>susy_QED::smuLsmuLsmuRsmuR=susy_QED::seLseLseRseR;
    std::complex<susy_QED::value_type>susy_QED::stauLstauLstauRstauR=susy_QED::seLseLseRseR;

    /* Gauge tag initialisation (to Feynman gauge): */

    int susy_QED::gauge=0;

    /* Constructor: */

    susy_QED::susy_QED()
    {
	/* Photon definition: */

	switch(gauge)
	{
	    case 0:
		add_vector<Feynman_gauge>("gamma",22);
		break;
	    case 1:
		add_vector<unitary_gauge>("gamma",22);
		break;
	    case 2:
		add_vector<R_vector_gauge>("gamma",22);
		break;
	    default:
		add_vector<Feynman_gauge>("gamma",22);
	}

	/* Lepton definitions: */

	(M_e>(value_type)0)?add_fermions("e-","e+",&M_e,11):add_fermions("e-","e+",11);
	(M_mu>(value_type)0)?add_fermions("mu-","mu+",&M_mu,13):add_fermions("mu-","mu+",13);
	(M_tau>(value_type)0)?add_fermions("tau-","tau+",&M_tau,15):add_fermions("tau-","tau+",15);

	/* Photino definition: */

	(M_chi>(value_type)0)?add_fermion("~chi_10",&M_chi,&W_chi,1000022):add_fermion("~chi_10",1000022);

	/* Slepton definitions: */

	(M_seL>(value_type)0)?add_scalars("~e_L-","~e_L+",&M_seL,&W_seL,1000011):add_scalars("~e_L-","~e_L+",1000011);
	(M_seR>(value_type)0)?add_scalars("~e_R-","~e_R+",&M_seR,&W_seR,2000011):add_scalars("~e_R-","~e_R+",2000011);
	(M_smuL>(value_type)0)?add_scalars("~mu_L-","~mu_L+",&M_smuL,&W_smuL,1000013):add_scalars("~mu_L-","~mu_L+",1000013);
	(M_smuR>(value_type)0)?add_scalars("~mu_R-","~mu_R+",&M_smuR,&W_smuR,2000013):add_scalars("~mu_R-","~mu_R+",2000013);
	(M_stauL>(value_type)0)?add_scalars("~tau_L-","~tau_L+",&M_stauL,&W_stauL,1000015):add_scalars("~tau_L-","~tau_L+",1000015);
	(M_stauR>(value_type)0)?add_scalars("~tau_R-","~tau_R+",&M_stauR,&W_stauR,2000015):add_scalars("~tau_R-","~tau_R+",2000015);

	/* Vertices: */

	add_vertex<vff>("gamma","e+","e-",&gammaee);
	add_vertex<vss>("gamma","~e_L-","~e_L+",&gammaseLseL);
	add_vertex<vss>("gamma","~e_R-","~e_R+",&gammaseRseR);
	add_vertex<sffR>("~e_L+","~chi_10","e-",&seLchie);
	add_vertex<sffL>("~e_L-","e+","~chi_10",&seLechi);
	add_vertex<sffL>("~e_R+","~chi_10","e-",&seRchie);
	add_vertex<sffR>("~e_R-","e+","~chi_10",&seRechi);
	add_vertex<ssvv>("~e_L-","~e_L+","gamma","gamma",&seLseLgammagamma);
	add_vertex<ssvv>("~e_R+","~e_R-","gamma","gamma",&seRseRgammagamma);
	add_vertex<ssss>("~e_L-","~e_L+","~e_L-","~e_L+",&seLseLseLseL);
	add_vertex<ssss>("~e_R+","~e_R-","~e_R+","~e_R-",&seRseRseRseR);
	add_vertex<ssss>("~e_L-","~e_L+","~e_R+","~e_R-",&seLseLseRseR);
	
	add_vertex<vff>("gamma","mu+","mu-",&gammamumu);
	add_vertex<vss>("gamma","~mu_L-","~mu_L+",&gammasmuLsmuL);
	add_vertex<vss>("gamma","~mu_R-","~mu_R+",&gammasmuRsmuR);
	add_vertex<sffR>("~mu_L+","~chi_10","mu-",&smuLchimu);
	add_vertex<sffL>("~mu_L-","mu+","~chi_10",&smuLmuchi);
	add_vertex<sffL>("~mu_R+","~chi_10","mu-",&smuRchimu);
	add_vertex<sffR>("~mu_R-","mu+","~chi_10",&smuRmuchi);
	add_vertex<ssvv>("~mu_L-","~mu_L+","gamma","gamma",&smuLsmuLgammagamma);
	add_vertex<ssvv>("~mu_R+","~mu_R-","gamma","gamma",&smuRsmuRgammagamma);
	add_vertex<ssss>("~mu_L-","~mu_L+","~mu_L-","~mu_L+",&smuLsmuLsmuLsmuL);
	add_vertex<ssss>("~mu_R+","~mu_R-","~mu_R+","~mu_R-",&smuRsmuRsmuRsmuR);
	add_vertex<ssss>("~mu_L-","~mu_L+","~mu_R+","~mu_R-",&smuLsmuLsmuRsmuR);

	add_vertex<vff>("gamma","tau+","tau-",&gammatautau);
	add_vertex<vss>("gamma","~tau_L-","~tau_L+",&gammastauLstauL);
	add_vertex<vss>("gamma","~tau_R-","~tau_R+",&gammastauRstauR);
	add_vertex<sffR>("~tau_L+","~chi_10","tau-",&stauLchitau);
	add_vertex<sffL>("~tau_L-","tau+","~chi_10",&stauLtauchi);
	add_vertex<sffL>("~tau_R+","~chi_10","tau-",&stauRchitau);
	add_vertex<sffR>("~tau_R-","tau+","~chi_10",&stauRtauchi);
	add_vertex<ssvv>("~tau_L-","~tau_L+","gamma","gamma",&stauLstauLgammagamma);
	add_vertex<ssvv>("~tau_R+","~tau_R-","gamma","gamma",&stauRstauRgammagamma);
	add_vertex<ssss>("~tau_L-","~tau_L+","~tau_L-","~tau_L+",&stauLstauLstauLstauL);
	add_vertex<ssss>("~tau_R+","~tau_R-","~tau_R+","~tau_R-",&stauRstauRstauRstauR);
	add_vertex<ssss>("~tau_L-","~tau_L+","~tau_R+","~tau_R-",&stauLstauLstauRstauR);

	/* Positively-charged lepton family definition: */

	construct_family("l+","e+,mu+");

	/* Positively-charged lepton family definition, including taus: */

	construct_family("L+","e+,mu+,tau+");

	/* Negatively-charged lepton family definition: */

	construct_family("l-","e-,mu-");

	/* Negatively-charged lepton family definition, including taus: */

	construct_family("L-","e-,mu-,tau-");

	/* Left-handed positive slepton family definitions: */

	construct_family("~L_L+","~e_L+,~mu_L+,~tau_L+");

	/* Left-handed negative slepton family definitions: */

	construct_family("~L_L-","~e_L-,~mu_L-,~tau_L-");

	/* Right-handed positive slepton family definitions: */

	construct_family("~L_R+","~e_R+,~mu_R+,~tau_R+");

	/* Right-handed negative slepton family definitions: */

	construct_family("~L_R-","~e_R-,~mu_R-,~tau_R-");

	/* Positive slepton family definitions: */

	construct_family("~L+","~e_L+,~e_R+,~mu_L+,~mu_R+,~tau_L+,~tau_R+");

	/* Negative slepton family definitions: */

	construct_family("~L-","~e_L-,~e_R-,~mu_L-,~mu_R-,~tau_L-,~tau_R-");

	/* Left-handed slepton family definitions: */

	construct_family("~L_L","~e_L-,~e_L+,~mu_L-,~mu_L+,~tau_L-,~tau_L+");

	/* Right-handed slepton family definitions: */

	construct_family("~L_R","~e_R-,~e_R+,~mu_R-,~mu_R+,~tau_R-,~tau_R+");
    }

    /* Function computing the couplings from the input parameters: */

    void susy_QED::refresh_couplings()
    {
	Q_e=std::sqrt((value_type)4*pi*alpha);

	gammaee=std::complex<value_type>(0,susy_QED::Q_e);
	gammamumu=gammaee;
	gammatautau=gammaee;

	gammaseLseL=std::complex<value_type>(0,susy_QED::Q_e);
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

    void susy_QED::refresh_masses()
    {
	if(initialised())
	{
	    (M_e>(value_type)0)?set_mass("e-",&M_e):set_massless("e-");
	    (M_mu>(value_type)0)?set_mass("mu-",&M_mu):set_massless("mu-");
	    (M_tau>(value_type)0)?set_mass("tau-",&M_tau):set_massless("tau-");
	    (M_chi>(value_type)0)?set_mass("~chi_10",&M_chi):set_massless("~chi_10");
	    (M_seL>(value_type)0)?set_mass("~e_L-",&M_seL):set_massless("~e_L-");
	    (M_seR>(value_type)0)?set_mass("~e_R-",&M_seR):set_massless("~e_R-");
	    (M_smuL>(value_type)0)?set_mass("~mu_L-",&M_smuL):set_massless("~mu_L-");
	    (M_smuR>(value_type)0)?set_mass("~mu_R-",&M_smuR):set_massless("~mu_R-");
	    (M_stauL>(value_type)0)?set_mass("~tau_L-",&M_stauL):set_massless("~tau_L-");
	    (M_stauR>(value_type)0)?set_mass("~tau_R-",&M_stauR):set_massless("~tau_R-");
	}
    }

    /* Function refreshing the particle decay widths: */

    void susy_QED::refresh_widths()
    {
	if(initialised())
	{
	    (W_chi>(value_type)0)?set_width("~chi_10",&M_chi):set_widthless("~chi_10");
	    (W_seL>(value_type)0)?set_width("~e_L-",&M_seL):set_widthless("~e_L-");
	    (W_seR>(value_type)0)?set_width("~e_R-",&M_seR):set_widthless("~e_R-");
	    (W_smuL>(value_type)0)?set_width("~mu_L-",&M_smuL):set_widthless("~mu_L-");
	    (W_smuR>(value_type)0)?set_width("~mu_R-",&M_smuR):set_widthless("~mu_R-");
	    (W_stauL>(value_type)0)?set_width("~tau_L-",&M_stauL):set_widthless("~tau_L-");
	    (W_stauR>(value_type)0)?set_width("~tau_R-",&M_stauR):set_widthless("~tau_R-");
	}
    }

    /* Function returning the QCD colours: */

    std::size_t susy_QED::QCD_colours()
    {
	return 1;
    }

    /* EM fine structure input: */

    void susy_QED::set_alpha(const value_type& a)
    {
	alpha=a;
	refresh_couplings();
    }

    /* Function setting the strong coupling: */

    void susy_QED::set_alpha_s(const value_type& a){}

    /* Function setting the QCD scale: */

    void susy_QED::set_QCD_scale(const value_type& mu){}

    /* Function setting the photon propagator to the unitary gauge: */

    void susy_QED::set_unitary_gauge()
    {
	if(initialised() and gauge!=1)
	{
	    set_propagator<unitary_gauge>("gamma");
	}
	gauge=1;
    }

    /* Function setting the photon propagator to the Feynman gauge: */

    void susy_QED::set_Feynman_gauge()
    {
	if(initialised() and gauge!=0)
	{
	    set_propagator<Feynman_gauge>("gamma");
	}
	gauge=0;
    }

    /* Function setting the photon propagator to the R-xi gauge: */

    void susy_QED::set_R_xi_gauge()
    {
	if(initialised() and gauge!=2)
	{
	    set_propagator<R_vector_gauge>("gamma");
	}
	gauge=2;
    }

    /* Function setting the gauge parameter: */

    void susy_QED::set_xi(const value_type& x)
    {
	R_gauge<susy_QED>::xi=x;
    }
}


