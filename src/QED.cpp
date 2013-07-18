//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/QED.h>
#include <Camgen/vector_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/vff.h>
#include <Camgen/SM_params.h>

namespace Camgen
{
    /* Compile time constants definitions: */

    const std::size_t QED::dimension;
    const bool QED::coloured;
    const bool QED::continuous_helicities;
    const int QED::beam_direction;
    
    /* Default alpha-value definition: */
    
    QED::value_type QED::alpha=(QED::value_type)SM_params::alpha;

    /* Alpha strong (not used): */

    const QED::value_type QED::alpha_s=-(QED::value_type)1;

    /* QCD scale (not used): */

    const QED::value_type QED::QCD_scale=-(QED::value_type)1;

    /* Pi computation: */

    const QED::value_type QED::pi=std::acos(-(QED::value_type)1);
    
    /* Electron charge initialisation: */
    
    QED::value_type QED::Q_e=std::sqrt((QED::value_type)4*QED::pi*QED::alpha);

    /* QED coupling initialisation: */

    std::complex<QED::value_type>QED::gammaee((QED::value_type)0,QED::Q_e);
    std::complex<QED::value_type>QED::gammamumu=QED::gammaee;
    std::complex<QED::value_type>QED::gammatautau=QED::gammaee;

    /* Lepton mass initialisations: */

    QED::value_type QED::M_e=(QED::value_type)0;
    QED::value_type QED::M_mu=(QED::value_type)0;
    QED::value_type QED::M_tau=(QED::value_type)SM_params::M_tau;

    /* Gauge parameter: */

    int QED::gauge=0;

    /* Constructor: */

    QED::QED()
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
	
	/* Electromagnetic interactions: */

	add_vertex<vff>("gamma","e+","e-",&gammaee);
	add_vertex<vff>("gamma","mu+","mu-",&gammamumu);
	add_vertex<vff>("gamma","tau+","tau-",&gammatautau);

	/* Positively-charged lepton family definition: */

	construct_family("l+","e+,mu+");

	/* Positively-charged lepton family definition, including taus: */

	construct_family("L+","e+,mu+,tau+");

	/* Negatively-charged lepton family definition: */

	construct_family("l-","e-,mu-");

	/* Negatively-charged lepton family definition, including taus: */

	construct_family("L-","e-,mu-,tau-");
    }

    /* Function computing the electron charge and vertex coupling from alpha: */

    void QED::refresh_couplings()
    {
	Q_e=std::sqrt((value_type)4*pi*alpha);
	gammaee=std::complex<value_type>(0,Q_e);
	gammamumu=gammaee;
	gammatautau=gammaee;
    }

    /* Refresh particle masses (to use when along the run massless fermions are
     * to be made massive): */

    void QED::refresh_fermion_masses()
    {
	if(initialised())
	{
	    (M_e>(value_type)0)?set_mass("e-",&M_e):set_massless("e-");
	    (M_mu>(value_type)0)?set_mass("mu-",&M_mu):set_massless("mu-");
	    (M_tau>(value_type)0)?set_mass("tau-",&M_tau):set_massless("tau-");
	}
    }

    /* Function returning the QCD colours: */

    std::size_t QED::QCD_colours()
    {
	return 1;
    }

    /* Input of alpha and computation of the couplings: */

    void QED::set_alpha(const QED::value_type& a)
    {
	alpha=a;
	refresh_couplings();
    }

    /* Dummy function setting the strong coupling: */

    void QED::set_alpha_s(const QED::value_type& a){}

    /* Dummy function setting the QCD scale: */

    void QED::set_QCD_scale(const QED::value_type& mu){}

    /* Function setting the photon propagator to the Feynman gauge: */

    void QED::set_unitary_gauge()
    {
	if(initialised() and gauge!=1)
	{
	    set_propagator<unitary_gauge>("gamma");
	}
	gauge=1;
    }

    /* Function setting the photon propagator to the unitary gauge: */

    void QED::set_Feynman_gauge()
    {
	if(initialised() and gauge!=0)
	{
	    set_propagator<Feynman_gauge>("gamma");
	}
	gauge=0;
    }

    /* Function setting the photon propagator to the R-xi gauge: */

    void QED::set_R_xi_gauge()
    {
	if(initialised() and gauge!=2)
	{
	    set_propagator<R_vector_gauge>("gamma");
	}
	gauge=2;
    }

    /* Function setting the xi-parameter: */

    void QED::set_xi(const value_type& x)
    {
	R_gauge<QED>::xi=x;
    }

    /* Reset all parameters using the default input parameters in the
     * SM_params.h file: */

    void QED::set_default_params()
    {
	M_e=(value_type)0;
	M_mu=(value_type)0;
	M_tau=(value_type)SM_params::M_tau;
	refresh_fermion_masses();
	alpha=(value_type)SM_params::alpha;
	refresh_couplings();
	set_Feynman_gauge();
    }
}

