//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_QED_BASE_H_
#define CAMGEN_QED_BASE_H_

#include <Camgen/model.h>
#include <Camgen/vector_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/vff.h>
#include <Camgen/SM_params.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the QED_base class template in Camgen.       *
 * Publicly deriving a model class model_t from QED_base<model_t> amounts to   *
 * adding QED to the model. It should be noted that the value type of base_QED *
 * has to be castable to that of the derived model.                            *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */                                                                            

namespace Camgen
{
    template<class model_t,class value_t>class QED_base: public model<model_t>
    {
	public:

	    /* Value type definition (is always casted to model_t's value type):
	     * */

	    typedef value_t value_type;

	    /* Fine structure constant: */

	    static value_type alpha;

	    /*  Strong coupling constant (only relevant for LH output): */

	    static const value_type alpha_s;

	    /* QCD scale (only relevant for LH output): */

	    static const value_type QCD_scale;
	    
	    /* Eelectron charge: */
	    
	    static value_type Q_e;

	    /* Pi: */

	    static const value_type pi;
	    
	    /* QED vertex coupling constants: */
	    
	    static std::complex<value_type>gammaee,gammamumu,gammatautau;

	    /* Lepton masses: */

	    static value_type M_e,M_mu,M_tau;

	    /* Constructor: */
	    
	    QED_base()
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
	
		/* Electromagnetic interactions: */

		model<model_t>::template add_vertex<vff>("gamma","e+","e-",&gammaee);
		model<model_t>::template add_vertex<vff>("gamma","mu+","mu-",&gammamumu);
		model<model_t>::template add_vertex<vff>("gamma","tau+","tau-",&gammatautau);
		
		/* Positively-charged lepton family definition: */

		model<model_t>::construct_family("l+","e+,mu+");

		/* Positively-charged lepton family definition, including taus: */

		model<model_t>::construct_family("L+","e+,mu+,tau+");

		/* Negatively-charged lepton family definition: */

		model<model_t>::construct_family("l-","e-,mu-");

		/* Negatively-charged lepton family definition, including taus: */

		model<model_t>::construct_family("L-","e-,mu-,tau-");
	    }

	    /* Function computing the electron charge and vertex coupling from alpha: */

	    static void refresh_couplings()
	    {
		Q_e=std::sqrt((value_type)4*pi*alpha);
		gammaee=std::complex<value_type>(0,Q_e);
		gammamumu=gammaee;
		gammatautau=gammaee;
	    }

	    /* Refresh particle masses (to use when along the run massless fermions are
	     * to be made massive): */

	    void refresh_fermion_masses()
	    {
		if(model<model_t>::initialised())
		{
		    (M_e>(value_type)0)?model<model_t>::set_mass("e-",&M_e):model<model_t>::set_massless("e-");
		    (M_mu>(value_type)0)?model<model_t>::set_mass("mu-",&M_mu):model<model_t>::set_massless("mu-");
		    (M_tau>(value_type)0)?model<model_t>::set_mass("tau-",&M_tau):model<model_t>::set_massless("tau-");
		}
	    }
	    
	    /* Output of the number of QCD colours (returns zero): */
	    
	    static std::size_t QCD_colours()
	    {
		return 0;
	    }
	    
	    /* Input of alpha and computation of the couplings: */
	    
	    static void set_alpha(const value_type& a)
	    {
		alpha=a;
		refresh_couplings();
	    }

	    /* Input of the strong coupling (not used): */

	    static void set_alpha_s(const value_type& a){}

	    /* Input of the QCD scale (not used): */

	    static void set_QCD_scale(const value_type& a){}

	    /* Function setting the photon propagator to the Feynman gauge: */

	    static void set_Feynman_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=0)
		{
		    model<model_t>::template set_propagator<Feynman_gauge>("gamma");
		}
		gauge=1;
	    }	

	    /* Function setting the photon propagator to the unitary gauge: */

	    static void set_unitary_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=1)
		{
		    model<model_t>::template set_propagator<unitary_gauge>("gamma");
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

	    /* Function setting the xi-parameter: */

	    static void set_xi(const value_type& x)
	    {
		R_gauge<model_t>::xi=x;
	    }

	private:

	    /* Gauge tag: */

	    static int gauge;
    };
    
    /* Default alpha-value definition: */
    
    template<class model_t,class value_t>typename QED_base<model_t,value_t>::value_type QED_base<model_t,value_t>::alpha=(value_t)SM_params::alpha;

    /* Alpha strong (not used): */

    template<class model_t,class value_t>const typename QED_base<model_t,value_t>::value_type QED_base<model_t,value_t>::alpha_s=-(value_t)1;

    /* QCD scale (not used): */

    template<class model_t,class value_t>const typename QED_base<model_t,value_t>::value_type QED_base<model_t,value_t>::QCD_scale=-(value_t)1;

    /* Pi computation: */

    template<class model_t,class value_t>const typename QED_base<model_t,value_t>::value_type QED_base<model_t,value_t>::pi=std::acos(-(value_t)1);
    
    /* Electron charge and coupling: */
    
    template<class model_t,class value_t>typename QED_base<model_t,value_t>::value_type QED_base<model_t,value_t>::Q_e=std::sqrt((value_t)4*QED_base<model_t,value_t>::pi*QED_base<model_t,value_t>::alpha);
    template<class model_t,class value_t>std::complex<value_t>QED_base<model_t,value_t>::gammaee(0,QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>QED_base<model_t,value_t>::gammamumu(0,QED_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>QED_base<model_t,value_t>::gammatautau(0,QED_base<model_t,value_t>::Q_e);

    /* Lepton mass initialisation: */

    template<class model_t,class value_t>value_t QED_base<model_t,value_t>::M_e=(value_t)0;
    template<class model_t,class value_t>value_t QED_base<model_t,value_t>::M_mu=(value_t)0;
    template<class model_t,class value_t>value_t QED_base<model_t,value_t>::M_tau=(value_t)SM_params::M_tau;
    
    /* Gauge parameter initialisation: */

    template<class model_t,class value_t>int QED_base<model_t,value_t>::gauge=0;
}

#endif /*CAMGEN_QED_BASE_H_*/

