//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/scalar_QED.h>
#include <Camgen/vector_particle.h>
#include <Camgen/scalar_particle.h>
#include <Camgen/vss.h>
#include <Camgen/ssvv.h>
#include <Camgen/SM_params.h>

namespace Camgen
{
    /* Static constant integer member definitions: */

    const std::size_t scalar_QED::dimension;
    const bool scalar_QED::coloured;
    const bool scalar_QED::continuous_helicities;
    const int scalar_QED::beam_direction;
    
    /* Initialisation of static floating-point data members: */
    
    scalar_QED::value_type scalar_QED::alpha=(scalar_QED::value_type)SM_params::alpha;
    const scalar_QED::value_type scalar_QED::pi=std::acos(-(scalar_QED::value_type)1);
    scalar_QED::value_type scalar_QED::Q_e=std::sqrt((scalar_QED::value_type)4*scalar_QED::pi*scalar_QED::alpha);
    std::complex<scalar_QED::value_type>scalar_QED::gammaee(0,scalar_QED::Q_e);
    std::complex<scalar_QED::value_type>scalar_QED::eegammagamma(0,(scalar_QED::value_type)2*scalar_QED::Q_e*scalar_QED::Q_e);
    
    /* Gauge information initialisation (initialised as Feynman gauge): */
    
    int scalar_QED::gauge=0;

    /* Constructor: */

    scalar_QED::scalar_QED()
    {
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
	add_scalars("e-","e+");
	
	add_vertex<vss>("gamma","e+","e-",&gammaee);
	add_vertex<ssvv>("e+","e-","gamma","gamma",&eegammagamma);
    }

    /* Recompute couplings from alpha: */

    void scalar_QED::refresh_couplings()
    {
	Q_e=std::sqrt((value_type)4*pi*alpha);
	gammaee=std::complex<value_type>(0,Q_e);
	eegammagamma=std::complex<value_type>(0,(value_type)2*Q_e*Q_e);
    }

    /* Input value of fine-structure constant: */

    void scalar_QED::set_alpha(const value_type& a)
    {
	alpha=a;
	refresh_couplings();
    }

    /* Switch to unitary gauge: */

    void scalar_QED::set_unitary_gauge()
    {
	if(initialised() and gauge!=1)
	{
	    set_propagator<unitary_gauge>("gamma");
	}
	gauge=1;
    }

    /* Switch to Feynman gauge: */

    void scalar_QED::set_Feynman_gauge()
    {
	if(initialised() and gauge!=0)
	{
	    set_propagator<Feynman_gauge>("gamma");
	}
	gauge=0;
    }

    /* Switch to R_xi-gauge: */

    void scalar_QED::set_R_xi_gauge()
    {
	if(initialised() and gauge!=2)
	{
	    set_propagator<R_vector_gauge>("gamma");
	}
	gauge=2;
    }

    /* Input value for xi: */

    void scalar_QED::set_xi(const value_type& x)
    {
	R_gauge<scalar_QED>::xi=x;
    }
}


