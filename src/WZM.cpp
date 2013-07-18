//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/WZM.h>
#include <Camgen/sss.h>
#include <Camgen/ssss.h>
#include <Camgen/sffL.h>
#include <Camgen/sffR.h>
#include <Camgen/fermion.h>
#include <Camgen/scalar_particle.h>

#define DEFAULT_MASS 120.00
#define DEFAULT_COUPLING 0.25

namespace Camgen
{
    /* Compile-time constants: */

    const std::size_t WZM::dimension;
    const bool WZM::coloured;
    const bool WZM::continuous_helicities;
    const int WZM::beam_direction;

    /* Initialisation of couplings: */

    WZM::value_type WZM::m=(WZM::value_type)DEFAULT_MASS;
    WZM::value_type WZM::g=(WZM::value_type)DEFAULT_COUPLING;
    std::complex<WZM::value_type>WZM::phi3(0,-WZM::m*WZM::g);
    std::complex<WZM::value_type>WZM::phi4(0,-WZM::g*WZM::g);
    std::complex<WZM::value_type>WZM::phipsipsi(0,-WZM::g);

    /* Constructor: */

    WZM::WZM()
    {
	add_scalars("phiL","phiR",&m);
	add_fermion("psi",&m);
	add_vertex<sss>("phiL","phiR","phiR",&phi3);
	add_vertex<sss>("phiR","phiL","phiL",&phi3);
	add_vertex<ssss>("phiL","phiL","phiR","phiR",&phi4);
	add_vertex<sffL>("phiL","psi","psi",&phipsipsi);
	add_vertex<sffR>("phiR","psi","psi",&phipsipsi);
    }

    /* Mass input: */

    void WZM::set_mass(const value_type& mass)
    {
	if(m==(value_type)0)
	{
	    if(mass>m)
	    {
		m=mass;
		phi3=std::complex<value_type>(0,-m*g);
		add_vertex<sss>("phiL","phiR","phiR",&phi3);
		add_vertex<sss>("phiR","phiL","phiL",&phi3);
	    }
	}
	else
	{
	    if(mass==(value_type)0)
	    {
		m=mass;
		phi3=0;
		erase_vertex("phiL","phiR","phiR");
		erase_vertex("phiR","phiL","phiL");
	    }
	    else
	    {
		m=mass;
		phi3.imag()=-m*g;
	    }
	}
    }
	    
    /* Coupling input: */
	    
    void WZM::set_coupling(const value_type& c)
    {
	g=c;
	refresh_couplings();
    }

    /* Function computing the couplings from g and m: */

    void WZM::refresh_couplings()
    {
	std::complex<WZM::value_type>phi3(0,-m*g);
	std::complex<WZM::value_type>phi4(0,-g*g);
	std::complex<WZM::value_type>phipsipsi(0,-g);
    }
}

#undef DEFAULT_MASS
#undef DEFAULT_COUPLING

