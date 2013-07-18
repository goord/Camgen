//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_QED_H_
#define CAMGEN_QED_H_

#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/hel_type.h>
#include <Camgen/model.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the QED model class in Camgen. The coupling constant values *
 * may be inserted manually, or computed from the basic input variable alpha.  *
 * Moreover, photon gauge can be changed at runtime to either Feynman, unitary *
 * or R_xi gauge.                                                              *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class QED: public model<QED>
    {
	public:

	    /* The usual compile-time data: */

	    typedef double value_type;
	    static const std::size_t dimension=4;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    static const bool coloured=false;
	    static const bool continuous_helicities=false;
	    typedef helicity_type spin_vector_type;
	    static const int beam_direction=3;

	    /* Fine structure constant: */

	    static value_type alpha;

	    /*  Strong coupling constant (only relevant for LH output): */

	    static const value_type alpha_s;

	    /* QCD scale (only relevant for LH output): */

	    static const value_type QCD_scale;
	    
	    /* Eelectron charge: */
	    
	    static value_type Q_e;

	    /* Lepton masses: */

	    static value_type M_e,M_mu,M_tau;

	    /* Pi: */

	    static const value_type pi;
	    
	    /* QED vertices: */
	    
	    static std::complex<value_type>gammaee,gammamumu,gammatautau;

	    /* Constructor: */

	    QED();

	    /* Function computing the electron charge and vertex coupling from alpha: */

	    static void refresh_couplings();

	    /* Function inserting the nonzero fermion mass variables in the
	     * corresponding particle objects, if instantiated: */

	    static void refresh_fermion_masses();
	    
	    /* Output of the number of QCD colours (returns zero): */
	    
	    static std::size_t QCD_colours();
	    
	    /* Input of alpha and computation of the couplings: */
	    
	    static void set_alpha(const value_type&);

	    /* Input of the strong coupling (not used): */

	    static void set_alpha_s(const value_type&);

	    /* Input of the QCD scale (not used): */

	    static void set_QCD_scale(const value_type&);

	    /* Function setting the photon propagator to the Feynman gauge: */

	    static void set_Feynman_gauge();   

	    /* Function setting the photon propagator to the unitary gauge: */

	    static void set_unitary_gauge();

	    /* Function setting the photon propagator to the R-xi gauge: */

	    static void set_R_xi_gauge();

	    /* Function setting the xi-parameter: */

	    static void set_xi(const value_type&);

	    /* Function setting the parameters to the initial configuration: */

	    static void set_default_params();

	private:

	    /* Gauge tag: */

	    static int gauge;
    };
}

#endif /*CAMGEN_QED_H_*/

