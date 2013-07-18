//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SUSY_QED_H_
#define CAMGEN_SUSY_QED_H_

#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/hel_type.h>
#include <Camgen/model.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the supersymmetric QED model class in Camgen. The coupling    *
 * constants may be inserted manually,or computed from the basic input variable  *
 * alpha. Unless manually modified, the sleptons all have a common mass defined  *
 * in the susy_params header.                                                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class susy_QED: public model<susy_QED>
    {
	public:

	    /* Compile-time data: */

	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    static const std::size_t dimension=4;
	    static const bool coloured=false;
	    typedef Pauli_basis Dirac_algebra_type;
	    static const bool continuous_helicities=true;
	    typedef helicity_type spin_vector_type;
	    static const int beam_direction=3;

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

	    susy_QED();

	    /* Function computing the electron charge and vertex coupling from alpha: */

	    static void refresh_couplings();

	    /* Function refreshing the particle masses: */

	    static void refresh_masses();

	    /* Function refreshing the particle decay widths: */

	    static void refresh_widths();
	    
	    /* Output of the number of QCD colours (returns one): */
	    
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

	private:

	    /* Gauge tag: */

	    static int gauge;
    };
}

#endif /*CAMGEN_SUSY_QED_H_*/

