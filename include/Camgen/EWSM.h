//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_EWSM_H_
#define CAMGEN_EWSM_H_

#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/hel_type.h>
#include <Camgen/model.h>
#include <Camgen/SM_params.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the electroweak colourless standard model class in Camgen.  *
 * The coupling constant values may be inserted manually, or computed by the   *
 * class from the basic input variables alpha, Fermi constant G_F, the Z-boson *
 * mass and the top and bottom masses. Gauge-invariant inclusion of particle   *
 * decay widths is implemented via the complex mass scheme, replacing all mass *
 * parameters in the computed couplings by the (complex) square root of        *
 *                                                                             *
 * 				M^2-iM\Gamma                                   *
 *                                                                             *
 * The decay width can again be manually inserted, or computed by the class.   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class EWSM: public model<EWSM>
    {
	public:

	    /* The usual compile-time model data: */

	    static const std::size_t dimension=4;
	    static const bool coloured=false;
	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    static const bool continuous_helicities=true;
	    typedef helicity_type spin_vector_type;
	    static const int beam_direction=3;
	    static const std::size_t N_c=SM_params::N_c;
	    static const value_type pi;

	    /* Particle masses: */

	    static value_type M_e,M_mu,M_tau,M_u,M_d,M_c,M_s,M_t,M_b,M_Z,M_W,M_h0;

	    /* Particle decay widths: */

	    static value_type W_t,W_Z,W_W,W_h0;
	    
	    /* Fine structure constant: */

	    static value_type alpha;

	    /* Alpha strong: */

	    static const value_type alpha_s;
	    
	    /* QCD scale: */
	    
	    static const value_type QCD_scale;

	    /* Fermi's constant: */

	    static value_type G_F;

	    /* CKM sines of Euler angles: */

	    static value_type s12,s23,s13;

	    /* CKM CP-violating phase: */

	    static value_type delta;
	    
	    /* Quadruple-vector couplings: */

	    static std::complex<value_type> WWWW,WZWZ,WgammaWZ,WgammaWgamma,gggg;

	    /* Triple-vector couplings: */

	    static std::complex<value_type> gammaWW,ZWW;

	    /* Quadruple-scalar couplings: */

	    static std::complex<value_type> hhhh,hhchichi,hhphiphi,chichichichi,chichiphiphi,phiphiphiphi;

	    /* Triple-scalar couplings: */

	    static std::complex<value_type> hhh,hchichi,hphiphi;

	    /* Vector-fermion-fermion couplings: */

	    static std::complex<value_type> gammaee,gammamumu,gammatautau,gammauu,gammacc,gammatt,gammadd,gammass,gammabb;

	    static std::complex<value_type> Znene,Znmunmu,Zntauntau,Zee_V,Zee_A,Zmumu_V,Zmumu_A,Ztautau_V,Ztautau_A,Zuu_V,Zuu_A,Zdd_V,Zdd_A,Zcc_V,Zcc_A,Zss_V,Zss_A,Ztt_V,Ztt_A,Zbb_V,Zbb_A;

	    static std::complex<value_type> Wnee,Wene,Wnmumu,Wmunmu,Wntautau,Wtauntau,Wud,Wdu,Wus,Wsu,Wub,Wbu,Wcd,Wdc,Wcs,Wsc,Wcb,Wbc,Wtd,Wdt,Wts,Wst,Wtb,Wbt;

	    /* Scalar-fermion-fermion couplings: */

	    static std::complex<value_type> hee,hmumu,htautau,huu,hdd,hcc,hss,htt,hbb;

	    static std::complex<value_type> chiee,chimumu,chitautau,chiuu,chidd,chicc,chiss,chitt,chibb;

	    static std::complex<value_type> phinee,phiene,phinmumu,phimunmu,phintautau,phitauntau,phiud_S,phiud_A,phidu_S,phidu_A,phius_S,phius_A,phisu_S,phisu_A,phiub_S,phiub_A,phibu_S,phibu_A,phicd_S,phicd_A,phidc_S,phidc_A,phics_S,phics_A,phisc_S,phisc_A,phicb_S,phicb_A,phibc_S,phibc_A,phitd_S,phitd_A,phidt_S,phidt_A,phits_S,phits_A,phist_S,phist_A,phitb_S,phitb_A,phibt_S,phibt_A;

	    /* Vector-scalar-scalar couplings: */

	    static std::complex<value_type> gammachih,Zchih,gammaphiphi,Zphiphi,Wpphimh,Wmphiph,Wpphimchi,Wmphipchi;

	    /* Scalar-vector-vector couplings: */

	    static std::complex<value_type> hWW,hZZ,hZgamma,phipWmZ,phimWpZ,phipWmgamma,phimWpgamma;
	    
	    /* Scalar-scalar-vector-vector couplings: */

	    static std::complex<value_type> hhWW,chichiWW,phiphiWW,phiphiZZ,phiphigammaZ,phiphigammagamma,hhZZ,chichiZZ,phiphWmZ,phimhWpZ,phiphWmgamma,phimhWpgamma,phimchiWpZ,phipchiWmZ,phimchiWpgamma,phipchiWmgamma;

	    /* Constructor declaration: */

	    EWSM();
	    
	    /* Output of the number of QCD colours: */

	    static std::size_t QCD_colours();

	    /* Function computing the couplings from given parameters: */

	    static void refresh_couplings();

	    /* Function computing the couplings from given parameters: */

	    static void refresh_strong_couplings();

	    /* Function computing the couplings from given parameters: */

	    static void refresh_weak_couplings();

	    /* Function computing the couplings from given parameters: */

	    static void refresh_Yukawa_couplings();

	    /* Function computing the CKM matrix from given parameters: */

	    static void refresh_CKM_matrix();

	    /* Function putting possible nonzero fermion masses in the model
	     * during runtime: */

	    static void refresh_fermion_masses();

	    /* Function computing the LO standard model W,Z,h and t widths: */

	    static void refresh_widths();
	    
	    /* Function computing the LO standard model W width: */
	    
	    static void refresh_W_width();

	    /* Function computing the LO standard model Z width: */

	    static void refresh_Z_width();
	    
	    /* Function computing the LO standard model H width: */
	    
	    static void refresh_Higgs_width();

	    /* Function computing the LO standard model t width: */

	    static void refresh_top_width();
	    
	    /* Function setting the fine structure constant to a given value and
	     * computing all couplings depending on alpha: */
	    
	    static void set_alpha(const value_type&);
	    
	    /* Function setting alpha strong and computing the couplings
	     * depending on it: */
	    
	    static void set_alpha_s(const value_type&);
	    
	    /* Function setting the strong scale and computing alpha_s (at one
	     * loop) and the resulting strong-interaction vertex couplings: */
	    
	    static void set_QCD_scale(const value_type&);

	    /* Function setting the Higgs mass and computing the couplings
	     * depending on it: */

	    static void set_Higgs_mass(const value_type&);

	    /* Function setting all widths to zero: */

	    static void discard_widths();

	    /* Function setting all vector propagators to the Feynman gauge: */

	    static void set_Feynman_gauge();
	    
	    /* Function setting all vector propagators to the unitary gauge and
	     * deleting all would-be Goldstones: */
	    
	    static void set_unitary_gauge();

	    /* Function setting all vector propagators to the R-xi gauge: */

	    static void set_R_xi_gauge();
	    
	    /* Input of the xi-value in the R_xi gauge: */
	    
	    static void set_xi(const value_type&);
	    
	    /* Reset all physical parameters to their default values (defined in
	     * SM_params): */

	    static void set_default_params();

	private:
	    
	    /* Unit charge: */

	    static value_type Q_e;
	    
	    /* Complex gauge boson masses: */

	    static std::complex<value_type> cM_Z,cM_W;
	    
	    /* Complex cosine and sine of the Weinberg angle: */
	    
	    static std::complex<value_type> cos_W,sin_W;

	    /* Temporary internal parameters: */

	    static std::complex<value_type> tan_W,cos2_W,sin2_W;

	    /* CKM cosines: */

	    static value_type c12,c23,c13;

	    /* CKM matrix: */

	    static std::complex<value_type>V_CKM[3][3];

	    /* Gauge flag: */

	    static int gauge;
    };
}

#endif /*CAMGEN_EWSM_H_*/


