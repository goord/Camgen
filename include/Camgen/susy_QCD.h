//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SUSY_QCD_H_
#define CAMGEN_SUSY_QCD_H_

#include <Camgen/Minkowski.h>
#include <Camgen/col_flow.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/hel_type.h>
#include <Camgen/model.h>
#include <Camgen/SM_params.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the supersymmetric QCD model class in Camgen. The coupling    *
 * constants may be inserted manually, or computed from the basic input variable *
 * alpha_s or QCD_scale.                                                         *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class susy_QCD: public model<susy_QCD>
    {
	public:

	    /* Compile-time data: */

	    typedef double value_type;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    typedef colour_flow colour_treatment;
	    typedef helicity_type spin_vector_type;
	    
	    static const std::size_t dimension=4;
	    static const bool coloured=true;
	    static const bool continuous_colours=false;
	    static const bool continuous_helicities=true;
	    static const int beam_direction=3;
	    static const std::size_t N_c=SM_params::N_c;

	    /* Fine structure constant (not used): */

	    static const value_type alpha;

	    /*  Strong coupling constant: */

	    static value_type alpha_s;

	    /* QCD scale: */

	    static value_type QCD_scale;

	    /* QCD coupling: */

	    static value_type g_s;
	    
	    /* Pi: */
	    
	    static const value_type pi;
	    
	    /* Particle masses: */

	    static value_type M_u,M_d,M_c,M_s,M_t,M_b,M_suL,M_suR,M_sdL,M_sdR,M_scL,M_scR,M_ssL,M_ssR,M_stL,M_stR,M_sbL,M_sbR,M_sg;

	    /* Particle widths: */

	    static value_type W_t,W_suL,W_suR,W_sdL,W_sdR,W_scL,W_scR,W_ssL,W_ssR,W_stL,W_stR,W_sbL,W_sbR,W_sg;

	    /* Triple gluon coupling: */

	    static std::complex<value_type> ggg;
	    
	    /* 4-gluon coupling: */

	    static std::complex<value_type> gggg;
	    
	    /* Auxiliary tensor-gluon coupling: */

	    static std::complex<value_type> Tgg;

	    /* Gluon-quark-quark couplings: */

	    static std::complex<value_type> guu,gdd,gcc,gss,gtt,gbb;

	    /* Gluon-squark-squark couplings: */

	    static std::complex<value_type> gsuLsuL,gsuRsuR,gsdLsdL,gsdRsdR,gscLscL,gscRscR,gssLssL,gssRssR,gstLstL,gstRstR,gsbLsbL,gsbRsbR;

	    /* Squark-gluino-quark couplings: */
	    
	    static std::complex<value_type> suLsgu,suRsgu,sdLsgd,sdRsgd,scLsgc,scRsgc,ssLsgs,ssRsgs,stLsgt,stRsgt,sbLsgb,sbRsgb;

	    /* Antisquark-antiquark-gluino couplings: */

	    static std::complex<value_type> suLusg,suRusg,sdLdsg,sdRdsg,scLcsg,scRcsg,ssLssg,ssRssg,stLtsg,stRtsg,sbLbsg,sbRbsg;
	    
	    /* Squark-squark-gluon-gluon couplings: */
	    
	    static std::complex<value_type>suLsuLgg,suRsuRgg,sdLsdLgg,sdRsdRgg,scLscLgg,scRscRgg,ssLssLgg,ssRssRgg,stLstLgg,stRstRgg,sbLsbLgg,sbRsbRgg;

	    /* 4-squark couplings: */

	    static std::complex<value_type>suLsuLsuLsuL,suRsuRsuRsuR,sdLsdLsdLsdL,sdRsdRsdRsdR,scLscLscLscL,scRscRscRscR,ssLssLssLssL,ssRssRssRssR,stLstLstLstL,stRstRstRstR,sbLsbLsbLsbL,sbRsbRsbRsbR;
	    
	    /* 2-2-squark couplings: */
	    
	    static std::complex<value_type>suLsuLsuRsuR,sdLsdLsdRsdR,scLscLscRscR,ssLssLssRssR,stLstLstRstR,sbLsbLsbRsbR;

	    /* Constructor: */

	    susy_QCD();
	    
	    /* Computation of the couplings: */

	    static void refresh_couplings();
	    
	    /* Function refreshing the particle masses: */

	    static void refresh_masses();

	    /* Function refreshing the particle decay widths: */

	    static void refresh_widths();

	    /* Function returning the number of QCD colours: */
	    
	    static std::size_t QCD_colours();
	    
	    /* Input of fine structure constant (not used): */
	    
	    static void set_alpha(const value_type&);

	    /* Input of alpha-strong and computation of the couplings: */

	    static void set_alpha_s(const value_type&);
	    
	    /* Input of the QCD scale and (LO) computation of alpha-strong and
	     * the couplings: */
	    
	    static void set_QCD_scale(const value_type&);
	    
	    /* Function setting the gluon propagator to the Feynman gauge: */
	    
	    static void set_Feynman_gauge();

	    /* Function setting the gluon propagator to the unitary gauge: */

	    static void set_unitary_gauge();
	    
	    /* Function setting the gluon propagator to the R-gauge: */
	    
	    static void set_R_xi_gauge();

	    /* Input of the gauge parameter: */

	    static void set_xi(const value_type&);

	    /* Function switching to a 4-gluon vertex Lagrangian: */

	    static void set_4_gluon_vertex();

	    /* Function switching to an auxiliary tensor field description of
	     * the 4-gluon vertex: */

	    static void set_auxiliary_QCD_field();
	    
	private:

	    /* Gauge flag: */

	    static int gauge;

	    /* 4-gluon vertex flag (true: model constructs a 4-gluon vertex,
	     * false: model uses the antisymmetric tensor field): */

	    static bool four_gluon_vertex;
    };
}

#endif /*CAMGEN_SUSY_QCD_H_*/

