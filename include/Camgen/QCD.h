//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_QCD_H_
#define CAMGEN_QCD_H_

#include <Camgen/Minkowski.h>
#include <Camgen/col_flow.h>
#include <Camgen/adjoint.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/hel_type.h>
#include <Camgen/model.h>
#include <Camgen/SM_params.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the QCD model class  in Camgen. The couplings may be inserted *
 * manually, or be computed from the alpha_s input value, or from the QCD scale  *
 * input by a leading-order evolution of alpha_s. Furthermore, the gluon gauge   *
 * may be chosen at runtime, and the gluon 4-vertex may be replaced by an        *
 * antisymmetric auxiliary tensor field.                                         *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class QCD: public model<QCD>
    {
	public:

	    /* The usual compile-time data: */

	    typedef double value_type;
	    static const std::size_t dimension=4;
	    typedef Minkowski_type spacetime_type;
	    typedef Pauli_basis Dirac_algebra_type;
	    typedef colour_flow colour_treatment;
	    static const bool continuous_helicities=true;
	    static const bool coloured=true;
	    static const bool continuous_colours=false;
	    typedef helicity_type spin_vector_type;
	    static const int beam_direction=3;
	    const static std::size_t N_c=SM_params::N_c;
	    
	    /* EM fine structure constant (unused): */

	    static const value_type alpha;

	    /* alpha strong: */

	    static value_type alpha_s;
	    
	    /* QCD scale parameter: */
	    
	    static value_type QCD_scale;

	    /* QCD coupling: */

	    static value_type g_s;
	    
	    /* Pi: */
	    
	    static const value_type pi;

	    /* Quark masses: */

	    static value_type M_u,M_d,M_c,M_s,M_t,M_b;

	    /* Triple-gluon coupling: */

	    static std::complex<value_type>ggg;
	    
	    /* Quadruple-gluon coupling: */
	    
	    static std::complex<value_type>gggg;

	    /* Tensor-gluon coupling: */

	    static std::complex<value_type>Tgg;
	    
	    /* Gluon-quark couplings: */
	    
	    static std::complex<value_type> guu,gdd,gcc,gss,gtt,gbb;

	    /* Constructor: */

	    QCD();
	    
	    /* Output of the number of QCD colours: */

	    static std::size_t QCD_colours();

	    /* Dummy function setting the EM fine-structure constant: */

	    static void set_alpha(const value_type&);
	    
	    /* Function setting alpha strong and computing the couplings: */
	    
	    static void set_alpha_s(const value_type&);

	    /* Function setting the scale and computing alpha_s (at LO) and the
	     * couplings: */

	    static void set_QCD_scale(const value_type&);
	    
	    /* Function computing the couplings from the current alpha_s value:
	     * */
	    
	    static void refresh_couplings();

	    /* Function inserting nonzero or zero fermion masses during runtime:
	     * */

	    static void refresh_fermion_masses();

	    /* Function setting the gluon-propagator to the Feynman gauge: */

	    static void set_Feynman_gauge();
	    
	    /* Function setting the gluon-propagator to the unitary gauge: */
	    
	    static void set_unitary_gauge();

	    /* Function setting the gluon-propagator to the R-xi gauge: */

	    static void set_R_xi_gauge();
	    
	    /* Function setting the xi-parameter: */
	    
	    static void set_xi(const value_type&);

	    /* Function switching to a 4-gluon vertex Lagrangian: */

	    static void set_4_gluon_vertex();

	    /* Function switching to an auxiliary tensor field description of
	     * the 4-gluon vertex: */

	    static void set_auxiliary_QCD_field();
	    
	    /* Reset all physical parameters to their default values: */

	    static void set_default_params();
	
	private:

	    /* 4-gluon vertex tag (true: model constructs a 4-gluon vertex,
	     * false: model uses the antisymmetric tensor field): */

	    static bool four_gluon_vertex;

	    /* Gauge tag: */

	    static int gauge;
    };
}

#endif /*CAMGEN_QCD_H_*/

