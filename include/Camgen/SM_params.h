//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SM_PARAMS_H_
#define CAMGEN_SM_PARAMS_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Default standard-model masses and couplings. These numbers were taken from  *
 * the particle data group in 2013.                                            *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace SM_params
    {
	/* Lepton masses: */

	const double M_e=0.000511;
	const double M_mu=0.1057;
	const double M_tau=1.777;
	
	/* Quark masses: */
	
	const double M_u=0.00255;
	const double M_d=0.00504;
	const double M_c=1.25;
	const double M_s=0.12;
	const double M_t=174;
	const double M_b=4.20;

	/* Gauge boson masses: */

	const double M_Z=91.1876;
	const double M_W=80.419;

	/* Fine-structure constant: */

	const double alpha=(double)1/127.9;

	/* Fermi's constant: */

	const double G_F=1.16639e-05;

	/* Strong coupling constant, evaluated at scale: */

	const double QCD_scale=M_Z;
	const double alpha_s=0.1172;

	/* Cabibbo angle: */

	const double s12=0.226;

	/* ub-mixing angle: */

	const double s23=0.0415;

	/* cb-mixing angle: */

	const double s13=0.0035;

	/* CP-violating phase: */

	const double delta=0.0209;

	/* Number of QCD colours: */

	const unsigned N_c=3;
    }
}

#endif /*CAMGEN_SM_PARAMS_H_*/

