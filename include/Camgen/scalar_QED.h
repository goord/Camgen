//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SCALAR_QED_H_
#define CAMGEN_SCALAR_QED_H_

#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/model.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Scalar-QED model declaration header. The particle content is a a charged      *
 * scalar pair and a photon field. The vertices are a ssv vertex with coupling   *
 * e=sqrt(4*pi*alpha) and a ssvv 4-vertex with coupling 2ie^2. The class allows  *
 * to change the photon gauge at runtime.                                        *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class scalar_QED: public model<scalar_QED>
    {
	public:

	    /* Numerical type definition: */

	    typedef double value_type;
	    
	    /* Spacetime dimensions definition: */
	    
	    static const std::size_t dimension=4;

	    /* Number of distinct colours a particle can have: */

	    static const bool coloured=false;
	    
	    /* Spacetime signature definition: */
	    
	    typedef Minkowski_type spacetime_type;

	    /* Dirac algebra basis definition: */

	    typedef Pauli_basis Dirac_algebra_type;

	    /* Beam direction: */

	    static const int beam_direction=3;
	    
	    /* Boolean denoting whether to mix helicities: */
	    
	    static const bool continuous_helicities=true;

	    /* Electromagnetic fine-structure constant: */

	    static value_type alpha;
	    
	    /* Scalar charge: */
	    
	    static value_type Q_e;

	    /* Useful constant: */

	    static const value_type pi;
	    
	    /* 3- and 4-vertex couplings: */
	    
	    static std::complex<value_type>gammaee;
	    static std::complex<value_type>eegammagamma;

	    /* Constructor: */

	    scalar_QED();
	    
	    /* Recompute couplings from value of alpha: */
	    
	    static void refresh_couplings();

	    /* Input alpha-value: */

	    static void set_alpha(const value_type&);
	    
	    /* Gauge-changing functions: */
	    
	    static void set_Feynman_gauge();
	    static void set_unitary_gauge();
	    static void set_R_xi_gauge();
	    static void set_xi(const value_type&);

	private:

	    /* Gauge information (0=Feynman gauge,1=unitary gauge,2=R_xi-gauge): */

	    static int gauge;
    };
}

#endif /*CAMGEN_SQED_H_*/

