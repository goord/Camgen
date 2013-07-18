//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_KS_TYPE_H_
#define CAMGEN_KS_TYPE_H_

#include <Camgen/vector.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the Kleiss-Stirling type spin vector class,     *
 * needed for the construction of massive spinors. This choice of spin vector is *
 * parallel to the four-vector k_0 defining the basis of basic spinors. The      *
 * standard choices, the Pauli and Weyl basis in Minkowski spacetime and beam    *
 * direction the z-axis, are explicitly specialised for optimisation purposes.   *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Wrapper class for the spin vector type: */

    class KS_type
    {
	public:
	    
	    /* Implementation class template, taking the spacetem type and beam
	     * direction as template arguments: */
	    
	    template<class spacetime_t,int beam_dir>class implementation
	    {
		public:
		    typedef vector<typename spacetime_t::value_type,spacetime_t::base_type::index_range> vector_type;

		    /* Construction of the spin vector. Simply calls a static
		     * member fuction of the spacetime implementation class
		     * template. Takes as arguments a vector reference to fill,
		     * a pointer to the momentum and a pointer to the mass: */

		    static vector_type& make_spinvec(vector_type& s,const vector_type* p,const typename spacetime_t::value_type* m)
		    {
			return spacetime_t::template make_KS_spinvec<beam_dir>(s,*p,*m);
		    }
	    };
    };
}

/* If the massive_spinor_factory class has been declared already,
 * appropriate specialisations are defined, for he purpose of optmisation: */

#ifdef CAMGEN_M_SPINOR_FAC_H_

#ifdef CAMGEN_PAULI_BASIS_H_
#include <Camgen/KSspinPb.h>
#endif /*CAMGEN_PAULI_BASIS_H_*/

#ifdef CAMGEN_WEYL_BASIS_H_
#include <Camgen/KSspinWb.h>
#endif /*CAMGEN_WEYL_BASIS_H_*/

#endif /*CAMGEN_M_SPINOR_FAC_H_*/

#endif /*CAMGEN_KS_TYPE_H_*/


