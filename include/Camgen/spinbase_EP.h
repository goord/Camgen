//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_EUCL_PAULI_BASIS_H_
#define CAMGEN_EUCL_PAULI_BASIS_H_

#include <Camgen/Euclidean.h>
#include <Camgen/Dirac_alg.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Euclidean Pauli basis type definition. It provides a matrix representation  *
 * of the Clifford algebra in Euclidean spacetimes. It can be derived from the *
 * Pauli basis in Minkowski space by multiplying all spacelike matrices with   *
 * the imaginary unit.                                                         *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class Euclidean_Pauli_basis
    {
	public:

	    /* Type definition of the corresponding spacetime-signature: */

	    typedef Euclidean_type spacetime_type;
	     
	    /* Implementation metafunction of the matrix basis: by default, this
	     * class is empty: */

	    template<class value_t,std::size_t dim,bool q=(dim%2==0)>class implementation: public Dirac_algebra<value_t,implementation<value_t,dim>,dim,q>{};
    };

    /* Implementation metafunction of the matrix basis in the four-dimensional
     * case. It is derived from the base class template Dirac
     * algebra<val_type,derived_type,dimension,diagonal_spacetime>: */

    template<class value_t>class Euclidean_Pauli_basis::implementation<value_t,4,true>: public Dirac_algebra<value_t,implementation<value_t,4,true>,4,true>
    {
	public:

	    /* Type definition of the corresponding spacetime implementation
	     * class: */

	    typedef Euclidean_type::template implementation<value_t,4> spacetime_type;
	    
	    /* Definition of the base class: */
	    
	    typedef Dirac_algebra<value_t,implementation<value_t,4,true>,4> base_type;

	    /* Utility type definitions: */

	    typedef value_t r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    
	    /* Signs of gamma matrices upon charge conjugation: */
	    
	    static const int eta_g=-1;
	    static const int eta_g5=1;
	    static const int eta_gg5=1;
	    static const int eta_gg=-1;
	    
	    /* Implementation of the gamma matrices: */
	    
	    static void fill_gamma_matrices()
	    {

		base_type::g[0][0][0]=value_type(-1,0);
		base_type::g[0][1][1]=value_type(-1,0);
		base_type::g[0][2][2]=value_type(1,0);
		base_type::g[0][3][3]=value_type(1,0);

		base_type::g[1][0][3]=value_type(0,-1);
		base_type::g[1][1][2]=value_type(0,-1);
		base_type::g[1][2][1]=value_type(0,1);
		base_type::g[1][3][0]=value_type(0,1);

		base_type::g[2][0][3]=value_type(-1,0);
		base_type::g[2][1][2]=value_type(1,0);
		base_type::g[2][2][1]=value_type(1,0);
		base_type::g[2][3][0]=value_type(-1,0);

		base_type::g[3][0][2]=value_type(0,-1);
		base_type::g[3][1][3]=value_type(0,1);
		base_type::g[3][2][0]=value_type(0,1);
		base_type::g[3][3][1]=value_type(0,-1);

	    }

	    /* Implementation of the chiral gamma matrix: */
	    
	    static void fill_gamma_5()
	    {
		
		base_type::g_5[0][2]=value_type(1,0);
		base_type::g_5[1][3]=value_type(1,0);
		base_type::g_5[2][0]=value_type(1,0);
		base_type::g_5[3][1]=value_type(1,0);
		
	    }

	    /* Implementation of the charge conjugation matrix: */
	    
	    static void fill_C_matrix()
	    {
		
		base_type::C_matrix[0][3]=value_type(0,1);
		base_type::C_matrix[1][2]=value_type(0,-1);
		base_type::C_matrix[2][1]=value_type(0,1);
		base_type::C_matrix[3][0]=value_type(0,-1);
		
	    }
    };
    template<class value_t>const int Euclidean_Pauli_basis::implementation<value_t,4,true>::eta_g;
    template<class value_t>const int Euclidean_Pauli_basis::implementation<value_t,4,true>::eta_g5;
    template<class value_t>const int Euclidean_Pauli_basis::implementation<value_t,4,true>::eta_gg5;
    template<class value_t>const int Euclidean_Pauli_basis::implementation<value_t,4,true>::eta_gg;
}
#endif /*CAMGEN_EUCL_PAULI_BASIS_H_*/

