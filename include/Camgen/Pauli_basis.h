//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PAULI_BASIS_H_
#define CAMGEN_PAULI_BASIS_H_

#include <Camgen/Minkowski.h>
#include <Camgen/Dirac_alg.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the Pauli basis of Dirac gamma matrices. Aside  *
 * the usual contruction functions of the matrices, all the spinor bilinears are *
 * overloaded to optimise the performance.                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class Pauli_basis
    {
	public:

	    /* Type definition of the corresponding spacetime-signature: */

	    typedef Minkowski_type spacetime_type;

	    /* Implementation metafunction of the matrix basis: by default, this
	     * class is empty: */

	    template<class value_t,std::size_t dim,bool q=(dim%2==0)>class implementation: public Dirac_algebra<value_t,implementation<value_t,dim>,dim,q>{};
    };

    /* Implementation metafunction of the matrix basis in the four-dimensional
     * case. It is derived from the base class template Dirac
     * algebra<val_type,derived_type,dimension,diagonal_spacetime>: */
    
    template<class value_t>class Pauli_basis::implementation<value_t,4,true>: public Dirac_algebra<value_t,implementation<value_t,4,true>,4,true>
    {
	public:

	    /* Type definition of the corresponding spacetime implementation
	     * class: */
	    
	    typedef Minkowski_type::template implementation<value_t,4> spacetime_type;
	    
	    /* Definition of the base class: */
	    
	    typedef Dirac_algebra<value_t,implementation<value_t,4,true>,4> base_type;
	    
	    /* Some useful type definitions, derived from the base class: */
	    
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::r_value_type r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;

	    /* Signs of gamma matrices upon charge conjugation: */
	    
	    static const int eta_g=-1;
	    static const int eta_g5=1;
	    static const int eta_gg5=1;
	    static const int eta_gg=-1;

	    /* Implementation of the gamma matrices: */
	    
	    static void fill_gamma_matrices()
	    {

		base_type::g[0][0][0]=value_type(1,0);
		base_type::g[0][1][1]=value_type(1,0);
		base_type::g[0][2][2]=value_type(-1,0);
		base_type::g[0][3][3]=value_type(-1,0);

		base_type::g[1][0][3]=value_type(1,0);
		base_type::g[1][1][2]=value_type(1,0);
		base_type::g[1][2][1]=value_type(-1,0);
		base_type::g[1][3][0]=value_type(-1,0);

		base_type::g[2][0][3]=value_type(0,-1);
		base_type::g[2][1][2]=value_type(0,1);
		base_type::g[2][2][1]=value_type(0,1);
		base_type::g[2][3][0]=value_type(0,-1);

		base_type::g[3][0][2]=value_type(1,0);
		base_type::g[3][1][3]=value_type(-1,0);
		base_type::g[3][2][0]=value_type(-1,0);
		base_type::g[3][3][1]=value_type(1,0);


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

		base_type::C_matrix[0][3]=value_type(-1,0);
		base_type::C_matrix[1][2]=value_type(1,0);
		base_type::C_matrix[2][1]=value_type(-1,0);
		base_type::C_matrix[3][0]=value_type(1,0);

	    }

	    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	     * In the following comments to the overloaded spinor vertex     *
	     * recursive relations, phi shall denote a scalar field, psi1 and*
	     * psi2 the two spinors and V a vector field,                    *
	     *                                                               *
	     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	    /* Basic spinor bilinear without prefactor x += bar{psi1}.psi2: */

	    static void I_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		x+=(A_1[0]*A_2[0]+A_1[1]*A_2[1]+A_1[2]*A_2[2]+A_1[3]*A_2[3]);
		
	    }

	    /* Spinor vertex recursive relation phi += c*bar{psi1}.psi2: */

	    static void I_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(A_1[0]*A_2[0]+A_1[1]*A_2[1]+A_1[2]*A_2[2]+A_1[3]*A_2[3]);

	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*psi2 */

	    static void I_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		
		A_1[0]+=c*A_2[0];
		A_1[1]+=c*A_2[1];
		A_1[2]+=c*A_2[2];
		A_1[3]+=c*A_2[3];

	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*phi*bar{psi1}: */

	    static void I_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		
		A_2[0]+=c*A_1[0];
		A_2[1]+=c*A_1[1];
		A_2[2]+=c*A_1[2];
		A_2[3]+=c*A_1[3];

	    }

	    /* Basic conjugate spinor bilinear without prefactor x +=
	     * bar{psi1}C.psi2: */

	    static void C_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		x+=(-A_1[0]*A_2[3]+A_1[1]*A_2[2]-A_1[2]*A_2[1]+A_1[3]*A_2[0]);

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}C.psi2: */

	    static void C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(-A_1[0]*A_2[3]+A_1[1]*A_2[2]-A_1[2]*A_2[1]+A_1[3]*A_2[0]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*phi*C.psi2:
	     * */

	    static void C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		
		A_1[0]-=c*A_2[3];
		A_1[1]+=c*A_2[2];
		A_1[2]-=c*A_2[1];
		A_1[3]+=c*A_2[0];

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi2*bar{psi1}.C: */

	    static void C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		
		A_2[0]+=c*A_1[3];
		A_2[1]-=c*A_1[2];
		A_2[2]+=c*A_1[1];
		A_2[3]-=c*A_1[0];

	    }

	    /* Left multiplication by C , psi1 += C.psi2: */

	    static void C_left(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		
		A_0[0]-=A_1[3];
		A_0[1]+=A_1[2];
		A_0[2]-=A_1[1];
		A_0[3]+=A_1[0];
		
	    }

	    /* Right multiplication by C, psi1 += psi2.C: */

	    static void C_right(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		
		A_0[0]+=A_1[3];
		A_0[1]-=A_1[2];
		A_0[2]+=A_1[1];
		A_0[3]-=A_1[0];
		
	    }

	    /* Basic conjugate spinor bilinear without prefactor x +=
	     * bar{psi1}.C^*.psi2: */
	    
	    static void Cc_first(value_type& x,const_iterator A_1,const_iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		x+=(-A_1[0]*A_2[3]+A_1[1]*A_2[2]-A_1[2]*A_2[1]+A_1[3]*A_2[0]);

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.psi2: */

	    static void Cc_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(-A_1[0]*A_2[3]+A_1[1]*A_2[2]-A_1[2]*A_2[1]+A_1[3]*A_2[0]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.psi2: */

	    static void Cc_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		
		A_1[0]-=c*A_2[3];
		A_1[1]+=c*A_2[2];
		A_1[2]-=c*A_2[1];
		A_1[3]+=c*A_2[0];

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*: */

	    static void Cc_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		
		A_2[0]+=c*A_1[3];
		A_2[1]-=c*A_1[2];
		A_2[2]+=c*A_1[1];
		A_2[3]-=c*A_1[0];

	    }

	    /* Left multiplication by C*, psi1 += C*.psi2: */

	    static void Cc_left(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		
		A_0[0]-=A_1[3];
		A_0[1]+=A_1[2];
		A_0[2]-=A_1[1];
		A_0[3]+=A_1[0];
		
	    }

	    /* Right multiplication by C*, psi1 += psi2.C*: */

	    static void Cc_right(iterator A_0,iterator A_1)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		
		A_0[0]+=A_1[3];
		A_0[1]-=A_1[2];
		A_0[2]+=A_1[1];
		A_0[3]-=A_1[0];
		
	    }

	    /* Spinor vertex recursive relation V(mu) += c*bar{psi1}.g(mu).psi2:
	     * */

	    static void g_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = A_1[0]*A_2[3] - A_1[2]*A_2[1];
		value_type Qm = A_1[1]*A_2[2] - A_1[3]*A_2[0];

		A_0[0] += C_0*(A_1[0]*A_2[0] + A_1[1]*A_2[1] - A_1[2]*A_2[2] - A_1[3]*A_2[3]);

		A_0[1] += C_0*(Qm + Qp);

		A_0[2] += C_0*(value_type(Qp.imag()-Qm.imag(),Qm.real()-Qp.real()));

		A_0[3] += C_0*(A_1[0]*A_2[2] - A_1[1]*A_2[3] - A_1[2]*A_2[0] + A_1[3]*A_2[1]);

	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.psi2: */

	    static void g_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_1[0] += C_0*(A_0[0]*A_2[0] - A_0[3]*A_2[2] - Qm*A_2[3]);

		A_1[1] += C_0*(A_0[0]*A_2[1] - Qp*A_2[2] + A_0[3]*A_2[3]);

		A_1[2] += C_0*(A_0[3]*A_2[0] + Qm*A_2[1] - A_0[0]*A_2[2]);

		A_1[3] += C_0*(Qp*A_2[0] - A_0[3]*A_2[1] - A_0[0]*A_2[3]);

	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.psi2: */

	    static void g_second(const value_type& C_0,const vector<r_value_type,4>& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1],A_0[2]);

		A_1[0] += C_0*(A_0[0]*A_2[0] - A_0[3]*A_2[2] - std::conj(Qp)*A_2[3]);

		A_1[1] += C_0*(A_0[0]*A_2[1] - Qp*A_2[2] + A_0[3]*A_2[3]);

		A_1[2] += C_0*(A_0[3]*A_2[0] + std::conj(Qp)*A_2[1] - A_0[0]*A_2[2]);

		A_1[3] += C_0*(Qp*A_2[0] - A_0[3]*A_2[1] - A_0[0]*A_2[3]);

	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*bar{psi1}.Vslash:
	     * */
	    
	    static void g_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_2[0] += C_0*(A_0[0]*A_1[0] + A_0[3]*A_1[2] + Qp*A_1[3]);

		A_2[1] += C_0*(A_0[0]*A_1[1] + Qm*A_1[2] - A_0[3]*A_1[3]);

		A_2[2] += C_0*(-A_0[3]*A_1[0] - Qp*A_1[1] - A_0[0]*A_1[2]);

		A_2[3] += C_0*(-Qm*A_1[0] + A_0[3]*A_1[1] - A_0[0]*A_1[3]);

	    }
	    
	    /* Spinor vertex recursive relation bar{psi2} += c*bar{psi1}.Vslash:
	     * */
	    
	    static void g_third(const value_type& C_0,const vector<r_value_type,4>& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1],A_0[2]);

		A_2[0] += C_0*(A_0[0]*A_1[0] + A_0[3]*A_1[2] + Qp*A_1[3]);

		A_2[1] += C_0*(A_0[0]*A_1[1] + std::conj(Qp)*A_1[2] - A_0[3]*A_1[3]);

		A_2[2] += C_0*(-A_0[3]*A_1[0] - Qp*A_1[1] - A_0[0]*A_1[2]);

		A_2[3] += C_0*(-std::conj(Qp)*A_1[0] + A_0[3]*A_1[1] - A_0[0]*A_1[3]);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).C.psi2: */

	    static void g_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = A_1[0]*A_2[0] - A_1[2]*A_2[2];
		value_type Qm = -A_1[1]*A_2[1] + A_1[3]*A_2[3];

		A_0[0] += C_0*(-A_1[0]*A_2[3] + A_1[1]*A_2[2] + A_1[2]*A_2[1] - A_1[3]*A_2[0]);

		A_0[1] += C_0*(Qm + Qp);

		A_0[2] += C_0*(value_type(Qp.imag()-Qm.imag(),Qm.real()-Qp.real()));

		A_0[3] += C_0*(-A_1[0]*A_2[1] - A_1[1]*A_2[0] + A_1[2]*A_2[3] + A_1[3]*A_2[2]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*Vslash.C.psi2: */

	    static void g_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_1[0] += C_0*(-A_0[0]*A_2[3] + A_0[3]*A_2[1] - Qm*A_2[0]);

		A_1[1] += C_0*(A_0[0]*A_2[2] + Qp*A_2[1] + A_0[3]*A_2[0]);

		A_1[2] += C_0*(-A_0[3]*A_2[3] + Qm*A_2[2] + A_0[0]*A_2[1]);

		A_1[3] += C_0*(-Qp*A_2[3] - A_0[3]*A_2[2] - A_0[0]*A_2[0]);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.C: */

	    static void g_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_2[0] += C_0*(-Qm*A_1[0] + A_0[3]*A_1[1] - A_0[0]*A_1[3]);

		A_2[1] -= C_0*(-A_0[3]*A_1[0] - Qp*A_1[1] - A_0[0]*A_1[2]);

		A_2[2] += C_0*(A_0[0]*A_1[1] + Qm*A_1[2] - A_0[3]*A_1[3]);

		A_2[3] -= C_0*(A_0[0]*A_1[0] + A_0[3]*A_1[2] + Qp*A_1[3]);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).psi2 */

	    static void Cc_g_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = A_1[3]*A_2[3] - A_1[1]*A_2[1];
		value_type Qm = -A_1[2]*A_2[2] + A_1[0]*A_2[0];

		A_0[0] += C_0*(A_1[3]*A_2[0] - A_1[2]*A_2[1] - A_1[1]*A_2[2] + A_1[0]*A_2[3]);

		A_0[1] += C_0*(Qm + Qp);

		A_0[2] += C_0*(value_type(Qp.imag()-Qm.imag(),Qm.real()-Qp.real()));

		A_0[3] += C_0*(A_1[3]*A_2[2] + A_1[2]*A_2[3] - A_1[1]*A_2[0] - A_1[0]*A_2[1]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 += c*C^*.Vslash.psi2: */

	    static void Cc_g_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_1[0] -= C_0*(Qp*A_2[0] - A_0[3]*A_2[1] - A_0[0]*A_2[3]);

		A_1[1] += C_0*(A_0[3]*A_2[0] + Qm*A_2[1] - A_0[0]*A_2[2]);

		A_1[2] -= C_0*(A_0[0]*A_2[1] - Qp*A_2[2] + A_0[3]*A_2[3]);

		A_1[3] += C_0*(A_0[0]*A_2[0] - A_0[3]*A_2[2] - Qm*A_2[3]);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*.Vslash: */

	    static void Cc_g_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_2[0] += C_0*(A_0[0]*A_1[3] + A_0[3]*A_1[1] - Qp*A_1[0]);

		A_2[1] += C_0*(-A_0[0]*A_1[2] + Qm*A_1[1] + A_0[3]*A_1[0]);

		A_2[2] += C_0*(-A_0[3]*A_1[3] + Qp*A_1[2] - A_0[0]*A_1[1]);

		A_2[3] += C_0*(-Qm*A_1[3] - A_0[3]*A_1[2] + A_0[0]*A_1[0]);	

	    }

	    /* Spinor vertex recursive relation phi += c*bar{psi1}.g5.psi2: */

	    static void g5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(A_1[0]*A_2[2]+A_1[1]*A_2[3]+A_1[2]*A_2[0]+A_1[3]*A_2[1]);

	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*g5.psi2: */

	    static void g5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		A_1[0]+=c*A_2[2];
		A_1[1]+=c*A_2[3];
		A_1[2]+=c*A_2[0];
		A_1[3]+=c*A_2[1];

	    }

	    /* Spinor vertex recursive relation bar{psi2} += c*phi*bar{psi1}.g5:
	     * */

	    static void g5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		A_2[0]+=c*A_1[2];
		A_2[1]+=c*A_1[3];
		A_2[2]+=c*A_1[0];
		A_2[3]+=c*A_1[1];

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.g5.C.psi2: */

	    static void g5_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(-A_1[0]*A_2[1]+A_1[1]*A_2[0]-A_1[2]*A_2[3]+A_1[3]*A_2[2]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*g5.C.psi2: */

	    static void g5_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		A_1[0]-=c*A_2[1];
		A_1[1]+=c*A_2[0];
		A_1[2]-=c*A_2[3];
		A_1[3]+=c*A_2[2];

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.g5.C: */

	    static void g5_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		A_2[0]+=c*A_1[1];
		A_2[1]-=c*A_1[0];
		A_2[2]+=c*A_1[3];
		A_2[3]-=c*A_1[2];

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.g5.psi2: */

	    static void Cc_g5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(A_1[3]*A_2[2]-A_1[2]*A_2[3]+A_1[1]*A_2[0]-A_1[0]*A_2[1]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.g5.psi2: */

	    static void Cc_g5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		A_1[3]+=c*A_2[2];
		A_1[2]-=c*A_2[3];
		A_1[1]+=c*A_2[0];
		A_1[0]-=c*A_2[1];

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.g5: */

	    static void Cc_g5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c=C_0*A_0[0];
		A_2[0]+=c*A_1[1];
		A_2[1]-=c*A_1[0];
		A_2[2]+=c*A_1[3];
		A_2[3]-=c*A_1[2];

	    }

	    /* Spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(c0+c1*g5).psi2: */

	    static void gVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(A_1[0]*A_2[0]+A_1[1]*A_2[1]+A_1[2]*A_2[2]+A_1[3]*A_2[3]);
		A_0[0]+=C_1*(A_1[0]*A_2[2]+A_1[1]*A_2[3]+A_1[2]*A_2[0]+A_1[3]*A_2[1]);

	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*(c0+c1*g5).psi2:
	     * */

	    static void gVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=C_0*A_0[0];
		value_type c2=C_1*A_0[0];

		A_1[0]+=(c1*A_2[0] + c2*A_2[2]);
		A_1[1]+=(c1*A_2[1] + c2*A_2[3]);
		A_1[2]+=(c2*A_2[0] + c1*A_2[2]);
		A_1[3]+=(c2*A_2[1] + c1*A_2[3]);

	    }

	    /* spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(c0+c1*g5): */

	    static void gVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=C_0*A_0[0];
		value_type c2=C_1*A_0[0];

		A_2[0]+=(c1*A_1[0] + c2*A_1[2]);
		A_2[1]+=(c1*A_1[1] + c2*A_1[3]);
		A_2[2]+=(c2*A_1[0] + c1*A_1[2]);
		A_2[3]+=(c2*A_1[1] + c1*A_1[3]);

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(c0+c1*g5).C.psi2: */
	    
	    static void gVA_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(-A_1[0]*A_2[3]+A_1[1]*A_2[2]-A_1[2]*A_2[1]+A_1[3]*A_2[0]);
		A_0[0]+=C_1*(-A_1[0]*A_2[1]+A_1[1]*A_2[0]-A_1[2]*A_2[3]+A_1[3]*A_2[2]);

	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*(c0+c1*g5).C.psi2: */
	    
	    static void gVA_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=C_0*A_0[0];
		value_type c2=C_1*A_0[0];

		A_1[0]-=(c1*A_2[3] + c2*A_2[1]);
		A_1[1]+=(c1*A_2[2] + c2*A_2[0]);
		A_1[2]-=(c2*A_2[3] + c1*A_2[1]);
		A_1[3]+=(c2*A_2[2] + c1*A_2[0]);

	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(c0+c1*g5).C: */
	    
	    static void gVA_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=C_0*A_0[0];
		value_type c2=C_1*A_0[0];

		A_2[0]+=(c2*A_1[1] + c1*A_1[3]);
		A_2[1]-=(c2*A_1[0] + c1*A_1[2]);
		A_2[2]+=(c1*A_1[1] + c2*A_1[3]);
		A_2[3]-=(c1*A_1[0] + c2*A_1[2]);

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.(c0+c1*g5).psi2: */
	    
	    static void Cc_gVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=C_0*(A_1[3]*A_2[0]-A_1[2]*A_2[1]+A_1[1]*A_2[2]-A_1[0]*A_2[3]);
		A_0[0]+=C_1*(A_1[3]*A_2[2]-A_1[2]*A_2[3]+A_1[1]*A_2[0]-A_1[0]*A_2[1]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.(c0+c1*g5).psi2: */
	    
	    static void Cc_gVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=C_0*A_0[0];
		value_type c2=C_1*A_0[0];

		A_1[0]-=(c2*A_2[1] + c1*A_2[3]);
		A_1[1]+=(c2*A_2[0] + c1*A_2[2]);
		A_1[2]-=(c1*A_2[1] + c2*A_2[3]);
		A_1[3]+=(c1*A_2[0] + c2*A_2[2]);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.(c0+c1*g5): */
	    
	    static void Cc_gVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=C_0*A_0[0];
		value_type c2=C_1*A_0[0];

		A_2[0]+=(c1*A_1[3] + c2*A_1[1]);
		A_2[1]-=(c1*A_1[2] + c2*A_1[0]);
		A_2[2]+=(c2*A_1[3] + c1*A_1[1]);
		A_2[3]-=(c2*A_1[2] + c1*A_1[0]);

	    }

	    /* Spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1+g5/2).psi2: */
	    
	    static void gL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*C_0*(A_1[0]*(A_2[0]+A_2[2])+A_1[1]*(A_2[1]+A_2[3])+A_1[2]*(A_2[2]+A_2[0])+A_1[3]*(A_2[3]+A_2[1]));

	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*(1+g5/2).psi2: */
	    
	    static void gL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_2[0]+A_2[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_2[1]+A_2[3]);

		A_1[0]+=c1;
		A_1[1]+=c2;
		A_1[2]+=c1;
		A_1[3]+=c2;

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1+g5/2): */

	    static void gL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_1[0]+A_1[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_1[1]+A_1[3]);

		A_2[0]+=c1;
		A_2[1]+=c2;
		A_2[2]+=c1;
		A_2[3]+=c2;

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1+g5/2).C.psi2: */

	    static void gL_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*C_0*(-A_1[0]*(A_2[3]+A_2[1])+A_1[1]*(A_2[2]+A_2[0])-A_1[2]*(A_2[1]+A_2[3])+A_1[3]*(A_2[0]+A_2[2]));

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*(1+g5/2).C.psi2: */
	    
	    static void gL_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_2[1]+A_2[3]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_2[0]+A_2[2]);

		A_1[0]-=c1;
		A_1[1]+=c2;
		A_1[2]-=c1;
		A_1[3]+=c2;

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1+g5/2).C: */
	    
	    static void gL_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_1[0]+A_1[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_1[1]+A_1[3]);

		A_2[0]+=c2;
		A_2[1]-=c1;
		A_2[2]+=c2;
		A_2[3]-=c1;

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.(1+g5/2).psi2: */
	    
	    static void Cc_gL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*C_0*(A_1[3]*(A_2[0]+A_2[2])-A_1[2]*(A_2[1]+A_2[3])+A_1[1]*(A_2[2]+A_2[0])-A_1[0]*(A_2[3]+A_2[1]));

	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.(1+g5/2).psi2: */
	    
	    static void Cc_gL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_2[0]+A_2[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_2[1]+A_2[3]);

		A_1[0]-=c2;
		A_1[1]+=c1;
		A_1[2]-=c2;
		A_1[3]+=c1;

	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.(1+g5/2): */
	    
	    static void Cc_gL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_1[1]+A_1[3]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_1[0]+A_1[2]);

		A_2[0]+=c1;
		A_2[1]-=c2;
		A_2[2]+=c1;
		A_2[3]-=c2;

	    }

	    /* Spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1-g5/2).psi2: */
	    
	    static void gR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*C_0*(A_1[0]*(A_2[0]-A_2[2])+A_1[1]*(A_2[1]-A_2[3])+A_1[2]*(A_2[2]-A_2[0])+A_1[3]*(A_2[3]-A_2[1]));

	    }

	    /* Spinor vertex recursive relation psi1 += c*phi*(1-g5/2).psi2: */
	    
	    static void gR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_2[0]-A_2[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_2[1]-A_2[3]);

		A_1[0]+=c1;
		A_1[1]+=c2;
		A_1[2]-=c1;
		A_1[3]-=c2;

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1-g5/2): */
	    
	    static void gR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_1[0]-A_1[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_1[1]-A_1[3]);

		A_2[0]+=c1;
		A_2[1]+=c2;
		A_2[2]-=c1;
		A_2[3]-=c2;

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.(1-g5/2).C.psi2: */
	    
	    static void gR_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*C_0*(A_1[0]*(A_2[1]-A_2[3])+A_1[1]*(A_2[2]-A_2[0])+A_1[2]*(A_2[3]-A_2[1])+A_1[3]*(A_2[0]-A_2[2]));

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*(1-g5/2).C.psi2: */
	    
	    static void gR_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_2[1]-A_2[3]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_2[2]-A_2[0]);

		A_1[0]+=c1;
		A_1[1]+=c2;
		A_1[2]-=c1;
		A_1[3]-=c2;

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.(1-g5/2).C: */
	    
	    static void gR_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_1[0]-A_1[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_1[1]-A_1[3]);

		A_2[0]-=c2;
		A_2[1]+=c1;
		A_2[2]+=c2;
		A_2[3]-=c1;

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * c*bar{psi1}.C^*.(1-g5/2).psi2: */
	    
	    static void Cc_gR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*C_0*(A_1[3]*(A_2[0]-A_2[2])-A_1[2]*(A_2[1]-A_2[3])+A_1[1]*(A_2[2]-A_2[0])-A_1[0]*(A_2[3]-A_2[1]));

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*phi*C^*.(1-g5/2).psi2: */
	    
	    static void Cc_gR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_2[0]-A_2[2]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_2[1]-A_2[3]);

		A_1[0]+=c2;
		A_1[1]-=c1;
		A_1[2]-=c2;
		A_1[3]+=c1;

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*phi*bar{psi1}.C^*.(1-g5/2): */
	    
	    static void Cc_gR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*C_0*A_0[0]*(A_1[3]-A_1[1]);
		value_type c2=(r_value_type)0.5*C_0*A_0[0]*(A_1[0]-A_1[2]);

		A_2[0]+=c1;
		A_2[1]+=c2;
		A_2[2]-=c1;
		A_2[3]-=c2;

	    }

	    /* Spinor vertex recursive relation phi +=
	     * bar{psi1}.(c0*PL+c1*PR).psi2: */
	    
	    static void gLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*(C_0+C_1)*(A_1[0]*A_2[0]+A_1[1]*A_2[1]+A_1[2]*A_2[2]+A_1[3]*A_2[3]);
		A_0[0]+=(r_value_type)0.5*(C_0-C_1)*(A_1[0]*A_2[2]+A_1[1]*A_2[3]+A_1[2]*A_2[0]+A_1[3]*A_2[1]);

	    }

	    /* Spinor vertex recursive relation psi1 += phi*(c0*PL+c1*PR).psi2:
	     * */

	    static void gLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*(C_0+C_1)*A_0[0];
		value_type c2=(r_value_type)0.5*(C_0-C_1)*A_0[0];

		A_1[0]+=(c1*A_2[0] + c2*A_2[2]);
		A_1[1]+=(c1*A_2[1] + c2*A_2[3]);
		A_1[2]+=(c2*A_2[0] + c1*A_2[2]);
		A_1[3]+=(c2*A_2[1] + c1*A_2[3]);

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * phi*bar{psi1}.(c0*PL+c1*PR): */

	    static void gLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*(C_0+C_1)*A_0[0];
		value_type c2=(r_value_type)0.5*(C_0-C_1)*A_0[0];

		A_2[0]+=(c1*A_1[0] + c2*A_1[2]);
		A_2[1]+=(c1*A_1[1] + c2*A_1[3]);
		A_2[2]+=(c2*A_1[0] + c1*A_1[2]);
		A_2[3]+=(c2*A_1[1] + c1*A_1[3]);

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * bar{psi1}.(c0*PL+c1*PR).C.psi2: */
	    
	    static void gLR_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*(C_0+C_1)*(-A_1[0]*A_2[3]+A_1[1]*A_2[2]-A_1[2]*A_2[1]+A_1[3]*A_2[0]);
		A_0[0]+=(r_value_type)0.5*(C_0-C_1)*(-A_1[0]*A_2[1]+A_1[1]*A_2[0]-A_1[2]*A_2[3]+A_1[3]*A_2[2]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * phi*(c0*PL+c1*PR).C.psi2: */

	    static void gLR_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*(C_0+C_1)*A_0[0];
		value_type c2=(r_value_type)0.5*(C_0-C_1)*A_0[0];

		A_1[0]-=(c1*A_2[3] + c2*A_2[1]);
		A_1[1]+=(c1*A_2[2] + c2*A_2[0]);
		A_1[2]-=(c2*A_2[3] + c1*A_2[1]);
		A_1[3]+=(c2*A_2[2] + c1*A_2[0]);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * phi*bar{psi1}.(c0*PL+c1*PR).C: */

	    static void gLR_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*(C_0+C_1)*A_0[0];
		value_type c2=(r_value_type)0.5*(C_0-C_1)*A_0[0];

		A_2[0]+=(c2*A_1[1] + c1*A_1[3]);
		A_2[1]-=(c2*A_1[0] + c1*A_1[2]);
		A_2[2]+=(c1*A_1[1] + c2*A_1[3]);
		A_2[3]-=(c1*A_1[0] + c2*A_1[2]);

	    }

	    /* Conjugate spinor vertex recursive relation phi +=
	     * bar{psi1}.C^*.(c0*PL+c1*PR).psi2: */
	    
	    static void Cc_gLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		A_0[0]+=(r_value_type)0.5*(C_0+C_1)*(A_1[3]*A_2[0]-A_1[2]*A_2[1]+A_1[1]*A_2[2]-A_1[0]*A_2[3]);
		A_0[0]+=(r_value_type)0.5*(C_0-C_1)*(A_1[3]*A_2[2]-A_1[2]*A_2[3]+A_1[1]*A_2[0]-A_1[0]*A_2[1]);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * phi*C^*.(c0*PL+c1*PR).C.psi2: */

	    static void Cc_gLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*(C_0+C_1)*A_0[0];
		value_type c2=(r_value_type)0.5*(C_0-C_1)*A_0[0];

		A_1[0]-=(c2*A_2[1] + c1*A_2[3]);
		A_1[1]+=(c2*A_2[0] + c1*A_2[2]);
		A_1[2]-=(c1*A_2[1] + c2*A_2[3]);
		A_1[3]+=(c1*A_2[0] + c2*A_2[2]);

	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * phi*bar{psi1}.C^*.(c0*PL+c1*PR): */
	    
	    static void Cc_gLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<1),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type c1=(r_value_type)0.5*(C_0+C_1)*A_0[0];
		value_type c2=(r_value_type)0.5*(C_0-C_1)*A_0[0];

		A_2[0]+=(c1*A_1[3] + c2*A_1[1]);
		A_2[1]-=(c1*A_1[2] + c2*A_1[0]);
		A_2[2]+=(c2*A_1[3] + c1*A_1[1]);
		A_2[3]-=(c2*A_1[2] + c1*A_1[0]);

	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).g5.psi2: */

	    static void gg5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = A_1[0]*A_2[1] - A_1[2]*A_2[3];
		value_type Qm = A_1[1]*A_2[0] - A_1[3]*A_2[2];

		A_0[0] += C_0*(A_1[0]*A_2[2] + A_1[1]*A_2[3] - A_1[2]*A_2[0] - A_1[3]*A_2[1]);

		A_0[1] += C_0*(Qm + Qp);

		A_0[2] += C_0*value_type(Qp.imag()-Qm.imag(),Qm.real()-Qp.real());

		A_0[3] += C_0*(A_1[0]*A_2[0] - A_1[1]*A_2[1] - A_1[2]*A_2[2] + A_1[3]*A_2[3]);

	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.g5.psi2: */

	    static void gg5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_1[0] += C_0*(-A_0[3]*A_2[0] + A_0[0]*A_2[2] - Qm*A_2[1]);

		A_1[1] += C_0*(A_0[0]*A_2[3] - Qp*A_2[0] + A_0[3]*A_2[1]);

		A_1[2] += C_0*(-A_0[0]*A_2[0] + Qm*A_2[3] + A_0[3]*A_2[2]);

		A_1[3] += C_0*(Qp*A_2[2] - A_0[3]*A_2[3] - A_0[0]*A_2[1]);

	    }

	    /* Spinor vertex recursive relation psi1 += c*Vslash.g5.psi2: */

	    static void gg5_second(const value_type& C_0,const vector<r_value_type,4>& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1],A_0[2]);

		A_1[0] += C_0*(-A_0[3]*A_2[0] + A_0[0]*A_2[2] - std::conj(Qp)*A_2[1]);

		A_1[1] += C_0*(A_0[0]*A_2[3] - Qp*A_2[0] + A_0[3]*A_2[1]);

		A_1[2] += C_0*(-A_0[0]*A_2[0] + std::conj(Qp)*A_2[3] + A_0[3]*A_2[2]);

		A_1[3] += C_0*(Qp*A_2[2] - A_0[3]*A_2[3] - A_0[0]*A_2[1]);

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.g5: */

	    static void gg5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_2[0] += C_0*(-A_0[0]*A_1[2] - A_0[3]*A_1[0] - Qp*A_1[1]);

		A_2[1] += C_0*(-A_0[0]*A_1[3] - Qm*A_1[0] + A_0[3]*A_1[1]);

		A_2[2] += C_0*(A_0[0]*A_1[0] + Qp*A_1[3] + A_0[3]*A_1[2]);

		A_2[3] += C_0*(Qm*A_1[2] + A_0[0]*A_1[1] - A_0[3]*A_1[3]);

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.g5: */

	    static void gg5_third(const value_type& C_0,const vector<r_value_type,4>& A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1],A_0[2]);

		A_2[0] += C_0*(-A_0[0]*A_1[2] - A_0[3]*A_1[0] - Qp*A_1[1]);

		A_2[1] += C_0*(-A_0[0]*A_1[3] - std::conj(Qp)*A_1[0] + A_0[3]*A_1[1]);

		A_2[2] += C_0*(A_0[0]*A_1[0] + Qp*A_1[3] + A_0[3]*A_1[2]);

		A_2[3] += C_0*(std::conj(Qp)*A_1[2] + A_0[0]*A_1[1] - A_0[3]*A_1[3]);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).g5.C.psi2: */
	    
	    static void gg5_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = A_1[0]*A_2[2] - A_1[2]*A_2[0];
		value_type Qm = -A_1[1]*A_2[3] + A_1[3]*A_2[1];

		A_0[0] += C_0*(-A_1[0]*A_2[1] + A_1[1]*A_2[0] + A_1[2]*A_2[3] - A_1[3]*A_2[2]);

		A_0[1] += C_0*(Qm + Qp);

		A_0[2] += C_0*value_type(Qp.imag()-Qm.imag(),Qm.real()-Qp.real());

		A_0[3] += C_0*(-A_1[0]*A_2[3] - A_1[1]*A_2[2] + A_1[2]*A_2[1] + A_1[3]*A_2[0]);

	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.g5.C.psi2: */
	    
	    static void gg5_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_1[0] += C_0*(A_0[3]*A_2[3] - A_0[0]*A_2[1] - Qm*A_2[2]);

		A_1[1] += C_0*(A_0[0]*A_2[0] + Qp*A_2[3] + A_0[3]*A_2[2]);

		A_1[2] += C_0*(A_0[0]*A_2[3] + Qm*A_2[0] - A_0[3]*A_2[1]);

		A_1[3] += C_0*(-Qp*A_2[1] - A_0[3]*A_2[0] - A_0[0]*A_2[2]);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.g5.C: */

	    static void gg5_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_2[0] += C_0*(Qm*A_1[2] + A_0[0]*A_1[1] - A_0[3]*A_1[3]);

		A_2[1] -= C_0*(A_0[0]*A_1[0] + Qp*A_1[3] + A_0[3]*A_1[2]);

		A_2[2] += C_0*(-A_0[0]*A_1[3] - Qm*A_1[0] + A_0[3]*A_1[1]);

		A_2[3] -= C_0*(-A_0[0]*A_1[2] - A_0[3]*A_1[0] - Qp*A_1[1]);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).g5.psi2: */

	    static void Cc_gg5_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = A_1[3]*A_2[1] - A_1[1]*A_2[3];
		value_type Qm = -A_1[2]*A_2[0] + A_1[0]*A_2[2];

		A_0[0] += C_0*(A_1[3]*A_2[2] - A_1[2]*A_2[3] - A_1[1]*A_2[0] + A_1[0]*A_2[1]);

		A_0[1] += C_0*(Qm + Qp);

		A_0[2] += C_0*value_type(Qp.imag()-Qm.imag(),Qm.real()-Qp.real());

		A_0[3] += C_0*(A_1[3]*A_2[0] + A_1[2]*A_2[1] - A_1[1]*A_2[2] - A_1[0]*A_2[3]);

	    }
	    
	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*C^*.Vslash.g5.psi2: */
	    
	    static void Cc_gg5_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_1[0] -= C_0*(Qp*A_2[2] - A_0[3]*A_2[3] - A_0[0]*A_2[1]);

		A_1[1] += C_0*(-A_0[0]*A_2[0] + Qm*A_2[3] + A_0[3]*A_2[2]);

		A_1[2] -= C_0*(A_0[0]*A_2[3] - Qp*A_2[0] + A_0[3]*A_2[1]);

		A_1[3] += C_0*(-A_0[3]*A_2[0] + A_0[0]*A_2[2] - Qm*A_2[1]);

	    }
	    
	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.C^*.Vslash.g5: */
	    
	    static void Cc_gg5_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real());
		value_type Qm(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real());

		A_2[0] += C_0*(-A_0[0]*A_1[1] - A_0[3]*A_1[3] + Qp*A_1[2]);

		A_2[1] += C_0*(A_0[0]*A_1[0] - Qm*A_1[3] - A_0[3]*A_1[2]);

		A_2[2] += C_0*(A_0[0]*A_1[3] - Qp*A_1[0] + A_0[3]*A_1[1]);

		A_2[3] += C_0*(Qm*A_1[1] - A_0[0]*A_1[2] + A_0[3]*A_1[0]);

	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0+c1*g5).psi2: */

	    static void ggVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = C_0*A_2[0] + C_1*A_2[2];

		value_type psi20 = C_1*A_2[0] + C_0*A_2[2];

		value_type psi13 = C_0*A_2[1] + C_1*A_2[3];

		value_type psi31 = C_1*A_2[1] + C_0*A_2[3];

		A_0[0] += (A_1[0]*psi02 + A_1[1]*psi13 - A_1[2]*psi20 - A_1[3]*psi31);

		A_0[1] += (A_1[0]*psi31 + A_1[1]*psi20 - A_1[2]*psi13 - A_1[3]*psi02);

		A_0[2] += times_i(-A_1[0]*psi31 + A_1[1]*psi20 + A_1[2]*psi13 - A_1[3]*psi02);

		A_0[3] += (A_1[0]*psi20 - A_1[1]*psi31 - A_1[2]*psi02 + A_1[3]*psi13);

	    }

	    /* Spinor vertex recursive relation psi1 += Vslash.(c0+c1*g5).psi2:
	     * */

	    static void ggVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = C_0*A_2[0] + C_1*A_2[2];

		value_type psi20 = C_1*A_2[0] + C_0*A_2[2];

		value_type psi13 = C_0*A_2[1] + C_1*A_2[3];

		value_type psi31 = C_1*A_2[1] + C_0*A_2[3];

		A_1[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi31);

		A_1[1] += (-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_1[2] += (-A_0[0]*psi20 + A_0[3]*psi02 + value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi13);

		A_1[3] += (value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi02 - A_0[0]*psi31 - A_0[3]*psi13);

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0+c1*g5): */

	    static void ggVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = C_0*A_1[0] - C_1*A_1[2];

		value_type psi20 = C_1*A_1[0] - C_0*A_1[2];

		value_type psi13 = C_0*A_1[1] - C_1*A_1[3];

		value_type psi31 = C_1*A_1[1] - C_0*A_1[3];

		A_2[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi31);

		A_2[1] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_2[2] += (A_0[0]*psi20 - A_0[3]*psi02 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi13);

		A_2[3] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi02 + A_0[0]*psi31 + A_0[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0+c1*g5).C.psi2: */
	    
	    static void ggVA_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = -C_0*A_2[3] - C_1*A_2[1];

		value_type psi20 = -C_1*A_2[3] - C_0*A_2[1];

		value_type psi13 = C_0*A_2[2] + C_1*A_2[0];

		value_type psi31 = C_1*A_2[2] + C_0*A_2[0];

		A_0[0] += (A_1[0]*psi02 + A_1[1]*psi13 - A_1[2]*psi20 - A_1[3]*psi31);

		A_0[1] += (A_1[0]*psi31 + A_1[1]*psi20 - A_1[2]*psi13 - A_1[3]*psi02);

		A_0[2] += times_i(-A_1[0]*psi31 + A_1[1]*psi20 + A_1[2]*psi13 - A_1[3]*psi02);

		A_0[3] += (A_1[0]*psi20 - A_1[1]*psi31 - A_1[2]*psi02 + A_1[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * Vslash.(c0+c1*g5).C.psi2: */
	    
	    static void ggVA_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = -C_0*A_2[3] - C_1*A_2[1];

		value_type psi20 = -C_1*A_2[3] - C_0*A_2[1];

		value_type psi13 = C_0*A_2[2] + C_1*A_2[0];

		value_type psi31 = C_1*A_2[2] + C_0*A_2[0];

		A_1[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi31);

		A_1[1] += (-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_1[2] += (-A_0[0]*psi20 + A_0[3]*psi02 + value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi13);

		A_1[3] += (value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi02 - A_0[0]*psi31 - A_0[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0+c1*g5).C: */

	    static void ggVA_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = C_0*A_1[0] - C_1*A_1[2];

		value_type psi20 = C_1*A_1[0] - C_0*A_1[2];

		value_type psi13 = C_0*A_1[1] - C_1*A_1[3];

		value_type psi31 = C_1*A_1[1] - C_0*A_1[3];

		A_2[3] -= (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi31);

		A_2[2] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_2[1] -= (A_0[0]*psi20 - A_0[3]*psi02 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi13);

		A_2[0] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi02 + A_0[0]*psi31 + A_0[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.C^*.g(mu).(c0+c1*g5).psi2: */
	    
	    static void Cc_ggVA_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = C_0*A_2[0] + C_1*A_2[2];

		value_type psi20 = C_1*A_2[0] + C_0*A_2[2];

		value_type psi13 = C_0*A_2[1] + C_1*A_2[3];

		value_type psi31 = C_1*A_2[1] + C_0*A_2[3];

		A_0[0] += (A_1[3]*psi02 - A_1[2]*psi13 - A_1[1]*psi20 + A_1[0]*psi31);

		A_0[1] += (A_1[3]*psi31 - A_1[2]*psi20 - A_1[1]*psi13 + A_1[0]*psi02);

		A_0[2] += times_i(-A_1[3]*psi31 - A_1[2]*psi20 + A_1[1]*psi13 + A_1[0]*psi02);

		A_0[3] += (A_1[3]*psi20 + A_1[2]*psi31 - A_1[1]*psi02 - A_1[0]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * C^*.Vslash.(c0+c1*g5).psi2: */

	    static void Cc_ggVA_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = C_0*A_2[0] + C_1*A_2[2];

		value_type psi20 = C_1*A_2[0] + C_0*A_2[2];

		value_type psi13 = C_0*A_2[1] + C_1*A_2[3];

		value_type psi31 = C_1*A_2[1] + C_0*A_2[3];

		A_1[0] -= (value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi02 - A_0[0]*psi31 - A_0[3]*psi13);

		A_1[1] += (-A_0[0]*psi20 + A_0[3]*psi02 + value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi13);

		A_1[2] -= (-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_1[3] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi31);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.C^*.Vslash.(c0+c1*g5): */

	    static void Cc_ggVA_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = C_0*A_1[3] - C_1*A_1[1];

		value_type psi20 = C_1*A_1[3] - C_0*A_1[1];

		value_type psi13 = -C_0*A_1[2] + C_1*A_1[0];

		value_type psi31 = -C_1*A_1[2] + C_0*A_1[0];

		A_2[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi31);

		A_2[1] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_2[2] += (A_0[0]*psi20 - A_0[3]*psi02 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi13);

		A_2[3] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi02 + A_0[0]*psi31 + A_0[3]*psi13);

	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1+g5/2).psi2: */

	    static void ggL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = (r_value_type)0.5*C_0*((A_1[0]-A_1[2])*(A_2[0]+A_2[2]));

		value_type Qm = (r_value_type)0.5*C_0*((A_1[1]-A_1[3])*(A_2[1]+A_2[3]));

		value_type Pp = (r_value_type)0.5*C_0*((A_1[1]-A_1[3])*(A_2[0]+A_2[2]));

		value_type Pm = (r_value_type)0.5*C_0*((A_1[0]-A_1[2])*(A_2[1]+A_2[3]));

		A_0[0] += (Qp + Qm);

		A_0[1] += (Pp + Pm);

		A_0[2] += value_type(Pm.imag()-Pp.imag(),Pp.real()-Pm.real());

		A_0[3] += (Qp - Qm);

	    }

	    /* Spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1+g5/2).psi2: */

	    static void ggL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]-A_0[3])*(A_2[0]+A_2[2])-value_type(A_0[1].real() + A_0[2].imag(),A_0[1].imag()-A_0[2].real())*(A_2[1]+A_2[3]));

		value_type Pm = (r_value_type)0.5*C_0*(-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag() + A_0[2].real())*(A_2[0]+A_2[2])+(A_0[0]+A_0[3])*(A_2[1]+A_2[3]));

		A_1[0] += Pp;

		A_1[1] += Pm;

		A_1[2] -= Pp;

		A_1[3] -= Pm;

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1+g5/2): */

	    static void ggL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]-A_0[3])*(A_1[0]-A_1[2])-value_type(A_0[1].real() - A_0[2].imag(),A_0[1].imag()+A_0[2].real())*(A_1[1]-A_1[3]));

		value_type Pm = (r_value_type)0.5*C_0*(-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag() - A_0[2].real())*(A_1[0]-A_1[2])+(A_0[0]+A_0[3])*(A_1[1]-A_1[3]));

		A_2[0] += Pp;

		A_2[1] += Pm;

		A_2[2] += Pp;

		A_2[3] += Pm;

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1+g5/2).C.psi2: */

	    static void ggL_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = -(r_value_type)0.5*C_0*((A_1[0]-A_1[2])*(A_2[1]+A_2[3]));

		value_type Qm = (r_value_type)0.5*C_0*((A_1[1]-A_1[3])*(A_2[0]+A_2[2]));

		value_type Pp = -(r_value_type)0.5*C_0*((A_1[1]-A_1[3])*(A_2[1]+A_2[3]));

		value_type Pm = (r_value_type)0.5*C_0*((A_1[0]-A_1[2])*(A_2[0]+A_2[2]));

		A_0[0] += (Qp + Qm);

		A_0[1] += (Pp + Pm);

		A_0[2] += value_type(Pm.imag()-Pp.imag(),Pp.real()-Pm.real());

		A_0[3] += (Qp - Qm);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1+g5/2).C.psi2: */

	    static void ggL_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*(-(A_0[0]-A_0[3])*(A_2[1]+A_2[3])-value_type(A_0[1].real() + A_0[2].imag(),A_0[1].imag()-A_0[2].real())*(A_2[0]+A_2[2]));

		value_type Pm = (r_value_type)0.5*C_0*(value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag() + A_0[2].real())*(A_2[1]+A_2[3])+(A_0[0]+A_0[3])*(A_2[0]+A_2[2]));

		A_1[0] += Pp;

		A_1[1] += Pm;

		A_1[2] -= Pp;

		A_1[3] -= Pm;

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1+g5/2).C: */

	    static void ggL_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]-A_0[3])*(A_1[0]-A_1[2])-value_type(A_0[1].real() - A_0[2].imag(),A_0[1].imag()+A_0[2].real())*(A_1[1]-A_1[3]));

		value_type Pm = (r_value_type)0.5*C_0*(-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag() - A_0[2].real())*(A_1[0]-A_1[2])+(A_0[0]+A_0[3])*(A_1[1]-A_1[3]));

		A_2[0] += Pm;

		A_2[1] -= Pp;

		A_2[2] += Pm;

		A_2[3] -= Pp;

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).(1+g5/2).psi2: */
	    
	    static void Cc_ggL_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = (r_value_type)0.5*C_0*((A_1[3]-A_1[1])*(A_2[0]+A_2[2]));

		value_type Qm = (r_value_type)0.5*C_0*((A_1[0]-A_1[2])*(A_2[1]+A_2[3]));

		value_type Pp = (r_value_type)0.5*C_0*((A_1[0]-A_1[2])*(A_2[0]+A_2[2]));

		value_type Pm = (r_value_type)0.5*C_0*((A_1[3]-A_1[1])*(A_2[1]+A_2[3]));

		A_0[0] += (Qp + Qm);

		A_0[1] += (Pp + Pm);

		A_0[2] += value_type(Pm.imag()-Pp.imag(),Pp.real()-Pm.real());

		A_0[3] += (Qp - Qm);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.C^*.(1+g5/2).psi2: */

	    static void Cc_ggL_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]-A_0[3])*(A_2[0]+A_2[2])-value_type(A_0[1].real() + A_0[2].imag(),A_0[1].imag()-A_0[2].real())*(A_2[1]+A_2[3]));

		value_type Pm = (r_value_type)0.5*C_0*(-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag() + A_0[2].real())*(A_2[0]+A_2[2])+(A_0[0]+A_0[3])*(A_2[1]+A_2[3]));

		A_1[3] += Pp;

		A_1[2] -= Pm;

		A_1[1] -= Pp;

		A_1[0] += Pm;

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.C^*.(1+g5/2): */

	    static void Cc_ggL_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]-A_0[3])*(A_1[3]-A_1[1])+value_type(A_0[1].real() - A_0[2].imag(),A_0[1].imag()+A_0[2].real())*(A_1[2]-A_1[0]));

		value_type Pm = (r_value_type)0.5*C_0*(-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag() - A_0[2].real())*(A_1[3]-A_1[1])+(A_0[0]+A_0[3])*(A_1[0]-A_1[2]));

		A_2[0] += Pp;

		A_2[1] += Pm;

		A_2[2] += Pp;

		A_2[3] += Pm;

	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1-g5/2).psi2: */
	    
	    static void ggR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = (r_value_type)0.5*C_0*((A_1[0]+A_1[2])*(A_2[0]-A_2[2]));

		value_type Qm = (r_value_type)0.5*C_0*((A_1[1]+A_1[3])*(A_2[1]-A_2[3]));

		value_type Pp = (r_value_type)0.5*C_0*((A_1[1]+A_1[3])*(A_2[0]-A_2[2]));

		value_type Pm = (r_value_type)0.5*C_0*((A_1[0]+A_1[2])*(A_2[1]-A_2[3]));

		A_0[0] += (Qp + Qm);

		A_0[1] -= (Pp + Pm);

		A_0[2] += value_type(Pp.imag()-Pm.imag(),Pm.real()-Pp.real());

		A_0[3] += (Qm - Qp);

	    }
	    
	    /* Spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1-g5/2).psi2: */
	    
	    static void ggR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]+A_0[3])*(A_2[0]-A_2[2])+value_type(A_0[1].real() + A_0[2].imag(),A_0[1].imag()-A_0[2].real())*(A_2[1]-A_2[3]));

		value_type Pm = (r_value_type)0.5*C_0*(value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag() + A_0[2].real())*(A_2[0]-A_2[2])+(A_0[0]-A_0[3])*(A_2[1]-A_2[3]));

		A_1[0] += Pp;

		A_1[1] += Pm;

		A_1[2] += Pp;

		A_1[3] += Pm;

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1-g5/2): */

	    static void ggR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]+A_0[3])*(A_1[0]+A_1[2])+value_type(A_0[1].real() - A_0[2].imag(),A_0[1].imag()+A_0[2].real())*(A_1[1]+A_1[3]));

		value_type Pm = (r_value_type)0.5*C_0*(value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag() - A_0[2].real())*(A_1[0]+A_1[2])+(A_0[0]-A_0[3])*(A_1[1]+A_1[3]));

		A_2[0] += Pp;

		A_2[1] += Pm;

		A_2[2] -= Pp;

		A_2[3] -= Pm;

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.g(mu).(1-g5/2).C.psi2: */

	    static void ggR_C_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = (r_value_type)0.5*C_0*((A_1[0]+A_1[2])*(A_2[1]-A_2[3]));

		value_type Qm = (r_value_type)0.5*C_0*((A_1[1]+A_1[3])*(A_2[2]-A_2[0]));

		value_type Pp = (r_value_type)0.5*C_0*((A_1[1]+A_1[3])*(A_2[1]-A_2[3]));

		value_type Pm = (r_value_type)0.5*C_0*((A_1[0]+A_1[2])*(A_2[2]-A_2[0]));

		A_0[0] += (Qp + Qm);

		A_0[1] -= (Pp + Pm);

		A_0[2] += value_type(Pp.imag()-Pm.imag(),Pm.real()-Pp.real());

		A_0[3] += (Qm - Qp);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.(1-g5/2).C.psi2: */

	    static void ggR_C_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]+A_0[3])*(A_2[1]-A_2[3])+value_type(A_0[1].real() + A_0[2].imag(),A_0[1].imag()-A_0[2].real())*(A_2[2]-A_2[0]));

		value_type Pm = (r_value_type)0.5*C_0*(value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag() + A_0[2].real())*(A_2[1]-A_2[3])+(A_0[0]-A_0[3])*(A_2[2]-A_2[0]));

		A_1[0] += Pp;

		A_1[1] += Pm;

		A_1[2] += Pp;

		A_1[3] += Pm;

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.(1-g5/2).C: */

	    static void ggR_C_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]+A_0[3])*(A_1[0]+A_1[2])+value_type(A_0[1].real() - A_0[2].imag(),A_0[1].imag()+A_0[2].real())*(A_1[1]+A_1[3]));

		value_type Pm = (r_value_type)0.5*C_0*(value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag() - A_0[2].real())*(A_1[0]+A_1[2])+(A_0[0]-A_0[3])*(A_1[1]+A_1[3]));

		A_2[0] -= Pm;

		A_2[1] += Pp;

		A_2[2] += Pm;

		A_2[3] -= Pp;

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * c*bar{psi1}.C^*.g(mu).(1-g5/2).psi2: */
	    
	    static void Cc_ggR_first(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Qp = (r_value_type)0.5*C_0*((A_1[1]+A_1[3])*(A_2[0]-A_2[2]));

		value_type Qm = (r_value_type)0.5*C_0*((A_1[0]+A_1[2])*(A_2[3]-A_2[1]));

		value_type Pp = (r_value_type)0.5*C_0*((A_1[0]+A_1[2])*(A_2[2]-A_2[0]));

		value_type Pm = (r_value_type)0.5*C_0*((A_1[1]+A_1[3])*(A_2[1]-A_2[3]));

		A_0[0] += (Qp + Qm);

		A_0[1] -= (Pp + Pm);

		A_0[2] += value_type(Pp.imag()-Pm.imag(),Pm.real()-Pp.real());

		A_0[3] += (Qm - Qp);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * c*Vslash.C^*.(1-g5/2).psi2: */

	    static void Cc_ggR_second(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]+A_0[3])*(A_2[0]-A_2[2])+value_type(A_0[1].real() + A_0[2].imag(),A_0[1].imag()-A_0[2].real())*(A_2[1]-A_2[3]));

		value_type Pm = (r_value_type)0.5*C_0*(value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag() + A_0[2].real())*(A_2[0]-A_2[2])+(A_0[0]-A_0[3])*(A_2[1]-A_2[3]));

		A_1[0] -= Pm;

		A_1[1] += Pp;

		A_1[2] -= Pm;

		A_1[3] += Pp;

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * c*bar{psi1}.Vslash.C^*.(1-g5/2): */

	    static void Cc_ggR_third(const value_type& C_0,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type Pp = (r_value_type)0.5*C_0*((A_0[0]+A_0[3])*(A_1[1]+A_1[3])-value_type(A_0[1].real() - A_0[2].imag(),A_0[1].imag()+A_0[2].real())*(A_1[0]+A_1[2]));

		value_type Pm = (r_value_type)0.5*C_0*(value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag() - A_0[2].real())*(A_1[1]+A_1[3])-(A_0[0]-A_0[3])*(A_1[0]+A_1[2]));

		A_2[0] += Pp;

		A_2[1] += Pm;

		A_2[2] -= Pp;

		A_2[3] -= Pm;

	    }

	    /* Spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0*PL+c1*PR).psi2: */
	    
	    static void ggLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = (r_value_type)0.5*((C_0+C_1)*A_2[0] + (C_0-C_1)*A_2[2]);

		value_type psi20 = (r_value_type)0.5*((C_0-C_1)*A_2[0] + (C_0+C_1)*A_2[2]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_2[1] + (C_0-C_1)*A_2[3]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_2[1] + (C_0+C_1)*A_2[3]);

		A_0[0] += (A_1[0]*psi02 + A_1[1]*psi13 - A_1[2]*psi20 - A_1[3]*psi31);

		A_0[1] += (A_1[0]*psi31 + A_1[1]*psi20 - A_1[2]*psi13 - A_1[3]*psi02);

		A_0[2] += times_i(-A_1[0]*psi31 + A_1[1]*psi20 + A_1[2]*psi13 - A_1[3]*psi02);

		A_0[3] += (A_1[0]*psi20 - A_1[1]*psi31 - A_1[2]*psi02 + A_1[3]*psi13);

	    }

	    /* Spinor vertex recursive relation psi1 +=
	     * Vslash.(c0*PL+c1*PR).psi2: */

	    static void ggLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = (r_value_type)0.5*((C_0+C_1)*A_2[0] + (C_0-C_1)*A_2[2]);

		value_type psi20 = (r_value_type)0.5*((C_0-C_1)*A_2[0] + (C_0+C_1)*A_2[2]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_2[1] + (C_0-C_1)*A_2[3]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_2[1] + (C_0+C_1)*A_2[3]);

		A_1[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi31);

		A_1[1] += (-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_1[2] += (-A_0[0]*psi20 + A_0[3]*psi02 + value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi13);

		A_1[3] += (value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi02 - A_0[0]*psi31 - A_0[3]*psi13);

	    }

	    /* Spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0*PL+c1*PR): */

	    static void ggLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = (r_value_type)0.5*((C_0+C_1)*A_1[0] - (C_0-C_1)*A_1[2]);

		value_type psi20 = (r_value_type)0.5*((C_0-C_1)*A_1[0] - (C_0+C_1)*A_1[2]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_1[1] - (C_0-C_1)*A_1[3]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_1[1] - (C_0+C_1)*A_1[3]);

		A_2[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi31);

		A_2[1] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_2[2] += (A_0[0]*psi20 - A_0[3]*psi02 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi13);

		A_2[3] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi02 + A_0[0]*psi31 + A_0[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.g(mu).(c0*PL+c1*PR).C.psi2: */
	    
	    static void ggLR_C_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = -(r_value_type)0.5*((C_0+C_1)*A_2[3] + (C_0-C_1)*A_2[1]);

		value_type psi20 = -(r_value_type)0.5*((C_0-C_1)*A_2[3] + (C_0+C_1)*A_2[1]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_2[2] + (C_0-C_1)*A_2[0]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_2[2] + (C_0+C_1)*A_2[0]);

		A_0[0] += (A_1[0]*psi02 + A_1[1]*psi13 - A_1[2]*psi20 - A_1[3]*psi31);

		A_0[1] += (A_1[0]*psi31 + A_1[1]*psi20 - A_1[2]*psi13 - A_1[3]*psi02);

		A_0[2] += times_i(-A_1[0]*psi31 + A_1[1]*psi20 + A_1[2]*psi13 - A_1[3]*psi02);

		A_0[3] += (A_1[0]*psi20 - A_1[1]*psi31 - A_1[2]*psi02 + A_1[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * Vslash.(c0*PL+c1*PR).C.psi2: */

	    static void ggLR_C_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = -(r_value_type)0.5*((C_0+C_1)*A_2[3] + (C_0-C_1)*A_2[1]);

		value_type psi20 = -(r_value_type)0.5*((C_0-C_1)*A_2[3] + (C_0+C_1)*A_2[1]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_2[2] + (C_0-C_1)*A_2[0]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_2[2] + (C_0+C_1)*A_2[0]);

		A_1[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi31);

		A_1[1] += (-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_1[2] += (-A_0[0]*psi20 + A_0[3]*psi02 + value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi13);

		A_1[3] += (value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi02 - A_0[0]*psi31 - A_0[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.Vslash.(c0*PL+c1*PR).C: */

	    static void ggLR_C_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = (r_value_type)0.5*((C_0+C_1)*A_1[0] - (C_0-C_1)*A_1[2]);

		value_type psi20 = (r_value_type)0.5*((C_0-C_1)*A_1[0] - (C_0+C_1)*A_1[2]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_1[1] - (C_0-C_1)*A_1[3]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_1[1] - (C_0+C_1)*A_1[3]);

		A_2[3] -= (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi31);

		A_2[2] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_2[1] -= (A_0[0]*psi20 - A_0[3]*psi02 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi13);

		A_2[0] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi02 + A_0[0]*psi31 + A_0[3]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation V(mu) +=
	     * bar{psi1}.C^*.g(mu).(c0*PL+c1*PR).psi2: */
	    
	    static void Cc_ggLR_first(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = (r_value_type)0.5*((C_0+C_1)*A_2[0] + (C_0-C_1)*A_2[2]);

		value_type psi20 = (r_value_type)0.5*((C_0-C_1)*A_2[0] + (C_0+C_1)*A_2[2]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_2[1] + (C_0-C_1)*A_2[3]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_2[1] + (C_0+C_1)*A_2[3]);

		A_0[0] += (A_1[3]*psi02 - A_1[2]*psi13 - A_1[1]*psi20 + A_1[0]*psi31);

		A_0[1] += (A_1[3]*psi31 - A_1[2]*psi20 - A_1[1]*psi13 + A_1[0]*psi02);

		A_0[2] += times_i(-A_1[3]*psi31 - A_1[2]*psi20 + A_1[1]*psi13 + A_1[0]*psi02);

		A_0[3] += (A_1[3]*psi20 + A_1[2]*psi31 - A_1[1]*psi02 - A_1[0]*psi13);

	    }

	    /* Conjugate spinor vertex recursive relation psi1 +=
	     * C^*.Vslash.(c0*PL+c1*PR).psi2: */

	    static void Cc_ggLR_second(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = (r_value_type)0.5*((C_0+C_1)*A_2[0] + (C_0-C_1)*A_2[2]);

		value_type psi20 = (r_value_type)0.5*((C_0-C_1)*A_2[0] + (C_0+C_1)*A_2[2]);

		value_type psi13 = (r_value_type)0.5*((C_0+C_1)*A_2[1] + (C_0-C_1)*A_2[3]);

		value_type psi31 = (r_value_type)0.5*((C_0-C_1)*A_2[1] + (C_0+C_1)*A_2[3]);

		A_1[0] -= (value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi02 - A_0[0]*psi31 - A_0[3]*psi13);

		A_1[1] += (-A_0[0]*psi20 + A_0[3]*psi02 + value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi13);

		A_1[2] -= (-value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_1[3] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi31);

	    }

	    /* Conjugate spinor vertex recursive relation bar{psi2} +=
	     * bar{psi1}.C^*.Vslash.(c0*PL+c1*PR): */

	    static void Cc_ggLR_third(const value_type& C_0,const value_type& C_1,iterator A_0,iterator A_1,iterator A_2)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<4),"tensor iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<4),"tensor iterator 2 out of range");
		
		value_type psi02 = (r_value_type)0.5*((C_0+C_1)*A_1[3] - (C_0-C_1)*A_1[1]);

		value_type psi20 = (r_value_type)0.5*((C_0-C_1)*A_1[3] - (C_0+C_1)*A_1[1]);

		value_type psi13 = (r_value_type)0.5*(-(C_0+C_1)*A_1[2] + (C_0-C_1)*A_1[0]);

		value_type psi31 = (r_value_type)0.5*(-(C_0-C_1)*A_1[2] + (C_0+C_1)*A_1[0]);

		A_2[0] += (A_0[0]*psi02 - A_0[3]*psi20 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi31);

		A_2[1] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi20 + A_0[0]*psi13 + A_0[3]*psi31);

		A_2[2] += (A_0[0]*psi20 - A_0[3]*psi02 - value_type(A_0[1].real()-A_0[2].imag(),A_0[1].imag()+A_0[2].real())*psi13);

		A_2[3] += (-value_type(A_0[1].real()+A_0[2].imag(),A_0[1].imag()-A_0[2].real())*psi02 + A_0[0]*psi31 + A_0[3]*psi13);

	    }

	    /* Massless Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion subamplitude to be propagated and p the momentum. */

	    static void massless_prop(const value_type& C_0,iterator A_0,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];

		temp[0] = C_0*(p[0]*A_0[0] - p[3]*A_0[2] - std::conj(Qp)*A_0[3]);

		temp[1] = C_0*(p[0]*A_0[1] - Qp*A_0[2] + p[3]*A_0[3]);

		temp[2] = C_0*(p[3]*A_0[0] + std::conj(Qp)*A_0[1] - p[0]*A_0[2]);

		temp[3] = C_0*(Qp*A_0[0] - p[3]*A_0[1] - p[0]*A_0[3]);

		A_0[0] = temp[0];

		A_0[1] = temp[1];

		A_0[2] = temp[2];

		A_0[3] = temp[3];

	    }

	    /* Massless Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion-multiplet subamplitude to be propagated and p the momentum. */

	    static void massless_prop(const value_type& C_0,iterator A_0,iterator A_1,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];
		do
		{
		    temp[0] = C_0*(p[0]*A_0[0] - p[3]*A_0[2] - std::conj(Qp)*A_0[3]);

		    temp[1] = C_0*(p[0]*A_0[1] - Qp*A_0[2] + p[3]*A_0[3]);

		    temp[2] = C_0*(p[3]*A_0[0] + std::conj(Qp)*A_0[1] - p[0]*A_0[2]);

		    temp[3] = C_0*(Qp*A_0[0] - p[3]*A_0[1] - p[0]*A_0[3]);

		    A_0[0] = temp[0];

		    A_0[1] = temp[1];

		    A_0[2] = temp[2];

		    A_0[3] = temp[3];

		    A_0+=4;
		}
		while(A_0 != A_1);

	    }

	    /* Massive Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion subamplitude to be propagated, p the momentum and
	     * m the (complex) mass. */

	    static void massive_prop(const value_type& C_0,const value_type& m,iterator A_0,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];

		temp[0] = C_0*((p[0]+m)*A_0[0] - p[3]*A_0[2] - std::conj(Qp)*A_0[3]);

		temp[1] = C_0*((p[0]+m)*A_0[1] - Qp*A_0[2] + p[3]*A_0[3]);

		temp[2] = C_0*(p[3]*A_0[0] + std::conj(Qp)*A_0[1] - (p[0]-m)*A_0[2]);

		temp[3] = C_0*(Qp*A_0[0] - p[3]*A_0[1] - (p[0]-m)*A_0[3]);

		A_0[0] = temp[0];

		A_0[1] = temp[1];

		A_0[2] = temp[2];

		A_0[3] = temp[3];

	    }

	    /* Massive Dirac propagator. C_0 is 1 divided by the denominator,
	     * A_0 the fermion multiplet subamplitude to be propagated, p the
	     * momentum and m the (complex) mass. */

	    static void massive_prop(const value_type& C_0,const value_type& m,iterator A_0,iterator A_1,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];
		do
		{
		    temp[0] = C_0*((p[0]+m)*A_0[0] - p[3]*A_0[2] - std::conj(Qp)*A_0[3]);

		    temp[1] = C_0*((p[0]+m)*A_0[1] - Qp*A_0[2] + p[3]*A_0[3]);

		    temp[2] = C_0*(p[3]*A_0[0] + std::conj(Qp)*A_0[1] - (p[0]-m)*A_0[2]);

		    temp[3] = C_0*(Qp*A_0[0] - p[3]*A_0[1] - (p[0]-m)*A_0[3]);

		    A_0[0] = temp[0];

		    A_0[1] = temp[1];

		    A_0[2] = temp[2];

		    A_0[3] = temp[3];

		    A_0+=4;
		}
		while(A_0 != A_1);

	    }

	    /* Massless Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion subamplitude to be propagated
	     * and p the momentum. */

	    static void massless_anti_prop(const value_type& C_0,iterator A_0,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];

		temp[0] = C_0*(-p[0]*A_0[0] - p[3]*A_0[2] - Qp*A_0[3]);

		temp[1] = C_0*(-p[0]*A_0[1] - std::conj(Qp)*A_0[2] + p[3]*A_0[3]);

		temp[2] = C_0*(p[3]*A_0[0] + Qp*A_0[1] + p[0]*A_0[2]);

		temp[3] = C_0*(std::conj(Qp)*A_0[0] - p[3]*A_0[1] + p[0]*A_0[3]);

		A_0[0] = temp[0];

		A_0[1] = temp[1];

		A_0[2] = temp[2];

		A_0[3] = temp[3];

	    }

	    /* Massless Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated and p the momentum. */

	    static void massless_anti_prop(const value_type& C_0,iterator A_0,iterator A_1,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];
		do
		{
		    temp[0] = C_0*(-p[0]*A_0[0] - p[3]*A_0[2] - Qp*A_0[3]);

		    temp[1] = C_0*(-p[0]*A_0[1] - std::conj(Qp)*A_0[2] + p[3]*A_0[3]);

		    temp[2] = C_0*(p[3]*A_0[0] + Qp*A_0[1] + p[0]*A_0[2]);

		    temp[3] = C_0*(std::conj(Qp)*A_0[0] - p[3]*A_0[1] + p[0]*A_0[3]);

		    A_0[0] = temp[0];

		    A_0[1] = temp[1];

		    A_0[2] = temp[2];

		    A_0[3] = temp[3];

		    A_0+=4;
		}
		while(A_0 != A_1);

	    }

	    /* Massive Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated, p the momentum and m the (complex) mass. */

	    static void massive_anti_prop(const value_type& C_0,const value_type& m,iterator A_0,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];

		temp[0] = C_0*((-p[0]+m)*A_0[0] - p[3]*A_0[2] - Qp*A_0[3]);

		temp[1] = C_0*((-p[0]+m)*A_0[1] - std::conj(Qp)*A_0[2] + p[3]*A_0[3]);

		temp[2] = C_0*(p[3]*A_0[0] + Qp*A_0[1] + (p[0]+m)*A_0[2]);

		temp[3] = C_0*(std::conj(Qp)*A_0[0] - p[3]*A_0[1] + (p[0]+m)*A_0[3]);

		A_0[0] = temp[0];

		A_0[1] = temp[1];

		A_0[2] = temp[2];

		A_0[3] = temp[3];

	    }

	    /* Massive Dirac antiparticle propagator. C_0 is 1 divided by the
	     * denominator, A_0 the anti-fermion multiplet subamplitude to be
	     * propagated, p the momentum and m the (complex) mass. */

	    static void massive_anti_prop(const value_type& C_0,const value_type& m,iterator A_0,iterator A_1,const vector<r_value_type,4>& p)
	    {
		CAMGEN_ERROR_IF((A_0.range()<4),"tensor iterator 0 out of range");
		
		value_type Qp(p[1],p[2]);
		value_type temp[4];
		do
		{
		    temp[0] = C_0*((-p[0]+m)*A_0[0] - p[3]*A_0[2] - Qp*A_0[3]);

		    temp[1] = C_0*((-p[0]+m)*A_0[1] - std::conj(Qp)*A_0[2] + p[3]*A_0[3]);

		    temp[2] = C_0*(p[3]*A_0[0] + Qp*A_0[1] + (p[0]+m)*A_0[2]);

		    temp[3] = C_0*(std::conj(Qp)*A_0[0] - p[3]*A_0[1] + (p[0]+m)*A_0[3]);

		    A_0[0] = temp[0];

		    A_0[1] = temp[1];

		    A_0[2] = temp[2];

		    A_0[3] = temp[3];

		    A_0+=4;
		}
		while(A_0 != A_1);

	    }
    };

    template<class value_t>const int Pauli_basis::implementation<value_t,4,true>::eta_g;
    template<class value_t>const int Pauli_basis::implementation<value_t,4,true>::eta_g5;
    template<class value_t>const int Pauli_basis::implementation<value_t,4,true>::eta_gg5;
    template<class value_t>const int Pauli_basis::implementation<value_t,4,true>::eta_gg;
}

#ifdef CAMGEN_SPINOR_FAC_H_
#include <Camgen/spinfacPb.h>
#endif /*CAMGEN_SPINOR_FAC_H_*/

#ifdef CAMGEN_M_SPINOR_FAC_H_

#ifdef CAMGEN_KS_TYPE_H_
#include <Camgen/KSspinPb.h>
#endif /*CAMGEN_KS_TYPE_H_*/

#ifdef CAMGEN_HEL_TYPE_H_
#include <Camgen/helspinPb.h>
#endif /*CAMGEN_HEL_TYPE_H_*/

#ifdef CAMGEN_POL_TYPE_H_
#include <Camgen/polspinPb.h>
#endif /*CAMGEN_POL_TYPE_H_*/

#endif /*CAMGEN_M_SPINOR_FAC_H_*/

#endif /*CAMGEN_PAULI_BASIS_H_*/

