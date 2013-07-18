//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/license_print.h>
#include <Camgen/stdrand.h>
#include <Camgen/rn_strm.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/KS_type.h>
#include <Camgen/hel_type.h>
#include <Camgen/pol_type.h>
#include <Camgen/m_spinor_fac.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for the massless and massive fermion/vector wave functions   *
 * in Camgen. In particular the template specialisations for the Pauli basis of *
 * the Dirac algebra are tested, in combination with the KS_type, helicity_type  *
 * or polarised spin vectors for massive fermions.                               *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Definition of a clone of the Pauli basis class, which does not contain all the
 * overloaded contractions and will not by instantiated by the spinor factory
 * specialisations: */

namespace Camgen
{
    class Pauli_basis_clone
    {
	public:

	    /* Type definition of the corresponding spacetime-signature: */

	    typedef Minkowski_type spacetime_type;

	    /* Implementation metafunction of the matrix basis: by default, this
	     * class is empty: */

	    template<class value_t,std::size_t dim,bool q=(dim%2==0)>class implementation: public Dirac_algebra<value_t,implementation<value_t,dim>,dim,q>{};
    };

    /* Implementation metafunction, only filling up the Dirac matrices: */

    template<class value_t>class Pauli_basis_clone::implementation<value_t,4,true>: public Dirac_algebra<value_t,implementation<value_t,4,true>,4,true>
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
	    
	    const static int eta_g=-1;
	    const static int eta_g5=1;
	    const static int eta_gg5=1;
	    const static int eta_gg=-1;

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
    };
}

/* Energy scale for testing momenta: */

#define ESCALE 75

using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing spinor/vector wave functions in the Pauli basis..."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    /* Useful type definitions: */

    typedef double r_value_type;
    typedef std::complex<r_value_type> value_type;
    typedef Pauli_basis::implementation<r_value_type,4> Dirac_algebra_type;
    typedef Pauli_basis_clone::implementation<r_value_type,4> Dirac_algebra_clone_type;
    typedef Minkowski_type::implementation<r_value_type,4> spacetime_type;
    typedef vector<r_value_type,4> momentum_type;
    typedef tensor<value_type> tensor_type;
    typedef random_number_stream<r_value_type,std::random> generator_type;

    /* Initialisation phase: */

    massless_spinor_factory<Dirac_algebra_type,3,4>::initialise();
    massless_spinor_factory<Dirac_algebra_clone_type,3,4>::initialise();

    /* Declarations of testing spinors, vectors, numbers and momenta: */

    tensor_type col_spinor(1,4);
    tensor_type row_spinor(1,4);
    tensor_type hplus_vector(1,4);
    tensor_type hlong_vector(1,4);
    tensor_type hmin_vector(1,4);
    tensor_type temp_spinor(1,4);
    tensor_type temp_vector(1,4);
    tensor_type check;
    momentum_type p,q,s;
    r_value_type m=generator_type::throw_number(ESCALE);
    value_type factor(generator_type::throw_number(),generator_type::throw_number());
    p[1]=generator_type::throw_number(ESCALE);
    p[2]=generator_type::throw_number(ESCALE);
    p[3]=generator_type::throw_number(ESCALE);
    p[0]=std::sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
    q[1]=p[1];
    q[2]=p[2];
    q[3]=p[3];
    q[0]=std::sqrt(p[0]*p[0]+m*m);
    tensor_type pslash(2,4,4);
    tensor_type pslash_L(2,4,4);
    tensor_type pslash_R(2,4,4);
    tensor_type sslash_L(2,4,4);
    tensor_type sslash_R(2,4,4);
    tensor_type uubar(2,4,4);
    for(unsigned mu=0;mu<4;++mu)
    {
	for(unsigned i=0;i<4;++i)
	{
	    for(unsigned j=0;j<4;++j)
	    {
		pslash(i,j)+=p[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::gamma(mu,i,j);
		pslash_L(i,j)+=p[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_L(mu,i,j);
		pslash_R(i,j)+=p[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_R(mu,i,j);
	    }
	}
    }

#define FACTORY massless_spinor_factory<Dirac_algebra_type,3,4>
#define FACTORY_CHECKER massless_spinor_factory<Dirac_algebra_clone_type,3,4>

    std::cout<<"Checking u+ spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(factor,col_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_p(factor,check.begin(),&p,NULL);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(col_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_p(check.begin(),&p,NULL);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether gL.u+=u+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
		check(i)+=Dirac_algebra_type::gamma_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether gR.u+=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=Dirac_algebra_type::gamma_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether pslash.u+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking ubar+ spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(factor,row_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_p_bar(factor,check.begin(),&p,NULL);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(row_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_p_bar(check.begin(),&p,NULL);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether (ubar+)gR=ubar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
		check(i)+=row_spinor(j)*(Dirac_algebra_type::gamma_R(j,i));
	}
    }
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    std::cout<<"Checking whether (ubar+).gL=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*(Dirac_algebra_type::gamma_L(j,i));
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar+).pslash=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u+)(ubar+)=PR.pslash..........";
    std::cout.flush();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor(i)*row_spinor(j);
	}
    }
    if(!equal_sequences(uubar,pslash_R))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking u- spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(factor,col_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_m(factor,check.begin(),&p,NULL);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(col_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_m(check.begin(),&p,NULL);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether gR.u-=u-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
		check(i)+=Dirac_algebra_type::gamma_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether gL.u-=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=Dirac_algebra_type::gamma_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether pslash.u-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking ubar- spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(factor,row_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_m_bar(factor,check.begin(),&p,NULL);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(row_spinor.begin(),&p,NULL);
    FACTORY_CHECKER::make_u_m_bar(check.begin(),&p,NULL);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether (ubar-)gL=ubar-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
		check(i)+=row_spinor(j)*(Dirac_algebra_type::gamma_L(j,i));
	}
    }
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether (ubar-).gR=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*(Dirac_algebra_type::gamma_R(j,i));
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar-).pslash=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u-)(ubar-)=PL.pslash..........";
    std::cout.flush();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor(i)*row_spinor(j);
	}
    }
    if(!equal_sequences(uubar,pslash_L))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massless e+ vector construction..........";
    std::cout.flush();
    hplus_vector.reset();
    check=hplus_vector;
    FACTORY::make_e_p(factor,hplus_vector.begin(),&p,NULL);
    massless_spinor_factory<Dirac_algebra_clone_type,3,4>::make_e_p(factor,check.begin(),&p,NULL);
    if(!equal_sequences(hplus_vector,check))
    {
	return 1;
    }
    hplus_vector.reset();
    check=hplus_vector;
    FACTORY::make_e_p(hplus_vector.begin(),&p,NULL);
    FACTORY_CHECKER::make_e_p(check.begin(),&p,NULL);
    if(!equal_sequences(hplus_vector,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massless e+ is massless..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hplus_vector,hplus_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massless p.e+ is zero..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(p,hplus_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking massless e- vector construction..........";
    std::cout.flush();
    hmin_vector.reset();
    check=hmin_vector;
    FACTORY::make_e_m(factor,hmin_vector.begin(),&p,NULL);
    FACTORY_CHECKER::make_e_m(factor,check.begin(),&p,NULL);
    if(!equal_sequences(hmin_vector,check))
    {
	return 1;
    }
    hmin_vector.reset();
    check=hmin_vector;
    FACTORY::make_e_m(hmin_vector.begin(),&p,NULL);
    FACTORY_CHECKER::make_e_m(check.begin(),&p,NULL);
    if(!equal_sequences(hmin_vector,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massless e- is massless..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hmin_vector,hmin_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massless p.e- is zero..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(p,hmin_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massless e- = e+*..........";
    std::cout.flush();
    if(!equal_sequences(hmin_vector,std::conj(hplus_vector)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massless e+.e- = -1..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hmin_vector,hplus_vector),value_type(-1,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking massive e+ vector construction..........";
    std::cout.flush();
    hplus_vector.reset();
    check=hplus_vector;
    FACTORY::make_m_e_p(factor,hplus_vector.begin(),&q,&m);
    FACTORY_CHECKER::make_m_e_p(factor,check.begin(),&q,&m);
    if(!equal_sequences(hplus_vector,check))
    {
	return 1;
    }
    hplus_vector.reset();
    check=hplus_vector;
    FACTORY::make_m_e_p(hplus_vector.begin(),&q,&m);
    FACTORY_CHECKER::make_m_e_p(check.begin(),&q,&m);
    if(!equal_sequences(hplus_vector,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive e+ is massless..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hplus_vector,hplus_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive p.e+ is zero..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(q,hplus_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking massive e- vector construction..........";
    std::cout.flush();
    hmin_vector.reset();
    check=hmin_vector;
    FACTORY::make_m_e_m(factor,hmin_vector.begin(),&q,&m);
    FACTORY_CHECKER::make_m_e_m(factor,check.begin(),&q,&m);
    if(!equal_sequences(hmin_vector,check))
    {
	return 1;
    }
    hmin_vector.reset();
    check=hmin_vector;
    FACTORY::make_m_e_m(hmin_vector.begin(),&q,&m);
    FACTORY_CHECKER::make_m_e_m(check.begin(),&q,&m);
    if(!equal_sequences(hmin_vector,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive e- is massless..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hmin_vector,hmin_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive p.e- is zero..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(q,hmin_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive e- = e+*..........";
    std::cout.flush();
    if(!equal_sequences(hmin_vector,std::conj(hplus_vector)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive e+.e- = -1..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hmin_vector,hplus_vector),value_type(-1,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive e0 vector construction..........";
    std::cout.flush();
    hlong_vector.reset();
    check=hlong_vector;
    FACTORY::make_m_e_L(factor,hlong_vector.begin(),&q,&m);
    FACTORY_CHECKER::make_m_e_L(factor,check.begin(),&q,&m);
    if(!equal_sequences(hlong_vector,check))
    {
	return 1;
    }
    hlong_vector.reset();
    check=hlong_vector;
    FACTORY::make_m_e_L(hlong_vector.begin(),&q,&m);
    FACTORY_CHECKER::make_m_e_L(check.begin(),&q,&m);
    if(!equal_sequences(hlong_vector,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive e0.e0 = -1..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hlong_vector,hlong_vector),value_type(-1,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive p.e0 is zero..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(q,hlong_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive e0.e+ = 0..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hlong_vector,hplus_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether massive e0.e- = 0..........";
    std::cout.flush();
    if(!equals(spacetime_type::dot(hlong_vector,hmin_vector),value_type(0,0)))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether u- = C.(ubar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_m(col_spinor.begin(),&p,NULL);
    FACTORY::make_u_p_bar(row_spinor.begin(),&p,NULL);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether u- = C.(ubar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_m(col_spinor.begin(),&p,NULL);
    FACTORY::make_u_p_bar(row_spinor.begin(),&p,NULL);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    /* Filling the massive momentum slash matrix: */

    pslash.reset();
    for(unsigned mu=0;mu<4;++mu)
    {
	for(unsigned i=0;i<4;++i)
	{
	    for(unsigned j=0;j<4;++j)
	    {
		pslash(i,j)+=q[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::gamma(mu,i,j);
	    }
	}
    }
#define SPINVEC_T helicity_type
#undef FACTORY
#define FACTORY massive_spinor_factory<Dirac_algebra_type,SPINVEC_T,3,4>
#undef FACTORY_CHECKER
#define FACTORY_CHECKER massive_spinor_factory<Dirac_algebra_clone_type,SPINVEC_T,3,4>
#define SPINVEC_PRINT "helicity"

    s.assign(0);
    SPINVEC_T::implementation<spacetime_type,3>::make_spinvec(s,&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    for(unsigned mu=0;mu<4;++mu)
	    {
		sslash_L(i,j)-=(r_value_type)0.5*s[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_5(mu,i,j);
		sslash_R(i,j)+=(r_value_type)0.5*s[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_5(mu,i,j);
	    }
	}
	sslash_L(i,i)+=(r_value_type)0.5;
	sslash_R(i,i)+=(r_value_type)0.5;
    }

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" u+ spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1+g5.$)u+=u+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)u+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash-m).u+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" ubar+ spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	std::cout<<"verkeerd 1: "<<row_spinor<<" neq "<<check<<std::endl;
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	std::cout<<"verkeerd 2: "<<row_spinor<<" neq "<<check<<std::endl;
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(ubar+)(1+g5.$)=ubar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(ubar+)(1-g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar+)(pslash-m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u+)(ubar+)=1/2(1+g5.sslash)(pslash+m)..........";
    std::cout.flush();
    uubar.reset();
    tensor_type rhs(2,4,4);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=m*sslash_L(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_L(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" u- spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash-m).u-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" ubar- spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1-g5.$)u-=u-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1+g5.$)u-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar-)(pslash-m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u-)(ubar-)=1/2(1-g5.sslash)(pslash+m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=m*sslash_R(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_R(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" v+ spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_p(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_p(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(ubar-)(1-g5.$)=ubar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(ubar-)(1+g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1+g5.$)v+=v+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)v+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash+m).v+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,-m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" vbar+ spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_p_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_p_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar+)(1+g5.$)=vbar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar+)(1-g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (vbar+)(pslash+m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,-m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (v+)(vbar+)=1/2(1+g5.sslash)(pslash-m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=-m*sslash_L(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_L(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" v- spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_m(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_m(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)v-=v-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1+g5.$)v-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash+m).v-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,-m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" vbar- spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_m_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_m_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(vbar-)(1-g5.$)=vbar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar-)(1+g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (vbar-)(pslash+m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,-m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (v-)(vbar-)=1/2(1-g5.sslash)(pslash-m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=-m*sslash_R(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_R(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;

    std::cout<<"Checking whether u+ = C.(vbar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_p(col_spinor.begin(),&q,&m);
    FACTORY::make_v_p_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether u- = C.(vbar-)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_m(col_spinor.begin(),&q,&m);
    FACTORY::make_v_m_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether v+ = C.(ubar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_v_p(col_spinor.begin(),&q,&m);
    FACTORY::make_u_p_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether v- = C.(ubar-)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_v_m(col_spinor.begin(),&q,&m);
    FACTORY::make_u_m_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

#undef SPINVEC_T
#define SPINVEC_T KS_type
#undef FACTORY
#define FACTORY massive_spinor_factory<Dirac_algebra_type,SPINVEC_T,3,4>
#undef FACTORY_CHECKER
#define FACTORY_CHECKER massive_spinor_factory<Dirac_algebra_clone_type,SPINVEC_T,3,4>
#undef SPINVEC_PRINT
#define SPINVEC_PRINT "KS_type"

    sslash_L.reset();
    sslash_R.reset();
    s.assign(0);
    SPINVEC_T::implementation<spacetime_type,3>::make_spinvec(s,&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    for(unsigned mu=0;mu<4;++mu)
	    {
		sslash_L(i,j)-=(r_value_type)0.5*s[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_5(mu,i,j);
		sslash_R(i,j)+=(r_value_type)0.5*s[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_5(mu,i,j);
	    }
	}
	sslash_L(i,i)+=(r_value_type)0.5;
	sslash_R(i,i)+=(r_value_type)0.5;
    }

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" u+ spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1+g5.$)u+=u+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)u+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash-m).u+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" ubar+ spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(ubar+)(1+g5.$)=ubar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(ubar+)(1-g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar+)(pslash-m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u+)(ubar+)=1/2(1+g5.sslash)(pslash+m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=m*sslash_L(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_L(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" u- spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash-m).u-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" ubar- spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1-g5.$)u-=u-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1+g5.$)u-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar-)(pslash-m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u-)(ubar-)=1/2(1-g5.sslash)(pslash+m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=m*sslash_R(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_R(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" v+ spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_p(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_p(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(ubar-)(1-g5.$)=ubar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(ubar-)(1+g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1+g5.$)v+=v+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)v+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash+m).v+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,-m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" vbar+ spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_p_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_p_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar+)(1+g5.$)=vbar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar+)(1-g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (vbar+)(pslash+m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,-m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (v+)(vbar+)=1/2(1+g5.sslash)(pslash-m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=-m*sslash_L(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_L(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" v- spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_m(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_m(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)v-=v-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1+g5.$)v-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash+m).v-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,-m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" vbar- spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_m_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_m_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(vbar-)(1-g5.$)=vbar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar-)(1+g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (vbar-)(pslash+m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,-m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (v-)(vbar-)=1/2(1-g5.sslash)(pslash-m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=-m*sslash_R(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_R(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking whether u+ = C.(vbar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_p(col_spinor.begin(),&q,&m);
    FACTORY::make_v_p_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether u- = C.(vbar-)^T..........";
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_m(col_spinor.begin(),&q,&m);
    FACTORY::make_v_m_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether v+ = C.(ubar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_v_p(col_spinor.begin(),&q,&m);
    FACTORY::make_u_p_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether v- = C.(ubar-)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_v_m(col_spinor.begin(),&q,&m);
    FACTORY::make_u_m_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

#undef SPINVEC_T
#define SPINVEC_T polarised_type
#undef FACTORY
#define FACTORY massive_spinor_factory<Dirac_algebra_type,SPINVEC_T,3,4>
#undef FACTORY_CHECKER
#define FACTORY_CHECKER massive_spinor_factory<Dirac_algebra_clone_type,SPINVEC_T,3,4>
#undef SPINVEC_PRINT
#define SPINVEC_PRINT "z-polarised"

    sslash_L.reset();
    sslash_R.reset();
    s.assign(0);
    SPINVEC_T::implementation<spacetime_type,3>::make_spinvec(s,&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    for(unsigned mu=0;mu<4;++mu)
	    {
		sslash_L(i,j)-=(r_value_type)0.5*s[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_5(mu,i,j);
		sslash_R(i,j)+=(r_value_type)0.5*s[mu]*spacetime_type::metric(mu,mu)*Dirac_algebra_type::ggamma_5(mu,i,j);
	    }
	}
	sslash_L(i,i)+=(r_value_type)0.5;
	sslash_R(i,i)+=(r_value_type)0.5;
    }

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" u+ spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_p(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1+g5.$)u+=u+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)u+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash-m).u+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" ubar+ spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_p_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_p_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(ubar+)(1+g5.$)=ubar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(ubar+)(1-g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar+)(pslash-m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u+)(ubar+)=1/2(1+g5.sslash)(pslash+m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=m*sslash_L(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_L(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" u- spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_u_m(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash-m).u-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" ubar- spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_u_m_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_u_m_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1-g5.$)u-=u-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1+g5.$)u-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (ubar-)(pslash-m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (u-)(ubar-)=1/2(1-g5.sslash)(pslash+m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=m*sslash_R(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_R(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" v+ spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_p(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_p(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(ubar-)(1-g5.$)=ubar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(ubar-)(1+g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(1+g5.$)v+=v+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)v+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash+m).v+=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,-m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" vbar+ spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_p_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_p_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_p_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar+)(1+g5.$)=vbar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar+)(1-g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (vbar+)(pslash+m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,-m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (v+)(vbar+)=1/2(1+g5.sslash)(pslash-m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=-m*sslash_L(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_L(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" v- spinor construction..........";
    std::cout.flush();
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_m(factor,col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m(factor,check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    col_spinor.reset();
    check=col_spinor;
    FACTORY::make_v_m(col_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m(check.begin(),&q,&m);
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1-g5.$)v-=v-..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_R(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(1+g5.$)v-=0..........";
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=sslash_L(i,j)*col_spinor(j);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (pslash+m).v-=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=pslash(i,j)*col_spinor[j];
	}
    }
    if(!equal_sequences(temp_spinor,-m*col_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking massive "<<SPINVEC_PRINT<<" vbar- spinor construction..........";
    std::cout.flush();
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_m_bar(factor,row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m_bar(factor,check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    row_spinor.reset();
    check=row_spinor;
    FACTORY::make_v_m_bar(row_spinor.begin(),&q,&m);
    FACTORY_CHECKER::make_v_m_bar(check.begin(),&q,&m);
    if(!equal_sequences(row_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 1/2(vbar-)(1-g5.$)=vbar+..........";
    std::cout.flush();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_R(j,i);
	}
    }
    if(!equal_sequences(check,row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 1/2(vbar-)(1+g5.$)=0..........";
    std::cout.flush();
    temp_spinor.reset();
    check.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check(i)+=row_spinor(j)*sslash_L(j,i);
	}
    }
    if(!equal_sequences(check,temp_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (vbar-)(pslash+m)=0..........";
    std::cout.flush();
    check.reset();
    temp_spinor.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    temp_spinor[i]+=row_spinor[j]*pslash(j,i);
	}
    }
    if(!equal_sequences(temp_spinor,-m*row_spinor))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether (v-)(vbar-)=1/2(1-g5.sslash)(pslash-m)..........";
    std::cout.flush();
    uubar.reset();
    rhs.reset();
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    uubar(i,j)=col_spinor[i]*row_spinor[j];
	    rhs(i,j)=-m*sslash_R(i,j);
	    for(unsigned k=0;k<4;++k)
	    {
		rhs(i,j)+=sslash_R(i,k)*pslash(k,j);
	    }
	}
    }
    if(!equal_sequences(uubar,rhs))
    {
	return 1;
    }
    std::cout<<"..........done"<<std::endl;
    
    std::cout<<"Checking whether u+ = C.(vbar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_p(col_spinor.begin(),&q,&m);
    FACTORY::make_v_p_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether u- = C.(vbar-)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_u_m(col_spinor.begin(),&q,&m);
    FACTORY::make_v_m_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether v+ = C.(ubar+)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_v_p(col_spinor.begin(),&q,&m);
    FACTORY::make_u_p_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether v- = C.(ubar-)^T..........";
    std::cout.flush();
    check.reset();
    col_spinor.reset();
    row_spinor.reset();
    FACTORY::make_v_m(col_spinor.begin(),&q,&m);
    FACTORY::make_u_m_bar(row_spinor.begin(),&q,&m);
    for(unsigned i=0;i<4;++i)
    {
	for(unsigned j=0;j<4;++j)
	{
	    check[i]+=(Dirac_algebra_type::C(i,j)*row_spinor[j]);
	}
    }
    if(!equal_sequences(col_spinor,check))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    return 0;
}











