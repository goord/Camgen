//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/license_print.h>
#include <Camgen/stdrand.h>
#include <Camgen/rn_strm.h>
#include <Camgen/Dirac_alg.h>
#include <Camgen/Pauli_basis.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Testing facility for the Pauli basis implementation of the 4-dimensional      *
 * Dirac algebra. The overridden recursive relations in the Pauli basis class    *
 * are tested against the default full contractions in the base class by acting  *
 * on random spinors, vectors and scalars.                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Testing program: */

int main()
{
    Camgen::license_print::disable();
    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing gamma matrix algebra in the Pauli basis..........."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    /* Useful type definitions: */

    typedef double value_type;
    typedef Camgen::Pauli_basis::implementation<value_type,4> Dirac_algebra_type;
    typedef Camgen::random_number_stream<value_type,std::random> generator_type;

    /* Dirac algebra initialisation: */

    Dirac_algebra_type::initialise();

    std::cerr<<"Checking Dirac algebra...........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_anticommutators())
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gamma 5...........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gamma_5())
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking Dirac conjugation matrix...........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_D_matrix())
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking charge conjugation matrix...........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_C_matrix())
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;

    /* Creating a random number generator object: */

    generator_type gen;
    
    std::cerr<<"Checking Id overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Id(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C* overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking g5 overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_g5(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking g5.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_g5_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.g5 overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_g5(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gL overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gL(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gL.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gL_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.gL overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_gL(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gR.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gR_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.gR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_gR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gVA overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gVA(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gVA.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gVA_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.gVA overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_gVA(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gLR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gLR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gLR.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gLR_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.gLR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_gLR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;

    std::cerr<<"Checking g overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_g(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking g.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_g_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.g overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_g(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gg5 overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gg5(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking gg5.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_gg5_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.gg5 overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_gg5(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggL overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggL(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggL.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggL_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.ggL overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_ggL(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggR.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggR_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.ggR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_ggR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggVA overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggVA(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggVA.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggVA_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.ggVA overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_ggVA(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggLR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggLR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking ggLR.C overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_ggLR_C(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;
    
    std::cerr<<"Checking C*.ggLR overridden recursive relations..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_Cc_ggLR(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;

    std::cerr<<"Checking overridden massless propagator..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_massless_prop(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;

    std::cerr<<"Checking overridden massive propagator..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_massive_prop(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;

    std::cerr<<"Checking overridden massless antiparticle propagator..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_massless_anti_prop(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;

    std::cerr<<"Checking overridden massive antiparticle propagator..........";
    std::cerr.flush();
    if(!Dirac_algebra_type::check_massive_anti_prop(gen))
    {
	return 1;
    }
    std::cerr<<"..........done."<<std::endl;

    return 0;
}
