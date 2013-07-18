//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SFFLR_H_
#define CAMGEN_SFFLR_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the left-handed and right_handed linear *
 * combination scalar-fermion-fermion vertex.                            *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition of the Feynman rule class: */

    template<class model_t>class sffLR
    {
	public:

	    /* Reference type definition to the model type: */

	    typedef model_t model_type;
	    
	    /* Dirac algebra type definition: */
	    
	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	private:

	    /* Definition of the spinor size: */

	    static const std::size_t spinor_size=Dirac_algebra_type::base_type::index_range;

	public:

	    /* Vertex rank definition: */

	    static const std::size_t rank=3;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=2;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=spinor_size*spinor_size;
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=false;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=true;

	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];
	    
	    /* Function returning the Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t sffLR<model_t>::spinor_size;
    template<class model_t>const std::size_t sffLR<model_t>::rank;
    template<class model_t>const std::size_t sffLR<model_t>::params;
    template<class model_t>const std::size_t sffLR<model_t>::tensor_size;
    template<class model_t>const bool sffLR<model_t>::p_dependent;
    template<class model_t>const bool sffLR<model_t>::fermionic;
    template<class model_t>const std::size_t sffLR<model_t>::sizes[4]={1,sffLR<model_t>::spinor_size,sffLR<model_t>::spinor_size,0};
    template<class model_t>const std::string sffLR<model_t>::formula="(c1*g_L+c2*g_R)(i,j)";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate<sffLR<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef sffLR<model_t> vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Dirac algebra type definition: */

	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	    /* Definition of the spinor size: */

	    static const std::size_t spinor_size=Dirac_algebra_type::base_type::index_range;

	    /* Sign of the chiral part of the vertex under fermion flow reversal: */

	    static const int sign=Dirac_algebra_type::eta_g5;
	
	public:	
	    
	    /* Initialisation phase: */
	    
	    static void initialise()
	    {
		Dirac_algebra_type::initialise();
	    }

	    /* function returning the index ranges of the interacting
	     * subamplitudes: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n==1 or n==2)
		{
		    v.push_back(spinor_size);
		}
		return v;
	    }

	    /* Recursive relations, calling the gLR-recursive relations in the
	     * Dirac algebra class: */

	    static void first(ARG_LIST)
	    {
		Dirac_algebra_type::gLR_first(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void second(ARG_LIST)
	    {
		Dirac_algebra_type::gLR_second(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void third(ARG_LIST)
	    {
		Dirac_algebra_type::gLR_third(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void fourth(ARG_LIST){}

	    /* Right charge conjugate recursive relations, calling the
	     * gLR_C-recursive relations in the Dirac algebra class: */

	    static void first_C(ARG_LIST)
	    {
		Dirac_algebra_type::gLR_C_first(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void second_C(ARG_LIST)
	    {
		Dirac_algebra_type::gLR_C_second(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void third_C(ARG_LIST)
	    {
		Dirac_algebra_type::gLR_C_third(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void fourth_C(ARG_LIST){}

	    /* Left charge conjugate recursive relations, calling the
	     * Cc_gLR-recursive relations in the Dirac algebra class: */

	    static void Cc_first(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_gLR_first(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void Cc_second(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_gLR_second(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void Cc_third(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_gLR_third(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void Cc_fourth(ARG_LIST){}

	    /* Reversed recursive relations, swapping the couplings if the
	     * reversal sign is -1: */

	    static void reversed_first(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::gLR_first(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::gLR_first(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void reversed_second(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::gLR_third(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::gLR_third(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void reversed_third(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::gLR_second(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::gLR_second(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void reversed_fourth(ARG_LIST){}

	    /* Right charge-conjugate reversed recursive relations, swapping the
	     * couplings and applying gLR_C if the reversal sign is -1: */

	    static void reversed_first_C(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::gLR_C_first(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::gLR_C_first(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void reversed_second_C(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::gLR_C_third(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::gLR_C_third(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void reversed_third_C(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::gLR_C_second(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::gLR_C_second(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void reversed_fourth_C(ARG_LIST){}
	    
	    /* Left charge-conjugate reversed recursive relations, swapping the
	     * couplings and applying Cc_gLR if the reversal sign is -1: */
	    
	    static void Cc_reversed_first(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::Cc_gLR_first(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::Cc_gLR_first(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void Cc_reversed_second(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::Cc_gLR_third(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::Cc_gLR_third(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void Cc_reversed_third(ARG_LIST)
	    {
		if(sign==1)
		{
		    Dirac_algebra_type::Cc_gLR_second(factor*C_0,factor*C_1,A_0,A_2,A_1);
		    return;
		}
		if(sign==-1)
		{
		    Dirac_algebra_type::Cc_gLR_second(factor*C_1,factor*C_0,A_0,A_2,A_1);
		    return;
		}
	    }
	    static void Cc_reversed_fourth(ARG_LIST){}
    };
    template<class model_t>const std::size_t evaluate< sffLR<model_t> >::spinor_size;
    template<class model_t>const int evaluate< sffLR<model_t> >::sign;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SFFLR_H_*/


