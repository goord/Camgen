//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VFF_H_
#define CAMGEN_VFF_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Vector-fermion-fermion Feynman rule declaration and definition. * 
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition of the Feynman rule class: */
    
    template<class model_t>class vff
    {
	public:

	    /* Reference type definition to the model type: */

	    typedef model_t model_type;
	    
	    /* Dirac algebra type definition: */
	    
	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);
	
	private:

	    /* Definition of the vector size: */

	    static const std::size_t vector_size=model_t::dimension;

	    /* Definition of the spinor size: */

	    static const std::size_t spinor_size=Dirac_algebra_type::base_type::index_range;

	public:

	    /* Vertex rank definition: */

	    static const std::size_t rank=3;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=vector_size*spinor_size*spinor_size;
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=false;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=true;

	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];

	    /* Function returning the Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t vff<model_t>::vector_size;
    template<class model_t>const std::size_t vff<model_t>::spinor_size;
    template<class model_t>const std::size_t vff<model_t>::rank;
    template<class model_t>const std::size_t vff<model_t>::params;
    template<class model_t>const std::size_t vff<model_t>::tensor_size;
    template<class model_t>const bool vff<model_t>::p_dependent;
    template<class model_t>const bool vff<model_t>::fermionic;
    template<class model_t>const std::size_t vff<model_t>::sizes[4]={vff<model_t>::vector_size,vff<model_t>::spinor_size,vff<model_t>::spinor_size,0};
    template<class model_t>const std::string vff<model_t>::formula="g_mu1(i2,j3)";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< vff<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef vff<model_t> vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Dirac algebra type definition: */

	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	    /* Definition of the vector size: */

	    static const std::size_t vector_size=model_t::dimension;

	    /* Definition of the spinor size: */

	    static const std::size_t spinor_size=Dirac_algebra_type::base_type::index_range;
	    
	    /* Sign of the vertex under fermion flow reversal: */
	    
	    static const r_value_type sign;
	public:
	    
	    /* Initialisation phase: */
	    
	    static void initialise()
	    {
		Dirac_algebra_type::initialise();
	    }

	    /* Function returning the index ranges of the interacting
	     * subamplitudes: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n==0)
		{
		    v.push_back(vector_size);
		}
		if(n==1 or n==2)
		{
		    v.push_back(spinor_size);
		}
		return v;
	    }

	    /* Recursive relations, calling the g-recursive relations in the
	     * Dirac algebra class: */

	    static void first(ARG_LIST)
	    {
		Dirac_algebra_type::g_first(factor*C_0,A_0,A_1,A_2);
	    }
	    static void second(ARG_LIST)
	    {
		Dirac_algebra_type::g_second(factor*C_0,A_0,A_1,A_2);
	    }
	    static void third(ARG_LIST)
	    {
		Dirac_algebra_type::g_third(factor*C_0,A_0,A_1,A_2);
	    }
	    static void fourth(ARG_LIST){}

	    /* Right charge-conjugate recursive relations, calling the g_C-recursive
	     * relations in the Dirac algebra class: */

	    static void first_C(ARG_LIST)
	    {
		Dirac_algebra_type::g_C_first(factor*C_0,A_0,A_1,A_2);
	    }
	    static void second_C(ARG_LIST)
	    {
		Dirac_algebra_type::g_C_second(factor*C_0,A_0,A_1,A_2);
	    }
	    static void third_C(ARG_LIST)
	    {
		Dirac_algebra_type::g_C_third(factor*C_0,A_0,A_1,A_2);
	    }
	    static void fourth_C(ARG_LIST){}

	    static void Cc_first(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_g_first(factor*C_0,A_0,A_1,A_2);
	    }
	    static void Cc_second(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_g_second(factor*C_0,A_0,A_1,A_2);
	    }
	    static void Cc_third(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_g_third(factor*C_0,A_0,A_1,A_2);
	    }
	    static void Cc_fourth(ARG_LIST){}

	    /* Reversed recursive relations, multiplying the g-recursive
	     * relations with the reversal sign: */

	    static void reversed_first(ARG_LIST)
	    {
		Dirac_algebra_type::g_first(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void reversed_second(ARG_LIST)
	    {
		Dirac_algebra_type::g_third(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void reversed_third(ARG_LIST)
	    {
		Dirac_algebra_type::g_second(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void reversed_fourth(ARG_LIST){}

	    /* Right charge-conjugated reversed recursive relations, multiplying
	     * the g_C-recursive relations with the reversal sign: */

	    static void reversed_first_C(ARG_LIST)
	    {
		Dirac_algebra_type::g_C_first(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void reversed_second_C(ARG_LIST)
	    {
		Dirac_algebra_type::g_C_third(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void reversed_third_C(ARG_LIST)
	    {
		Dirac_algebra_type::g_C_second(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void reversed_fourth_C(ARG_LIST){}

	    /* Left charge-conjugated reversed recursive relations, multiplying
	     * the g_C-recursive relations with the reversal sign: */

	    static void Cc_reversed_first(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_g_first(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void Cc_reversed_second(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_g_third(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void Cc_reversed_third(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_g_second(sign*factor*C_0,A_0,A_2,A_1);
	    }
	    static void Cc_reversed_fourth(ARG_LIST){}
    };
    template<class model_t>const std::size_t evaluate< vff<model_t> >::vector_size;
    template<class model_t>const std::size_t evaluate< vff<model_t> >::spinor_size;
    template<class model_t>const typename evaluate< vff<model_t> >::r_value_type evaluate< vff<model_t> >::sign=evaluate< vff<model_t> >::Dirac_algebra_type::eta_g;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VFF_H_*/


