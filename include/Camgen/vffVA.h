//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VFFVA_H_
#define CAMGEN_VFFVA_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * V-A-type vector-fermion-fermion vertex declaration and definition.  *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition of the Feynman rule class: */
    
    template<class model_t>class vffVA
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
	    
	    static const std::size_t params=2;

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
    template<class model_t>const std::size_t vffVA<model_t>::vector_size;
    template<class model_t>const std::size_t vffVA<model_t>::spinor_size;
    template<class model_t>const std::size_t vffVA<model_t>::rank;
    template<class model_t>const std::size_t vffVA<model_t>::params;
    template<class model_t>const std::size_t vffVA<model_t>::tensor_size;
    template<class model_t>const bool vffVA<model_t>::p_dependent;
    template<class model_t>const bool vffVA<model_t>::fermionic;
    template<class model_t>const std::size_t vffVA<model_t>::sizes[4]={vffVA<model_t>::vector_size,vffVA<model_t>::spinor_size,vffVA<model_t>::spinor_size,0};
    template<class model_t>const std::string vffVA<model_t>::formula="g_(mu1)(c1+c2*g_5)(i2,j3)";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< vffVA<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef vffVA<model_t> vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Dirac algebra type definition: */

	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	    /* Definition of the vector size: */

	    static const std::size_t vector_size=model_t::dimension;

	    /* Definition of the spinor size: */

	    static const std::size_t spinor_size=Dirac_algebra_type::base_type::index_range;
	    
	    /* Sign of the vector part under fermion flow reversal: */
	    
	    static const r_value_type Vsign;

	    /* Sign of the axial part under fermion flow reversal: */

	    static const r_value_type Asign;
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

	    /* Recursive relations, calling the ggVA-recursive relations in the
	     * Dirac algebra class: */

	    static void first(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_first(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void second(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_second(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void third(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_third(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void fourth(ARG_LIST){}

	    /* Right charge-conjugate recursive relations, calling the ggVA_C-recursive
	     * relations in the Dirac algebra class: */

	    static void first_C(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_C_first(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void second_C(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_C_second(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void third_C(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_C_third(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void fourth_C(ARG_LIST){}

	    /* Left charge-conjugate recursive relations, calling the Cc_ggVA-recursive
	     * relations in the Dirac algebra class: */

	    static void Cc_first(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_ggVA_first(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void Cc_second(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_ggVA_second(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void Cc_third(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_ggVA_third(factor*C_0,factor*C_1,A_0,A_1,A_2);
	    }
	    static void Cc_fourth(ARG_LIST){}

	    /* Reversed recursive relations, cllaing the ggVA-relations and
	     * multiplying the V- and A- couplings by the V- and A-reversal
	     * signs respectively: */

	    static void reversed_first(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_first(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void reversed_second(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_third(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void reversed_third(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_second(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void reversed_fourth(ARG_LIST){}

	    /* Right charge conjugate reversed recursive relations, cllaing the
	     * ggVA_C-relations and multiplying the V- and A- couplings by the
	     * V- and A-reversal signs respectively: */

	    static void reversed_first_C(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_C_first(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void reversed_second_C(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_C_third(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void reversed_third_C(ARG_LIST)
	    {
		Dirac_algebra_type::ggVA_C_second(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void reversed_fourth_C(ARG_LIST){}

	    /* Left charge conjugate reversed recursive relations, cllaing the
	     * Cc_ggVA-relations and multiplying the V- and A- couplings by the
	     * V- and A-reversal signs respectively: */

	    static void Cc_reversed_first(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_ggVA_first(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void Cc_reversed_second(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_ggVA_third(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void Cc_reversed_third(ARG_LIST)
	    {
		Dirac_algebra_type::Cc_ggVA_second(Vsign*factor*C_0,Asign*factor*C_1,A_0,A_2,A_1);
	    }
	    static void Cc_reversed_fourth(ARG_LIST){}
    };
    template<class model_t>const std::size_t evaluate< vffVA<model_t> >::vector_size;
    template<class model_t>const std::size_t evaluate< vffVA<model_t> >::spinor_size;
    template<class model_t>const typename evaluate< vffVA<model_t> >::r_value_type evaluate< vffVA<model_t> >::Vsign=evaluate< vffVA<model_t> >::Dirac_algebra_type::eta_g;
    template<class model_t>const typename evaluate< vffVA<model_t> >::r_value_type evaluate< vffVA<model_t> >::Asign=evaluate< vffVA<model_t> >::Dirac_algebra_type::eta_gg5;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VFFLR_H_*/


