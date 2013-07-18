//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VSS_H_
#define CAMGEN_VSS_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the vector-scalar-scalar vertex, with Feynman   *
 * rule                                                                          *
 *                                                                               *
 * 				p1^{mu}-p2*{mu}                                  *
 *                                                                               *
 * where p1 and p2 are the momenta of the respective scalar fields and mu is the *
 * Lorentz index of the vector field.                                            *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class vss
    {
	public:

	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=3;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=model_t::dimension;
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=true;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=false;

	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];

	    /* Function returning the Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t vss<model_t>::rank;
    template<class model_t>const std::size_t vss<model_t>::params;
    template<class model_t>const std::size_t vss<model_t>::tensor_size;
    template<class model_t>const bool vss<model_t>::p_dependent;
    template<class model_t>const bool vss<model_t>::fermionic;
    template<class model_t>const std::size_t vss<model_t>::sizes[4]={model_t::dimension,1,1,0};
    template<class model_t>const std::string vss<model_t>::formula="(p2(mu1)-p3(mu1))";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< vss<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef vss<model_t> vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);

	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
	    }

	    /* function returning the index ranges of the interacting
	     * subamplitudes: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n==0)
		{
		    v.push_back(model_t::dimension);
		}
		return v;
	    }

	    /* Recursion relation evaluating the vector subamplitude: */

	    static void first(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");

		value_type c=factor*C_0*(*A_1)*(*A_2);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_0[mu]+=c*(P_1[mu]-P_2[mu]);
		}
	    }

	    /* Recursion relation evaluating the first scalar subamplitude: */

	    static void second(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");

		momentum_type Q=P_0+(r_value_type)2*P_2;
		(*A_1)-=factor*C_0*(*A_2)*(spacetime_type::dot(A_0,Q));
	    }

	    /* Recursion relation evaluating the second scalar subamplitude: */

	    static void third(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");

		momentum_type Q=P_0+(r_value_type)2*P_1;
		(*A_2)+=factor*C_0*(*A_1)*(spacetime_type::dot(A_0,Q));
	    }
	    static void fourth(ARG_LIST){}

    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VSS_H_*/


