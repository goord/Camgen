//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SSSS_H_
#define CAMGEN_SSSS_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the 4-scalar vertex.  *
 *                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Feynman rule class definition: */

    template<class model_t>class ssss
    {
	public:

	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=4;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=1;
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=false;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=false;
	    
	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];

	    /* Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t ssss<model_t>::rank;
    template<class model_t>const std::size_t ssss<model_t>::params;
    template<class model_t>const std::size_t ssss<model_t>::tensor_size;
    template<class model_t>const bool ssss<model_t>::p_dependent;
    template<class model_t>const bool ssss<model_t>::fermionic;
    template<class model_t>const std::size_t ssss<model_t>::sizes[4]={1,1,1,1};
    template<class model_t>const std::string ssss<model_t>::formula="1";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< ssss<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef ssss<model_t> vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	    /* Initialisation phase: */
	    
	    static void initialise(){}

	    /* function returning the index ranges of the interacting
	     * subamplitudes: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		return v;
	    }

	    /* Recursive relations: */

	    static void first(ARG_LIST)
	    {
		(*A_0)+=factor*C_0*(*A_1)*(*A_2)*(*A_3);
	    }
	    static void second(ARG_LIST)
	    {
		(*A_1)+=factor*C_0*(*A_0)*(*A_2)*(*A_3);
	    }
	    static void third(ARG_LIST)
	    {
		(*A_2)+=factor*C_0*(*A_0)*(*A_1)*(*A_3);
	    }
	    static void fourth(ARG_LIST)
	    {
		(*A_3)+=factor*C_0*(*A_0)*(*A_1)*(*A_2);
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SSSS_H_*/

