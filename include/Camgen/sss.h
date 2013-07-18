//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SSS_H_
#define CAMGEN_SSS_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the 3-scalar vertex.  *
 *                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Feynman rule class definition: */

    template<class model_t>class sss
    {
	public:

	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=3;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=1;
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=false;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=false;

	    /* Interacting subamplitude sizes: */

	    static const std::size_t sizes[4];

	    /* Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t sss<model_t>::rank;
    template<class model_t>const std::size_t sss<model_t>::params;
    template<class model_t>const std::size_t sss<model_t>::tensor_size;
    template<class model_t>const bool sss<model_t>::p_dependent;
    template<class model_t>const bool sss<model_t>::fermionic;
    template<class model_t>const std::size_t sss<model_t>::sizes[4]={1,1,1,0};
    template<class model_t>const std::string sss<model_t>::formula="1";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< sss<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef sss<model_t> vertex_type;
	    
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
		(*A_0)+=C_0*factor*(*A_1)*(*A_2);
	    }
	    static void second(ARG_LIST)
	    {
		(*A_1)+=C_0*factor*(*A_0)*(*A_2);
	    }
	    static void third(ARG_LIST)
	    {
		(*A_2)+=C_0*factor*(*A_0)*(*A_1);
	    }
	    static void fourth(ARG_LIST){}
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SSS_H_*/

