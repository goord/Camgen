//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SVV3_H_
#define CAMGEN_SVV3_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Decalaration and definition of the scalar-vector-vector vertex. *
 * Mirroring the VVS3 Lorentz structure from MadGraph 5            *
 *        p1^{mu}*p2^{nu} - p1.p2*g^{mu,nu}			   *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Feynman rule class definition: */

    template<class model_t>class svv3
    {
	public:

	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=3;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=(model_t::dimension)*(model_t::dimension);
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=true;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=false;
	    
	    /* Initialisation phase: */
	    
	    static void initialise(){}

	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];

	    /* Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t svv3<model_t>::rank;
    template<class model_t>const std::size_t svv3<model_t>::params;
    template<class model_t>const std::size_t svv3<model_t>::tensor_size;
    template<class model_t>const bool svv3<model_t>::p_dependent;
    template<class model_t>const bool svv3<model_t>::fermionic;
    template<class model_t>const std::size_t svv3<model_t>::sizes[4]={1,model_t::dimension,model_t::dimension};
    template<class model_t>const std::string svv3<model_t>::formula="(p1(mu1)*p2(mu2)-g(mu1,mu2)*g(p1,p2))";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< svv3<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef svv3<model_t> vertex_type;
	    
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
		if(n==1 or n==2)
		{
		    v.push_back(model_t::dimension);
		}
		return v;
	    }

	    /* Scalar recursive relation: */
	    
	    static void first(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		(*A_0)+=factor*C_0*(spacetime_type::dot(P_1,A_1)*spacetime_type::dot(P_2,A_2) - spacetime_type::dot(P_1,P_2)*spacetime_type::dot(A_1,A_2));
	    }

	    /* First vector recursive relation: */
	    
	    static void second(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		value_type c=factor*C_0*(*A_0);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_1[mu]+= -c*(P_0[mu]+P_2[mu])*spacetime_type::dot(P_2,A_2);
		    A_1[mu]+= c*spacetime_type::dot(P_0,P_2)*A_2[mu];
		    A_1[mu]+= c*spacetime_type::dot(P_2,P_2)*A_2[mu];
		}
	    }

	    /* Second vector recursive relation: */
	    
	    static void third(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		value_type c=factor*C_0*(*A_0);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_2[mu]+= -c*(P_0[mu]+P_1[mu])*spacetime_type::dot(P_1,A_1);
		    A_2[mu]+= c*spacetime_type::dot(P_0,P_1)*A_1[mu];
		    A_2[mu]+= c*spacetime_type::dot(P_1,P_1)*A_1[mu];
		}
	    }
	    static void fourth(ARG_LIST){}
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SVV3_H_*/

