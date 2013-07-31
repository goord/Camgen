//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SSVV3_H_
#define CAMGEN_SSVV3_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Decalaration and definition of the scalar-scalar-vector-vector vertex.  *
 * Mirroring the VVSS3 Lorentz structure from MadGraph 5                   *
 *        p1^{mu}*p2^{nu} - p1.p2*g^{mu,nu}			           *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Feynman rule class definition: */

    template<class model_t>class ssvv3
    {
	public:

	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=4;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=(model_t::dimension)*(model_t::dimension);
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=true;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=false;

	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];

	    /* Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t ssvv3<model_t>::rank;
    template<class model_t>const std::size_t ssvv3<model_t>::params;
    template<class model_t>const std::size_t ssvv3<model_t>::tensor_size;
    template<class model_t>const bool ssvv3<model_t>::p_dependent;
    template<class model_t>const bool ssvv3<model_t>::fermionic;
    template<class model_t>const std::size_t ssvv3<model_t>::sizes[4]={1,1,model_t::dimension,model_t::dimension};
    template<class model_t>const std::string ssvv3<model_t>::formula="(p2(mu2)*p3(mu3)-g(mu2,mu3)*g(p2,p3))";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< ssvv3<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef ssvv3<model_t> vertex_type;
	    
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
		if(n==2 or n==3)
		{
		    v.push_back(model_t::dimension);
		}
		return v;
	    }

	    /* First scalar recursive relation: */
	    
	    static void first(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		(*A_0)+=factor*C_0*(*A_1)*(spacetime_type::dot(P_2,A_2)*spacetime_type::dot(P_3,A_3) - spacetime_type::dot(P_2,P_3)*spacetime_type::dot(A_2,A_3));
	    }

	    /* Second scalar recursive relation: */
	    
	    static void second(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		(*A_1)+=factor*C_0*(*A_0)*(spacetime_type::dot(P_2,A_2)*spacetime_type::dot(P_3,A_3) - spacetime_type::dot(P_2,P_3)*spacetime_type::dot(A_2,A_3));
	    }

	    /* First vector recursive relation: */
	    
	    static void third(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		value_type c=factor*C_0*(*A_0)*(*A_1);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_2[mu]+= -c*(P_0[mu]+P_1[mu]+P_3[mu])*spacetime_type::dot(P_3,A_3);
		    A_2[mu]+= c*spacetime_type::dot(P_0,P_3)*A_3[mu];
		    A_2[mu]+= c*spacetime_type::dot(P_1,P_3)*A_3[mu];
		    A_2[mu]+= c*spacetime_type::dot(P_3,P_3)*A_3[mu];
		}
	    }

	    /* Second vector recursive relation: */
	    
	    static void fourth(ARG_LIST)
	    {
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		value_type c=factor*C_0*(*A_0)*(*A_1);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_3[mu]+= -c*(P_0[mu]+P_1[mu]+P_2[mu])*spacetime_type::dot(P_2,A_2);
		    A_3[mu]+= c*spacetime_type::dot(P_0,P_2)*A_2[mu];
		    A_3[mu]+= c*spacetime_type::dot(P_1,P_2)*A_2[mu];
		    A_3[mu]+= c*spacetime_type::dot(P_2,P_2)*A_2[mu];
		}
	    }
    };	
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SSVV3_H_*/

