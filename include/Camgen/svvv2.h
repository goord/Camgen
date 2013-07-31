//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SVVV2_H_
#define CAMGEN_SVVV2_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition and declaration of the Yang-Mills 3-vertex together with a     *
 * scalar current:                                                           *
 *                                                                           *
 * g^{mu,nu}(p2-p3)^{rho} + g^{mu,rho}(p4-p2)^{nu} + g^{nu,rho}(p3-p4)^{mu}  *
 *                                                                           *
 * where p2, p3 and p4 are the momenta of interacting vector particles       *
 * components V2^{mu}, V3^{nu} and V4^{rho}.                                 *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Feynman rule class definition: */

    template<class model_t>class svvv2
    {
	public:
	    
	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=4;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=(model_t::dimension)*(model_t::dimension)*(model_t::dimension);
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=true;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=false;

	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];

	    /* Function returning the Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t svvv2<model_t>::rank;
    template<class model_t>const std::size_t svvv2<model_t>::params;
    template<class model_t>const std::size_t svvv2<model_t>::tensor_size;
    template<class model_t>const bool svvv2<model_t>::p_dependent;
    template<class model_t>const bool svvv2<model_t>::fermionic;
    template<class model_t>const std::size_t svvv2<model_t>::sizes[4]={1,model_t::dimension,model_t::dimension,model_t::dimension};
    template<class model_t>const std::string svvv2<model_t>::formula="(g(mu2,mu3)(p2-p3)(mu4)+g(mu2,mu3)(p4-p2)(mu3)+g(mu3,mu4)(p3-p4)(mu2))";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< svvv2<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef svvv2<model_t> vertex_type;
	    
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
		if((n == 1) || (n == 2) || (n == 3))
		{
		    v.push_back(model_t::dimension);
		}
		return v;
	    }

	    /* Recursive relations: */

	    static void first(ARG_LIST)
	    {
	        CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
	        CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
	        CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");

		value_type c0=factor*C_0;
		value_type c1=spacetime_type::dot(A_1,A_2)*(spacetime_type::dot(P_1,A_3)-spacetime_type::dot(P_2,A_3));
		value_type c2=spacetime_type::dot(A_1,A_3)*(spacetime_type::dot(P_3,A_2)-spacetime_type::dot(P_1,A_2));
		value_type c3=spacetime_type::dot(A_2,A_3)*(spacetime_type::dot(P_2,A_1)-spacetime_type::dot(P_3,A_1));

		(*A_0)=c0*(c1+c2+c3);

	    }
	    static void second(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		
		value_type c0=factor*C_0*(*A_0);
		value_type c1=(r_value_type)2*spacetime_type::dot(P_2,A_3)+spacetime_type::dot(P_3,A_3)+spacetime_type::dot(P_0,A_3);
		value_type c2=(r_value_type)2*spacetime_type::dot(P_3,A_2)+spacetime_type::dot(P_2,A_2)+spacetime_type::dot(P_0,A_2);
		value_type c3=spacetime_type::dot(A_2,A_3);

		for(std::size_t mu=0;mu<model_t::dimension;++mu)
		{
		    A_1[mu]+=c0*(-c1*A_2[mu]+c3*(P_2[mu]-P_3[mu])+c2*A_3[mu]);
		}
		
	    }
	    static void third(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		
		value_type c0=factor*C_0*(*A_0);
		value_type c1=(r_value_type)2*spacetime_type::dot(P_1,A_3)+spacetime_type::dot(P_3,A_3)+spacetime_type::dot(P_0,A_3);
		value_type c2=(r_value_type)2*spacetime_type::dot(P_3,A_1)+spacetime_type::dot(P_1,A_1)+spacetime_type::dot(P_0,A_1);
		value_type c3=spacetime_type::dot(A_1,A_3);
		
		for(std::size_t mu=0;mu<model_t::dimension;++mu)
		{
		    A_2[mu]-=c0*(-c1*A_1[mu]+c3*(P_1[mu]-P_3[mu])+c2*A_3[mu]);
		}
		
	    }
	    static void fourth(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		
		value_type c0=factor*C_0*(*A_0);
		value_type c1=(r_value_type)2*spacetime_type::dot(P_1,A_2)+spacetime_type::dot(P_2,A_2)+spacetime_type::dot(P_0,A_2);
		value_type c2=(r_value_type)2*spacetime_type::dot(P_2,A_1)+spacetime_type::dot(P_1,A_1)+spacetime_type::dot(P_0,A_1);
		value_type c3=spacetime_type::dot(A_1,A_2);
		
		for(std::size_t mu=0;mu<model_t::dimension;++mu)
		{
		    A_3[mu]+=c0*(-c1*A_1[mu]+c3*(P_1[mu]-P_2[mu])+c2*A_2[mu]);
		}
		
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SVVV2_H_*/

