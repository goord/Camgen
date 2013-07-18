//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VVV_H_
#define CAMGEN_VVV_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition and declaration of the Yang-Mills 3-vertex:                    *
 *                                                                           *
 * g^{mu,nu}(p1-p2)^{rho} + g^{mu,rho}(p3-p1)^{nu} + g^{nu,rho}(p2-p3)^{mu}  *
 *                                                                           *
 * where p1, p2 and p3 are the momenta of interacting vector particles       *
 * components V1^{mu}, V2^{nu} and V3^{rho}.                                 *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Feynman rule class definition: */

    template<class model_t>class vvv
    {
	public:
	    
	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=3;
	    
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
    template<class model_t>const std::size_t vvv<model_t>::rank;
    template<class model_t>const std::size_t vvv<model_t>::params;
    template<class model_t>const std::size_t vvv<model_t>::tensor_size;
    template<class model_t>const bool vvv<model_t>::p_dependent;
    template<class model_t>const bool vvv<model_t>::fermionic;
    template<class model_t>const std::size_t vvv<model_t>::sizes[4]={model_t::dimension,model_t::dimension,model_t::dimension,0};
    template<class model_t>const std::string vvv<model_t>::formula="(g(mu1,mu2)(p1-p2)(mu3)+g(mu1,mu3)(p3-p1)(mu2)+g(mu2,mu3)(p2-p3)(mu1))";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< vvv<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef vvv<model_t> vertex_type;
	    
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
		if(n<3)
		{
		    v.push_back(model_t::dimension);
		}
		return v;
	    }

	    /* Recursive relations: */

	    static void first(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		
		value_type c0=factor*C_0;
		value_type c1=(r_value_type)2*spacetime_type::dot(P_1,A_2)+spacetime_type::dot(P_2,A_2);
		value_type c2=(r_value_type)2*spacetime_type::dot(P_2,A_1)+spacetime_type::dot(P_1,A_1);
		value_type c3=spacetime_type::dot(A_1,A_2);

		for(std::size_t mu=0;mu<model_t::dimension;++mu)
		{
		    A_0[mu]+=c0*(-c1*A_1[mu]+c3*(P_1[mu]-P_2[mu])+c2*A_2[mu]);
		}
		
	    }
	    static void second(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		
		value_type c0=factor*C_0;
		value_type c1=(r_value_type)2*spacetime_type::dot(P_0,A_2)+spacetime_type::dot(P_2,A_2);
		value_type c2=(r_value_type)2*spacetime_type::dot(P_2,A_0)+spacetime_type::dot(P_0,A_0);
		value_type c3=spacetime_type::dot(A_0,A_2);
		
		for(std::size_t mu=0;mu<model_t::dimension;++mu)
		{
		    A_1[mu]-=c0*(-c1*A_0[mu]+c3*(P_0[mu]-P_2[mu])+c2*A_2[mu]);
		}
		
	    }
	    static void third(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		
		value_type c0=factor*C_0;
		value_type c1=(r_value_type)2*spacetime_type::dot(P_0,A_1)+spacetime_type::dot(P_1,A_1);
		value_type c2=(r_value_type)2*spacetime_type::dot(P_1,A_0)+spacetime_type::dot(P_0,A_0);
		value_type c3=spacetime_type::dot(A_0,A_1);
		
		for(std::size_t mu=0;mu<model_t::dimension;++mu)
		{
		    A_2[mu]+=c0*(-c1*A_0[mu]+c3*(P_0[mu]-P_1[mu])+c2*A_1[mu]);
		}
		
	    }
	    static void fourth(ARG_LIST){}
    };
}

/* If the structure constant colour structure is included by the user, we include the
 * specialisation of the Feynman rule composition of the Yang-Mills 3-vertex with a
 * structure constants tensor: */

#ifdef CAMGEN_F_H_
#include <Camgen/ggg.h>
#endif /*CAMGEN_F_H_*/

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VVV_H_*/

