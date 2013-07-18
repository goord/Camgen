//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_VVVV_H_
#define CAMGEN_VVVV_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Declaration and definition of the Yang-Mills four-vertex Feynman rule       *
 *                                                                             *
 *  2g^{mu,rho}g^{nu,sigma} - g*{mu,nu}g^{rho,sigma} - g^{mu,sigma}g^{nu,rho}  *
 *                                                                             *
 * where mu, nu, rho and sigma are the Lorentz indices of the respective first,*
 * second third and fourth interacting vector subamplitudes.                   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Feynman rule class definition: */

    template<class model_t>class vvvv
    {
	public:
	    
	    /* Reference type definition to the model type: */

	    typedef model_t model_type;

	    /* Vertex rank definition: */

	    static const std::size_t rank=4;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=(model_t::dimension)*(model_t::dimension)*(model_t::dimension)*(model_t::dimension);
	    
	    /* Boolean denoting whether the vertex tensor depends on momenta: */
	    
	    static const bool p_dependent=false;

	    /* Boolean denoting whether the vertex interacts fermions: */

	    static const bool fermionic=false;

	    /* Tensor sizes of interacting amplitudes: */

	    static const std::size_t sizes[4];

	    /* Function returning the Feynman rule formula: */

	    static const std::string formula;
    };
    template<class model_t>const std::size_t vvvv<model_t>::rank;
    template<class model_t>const std::size_t vvvv<model_t>::params;
    template<class model_t>const std::size_t vvvv<model_t>::tensor_size;
    template<class model_t>const bool vvvv<model_t>::p_dependent;
    template<class model_t>const bool vvvv<model_t>::fermionic;
    template<class model_t>const std::size_t vvvv<model_t>::sizes[4]={model_t::dimension,model_t::dimension,model_t::dimension,model_t::dimension};
    template<class model_t>const std::string vvvv<model_t>::formula="(2g(mu1,mu3)g(mu2,mu4)-g(mu1,mu2)g(mu3,mu4)-g(mu1,mu4)g(mu2,mu3))";

    /* Specialisation of the evaluate class template: */

    template<class model_t>class evaluate< vvvv<model_t> >
    {
	public:

	    /* Reference type definition of the vertex type: */

	    typedef vvvv<model_t> vertex_type;
	    
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

	    /* Function returning the index ranges of the interacting
	     * subamplitudes: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n<4)
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
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		
		value_type c0=factor*C_0;
		value_type c1=spacetime_type::dot(A_1,A_2);
		value_type c2=spacetime_type::dot(A_1,A_3);
		value_type c3=spacetime_type::dot(A_2,A_3);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_0[mu]+=c0*((r_value_type)2*c2*A_2[mu]-c1*A_3[mu]-c3*A_1[mu]);
		}
	    }
	    static void second(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		
		value_type c0=factor*C_0;
		value_type c1=spacetime_type::dot(A_0,A_2);
		value_type c2=spacetime_type::dot(A_0,A_3);
		value_type c3=spacetime_type::dot(A_2,A_3);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_1[mu]+=c0*((r_value_type)2*c1*A_3[mu]-c3*A_0[mu]-c2*A_2[mu]);
		}
	    }
	    static void third(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		
		value_type c0=factor*C_0;
		value_type c1=spacetime_type::dot(A_0,A_1);
		value_type c2=spacetime_type::dot(A_0,A_3);
		value_type c3=spacetime_type::dot(A_1,A_3);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_2[mu]+=c0*((r_value_type)2*c3*A_0[mu]-c2*A_1[mu]-c1*A_3[mu]);
		}
	    }
	    static void fourth(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<model_t::dimension),"tensor 0 iterator out of range");
		CAMGEN_ERROR_IF((A_1.range()<model_t::dimension),"tensor 1 iterator out of range");
		CAMGEN_ERROR_IF((A_2.range()<model_t::dimension),"tensor 2 iterator out of range");
		CAMGEN_ERROR_IF((A_3.range()<model_t::dimension),"tensor 3 iterator out of range");
		
		value_type c0=factor*C_0;
		value_type c1=spacetime_type::dot(A_0,A_1);
		value_type c2=spacetime_type::dot(A_0,A_2);
		value_type c3=spacetime_type::dot(A_1,A_2);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    A_3[mu]+=c0*((r_value_type)2*c2*A_1[mu]-c1*A_2[mu]-c3*A_0[mu]);
		}
	    }
    };
}


#ifdef CAMGEN_FF_CONTR_H_
#include <Camgen/gggg.h>
#endif /*CAMGEN_FF_CONTR_H_*/

#include <Camgen/undef_args.h>

#endif /*CAMGEN_VVVV_H_*/

