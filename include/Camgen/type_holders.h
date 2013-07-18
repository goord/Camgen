//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_TYPE_HOLDERS_H_
#define CAMGEN_TYPE_HOLDERS_H_

#include <set>
#include <Camgen/unused.h>
#include <Camgen/tensor.h>
#include <Camgen/forward_decs.h>
#include <Camgen/phase_space.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility class templates with the purpose of defining implementation classes   *
 * from the model classes static data. Every model class should at least contain *
 * a type definition of the numerical type called 'value_type', for which the    *
 * usual arithmetic operations and std::sqrt should be defined, a const static   *
 * int called 'dimension', denoting the spacetime dimension, and a const static  *
 * bool called 'coloured' denoting whether the model contains unbroken           *
 * nonabelian symmetries. If the dimension is bigger than zero, a spacetime type *
 * must by given, and a static boolean called 'continuous_helicities' denoting   *
 * whether the external states are generally not helicity eigenstates. Moreover, *
 * if fermions are present, a definition of the Dirac algebra basis should be    *
 * contained in the model as well as a const static integer 'beam_direction',    *
 * and optionally if these are massive, a definition of the spin vector type. On *
 * the other hand, if the 'coloured' boolean is true, a const static boolean     *
 * called 'continuous_colours', denoting whether the colour states are           *
 * eigenstates or not. Finally, one should then also define a type called        *
 * 'colour treatment', which determines what colour basis is used by Camgen.    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Class template defining the spacetime implementation class: */

    template<class model_t,std::size_t dim=model_t::dimension>class get_spacetime_type
    {
	public:
	    typedef typename model_t::spacetime_type signature;
	    typedef typename signature::template implementation<typename model_t::value_type,model_t::dimension> type;

	    static const bool continuous_helicities=model_t::continuous_helicities;
    };
    template<class model_t,std::size_t dim>const bool get_spacetime_type<model_t,dim>::continuous_helicities;

    /* Trivial specialisation for zero dimensions: */

    template<class model_t>class get_spacetime_type<model_t,0>
    {
	public:
	    typedef unused signature;
	    typedef unused type;

	    static const bool continuous_helicities=false;
    };
    template<class model_t>const bool get_spacetime_type<model_t,0>::continuous_helicities;

    /* Class template defining the Dirac algebra implementation class: */

    template<class model_t>class get_Dirac_algebra_type
    {
	public:
	    typedef typename model_t::Dirac_algebra_type basis_type;
	    typedef typename basis_type::template implementation<typename model_t::value_type,model_t::dimension> type;
    };

    /* Class template defining the colour treatment class: */

    template<class model_t,bool cols=model_t::coloured>class get_colour_treatment;
    
    template<class model_t>class get_colour_treatment<model_t,true>
    {
	public:
	    typedef typename model_t::colour_treatment type;

	    static const bool continuous_colours=model_t::continuous_colours;
	    static const bool decomposes=type::template decomposes<continuous_colours>::value;
    };
    template<class model_t>const bool get_colour_treatment<model_t,true>::continuous_colours;
    template<class model_t>const bool get_colour_treatment<model_t,true>::decomposes;

    /* Compile-time error specialisation for uncoloured models: */

    template<class model_t>class get_colour_treatment<model_t,false>
    {
	public:
	    typedef unused type;

	    static const bool continuous_colours=false;
	    static const bool decomposes=false;
    };
    template<class model_t>const bool get_colour_treatment<model_t,false>::continuous_colours;
    template<class model_t>const bool get_colour_treatment<model_t,false>::decomposes;

    /* Class template defining the particle phase space type: */

    template<class model_t>class get_ps_type
    {
	public:
	    typedef particle_ps<model_t,model_t::dimension,model_t::coloured> type;
    };

    /* Class template defining the implementation of a group or representation
     * wrapper class: */

    template<class model_t,class group_t>class get_implementation
    {
	public:
	    typedef typename group_t::template implementation<typename model_t::value_type> type;
    };

    /* Class template defining the implementation of vertex Feynman rules: */

    template<class model_t,bool decomp=get_colour_treatment<model_t,model_t::coloured>::decomposes>class get_recursive_relation;

    /* Specialisation for models with no colour decomposition: */

    template<class model_t>class get_recursive_relation<model_t,false>
    {
	private:
	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef typename tensor<value_type>::iterator iterator;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;
	public:
	    typedef void(*vert_func)(const value_type&,const std::vector<const value_type*>&,std::vector<iterator>&,const std::vector<const momentum_type*>&);
	    
	    template<class Feynrule_t>class apply
	    {
		public:
		    typedef evaluate<Feynrule_t> type;
	    };
    };

    /* Specialisation for models with colour decomposition: */

    template<class model_t>class get_recursive_relation<model_t,true>
    {
	private:
	    typedef typename model_t::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef typename tensor<value_type>::iterator iterator;
	    typedef vector<r_value_type,model_t::dimension> momentum_type;
	public:
	    typedef void(*vert_func)(const value_type&,const std::vector<const value_type*>&,std::vector<iterator>&,const std::vector<const momentum_type*>&,std::set<iterator>&);	    
	    template<class Feynrule_t>class apply
	    {
		public:
		    typedef cfd_evaluate<Feynrule_t> type;
	    };
    };

    /* Class template defining all useful types corresponding to a model: */

    template<class model_t>class get_basic_types
    {
	public:
	    typedef model_t model_type;
	    typedef typename model_type::value_type r_value_type;
	    typedef std::complex<r_value_type> value_type;
	    typedef tensor<value_type> tensor_type;
	    typedef typename tensor_type::size_type size_type;
	    typedef typename tensor_type::iterator iterator;
	    typedef typename tensor_type::const_iterator const_iterator;
	    typedef vector<value_type,model_type::dimension> vector_type;
	    typedef vector<r_value_type,model_type::dimension> momentum_type;
	    typedef particle<model_type> particle_type;
	    typedef vertex<model_type> vertex_type;
	    typedef particle_ps<model_type,model_type::dimension,model_type::coloured> phase_space_type;
	    typedef typename phase_space_type::wave_func wave_func;
	    typedef typename phase_space_type::contr_func contr_func;
	    typedef void(*prop_func)(iterator,iterator);
	    typedef void(*prop_refresher)(const momentum_type*,const r_value_type*,const r_value_type*);
	    typedef typename get_recursive_relation<model_t,get_colour_treatment<model_t,model_t::coloured>::decomposes>::vert_func vert_func;
    };
}

#endif /*CAMGEN_TYPE_HOLDERS_H_*/

