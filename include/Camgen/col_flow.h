//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_COL_FLOW_H_
#define CAMGEN_COL_FLOW_H_

#include <Camgen/su(n).h>
#include <Camgen/adj_rep.h>
#include <Camgen/comp_vert.h>
#include <Camgen/comp_contr.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Colour_flow colour treatment type definition. When dealing with coloured      *
 * particles, a colour treatment type has to be defined in the model class. The  *
 * colour flow type expands SU(N)-adjoint colour indices into doublets of        *
 * colour-anti-colour indices in the fundamental representation. The colour      *
 * structures in Feynman rules are replaced with their 'CF_type', replacing      *
 * SU(N) structure constants, generators and delta functions with their          *
 * corresponding expressions in the colour flow representation.                  *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class colour_flow
    {
	public:

	    /* The apply metafunction maps representations to representations. In the
	     * colour flow treatment it maps by default a representation to itself. */

	    template<class rep_t>class apply
	    {
		public:

		    /* The type 'type' denotes the mapped representation. In the
		     * default case the representation itself. */

		    typedef rep_t type;

		    /* The integer 'factor' denotes the number of indices a rep-type
		     * index has to be replaced with. In the default case one. */

		    static const std::size_t factor=1;

		    /* The boolean 'decomposes' denotes whether for discrete colours
		     * Camgen has to keep track of internal generated colour singlet
		     * tensors. In the default case not. */

		    static const bool decomposes=false;

		    /* The static function 'fill_colour_numbers' denotes if
		     * anti-colours are assigned to indices, in which case -1 will be
		     * inserted in the argument vector. Otherwise, +1 is inserted. In
		     * the default case no anticolours are assigned. */

		    static void fill_colour_numbers(std::vector<int>& v)
		    {
			v.insert(v.begin(),1);
		    }
	    };

	    /* The 'make_Feynrule' metafunction class converts a Feynman rule
	     * (template parameter) to 'type'. In the default case the Feynman rule
	     * is passed as it is. */
	    
	    template<class Feynrule_t>class make_Feynrule
	    {
		public:
		    typedef Feynrule_t type;
	    };

	    /* The 'make_contraction' metafunction class converts a colour index
	     * contraction (template parameter) to 'type'. In the default case the
	     * colour contraction is passed as it is. */

	    template<class contr_t>class make_contraction
	    {
		public:
		    typedef contr_t type;
	    };

	    /* The 'decomposes' metafunction class tells Camgen whether it should
	     * keep track of internal generate colour singlets for discrete colour
	     * computations. Not in the default case. */

	    template<bool cont_cols>class decomposes
	    {
		public:
		    static const bool value=!cont_cols;
	    };
    };
    template<class rep_t>const std::size_t colour_flow::template apply<rep_t>::factor;
    template<class rep_t>const bool colour_flow::template apply<rep_t>::decomposes;
    template<bool cont_cols>const bool colour_flow::template decomposes<cont_cols>::value;

    /* Specialisation of the apply metafunction: the SU(N) adjoint representation is
     * mapped to a doublet of fundamental representations. In the discrete colour
     * case, Camgen will decompose these doublets and keep track of propagating
     * colour-anticolour modes.*/

    template<std::size_t N>class colour_flow::apply< adjoint_rep< SU<N> > >
    {
	public:
	    typedef fundamental_rep< SU<N> > type;
	    static const std::size_t factor=2;
	    static const bool decomposes=true;
	    static void fill_colour_numbers(std::vector<int>& v)
	    {
		v.insert(v.begin(),-1);
		v.insert(v.begin(),1);
	    }
    };
    template<std::size_t N>const std::size_t colour_flow::template apply< adjoint_rep< SU<N> > >::factor;
    template<std::size_t N>const bool colour_flow::template apply< adjoint_rep< SU<N> > >::decomposes;

    /* Specialisation of the apply metafunction: the SU(N) fundamental representation
     * decomposes in the discrete colours case. */

    template<std::size_t N>class colour_flow::apply< fundamental_rep< SU<N> > >
    {
	public:
	    typedef fundamental_rep< SU<N> > type;
	    static const std::size_t factor=1;
	    static const bool decomposes=true;
	    static void fill_colour_numbers(std::vector<int>& v)
	    {
		v.insert(v.begin(),1);
	    }
    };
    template<std::size_t N>const std::size_t colour_flow::template apply< fundamental_rep< SU<N> > >::factor;
    template<std::size_t N>const bool colour_flow::template apply< fundamental_rep< SU<N> > >::decomposes;

    /* Specialisation of the make_Feynrule metafunction: the colour part (first
     * template argument in compose_vertices) is replaced with its CF-type. */

    template<class Feynrule_t1,class Feynrule_t2>class colour_flow::make_Feynrule< compose_vertices<Feynrule_t1,Feynrule_t2> >
    {
	public:
	    typedef compose_vertices<typename Feynrule_t1::CF_type,typename colour_flow::template make_Feynrule<Feynrule_t2>::type> type;
    };

    /* Specialisation of the make_contraction metafunction: the colour part (first
     * template argument in compose_vertices) is replaced with its CF-type. */

    template<class contr_t1,class contr_t2>class colour_flow::make_contraction< compose_contraction<contr_t1,contr_t2> >
    {
	public:
	    typedef compose_contraction<typename contr_t1::CF_type,typename colour_flow::template make_contraction<contr_t2>::type> type;
    };
}

#endif /*CAMGEN_COL_FLOW_H_*/

