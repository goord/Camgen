//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_ADJOINT_H_
#define CAMGEN_ADJOINT_H_

#include <vector>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Adjoint colour treatment type definition. When dealing with coloured          *
 * particles, a colour treatment type has to be defined in the model class. The  *
 * adjoint type is the 'default' type in the sense that adjoint                  *
 * representations are not expanded into fundamental representation doublets.    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class adjoint
    {
	public:

	    /* The apply metafunction maps representations to representations. In the
	     * adjoint colour treatment it maps every representation to itself. */

	    template<class rep_t>class apply
	    {
		public:

		    /* The type 'type' denotes the mapped representation. In this
		     * case the representation itself. */

		    typedef rep_t type;

		    /* The integer 'factor' denotes the number of indices a rep-type
		     * index has to be replaced with. */

		    static const std::size_t factor=1;

		    /* The boolean 'decomposes' denotes whether for discrete colours
		     * Camgen has to keep track of internal generated colour singlet
		     * tensors. */

		    static const bool decomposes=false;

		    /* The static function 'fill_colour_numbers' denotes if
		     * anti-colours are assigned to indices, in which case a -1 will
		     * be inserted in the argument vector. Otherwise, +1 is inserted.
		     * In the adjoint colour treatment no anticolours are assigned.
		     * */

		    static void fill_colour_numbers(std::vector<int>& v)
		    {
			v.insert(v.begin(),1);
		    }
	    };

	    /* The 'make_Feynrule' metafunction class converts a Feynman rule
	     * (template parameter) to 'type'. In the adjoint colour treatment the
	     * Feynman rule is passed as it is. */

	    template<class Feynrule_t>class make_Feynrule
	    {
		public:
		    typedef Feynrule_t type;
	    };

	    /* The 'make_contraction' metafunction class converts a colour index
	     * contraction (template parameter) to 'type'. In the adjoint colour
	     * treatment the colour contraction is passed as it is. */

	    template<class contr_t>class make_contraction
	    {
		public:
		    typedef contr_t type;
	    };

	    /* The 'decomposes' metafunction class tells Camgen whether it should
	     * keep track of internal generated colour singlets for discrete colour
	     * computations. Never in the djoint colour treatment. */

	    template<bool cont_cols>class decomposes
	    {
		public:
		    static const bool value=false;
	    };
    };
    template<class rep_t>const std::size_t adjoint::template apply<rep_t>::factor;
    template<class rep_t>const bool adjoint::template apply<rep_t>::decomposes;
    template<bool cont_cols>const bool adjoint::template decomposes<cont_cols>::value;
}

#endif /*CAMGEN_ADJOINT_H_*/

