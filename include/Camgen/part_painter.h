//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PART_PAINTER_H_
#define CAMGEN_PART_PAINTER_H_

#include <cstdlib>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Utility class template definitions for the implementation of coloured *
 * particles in Camgen.                                                  *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Utility class template to define coloured particles in Camgen: */

    template<class rep_t,std::size_t N=rep_t::multiplicity>class particle_painter
    {
	public:

	    template<class particle_t>static void apply(particle_t* phi)
	    {
		particle_painter<typename rep_t::back,1>::template apply<particle_t>(phi);
		particle_painter<typename rep_t::min_back,N-1>::template apply<particle_t>(phi);
	    }
    };

    /* Specialisation for representation multiplicity 1: */

    template<class rep_t>class particle_painter<rep_t,1>
    {
	public:

	    template<class particle_t>static void apply(particle_t* phi)
	    {
		/* Determining the properties of the representation under
		 * model_t's colour treatment: */

		typename particle_t::size_type n=particle_t::model_type::colour_treatment::template apply<rep_t>::factor;
		typename particle_t::size_type d=particle_t::model_type::colour_treatment::template apply<rep_t>::type::dimension;

		(phi->col_dim)*=rep_t::dimension;
		
		/* If the colour type decomposes in model_t, insert the index
		 * and raise the decomposed_colours variable: */
		
		if(particle_t::model_type::colour_treatment::template apply<rep_t>::decomposes)
		{
		    for(typename particle_t::size_type i=0;i<n;++i)
		    {
			phi->colour_index_ranges.insert(phi->colour_index_ranges.begin(),d);
			phi->index_ranges.insert(phi->index_ranges.begin(),d);
			++(phi->decomposed_colours);
		    }
		}

		/* If not, just insert the index: */

		else
		{
		    for(typename particle_t::size_type i=0;i<n;++i)
		    {
			phi->colour_index_ranges.insert(phi->colour_index_ranges.begin(),d);
			phi->index_ranges.insert(phi->index_ranges.begin(),d);
			(phi->block_size)*=d;
		    }
		}

		/* Let the colour treatment determine the colour/anti-colour
		 * assignment of the representation: */

		particle_t::model_type::colour_treatment::template apply<rep_t>::fill_colour_numbers(phi->colour_numbers);
		
		/* Determine the corresponding swapped colours; by default,
		 * rep_t doesn't induce colour swapping so all swapped indices
		 * move n places (remember the painting corresponds to inserting
		 * indices in the front): */
		
		for(typename particle_t::size_type i=0;i<phi->swapped_colours.size();++i)
		{
		    phi->swapped_colours[i]+=n;
		    phi->swapped_anti_colours[i]+=n;
		}

		/* For colour doublets with opposite colour charge, a swapping
		 * procedure is inserted: */

		if(n==2 and phi->colour_numbers[0]==1 and phi->colour_numbers[1]==-1)
		{
		    phi->swapped_colours.insert(phi->swapped_colours.begin(),0);
		    phi->swapped_anti_colours.insert(phi->swapped_anti_colours.begin(),1);
		}
	    }
    };
}

#endif /*CAMGEN_PART_PAINTER_H_*/

