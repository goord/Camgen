//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_COMP_CONTR_H_
#define CAMGEN_COMP_CONTR_H_

#include <vector>
#include <Camgen/debug.h>
#include <Camgen/eval.h>
#include <Camgen/d.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition and declaration of the compose_contraction class                 *
 * template.Determines the contraction of the final off-shell current in the   *
 * tree with the corresponding external wave function in the case of coloured  *
 * particles. In particular in the colour-flow treatment of gluons, this       *
 * requires a subtraction of the external trace.                               *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Declaration and definition of the contraction composer: the first template
     * parameter is the outer colour delta function, the second parameter the inner colour
     * structure or spacetime index contraction: */

    template<class contr_t1,class contr_t2>class compose_contraction
    {
	public:
	    /* Fermionic behaviour determined by inner contraction: */

	    static const bool fermionic=contr_t2::fermionic;
	    
	    /* Size of the total Kronecker delta-tensor is the product: */
	    
	    static const std::size_t tensor_size=contr_t1::tensor_size*contr_t2::tensor_size;

	    /* Size of the tensors to be contracted: */

	    static const std::size_t size;
    };
    template<class contr_t1,class contr_t2>const bool compose_contraction<contr_t1,contr_t2>::fermionic;
    template<class contr_t1,class contr_t2>const std::size_t compose_contraction<contr_t1,contr_t2>::tensor_size;
    template<class contr_t1,class contr_t2>const std::size_t compose_contraction<contr_t1,contr_t2>::size=contr_t1::ranges[0]*contr_t2::size;

    /* Overloaded evaluate class template: */

    template<class rep_t,class contr_t>class evaluate<compose_contraction<colour_tensor::d<rep_t,0,1>,contr_t> >
    {
	public:

	    /* Reference typedef to the contraction class: */

	    typedef compose_contraction< colour_tensor::d<rep_t,0,1>,contr_t> contraction_type;
	    
	    /* Standard typedefs determined by the inner contraction: */
	    
	    typedef typename evaluate<contr_t>::model_type model_type;
	    typedef typename evaluate<contr_t>::value_type value_type;
	    typedef typename evaluate<contr_t>::r_value_type r_value_type;
	    typedef typename evaluate<contr_t>::size_type size_type;
	    typedef typename evaluate<contr_t>::tensor_type tensor_type;
	    typedef typename evaluate<contr_t>::iterator iterator;
	    typedef typename evaluate<contr_t>::const_iterator const_iterator;
	    typedef typename evaluate<contr_t>::momentum_type momentum_type;

	    /* Intialisation function. Fills the stepsize with the size of the inner
	     * tensors: */
	
	    static void initialise()
	    {
		evaluate<contr_t>::initialise();
	    }

	    /* Rank-vector filler for checking purposes: */

	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		evaluate<contr_t>::fill_rank_vector(r);
		r.insert(r.begin(),rep_t::dimension);
		return r;
	    }

	    /* Implementation of the contraction for discrete colours; it2, the external
	     * amplitude iterator is supposed point to the tensor block determined by
	     * colours[col]. The implementation therefore just moves it1, the internal
	     * amplitude iterator to that position and finally performs the spacetime
	     * contraction: */

	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type col)
	    {
		CAMGEN_ERROR_IF((it1.range()<(rep_t::dimension*contr_t::size)),"iterator 1 out of range");
		CAMGEN_ERROR_IF((it2.range()<(rep_t::dimension*contr_t::size)),"iterator 2 out of range");
		it1+=colours[col]*contr_t::size;
		it2+=colours[col]*contr_t::size;
		return evaluate<contr_t>::apply(it1,it2,colours,col+1);
	    }

	    /* Implementation of the contraction for continuous colours; this function
	     * traverses both tensors by the stepsize, and each time performs the inner
	     * contraction, adding the result to the first argument: */

	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		CAMGEN_ERROR_IF((it1.range()<(rep_t::dimension*contr_t::size)),"iterator 1 out of range");
		CAMGEN_ERROR_IF((it2.range()<(rep_t::dimension*contr_t::size)),"iterator 2 out of range");
		for(std::size_t a=0;a<rep_t::dimension;++a)
		{
		    evaluate<contr_t>::apply(result,it1,it2);
		    it1+=contr_t::size;
		    it2+=contr_t::size;
		}
		it1-=(rep_t::dimension)*contr_t::size;
		it2-=(rep_t::dimension)*contr_t::size;
	    }
    };

    /* Specialisation for the case of gluons in the colour-flow treatment: */

    template<std::size_t N,class contr_t>class evaluate<compose_contraction< colour_tensor::CF_d<N,0,1>,contr_t> >
    {
	public:

	    /* Reference typedef to the contraction class: */

	    typedef compose_contraction<colour_tensor::CF_d<N,0,1>,contr_t> contraction_type;
	    
	    /* Standard typedefs determined by the inner contraction: */
	    
	    typedef typename evaluate<contr_t>::model_type model_type;
	    typedef typename evaluate<contr_t>::value_type value_type;
	    typedef typename evaluate<contr_t>::r_value_type r_value_type;
	    typedef typename evaluate<contr_t>::size_type size_type;
	    typedef typename evaluate<contr_t>::tensor_type tensor_type;
	    typedef typename evaluate<contr_t>::iterator iterator;
	    typedef typename evaluate<contr_t>::const_iterator const_iterator;
	    typedef typename evaluate<contr_t>::momentum_type momentum_type;

	    /* Intialisation function. Fills the smallstep with the size of the inner
	     * tensors and the medium- and bigstep by recursive multiplications with N: */
	
	    static void initialise()
	    {
		evaluate<contr_t>::initialise();
	    }

	    /* Rank-vector filler for checking purposes: */

	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		evaluate<contr_t>::fill_rank_vector(r);
		r.insert(r.begin(),2,N);
		return r;
	    }

	    /* Implementation of the contraction for discrete colours: */

	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		CAMGEN_ERROR_IF((it1.range()<(bigstep)),"iterator 1 out of range");
		CAMGEN_ERROR_IF((it2.range()<(bigstep)),"iterator 2 out of range");
		return value_type(2,0)*(evaluate<contr_t>::apply(it1+(colours[c]+colours[c+1]*N)*smallstep,it2+(colours[c]*N+colours[c+1])*smallstep,colours,c+2));
	    }

	    /* Implementation of the contraction in the case of continuous colours, where
	     * the external colour wave functions are assumed to be traceless, simply 
	     * returns 2*Gint_{AB}Gext_{BA}: */

	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		CAMGEN_ERROR_IF((it1.range()<(bigstep)),"iterator 1 out of range");
		CAMGEN_ERROR_IF((it2.range()<(bigstep)),"iterator 2 out of range");
		
		value_type temp(0,0);
		for(std::size_t a=0;a<N;++a)
		{
		    for(std::size_t b=0;b<N;++b)
		    {
			evaluate<contr_t>::apply(temp,it1,it2);
			it1+=mediumstep;
			it2+=smallstep;
		    }
		    it1-=bigstep;
		    it1+=smallstep;
		}
		it2-=bigstep;
		it1-=mediumstep;

		result+=(value_type(2,0)*temp);

	    }

	private:

	    /* Three useful step sizes: */
	    
	    static const std::size_t smallstep;
	    static const std::size_t mediumstep;
	    static const std::size_t bigstep;
    };
    template<std::size_t N,class contr_t>const std::size_t evaluate<compose_contraction<colour_tensor::CF_d<N,0,1>,contr_t> >::smallstep=contr_t::size;
    template<std::size_t N,class contr_t>const std::size_t evaluate<compose_contraction<colour_tensor::CF_d<N,0,1>,contr_t> >::mediumstep=N*contr_t::size;
    template<std::size_t N,class contr_t>const std::size_t evaluate<compose_contraction<colour_tensor::CF_d<N,0,1>,contr_t> >::bigstep=N*N*contr_t::size;
}

#endif /*CAMGEN_COMP_CONTR_H_*/

