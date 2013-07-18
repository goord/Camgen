//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SU_N_H_
#define CAMGEN_SU_N_H_

#include <Camgen/group.h>
#include <Camgen/fundam_rep.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the SU(N) fundamental representation class and the SU(N) group class. The *
 * latter's structure constants are derived from the fundamental representation by         *
 *                                                                                         *
 * 				f_{abc}=-2*i*Tr([T_a,T_b]T_c)                             *
 *                                                                                         *
 * where the T_a denote the generators in the fundamental representation. These are        *
 * constructed as usual: starting in the upper-left corner, embed the SU(2) matrices, then *
 * move to the lower right, adding each time all possible hermitian offdiagonal matrices   *
 * and the traceless diagonal generator.                                                   *
 *                                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{   
    /* Declaration of the group class: */

    template<std::size_t N>class SU;

    /* Decaration and definition of the fundamental representation class: */
    
    template<std::size_t N>class fundamental_rep< SU<N> >
    {
	public:
	    
	    /* Group wrapper class type definion: */

	    typedef SU<N> group_type;
	    
	    /* Dimension of the representation: */

	    static const std::size_t dimension=N;

	    /* Defining the multiplicity w.r.t. representation composition: */

	    static const std::size_t multiplicity=1;

	    /* Implementation class template: */

	    template<class value_t>class implementation: public representation<value_t,implementation<value_t>,N,N*N-1>
	    {
		public:

		    /* Type reference to the base class template: */
		    
		    typedef representation<value_t,implementation<value_t>,N,N*N-1> base_type;
		    
		    /* Group implementation class type definition: */

		    typedef typename SU<N>::template implementation<value_t> group_implementation;
		    
		    /* Hermiticity tag set true: */
		    
		    static const bool hermitian=true;

		    /* Instead of providing a fill_generators() function, the
		     * initialisation function is overloaded (omitting the
		     * checks usually performed): */

		    static void initialise()
		    {
			if(!initialised)
			{

			    std::size_t n=0;
			    for(std::size_t i=1;i<N;++i)
			    {
				/* Constructing the off-diagonal generators: */

				for(std::size_t j=0;j<i;++j)
				{
				    base_type::T[n][i][j]=std::complex<value_t>(0.5,0);
				    base_type::T[n][j][i]=std::complex<value_t>(0.5,0);
				    base_type::T[n+1][i][j]=std::complex<value_t>(0,0.5);
				    base_type::T[n+1][j][i]=-std::complex<value_t>(0,0.5);
				    n+=2;
				}

				/* Constructing the diagonal generator: */

				std::complex<value_t>factor=std::complex<value_t>(1,0)/std::sqrt((std::complex<value_t>)(2*i*(i+1)));
				for(std::size_t k=0;k<i;++k)
				{
				    base_type::T[n][k][k]=factor;
				}
				base_type::T[n][i][i]=-((std::complex<value_t>)i)*factor;
				++n;
			    }
			    initialised=true;

			}
		    }
		private:
		    static bool initialised;
	    };
    };
    template<std::size_t N>const std::size_t fundamental_rep< SU<N> >::dimension;
    template<std::size_t N>template<class value_t>const bool fundamental_rep< SU<N> >::implementation<value_t>::hermitian;
    template<std::size_t N>template<class value_t>bool fundamental_rep< SU<N> >::implementation<value_t>::initialised=false;

    /* Definition of the SU(N) group class: */

    template<std::size_t N>class SU
    {
	public:

	    /* Rank of the group: */

	    static const std::size_t rank=N*N-1;

	    /* Implementation class template: */

	    template<class value_t>class implementation: public group<value_t,implementation<value_t>,N*N-1>
	    {
		private:

		    /* Helper class to determine the structure constants: the
		     * fundamental representation: */

		    typedef typename fundamental_rep< SU<N> >::template implementation<value_t> helper_rep;

		public:

		    /* Type reference to the base class template: */

		    typedef group<value_t,implementation<value_t>,N*N-1> base_type;

		    /* Antisymmetry tag is set true: */

		    static const bool antisymmetric=true;

		    /* Implementation of the structure constants: */

		    static void fill_structure_constants()
		    {

			helper_rep::initialise();

			/* The structure constants are implemented as f_{abc} =
			 * -2*i*Tr([T_a,T_b]T_b): */

			for(std::size_t a=0;a<N*N-1;++a)
			{
			    for(std::size_t b=a+1;b<N*N-1;++b)
			    {
				for(std::size_t c=b+1;c<N*N-1;++c)
				{
				    std::complex<value_t>x=0;
				    for(std::size_t i=0;i<N;++i)
				    {
					for(std::size_t j=0;j<N;++j)
					{
					    for(std::size_t k=0;k<N;++k)
					    {
						x+=(helper_rep::generator(a,i,j)*helper_rep::generator(b,j,k)*helper_rep::generator(c,k,i));
						x-=(helper_rep::generator(b,i,j)*helper_rep::generator(a,j,k)*helper_rep::generator(c,k,i));
					    }
					}
				    }
				    base_type::f[a][b][c]=-std::complex<value_t>(0,2)*x;
				}
			    }
			}

		    }
	    };
    };
    template<std::size_t N>const std::size_t SU<N>::rank;
    template<std::size_t N>template<class value_t>const bool SU<N>::implementation<value_t>::antisymmetric;
}

#endif /*CAMGEN_SU_N_H_*/

