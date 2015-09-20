//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file parni_int.h
    \brief Monte-Carlo integrator class using an adaptive, unequal-spaced grid.
 */

#ifndef CAMGEN_PARNI_INT_H_
#define CAMGEN_PARNI_INT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Parni is a Monte-Carlo integrating class, which uses an adaptive grid with  *
 * unequal spacings. When calling the update method, the bin weights are       *
 * updated, and when the adapt method is called, the highest-weight bin is     *
 * split, approximating the integrand optimally.                               *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/parni_gen.h>
#include <Camgen/MC_int.h>

namespace Camgen
{
    /// Parni adaptive grid for multi-dimensional integration.

    template<class value_t,std::size_t D,class rng_t,class key_t>class parni_integrator: public object_allocator<typename parni_bin<value_t,D,rng_t,key_t>::point_type>, public MC_integrator<value_t>
    {
	/* Type definitions: */

	typedef object_allocator<typename parni_bin<value_t,D,rng_t>::point_type> base_type;

	public:

	    typedef value_t value_type;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef std::size_t size_type;
	    typedef parni_bin<value_t,D,rng_t,key_t> bin_type;
	    typedef typename bin_type::point_type point_type;
	    typedef typename bin_type::bin_iterator bin_iterator;
	    typedef typename bin_type::const_bin_iterator const_bin_iterator;
	    typedef typename bin_type::reverse_bin_iterator reverse_bin_iterator;
	    typedef typename bin_type::const_reverse_bin_iterator const_reverse_bin_iterator;
	    typedef typename bin_type::leaf_iterator leaf_iterator;
	    typedef typename bin_type::const_leaf_iterator const_leaf_iterator;
	    typedef typename bin_type::reverse_leaf_iterator reverse_leaf_iterator;
	    typedef typename bin_type::const_reverse_leaf_iterator const_reverse_leaf_iterator;

	    /* Public data members: */
	    /*----------------------*/

	    /// Maximal number of bins.

	    const size_type n_bins;

	    /// Splitting criterium
	    
	    const grid_modes::type mode;

	    /// Constructor with lower-left and upper-right integration domain
	    /// corners and the maximal number of bins as arguments.

	    /* Public constructors: */
	    /*----------------------*/
	    
	    parni_integrator(const point_type& xmin_,const point_type& xmax_,size_type n_bins_=500,grid_modes::type mode_=grid_modes::cumulant_weights):base_type(1),n_bins(n_bins_),mode(mode_),generator(&(this->object(0),xmin_,xmax_,n_bins_,mode_)){}

	    /// Constructor with an external point address (to be filled) and
	    /// lower-left and upper-right integration domain corners and the
	    /// maximal number of bins as arguments.
	    
	    parni_integrator(point_type* point_,const point_type& xmin_,const point_type& xmax_,size_type n_bins_=500,grid_modes::type mode_=grid_modes::cumulant_weights):base_type(point_),n_bins(n_bins_),mode(mode_),generator(point_,xmin_,xmax_,n_bins_,mode_){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Generation method implementation.

	    bool generate()
	    {
		bool result=generator.generate();
		this->weight()=generator.weight();
		return result;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		bool result=generator.evaluate_weight();
		this->weight()=generator.weight();
		return result;
	    }

	    /// Overridden updating method.

	    void update()
	    {
		generator.integrand()=this->integrand();
		generator.update();
	    }

	    /// Overridden adaptation method.

	    void adapt()
	    {
		generator.adapt();
	    }

	    /// Resets the cross section and re-initialises the grid to a single
	    /// bin.

	    void reset()
	    {
		generator.reset();
	    }
	    
	    /// Creates a new sub grid generator with the current instance as parent.

	    parni_sub_grid<value_t,D,rng_t,key_t>* create_sub_grid(const point_type& xmin,const point_type& xmax)
	    {
		return generator.create_sub_grid(xmin,xmax);
	    }

	    /* Public readout functions: */
	    /*---------------------------*/

	    /// Returns whether the maximal number of bins has been reached.

	    bool final() const
	    {
		return generator.final();
	    }

	    /// Returns the number of bins in the tree.

	    size_type count_bins() const
	    {
		return generator.count_bins();
	    }

	    /// Returns the number of leaves in the tree.

	    size_type count_leaves() const
	    {
		return generator.count_leaves();
	    }

	    /// Weight normalisation.

	    const value_type& norm() const
	    {
		return generator.norm();
	    }

	    /// Printing method.

	    std::ostream& print(std::ostream& os) const
	    {
		return generator.print(os);
	    }

	    /// Method printing only child bins.

	    std::ostream& print_leaves(std::ostream& os) const
	    {
		return generator.print_leaves(os);
	    }

	    /// Adds the rectangles to the plotscript pointer

	    plot_script* print_rectangles(plot_script* ps) const
	    {
		return generator.print_rectangles(ps);
	    }

	private:

	    /* Engine generator: */

	    parni_generator<value_t,D,rng_t,key_t> generator;
    };
}

#endif /*CAMGEN_PARNI_INT_H_*/

