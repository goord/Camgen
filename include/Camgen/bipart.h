//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_BI_PART_H_
#define CAMGEN_BI_PART_H_

#include <vector>
#include <set>
#include <utility>
#include <iostream>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Partitioning algorithm for contributing invariant mass channels into pairs    *
 * for determining the minimal invariant mass of a timelike phase-space channel. *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class bi_partition: public std::vector< std::vector< std::pair<std::size_t,std::size_t> > >
    {
	public:

	    typedef std::size_t index_type;
	    typedef std::size_t size_type;

	    /// Constructor.

	    bi_partition(const std::set<index_type>&);

	private:
	    
	    /* Partition construction iteration: */

	    void add_partition(const std::set<size_type>&);

	    /* Partition construction iteration helper: */

	    bool next_index_set(size_type);

	    size_type n;
	    std::vector<size_type>kappa;
	    std::vector<size_type>maxima;
    };

    /* Partition printing facility: */

    std::ostream& operator << (std::ostream&,const bi_partition&);
}

#endif /*CAMGEN_BI_PART_H_*/

