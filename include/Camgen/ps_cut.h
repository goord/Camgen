//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_cut.h
    \brief Abstract interface base class for phase space cut objects.
 */

#ifndef CAMGEN_PS_CUT_H_
#define CAMGEN_PS_CUT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for phase space cut objects. Derived classes should *
 * implement the pass() method.                                            *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

namespace Camgen
{
    //TODO: Refactor, pass needs to get an ps-generator argument, ps generators should not be cuts, single
    // responsability!
    /// Abstract base class for phase space cut objects.

    class phase_space_cut
    {
	public:

	    /// Constructor.

	    phase_space_cut();

	    /// Destructor.

	    virtual ~phase_space_cut();

	    /// Abstract method determining whether cuts are passed.

	    virtual bool pass()=0;

	    /// Virtual printing method.

	    virtual std::ostream& print(std::ostream&) const;
    };

    /* Phase space cut class negating a cut. */

    class phase_space_cut_not: public phase_space_cut
    {
	public:

	    /* Argument cut: */

	    phase_space_cut& cut;

	    /* Constructor: */
	    
	    phase_space_cut_not(phase_space_cut&);

	    /* Destructor: */

	    ~phase_space_cut_not();

	    /// Implementation of the cut passing method.

	    bool pass();

	    /// Virtual printing method.

	    std::ostream& print(std::ostream&) const;
    };

    /// Operator producing a cut that negates the argument cut.

    phase_space_cut_not operator !(phase_space_cut&);

    /* Phase space cut class merging 2 cuts with the 'and' operator. */

    class phase_space_cut_and: public phase_space_cut
    {
	public:
	    
	    /* Argument cuts: */
	    
	    phase_space_cut& first;
	    phase_space_cut& second;

	    /* Constructor: */
	    
	    phase_space_cut_and(phase_space_cut&,phase_space_cut&);

	    /* Destructor: */

	    ~phase_space_cut_and();

	    /// Implementation of the cut passing method.

	    bool pass();

	    /// Virtual printing method.

	    std::ostream& print(std::ostream&) const;
    };

    /// Operator producing a cut that implements the logical 'and' of the arguments.

    phase_space_cut_and operator && (phase_space_cut&,phase_space_cut&);

    /* Phase space cut class merging 2 cuts with the 'or' operator. */

    class phase_space_cut_or: public phase_space_cut
    {
	public:

	    /* Argument cuts: */

	    phase_space_cut& first;
	    phase_space_cut& second;

	    /* Constructor: */
	    
	    phase_space_cut_or(phase_space_cut&,phase_space_cut&);

	    /* Destructor: */

	    ~phase_space_cut_or();

	    /// Implementation of the cut passing method.

	    bool pass();

	    /// Virtual printing method.

	    std::ostream& print(std::ostream&) const;
    };

    /// Operator producing a cut that implements the logical 'or' of the arguments.

    phase_space_cut_or operator || (phase_space_cut&,phase_space_cut&);
}

#endif /*CAMGEN_PS_CUT_H_*/


