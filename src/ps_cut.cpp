//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/ps_cut.h>

namespace Camgen
{
    /* Phase space cut constructor: */

    phase_space_cut::phase_space_cut(){}

    /* Phase space cut destructor: */

    phase_space_cut::~phase_space_cut(){}

    /* Printing method: */

    std::ostream& phase_space_cut::print(std::ostream& os) const
    {
	return os;
    }
    

    /* Phase space cut negation constructor: */

    phase_space_cut_not::phase_space_cut_not(phase_space_cut& cut_):cut(cut_){}

    /* Phase space cut negation destructor: */

    phase_space_cut_not::~phase_space_cut_not(){}

    /* Phase space cut negation method implementation: */

    bool phase_space_cut_not::pass()
    {
	return !cut.pass();
    }

    /* Printing method: */

    std::ostream& phase_space_cut_not::print(std::ostream& os) const
    {
	os<<"!(";
	cut.print(os);
	os<<')';
	return os;
    }

    /* Phase space cut negation unary operator: */

    phase_space_cut_not operator ! (phase_space_cut& cut_)
    {
	return phase_space_cut_not(cut_);
    }


    /* Phase space cut logical 'and' constructor: */

    phase_space_cut_and::phase_space_cut_and(phase_space_cut& first_,phase_space_cut& second_):first(first_),second(second_){}

    /* Phase space cut logical 'and' destructor: */

    phase_space_cut_and::~phase_space_cut_and(){}

    /* Phase space cut logical 'and' method implementation: */

    bool phase_space_cut_and::pass()
    {
	return (first.pass() and second.pass());
    }

    /* Printing method: */

    std::ostream& phase_space_cut_and::print(std::ostream& os) const
    {
	first.print(os);
	second.print(os);
	return os;
    }

    /* Phase space cut logical 'and' binary operator: */

    phase_space_cut_and operator && (phase_space_cut& first_,phase_space_cut& second_)
    {
	return phase_space_cut_and(first_,second_);
    }


    /* Phase space cut logical 'or' constructor: */

    phase_space_cut_or::phase_space_cut_or(phase_space_cut& first_,phase_space_cut& second_):first(first_),second(second_){}

    /* Phase space cut logical 'and' destructor: */

    phase_space_cut_or::~phase_space_cut_or(){}

    /* Phase space cut logical 'and' method implementation: */

    bool phase_space_cut_or::pass()
    {
	return (first.pass() or second.pass());
    }

    /* Printing method: */

    std::ostream& phase_space_cut_or::print(std::ostream& os) const
    {
	os<<'(';
	first.print(os);
	os<<")\n or \n(";
	second.print(os);
	os<<')';
	return os;
    }

    /* Phase space cut logical 'and' binary operator: */

    phase_space_cut_or operator || (phase_space_cut& first_,phase_space_cut& second_)
    {
	return phase_space_cut_or(first_,second_);
    }
}
