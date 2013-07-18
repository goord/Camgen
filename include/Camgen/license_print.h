//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_LICENSE_PRINT_H_
#define CAMGEN_LICENSE_PRINT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the license_print class, whose single purpose is to print a  *
 * block of text to the standard output stream containing the author name and  *
 * copyright/license information.                                              *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class license_print
    {
	public:

	    /* Static initialiser method, printing the information: */

	    static bool initialise();
	    
	    /* Switch to disable printing: */
	    
	    static void disable();

	    /* Switch to enable printing: */

	    static void enable();
	
	private:
	    
	    /* Static initialisation flag, so that a program only prints once: */
	    
	    static bool initialised;

	    /* External control switch: */

	    static bool enabled;
	    
	    /* Private constructors and assignment operator: */
	    
	    license_print();
	    license_print(const license_print&);
	    license_print& operator = (const license_print&);
    };
}

#endif /*CAMGEN_LICENSE_PRINT_H_*/

