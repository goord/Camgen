//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PDF_WRAPPER_H_
#define CAMGEN_PDF_WRAPPER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This class wraps the LHAPDF dependence, such that we do not have direct *
 * dependence on the configuration header in public headers, but only      *
 * implementation.                                                         *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <string>

namespace Camgen
{
    /// Wrapper class for LHAPDF dependence

    class pdf_wrapper
    {
	public:

	    /* Initialises pdfs: */

	    static void initialise(const char*,int);

	    /* Resets pdfs: */

	    static void reset();

	    /* Returns the minimal x-value: */

	    static double xmin();

	    /* Returns the maximal x-value: */

	    static double xmax();

	    /* Returns the minimal y-value: */

	    static double ymin();

	    /* Returns the maximal y-value: */

	    static double ymax();

	    /* Returns the pdf set name: */

	    static std::string name();

	    /* Returns the pdf set number: */

	    static int number();

	    /* Returns the pdf group name: */

	    static int group();

	    /* Returns the pdf set name: */

	    static int set();

	    /* Returns x*f(x) evaluated at momentum fraction x, for parton q at
	     * qcd scale mu: */

	    static double xf(double x,int q,double mu);

	    /* Returns the evaluated pdf at momentum fraction x, for parton q at
	     * qcd scale mu: */

	    static double f(double x,int q,double mu);

	    /* Returns the flux factor corresponding to the arguments: */

	    static double ff(double x1,double x2,int q1,int q2,double mu);

	    /* Returns the strong couplings constant: */

	    static double alpha_s(double mu);

	private:

	    /* Initialisation flag: */

	    static bool init;

	    /* PDF set name: */

	    static std::string setname;

	    /* PDF set number: */

	    static int setnr;

    };
}

#endif /*CAMGEN_PDF_WRAPPER_H_*/

