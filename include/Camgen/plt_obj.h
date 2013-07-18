//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file pltobj.h
    \brief Plot object definitions for plot-scripting class.
 */

#ifndef CAMGEN_PLT_OBJ_H_
#define CAMGEN_PLT_OBJ_H_

#include <string>
#include <iostream>

namespace Camgen
{
    /// Abstract base type of the plot objects classes.

    class plot_object
    {
	public:

	    /// RGB fill colour specification.

	    std::string fill_colour;

	    /// Fill style specification.
	    
	    std::string fill_style;

	    /// Line type specification.
	    
	    int line_type;

	    /// Line width specification.

	    double line_width;

	    /// RGB line colour specification.

	    std::string line_colour;

	    /// No border flag.

	    bool no_border;

	    /// Trivial constructor.

	    plot_object();
	    
	    /// Trivial destructor.
	    
	    virtual ~plot_object();

	    /// Pure abstract method outputting the specific object.

	    virtual void write_shape(std::ostream&) const=0;

	    /// Clone method.
	    
	    virtual plot_object* clone() const=0;

	    /// Sets the object to the front of the rest.

	    void set_front();

	    /// Sets the object to the back of the rest.

	    void set_back();

	    /// Sets the object as backgound to the rest.

	    void set_behind();

	    /// Sets the fill pattern

	    void set_pattern(int);

	    /// Plot command writing method.

	    void write(std::ostream&,int) const;

	private:

	    std::string visibility;
	    int pattern;
    };

    /// Object plot class for circles.

    class plot_circle: public plot_object
    {
	public:

	    /// Midpoint positions and radius.

	    double x,y,r;

	    /// Constructor.

	    plot_circle(const double&,const double&,const double&);

	    /// Shape writing method implementation.

	    void write_shape(std::ostream&) const;

	    /// Clone implementation.
	    
	    plot_circle* clone() const;
    };

    /// Object plot class for rectangles.

    class plot_rectangle: public plot_object
    {
	public:

	    /// Lower left point and upper right point.

	    double xmin,ymin,xmax,ymax;

	    /// Constructor.

	    plot_rectangle(const double&,const double&,const double&,const double&);

	    /// Shape writing method implementation.

	    void write_shape(std::ostream&) const;

	    /// Clone implementation.
	    
	    plot_rectangle* clone() const;
    };

    /// Object plot class for ellipses.

    class plot_ellipsis: public plot_object
    {
	public:

	    /// Midpoint position and radii.

	    double x,y,r1,r2;

	    /// Constructor.

	    plot_ellipsis(const double&,const double&,const double&,const double&);

	    /// Shape writing method implementation.

	    void write_shape(std::ostream&) const;

	    /// Clone implementation.
	    
	    plot_ellipsis* clone() const;
    };
}

#endif /*CAMGEN_PLT_OBJ_H_*/

