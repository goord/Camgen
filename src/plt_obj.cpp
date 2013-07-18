//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/plt_obj.h>

namespace Camgen
{
    /* Default base class constructor: */

    plot_object::plot_object():fill_style("empty"),line_type(1),line_width(1),no_border(false),visibility("front"),pattern(0){}

    /* Base class destructor (trivial): */

    plot_object::~plot_object(){}

    /* Puts the object on the foreground: */
    
    void plot_object::set_front()
    {
	visibility="front";
    }

    /* Puts the object on the background: */

    void plot_object::set_back()
    {
	visibility="back";
    }

    /* Puts the object as a background: */

    void plot_object::set_behind()
    {
	visibility="behind";
    }

    /* Sets the fill pattern: */

    void plot_object::set_pattern(int p)
    {
	pattern=p;
    }

    /* Script-writing method: */

    void plot_object::write(std::ostream& os,int n) const
    {
	os<<"set object "<<n<<" ";
	write_shape(os);
	if(fill_colour.size()!=0)
	{
	    os<<" fc rgb \""<<fill_colour<<"\"";
	}
	if(fill_style.size()!=0)
	{
	    if(fill_style=="pattern")
	    {
		os<<" fs pattern "<<pattern;
	    }
	    else
	    {
		os<<" fs "<<fill_style;
	    }
	}
	if(line_type!=1)
	{
	    os<<" lt "<<line_type;
	}
	if(line_width!=1)
	{
	    os<<" lw "<<line_width;
	}
	if(line_colour.size()!=0)
	{
	    os<<" lc rgb \""<<line_colour<<"\"";
	}
	if(no_border)
	{
	    os<<" noborder";
	}
	os<<" "<<visibility;
    }

    /* Plot_circle constructor: */

    plot_circle::plot_circle(const double& x_,const double& y_,const double& r_):x(x_),y(y_),r(r_){}

    /* Shape-writing method implementation: */

    void plot_circle::write_shape(std::ostream& os) const
    {
	os<<"circle at "<<x<<","<<y<<" size "<<r;
    }

    /* Clone method implementation: */

    plot_circle* plot_circle::clone() const
    {
	return new plot_circle(*this);
    }

    /* Plot_rectangle constructor: */

    plot_rectangle::plot_rectangle(const double& xmin_,const double& ymin_,const double& xmax_,const double& ymax_):xmin(xmin_),ymin(ymin_),xmax(xmax_),ymax(ymax_){}

    /* Shape-writing method implementation: */

    void plot_rectangle::write_shape(std::ostream& os) const
    {
	os<<"rect from "<<xmin<<","<<ymin<<" to "<<xmax<<","<<ymax;
    }

    /* Clone method implementation: */

    plot_rectangle* plot_rectangle::clone() const
    {
	return new plot_rectangle(*this);
    }

    /* Plot_ellipsis constructor: */

    plot_ellipsis::plot_ellipsis(const double& x_,const double& y_,const double& r1_,const double& r2_):x(x_),y(y_),r1(r1_),r2(r2_){}

    /* Shape-writing method implementation: */

    void plot_ellipsis::write_shape(std::ostream& os) const
    {
	os<<"ellipse at "<<x<<","<<y<<" size "<<r1<<","<<r2;
    }

    /* Clone method implementation: */

    plot_ellipsis* plot_ellipsis::clone() const
    {
	return new plot_ellipsis(*this);
    }
}

