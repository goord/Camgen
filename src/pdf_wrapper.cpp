//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <cmath>
#include <limits>
#include <config.h>
#include <Camgen/pdf_wrapper.h>

#if HAVE_LHAPDF_H_

#include <LHAPDF/LHAPDF.h>

namespace Camgen
{
    bool pdf_wrapper::init=false;
    std::string pdf_wrapper::setname;
    int pdf_wrapper::setnr=0;

    void pdf_wrapper::initialise(const char* setname_, int setnr_)
    {
	std::string namestr(setname_);
	if(setname!=namestr or setnr!=setnr_)
	{
	    setname=namestr;
	    setnr=setnr_;
	    LHAPDF::initPDFSet(setname.c_str());
	    LHAPDF::initPDF(setnr);
	}
	init=true;
    }

    void pdf_wrapper::reset()
    {
	setname.clear();
	setnr=0;
	init=false;
    }

    double pdf_wrapper::xmin()
    {
	return init?LHAPDF::getXmin(setnr):(double)0;
    }

    double pdf_wrapper::xmax()
    {
	return init?LHAPDF::getXmax(setnr):(double)1;
    }

    double pdf_wrapper::ymin()
    {
	return -ymax();
    }

    double pdf_wrapper::ymax()
    {
	return init?0.5*std::log(LHAPDF::getXmax(setnr)/LHAPDF::getXmin(setnr)):std::numeric_limits<double>::infinity();
    }

    std::string pdf_wrapper::name()
    {
	return init?setname:"pdf set not initialised";
    }

    int pdf_wrapper::number()
    {
	return setnr;
    }

    int pdf_wrapper::group()
    {
#ifdef __APPLE__
	return 0;
#else
	return init?LHAPDF::getPDFSetInfo(setname.c_str(),setnr).pdflibNGroup:0;
#endif
    }

    int pdf_wrapper::set()
    {
#ifdef __APPLE__
	return 0;
#else
	return init?LHAPDF::getPDFSetInfo(setname.c_str(),setnr).pdflibNSet:0;
#endif
    }

    double pdf_wrapper::xf(double x,int q,double mu)
    {
	return init?LHAPDF::xfx(x,mu,q):x;
    }

    double pdf_wrapper::f(double x,int q,double mu)
    {
	return init?(LHAPDF::xfx(x,mu,q)/x):(double)1;
    }

    double pdf_wrapper::ff(double x1,double x2,int q1,int q2,double mu)
    {
	if(init)
	{
	    double xfx1=LHAPDF::xfx(x1,mu,q1);
	    double xfx2=LHAPDF::xfx(x2,mu,q2);
	    double symm=(q1==q2)?double(0.5):double(1);
	    return (xfx1>(double)0 and xfx2>(double)0)?(symm*xfx1*xfx2/(x1*x2)):(double)0;
	}
	return (double)0.5;
    }

    double pdf_wrapper::alpha_s(double mu)
    {
	return init?LHAPDF::alphasPDF(mu):double(-1);
    }
}

#else

namespace Camgen
{
    bool pdf_wrapper::init=true;
    std::string pdf_wrapper::setname;
    int pdf_wrapper::setnr=0;

    void pdf_wrapper::initialise(const char* setname_, int setnr_){}

    void pdf_wrapper::reset(){}

    double pdf_wrapper::xmin()
    {
	return (double)0;
    }

    double pdf_wrapper::xmax()
    {
	return (double)1;
    }

    double pdf_wrapper::ymin()
    {
	return -std::numeric_limits<double>::infinity();
    }

    double pdf_wrapper::ymax()
    {
	return std::numeric_limits<double>::infinity();
    }

    std::string pdf_wrapper::name()
    {
	return "pdf set not initialised";
    }

    int pdf_wrapper::number()
    {
	return 0;
    }

    int pdf_wrapper::group()
    {
	return 0;
    }

    int pdf_wrapper::set()
    {
	return 0;
    }

    double pdf_wrapper::xf(double x,int q,double mu)
    {
	return x;
    }

    double pdf_wrapper::f(double x,int q,double mu)
    {
	return (double)1;
    }

    double pdf_wrapper::ff(double x1,double x2,int q1,int q2,double mu)
    {
	return (double)1;
    }

    double pdf_wrapper::alpha_s(double mu)
    {
	return double(-1);
    }
}

#endif /*HAVE_LHAPDF_H_*/

