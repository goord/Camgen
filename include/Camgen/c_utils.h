//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_C_UTILS_H_
#define CAMGEN_C_UTILS_H_

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <Camgen/num_utils.h>

/* * * * * * * * * * * * * * * * * * * * * *
 * Utilities for complex numbers in Camgen *
 *                                         *
 * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{

    /* Optimised multiplication by the imaginary unit: */

    template<class T>std::complex<T>times_i(const std::complex<T>& z)
    {
	return std::complex<T>(-z.imag(),z.real());
    }

    /* Compute z1+i*z2: */

    template<class T>std::complex<T>make_z(const T& z1,const T& z2)
    {
	return std::complex<T>(z1,z2);
    }

    /* Compute z1+i*z2, where z1 and z2 are complex numbers: */

    template<class T>std::complex<T>make_z(const std::complex<T>& z1,const std::complex<T>& z2)
    {
	return std::complex<T>(z1.real()-z2.imag(),z1.imag()+z2.real());
    }

    /* Compute z1-i*z2: */

    template<class T>std::complex<T>make_zbar(const T& z1,const T& z2)
    {
	return std::complex<T>(z1,-z2);
    }
    /* Compute z1-i*z2, where z1 and z2 are complex numbers: */

    template<class T>std::complex<T>make_zbar(const std::complex<T>& z1,const std::complex<T>& z2)
    {
	return std::complex<T>(z1.real()+z2.imag(),z1.imag()-z2.real());
    }

    /* Output utilities for simple complex numbers appearing in algebraic
     * expressions: */

    template<class T>std::ostream& shortprint(std::ostream& os,const std::complex<T>& z)
    {
	if(z.imag()==(T)0)
	{
	    os<<z.real();
	}
	else if(z.real()==(T)0)
	{
	    if(z.imag()==(T)1)
	    {
		os<<"i";
	    }
	    else if(z.imag()==(T)-1)
	    {
		os<<"-i";
	    }
	    else
	    {
		os<<z.imag()<<"i";
	    }
	}
	else
	{
	    if(z.imag()==(T)1)
	    {
		os<<z.real()<<"+i";
	    }
	    else if(z.imag()==(T)-1)
	    {
		os<<z.real()<<"-i";
	    }
	    else if(z.imag()<(T)0)
	    {
		os<<z.real()<<"-"<<-z.imag()<<"i";
	    }
	    else
	    {
		os<<z.real()<<"+"<<z.imag()<<"i";
	    }
	}
	return os;
    }

    template<class T>std::ostream& shortprint_prefactor(std::ostream& os,const std::complex<T>& z)
    {
	if(z.imag()==(T)0)
	{
	    if(z.real()==(T)1)
	    {
		os<<" + ";
		return os;
	    }
	    else if(z.real()==(T)-1)
	    {
		os<<" - ";
		return os;
	    }
	    else if(z.real()<0)
	    {
		os<<" - "<<-z.real();
		return os;
	    }
	    else
	    {
		os<<" + "<<z.real();
		return os;
	    }
	}
	else if(z.real()==(T)0)
	{
	    if(z.imag()==(T)1)
	    {
		os<<" + i";
		return os;
	    }
	    else if(z.imag()==(T)-1)
	    {
		os<<" - i";
		return os;
	    }
	    else if(z.imag()<0)
	    {
		os<<" - "<<-z.imag()<<"i";
		return os;
	    }
	    else
	    {
		os<<" + "<<z.imag()<<"i";
	    }
	}
	else
	{
	    os<<" + (";
	    shortprint(os,z);
	    os<<")";
	    return os;
	}
    }

    template<class T>std::ostream& shortprint_first_prefactor(std::ostream& os,const std::complex<T>& z)
    {
	if(z.imag()==(T)0)
	{
	    if(z.real()==(T)1)
	    {
		return os;
	    }
	    else if(z.real()==(T)-1)
	    {
		os<<"-";
		return os;
	    }
	    else if(z.real()<0)
	    {
		os<<"-"<<-z.real();
		return os;
	    }
	    else
	    {
		os<<z.real();
		return os;
	    }
	}
	else if(z.real()==(T)0)
	{
	    if(z.imag()==(T)1)
	    {
		os<<"i";
		return os;
	    }
	    else if(z.imag()==(T)-1)
	    {
		os<<"-i";
		return os;
	    }
	    else if(z.imag()<0)
	    {
		os<<"-"<<-z.imag()<<"i";
		return os;
	    }
	    else
	    {
		os<<z.imag()<<"i";
	    }
	}
	else
	{
	    os<<"(";
	    shortprint(os,z);
	    os<<")";
	    return os;
	}
    }

    /* String-producing function for object types for which the streaming
     * operator is overloaded: */

    //TODO: move to general utils...
    template<class T>std::string make_string(const T& object)
    {
	std::stringstream ss;
	ss<<object;
	return ss.str();
    }
}

#endif /*CAMGEN_C_UTILS_H_*/

