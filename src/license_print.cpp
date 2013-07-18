//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/license_print.h>
#include <iostream>
#include <config.h>

namespace Camgen
{
    bool license_print::initialised=false;
    bool license_print::enabled=true;

    license_print::license_print(){}

    license_print::license_print(const license_print& obj){}

    license_print& license_print::operator=(const license_print& obj)
    {
	return *this;
    }

    bool license_print::initialise()
    {
	if(enabled and !initialised)
	{
	    std::cout<<std::endl;
	    std::cout<<"================================================================"<<std::endl<<std::endl;
	    std::cout<<"                 This is "<<PACKAGE_NAME<<", version "<<PACKAGE_VERSION<<std::endl;
	    std::cout<<"              Copyright (C) 2010 Gijs van den Oord,"<<std::endl;
	    std::cout<<"              licensed under the GNU GPL version 2                "<<std::endl;
	    std::cout<<" Send your bug reports to "<<PACKAGE_BUGREPORT<<std::endl<<std::endl;
	    std::cout<<"================================================================"<<std::endl<<std::endl;
	    initialised=true;
	}
	return initialised;
    }

    void license_print::disable()
    {
	enabled=false;
    }
    
    void license_print::enable()
    {
	enabled=true;
    }
}

