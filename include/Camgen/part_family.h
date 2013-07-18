//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PART_FAMILY_H_
#define CAMGEN_PART_FAMILY_H_

#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <Camgen/forward_decs.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Jet definition class. It is derived from the vector of particle pointer type, *
 * with some functionality and a name added.                                     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Definition and declaration of the class template: */

    template<class model_t>class particle_family: public std::vector<const particle<model_t>*>
    {
	private:
	    
	    /* Base type declaration: */

	    typedef std::vector<const particle<model_t>*> base_type;
	
	public:
	    
	    /* STL-container compliant type definitions, inherited from the base type: */
	    
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::reference reference;
	    typedef typename base_type::const_reference const_reference;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;
	    typedef typename base_type::reverse_iterator reverse_iterator;
	    typedef typename base_type::const_reverse_iterator const_reverse_iterator;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::difference_type difference_type;
	    typedef typename base_type::allocator_type allocator_type;
	    typedef typename base_type::pointer pointer;
	    typedef typename base_type::const_pointer const_pointer;

	    /* Trivial constructor: */

	    particle_family(){}
	    
	    /* Named constructor: */
	    
	    particle_family(const std::string& str):name(str){}

	    /* Jet name readout: */

	    std::string get_name() const
	    {
		return name;
	    }

	    /* Particle insertion function: */

	    particle_family<model_t>& add_particle(const particle<model_t>* phi)
	    {
		if(phi!=NULL and !std::binary_search(this->begin(),this->end(),phi))
		{
		    this->push_back(phi);
		}
		return *this;
	    }

	    /* Particle removal function: */

	    particle_family<model_t>& remove_particle(const particle<model_t>* phi)
	    {
		std::remove(this->begin(),this->end(),phi);
		return *this;
	    }

	private:

	    /* Name of the family object: */

	    std::string name;
    };

    /* Outstream operator for particle families: */

    template<class model_t>std::ostream& operator << (std::ostream& os,const particle_family<model_t>& j)
    {
	os<<j.get_name()<<"=";
	for(std::size_t i=0;i<j.size()-1;++i)
	{
	    os<<j[i]->get_name()<<",";
	}
	os<<j.back()->get_name()<<std::endl;
	return os;
    }
}



#endif /*CAMGEN_PART_FAMILY_H_*/

