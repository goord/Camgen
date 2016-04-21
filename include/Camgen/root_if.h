//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file root_if.h
    \brief Root TTree output interface for phase space generators.
 */

#ifndef CAMGEN_ROOT_IF_H_
#define CAMGEN_ROOT_IF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Root TTree output interface class implementation. It creates a root file and  *
 * implements the branching members as creators of branchings in the TTree.      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/if_output.h>
#include <Camgen/root_tree.h>

namespace Camgen
{
    
    /// ROOT TTree output interface class.
    
    template<class model_t>class root_interface: public interface_output<model_t>
    {
	typedef interface_output<model_t> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;

	    /// ROOT tree name.

	    const std::string tree_name;

	    /// Constructor with file name, tree name and tree description arguments.

	    root_interface(const std::string& file_name_="output",
		           const std::string tree_name_="output_tree",
		           const std::string description_="no description available"):base_type(file_name_,description_),tree_name(tree_name_){}

	    /* Factory method implementation: */

	    base_type* create(const std::string& file_name_) const
	    {
		return new root_interface<model_t>(file_name_,tree_name,this->description);
	    }

	    /* Opens the root file. */

	    bool open_file()
	    {
                return tree.open(this->file_name,tree_name,this->description);
	    }

	    /* Closes the .root file. */

	    bool close_file()
	    {
                return tree.close();
	    }

	    /* Inserts a branch holding a momentum: */

	    bool branch(const momentum_type* p,const std::string& name)
	    {
                return tree.branch(p->data,p->size(),name);
	    }

	    /* Inserts a branch holding a floating-point number: */

	    bool branch(const value_type* x,const std::string& name)
	    {
                return tree.branch(x,name);
	    }

	    /* Inserts a branch holding an integer: */

	    bool branch(const int* n,const std::string& name)
	    {
                return tree.branch(n,name);
	    }

	    /* Inserts a branch holding a boolean: */

	    bool branch(const bool* n,const std::string& name)
	    {
                return tree.branch(n,name);
	    }

	    /* Writes the event to the ROOT tree. */

	    bool write_event()
	    {
                return tree.fill();
	    }

	private:

            root_tree tree;
    };
}

#endif /*CAMGEN_ROOT_TTREE*/

