//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file root_file.h
    \brief Root TTree output interface for phase space generators.
 */

#ifndef CAMGEN_ROOT_FILE_H_
#define CAMGEN_ROOT_FILE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Root TTree output interface class implementation. It creates a root file and  *
 * implements the branching members as creators of branchings in the TTree.      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/evt_output.h>
#include <Camgen/root_tree.h>

namespace Camgen
{
    
    /// ROOT TTree output interface class.
    
    template<class model_t,std::size_t N_in,std::size_t N_out>class root_file: public event_output<model_t,N_in,N_out>
    {
	typedef event_output<model_t,N_in,N_out> base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::event_type event_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;

	    /// ROOT tree name.

	    const std::string tree_name;

	    /// Constructor with file name, tree name and tree description arguments.

	    root_file(const std::string& file_name_="output",
		      const std::string tree_name_="output_tree",
		      const std::string description_="no description available"):base_type(file_name_,description_),tree_name(tree_name_){}

	    /* Factory method implementation: */

	    base_type* create(const std::string& file_name_) const
	    {
		return new root_file<model_t,N_in,N_out>(file_name_,tree_name,this->description);
	    }

	    /* Opens the root file. */

	    bool open_file()
	    {
                return tree.open(this->file_name,tree_name,this->description);
	    }

	    /* Closes the .root file. */

	    bool close_file()
	    {
                tree.close();
                return true;
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

