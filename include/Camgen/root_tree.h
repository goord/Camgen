//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file root_tree.h
    \brief Root TTree wrapper.
 */

#ifndef CAMGEN_ROOT_TREE_H_
#define CAMGEN_ROOT_TREE_H_

#include <string>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Root TTree wrapper. This is the correct way of doing things because only the source file will depend upon the *
 * configuration header. We therefore use void-pointers to the root-specific classes.                            *
 *                                                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// ROOT TTree wrapper class.

    class root_tree
    {
        public:

            typedef std::size_t size_type;

            /// Boolean type identifyer for tree branches.

            static const char* id_bool;

            /// Integer type identifyer for tree branches.
            
            static const char* id_int;

            /// Long int type identifyer for tree branches.
            
            static const char* id_long;

            /// Float type identifyer for tree branches.
            
            static const char* id_double;

            /// Double-precision type identifyer for tree branches.
            
            static const char* id_float;

            /// Constructor. Empty.
            
            root_tree();

            /// Destructor. Writes buffer and closes file if necessary.

            ~root_tree();

            /// Opens a root file with name first argument and a TTree with name second argument.

            bool open(const std::string&,const std::string&);

            /// Opens a root file with name first argument, a TTree with name second argument and description third
            /// argument.

            bool open(const std::string&,const std::string&,const std::string&);

            /// Writes buffer and closes file.
            
            void close();

            /// Returns whether the file is open.
            
            bool is_open() const;

            /// Creates a branch in the tree containing booleans, named last argument.

            bool branch(const bool*,const std::string&);

            /// Creates a branch in the tree containing boolean arrays of length second argument, named last argument.

            bool branch(const bool*,size_type,const std::string&);

            /// Creates a branch in the tree containing integers, named last argument.
            
            bool branch(const int*,const std::string&);

            /// Creates a branch in the tree containing integer arrays of length second argument, named last argument.
            
            bool branch(const int*,size_type,const std::string&);

            /// Creates a branch in the tree containing long integers, named last argument.
            
            bool branch(const long int*,const std::string&);

            /// Creates a branch in the tree containing long integer arrays of length second argument, named last argument.
            
            bool branch(const long int*,size_type,const std::string&);

            /// Creates a branch in the tree containing floats, named last argument.
            
            bool branch(const float*,const std::string&);

            /// Creates a branch in the tree containing float arrays of length second argument, named last argument.
            
            bool branch(const float*,size_type,const std::string&);

            /// Creates a branch in the tree containing doubles, named last argument.
            
            bool branch(const double*,const std::string&);

            /// Creates a branch in the tree containing double arrays of length second argument, named last argument.
            
            bool branch(const double*,size_type,const std::string&);

            /// Streams all pointer values to the event in the tree.

            bool fill();

        private:

            /* Helper function for opening the file: */

            bool open(const char*,const char*,const char*);

            /* ROOT TTree, or NULL: */

            void* roottree;

            /* ROOT TFile, or NULL: */

            void* rootfile;
    };
}

#endif /*CAMGEN_ROOT_TREE_H_*/
