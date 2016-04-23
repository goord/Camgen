//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file root_tree.h
    \brief Root TTree wrapper implementation.
 */

#include <Camgen/root_tree.h>
#include <config.h>

#if HAVE_ROOT_H
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#endif

namespace Camgen
{
    const char* root_tree::id_bool="/O";
    const char* root_tree::id_int=sizeof(int)==8?"/L":"/I";
    const char* root_tree::id_long=sizeof(long int)==8?"/L":"/I";
    const char* root_tree::id_float="/F";
    const char* root_tree::id_double="/D";

    root_tree::root_tree():roottree(NULL),rootfile(NULL){}

    root_tree::~root_tree()
    {
        if(is_open())
        {
            close();
        }
    }

    bool root_tree::open(const std::string& fname,const std::string& tname)
    {
        return open_core(fname.c_str(),tname.c_str(),NULL);
    }

    bool root_tree::open(const std::string& fname,const std::string& tname, const std::string& descr)
    {
        return open_core(fname.c_str(),tname.c_str(),descr.c_str());
    }

    bool root_tree::open_core(const char* fname,const char* tname,const char* descr)
    {
        if(is_open())
        {
            close();
        }
#if HAVE_ROOT_H
        std::string filename(fname);
        std::string::size_type n=filename.find_last_of('.');
        if(n!=filename.size()-5 or filename.substr(n+1,5)!="root")
        {
            filename+=".root";
        }

        rootfile=new TFile(filename.c_str(),"RECREATE");
        roottree=new TTree(tname,descr);
        return static_cast<TFile*>(rootfile)->IsOpen();
#else
        return false;
#endif
    }

    bool root_tree::is_open() const
    {
#if HAVE_ROOT_H
        if(rootfile!=NULL)
        {
            return static_cast<TFile*>(rootfile)->IsOpen();
        }
#endif
        return false;
    }

    void root_tree::close()
    {
#if HAVE_ROOT_H
        if(rootfile!=NULL)
        {
            static_cast<TFile*>(rootfile)->Write();
            static_cast<TFile*>(rootfile)->Close();
//            delete static_cast<TTree*>(roottree);
            delete static_cast<TFile*>(rootfile);
            rootfile=NULL;
            roottree=NULL;
        }
#endif
    }

    bool root_tree::branch(const bool* data,const std::string& varname)
    {
        return branch(data,1,varname);
    }

    bool root_tree::branch(const bool* data,size_type n,const std::string& varname)
    {
       if(!is_open() or n<1)
       {
           return false;
       }
#if HAVE_ROOT_H
       std::string id;
       if(n==1)
       {
           id=std::string(root_tree::id_bool);
       }
       else
       {
           std::stringstream ss;
           ss<<varname<<'['<<n<<']'<<root_tree::id_bool;
           id=ss.str();
       }
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<bool*>(data),create_var_id(varname,n,root_tree::id_bool).c_str());
           return true;
       }
#endif
       return false;
    }

    bool root_tree::branch(const int* data,const std::string& varname)
    {
        return branch(data,1,varname);
    }

    bool root_tree::branch(const int* data,size_type n,const std::string& varname)
    {
       if(!is_open() or n<1)
       {
           return false;
       }
#if HAVE_ROOT_H
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<int*>(data),create_var_id(varname,n,root_tree::id_int).c_str());
           return true;
       }
#endif
       return false;
    }

    bool root_tree::branch(const long int* data,const std::string& varname)
    {
        return branch(data,1,varname);
    }

    bool root_tree::branch(const long int* data,size_type n,const std::string& varname)
    {
       if(!is_open() or n<1)
       {
           return false;
       }
#if HAVE_ROOT_H
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<long int*>(data),create_var_id(varname,n,root_tree::id_long).c_str());
           return true;
       }
#endif
       return false;
    }

    bool root_tree::branch(const float* data,const std::string& varname)
    {
        return branch(data,1,varname);
    }

    bool root_tree::branch(const float* data,size_type n,const std::string& varname)
    {
       if(!is_open() or n<1)
       {
           return false;
       }
#if HAVE_ROOT_H
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<float*>(data),create_var_id(varname,n,root_tree::id_float).c_str());
           return true;
       }
#endif
       return false;
    }

    bool root_tree::branch(const double* data,const std::string& varname)
    {
        return branch(data,1,varname);
    }

    bool root_tree::branch(const double* data,size_type n,const std::string& varname)
    {
       if(!is_open() or n<1)
       {
           return false;
       }
#if HAVE_ROOT_H
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<double*>(data),create_var_id(varname,n,root_tree::id_double).c_str());
           return true;
       }
#endif
       return false;
    }

    bool root_tree::fill()
    {
        if(!is_open())
        {
            return false;
        }
#if HAVE_ROOT_H
        if(roottree!=NULL)
        {
            static_cast<TTree*>(roottree)->Fill();
            return true;
        }
#endif
        return false;
    }

    std::string root_tree::create_var_id(const std::string& name,size_type n,const char* type)
    {
#if HAVE_ROOT_H
        std::stringstream ss;
        if(n==1)
        {
            ss<<name<<type;
        }
        else
        {
            ss<<name<<'['<<n<<']'<<type;
        }
        return ss.str();
#else
        return "";
#endif
    }
}
