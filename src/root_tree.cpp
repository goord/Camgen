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

#if HAVE_ROOT_H_
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
        return open(fname.c_str(),tname.c_str(),NULL);
    }

    bool root_tree::open(const std::string& fname,const std::string& tname, const std::string& descr)
    {
        return open(fname.c_str(),tname.c_str(),descr.c_str());
    }

    bool root_tree::open(const char* fname,const char* tname,const char* descr)
    {
        if(is_open())
        {
            close();
        }
#if HAVE_ROOT_H_
        std::string filename(fname);
        std::size_t n=filename.find_last_of('.');
        if(n!=filename.size()-5 and filename.substr(n+1,5)!="root")
        {
            filename+=".root";
        }
        rootfile=new TFile(filename.c_str(),"RECREATE");
        roottree=new TTree(tname,descr);
        return rootfile->IsOpen();
#else
        return false;
#endif
    }

    bool root_tree::is_open() const
    {
#if HAVE_ROOT_H_
        if(rootfile!=NULL)
        {
            return static_cast<TFile*>(rootfile)->IsOpen();
        }
#endif
        return false;
    }

    void root_tree::close()
    {
#if HAVE_ROOT_H_
        if(rootfile!=NULL)
        {
            static_cast<TFile*>(rootfile)->Write();
            static_cast<TFile*>(rootfile)->Close();
            delete rootfile;
            rootfile=NULL;
            delete roottree;
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
#if HAVE_ROOT_H_
       std::string id;
       if(n==1)
       {
           id=std::string(root_tree::id_bool);
       }
       else
       {
           std::stringstream ss;
           ss<<name<<'['<<n<<']'<<root_tree::id_bool;
           id=ss.str();
       }
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<bool*>(data),id.c_str());
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
#if HAVE_ROOT_H_
       std::string id;
       if(n==1)
       {
           id=std::string(root_tree::id_int);
       }
       else
       {
           std::stringstream ss;
           ss<<name<<'['<<n<<']'<<root_tree::id_int;
           id=ss.str();
       }
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<int*>(data),id.c_str());
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
#if HAVE_ROOT_H_
       std::string id;
       if(n==1)
       {
           id=std::string(root_tree::id_long);
       }
       else
       {
           std::stringstream ss;
           ss<<name<<'['<<n<<']'<<root_tree::id_long;
           id=ss.str();
       }
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<long int*>(data),id.c_str());
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
#if HAVE_ROOT_H_
       std::string id;
       if(n==1)
       {
           id=std::string(root_tree::id_float);
       }
       else
       {
           std::stringstream ss;
           ss<<name<<'['<<n<<']'<<root_tree::id_float;
           id=ss.str();
       }
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<float*>(data),id.c_str());
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
#if HAVE_ROOT_H_
       std::string id;
       if(n==1)
       {
           id=std::string(root_tree::id_double);
       }
       else
       {
           std::stringstream ss;
           ss<<name<<'['<<n<<']'<<root_tree::id_double;
           id=ss.str();
       }
       if(roottree!=NULL)
       {
           static_cast<TTree*>(roottree)->Branch(varname.c_str(),const_cast<double*>(data),id.c_str());
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
#if HAVE_ROOT_H_
        if(roottree!=NULL)
        {
            static_cast<TTree*>(roottree)->Fill();
            return true;
        }
#endif
        return false;
    }
}
