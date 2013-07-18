//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ascii_if.h
    \brief ASCII interface output for phase space generators.
 */

#ifndef CAMGEN_ASCII_IF_H_
#define CAMGEN_ASCII_IF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * ASCII datafile output interface class implementation. Creates a datafile  *
 * where each row represents an event.                                       *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <fstream>
#include <Camgen/if_output.h>

namespace Camgen
{
    /// ACII-file output interface class.
    
    template<class model_t>class ascii_file: public interface_output<model_t>
    {
	typedef interface_output<model_t> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_type::value_type value_type;
	    typedef vector<value_type,model_t::dimension> momentum_type;
	    typedef typename std::vector<value_type>::iterator weight_iterator;
	    typedef typename std::vector<value_type>::const_iterator const_weight_iterator;
	    
	    typedef typename std::map<std::string,const momentum_type*>::iterator vector_iterator;
	    typedef typename std::map<std::string,const value_type*>::iterator value_iterator;
	    typedef typename std::map<std::string,const int*>::iterator integer_iterator;
	    typedef typename std::map<std::string,const bool*>::iterator boolean_iterator;

	    /* Constructor with file name argument */

	    ascii_file(const std::string& file_name_):base_type(file_name_),line(0){}

	    /* Constructor with file name, tree name and tree description arguments: */

	    ascii_file(const std::string& file_name_,const std::string description_):base_type(file_name_,description_),line(0){}

	    /* Destructor: */

	    ~ascii_file(){}

	    /* Creation method implementation: */

	    interface_output<model_t>* create(const std::string& file_name_) const
	    {
		return new ascii_file<model_t>(file_name_,this->description);
	    }

	    /* Opens the (temporary) datafile. */

	    bool open_file()
	    {
		std::string fname=this->file_name+".dat";
		ofs.open(fname.c_str());
		if(ofs.is_open())
		{
		    return true;
		}
		return false;
	    }

	    /* Closes the datafile. */

	    bool close_file()
	    {
		if(ofs.is_open())
		{
		    ofs.close();
		}
		return !(ofs.is_open());
	    }

	    /* Virtual method adding a branch holding a Lorentz vector: */

	    bool branch(const momentum_type* p,const std::string& varname)
	    {
		if(line!=0)
		{
		    return false;
		}
		vectors[varname]=p;
		values.erase(varname);
		integers.erase(varname);
		booleans.erase(varname);
		return true;
	    }

	    /* Virtual method adding a branch holding a floating-point number:
	     * */

	    bool branch(const value_type* x,const std::string& varname)
	    {
		if(line!=0)
		{
		    return false;
		}
		vectors.erase(varname);
		values[varname]=x;
		integers.erase(varname);
		booleans.erase(varname);
		return true;
	    }

	    /* Virtual method adding a branch holding an integer. */

	    bool branch(const int* n,const std::string& varname)
	    {
		if(line!=0)
		{
		    return false;
		}
		vectors.erase(varname);
		values.erase(varname);
		integers[varname]=n;
		booleans.erase(varname);
		return true;
	    }

	    /* Virtual method adding a branch holding a boolean. */

	    bool branch(const bool* b,const std::string& varname)
	    {
		if(line!=0)
		{
		    return false;
		}
		vectors.erase(varname);
		values.erase(varname);
		integers.erase(varname);
		booleans[varname]=b;
		return true;
	    }

	    /* Writes the event to the output file. */

	    bool write_event()
	    {
		if(!ofs.is_open())
		{
		    return false;
		}
		if(line==0)
		{
		    ofs<<'#';
		    for(vector_iterator it=vectors.begin();it!=vectors.end();++it)
		    {
			for(typename momentum_type::size_type i=0;i<model_type::dimension;++i)
			{
			    std::stringstream ss;
			    ss<<it->first<<'['<<i<<']';
			    ofs<<std::setw(20)<<ss.str();
			}
		    }
		    for(value_iterator it=values.begin();it!=values.end();++it)
		    {
			ofs<<std::setw(20)<<it->first;
		    }
		    for(integer_iterator it=integers.begin();it!=integers.end();++it)
		    {
			ofs<<std::setw(10)<<it->first;
		    }
		    for(boolean_iterator it=booleans.begin();it!=booleans.end();++it)
		    {
			ofs<<std::setw(10)<<it->first;
		    }
		    ofs<<std::endl;
		    line=1;
		}
		for(vector_iterator it=vectors.begin();it!=vectors.end();++it)
		{
		    for(typename momentum_type::size_type i=0;i<model_type::dimension;++i)
		    {
			ofs<<std::setw(20)<<(*(it->second))[i];
		    }
		}
		for(value_iterator it=values.begin();it!=values.end();++it)
		{
		    ofs<<std::setw(20)<<*(it->second);
		}
		for(integer_iterator it=integers.begin();it!=integers.end();++it)
		{
		    ofs<<std::setw(10)<<*(it->second);
		}
		for(boolean_iterator it=booleans.begin();it!=booleans.end();++it)
		{
		    ofs<<std::setw(10)<<*(it->second);
		}
		ofs<<std::endl;
		++line;
		return true;
	    }

	private:

	    std::ofstream ofs;

	    int line;

	    std::map<std::string,const momentum_type*> vectors;
	    std::map<std::string,const value_type*> values;
	    std::map<std::string,const int*> integers;
	    std::map<std::string,const bool*> booleans;
    };
}

#endif /*CAMGEN_ASCII_IF_H_*/

