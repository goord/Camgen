//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <config.h>
#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>
#include <Camgen/root_if.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Tests for ROOT interface to event/process generators  *
 *                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
#if HAVE_ROOT_H

    typedef SM model_type;
    typedef std::random rn_engine;
    typedef typename model_type::value_type value_type;
    typedef std::size_t size_type;

    license_print::disable();

    file_utils::create_directory("test_output/root_test");
    size_type n_evts=10000;

    {
        Camgen::log.enable_level=log_level::error;
        std::cerr<<"Checking root tree wrapper for toy distributions............";
        std::cerr.flush();

        vector<double,3>p;
        float f;
        int n;
        bool q;

        uniform_sphere<double,2,rn_engine>gen(&p);
        root_tree tree;
        
        if(!tree.open("test_output/root_test/uni_sphere","sphere","test tree"))
        {
            return 1;
        }
        if(!tree.is_open())
        {
            return 1;
        }
        if(!tree.branch(p.data,3,"vector"))
        {
            return 1;
        }
        if(!tree.branch(&f,"sum"))
        {
            return 1;
        }
        if(!tree.branch(&n,"diff"))
        {
            return 1;
        }
        if(!tree.branch(p.data,"sign"))
        {
            return 1;
        }
        for(size_type i=0;i<n_evts;++i)
        {
            gen.generate();
            f=(float)(p[0]+p[1]+p[2])/3;
            n=(int)(p[0]-p[1]);
            q=(p[0]>0);
            if(!tree.fill())
            {
                return 1;
            }
        }
        tree.close();
        std::cerr<<"done."<<std::endl;
        Camgen::log.enable_level=log_level::warning;
    }

#endif

    return 0;
}
