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
#include <Camgen/if_engine.h>
#include <Camgen/gen_if.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Tests for ROOT interface to event/process generators  *
 *                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

template<class model_t>class test_output: public interface_engine<model_t>
{
    public:

        typedef interface_engine<model_t> base_type;
        typedef typename base_type::value_type value_type;
        typedef typename base_type::momentum_type momentum_type;

        momentum_type p1;
        value_type pT2;
        value_type eta2;
        value_type M12;

        interface_engine<model_t>* clone() const
        {
            return new test_output(*this);
        }

        void fill()
        {
            p1=this->p(1);
            pT2=this->pT(2);
            eta2=this->eta(2);
            M12=this->s(1,2);
        }

    protected:

        void add_variables()
        {
            this->add_variable(p1,"p1");
            this->add_variable(pT2,"pT2");
            this->add_variable(eta2,"eta2");
            this->add_variable(M12,"M12");
        }
};

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

    {
	model_type::M_h0=200;
	model_type::refresh_widths();
	std::string process("h0 > e-,nu_ebar,mu+,nu_mu");
	std::string fname("test_output/root_test/h_WW_2l2n");
	std::cerr<<"Checking root interface for "<<process<<"............";
	std::cerr.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
        process_generator_factory<model_type,1,4,rn_engine> factory;
        process_generator<model_type,1,4,rn_engine>* proc_gen=factory.create_generator(algo.get_tree_iterator());
        proc_gen->refresh_Ecm();
        generator_interface<model_type>* gen_if=new generator_interface<model_type>(proc_gen,new root_interface<model_type>(fname,"test-tree"),new test_output<model_type>());
        for(size_type i=0;i<n_evts;++i)
        {
            proc_gen->generate();
            gen_if->fill();
        }
        gen_if->write_statistics();
        gen_if->write();
    }

#endif

    return 0;
}
