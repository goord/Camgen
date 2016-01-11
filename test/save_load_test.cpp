//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/parni_gen.h>
#include <Camgen/stdrand.h>
#include <Camgen/SM.h>
#include <Camgen/proc_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Facility testing the serialization of event generators. *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"testing saving/loading MC generators..................................."<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;

//    typedef SM model_type;
//    typedef SM::value_type value_type;
//    typedef std::size_t size_type;
//    typedef std::random rng_type;
    
    //////////////////////////////////////////////////////////////////

//    size_type N_events = 100000;
//    size_type N_batch  = 100;
//    size_type N_bins   = 100;
//    grid_modes::type mode = grid_modes::maximum_weights;

    //////////////////////////////////////////////////////////////////


    std::cerr<<"...serialization not implemented yet..."<<std::endl;

	
/*    size_type offset=N_events/N_bins;
    value_type pi=std::acos(-(value_type)1);
    std::string filename="test_output.dat";

    {
	value_type x=1,xmin=0,xmax=10,m=5,w=0.75;
	std::cerr<<"Checking 1D parni save/load on Cauchy distribution..........";
        std::cerr.flush();
	parni<value_type,1,rng_type>* gen=new parni<value_type,1,rng_type>(&x,xmin,xmax,N_bins,mode);
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    gen->integrand()=(value_type)1/(std::pow(x-m,(int)2)+w*w);
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	}
	std::ofstream output_file(filename.c_str());
	gen->save(output_file);
	random_number_stream<value_type,rng_type>::reset_engine();
	std::vector<value_type>xvals(100);
	std::vector<value_type>wvals(100);
	for(int n=0;n<100;++n)
	{
	    gen->generate();
	    xvals[n]=x;
	    wvals[n]=gen->weight();
	}
	delete gen;
	gen=new parni<value_type,1,rng_type>(&x,xmin,xmax,N_bins,mode);
	std::ifstream input_file(filename.c_str());
	gen->load(input_file);
	random_number_stream<value_type,rng_type>::reset_engine();
	for(int n=0;n<100;++n)
	{
	    gen->generate();
	    if(!equals(x,xvals[n]))
	    {
		std::cerr<<"Different value detected after save/load for "<<n<<" throws: "<<x<<" not equal to original value "<<xvals[n]<<'.'<<std::endl;
		return 1;
	    }
	    if(!equals(gen->weight(),wvals[n]))
	    {
		std::cerr<<"Different weight detected after save/load for "<<n<<" throws: "<<gen->weight()<<" not equal to original weight "<<wvals[n]<<'.'<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::remove(filename.c_str());
	std::cerr<<"...........done."<<std::endl;
    }
    {
	vector<value_type,2>x;
	vector<value_type,2>xmin={0,0};
	vector<value_type,2>xmax={10,10};
	value_type m=5,w=0.75;
	std::cerr<<"Checking 2D parni save/load on Cauchy distribution..........";
        std::cerr.flush();
	parni<value_type,2,rng_type>* gen=new parni<value_type,2,rng_type>(&x,xmin,xmax,N_bins,mode);
	for(size_type n=0;n<N_events;++n)
	{
	    gen->generate();
	    gen->integrand()=(value_type)1/(std::pow(std::sqrt(x[0]*x[0]+x[1]*x[1])-m,(int)2)+w*w);
	    gen->update();
	    gen->refresh_cross_section();
	    if(n%N_batch==0 and n!=0)
	    {
		gen->adapt();
	    }
	}
	std::ofstream output_file(filename.c_str());
	gen->save(output_file);
	random_number_stream<value_type,rng_type>::reset_engine();
	std::vector<vector<value_type,2> >xvals(100);
	std::vector<value_type>wvals(100);
	for(int n=0;n<100;++n)
	{
	    gen->generate();
	    xvals[n][0]=x[0];
	    xvals[n][1]=x[1];
	    wvals[n]=gen->weight();
	}
	delete gen;
	gen=new parni<value_type,2,rng_type>(&x,xmin,xmax,N_bins,mode);
	std::ifstream input_file(filename.c_str());
	gen->load(input_file);
	random_number_stream<value_type,rng_type>::reset_engine();
	for(int n=0;n<100;++n)
	{
	    gen->generate();
	    if(!equals(x[0],xvals[n][0]))
	    {
		std::cerr<<"Different value detected after save/load for "<<n<<" throws: "<<x<<" not equal to original value "<<xvals[n]<<'.'<<std::endl;
		return 1;
	    }
	    if(!equals(x[1],xvals[n][1]))
	    {
		std::cerr<<"Different value detected after save/load for "<<n<<" throws: "<<x<<" not equal to original value "<<xvals[n]<<'.'<<std::endl;
		return 1;
	    }
	    if(!equals(gen->weight(),wvals[n]))
	    {
		std::cerr<<"Different weight detected after save/load for "<<n<<" throws: "<<gen->weight()<<" not equal to original weight "<<wvals[n]<<'.'<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::remove(filename.c_str());
	std::cerr<<"...........done."<<std::endl;
    }

    set_initial_state_type(initial_states::partonic);
    set_phase_space_generator_type(phase_space_generators::uniform);
    
    {
	std::string process("h0 > mu-,nu_mubar,mu+,nu_mu");
	std::cout<<"Checking recursive MC generator save/load for "<<process<<"............";
	std::cout.flush();
	CM_algorithm<model_type,1,4>algo(process);
	algo.load();
	algo.construct();
	process_generator<model_type,1,4,rng_type>* gen=new process_generator<model_type,1,4,rng_type>(algo.get_tree_iterator());
	gen->initialise(true);
	std::ofstream output_file(filename.c_str());
	gen->save(output_file);
	random_number_stream<value_type,rng_type>::reset_engine();
	std::vector<value_type>mvals(100);
	std::vector<value_type>wvals(100);
	for(int n=0;n<100;++n)
	{
	    gen->generate();
	    mvals[n]=gen->integrand();
	    wvals[n]=gen->weight();
	}
	delete gen;
	gen=new process_generator<model_type,1,4,rng_type>(algo.get_tree_iterator());
	std::ifstream input_file(filename.c_str());
	gen->load(input_file);
	random_number_stream<value_type,rng_type>::reset_engine();
	for(int n=0;n<100;++n)
	{
	    gen->generate();
	    if(!equals(gen->integrand(),mvals[n]))
	    {
		std::cerr<<"Different integrand detected after save/load for "<<n<<" throws: "<<gen->integrand()<<" not equal to original weight "<<mvals[n]<<'.'<<std::endl;
		return 1;
	    }
	    if(!equals(gen->weight(),wvals[n]))
	    {
		std::cerr<<"Different weight detected after save/load for "<<n<<" throws: "<<gen->weight()<<" not equal to original weight "<<wvals[n]<<'.'<<std::endl;
		return 1;
	    }
	}
	delete gen;
	std::remove(filename.c_str());
	std::cerr<<"...........done."<<std::endl;
    }*/
}
