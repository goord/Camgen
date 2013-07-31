#include "SMseesaw.h"
#include <Camgen/CM_algo.h>
#include <Camgen/proc_gen.h>
#include <Camgen/rcarry.h>

using namespace Camgen;

std::ostream& print_colours_momenta(const CM_algorithm<SMseesaw,2,4>& algo,std::ostream& os)
{
    os<<"particle colours and momenta:"<<std::endl<<std::endl;
    os<<"     cols  \tE                  \tpx                   \tpy                    \tpz                 "<<std::endl;
    os<<"---------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    os.precision(10);
    os.setf(std::ios::scientific,std::ios::floatfield);
    for(unsigned i=0;i<6;++i)
    {
	os<<i;
	if(i<4)
	{
	    os<<"\t("<<algo.get_phase_space(i)->colour(0)<<")";
	}
	else
	{
	    os<<"\t( )";
	}
	os<<"\t"<<algo.get_phase_space(i)->momentum(0);
	os<<"\t"<<algo.get_phase_space(i)->momentum(1);
	os<<"\t"<<algo.get_phase_space(i)->momentum(2);
	os<<"\t"<<algo.get_phase_space(i)->momentum(3)<<std::endl;
    }
    return os;
}

int main(int argc,char *argv[])
{
    bool publisher_mode=false;

    if(argc>1)
    {
	std::string option=argv[1];
	if(option=="-print")
	{
	    publisher_mode=true;
	}
    }

    CM_algorithm<SMseesaw,2,4>algo("d,ubar > u,dbar,e-,e-");
    algo.load();
    algo.construct();
    algo.print_tree(std::cout);
    std::cout<<"Nr of Feynman diagrams:"<<algo.count_diagrams()<<std::endl;

    double Ecm=500;
    set_beam_energy(1,0.5*Ecm);
    set_beam_energy(2,0.5*Ecm);

    set_initial_state_type(initial_states::partonic);
    set_phase_space_generator_type(phase_space_generators::uniform);
    set_colour_generator_type(colour_generators::flow_sampling);
    set_helicity_generator_type(helicity_generators::summation);

    process_generator<SMseesaw,2,4,rcarry>gen(algo.get_tree_iterator());
    gen.generate();

    std::ofstream os;
    if(publisher_mode)
    {
	os.open("example2.txt");
	os<<"This is a file containing the checking result of a spin-summed d,ubar->u,dbar+2e- amplitude in example2.cpp. "<<std::endl<<std::endl;
	print_colours_momenta(algo,os);
    }
    print_colours_momenta(algo,std::cout);
    
    std::complex<double>amplitude=algo.evaluate();

    if(publisher_mode)
    {
	os<<std::endl<<"The spin-summed amplitude is "<<amplitude<<std::endl<<std::endl;
    }
    std::cout<<std::endl<<"The spin-summed amplitude is "<<amplitude<<std::endl<<std::endl;
}
