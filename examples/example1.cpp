#include <Camgen/QCD.h>
#include <Camgen/CM_algo.h>
#include <Camgen/proc_gen.h>
#include <Camgen/rcarry.h>

using namespace Camgen;

std::ostream& print_colours_momenta(const CM_algorithm<QCD,2,10>& algo,std::ostream& os)
{
    os<<"gluon colours and momenta:"<<std::endl<<std::endl;
    os<<"g    cols  \tE                  \tpx                   \tpy                    \tpz                 "<<std::endl;
    os<<"---------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    os.precision(10);
    os.setf(std::ios::scientific,std::ios::floatfield);
    for(unsigned i=0;i<12;++i)
    {
	os<<i;
	os<<"\t("<<algo.get_phase_space(i)->colour(0)<<","<<algo.get_phase_space(i)->colour(1)<<")";
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

    CM_algorithm<QCD,2,10>algo("g,g > g,g,g,g,g,g,g,g,g,g");
    algo.load();
    algo.construct();
    std::cout<<"Nr of Feynman diagrams:"<<algo.count_diagrams()<<std::endl;

    double Ecm=1000;
    set_beam_energy(1,0.5*Ecm);
    set_beam_energy(2,0.5*Ecm);

    set_initial_state_type(initial_states::partonic);
    set_phase_space_generator_type(phase_space_generators::uniform);
    set_colour_generator_type(colour_generators::flow_sampling);

    process_generator<QCD,2,10,rcarry>gen(algo.get_tree_iterator());
    gen.generate();
    for(unsigned i=0;i<12;++i)
    {
	algo.get_phase_space(i)->helicity_phase(-1)=1.0;
	algo.get_phase_space(i)->helicity_phase(1)=0.0;
    }

    std::ofstream os;
    if(publisher_mode)
    {
	os.open("example1.txt");
	os<<"This is a file containing the checking results of some gg->10g amplitudes in example1.cpp. "<<std::endl<<std::endl;
	print_colours_momenta(algo,os);
    }
    print_colours_momenta(algo,std::cout);
    
    std::complex<double>amplitude=algo.evaluate();

    if(publisher_mode)
    {
	os<<std::endl<<"The (------------) helicity amplitude is "<<amplitude<<std::endl<<std::endl;
    }
    std::cout<<std::endl<<"The (------------) helicity amplitude is "<<amplitude<<std::endl<<std::endl;
    
    algo.get_phase_space(0)->helicity_phase(-1)=0.0;
    algo.get_phase_space(0)->helicity_phase(1)=1.0;
    
    amplitude=algo.evaluate();

    if(publisher_mode)
    {
	os<<std::endl<<"The (+-----------) helicity amplitude is "<<amplitude<<std::endl<<std::endl;
    }
    std::cout<<std::endl<<"The (+-----------) helicity amplitude is "<<amplitude<<std::endl<<std::endl;
    
    algo.get_phase_space(0)->helicity_phase(-1)=1.0;
    algo.get_phase_space(0)->helicity_phase(1)=0.0;
    algo.get_phase_space(2)->helicity_phase(-1)=0.0;
    algo.get_phase_space(2)->helicity_phase(1)=1.0;
    
    amplitude=algo.evaluate();

    if(publisher_mode)
    {
	os<<std::endl<<"The (--+---------) helicity amplitude is "<<amplitude<<std::endl<<std::endl;
    }
    std::cout<<std::endl<<"The (--+---------) helicity amplitude is "<<amplitude<<std::endl<<std::endl;

    if(publisher_mode)
    {
	os.close();
    }
}
