//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file LH_if.h
    \brief LH_interface class template interface and implementation header.
 */

#ifndef CAMGEN_LH_INTERFACE_H_
#define CAMGEN_LH_INTERFACE_H_

#include <cxxabi.h>
#include <cmath>
#include <Camgen/CM_algo.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Interface class for reading and writing Les-Houches-Accord event files. Can   *
 * be used to feed Camgen's output to parton shower Monte Carlo program, or to  *
 * input hard process events into Camgen, to obtain the matrix element          *
 * correction.                                                                   *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// CM_algorithm subclass with Les Houches IO support.
    
    template<class model_t,std::size_t N_out>class LH_interface: public CM_algorithm<model_t,2,N_out>
    {
	public:

	    /* Base type definition: */

	    typedef CM_algorithm<model_t,2,N_out> base_type;
	    
	    /* Model type definition: */
	    
	    typedef model_t model_type;

	    /* Spacetime type definition: */

	    typedef typename get_spacetime_type<model_t,model_t::dimension>::type spacetime_type;

	    /* Inherited type definitions: */

	    typedef typename base_type::r_value_type r_value_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::tensor_type tensor_type;
	    typedef typename base_type::iterator iterator;
	    typedef typename base_type::const_iterator const_iterator;
	    typedef typename base_type::tree_type tree_type;
	    typedef typename base_type::process_type process_type;
	    typedef typename base_type::phase_space_type phase_space_type;
	    typedef typename base_type::tree_iterator tree_iterator;
	    typedef typename base_type::const_tree_iterator const_tree_iterator;
	    typedef typename base_type::process_iterator process_iterator;
	    typedef typename base_type::const_process_iterator const_process_iterator;
	    
	    /* Number of incoming and outgoing particles: */

	    static const std::size_t N_incoming=2;
	    static const std::size_t N_outgoing=N_out;
	    static const std::size_t N_external=2+N_out;

	    /// Trivial constructor.
	    /// The optional integer argument denotes the external particle
	    /// chosen by the algorithm as the final current in the recursion.
	    /// The amplitudes should not depend on this choice, unless there is
	    /// a bug in Camgen, or the model violates CPT invariance.  

	    LH_interface(std::size_t n=0):base_type(n),beam_id1(0),beam_id2(0),beam_E1(0),beam_E2(0),pdf_authors1(-1),pdf_authors2(-1),pdf_cernlib_code1(-1),pdf_cernlib_code2(-1),weight_switch(1),N_cols(model_t::QCD_colours()),NR_iterations(10)
	    {
		is.precision(10);
		is.setf(std::ios::scientific,std::ios::floatfield);
		os.precision(10);
		os.setf(std::ios::scientific,std::ios::floatfield);
	    }
	    
	    /// Constructor specifying the process.
	    /// The first argument should be of the form "phi1,...,phiN_in >
	    /// psi1,...,psiN_out", and the second argument is the optional final
	    /// current (see trivial constructor).
	    
	    LH_interface(const std::string& proc,std::size_t n=0):base_type(proc,n),beam_id1(0),beam_id2(0),beam_E1(0),beam_E2(0),pdf_authors1(-1),pdf_authors2(-1),pdf_cernlib_code1(-1),pdf_cernlib_code2(-1),weight_switch(1),N_cols(model_t::QCD_colours()),NR_iterations(10)
	    {
		is.precision(10);
		is.setf(std::ios::scientific,std::ios::floatfield);
		os.precision(10);
		os.setf(std::ios::scientific,std::ios::floatfield);
	    }

	    /// Beam PDG id input (default values: 0).

	    void set_beams_ids(int i1,int i2)
	    {
		beam_id1=i1;
		beam_id2=i2;
	    }

	    /// Proton-proton beam input.
	    
	    void set_LHC_beams()
	    {
		beam_id1=2212;
		beam_id2=2212;
	    }

	    /// Proton-antiproton beam input.

	    void set_tevatron_beams()
	    {
		beam_id1=2212;
		beam_id2=-2212;
	    }

	    /// Electron-positron beam input.

	    void set_LEP_beams()
	    {
		beam_id1=11;
		beam_id2=-11;
	    }

	    /// Beam energies input (default values: 0).

	    void set_beam_energies(const r_value_type& E1,const r_value_type& E2)
	    {
		beam_E1=E1;
		beam_E2=E2;
	    }

	    /// pdf authors input (default values: -1).

	    void set_pdf_authors(int n1,int n2)
	    {
		pdf_authors1=n1;
		pdf_authors2=n2;
	    }

	    /// pdf authors cernlib code input (default values: -1).

	    void set_pdf_cernlib_codes(int n1,int n2)
	    {
		pdf_cernlib_code1=n1;
		pdf_cernlib_code2=n2;
	    }

	    /// Destructor, closing open files.

	    ~LH_interface()
	    {
		if(is.is_open())
		{
		    is.close();
		}
	    }

	    /// Loads the inserted processes and sorts them by PDG-id.
	    
	    void load()
	    {
		this->base_type::load();
		this->sort_by_pdg_id();
	    }

	    /// Starts reading a Les-Houches event file.

	    void read(const char* fname)
	    {
		is.open(fname,std::ios::in);
		if(is.is_open())
		{
		    log(log_level::message)<<"Reading Les Houches interface file "<<fname<<".........."<<endlog;
		    char firstline[256];
		    do
		    {
			is.getline(firstline,256);
		    }
		    while(std::string(firstline)!=initbegin and !is.eof());
		    if(is.eof())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached without initialisation block"<<endlog;
			return;
		    }
		    is>>beam_id1>>beam_id2>>beam_E1>>beam_E2>>pdf_authors1>>pdf_authors2>>pdf_cernlib_code1>>pdf_cernlib_code2;
		    int n=0;
		    is>>weight_switch>>n>>this->Xsec>>this->sigma_Xsec>>this->max_wght>>n;
		    is.getline(firstline,256);
		}
		else
		{
		    log(log_level::error)<<CAMGEN_STREAMLOC<<"Error opening file "<<fname<<endlog;
		}
	    }

	    /// Starts writing a Les-Houches event file.

	    void write(const char* fname)
	    {
		os.open(fname,std::ios::out);
		if(os.is_open())
		{
		    std::cout<<"Writing Les Houches interface file "<<fname<<".........."<<std::endl;
		    os<<"<LesHouchesEvents version=\"1.0\">"<<std::endl;
		    os<<"<!--"<<std::endl;
		    os<<"File written by "<<abi::__cxa_demangle(typeid(*this).name(),0,0,NULL)<<" on "<<__DATE__<<" at "<<__TIME__<<std::endl;
		    os<<"-->"<<std::endl;
		    os<<"<init>"<<std::endl;
		    os<<"\t"<<beam_id1<<"\t"<<beam_id2<<"\t"<<beam_E1<<"\t"<<beam_E2<<"\t"<<pdf_authors1<<"\t"<<pdf_authors2<<"\t"<<pdf_cernlib_code1<<"\t"<<pdf_cernlib_code2<<std::endl;
		    os<<"\t"<<weight_switch<<"\t 1 \t"<<this->Xsec<<"\t"<<this->sigma_Xsec<<"\t"<<this->max_wght<<"\t9999"<<std::endl;
		    os<<"</init>"<<std::endl;
		}
		else
		{
		    log(log_level::error)<<CAMGEN_STREAMLOC<<"error opening file "<<fname<<endlog;
		}
	    }

	    /// Reads the next event in the file.
	    /// Returns a boolean denoting whether the read-in was succesful.

	    bool read_event()
	    {
		if(is.is_open())
		{
		    char firstline[256];
		    do
		    {
			is.getline(firstline,256);
		    }
		    while(std::string(firstline)!=evtbegin and !is.eof());
		    if(is.eof())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached"<<endlog;
			return false;
		    }
		    is>>n_particles>>process_id>>event_weight>>qcd_scale>>alpha_em>>alpha_s;
		    if(n_particles<N_external)
		    {
			return false;
		    }
		    unsigned k=0;
		    CM_momentum.assign(0);
		    int pdgid,dir,dummy,col,anticol;
		    momentum_type p;
		    r_value_type m,s;
		    for(std::size_t n=0;n<n_particles;++n)
		    {
			is>>pdgid>>dir>>dummy>>dummy>>col>>anticol;
			for(std::size_t j=1;j<p.static_size;++j)
			{
			    is>>p[j];
			}
			is>>p[0]>>m>>s>>s;
			if(dir==-1)
			{
			    if(k==N_external)
			    {
				return false;
			    }
			    pdg_id_vec[k]=pdgid;
			    col_vec[k]=col;
			    anti_col_vec[k]=anticol;
			    masses[k]=m;
			    momenta[k]=p;
			    CM_momentum+=p;
			    ++k;
			}
			if(dir==1)
			{
			    if(k==N_external)
			    {
				return false;
			    }
			    pdg_id_vec[k]=pdgid;
			    col_vec[k]=col;
			    anti_col_vec[k]=anticol;
			    masses[k]=m;
			    momenta[k]=p;
			    ++k;
			}
		    }
		    if(k!=N_external)
		    {
			return false;
		    }
		    Ecm=std::sqrt(spacetime_type::dot(CM_momentum,CM_momentum));
		    return true;
		}
		return false;
	    }

	    /// Writes the current event to the Les-Houches event file.
	    /// Returns a boolean denoting whether the output was succesful.

	    bool write_event()
	    {
		if(os.is_open())
		{
		    vector<int,2+N_out>cols,anticols;
		    if(this->get_colour_generator()!=NULL)
		    {
			this->get_colour_generator()->LH_output(cols,anticols);
			cols.assign(0);
			anticols.assign(0);
		    }
		    else if(is.is_open())
		    {
			cols=col_vec;
			anticols=anti_col_vec;
		    }
		    os<<"<event>"<<std::endl;
		    os<<(2+N_out)<<"\t9999\t"<<this->wght<<"\t"<<model_t::QCD_scale<<"\t"<<model_t::alpha<<"\t"<<model_t::alpha_s<<std::endl;
		    for(unsigned i=0;i<2;++i)
		    {
			os<<std::setw(4)<<this->get_pdg_id(i);
			os<<std::setw(6)<<-1;
			os<<std::setw(6)<<0;
			os<<std::setw(6)<<0;
			os<<std::setw(6)<<cols[i];
			os<<std::setw(6)<<anticols[i];
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(1);
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(2);
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(3);
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(0);
			os<<std::setw(20)<<std::sqrt(std::abs(spacetime_type::dot(this->get_phase_space(i)->momentum(),this->get_phase_space(i)->momentum())));
			os<<std::setw(6)<<"0.";
			os<<std::setw(6)<<"9.";
			os<<std::endl;
		    }
		    for(unsigned i=2;i<N_out+2;++i)
		    {
			os<<std::setw(4)<<this->get_pdg_id(i);
			os<<std::setw(6)<<1;
			os<<std::setw(6)<<1;
			os<<std::setw(6)<<2;
			os<<std::setw(6)<<cols[i];
			os<<std::setw(6)<<anticols[i];
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(1);
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(2);
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(3);
			os<<std::setw(20)<<this->get_phase_space(i)->momentum(0);
			os<<std::setw(20)<<std::sqrt(std::abs(spacetime_type::dot(this->get_phase_space(i)->momentum(),this->get_phase_space(i)->momentum())));
			os<<std::setw(6)<<"0.";
			os<<std::setw(6)<<"9.";
			os<<std::endl;
		    }
		    os<<"</event>"<<std::endl;
		    return true;
		}
		return false;
	    }

	    /* Utility printing facilities: */

	    void print_info()
	    {
		std::cout<<"reading process "<<pdg_id_vec<<".....";
		this->set_process_pdg_ids(pdg_id_vec);
		std::cout<<std::endl<<"interfaced flavour event:"<<std::endl;
		for(unsigned i=0;i<2;++i)
		{
		    std::cout<<this->get_pdg_id(i)<<"  ";
		}
		std::cout<<">";
		for(unsigned i=2;i<N_external;++i)
		{
		    std::cout<<"  "<<this->get_pdg_id(i);
		}
		std::cout<<std::endl<<std::endl;
	    }

	    void print_momenta() const
	    {
		std::cout.precision(10);
		std::cout.setf(std::ios::fixed,std::ios::floatfield);
		momentum_type k;
		k.assign(0);
		for(size_type i=0;i<2;++i)
		{
		    std::cout<<-1<<"\t"<<momenta[i]<<"\t"<<std::sqrt(std::fabs(momenta[i].s))<<"\t"<<masses[i]<<std::endl;
		    k+=momenta[i];
		}
		std::cout<<"sum: "<<k<<std::endl<<std::endl;
		k.assign(0);
		for(size_type i=2;i<N_external;++i)
		{
		    std::cout<<1 <<"\t"<<momenta[i]<<"\t"<<std::sqrt(std::fabs(momenta[i].s))<<"\t"<<masses[i]<<std::endl;
		    k+=momenta[i];
		}
		std::cout<<"sum: "<<k<<std::endl<<std::endl;
	    }

	    void find_process()
	    {
		this->set_process_pdg_ids(pdg_id_vec);
	    }

	    /// Inserts the event in the vertex tree.
	    /// Starts by searching the internal process list for the most
	    /// recently read process from the LH file. If this process is not
	    /// found, will automatically an entry in the process list and
	    /// construct a corresponding vertex tree. Then, inserts the momenta
	    /// and colours in the phase space instances of the external
	    /// particles. 

	    void input_event()
	    {
		if(!this->tree_it->is_empty())
		{
		    this->set_process_pdg_ids(pdg_id_vec);
		    model_t::set_QCD_scale(qcd_scale);
		    model_t::set_alpha(alpha_em);
		    model_t::set_alpha_s(alpha_s);
		    for(std::size_t n=0;n<N_external;++n)
		    {
			this->get_phase_space(n)->momentum()=momenta[n];
			if(col_vec[n]!=0)
			{
			    if(anti_col_vec[n]!=0)
			    {
				this->get_phase_space(n)->colour(0)=col_vec[n]%N_cols;
				this->get_phase_space(n)->colour(1)=anti_col_vec[n]%N_cols;
			    }
			    else
			    {
				this->get_phase_space(n)->colour(0)=col_vec[n]%N_cols;
			    }
			}
			else
			{
			    if(anti_col_vec[n]!=0)
			    {
				this->get_phase_space(n)->colour(0)=anti_col_vec[n]%N_cols;
			    }
			}
		    }
		}
	    }

	    /* Function rescaling momenta to obtain the correct masses: */

	    void rescale_outgoing_momenta()
	    {
		if(!this->tree_it->is_empty())
		{
		    r_value_type msq[N_out],f,df,pvec[N_out],d;
		    for(size_type n=0;n<N_out;++n)
		    {
			msq[n]=(this->get_phase_space(n+2)->mass())*(this->get_phase_space(n+2)->mass());
			pvec[n]=spacetime_type::space_dot(momenta[n+2],momenta[n+2]);
		    }
		    r_value_type xi=1;
		    for(size_type i=0;i<NR_iterations;++i)
		    {
			f=-Ecm;
			df=0;
			for(size_type n=0;n<N_out;++n)
			{
			    d=std::sqrt(xi*xi*pvec[n]+msq[n]);
			    f+=d;
			    df+=pvec[n]/d;
			}
			xi-=f/(xi*df);
		    }
		    r_value_type wnum=0;
		    r_value_type wdenom=0;
		    r_value_type wfac=1;
		    for(size_type i=2;i<N_external;++i)
		    {
			wnum+=(spacetime_type::space_dot(momenta[i],momenta[i])/momenta[i][0]);
			wfac*=momenta[i][0];
			for(size_type mu=1;mu<model_t::dimension;++mu)
			{
			    momenta[i][mu]*=xi;
			}
			momenta[i][0]=std::sqrt(xi*xi*pvec[i-2]+msq[i-2]);
			wdenom+=(spacetime_type::space_dot(momenta[i],momenta[i])/momenta[i][0]);
			wfac/=momenta[i][0];
		    }
		    int pw=(model_t::dimension-2)*N_out+1-model_t::dimension;
		    scaling_weight=wfac*(wnum/wdenom)*std::pow(xi,pw);
		}
	    }

	    /* Function boosting the external momenta to the CM frame: */

	    void boost_to_CM_frame()
	    {
		r_value_type kq;
		for(size_type i=0;i<N_external;++i)
		{
		    kq=spacetime_type::dot(CM_momentum,momenta[i]);
		    for(size_type mu=1;mu<model_t::dimension;++mu)
		    {
			momenta[i][mu]-=((kq/Ecm+momenta[i][0])*CM_momentum[mu]/(Ecm+CM_momentum[0]));
		    }
		    momenta[i][0]=kq/Ecm;
		}
	    }

	    /* Function boosting the external momenta from the CM frame: */

	    void boost_from_CM_frame()
	    {
		momentum_type k=-CM_momentum;
		k[0]=-k[0];
		r_value_type kq;
		for(size_type i=0;i<N_external;++i)
		{
		    kq=spacetime_type::dot(k,momenta[i]);
		    for(size_type mu=1;mu<model_t::dimension;++mu)
		    {
			momenta[i][mu]-=((kq/Ecm+momenta[i][0])*k[mu]/(Ecm+k[0]));
		    }
		    momenta[i][0]=kq/Ecm;
		}
	    }

	    /// Function boosting and rescaling the inputted momenta to make them on-shell.
	    /// For events with effective quark masses etc., this routine
	    /// recales the input momenta to match with Camgen's internal masses
	    /// while retaining momentum conservation. As a result, the event
	    /// attains an extra weight, which is accessed by the
	    /// get_scaling_weight function.  

	    void rescale_phase_space()
	    {
		find_process();
		boost_to_CM_frame();
		rescale_outgoing_momenta();
		boost_from_CM_frame();
	    }

	    /// Determines the number of Newton-Raphson iterations for obtaining momentum scale factor.
	    
	    void set_NR_iterations(size_type n)
	    {
		NR_iterations=n;
	    }

	    /// Returns the extra weight due the rescale_phase_space call.

	    const r_value_type& get_scaling_weight() const
	    {
		return scaling_weight;
	    }

	private:

	    /* Input file stream object: */

	    std::ifstream is;

	    /* Output file stream object: */

	    std::ofstream os;

	    /* Beam particle id's: */

	    int beam_id1,beam_id2;

	    /* Beam particle energies: */

	    r_value_type beam_E1,beam_E2;

	    /* PDF authors and cernlib codes: */

	    int pdf_authors1,pdf_authors2,pdf_cernlib_code1,pdf_cernlib_code2;

	    /* Les Houches output weight switch: */

	    int weight_switch;
	    
	    /* Number of external particles described by the input file: */
	    
	    std::size_t n_particles;

	    /* Event weight given by the input file: */

	    r_value_type event_weight;
	    
	    /* QCD scale, alpha weak and alpha strong used in the event: */
	    
	    r_value_type qcd_scale,alpha_s,alpha_em;

	    /* Number of QCD (quark-) colours in the model" */

	    std::size_t N_cols;
	    
	    /* LH-file process id: */
	    
	    int process_id;

	    /* Vector of pdg-codes to be read from the input file: */

	    vector<int,N_external>pdg_id_vec;
	    
	    /* Colour and anti-colour vectors: */
	    
	    vector<int,N_external>col_vec;
	    vector<int,N_external>anti_col_vec;

	    /* Vector of momenta of external particles: */

	    vector<momentum_type,N_external>momenta;

	    /* Center-of-momentum vector: */

	    momentum_type CM_momentum;

	    /* Center-of-mass energy: */

	    r_value_type Ecm;
	    
	    /* External particle masses: */
	    
	    vector<r_value_type,N_external>masses;

	    /* Utility search string: */

	    static const std::string evtbegin,initbegin;

	    /* Number of Newton-Raphson iterations to obtain the scaling factor:
	     * */

	    size_type NR_iterations;

	    /* Extra event weight due to rescaling: */

	    r_value_type scaling_weight;
    };
    template<class model_t,std::size_t N_out>const std::size_t LH_interface<model_t,N_out>::N_incoming;
    template<class model_t,std::size_t N_out>const std::size_t LH_interface<model_t,N_out>::N_outgoing;
    template<class model_t,std::size_t N_out>const std::size_t LH_interface<model_t,N_out>::N_external;
    template<class model_t,std::size_t N_out>const std::string LH_interface<model_t,N_out>::evtbegin="<event>";
    template<class model_t,std::size_t N_out>const std::string LH_interface<model_t,N_out>::initbegin="<init>";
}

#endif /*CAMGEN_LH_INTERFACE_H_*/

