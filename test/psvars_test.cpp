//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/QED.h>
#include <Camgen/ps_gen_base.h>
#include <Camgen/stdrand.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Program checking the calculation of phase space variables by the ps_gen_base  *
 * class in Camgen.                                                              *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

using namespace Camgen;

/* Dummy model containing spacetime info and beam direction: */

class model_mock
{
    public:

	typedef double value_type;
	typedef Minkowski_type spacetime_type;

	static const std::size_t dimension=4;
	static const int beam_direction=2;
};

/* Method generating 3 momenta at the given CM momentum scale: */

template<std::size_t N>void generate_momenta(vector<double,N>& p1,vector<double,N>& p2,vector<double,N>& p3, double scale)
{
    random_number_stream<double,std::random>rng;
    double P1=0,P2=0,P3=0;
    vector<double,N>ptot;
    double Etot=0;
    do
    {
	for(unsigned j=1;j<N;++j)
	{
	    p1[j]=rng(-1,1);
	    P1+=p1[j]*p1[j];

	    p2[j]=rng(-1,1);
	    P2+=p2[j]*p2[j];

	    p3[j]=rng(-1,1);
	    P3+=p3[j]*p3[j];

	    ptot[j]=p1[j]+p2[j]+p3[j];
	}
	p1[0]=std::sqrt(P1)+rng();
	p2[0]=std::sqrt(P2)+rng();
	p3[0]=std::sqrt(P3)+rng();

	ptot[0]=p1[0]+p2[0]+p3[0];

	Etot=ptot[0]*ptot[0];
	for(unsigned j=1;j<N;++j)
	{
	    Etot-=ptot[j]*ptot[j];
	}
    }
    while(Etot<=0);
    if(scale!=0)
    {
	double factor=scale/std::sqrt(Etot);
	p1*=factor;
	p2*=factor;
	p3*=factor;
    }
}

int main()
{
    typedef model_mock model_type;
    typedef model_type::value_type value_type;
    typedef vector<value_type,model_type::dimension> momentum_type;

    std::cout<<"-----------------------------------------------------"<<std::endl;
    std::cout<<"testing computed phase space variables..............."<<std::endl;
    std::cout<<"-----------------------------------------------------"<<std::endl;

    momentum_collection<model_type>ps(3);
    const int N_events=100000;

    {
	std::cerr<<"Checking invariant particle mass calculations..........";
	random_number_stream<value_type,std::random>rng;
	for(int i=0;i<N_events;++i)
	{
	    momentum_type& p1=ps.p_out(0);
	    momentum_type& p2=ps.p_out(1);
	    momentum_type& p3=ps.p_out(2);

	    value_type scale=rng(0.0001,10000);
	    generate_momenta(p1,p2,p3,scale);

	    ps.next_event();

	    value_type s1=ps.s_out(0);
	    value_type s1_check=p1[0]*p1[0]-p1[1]*p1[1]-p1[2]*p1[2]-p1[3]*p1[3];
	    if(!equals(s1,s1_check))
	    {
		std::cerr<<"incorrect invariant mass calculated for p1 "<<p1<<": "<<s1<<" not equal to checked value "<<s1_check<<std::endl;
		return 1;
	    }

	    value_type s2=ps.s_out(1);
	    value_type s2_check=p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]-p2[3]*p2[3];
	    if(!equals(s2,s2_check))
	    {
		std::cerr<<"incorrect invariant mass calculated for p2 "<<p2<<": "<<s2<<" not equal to checked value "<<s2_check<<std::endl;
		return 1;
	    }

	    value_type s3=ps.s_out(2);
	    value_type s3_check=p3[0]*p3[0]-p3[1]*p3[1]-p3[2]*p3[2]-p3[3]*p3[3];
	    if(!equals(s3,s3_check))
	    {
		std::cerr<<"incorrect invariant mass calculated for p3 "<<p3<<": "<<s3<<" not equal to checked value "<<s3_check<<std::endl;
		return 1;
	    }

	    value_type m1=ps.m_out(0);
	    value_type m1_check=std::sqrt(s1_check);
	    if(!equals(m1,m1_check))
	    {
		std::cerr<<"incorrect invariant mass calculated for p1 "<<p1<<": "<<m1<<" not equal to checked value "<<m1_check<<std::endl;
		return 1;
	    }

	    value_type m2=ps.m_out(1);
	    value_type m2_check=std::sqrt(s2_check);
	    if(!equals(m2,m2_check))
	    {
		std::cerr<<"incorrect invariant mass calculated for p2 "<<p2<<": "<<m2<<" not equal to checked value "<<m2_check<<std::endl;
		return 1;
	    }

	    value_type m3=ps.m_out(2);
	    value_type m3_check=std::sqrt(s3_check);
	    if(!equals(m3,m3_check))
	    {
		std::cerr<<"incorrect invariant mass calculated for p3 "<<p3<<": "<<m3<<" not equal to checked value "<<m3_check<<std::endl;
		return 1;
	    }

	    value_type msum=ps.ps_generator_base<model_type>::M_out();
	    value_type msum_check=m1+m2+m3;
	    if(!equals(msum,msum_check))
	    {
		std::cerr<<"incorrect total mass calculated for momenta "<<p1<<','<<p2<<','<<p3<<": value "<<msum<<" not equal to checked value "<<msum_check<<std::endl;
	    }

	    if(!equals(ps.s_tot(),scale*scale))
	    {
		std::cerr<<"incorrect total invariant mass calculated for momenta "<<p1<<','<<p2<<','<<p3<<": value "<<ps.s_tot()<<" not equal to checked value "<<scale*scale<<std::endl;
		return 1;
	    }

	    if(!equals(ps.s_hat(),scale*scale))
	    {
		std::cerr<<"incorrect partonic invariant mass calculated for momenta "<<p1<<','<<p2<<','<<p3<<": value "<<ps.s_hat()<<" not equal to checked value "<<scale*scale<<std::endl;
		return 1;
	    }

	    if(!equals(ps.Ecm_hat(),scale))
	    {
		std::cerr<<"incorrect total CM-energy calculated for momenta "<<p1<<','<<p2<<','<<p3<<": value "<<ps.Ecm_hat()<<" not equal to checked value "<<scale<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    {
	std::cerr<<"Checking particle momentum calculations..........";
	random_number_stream<value_type,std::random>rng;
	for(int i=0;i<N_events;++i)
	{
	    momentum_type& p1=ps.p_out(0);
	    momentum_type& p2=ps.p_out(1);
	    momentum_type& p3=ps.p_out(2);

	    value_type scale=rng(0.0001,10000);
	    generate_momenta(p1,p2,p3,scale);

	    ps.next_event();

	    value_type P12=p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3];
	    if(!equals(ps.P2(1),P12))
	    {
		std::cerr<<"incorrect momentum-squared calculated for momentum "<<p1<<": value "<<ps.P2(1)<<" not equal to checked value "<<P12<<std::endl;
		return 1;
	    }

	    value_type P1=std::sqrt(P12);
	    if(!equals(ps.P(1),P1))
	    {
		std::cerr<<"incorrect momentum-squared calculated for momentum "<<p1<<": value "<<ps.P(1)<<" not equal to checked value "<<P1<<std::endl;
		return 1;
	    }

	    value_type P22=p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3];
	    if(!equals(ps.P2(2),P22))
	    {
		std::cerr<<"incorrect momentum-squared calculated for momentum "<<p2<<": value "<<ps.P2(2)<<" not equal to checked value "<<P22<<std::endl;
		return 1;
	    }

	    value_type P2=std::sqrt(P22);
	    if(!equals(ps.P(2),P2))
	    {
		std::cerr<<"incorrect momentum-squared calculated for momentum "<<p2<<": value "<<ps.P(2)<<" not equal to checked value "<<P2<<std::endl;
		return 1;
	    }
	    
	    value_type P32=p3[1]*p3[1]+p3[2]*p3[2]+p3[3]*p3[3];
	    if(!equals(ps.P2(3),P32))
	    {
		std::cerr<<"incorrect momentum-squared calculated for momentum "<<p3<<": value "<<ps.P2(3)<<" not equal to checked value "<<P32<<std::endl;
		return 1;
	    }

	    value_type P3=std::sqrt(P32);
	    if(!equals(ps.P(3),P3))
	    {
		std::cerr<<"incorrect momentum-squared calculated for momentum "<<p3<<": value "<<ps.P(3)<<" not equal to checked value "<<P3<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    {
	std::cerr<<"Checking particle dimass calculations..........";
	random_number_stream<value_type,std::random>rng;
	for(int i=0;i<N_events;++i)
	{
	    momentum_type& p1=ps.p_out(0);
	    momentum_type& p2=ps.p_out(1);
	    momentum_type& p3=ps.p_out(2);

	    value_type scale=rng(0.0001,10000);
	    generate_momenta(p1,p2,p3,scale);

	    ps.next_event();

	    value_type s12=ps.s(1,2);
	    momentum_type p12=p1+p2;
	    value_type s12_check=p12[0]*p12[0]-p12[1]*p12[1]-p12[2]*p12[2]-p12[3]*p12[3];
	    if(!equals(s12,s12_check))
	    {
		std::cerr<<"incorrect dimass-squared calculated for momenta "<<p1<<','<<p2<<": value "<<s12<<" not equal to checked value "<<s12_check<<std::endl;
	    }

	    value_type m12=ps.m(1,2);
	    value_type m12_check=std::sqrt(s12_check);
	    if(!equals(m12,m12_check))
	    {
		std::cerr<<"incorrect dimass calculated for momenta "<<p1<<','<<p2<<": value "<<m12<<" not equal to checked value "<<m12_check<<std::endl;
	    }

	    value_type p1p2=ps.dot(1,2);
	    value_type p1p2_check=p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]-p1[3]*p2[3];
	    if(!equals(p1p2,p1p2_check))
	    {
		std::cerr<<"incorrect dot product calculated for momenta "<<p1<<','<<p2<<": value "<<p1p2<<" not equal to checked value "<<p1p2_check<<std::endl;
	    }

	    value_type v1v2=ps.space_dot(1,2);
	    value_type v1v2_check=p1[1]*p2[1]+p1[2]*p2[2]+p1[3]*p2[3];
	    if(!equals(v1v2,v1v2_check))
	    {
		std::cerr<<"incorrect dot product calculated for momenta "<<p1<<','<<p2<<": value "<<v1v2<<" not equal to checked value "<<v1v2_check<<std::endl;
	    }
	    
	    value_type s13=ps.s(1,3);
	    momentum_type p13=p1+p3;
	    value_type s13_check=p13[0]*p13[0]-p13[1]*p13[1]-p13[2]*p13[2]-p13[3]*p13[3];
	    if(!equals(s13,s13_check))
	    {
		std::cerr<<"incorrect dimass-squared calculated for momenta "<<p1<<','<<p3<<": value "<<s13<<" not equal to checked value "<<s13_check<<std::endl;
	    }

	    value_type m13=ps.m(1,3);
	    value_type m13_check=std::sqrt(s13_check);
	    if(!equals(m13,m13_check))
	    {
		std::cerr<<"incorrect dimass calculated for momenta "<<p1<<','<<p3<<": value "<<m13<<" not equal to checked value "<<m13_check<<std::endl;
	    }

	    value_type p1p3=ps.dot(1,3);
	    value_type p1p3_check=p1[0]*p3[0]-p1[1]*p3[1]-p1[2]*p3[2]-p1[3]*p3[3];
	    if(!equals(p1p3,p1p3_check))
	    {
		std::cerr<<"incorrect dot product calculated for momenta "<<p1<<','<<p3<<": value "<<p1p3<<" not equal to checked value "<<p1p3_check<<std::endl;
	    }

	    value_type v1v3=ps.space_dot(1,3);
	    value_type v1v3_check=p1[1]*p3[1]+p1[2]*p3[2]+p1[3]*p3[3];
	    if(!equals(v1v3,v1v3_check))
	    {
		std::cerr<<"incorrect dot product calculated for momenta "<<p1<<','<<p3<<": value "<<v1v3<<" not equal to checked value "<<v1v3_check<<std::endl;
	    }
	    
	    value_type s23=ps.s(2,3);
	    momentum_type p23=p2+p3;
	    value_type s23_check=p23[0]*p23[0]-p23[1]*p23[1]-p23[2]*p23[2]-p23[3]*p23[3];
	    if(!equals(s23,s23_check))
	    {
		std::cerr<<"incorrect dimass-squared calculated for momenta "<<p2<<','<<p3<<": value "<<s23<<" not equal to checked value "<<s23_check<<std::endl;
	    }

	    value_type m23=ps.m(2,3);
	    value_type m23_check=std::sqrt(s23_check);
	    if(!equals(m23,m23_check))
	    {
		std::cerr<<"incorrect dimass calculated for momenta "<<p2<<','<<p3<<": value "<<m23<<" not equal to checked value "<<m23_check<<std::endl;
	    }

	    value_type p2p3=ps.dot(2,3);
	    value_type p2p3_check=p2[0]*p3[0]-p2[1]*p3[1]-p2[2]*p3[2]-p2[3]*p3[3];
	    if(!equals(p2p3,p2p3_check))
	    {
		std::cerr<<"incorrect dot product calculated for momenta "<<p2<<','<<p3<<": value "<<p2p3<<" not equal to checked value "<<p2p3_check<<std::endl;
	    }

	    value_type v2v3=ps.space_dot(2,3);
	    value_type v2v3_check=p2[1]*p3[1]+p2[2]*p3[2]+p2[3]*p3[3];
	    if(!equals(v2v3,v2v3_check))
	    {
		std::cerr<<"incorrect dot product calculated for momenta "<<p2<<','<<p3<<": value "<<v2v3<<" not equal to checked value "<<v2v3_check<<std::endl;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    {
	std::cerr<<"Checking particle transverse mass calculations..........";
	random_number_stream<value_type,std::random>rng;
	for(int i=0;i<N_events;++i)
	{
	    momentum_type& p1=ps.p_out(0);
	    momentum_type& p2=ps.p_out(1);
	    momentum_type& p3=ps.p_out(2);

	    value_type scale=rng(0.0001,10000);
	    generate_momenta(p1,p2,p3,scale);

	    ps.next_event();

	    value_type pT12=ps.pT2(1);
	    value_type pT12_check=p1[1]*p1[1]+p1[3]*p1[3];
	    if(!equals(pT12,pT12_check))
	    {
		std::cerr<<"incorrect transverse mass-squared calculated for p1 "<<p1<<": "<<pT12<<" not equal to checked value "<<pT12_check<<std::endl;
		return 1;
	    }

	    value_type pT1=ps.pT(1);
	    value_type pT1_check=std::sqrt(pT12_check);
	    if(!equals(pT1,pT1_check))
	    {
		std::cerr<<"incorrect transverse mass calculated for p1 "<<p1<<": "<<pT1<<" not equal to checked value "<<pT1_check<<std::endl;
		return 1;
	    }

	    value_type pT22=ps.pT2(2);
	    value_type pT22_check=p2[1]*p2[1]+p2[3]*p2[3];
	    if(!equals(pT22,pT22_check))
	    {
		std::cerr<<"incorrect transverse mass-squared calculated for p2 "<<p2<<": "<<pT22<<" not equal to checked value "<<pT22_check<<std::endl;
		return 1;
	    }

	    value_type pT2=ps.pT(2);
	    value_type pT2_check=std::sqrt(pT22_check);
	    if(!equals(pT2,pT2_check))
	    {
		std::cerr<<"incorrect transverse mass calculated for p2 "<<p2<<": "<<pT2<<" not equal to checked value "<<pT2_check<<std::endl;
		return 1;
	    }

	    value_type pT32=ps.pT2(3);
	    value_type pT32_check=p3[1]*p3[1]+p3[3]*p3[3];
	    if(!equals(pT32,pT32_check))
	    {
		std::cerr<<"incorrect transverse mass-squared calculated for p3 "<<p3<<": "<<pT32<<" not equal to checked value "<<pT32_check<<std::endl;
		return 1;
	    }

	    value_type pT3=ps.pT(3);
	    value_type pT3_check=std::sqrt(pT32_check);
	    if(!equals(pT3,pT3_check))
	    {
		std::cerr<<"incorrect transverse mass calculated for p3 "<<p3<<": "<<pT3<<" not equal to checked value "<<pT3_check<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    {
	std::cerr<<"Checking particle (pseudo-)rapidity calculations..........";
	random_number_stream<value_type,std::random>rng;
	for(int i=0;i<N_events;++i)
	{
	    momentum_type& p1=ps.p_out(0);
	    momentum_type& p2=ps.p_out(1);
	    momentum_type& p3=ps.p_out(2);

	    value_type scale=rng(0.0001,10000);
	    generate_momenta(p1,p2,p3,scale);

	    ps.next_event();

	    value_type y1=ps.y(1);
	    value_type y1_check=0.5*std::log((p1[0]+p1[2])/(p1[0]-p1[2]));
	    if(!equals(y1,y1_check))
	    {
		std::cerr<<"incorrect rapidity calculated for p1 "<<p1<<": "<<y1<<" not equal to checked value "<<y1_check<<std::endl;
		return 1;
	    }

	    value_type eta1=ps.eta(1);
	    value_type pvec1=std::sqrt(p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3]);
	    value_type eta1_check=0.5*std::log((pvec1+p1[2])/(pvec1-p1[2]));
	    if(!equals(eta1,eta1_check))
	    {
		std::cerr<<"incorrect pseudo-rapidity calculated for p1 "<<p1<<": "<<eta1<<" not equal to checked value "<<eta1_check<<std::endl;
		return 1;
	    }

	    value_type y2=ps.y(2);
	    value_type y2_check=0.5*std::log((p2[0]+p2[2])/(p2[0]-p2[2]));
	    if(!equals(y2,y2_check))
	    {
		std::cerr<<"incorrect rapidity calculated for p2 "<<p2<<": "<<y2<<" not equal to checked value "<<y2_check<<std::endl;
		return 1;
	    }

	    value_type eta2=ps.eta(2);
	    value_type pvec2=std::sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3]);
	    value_type eta2_check=0.5*std::log((pvec2+p2[2])/(pvec2-p2[2]));
	    if(!equals(eta2,eta2_check))
	    {
		std::cerr<<"incorrect pseudo-rapidity calculated for p2 "<<p2<<": "<<eta2<<" not equal to checked value "<<eta2_check<<std::endl;
		return 1;
	    }

	    value_type y3=ps.y(3);
	    value_type y3_check=0.5*std::log((p3[0]+p3[2])/(p3[0]-p3[2]));
	    if(!equals(y3,y3_check))
	    {
		std::cerr<<"incorrect rapidity calculated for p3 "<<p3<<": "<<y3<<" not equal to checked value "<<y3_check<<std::endl;
		return 1;
	    }

	    value_type eta3=ps.eta(3);
	    value_type pvec3=std::sqrt(p3[1]*p3[1]+p3[2]*p3[2]+p3[3]*p3[3]);
	    value_type eta3_check=0.5*std::log((pvec3+p3[2])/(pvec3-p3[2]));
	    if(!equals(eta3,eta3_check))
	    {
		std::cerr<<"incorrect pseudo-rapidity calculated for p2 "<<p2<<": "<<eta2<<" not equal to checked value "<<eta2_check<<std::endl;
		return 1;
	    }

	    value_type dy12=ps.d_y(1,2);
	    value_type dy12_check=y1-y2;
	    if(!equals(dy12,dy12_check))
	    {
		std::cerr<<"incorrect rapidity difference calculated for p1,p2 "<<p1<<','<<p2<<": "<<dy12<<" not equal to checked value "<<dy12_check<<std::endl;
		return 1;
	    }

	    value_type deta12=ps.d_eta(1,2);
	    value_type deta12_check=eta1-eta2;
	    if(!equals(deta12,deta12_check))
	    {
		std::cerr<<"incorrect pseudo-rapidity difference calculated for p1,p2 "<<p1<<','<<p2<<": "<<deta12<<" not equal to checked value "<<deta12_check<<std::endl;
		return 1;
	    }

	    value_type dy13=ps.d_y(1,3);
	    value_type dy13_check=y1-y3;
	    if(!equals(dy13,dy13_check))
	    {
		std::cerr<<"incorrect rapidity difference calculated for p1,p3 "<<p1<<','<<p3<<": "<<dy13<<" not equal to checked value "<<dy13_check<<std::endl;
		return 1;
	    }

	    value_type deta13=ps.d_eta(1,3);
	    value_type deta13_check=eta1-eta3;
	    if(!equals(deta13,deta13_check))
	    {
		std::cerr<<"incorrect pseudo-rapidity difference calculated for p1,p3 "<<p1<<','<<p3<<": "<<deta13<<" not equal to checked value "<<deta13_check<<std::endl;
		return 1;
	    }

	    value_type dy23=ps.d_y(2,3);
	    value_type dy23_check=y2-y3;
	    if(!equals(dy23,dy23_check))
	    {
		std::cerr<<"incorrect rapidity difference calculated for p2,p3 "<<p2<<','<<p3<<": "<<dy23<<" not equal to checked value "<<dy23_check<<std::endl;
		return 1;
	    }

	    value_type deta23=ps.d_eta(2,3);
	    value_type deta23_check=eta2-eta3;
	    if(!equals(deta23,deta23_check))
	    {
		std::cerr<<"incorrect pseudo-rapidity difference calculated for p2,p3 "<<p2<<','<<p3<<": "<<deta23<<" not equal to checked value "<<deta23_check<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }

    {
	std::cerr<<"Checking particle angle calculations..........";
	random_number_stream<value_type,std::random>rng;
	for(int i=0;i<N_events;++i)
	{
	    momentum_type& p1=ps.p_out(0);
	    momentum_type& p2=ps.p_out(1);
	    momentum_type& p3=ps.p_out(2);

	    value_type scale=rng(0.0001,10000);
	    generate_momenta(p1,p2,p3,scale);

	    ps.next_event();

	    value_type ct1=ps.cos_theta(1);
	    value_type p1vec=p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3];
	    value_type ct1_check=p1[2]/std::sqrt(p1vec);
	    if(!equals(ct1,ct1_check))
	    {
		std::cerr<<"incorrect cosine theta calculated for p1 "<<p1<<": "<<ct1<<" not equal to checked value "<<ct1_check<<std::endl;
		return 1;
	    }

	    value_type theta1=ps.theta(1);
	    value_type theta1_check=std::acos(ct1_check);
	    if(!equals(theta1,theta1_check))
	    {
		std::cerr<<"incorrect theta calculated for p1 "<<p1<<": "<<theta1<<" not equal to checked value "<<theta1_check<<std::endl;
		return 1;
	    }

	    value_type phi1=ps.phi(1);
	    value_type phi1_check=std::atan2(p1[1],p1[3]);
	    if(!equals(phi1,phi1_check))
	    {
		std::cerr<<"incorrect phi calculated for p1 "<<p1<<": "<<phi1<<" not equal to checked value "<<phi1_check<<std::endl;
		return 1;
	    }

	    value_type ct2=ps.cos_theta(2);
	    value_type p2vec=p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3];
	    value_type ct2_check=p2[2]/std::sqrt(p2vec);
	    if(!equals(ct2,ct2_check))
	    {
		std::cerr<<"incorrect cosine theta calculated for p2 "<<p2<<": "<<ct2<<" not equal to checked value "<<ct2_check<<std::endl;
		return 1;
	    }

	    value_type theta2=ps.theta(2);
	    value_type theta2_check=std::acos(ct2_check);
	    if(!equals(theta2,theta2_check))
	    {
		std::cerr<<"incorrect theta calculated for p2 "<<p2<<": "<<theta2<<" not equal to checked value "<<theta2_check<<std::endl;
		return 1;
	    }

	    value_type phi2=ps.phi(2);
	    value_type phi2_check=std::atan2(p2[1],p2[3]);
	    if(!equals(phi2,phi2_check))
	    {
		std::cerr<<"incorrect phi calculated for p2 "<<p2<<": "<<phi2<<" not equal to checked value "<<phi2_check<<std::endl;
		return 1;
	    }

	    value_type ct3=ps.cos_theta(3);
	    value_type p3vec=p3[1]*p3[1]+p3[2]*p3[2]+p3[3]*p3[3];
	    value_type ct3_check=p3[2]/std::sqrt(p3vec);
	    if(!equals(ct3,ct3_check))
	    {
		std::cerr<<"incorrect cosine theta calculated for p3 "<<p3<<": "<<ct3<<" not equal to checked value "<<ct3_check<<std::endl;
		return 1;
	    }

	    value_type theta3=ps.theta(3);
	    value_type theta3_check=std::acos(ct3_check);
	    if(!equals(theta3,theta3_check))
	    {
		std::cerr<<"incorrect theta calculated for p3 "<<p3<<": "<<theta3<<" not equal to checked value "<<theta3_check<<std::endl;
		return 1;
	    }

	    value_type phi3=ps.phi(3);
	    value_type phi3_check=std::atan2(p3[1],p3[3]);
	    if(!equals(phi3,phi3_check))
	    {
		std::cerr<<"incorrect phi calculated for p3 "<<p3<<": "<<phi3<<" not equal to checked value "<<phi3_check<<std::endl;
		return 1;
	    }

	    value_type P1=p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3];
	    value_type P2=p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3];
	    value_type P3=p3[1]*p3[1]+p3[2]*p3[2]+p3[3]*p3[3];

	    value_type cos_alpha12=ps.cos_alpha(1,2);
	    value_type cos_alpha12_check=(p1[1]*p2[1]+p1[2]*p2[2]+p1[3]*p2[3])/std::sqrt(P1*P2);
	    if(!equals(cos_alpha12,cos_alpha12_check))
	    {
		std::cerr<<"incorrect cosine angle calculated for p1 "<<p1<<", p2 "<<p2<<": "<<cos_alpha12<<" not equal to checked value "<<cos_alpha12_check<<std::endl;
		return 1;
	    }

	    value_type alpha12=ps.alpha(1,2);
	    value_type alpha12_check=std::acos(cos_alpha12_check);
	    if(!equals(alpha12,alpha12_check))
	    {
		std::cerr<<"incorrect angle calculated for p1 "<<p1<<", p2 "<<p2<<": "<<alpha12<<" not equal to checked value "<<alpha12_check<<std::endl;
		return 1;
	    }

	    value_type cos_alpha13=ps.cos_alpha(1,3);
	    value_type cos_alpha13_check=(p1[1]*p3[1]+p1[2]*p3[2]+p1[3]*p3[3])/std::sqrt(P1*P3);
	    if(!equals(cos_alpha13,cos_alpha13_check))
	    {
		std::cerr<<"incorrect cosine angle calculated for p1 "<<p1<<", p3 "<<p3<<": "<<cos_alpha13<<" not equal to checked value "<<cos_alpha13_check<<std::endl;
		return 1;
	    }

	    value_type alpha13=ps.alpha(1,3);
	    value_type alpha13_check=std::acos(cos_alpha13_check);
	    if(!equals(alpha13,alpha13_check))
	    {
		std::cerr<<"incorrect angle calculated for p1 "<<p1<<", p3 "<<p3<<": "<<alpha13<<" not equal to checked value "<<alpha13_check<<std::endl;
		return 1;
	    }

	    value_type cos_alpha23=ps.cos_alpha(2,3);
	    value_type cos_alpha23_check=(p2[1]*p3[1]+p2[2]*p3[2]+p2[3]*p3[3])/std::sqrt(P2*P3);
	    if(!equals(cos_alpha23,cos_alpha23_check))
	    {
		std::cerr<<"incorrect cosine angle calculated for p2 "<<p2<<", p3 "<<p3<<": "<<cos_alpha23<<" not equal to checked value "<<cos_alpha23_check<<std::endl;
		return 1;
	    }

	    value_type alpha23=ps.alpha(2,3);
	    value_type alpha23_check=std::acos(cos_alpha23_check);
	    if(!equals(alpha23,alpha23_check))
	    {
		std::cerr<<"incorrect angle calculated for p2 "<<p1<<", p3 "<<p3<<": "<<alpha23<<" not equal to checked value "<<alpha23_check<<std::endl;
		return 1;
	    }
	}
	std::cerr<<".........done."<<std::endl;
    }
}

