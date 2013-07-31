//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SM_BASE_H_ 
#define CAMGEN_SM_BASE_H_

#include <Camgen/scalar_particle.h>
#include <Camgen/vector_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/vff.h>
#include <Camgen/sff.h>
#include <Camgen/vffVA.h>
#include <Camgen/vffR.h>
#include <Camgen/svv.h>
#include <Camgen/svv3.h>
#include <Camgen/vss.h>
#include <Camgen/ssvv.h>
#include <Camgen/ssvv3.h>
#include <Camgen/svvv3.h>
#include <Camgen/sss.h>
#include <Camgen/ssss.h>
#include <Camgen/sffL.h>
#include <Camgen/sffR.h>
#include <Camgen/sff5.h>
#include <Camgen/sffVA.h>
#include <Camgen/vvv.h>
#include <Camgen/vvvv.h>
#include <Camgen/su(n).h>
#include <Camgen/f.h>
#include <Camgen/ff_contr.h>
#include <Camgen/T.h>
#include <Camgen/d.h>
#include <Camgen/model.h>
#include <Camgen/SM_params.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the standard model base class in Camgen.  Deriving a class    *
 * model_t with value_type value_t from SM_base<model_t,value_t> inserts all the *
 * standard model particles and couplings automatically in your model. The       *
 * coupling constant values may be values may be inserted manually, or computed  *
 * by the class from the basic input variables alpha, alpha_s, QCD scale, Fermi  *
 * constant G_F, the Z-boson mass and the top and bottom masses. Gauge-invariant *
 * inclusion of particle decay widths is implemented via the complex mass        *
 * scheme, replacing all mass parameters in the computed couplings by the        *
 * (complex) square root of                                                      *
 *                                                                               *
 * 				M^2-iM\Gamma                                     *
 *                                                                               *
 * The decay width can again be manually inserted, or computed by the class      *
 * methods. Note that all computations are of leading order in electroweak and   *
 * strong couplings.                                                             *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,class value_t>class SM_base: public model<model_t>
    {
	public:

	    /* Value type definition: */

	    typedef value_t value_type;

	    /* Pi definition: */

	    static const value_type pi;

	    /* Particle masses: */

	    static value_type M_e,M_mu,M_tau,M_u,M_d,M_c,M_s,M_t,M_b,M_Z,M_W,M_h0;

	    /* Particle decay widths: */

	    static value_type W_t,W_Z,W_W,W_h0;
	    
	    /* Fine structure constant: */

	    static value_type alpha;

	    /* Alpha strong: */

	    static value_type alpha_s;
	    
	    /* QCD scale: */
	    
	    static value_type QCD_scale;

	    /* Fermi's constant: */

	    static value_type G_F;

	    /* CKM sines of Euler angles: */

	    static value_type s12,s23,s13;

	    /* CP-violating phase: */

	    static value_type delta;
	    
	    /* Quadruple-vector couplings: */

	    static std::complex<value_type> WWWW,WZWZ,WgammaWZ,WgammaWgamma,gggg;

	    /* Triple-vector couplings: */

	    static std::complex<value_type> gammaWW,ZWW,ggg;

	    /* Quadruple-scalar couplings: */

	    static std::complex<value_type> hhhh,hhchichi,hhphiphi,chichichichi,chichiphiphi,phiphiphiphi;

	    /* Triple-scalar couplings: */

	    static std::complex<value_type> hhh,hchichi,hphiphi;

	    /* Vector-fermion-fermion couplings: */

	    static std::complex<value_type> gammaee,gammamumu,gammatautau,gammauu,gammacc,gammatt,gammadd,gammass,gammabb;

	    static std::complex<value_type> Znene,Znmunmu,Zntauntau,Zee_V,Zee_A,Zmumu_V,Zmumu_A,Ztautau_V,Ztautau_A,Zuu_V,Zuu_A,Zdd_V,Zdd_A,Zcc_V,Zcc_A,Zss_V,Zss_A,Ztt_V,Ztt_A,Zbb_V,Zbb_A;

	    static std::complex<value_type> Wnee,Wene,Wnmumu,Wmunmu,Wntautau,Wtauntau,Wud,Wdu,Wus,Wsu,Wub,Wbu,Wcd,Wdc,Wcs,Wsc,Wcb,Wbc,Wtd,Wdt,Wts,Wst,Wtb,Wbt;
	    
	    static std::complex<value_type> guu,gdd,gcc,gss,gtt,gbb;

	    /* Scalar-fermion-fermion couplings: */

	    static std::complex<value_type> hee,hmumu,htautau,huu,hdd,hcc,hss,htt,hbb;

	    static std::complex<value_type> chiee,chimumu,chitautau,chiuu,chidd,chicc,chiss,chitt,chibb;

	    static std::complex<value_type> phinee,phiene,phinmumu,phimunmu,phintautau,phitauntau,phiud_S,phiud_A,phidu_S,phidu_A,phius_S,phius_A,phisu_S,phisu_A,phiub_S,phiub_A,phibu_S,phibu_A,phicd_S,phicd_A,phidc_S,phidc_A,phics_S,phics_A,phisc_S,phisc_A,phicb_S,phicb_A,phibc_S,phibc_A,phitd_S,phitd_A,phidt_S,phidt_A,phits_S,phits_A,phist_S,phist_A,phitb_S,phitb_A,phibt_S,phibt_A;

	    /* Vector-scalar-scalar couplings: */

	    static std::complex<value_type> gammachih,Zchih,gammaphiphi,Zphiphi,Wpphimh,Wmphiph,Wpphimchi,Wmphipchi;

	    /* Scalar-vector-vector couplings: */

	    static std::complex<value_type> hWW,hZZ,hZgamma,phipWmZ,phimWpZ,phipWmgamma,phimWpgamma;
	    
	    /* Scalar-scalar-vector-vector couplings: */

	    static std::complex<value_type> hhWW,chichiWW,phiphiWW,phiphiZZ,phiphigammaZ,phiphigammagamma,hhZZ,chichiZZ,phiphWmZ,phimhWpZ,phiphWmgamma,phimhWpgamma,phimchiWpZ,phipchiWmZ,phimchiWpgamma,phipchiWmgamma;

	    /* Tensor-vector-vector couplings: */

	    static std::complex<value_type>Tgg;

	    /* Constructor declaration: */

	    SM_base()
	    {
		/* Lepton definitions: */

		(M_e>(value_type)0)?model<model_t>::add_fermions("e-","e+",&M_e,11):model<model_t>::add_fermions("e-","e+",11);
		(M_mu>(value_type)0)?model<model_t>::add_fermions("mu-","mu+",&M_mu,13):model<model_t>::add_fermions("mu-","mu+",13);
		(M_tau>(value_type)0)?model<model_t>::add_fermions("tau-","tau+",&M_tau,15):model<model_t>::add_fermions("tau-","tau+",15);
		model<model_t>::add_fermions("nu_e","nu_ebar",12);
		model<model_t>::add_fermions("nu_mu","nu_mubar",14);
		model<model_t>::add_fermions("nu_tau","nu_taubar",16);

		/* Quark definitions: */

		(M_d>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("d","dbar",&M_d,1):model<model_t>::template add_quarks< SU<model_t::N_c> >("d","dbar",1);
		(M_u>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("u","ubar",&M_u,2):model<model_t>::template add_quarks< SU<model_t::N_c> >("u","ubar",2);
		(M_s>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("s","sbar",&M_s,3):model<model_t>::template add_quarks< SU<model_t::N_c> >("s","sbar",3);
		(M_c>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("c","cbar",&M_c,4):model<model_t>::template add_quarks< SU<model_t::N_c> >("c","cbar",4);
		(M_b>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("b","bbar",&M_b,5):model<model_t>::template add_quarks< SU<model_t::N_c> >("b","bbar",5);
		(M_t>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("t","tbar",&M_t,&W_t,6):model<model_t>::template add_quarks< SU<model_t::N_c> >("t","tbar",6);

		/* Gauge field insertions: */

		model<model_t>::template add_vector<Feynman_gauge>("gamma",22);
		model<model_t>::template add_vectors<Feynman_gauge>("W+","W-",&M_W,&W_W,24);
		model<model_t>::template add_vector<Feynman_gauge>("Z",&M_Z,&W_Z,23);
		model<model_t>::template add_gluon<Feynman_gauge,SU<model_t::N_c> >();

		/* Scalar particle insertions: */

		model<model_t>::add_scalar("h0",&M_h0,&W_h0,25);
		model<model_t>::add_scalar("chi",&M_Z,&W_Z);
		model<model_t>::add_scalars("phi+","phi-",&M_W,&W_W);

		/* 4-vector vertices: */

		model<model_t>::template add_vertex<vvvv>("W+","W-","W+","W-",&WWWW);
		model<model_t>::template add_vertex<vvvv>("W+","Z","W-","Z",&WZWZ);
		model<model_t>::template add_vertex<vvvv>("W+","gamma","W-","Z",&WgammaWZ);
		model<model_t>::template add_vertex<vvvv>("W+","gamma","W-","gamma",&WgammaWgamma);

		if(four_gluon_vertex)
		{
		    model<model_t>::template add_vertex<colour_tensor::ff_contr< SU<model_t::N_c> >,vvvv>("g","g","g","g",&gggg);
		}
		else
		{
		    model<model_t>::template add_fast_4g_vertex< SU<model_t::N_c> >("g",&Tgg);
		}

		/* 3-vector vertices: */

		model<model_t>::template add_vertex<vvv>("gamma","W+","W-",&gammaee);
		model<model_t>::template add_vertex<vvv>("Z","W+","W-",&ZWW);
		model<model_t>::template add_vertex<colour_tensor::f< SU<model_t::N_c> >,vvv>("g","g","g",&ggg);

		/* 4-scalar vertices: */

		model<model_t>::template add_vertex<ssss>("h0","h0","h0","h0",&hhhh);
		model<model_t>::template add_vertex<ssss>("h0","h0","chi","chi",&hhchichi);
		model<model_t>::template add_vertex<ssss>("h0","h0","phi+","phi-",&hhphiphi);
		model<model_t>::template add_vertex<ssss>("chi","chi","chi","chi",&chichichichi);
		model<model_t>::template add_vertex<ssss>("chi","chi","phi+","phi-",&chichiphiphi);
		model<model_t>::template add_vertex<ssss>("phi+","phi-","phi+","phi-",&phiphiphiphi);

		/* 3-scalar vertices: */

		model<model_t>::template add_vertex<sss>("h0","h0","h0",&hhh);
		model<model_t>::template add_vertex<sss>("h0","chi","chi",&hchichi);
		model<model_t>::template add_vertex<sss>("h0","phi+","phi-",&hphiphi);

		/* Vector-fermion-fermion vertices: */

		model<model_t>::template add_vertex<vff>("gamma","e+","e-",&gammaee);
		model<model_t>::template add_vertex<vff>("gamma","mu+","mu-",&gammamumu);
		model<model_t>::template add_vertex<vff>("gamma","tau+","tau-",&gammatautau);

		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vff>("gamma","ubar","u",&gammauu);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vff>("gamma","dbar","d",&gammadd);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vff>("gamma","cbar","c",&gammacc);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vff>("gamma","sbar","s",&gammass);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vff>("gamma","tbar","t",&gammatt);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vff>("gamma","bbar","b",&gammabb);

		model<model_t>::template add_vertex<vffVA>("Z","e+","e-",&Zee_V,&Zee_A);
		model<model_t>::template add_vertex<vffR>("Z","nu_ebar","nu_e",&Znene);
		model<model_t>::template add_vertex<vffVA>("Z","mu+","mu-",&Zmumu_V,&Zmumu_A);
		model<model_t>::template add_vertex<vffR>("Z","nu_mubar","nu_mu",&Znmunmu);
		model<model_t>::template add_vertex<vffVA>("Z","tau+","tau-",&Ztautau_V,&Ztautau_A);
		model<model_t>::template add_vertex<vffR>("Z","nu_taubar","nu_tau",&Zntauntau);

		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffVA>("Z","ubar","u",&Zuu_V,&Zuu_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffVA>("Z","dbar","d",&Zdd_V,&Zdd_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffVA>("Z","cbar","c",&Zcc_V,&Zcc_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffVA>("Z","sbar","s",&Zss_V,&Zss_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffVA>("Z","tbar","t",&Ztt_V,&Ztt_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffVA>("Z","bbar","b",&Zbb_V,&Zbb_A);

		model<model_t>::template add_vertex<vffR>("W+","nu_ebar","e-",&Wnee);
		model<model_t>::template add_vertex<vffR>("W-","e+","nu_e",&Wene);
		model<model_t>::template add_vertex<vffR>("W+","nu_mubar","mu-",&Wnmumu);
		model<model_t>::template add_vertex<vffR>("W-","mu+","nu_mu",&Wmunmu);
		model<model_t>::template add_vertex<vffR>("W+","nu_taubar","tau-",&Wntautau);
		model<model_t>::template add_vertex<vffR>("W-","tau+","nu_tau",&Wtauntau);

		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","ubar","d",&Wud);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","dbar","u",&Wdu);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","ubar","s",&Wus);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","sbar","u",&Wsu);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","ubar","b",&Wub);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","bbar","u",&Wbu);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","cbar","d",&Wcd);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","dbar","c",&Wdc);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","cbar","s",&Wcs);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","sbar","c",&Wsc);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","cbar","b",&Wcb);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","bbar","c",&Wbc);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","tbar","d",&Wtd);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","dbar","t",&Wdt);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","tbar","s",&Wts);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","sbar","t",&Wst);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W+","tbar","b",&Wtb);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,vffR>("W-","bbar","t",&Wbt);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","ubar","u",&guu);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","dbar","d",&gdd);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","cbar","c",&gcc);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","sbar","s",&gss);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","tbar","t",&gtt);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","bbar","b",&gbb);

		/* Scalar-fermion-fermion vertices: */

		model<model_t>::template add_vertex<sff>("h0","e+","e-",&hee);
		model<model_t>::template add_vertex<sff>("h0","mu+","mu-",&hmumu);
		model<model_t>::template add_vertex<sff>("h0","tau+","tau-",&htautau);

		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff>("h0","ubar","u",&huu);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff>("h0","dbar","d",&hdd);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff>("h0","cbar","c",&hcc);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff>("h0","sbar","s",&hss);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff>("h0","tbar","t",&htt);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff>("h0","bbar","b",&hbb);

		model<model_t>::template add_vertex<sff5>("chi","e+","e-",&chiee);
		model<model_t>::template add_vertex<sff5>("chi","mu+","mu-",&chimumu);
		model<model_t>::template add_vertex<sff5>("chi","tau+","tau-",&chitautau);

		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff5>("chi","ubar","u",&chiuu);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff5>("chi","dbar","d",&chidd);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff5>("chi","cbar","c",&chicc);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff5>("chi","sbar","s",&chiss);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff5>("chi","tbar","t",&chitt);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sff5>("chi","bbar","b",&chibb);

		model<model_t>::template add_vertex<sffL>("phi+","nu_ebar","e-",&phinee);
		model<model_t>::template add_vertex<sffR>("phi-","e+","nu_e",&phiene);
		model<model_t>::template add_vertex<sffL>("phi+","nu_mubar","mu-",&phinmumu);
		model<model_t>::template add_vertex<sffR>("phi-","mu+","nu_mu",&phimunmu);
		model<model_t>::template add_vertex<sffL>("phi+","nu_taubar","tau-",&phintautau);
		model<model_t>::template add_vertex<sffR>("phi-","tau+","nu_tau",&phitauntau);

		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","ubar","d",&phiud_S,&phiud_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","dbar","u",&phidu_S,&phidu_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","ubar","s",&phius_S,&phius_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","sbar","u",&phisu_S,&phisu_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","ubar","b",&phiub_S,&phiub_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","bbar","u",&phibu_S,&phibu_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","cbar","d",&phicd_S,&phicd_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","dbar","c",&phidc_S,&phidc_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","cbar","s",&phics_S,&phics_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","sbar","c",&phisc_S,&phisc_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","cbar","b",&phicb_S,&phicb_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","bbar","c",&phibc_S,&phibc_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","tbar","d",&phitd_S,&phitd_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","dbar","t",&phidt_S,&phidt_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","tbar","s",&phits_S,&phits_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","sbar","t",&phist_S,&phist_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi+","tbar","b",&phitb_S,&phitb_A);
		model<model_t>::template add_vertex<colour_tensor::d<fundamental_rep< SU<model_t::N_c> >,1,2>,sffVA>("phi-","bbar","t",&phibt_S,&phibt_A);

		/* Vector-scalar-scalar vertices: */

		model<model_t>::template add_vertex<vss>("Z","chi","h0",&Zchih);
		model<model_t>::template add_vertex<vss>("gamma","phi+","phi-",&gammaphiphi);
		model<model_t>::template add_vertex<vss>("Z","phi+","phi-",&Zphiphi);
		model<model_t>::template add_vertex<vss>("W+","phi-","h0",&Wpphimh);
		model<model_t>::template add_vertex<vss>("W-","phi+","h0",&Wmphiph);
		model<model_t>::template add_vertex<vss>("W+","phi-","chi",&Wpphimchi);
		model<model_t>::template add_vertex<vss>("W-","phi+","chi",&Wmphipchi);

		/* Scalar-vector-vector vertex insertions: */

		model<model_t>::template add_vertex<svv>("h0","W+","W-",&hWW);
		model<model_t>::template add_vertex<svv>("h0","Z","Z",&hZZ);
		model<model_t>::template add_vertex<svv>("phi+","W-","Z",&phipWmZ);
		model<model_t>::template add_vertex<svv>("phi-","W+","Z",&phimWpZ);
		model<model_t>::template add_vertex<svv>("phi+","W-","gamma",&phipWmgamma);
		model<model_t>::template add_vertex<svv>("phi-","W+","gamma",&phimWpgamma);
		
		/* Scalar-scalar-vector-vector vertex insertions: */

		model<model_t>::template add_vertex<ssvv>("h0","h0","W+","W-",&hhWW);
		model<model_t>::template add_vertex<ssvv>("h0","h0","Z","Z",&hhZZ);
		model<model_t>::template add_vertex<ssvv>("phi+","phi-","W+","W-",&phiphiWW);
		model<model_t>::template add_vertex<ssvv>("phi+","phi-","Z","Z",&phiphiZZ);
		model<model_t>::template add_vertex<ssvv>("chi","chi","W+","W-",&chichiWW);
		model<model_t>::template add_vertex<ssvv>("chi","chi","Z","Z",&chichiZZ);
		model<model_t>::template add_vertex<ssvv>("phi+","phi-","gamma","Z",&phiphigammaZ);
		model<model_t>::template add_vertex<ssvv>("phi+","phi-","gamma","gamma",&phiphigammagamma);
		model<model_t>::template add_vertex<ssvv>("phi+","h0","W-","Z",&phiphWmZ);
		model<model_t>::template add_vertex<ssvv>("phi-","h0","W+","Z",&phimhWpZ);
		model<model_t>::template add_vertex<ssvv>("phi+","h0","W-","gamma",&phiphWmgamma);
		model<model_t>::template add_vertex<ssvv>("phi-","h0","W+","gamma",&phimhWpgamma);
		model<model_t>::template add_vertex<ssvv>("phi-","chi","W+","Z",&phimchiWpZ);
		model<model_t>::template add_vertex<ssvv>("phi+","chi","W-","Z",&phipchiWmZ);
		model<model_t>::template add_vertex<ssvv>("phi-","chi","W+","gamma",&phimchiWpgamma);
		model<model_t>::template add_vertex<ssvv>("phi+","chi","W-","gamma",&phipchiWmgamma);

		model<model_t>::decouple_zero_vertices();

		if(gauge==1)
		{
		    model<model_t>::decouple_particle("phi+");
		    model<model_t>::decouple_particle("chi");
		    model<model_t>::template set_propagator<unitary_gauge>("W+");
		    model<model_t>::template set_propagator<unitary_gauge>("W-");
		    model<model_t>::template set_propagator<unitary_gauge>("Z");
		    model<model_t>::template set_propagator<unitary_gauge>("gamma");
		    model<model_t>::template set_propagator<unitary_gauge>("g");
		}
		if(gauge==2)
		{
		    model<model_t>::template set_propagator<R_vector_gauge>("W+");
		    model<model_t>::template set_propagator<R_vector_gauge>("W-");
		    model<model_t>::template set_propagator<R_vector_gauge>("Z");
		    model<model_t>::template set_propagator<R_vector_gauge>("gamma");
		    model<model_t>::template set_propagator<R_vector_gauge>("g");
		    model<model_t>::template set_propagator<R_scalar_gauge>("phi+");
		    model<model_t>::template set_propagator<R_scalar_gauge>("phi-");
		    model<model_t>::template set_propagator<R_scalar_gauge>("chi");
		}

		/* Positively-charged lepton family definition: */

		model<model_t>::construct_family("l+","e+,mu+");

		/* Positively-charged lepton family definition, including taus: */

		model<model_t>::construct_family("L+","e+,mu+,tau+");

		/* Negatively-charged lepton family definition: */

		model<model_t>::construct_family("l-","e-,mu-");

		/* Negatively-charged lepton family definition, including taus: */

		model<model_t>::construct_family("L-","e-,mu-,tau-");

		/* Neutrino family definition: */

		model<model_t>::construct_family("nu","nu_e,nu_mu,nu_tau");

		/* Antineutrino family definition: */

		model<model_t>::construct_family("nubar","nu_ebar,nu_mubar,nu_taubar");

		/* Massive gauge boson definition: */

		model<model_t>::construct_family("V","Z,W+,W-");

		/* Quark family definition: */

		model<model_t>::construct_family("q","u,d,c,s");

		/* Quark family definition including b's: */

		model<model_t>::construct_family("Q","u,d,c,s,b");

		/* Quark family definition including b's and tops: */

		model<model_t>::construct_family("Qt","u,d,c,s,b,t");

		/* Antiquark family definition: */

		model<model_t>::construct_family("qbar","dbar,ubar,sbar,cbar");

		/* Antiquark family definition, including b's: */

		model<model_t>::construct_family("Qbar","dbar,ubar,sbar,cbar,bbar");

		/* Antiquark family definition, including b's and tops: */

		model<model_t>::construct_family("Qtbar","dbar,ubar,sbar,cbar,bbar,tbar");

		/* Up-type quark family definition: */

		model<model_t>::construct_family("q_up","u,c");

		/* Up-type quark family definition, including b's: */

		model<model_t>::construct_family("Q_up","u,c");

		/* Up-type quark family definition, including b's and tops: */

		model<model_t>::construct_family("Qt_up","u,c,t");

		/* Up-type antiquark family definition: */

		model<model_t>::construct_family("qbar_up","ubar,cbar");

		/* Up-type antiquark family definition, including b's: */

		model<model_t>::construct_family("Qbar_up","ubar,cbar");

		/* Up-type antiquark family definition, including b's and tops: */

		model<model_t>::construct_family("Qtbar_up","ubar,cbar,tbar");

		/* Down-type quark family definition: */

		model<model_t>::construct_family("q_down","d,s");

		/* Down-type quark family definition, including b's: */

		model<model_t>::construct_family("Q_down","d,s,b");

		/* Down-type quark family definition, including b's and tops: */

		model<model_t>::construct_family("Qt_down","d,s,b");

		/* Down-type antiquark family definition: */

		model<model_t>::construct_family("qbar_down","dbar,sbar");

		/* Down-type antiquark family definition, including b's: */

		model<model_t>::construct_family("Qbar_down","dbar,sbar,bbar");

		/* Down-type antiquark family definition, including b's and tops: */

		model<model_t>::construct_family("Qtbar_down","dbar,sbar,bbar");
		
		/* Positively-charged quark family definition: */

		model<model_t>::construct_family("q+","u,dbar,c,sbar");
		
		/* Positively-charged quark family definition including b's: */

		model<model_t>::construct_family("Q+","u,dbar,c,sbar,bbar");
		
		/* Positively-charged quark family definition including b's and
		 * tops: */

		model<model_t>::construct_family("Qt+","u,dbar,c,sbar,bbar,t");

		/* Negatively-charged quark family definition: */

		model<model_t>::construct_family("q-","d,ubar,s,cbar");

		/* Negatively-charged quark family definition, including b's: */

		model<model_t>::construct_family("Q-","d,ubar,s,cbar,b");

		/* Negatively-charged quark family definition, including b's and
		 * tops: */

		model<model_t>::construct_family("Qt-","d,ubar,s,cbar,b,tbar");

		/* Proton constituent parton family definition: */

		model<model_t>::construct_family("p","u,ubar,d,dbar,c,cbar,s,sbar,g");

		/* Proton constituent parton family definition, including b's: */

		model<model_t>::construct_family("P","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

		/* QCD jet parton family definition: */

		model<model_t>::construct_family("j","u,ubar,d,dbar,c,cbar,s,sbar,g");

		/* QCD jet parton family definition, including b's: */

		model<model_t>::construct_family("J","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

		/* QCD jet parton family definition, including b's and tops: */

		model<model_t>::construct_family("Jt","u,ubar,d,dbar,c,cbar,s,sbar,t,tbar,b,bbar,g");
	    }
	    
	    /* Output of the number of QCD colours: */

	    static std::size_t QCD_colours()
	    {
		return model_t::N_c;
	    }

	    /* Function computing the couplings from given parameters: */

	    static void refresh_couplings()
	    {
		refresh_weak_couplings();
		refresh_strong_couplings();
	    }

	    /* Function computing the couplings from given parameters: */

	    static void refresh_strong_couplings()
	    {
		g_s=std::sqrt((value_type)4*pi*alpha_s);

		ggg=std::complex<value_type>(g_s,0);
		gggg=std::complex<value_type>(0,-g_s*g_s);
		Tgg=std::complex<value_type>(0,std::sqrt((value_type)0.5)*g_s);
		guu=std::complex<value_type>(0,g_s);
		gdd=guu;
		gcc=guu;
		gss=guu;
		gtt=guu;
		gbb=guu;
	    }

	    /* Function computing the couplings from given parameters: */

	    static void refresh_weak_couplings()
	    {
		Q_e=std::sqrt((value_type)4*pi*alpha);
		cM_Z=std::sqrt(std::complex<value_type>(M_Z*M_Z,-M_Z*W_Z));
		cM_W=std::sqrt(std::complex<value_type>(M_W*M_W,-M_W*W_W));
		cos_W=cM_W/cM_Z;
		sin_W=std::sqrt(std::complex<value_type>(1,0)-cos_W*cos_W);
		tan_W=sin_W/cos_W;
		cos2_W=cos_W*cos_W-sin_W*sin_W;
		sin2_W=(value_type)2*sin_W*cos_W;
		refresh_CKM_matrix();

		/* VVVV couplings: */

		WWWW=std::complex<value_type>(0,Q_e*Q_e)/(sin_W*sin_W);
		WZWZ=-WWWW*cos_W*cos_W;
		WgammaWZ=WWWW*cos_W*sin_W;
		WgammaWgamma=-WWWW*sin_W*sin_W;

		/* VVV couplings: */

		gammaWW=std::complex<value_type>(0,Q_e);
		ZWW=-gammaWW/tan_W;

		/* SSSS couplings: */

		hhhh=-(value_type)0.75*std::complex<value_type>(0,Q_e*Q_e)*std::pow(M_h0/(sin_W*cM_W),2);
		hhchichi=hhhh/(value_type)3;
		hhphiphi=hhchichi;
		chichichichi=hhhh;
		chichiphiphi=hhphiphi;
		phiphiphiphi=(value_type)2*hhchichi;

		/* SSS couplings: */

		hhh=-(value_type)1.5*std::complex<value_type>(0,Q_e)*M_h0*M_h0/(sin_W*cM_W);
		hchichi=hhh/(value_type)3;
		hphiphi=hchichi;

		/* VFF couplings: */

		gammaee=std::complex<value_type>(0,Q_e);
		gammamumu=gammaee;
		gammatautau=gammaee;
		gammauu=std::complex<value_type>(0,-(value_type)2*Q_e/(value_type)3);
		gammadd=std::complex<value_type>(0,Q_e/(value_type)3);
		gammacc=gammauu;
		gammass=gammadd;
		gammatt=gammauu;
		gammabb=gammadd;

		Znene=gammaee/sin2_W;
		Zee_V=Znene*((value_type)2*sin_W*sin_W-(value_type)0.5);
		Zee_A=(value_type)0.5*Znene;
		Znmunmu=Znene;
		Zmumu_V=Zee_V;
		Zmumu_A=Zee_A;
		Zntauntau=Znene;
		Ztautau_V=Zee_V;
		Ztautau_A=Zee_A;
		Zuu_V=Znene*((value_type)0.5-(value_type)4*sin_W*sin_W/(value_type)3);
		Zuu_A=-Zee_A;
		Zdd_V=Znene*(-(value_type)0.5+(value_type)2*sin_W*sin_W/(value_type)3);
		Zdd_A=Zee_A;
		Zcc_V=Zuu_V;
		Zcc_A=Zuu_A;
		Zss_V=Zdd_V;
		Zss_A=Zdd_A;
		Ztt_V=Zuu_V;
		Ztt_A=Zuu_A;
		Zbb_V=Zdd_V;
		Zbb_A=Zdd_A;

		Wnee=gammaee/(std::sqrt((value_type)2)*sin_W);
		Wene=Wnee;
		Wnmumu=Wnee;
		Wmunmu=Wnmumu;
		Wntautau=Wnee;
		Wtauntau=Wntautau;

		Wud=Wnee*V_CKM[0][0];
		Wdu=Wnee*std::conj(V_CKM[0][0]);
		Wus=Wnee*V_CKM[0][1];
		Wsu=Wnee*std::conj(V_CKM[0][1]);
		Wub=Wnee*V_CKM[0][2];
		Wbu=Wnee*std::conj(V_CKM[0][2]);

		Wcd=Wnee*V_CKM[1][0];
		Wdc=Wnee*std::conj(V_CKM[1][0]);
		Wcs=Wnee*V_CKM[1][1];
		Wsc=Wnee*std::conj(V_CKM[1][1]);
		Wcb=Wnee*V_CKM[1][2];
		Wbc=Wnee*std::conj(V_CKM[1][2]);

		Wtd=Wnee*V_CKM[2][0];
		Wdt=Wnee*std::conj(V_CKM[2][0]);
		Wts=Wnee*V_CKM[2][1];
		Wst=Wnee*std::conj(V_CKM[2][1]);
		Wtb=Wnee*V_CKM[2][2];
		Wbt=Wnee*std::conj(V_CKM[2][2]);

		/* SFF couplings: */

		refresh_Yukawa_couplings();

		/* VSS couplings: */

		Zchih=Q_e/sin2_W;
		gammaphiphi=-gammaee;
		Zphiphi=gammaee*cos2_W/sin2_W;
		Wpphimh=-(value_type)0.5*gammaee/sin_W;
		Wmphiph=-Wpphimh;
		Wpphimchi=(value_type)0.5*Q_e/sin_W;
		Wmphipchi=Wpphimchi;

		/* SVV couplings: */

		hWW=gammaee*cM_W/sin_W;
		hZZ=(value_type)2*gammaee*cM_Z/sin2_W;
		phipWmZ=-gammaee*cM_W*tan_W;
		phimWpZ=phipWmZ;
		phipWmgamma=-gammaee*cM_W;
		phimWpgamma=phipWmgamma;

		/* SSVV couplings: */

		hhWW=std::complex<value_type>((value_type)0,(value_type)0.5)*std::pow(Q_e/sin_W,2);
		chichiWW=hhWW;
		phiphiWW=hhWW;
		hhZZ=hhWW/(cos_W*cos_W);
		chichiZZ=hhZZ;
		phiphiZZ=hhZZ*cos2_W*cos2_W;
		phiphigammaZ=-(value_type)2*hhWW*cos2_W*tan_W;
		phiphigammagamma=(value_type)2*Q_e*gammaee;
		phiphWmZ=-(value_type)0.5*Q_e*gammaee/cos_W;
		phimhWpZ=phiphWmZ;
		phiphWmgamma=phiphWmZ/tan_W;
		phimhWpgamma=phiphWmgamma;
		phimchiWpZ=std::complex<value_type>((value_type)0,(value_type)1)*phiphWmZ;
		phipchiWmZ=-phimchiWpZ;
		phimchiWpgamma=std::complex<value_type>((value_type)0,(value_type)1)*phiphWmgamma;
		phipchiWmgamma=-phimchiWpgamma;
	    }

	    /* Function computing the couplings from given parameters: */

	    static void refresh_Yukawa_couplings()
	    {
		std::complex<typename model_t::value_type>prefactor=-(value_type)0.5*gammaee/(sin_W*cM_W);

		hee=prefactor*M_e;
		hmumu=prefactor*M_mu;
		htautau=prefactor*M_tau;
		huu=prefactor*M_u;
		hdd=prefactor*M_d;
		hcc=prefactor*M_c;
		hss=prefactor*M_s;
		htt=prefactor*M_t;
		hbb=prefactor*M_b;

		prefactor=(value_type)0.5*Q_e/(sin_W*cM_W);

		chiee=prefactor*M_e;
		chimumu=prefactor*M_mu;
		chitautau=prefactor*M_tau;
		chiuu=-prefactor*M_u;
		chidd=prefactor*M_d;
		chicc=-prefactor*M_c;
		chiss=prefactor*M_s;
		chitt=-prefactor*M_t;
		chibb=prefactor*M_b;

		prefactor=-gammaee/(std::sqrt((value_type)2)*sin_W*cM_W);

		phinee=prefactor*M_e;
		phiene=phinee;
		phinmumu=prefactor*M_mu;
		phimunmu=phinmumu;
		phintautau=prefactor*M_tau;
		phitauntau=phintautau;

		prefactor=(value_type)0.5*gammaee/(std::sqrt((value_type)2)*sin_W*cM_W);

		phiud_S=prefactor*V_CKM[0][0]*(M_u-M_d);
		phiud_A=-prefactor*V_CKM[0][0]*(M_u+M_d);
		phidu_S=prefactor*std::conj(V_CKM[0][0])*(M_u-M_d);
		phidu_A=prefactor*std::conj(V_CKM[0][0])*(M_u+M_d);
		phius_S=prefactor*V_CKM[0][1]*(M_u-M_s);
		phius_A=-prefactor*V_CKM[0][1]*(M_u+M_s);
		phisu_S=prefactor*std::conj(V_CKM[0][1])*(M_u-M_s);
		phisu_A=prefactor*std::conj(V_CKM[0][1])*(M_u+M_s);
		phiub_S=prefactor*V_CKM[0][2]*(M_u-M_b);
		phiub_A=-prefactor*V_CKM[0][2]*(M_u+M_b);
		phibu_S=prefactor*std::conj(V_CKM[0][2])*(M_u-M_b);
		phibu_A=prefactor*std::conj(V_CKM[0][2])*(M_u+M_b);

		phicd_S=prefactor*V_CKM[1][0]*(M_c-M_d);
		phicd_A=-prefactor*V_CKM[1][0]*(M_c+M_d);
		phidc_S=prefactor*std::conj(V_CKM[1][0])*(M_c-M_d);
		phidc_A=prefactor*std::conj(V_CKM[1][0])*(M_c+M_d);
		phics_S=prefactor*V_CKM[1][1]*(M_c-M_s);
		phics_A=-prefactor*V_CKM[1][1]*(M_c+M_s);
		phisc_S=prefactor*std::conj(V_CKM[1][1])*(M_c-M_s);
		phisc_A=prefactor*std::conj(V_CKM[1][1])*(M_c+M_s);
		phicb_S=prefactor*V_CKM[1][2]*(M_c-M_b);
		phicb_A=-prefactor*V_CKM[1][2]*(M_c+M_b);
		phibc_S=prefactor*std::conj(V_CKM[1][2])*(M_c-M_b);
		phibc_A=prefactor*std::conj(V_CKM[1][2])*(M_c+M_b);

		phitd_S=prefactor*V_CKM[2][0]*(M_t-M_d);
		phitd_A=-prefactor*V_CKM[2][0]*(M_t+M_d);
		phidt_S=prefactor*std::conj(V_CKM[2][0])*(M_t-M_d);
		phidt_A=prefactor*std::conj(V_CKM[2][0])*(M_t+M_d);
		phits_S=prefactor*V_CKM[2][1]*(M_t-M_s);
		phits_A=-prefactor*V_CKM[2][1]*(M_t+M_s);
		phist_S=prefactor*std::conj(V_CKM[2][1])*(M_t-M_s);
		phist_A=prefactor*std::conj(V_CKM[2][1])*(M_t+M_s);
		phitb_S=prefactor*V_CKM[2][2]*(M_t-M_b);
		phitb_A=-prefactor*V_CKM[2][2]*(M_t+M_b);
		phibt_S=prefactor*std::conj(V_CKM[2][2])*(M_t-M_b);
		phibt_A=prefactor*std::conj(V_CKM[2][2])*(M_t+M_b);
	    }

	    /* Function computing the CKM matrix from given parameters: */

	    static void refresh_CKM_matrix()
	    {
		c12=std::sqrt(((value_type)1+s12)*((value_type)1-s12));
		c23=std::sqrt(((value_type)1+s23)*((value_type)1-s23));
		c13=std::sqrt(((value_type)1+s13)*((value_type)1-s13));
		std::complex<value_type>d13=std::polar(s13,delta);

		V_CKM[0][0]=c12*c13;
		V_CKM[0][1]=s12*c13;
		V_CKM[0][2]=std::conj(d13);
		V_CKM[1][0]=-s12*c23-c12*s23*d13;
		V_CKM[1][1]=c12*c23-s12*s23*d13;
		V_CKM[1][2]=s23*c13;
		V_CKM[2][0]=s12*s23-c12*c23*d13;
		V_CKM[2][1]=-c12*s23-s12*c23*d13;
		V_CKM[2][2]=c23*c13;
	    }

	    /* Function putting possible nonzero fermion masses in the model
	     * during runtime: */

	    static void refresh_fermion_masses()
	    {
		if(model<model_t>::initialised())
		{
		    (M_e>(value_type)0)?model<model_t>::set_mass("e-",&M_e):model<model_t>::set_massless("e-");
		    (M_mu>(value_type)0)?model<model_t>::set_mass("mu-",&M_mu):model<model_t>::set_massless("mu-");
		    (M_tau>(value_type)0)?model<model_t>::set_mass("tau-",&M_tau):model<model_t>::set_massless("tau-");
		    (M_u>(value_type)0)?model<model_t>::set_mass("u",&M_u):model<model_t>::set_massless("u");
		    (M_d>(value_type)0)?model<model_t>::set_mass("d",&M_d):model<model_t>::set_massless("d");
		    (M_c>(value_type)0)?model<model_t>::set_mass("c",&M_c):model<model_t>::set_massless("c");
		    (M_s>(value_type)0)?model<model_t>::set_mass("s",&M_s):model<model_t>::set_massless("s");
		}
	    }

	    /* Function computing the LO standard model W,Z,h and t widths: */

	    static void refresh_widths()
	    {
		refresh_W_width();
		refresh_Z_width();
		refresh_Higgs_width();
		refresh_top_width();
	    }

	    /* Function computing the LO standard model W width: */

	    static void refresh_W_width()
	    {
		W_W=value_type(3+2*model_t::N_c)*G_F*std::pow(M_W,3)/((value_type)6*pi*std::sqrt((value_type)2));
	    }

	    /* Function computing the LO standard model Z width: */

	    static void refresh_Z_width()
	    {
		value_type cw2=(M_W*M_W)/(M_Z*M_Z);
		value_type factor=G_F*std::pow(M_Z,3)/(108*std::sqrt((value_type)2)*pi);
		W_Z=(value_type(81+49*model_t::N_c)-(value_type)(162+92*model_t::N_c)*cw2+(value_type)(108+88*model_t::N_c)*cw2*cw2)*factor;
	    }

	    /* Function computing the LO standard model H width: */

	    static void refresh_Higgs_width()
	    {
		W_h0=0;
		value_type ml[3]={M_e,M_mu,M_tau};
		value_type prefactor=G_F*M_h0/(std::sqrt((value_type)2)*(value_type)4*pi);
		for(int i=0;i<3;++i)
		{
		    if(ml[i]>(value_type)0 and ml[i]<=(value_type)0.5*M_h0)
		    {
			W_h0+=prefactor*ml[i]*ml[i]*std::pow((value_type)1-std::pow(2*ml[i]/M_h0,(int)2),(value_type)1.5);
		    }
		}
		value_type mq[6]={M_d,M_u,M_s,M_c,M_b,M_t};
		prefactor*=(model_t::N_c);
		for(int i=0;i<6;++i)
		{
		    if(mq[i]>(value_type)0 and mq[i]<=(value_type)0.5*M_h0)
		    {
			W_h0+=prefactor*mq[i]*mq[i]*std::pow((value_type)1-std::pow(2*mq[i]/M_h0,(int)2),(value_type)1.5);
		    }
		}
		if(M_W<=(value_type)0.5*M_h0)
		{
		    value_type cw2=M_W*M_W/(M_Z*M_Z);
		    value_type xZ=(value_type)2*M_Z/M_h0;
		    xZ*=xZ;
		    value_type xW=(value_type)2*M_W/M_h0;
		    xW*=xW;

		    W_h0+=Q_e*Q_e*M_h0/((value_type)16*pi*((value_type)1-cw2)*xW)*std::sqrt((value_type)1-xW)*((value_type)1-xW+(value_type)0.75*xW*xW);
		}
		if(M_Z<=(value_type)0.5*M_h0)
		{
		    value_type cw2=M_W*M_W/(M_Z*M_Z);
		    value_type xZ=(value_type)2*M_Z/M_h0;
		    xZ*=xZ;
		    value_type xW=(value_type)2*M_W/M_h0;
		    xW*=xW;

		    W_h0+=Q_e*Q_e*M_h0/((value_type)32*pi*((value_type)1-cw2)*xW)*std::sqrt((value_type)1-xZ)*((value_type)1-xZ+(value_type)0.75*xZ*xZ);
		}
	    }

	    /* Function computing the LO standard model t width: */

	    static void refresh_top_width()
	    {
		W_t=0;
		if(M_W+M_b<=M_t)
		{
		    value_type xb=M_b/M_t;
		    xb*=xb;
		    value_type xW=M_W/M_t;
		    xW*=xW;
		    W_t+=G_F*M_t*M_t*M_t/(std::sqrt((value_type)2)*(value_type)8*pi)*std::sqrt((value_type)1+xb*xb+xW*xW-(value_type)2*(xb+xW+xb*xW))*((value_type)1-(value_type)2*xb+xb*xb+xW-(value_type)2*xW*xW+xb*xW);
		}
	    }
	    
	    /* Function setting the fine structure constant to a given value and
	     * computing all couplings depending on alpha: */
	    
	    static void set_alpha(const value_type& a)
	    {
		alpha=a;
		refresh_weak_couplings();
	    }
	    
	    /* Function setting alpha strong and computing the couplings
	     * depending on it: */
	    
	    static void set_alpha_s(const value_type& a)
	    {
		alpha_s=a;
		refresh_strong_couplings();
	    }
	    
	    /* Function setting the strong scale and computing alpha_s (at one
	     * loop) and the resulting strong-interaction vertex couplings: */
	    
	    static void set_QCD_scale(const value_type& mu)
	    {
		QCD_scale=mu;
		value_type s=mu/(value_type)SM_params::QCD_scale;
		alpha_s=(value_type)SM_params::alpha_s/((value_type)1+(value_type)SM_params::alpha_s*(value_type)(11*model_t::N_c-12)*std::log(s)/((value_type)6*pi));
		refresh_strong_couplings();
	    }

	    /* Function setting the Higgs mass and computing the couplings
	     * depending on it: */

	    static void set_Higgs_mass(const value_type& m)
	    {
		M_h0=m;
		refresh_Higgs_width();
		hhh=-(value_type)1.5*gammaee*M_h0*M_h0/(sin_W*cM_W);
		hchichi=hhh/(value_type)3;
		hphiphi=hchichi;
		hhhh=-(value_type)0.75*gammaee*Q_e*M_h0*M_h0/(sin_W*sin_W*cM_W*cM_W);
		hhchichi=hhhh/(value_type)3;
		hhphiphi=hhchichi;
		chichichichi=hhhh;
		chichiphiphi=hhphiphi;
		phiphiphiphi=(value_type)2*hhchichi;
	    }

	    /* Function setting all widths to zero: */

	    static void discard_widths()
	    {
		W_Z=(value_type)0;
		W_W=(value_type)0;
		W_h0=(value_type)0;
		W_t=(value_type)0;
		refresh_weak_couplings();
	    }

	    /* Function setting all vector propagators to the Feynman gauge: */

	    static void set_Feynman_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=0)
		{
		    if(gauge==1)
		    {
			model<model_t>::couple_particle("chi");
			model<model_t>::couple_particle("phi+");
			model<model_t>::decouple_zero_vertices();
		    }
		    model<model_t>::template set_propagator<Feynman_gauge>("W+");
		    model<model_t>::template set_propagator<Feynman_gauge>("W-");
		    model<model_t>::template set_propagator<Feynman_gauge>("Z");
		    model<model_t>::template set_propagator<Feynman_gauge>("gamma");
		    model<model_t>::template set_propagator<Feynman_gauge>("g");
		    model<model_t>::template set_propagator<scalar_propagator>("chi");
		    model<model_t>::template set_propagator<scalar_propagator>("phi+");
		    model<model_t>::template set_propagator<scalar_propagator>("phi-");
		}
		gauge=0;
	    }

	    /* Function setting all vector propagators to the unitary gauge and
	     * deleting all would-be Goldstones: */

	    static void set_unitary_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=1)
		{
		    model<model_t>::decouple_particle("phi+");
		    model<model_t>::decouple_particle("chi");
		    model<model_t>::template set_propagator<unitary_gauge>("W+");
		    model<model_t>::template set_propagator<unitary_gauge>("W-");
		    model<model_t>::template set_propagator<unitary_gauge>("Z");
		    model<model_t>::template set_propagator<unitary_gauge>("gamma");
		    model<model_t>::template set_propagator<unitary_gauge>("g");
		}
		gauge=1;
	    }

	    /* Function setting all vector propagators to the R-xi gauge: */

	    static void set_R_xi_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=2)
		{
		    if(gauge==1)
		    {
			model<model_t>::couple_particle("chi");
			model<model_t>::couple_particle("phi+");
			model<model_t>::decouple_zero_vertices();
		    }
		    model<model_t>::template set_propagator<R_vector_gauge>("W+");
		    model<model_t>::template set_propagator<R_vector_gauge>("W-");
		    model<model_t>::template set_propagator<R_vector_gauge>("Z");
		    model<model_t>::template set_propagator<R_vector_gauge>("gamma");
		    model<model_t>::template set_propagator<R_vector_gauge>("g");
		    model<model_t>::template set_propagator<R_scalar_gauge>("phi+");
		    model<model_t>::template set_propagator<R_scalar_gauge>("phi-");
		    model<model_t>::template set_propagator<R_scalar_gauge>("chi");
		}
		gauge=2;
	    }

	    /* Input of the xi-value in the R_xi gauge: */

	    static void set_xi(const value_type& x)
	    {
		R_gauge<model_t>::xi=x;
	    }

	    /* Function switching to a 4-gluon vertex Lagrangian: */

	    static void set_4_gluon_vertex()
	    {
		if(model<model_t>::initialised() and !four_gluon_vertex)
		{
		    model<model_t>::erase_vertex("H_qcd","g","g");
		    model<model_t>::template add_vertex<colour_tensor::ff_contr< SU<model_t::N_c> >,vvvv>("g","g","g","g",&gggg);
		}
		four_gluon_vertex=true;
	    }

	    /* Function switching to an auxiliary tensor field description of
	     * the 4-gluon vertex: */

	    static void set_auxiliary_QCD_field()
	    {
		if(model<model_t>::initialised() and four_gluon_vertex)
		{
		    model<model_t>::erase_vertex("g","g","g","g");
		    model<model_t>::template add_fast_4g_vertex< SU<model_t::N_c> >("g",&Tgg);
		}
		four_gluon_vertex=false;
	    }

	    /* Reset all physical parameters to their default values (defined in
	     * SM_params): */

	    static void set_default_params()
	    {
		M_e=(value_type)0;
		M_mu=(value_type)0;
		M_tau=(value_type)SM_params::M_tau;
		M_u=(value_type)0;
		M_d=(value_type)0;
		M_c=(value_type)SM_params::M_c;
		M_s=(value_type)0;
		M_t=(value_type)SM_params::M_t;
		M_b=(value_type)SM_params::M_b;
		M_Z=(value_type)SM_params::M_Z;
		refresh_fermion_masses();
		M_h0=(value_type)120;
		alpha=(value_type)SM_params::alpha;
		G_F=(value_type)SM_params::G_F;
		QCD_scale=(value_type)SM_params::QCD_scale;
		alpha_s=(value_type)SM_params::alpha_s;
		s12=(value_type)SM_params::s12;
		s23=(value_type)SM_params::s23;
		s13=(value_type)SM_params::s13;
		delta=(value_type)0;
		refresh_couplings();
		set_Feynman_gauge();
		set_auxiliary_QCD_field();
	    }

	private:
	    
	    /* QCD coupling: */

	    static value_type g_s;

	    /* Unit charge: */

	    static value_type Q_e;
	    
	    /* Complex gauge boson masses: */

	    static std::complex<value_type> cM_Z,cM_W;
	    
	    /* Complex cosine and sine of the Weinberg angle: */
	    
	    static std::complex<value_type> cos_W,sin_W;

	    /* Temporary internal parameters: */

	    static std::complex<value_type> tan_W,cos2_W,sin2_W;

	    /* CKM cosines: */

	    static value_type c12,c23,c13;

	    /* CKM matrix: */

	    static std::complex<value_type>V_CKM[3][3];

	    /* Gauge flag: */

	    static int gauge;
	    
	    /* Boolean denoting whether the 4-gluon vertex should be inserted,
	     * or the antisymmetric tensor field: */

	    static bool four_gluon_vertex;

    };
    
    template<class model_t,class value_t>const typename SM_base<model_t,value_t>::value_type SM_base<model_t,value_t>::pi=std::acos(-(value_t)1);
    
    /* Real input parameters: */

    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_e=(value_t)0;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_mu=(value_t)0;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_tau=(value_t)SM_params::M_tau;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_u=(value_t)0;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_d=(value_t)0;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_c=(value_t)SM_params::M_c;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_s=(value_t)0;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_t=(value_t)SM_params::M_t;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_b=(value_t)SM_params::M_b;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_Z=(value_t)SM_params::M_Z;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_h0=(value_t)120;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::alpha=(value_t)SM_params::alpha;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::G_F=(value_t)SM_params::G_F;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::QCD_scale=(value_t)SM_params::QCD_scale;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::alpha_s=(value_t)SM_params::alpha_s;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::s12=(value_t)SM_params::s12;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::s23=(value_t)SM_params::s23;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::s13=(value_t)SM_params::s13;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::delta=(value_t)0;

    /* Computed parameters: */
    
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::Q_e=std::sqrt((value_t)4*SM_base<model_t,value_t>::pi*SM_base<model_t,value_t>::alpha);
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::g_s=std::sqrt((value_t)4*SM_base<model_t,value_t>::pi*SM_base<model_t,value_t>::alpha_s);
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::M_W = SM_base<model_t,value_t>::M_Z*std::sqrt((value_t)0.5+(value_t)0.5*std::sqrt((value_t)1-(value_t)4*SM_base<model_t,value_t>::pi*SM_base<model_t,value_t>::alpha/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::G_F*SM_base<model_t,value_t>::M_Z*SM_base<model_t,value_t>::M_Z)));
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::W_W=(3+2*model_t::N_c)*SM_base<model_t,value_t>::G_F*std::pow(SM_base<model_t,value_t>::M_W,3)/(SM_base<model_t,value_t>::pi*std::sqrt((value_t)2))/(value_t)6;
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::W_Z = ((value_t)(162+49*model_t::N_c)-(value_t)(324+92*model_t::N_c)*std::pow(SM_base<model_t,value_t>::M_W/SM_base<model_t,value_t>::M_Z,2)+(value_t)(216+88*model_t::N_c)*std::pow(SM_base<model_t,value_t>::M_W/SM_base<model_t,value_t>::M_Z,4))*SM_base<model_t,value_t>::G_F*std::pow(SM_base<model_t,value_t>::M_Z,3)/((value_t)108*SM_base<model_t,value_t>::pi*std::sqrt((value_t)2));
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::W_h0 = (value_t)model_t::N_c*SM_base<model_t,value_t>::G_F*SM_base<model_t,value_t>::M_b*SM_base<model_t,value_t>::M_b*SM_base<model_t,value_t>::M_h0*std::pow((value_t)1-std::pow((value_t)2*SM_base<model_t,value_t>::M_b/SM_base<model_t,value_t>::M_h0,2),1.5)/(std::sqrt((value_t)2)*(value_t)4*SM_base<model_t,value_t>::pi);
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::W_t = SM_base<model_t,value_t>::G_F*std::pow(SM_base<model_t,value_t>::M_t,3)*((value_t)1-std::pow(SM_base<model_t,value_t>::M_W/SM_base<model_t,value_t>::M_t,2))*((value_t)1+std::pow(SM_base<model_t,value_t>::M_W/SM_base<model_t,value_t>::M_t,2)-2*std::pow(SM_base<model_t,value_t>::M_W/SM_base<model_t,value_t>::M_t,4))/(std::sqrt((value_t)128)*SM_base<model_t,value_t>::pi);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::cM_Z=std::sqrt(std::complex<value_t>(SM_base<model_t,value_t>::M_Z*SM_base<model_t,value_t>::M_Z,-SM_base<model_t,value_t>::M_Z*SM_base<model_t,value_t>::W_Z));
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::cM_W=std::sqrt(std::complex<value_t>(SM_base<model_t,value_t>::M_W*SM_base<model_t,value_t>::M_W,-SM_base<model_t,value_t>::M_W*SM_base<model_t,value_t>::W_W));
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::cos_W=SM_base<model_t,value_t>::cM_W/SM_base<model_t,value_t>::cM_Z;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::sin_W=std::sqrt(std::complex<value_t>(1,0)-SM_base<model_t,value_t>::cos_W*SM_base<model_t,value_t>::cos_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::tan_W=SM_base<model_t,value_t>::sin_W/SM_base<model_t,value_t>::cos_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::cos2_W=SM_base<model_t,value_t>::cos_W*SM_base<model_t,value_t>::cos_W-SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::sin_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::sin2_W=(value_t)2*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cos_W;
    
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::c12=std::sqrt(((value_t)1+SM_base<model_t,value_t>::s12)*((value_t)1-SM_base<model_t,value_t>::s12));
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::c23=std::sqrt(((value_t)1+SM_base<model_t,value_t>::s23)*((value_t)1-SM_base<model_t,value_t>::s23));
    template<class model_t,class value_t>value_t SM_base<model_t,value_t>::c13=std::sqrt(((value_t)1+SM_base<model_t,value_t>::s13)*((value_t)1-SM_base<model_t,value_t>::s13));
    
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::V_CKM[3][3]=
    {{ SM_base<model_t,value_t>::c12*SM_base<model_t,value_t>::c13,
       SM_base<model_t,value_t>::s12*SM_base<model_t,value_t>::c13,
       SM_base<model_t,value_t>::s13},
     {-SM_base<model_t,value_t>::s12*SM_base<model_t,value_t>::c23-SM_base<model_t,value_t>::c12*SM_base<model_t,value_t>::s23*SM_base<model_t,value_t>::s13,
       SM_base<model_t,value_t>::c12*SM_base<model_t,value_t>::c23-SM_base<model_t,value_t>::s12*SM_base<model_t,value_t>::s23*SM_base<model_t,value_t>::s13,
       SM_base<model_t,value_t>::s23*SM_base<model_t,value_t>::c13},
     { SM_base<model_t,value_t>::s12*SM_base<model_t,value_t>::s23-SM_base<model_t,value_t>::c12*SM_base<model_t,value_t>::c23*SM_base<model_t,value_t>::s13,
      -SM_base<model_t,value_t>::c12*SM_base<model_t,value_t>::s23-SM_base<model_t,value_t>::s12*SM_base<model_t,value_t>::c23*SM_base<model_t,value_t>::s13,
       SM_base<model_t,value_t>::c23*SM_base<model_t,value_t>::c13}};
    
    /* Gauge switch initialisation: */

    template<class model_t,class value_t>int SM_base<model_t,value_t>::gauge=0;

    /* Four-gluon vertex switch initialisation: */

    template<class model_t,class value_t>bool SM_base<model_t,value_t>::four_gluon_vertex=false;

    /* VVVV couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::WWWW=std::complex<value_t>(0,SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::Q_e)/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::sin_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::WZWZ=-SM_base<model_t,value_t>::WWWW*SM_base<model_t,value_t>::cos_W*SM_base<model_t,value_t>::cos_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::WgammaWZ=SM_base<model_t,value_t>::WWWW*SM_base<model_t,value_t>::cos_W*SM_base<model_t,value_t>::sin_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::WgammaWgamma=-SM_base<model_t,value_t>::WWWW*SM_base<model_t,value_t>::sin_W*sin_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gggg(0,-SM_base<model_t,value_t>::g_s*SM_base<model_t,value_t>::g_s);

    /* VVV couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammaWW(0,SM_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::ZWW=-SM_base<model_t,value_t>::gammaWW/SM_base<model_t,value_t>::tan_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::ggg(SM_base<model_t,value_t>::g_s,0);

    /* SSSS couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hhhh=-(value_t)0.75*std::complex<value_t>(0,SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::Q_e)*std::pow(SM_base<model_t,value_t>::M_h0/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W),2);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hhchichi=SM_base<model_t,value_t>::hhhh/(value_t)3;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hhphiphi=SM_base<model_t,value_t>::hhchichi;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chichichichi=SM_base<model_t,value_t>::hhhh;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chichiphiphi=SM_base<model_t,value_t>::hhphiphi;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiphiphiphi=(value_t)2*SM_base<model_t,value_t>::hhchichi;

    /* SSS couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hhh=-(value_t)1.5*std::complex<value_t>(0,SM_base<model_t,value_t>::Q_e)*SM_base<model_t,value_t>::M_h0*SM_base<model_t,value_t>::M_h0/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hchichi=SM_base<model_t,value_t>::hhh/(value_t)3;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hphiphi=SM_base<model_t,value_t>::hchichi;
    
    /* VFF couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammaee(0,SM_base<model_t,value_t>::Q_e);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammamumu=SM_base<model_t,value_t>::gammaee;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammatautau=SM_base<model_t,value_t>::gammaee;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammauu(0,-(value_t)2*SM_base<model_t,value_t>::Q_e/(value_t)3);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammadd(0,SM_base<model_t,value_t>::Q_e/(value_t)3);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammacc=SM_base<model_t,value_t>::gammauu;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammass=SM_base<model_t,value_t>::gammadd;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammatt=SM_base<model_t,value_t>::gammauu;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammabb=SM_base<model_t,value_t>::gammadd;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Znene=SM_base<model_t,value_t>::gammaee/SM_base<model_t,value_t>::sin2_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zee_V=SM_base<model_t,value_t>::Znene*((value_t)2*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::sin_W-(value_t)0.5);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zee_A=(value_t)0.5*SM_base<model_t,value_t>::Znene;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Znmunmu=SM_base<model_t,value_t>::Znene;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zmumu_V=SM_base<model_t,value_t>::Zee_V;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zmumu_A=SM_base<model_t,value_t>::Zee_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zntauntau=SM_base<model_t,value_t>::Znene;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Ztautau_V=SM_base<model_t,value_t>::Zee_V;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Ztautau_A=SM_base<model_t,value_t>::Zee_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zuu_V=SM_base<model_t,value_t>::Znene*((value_t)0.5-(value_t)4*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::sin_W/(value_t)3);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zuu_A=-SM_base<model_t,value_t>::Zee_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zdd_V=SM_base<model_t,value_t>::Znene*(-(value_t)0.5+(value_t)2*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::sin_W/(value_t)3);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zdd_A=SM_base<model_t,value_t>::Zee_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zcc_V=SM_base<model_t,value_t>::Zuu_V;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zcc_A=SM_base<model_t,value_t>::Zuu_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zss_V=SM_base<model_t,value_t>::Zdd_V;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zss_A=SM_base<model_t,value_t>::Zdd_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Ztt_V=SM_base<model_t,value_t>::Zuu_V;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Ztt_A=SM_base<model_t,value_t>::Zuu_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zbb_V=SM_base<model_t,value_t>::Zdd_V;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zbb_A=SM_base<model_t,value_t>::Zdd_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wnee=SM_base<model_t,value_t>::gammaee/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wnmumu=SM_base<model_t,value_t>::Wnee;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wntautau=SM_base<model_t,value_t>::Wnee;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wene=SM_base<model_t,value_t>::Wnee;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wmunmu=SM_base<model_t,value_t>::Wnmumu;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wtauntau=SM_base<model_t,value_t>::Wntautau;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wud=SM_base<model_t,value_t>::Wnee*SM_base<model_t,value_t>::V_CKM[0][0];
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wdu=SM_base<model_t,value_t>::Wud;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wus=SM_base<model_t,value_t>::Wnee*SM_base<model_t,value_t>::V_CKM[0][1];
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wsu=SM_base<model_t,value_t>::Wus;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wub=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wbu=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wcd=SM_base<model_t,value_t>::Wnee*SM_base<model_t,value_t>::V_CKM[1][0];
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wdc=SM_base<model_t,value_t>::Wcd;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wcs=SM_base<model_t,value_t>::Wnee*SM_base<model_t,value_t>::V_CKM[1][1];
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wsc=SM_base<model_t,value_t>::Wcs;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wcb=SM_base<model_t,value_t>::Wnee*SM_base<model_t,value_t>::V_CKM[1][2];
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wbc=SM_base<model_t,value_t>::Wcb;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wtd=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wdt=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wts=SM_base<model_t,value_t>::Wnee*SM_base<model_t,value_t>::V_CKM[2][1];
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wst=SM_base<model_t,value_t>::Wts;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wtb=SM_base<model_t,value_t>::Wnee*SM_base<model_t,value_t>::V_CKM[2][2];
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wbt=SM_base<model_t,value_t>::Wtb;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::guu(0,SM_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gdd(0,SM_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gcc(0,SM_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gss(0,SM_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gtt(0,SM_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gbb(0,SM_base<model_t,value_t>::g_s);

    /* SFF couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hee=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hmumu=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::htautau=-(value_t)0.5*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_tau/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::huu=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hdd=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hcc=-(value_t)0.5*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_c/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hss=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::htt=-(value_t)0.5*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_t/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hbb=-(value_t)0.5*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_b/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chiee=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chimumu=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chitautau=(value_t)0.5*SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::M_tau/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chiuu=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chidd=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chicc=-(value_t)0.5*SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::M_t/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chiss=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chitt=-(value_t)0.5*SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::M_t/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chibb=(value_t)0.5*SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::M_b/(SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phinee=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiene=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phinmumu=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phimunmu=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phintautau=-SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_tau/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phitauntau=SM_base<model_t,value_t>::phintautau;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiud_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiud_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phidu_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phidu_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phius_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phius_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phisu_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phisu_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiub_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiub_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phibu_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phibu_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phicd_S=(value_t)0.5*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_c*SM_base<model_t,value_t>::V_CKM[1][0]/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phicd_A=-SM_base<model_t,value_t>::phicd_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phidc_S=SM_base<model_t,value_t>::phicd_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phidc_A=-SM_base<model_t,value_t>::phicd_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phics_S=(value_t)0.5*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_c*SM_base<model_t,value_t>::V_CKM[1][1]/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phics_A=-SM_base<model_t,value_t>::phics_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phisc_S=SM_base<model_t,value_t>::phics_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phisc_A=-SM_base<model_t,value_t>::phics_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phicb_S=(value_t)0.5*SM_base<model_t,value_t>::gammaee*(SM_base<model_t,value_t>::M_c-SM_base<model_t,value_t>::M_b)*SM_base<model_t,value_t>::V_CKM[1][2]/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phicb_A=-(value_t)0.5*SM_base<model_t,value_t>::gammaee*(SM_base<model_t,value_t>::M_c+SM_base<model_t,value_t>::M_b)*SM_base<model_t,value_t>::V_CKM[1][2]/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phibc_S=SM_base<model_t,value_t>::phicb_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phibc_A=-SM_base<model_t,value_t>::phicb_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phitd_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phitd_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phidt_S=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phidt_A=(value_t)0;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phits_S=(value_t)0.5*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::M_t*SM_base<model_t,value_t>::V_CKM[2][1]/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phits_A=-SM_base<model_t,value_t>::phits_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phist_S=SM_base<model_t,value_t>::phits_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phist_A=-SM_base<model_t,value_t>::phits_A;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phitb_S=(value_t)0.5*SM_base<model_t,value_t>::gammaee*(SM_base<model_t,value_t>::M_t-SM_base<model_t,value_t>::M_b)/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phitb_A=-(value_t)0.5*SM_base<model_t,value_t>::gammaee*(SM_base<model_t,value_t>::M_t+SM_base<model_t,value_t>::M_b)/(std::sqrt((value_t)2)*SM_base<model_t,value_t>::sin_W*SM_base<model_t,value_t>::cM_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phibt_S=SM_base<model_t,value_t>::phitb_S;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phibt_A=-SM_base<model_t,value_t>::phitb_A;

    /* VSS couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zchih=SM_base<model_t,value_t>::Q_e/SM_base<model_t,value_t>::sin2_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::gammaphiphi=-SM_base<model_t,value_t>::gammaee;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Zphiphi=SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::cos2_W/SM_base<model_t,value_t>::sin2_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wpphimh=-(value_t)0.5*SM_base<model_t,value_t>::gammaee/SM_base<model_t,value_t>::sin_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wmphiph=-SM_base<model_t,value_t>::Wpphimh;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wpphimchi=(value_t)0.5*SM_base<model_t,value_t>::Q_e/SM_base<model_t,value_t>::sin_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Wmphipchi=SM_base<model_t,value_t>::Wpphimchi;
    
    /* SVV couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hWW=SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::cM_W/SM_base<model_t,value_t>::sin_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hZZ=(value_t)2*SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::cM_Z/SM_base<model_t,value_t>::sin2_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phipWmZ=-SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::cM_W*SM_base<model_t,value_t>::tan_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phimWpZ=SM_base<model_t,value_t>::phipWmZ;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phipWmgamma=-SM_base<model_t,value_t>::gammaee*SM_base<model_t,value_t>::cM_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phimWpgamma=SM_base<model_t,value_t>::phipWmgamma;
    
    /* SSVV couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hhWW=std::complex<value_t>((value_t)0,(value_t)0.5)*std::pow(SM_base<model_t,value_t>::Q_e/SM_base<model_t,value_t>::sin_W,2);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chichiWW=SM_base<model_t,value_t>::hhWW;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiphiWW=SM_base<model_t,value_t>::hhWW;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::hhZZ=SM_base<model_t,value_t>::hhWW/(SM_base<model_t,value_t>::cos_W*SM_base<model_t,value_t>::cos_W);
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::chichiZZ=SM_base<model_t,value_t>::hhZZ;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiphiZZ=SM_base<model_t,value_t>::hhZZ*SM_base<model_t,value_t>::cos2_W*SM_base<model_t,value_t>::cos2_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiphigammaZ=-(value_t)2*SM_base<model_t,value_t>::hhWW*SM_base<model_t,value_t>::cos2_W*SM_base<model_t,value_t>::tan_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiphigammagamma=(value_t)2*SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::gammaee;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiphWmZ=-(value_t)0.5*SM_base<model_t,value_t>::Q_e*SM_base<model_t,value_t>::gammaee/SM_base<model_t,value_t>::cos_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phimhWpZ=SM_base<model_t,value_t>::phiphWmZ;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phiphWmgamma=SM_base<model_t,value_t>::phimhWpZ/SM_base<model_t,value_t>::tan_W;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phimhWpgamma=SM_base<model_t,value_t>::phiphWmgamma;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phimchiWpZ=std::complex<value_t>((value_t)0,(value_t)1)*SM_base<model_t,value_t>::phiphWmZ;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phipchiWmZ=-SM_base<model_t,value_t>::phimchiWpZ;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phimchiWpgamma=std::complex<value_t>((value_t)0,(value_t)1)*SM_base<model_t,value_t>::phiphWmgamma;
    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::phipchiWmgamma=-SM_base<model_t,value_t>::phimchiWpgamma;

    /* TVV couplings: */

    template<class model_t,class value_t>std::complex<value_t>SM_base<model_t,value_t>::Tgg(0,std::sqrt((value_t)0.5)*SM_base<model_t,value_t>::g_s);
}

#endif /*CAMGEN_SM_H_*/


