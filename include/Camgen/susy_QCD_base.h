//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SUSY_QCD_BASE_H_
#define CAMGEN_SUSY_QCD_BASE_H_

#include <Camgen/model.h>
#include <Camgen/SM_params.h>
#include <Camgen/vector_particle.h>
#include <Camgen/scalar_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/su(n).h>
#include <Camgen/vvv.h>
#include <Camgen/vvvv.h>
#include <Camgen/vff.h>
#include <Camgen/sff.h>
#include <Camgen/sffL.h>
#include <Camgen/sffR.h>
#include <Camgen/ssvv.h>
#include <Camgen/ssss.h>
#include <Camgen/vss.h>
#include <Camgen/f.h>
#include <Camgen/ff_contr.h>
#include <Camgen/T.h>
#include <Camgen/TT_contr.h>
#include <Camgen/TT_plus.h>
#include <Camgen/dd_plus.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration of the supersymmetric QCD base class in Camgen. Deriving from    *
 * susy_QCD<Derived,value_type> inserts automatically all susy QCD particles and *
 * vertices in the Derived model. Note that for this to work, the numerical type *
 * in Derived should be castable to the value_type in the base class, i.e. the   *
 * second template argument.                                                     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */                                                                              

namespace Camgen
{
    template<class model_t,class value_t>class susy_QCD_base: public model<model_t>
    {
	public:
	    typedef value_t value_type;

	    /* Fine structure constant (not used): */

	    static const value_type alpha;

	    /*  Strong coupling constant: */

	    static value_type alpha_s;

	    /* QCD scale: */

	    static value_type QCD_scale;

	    /* QCD coupling: */

	    static value_type g_s;
	    
	    /* Pi: */
	    
	    static const value_type pi;
	    
	    /* Particle masses: */

	    static value_type M_u,M_d,M_c,M_s,M_t,M_b,M_suL,M_suR,M_sdL,M_sdR,M_scL,M_scR,M_ssL,M_ssR,M_stL,M_stR,M_sbL,M_sbR,M_sg;

	    /* Particle widths: */

	    static value_type W_t,W_suL,W_suR,W_sdL,W_sdR,W_scL,W_scR,W_ssL,W_ssR,W_stL,W_stR,W_sbL,W_sbR,W_sg;

	    /* Triple gluon coupling: */

	    static std::complex<value_type> ggg;
	    
	    /* 4-gluon coupling: */

	    static std::complex<value_type> gggg;
	    
	    /* Auxiliary tensor-gluon coupling: */

	    static std::complex<value_type> Tgg;

	    /* Gluon-quark-quark couplings: */

	    static std::complex<value_type> guu,gdd,gcc,gss,gtt,gbb;

	    /* Gluon-squark-squark couplings: */

	    static std::complex<value_type> gsuLsuL,gsuRsuR,gsdLsdL,gsdRsdR,gscLscL,gscRscR,gssLssL,gssRssR,gstLstL,gstRstR,gsbLsbL,gsbRsbR;

	    /* Squark-gluino-quark couplings: */
	    
	    static std::complex<value_type> suLsgu,suRsgu,sdLsgd,sdRsgd,scLsgc,scRsgc,ssLsgs,ssRsgs,stLsgt,stRsgt,sbLsgb,sbRsgb;

	    /* Antisquark-antiquark-gluino couplings: */

	    static std::complex<value_type> suLusg,suRusg,sdLdsg,sdRdsg,scLcsg,scRcsg,ssLssg,ssRssg,stLtsg,stRtsg,sbLbsg,sbRbsg;
	    
	    /* Squark-squark-gluon-gluon couplings: */
	    
	    static std::complex<value_type>suLsuLgg,suRsuRgg,sdLsdLgg,sdRsdRgg,scLscLgg,scRscRgg,ssLssLgg,ssRssRgg,stLstLgg,stRstRgg,sbLsbLgg,sbRsbRgg;

	    /* 4-squark couplings: */

	    static std::complex<value_type>suLsuLsuLsuL,suRsuRsuRsuR,sdLsdLsdLsdL,sdRsdRsdRsdR,scLscLscLscL,scRscRscRscR,ssLssLssLssL,ssRssRssRssR,stLstLstLstL,stRstRstRstR,sbLsbLsbLsbL,sbRsbRsbRsbR;
	    
	    /* 2-2-squark couplings: */
	    
	    static std::complex<value_type>suLsuLsuRsuR,sdLsdLsdRsdR,scLscLscRscR,ssLssLssRssR,stLstLstRstR,sbLsbLsbRsbR;

	    /* Constructor: */

	    susy_QCD_base()
	    {
		/* Gluon definition: */

		switch(gauge)
		{
		    case 0:
			model<model_t>::template add_gluon<Feynman_gauge,SU<model_t::N_c> >();
			break;
		    case 1:
			model<model_t>::template add_gluon<unitary_gauge,SU<model_t::N_c> >();
			break;
		    case 2:
			model<model_t>::template add_gluon<R_vector_gauge,SU<model_t::N_c> >();
			break;
		    default:
			model<model_t>::template add_gluon<Feynman_gauge,SU<model_t::N_c> >();
		}

		/* Gluino definition: */

		(M_sg>(value_type)0)?model<model_t>::template add_fermion< adjoint_rep< SU<model_t::N_c> > >("~g",&M_sg,&W_sg,1000021):model<model_t>::template add_fermion< adjoint_rep< SU<model_t::N_c> > >("~g",1000021);

		/* Quark definitions: */

		(M_d>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("d","dbar",&M_d,1):model<model_t>::template add_quarks< SU<model_t::N_c> >("d","dbar",1);
		(M_u>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("u","ubar",&M_u,2):model<model_t>::template add_quarks< SU<model_t::N_c> >("u","ubar",2);
		(M_s>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("s","sbar",&M_s,3):model<model_t>::template add_quarks< SU<model_t::N_c> >("s","sbar",3);
		(M_c>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("c","cbar",&M_c,4):model<model_t>::template add_quarks< SU<model_t::N_c> >("c","cbar",4);
		(M_b>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("b","bbar",&M_b,5):model<model_t>::template add_quarks< SU<model_t::N_c> >("b","bbar",5);
		(M_t>(value_type)0)?model<model_t>::template add_quarks< SU<model_t::N_c> >("t","tbar",&M_t,&W_t,6):model<model_t>::template add_quarks< SU<model_t::N_c> >("t","tbar",6);

		/* Squark definitions: */

		(M_sdL>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~d_L+","~d_L-",&M_sdL,&W_sdL,1000001):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~d_L+","~d_L-",1000001);
		(M_sdR>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~d_R+","~d_R-",&M_sdR,&W_sdR,2000001):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~d_R+","~d_R-",2000001);
		(M_suL>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~u_L+","~u_L-",&M_suL,&W_suL,1000002):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~u_L+","~u_L-",1000002);
		(M_suR>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~u_R+","~u_R-",&M_suR,&W_suR,2000002):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~u_R+","~u_R-",2000002);
		(M_ssL>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~s_L+","~s_L-",&M_ssL,&W_ssL,1000003):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~s_L+","~s_L-",1000003);
		(M_ssR>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~s_R+","~s_R-",&M_ssR,&W_ssR,2000003):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~s_R+","~s_R-",2000003);
		(M_scL>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~c_L+","~c_L-",&M_scL,&W_scL,1000004):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~c_L+","~c_L-",1000004);
		(M_scR>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~c_R+","~c_R-",&M_scR,&W_scR,2000004):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~c_R+","~c_R-",2000004);
		(M_sbL>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~b_L+","~b_L-",&M_sbL,&W_sbL,1000005):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~b_L+","~b_L-",1000005);
		(M_sbR>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~b_R+","~b_R-",&M_sbR,&W_sbR,2000005):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~b_R+","~b_R-",2000005);
		(M_sbL>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~t_L+","~t_L-",&M_stL,&W_stL,1000006):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~t_L+","~t_L-",1000006);
		(M_sbR>(value_type)0)?model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~t_R+","~t_R-",&M_stR,&W_stR,2000006):model<model_t>::template add_scalars< fundamental_rep< SU<model_t::N_c> > >("~t_R+","~t_R-",2000006);

		/* Triple gluon vertex definition: */

		model<model_t>::template add_vertex<colour_tensor::f< SU<model_t::N_c> >,vvv>("g","g","g",&ggg);

		/* Gluon-squark-squark vertices: */

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~d_L+","~d_L-",&gsdLsdL);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~d_R+","~d_R-",&gsdRsdR);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~u_L+","~u_L-",&gsuLsuL);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~u_R+","~u_R-",&gsuRsuR);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~s_L+","~s_L-",&gssLssL);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~s_R+","~s_R-",&gssRssR);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~c_L+","~c_L-",&gscLscL);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~c_R+","~c_R-",&gscRscR);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~b_L+","~b_L-",&gsbLsbL);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~b_R+","~b_R-",&gsbRsbR);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~t_L+","~t_L-",&gstLstL);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,0,2,1>,vss>("g","~t_R+","~t_R-",&gstRstR);

		/* Gluon-quark-quark vertices: */

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","dbar","d",&gdd);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","ubar","u",&guu);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","sbar","s",&gss);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","cbar","c",&gcc);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","bbar","b",&gbb);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> > >,vff>("g","tbar","t",&gtt);

		/* Gluon-gluino-gluino vertex: */

		model<model_t>::template add_vertex<colour_tensor::f<SU<model_t::N_c>,0,1,2>,vff>("g","~g","~g",&ggg);

		/* Squark-gluino-quark vertices: */

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffL>("~d_L+","~g","d",&sdLsgd);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffR>("~d_L-","dbar","~g",&sdLdsg);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffR>("~d_R+","~g","d",&sdRsgd);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffL>("~d_R-","dbar","~g",&sdRdsg);

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffL>("~u_L-","~g","u",&suLsgu);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffR>("~u_L+","ubar","~g",&suLusg);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffR>("~u_R-","~g","u",&suRsgu);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffL>("~u_R+","ubar","~g",&suRusg);

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffL>("~s_L+","~g","s",&ssLsgs);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffR>("~s_L-","sbar","~g",&ssLssg);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffR>("~s_R+","~g","s",&ssRsgs);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffL>("~s_R-","sbar","~g",&ssRssg);

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffL>("~c_L-","~g","c",&scLsgc);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffR>("~c_L+","cbar","~g",&scLcsg);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffR>("~c_R-","~g","c",&scRsgc);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffL>("~c_R+","cbar","~g",&scRcsg);

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffL>("~b_L+","~g","b",&sbLsgb);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffR>("~b_L-","bbar","~g",&sbLbsg);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffR>("~b_R+","~g","b",&sbRsgb);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffL>("~b_R-","bbar","~g",&sbRbsg);

		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffL>("~t_L-","~g","t",&stLsgt);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffR>("~t_L+","tbar","~g",&stLtsg);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,1,0,2>,sffR>("~t_R-","~g","t",&stRsgt);
		model<model_t>::template add_vertex<colour_tensor::T<fundamental_rep< SU<model_t::N_c> >,2,1,0>,sffL>("~t_R+","tbar","~g",&stRtsg);

		/* 4-gluon vertex definition: */

		if(four_gluon_vertex)
		{
		    model<model_t>::template add_vertex<colour_tensor::ff_contr< SU<model_t::N_c> >,vvvv>("g","g","g","g",&gggg);
		}
		else
		{
		    model<model_t>::template add_fast_4g_vertex< SU<model_t::N_c> >("g",&Tgg);
		}

		/* Squark-squark-gluon-gluon vertices: */

		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~d_L-","~d_L+","g","g",&sdLsdLgg);
		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~d_R-","~d_R+","g","g",&sdRsdRgg);

		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~u_L-","~u_L+","g","g",&suLsuLgg);
		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~u_R-","~u_R+","g","g",&suRsuRgg);

		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~s_L-","~s_L+","g","g",&ssLssLgg);
		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~s_R-","~s_R+","g","g",&ssRssRgg);

		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~c_L-","~c_L+","g","g",&scLscLgg);
		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~c_R-","~c_R+","g","g",&scRscRgg);

		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~b_L-","~b_L+","g","g",&sbLsbLgg);
		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~b_R-","~b_R+","g","g",&sbRsbRgg);

		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~t_L-","~t_L+","g","g",&stLstLgg);
		model<model_t>::template add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<model_t::N_c> >,2,3,0,1>,ssvv>("~t_R-","~t_R+","g","g",&stRstRgg);

		/* 4-Squark vertices: */

		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~d_L+","~d_L-","~d_L+","~d_L-",&sdLsdLsdLsdL);
		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~d_R+","~d_R-","~d_R+","~d_R-",&sdRsdRsdRsdR);
		model<model_t>::template add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<model_t::N_c> > >,ssss>("~d_L+","~d_L-","~d_R+","~d_R-",&sdLsdLsdRsdR);

		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~u_L+","~u_L-","~u_L+","~u_L-",&suLsuLsuLsuL);
		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~u_R+","~u_R-","~u_R+","~u_R-",&suRsuRsuRsuR);
		model<model_t>::template add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<model_t::N_c> > >,ssss>("~u_L+","~u_L-","~u_R+","~u_R-",&suLsuLsuRsuR);

		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~s_L+","~s_L-","~s_L+","~s_L-",&ssLssLssLssL);
		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~s_R+","~s_R-","~s_R+","~s_R-",&ssRssRssRssR);
		model<model_t>::template add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<model_t::N_c> > >,ssss>("~s_L+","~s_L-","~s_R+","~s_R-",&ssLssLssRssR);

		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~c_L+","~c_L-","~c_L+","~c_L-",&scLscLscLscL);
		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~c_R+","~c_R-","~c_R+","~c_R-",&scRscRscRscR);
		model<model_t>::template add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<model_t::N_c> > >,ssss>("~c_L+","~c_L-","~c_R+","~c_R-",&scLscLscRscR);

		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~b_L+","~b_L-","~b_L+","~b_L-",&sbLsbLsbLsbL);
		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~b_R+","~b_R-","~b_R+","~b_R-",&sbRsbRsbRsbR);
		model<model_t>::template add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<model_t::N_c> > >,ssss>("~b_L+","~b_L-","~b_R+","~b_R-",&sbLsbLsbRsbR);

		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~t_L+","~t_L-","~t_L+","~t_L-",&stLstLstLstL);
		model<model_t>::template add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<model_t::N_c> > >,ssss>("~t_R+","~t_R-","~t_R+","~t_R-",&stRstRstRstR);
		model<model_t>::template add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<model_t::N_c> > >,ssss>("~t_L+","~t_L-","~t_R+","~t_R-",&stLstLstRstR);

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
		
		/* Left-handed positive squark family definitions: */

		model<model_t>::construct_family("~Q_L+","~d_L+,~u_L+,~s_L+,~c_L+,~b_L+,~t_L+");

		/* Left-handed negative squark family definitions: */

		model<model_t>::construct_family("~Q_L-","~d_L-,~u_L-,~s_L-,~c_L-,~b_L-,~t_L-");

		/* Right-handed positive squark family definitions: */

		model<model_t>::construct_family("~Q_R+","~d_R+,~u_R+,~s_R+,~c_R+,~b_R+,~t_R+");

		/* Right-handed negative squark family definitions: */

		model<model_t>::construct_family("~Q_R-","~d_R-,~u_R-,~s_R-,~c_R-,~b_R-,~t_R-");

		/* Left-handed squark family definitions: */

		model<model_t>::construct_family("~Q_L","~d_L+,~d_L-,~u_L+,~u_L-,~s_L+,~s_L-,~c_L+,~c_L-,~b_L+,~b_L-,~t_L+,~t_L-");

		/* Right-handed squark family definitions: */

		model<model_t>::construct_family("~Q_R","~d_R+,~d_R-,~u_R+,~u_R-,~s_R+,~s_R-,~c_R+,~c_R-,~b_R+,~b_R-,~t_R+,~t_R-");

		/* Positive squark family definitions: */

		model<model_t>::construct_family("~Q+","~d_L+,~d_R+,~u_L+,~u_R+,~s_L+,~s_R+,~c_L+,~c_R+,~b_L+,~b_R+,~t_L+,~t_R+");

		/* Negative squark family definitions: */

		model<model_t>::construct_family("~Q-","~d_L-,~d_R-,~u_L-,~u_R-,~s_L-,~s_R-,~c_L-,~c_R-,~b_L-,~b_R-,~t_L-,~t_R-");
	    }
	    
	    /* Computation of the couplings: */

	    static void refresh_couplings()
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

		gsuLsuL=std::complex<value_type>(0,g_s);
		gsdLsdL=gsuLsuL;
		gscLscL=gsuLsuL;
		gssLssL=gsuLsuL;
		gstLstL=gsuLsuL;
		gsbLsbL=gsuLsuL;

		gsuRsuR=gsuLsuL;
		gsdRsdR=gsdLsdL;
		gscRscR=gscLscL;
		gssRssR=gssLssL;
		gstRstR=gstLstL;
		gsbRsbR=gsbLsbL;

		suLsgu=std::complex<value_type>(0,-std::sqrt((value_type)2)*g_s);
		sdLsgd=suLsgu;
		scLsgc=suLsgu;
		ssLsgs=suLsgu;
		stLsgt=suLsgu;
		sbLsgb=suLsgu;

		suRsgu=suLsgu;
		sdRsgd=sdLsgd;
		scRsgc=scLsgc;
		ssRsgs=ssLsgs;
		stRsgt=stLsgt;
		sbRsgb=sbLsgb;

		suLusg=std::complex<value_type>(0,std::sqrt((value_type)2)*g_s);
		sdLdsg=suLusg;
		scLcsg=suLusg;
		ssLssg=suLusg;
		stLtsg=suLusg;
		sbLbsg=suLusg;

		suRusg=suLusg;
		sdRdsg=sdLdsg;
		scRcsg=scLcsg;
		ssRssg=ssLssg;
		stRtsg=stLtsg;
		sbRbsg=sbLbsg;

		suLsuLgg=std::complex<value_type>(0,g_s*g_s);
		sdLsdLgg=suLsuLgg;
		scLscLgg=suLsuLgg;
		ssLssLgg=suLsuLgg;
		stLstLgg=suLsuLgg;
		sbLsbLgg=suLsuLgg;

		suRsuRgg=suLsuLgg;
		sdRsdRgg=sdLsdLgg;
		scRscRgg=scLscLgg;
		ssRssRgg=ssLssLgg;
		stRstRgg=stLstLgg;
		sbRsbRgg=sbLsbLgg;

		suLsuLsuLsuL=std::complex<value_type>(0,-(value_type)(model_t::N_c-1)*g_s*g_s/value_type(2*model_t::N_c));
		sdLsdLsdLsdL=suLsuLsuLsuL;
		scLscLscLscL=suLsuLsuLsuL;
		ssLssLssLssL=suLsuLsuLsuL;
		stLstLstLstL=suLsuLsuLsuL;
		sbLsbLsbLsbL=suLsuLsuLsuL;

		suRsuRsuRsuR=suLsuLsuLsuL;
		sdRsdRsdRsdR=sdLsdLsdLsdL;
		scRscRscRscR=scLscLscLscL;
		ssRssRssRssR=ssLssLssLssL;
		stRstRstRstR=stLstLstLstL;
		sbRsbRsbRsbR=sbLsbLsbLsbL;

		suLsuLsuRsuR=std::complex<value_type>(0,g_s*g_s);
		sdLsdLsdRsdR=suLsuLsuRsuR;
		scLscLscRscR=suLsuLsuRsuR;
		ssLssLssRssR=suLsuLsuRsuR;
		stLstLstRstR=suLsuLsuRsuR;
		sbLsbLsbRsbR=suLsuLsuRsuR;
	    }
	    
	    /* Function refreshing the particle masses: */

	    static void refresh_masses()
	    {
		if(model<model_t>::initialised())
		{
		    (M_d>(value_type)0)?model<model_t>::set_mass("d",&M_d):model<model_t>::set_massless("d");
		    (M_u>(value_type)0)?model<model_t>::set_mass("u",&M_u):model<model_t>::set_massless("u");
		    (M_s>(value_type)0)?model<model_t>::set_mass("s",&M_s):model<model_t>::set_massless("s");
		    (M_c>(value_type)0)?model<model_t>::set_mass("c",&M_c):model<model_t>::set_massless("c");
		    (M_b>(value_type)0)?model<model_t>::set_mass("b",&M_b):model<model_t>::set_massless("b");
		    (M_t>(value_type)0)?model<model_t>::set_mass("t",&M_t):model<model_t>::set_massless("t");

		    (M_sdL>(value_type)0)?model<model_t>::set_mass("~d_L",&M_sdL):model<model_t>::set_massless("~d_L");
		    (M_suL>(value_type)0)?model<model_t>::set_mass("~u_L",&M_suL):model<model_t>::set_massless("~u_L");
		    (M_ssL>(value_type)0)?model<model_t>::set_mass("~s_L",&M_ssL):model<model_t>::set_massless("~s_L");
		    (M_scL>(value_type)0)?model<model_t>::set_mass("~c_L",&M_scL):model<model_t>::set_massless("~c_L");
		    (M_sbL>(value_type)0)?model<model_t>::set_mass("~b_L",&M_sbL):model<model_t>::set_massless("~b_L");
		    (M_stL>(value_type)0)?model<model_t>::set_mass("~t_L",&M_stL):model<model_t>::set_massless("~t_L");

		    (M_sdR>(value_type)0)?model<model_t>::set_mass("~d_R",&M_sdR):model<model_t>::set_massless("~d_R");
		    (M_suR>(value_type)0)?model<model_t>::set_mass("~u_R",&M_suR):model<model_t>::set_massless("~u_R");
		    (M_ssR>(value_type)0)?model<model_t>::set_mass("~s_R",&M_ssR):model<model_t>::set_massless("~s_R");
		    (M_scR>(value_type)0)?model<model_t>::set_mass("~c_R",&M_scR):model<model_t>::set_massless("~c_R");
		    (M_sbR>(value_type)0)?model<model_t>::set_mass("~b_R",&M_sbR):model<model_t>::set_massless("~b_R");
		    (M_stR>(value_type)0)?model<model_t>::set_mass("~t_R",&M_stR):model<model_t>::set_massless("~t_R");

		    (M_sg>(value_type)0)?model<model_t>::set_mass("~g",&M_sg):model<model_t>::set_massless("~g");
		}
	    }

	    /* Function refreshing the particle decay widths: */

	    static void refresh_widths()
	    {
		if(model<model_t>::initialised())
		{
		    (W_t>(value_type)0)?model<model_t>::set_width("t",&W_t):model<model_t>::set_widthless("t");
		    (W_sdL>(value_type)0)?model<model_t>::set_width("~d_L",&W_sdL):model<model_t>::set_widthless("~d_L");
		    (W_suL>(value_type)0)?model<model_t>::set_width("~u_L",&W_suL):model<model_t>::set_widthless("~u_L");
		    (W_ssL>(value_type)0)?model<model_t>::set_width("~s_L",&W_ssL):model<model_t>::set_widthless("~s_L");
		    (W_scL>(value_type)0)?model<model_t>::set_width("~c_L",&W_scL):model<model_t>::set_widthless("~c_L");
		    (W_sbL>(value_type)0)?model<model_t>::set_width("~b_L",&W_sbL):model<model_t>::set_widthless("~b_L");
		    (W_stL>(value_type)0)?model<model_t>::set_width("~t_L",&W_stL):model<model_t>::set_widthless("~t_L");
		    (W_sdR>(value_type)0)?model<model_t>::set_width("~d_R",&W_sdR):model<model_t>::set_widthless("~d_R");
		    (W_suR>(value_type)0)?model<model_t>::set_width("~u_R",&W_suR):model<model_t>::set_widthless("~u_R");
		    (W_ssR>(value_type)0)?model<model_t>::set_width("~s_R",&W_ssR):model<model_t>::set_widthless("~s_R");
		    (W_scR>(value_type)0)?model<model_t>::set_width("~c_R",&W_scR):model<model_t>::set_widthless("~c_R");
		    (W_sbR>(value_type)0)?model<model_t>::set_width("~b_R",&W_sbR):model<model_t>::set_widthless("~b_R");
		    (W_stR>(value_type)0)?model<model_t>::set_width("~t_R",&W_stR):model<model_t>::set_widthless("~t_R");
		    (W_sg>(value_type)0)?model<model_t>::set_width("~g",&W_sg):model<model_t>::set_widthless("~g");
		}
	    }

	    /* Function returning the number of QCD colours: */
	    
	    static std::size_t QCD_colours()
	    {
		return model_t::N_c;
	    }
	    
	    /* Input of fine structure constant (not used): */
	    
	    static void set_alpha(const value_type& a){}

	    /* Input of alpha-strong and computation of the couplings: */

	    static void set_alpha_s(const value_type& a)
	    {
		alpha_s=a;
		refresh_couplings();
	    }
	    
	    /* Input of the QCD scale and (LO) computation of alpha-strong and
	     * the couplings: */
	    
	    static void set_QCD_scale(const value_type& mu)
	    {
		QCD_scale=mu;
		value_type rho=mu/(value_type)SM_params::QCD_scale;
		alpha_s=SM_params::alpha_s/((value_type)1+(value_type)SM_params::alpha_s*value_type(11*model_t::N_c-12)*std::log(rho)/((value_type)6*pi));
		refresh_couplings();
	    }
	    
	    /* Function setting the gluon propagator to the Feynman gauge: */
	    
	    static void set_Feynman_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=0)
		{
		    model<model_t>::template set_propagator<Feynman_gauge>("g");
		}
		gauge=0;
	    }

	    /* Function setting the gluon propagator to the unitary gauge: */

	    static void set_unitary_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=1)
		{
		    model<model_t>::template set_propagator<unitary_gauge>("g");
		}
		gauge=1;
	    }
	    
	    /* Function setting the gluon propagator to the R-gauge: */
	    
	    static void set_R_xi_gauge()
	    {
		if(model<model_t>::initialised() and gauge!=2)
		{
		    model<model_t>::template set_propagator<R_vector_gauge>("g");
		}
		gauge=2;
	    }

	    /* Input of the gauge parameter: */

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
	    
	private:

	    /* Gauge flag: */

	    static int gauge;

	    /* 4-gluon vertex flag (true: model constructs a 4-gluon vertex,
	     * false: model uses the antisymmetric tensor field): */

	    static bool four_gluon_vertex;
    };
    template<class model_t,class value_t>const value_t susy_QCD_base<model_t,value_t>::pi=std::acos(-(value_t)1);
    template<class model_t,class value_t>const value_t susy_QCD_base<model_t,value_t>::alpha=-(value_t)1;
    
    /* Initialisation of real input parameters: */

    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::QCD_scale=(value_t)SM_params::QCD_scale;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::alpha_s=(value_t)SM_params::alpha_s;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::g_s=std::sqrt((value_t)4*susy_QCD_base<model_t,value_t>::pi*susy_QCD_base<model_t,value_t>::alpha_s);

    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_u=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_d=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_c=(value_t)SM_params::M_c;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_s=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_t=(value_t)SM_params::M_t;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_b=(value_t)SM_params::M_b;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_suL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_sdL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_scL=(value_t)SM_params::M_c;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_ssL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_stL=(value_t)SM_params::M_t;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_sbL=(value_t)SM_params::M_b;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_suR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_sdR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_scR=(value_t)SM_params::M_c;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_ssR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_stR=(value_t)SM_params::M_t;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_sbR=(value_t)SM_params::M_b;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::M_sg=(value_t)0;

    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_t=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_suL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_sdL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_scL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_ssL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_stL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_sbL=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_suR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_sdR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_scR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_ssR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_stR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_sbR=(value_t)0;
    template<class model_t,class value_t>value_t susy_QCD_base<model_t,value_t>::W_sg=(value_t)0;
    
    /* Initialisation of coupling constants: */

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ggg(susy_QCD_base<model_t,value_t>::g_s,0);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gggg(0,-susy_QCD_base<model_t,value_t>::g_s*susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::Tgg(0,std::sqrt((value_t)0.5)*susy_QCD_base<model_t,value_t>::g_s);
    
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::guu(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gdd(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gcc(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gss(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gtt(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gbb(0,susy_QCD_base<model_t,value_t>::g_s);
    
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gsuLsuL(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gsdLsdL(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gscLscL(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gssLssL(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gstLstL(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gsbLsbL(0,susy_QCD_base<model_t,value_t>::g_s);
    
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gsuRsuR(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gsdRsdR(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gscRscR(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gssRssR(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gstRstR(0,susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::gsbRsbR(0,susy_QCD_base<model_t,value_t>::g_s);

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suLsgu(0,-std::sqrt((value_t)2)*susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdLsgd=susy_QCD_base<model_t,value_t>::suLsgu;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scLsgc=susy_QCD_base<model_t,value_t>::suLsgu;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssLsgs=susy_QCD_base<model_t,value_t>::suLsgu;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stLsgt=susy_QCD_base<model_t,value_t>::suLsgu;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbLsgb=susy_QCD_base<model_t,value_t>::suLsgu;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suRsgu=susy_QCD_base<model_t,value_t>::suLsgu;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdRsgd=susy_QCD_base<model_t,value_t>::sdLsgd;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scRsgc=susy_QCD_base<model_t,value_t>::scLsgc;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssRsgs=susy_QCD_base<model_t,value_t>::ssLsgs;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stRsgt=susy_QCD_base<model_t,value_t>::stLsgt;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbRsgb=susy_QCD_base<model_t,value_t>::sbLsgb;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suLusg(0,std::sqrt((value_t)2)*susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdLdsg=susy_QCD_base<model_t,value_t>::suLusg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scLcsg=susy_QCD_base<model_t,value_t>::suLusg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssLssg=susy_QCD_base<model_t,value_t>::suLusg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stLtsg=susy_QCD_base<model_t,value_t>::suLusg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbLbsg=susy_QCD_base<model_t,value_t>::suLusg;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suRusg=susy_QCD_base<model_t,value_t>::suLusg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdRdsg=susy_QCD_base<model_t,value_t>::sdLdsg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scRcsg=susy_QCD_base<model_t,value_t>::scLcsg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssRssg=susy_QCD_base<model_t,value_t>::ssLssg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stRtsg=susy_QCD_base<model_t,value_t>::stLtsg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbRbsg=susy_QCD_base<model_t,value_t>::sbLbsg;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suLsuLgg(0,susy_QCD_base<model_t,value_t>::g_s*susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdLsdLgg=susy_QCD_base<model_t,value_t>::suLsuLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scLscLgg=susy_QCD_base<model_t,value_t>::suLsuLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssLssLgg=susy_QCD_base<model_t,value_t>::suLsuLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stLstLgg=susy_QCD_base<model_t,value_t>::suLsuLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbLsbLgg=susy_QCD_base<model_t,value_t>::suLsuLgg;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suRsuRgg=susy_QCD_base<model_t,value_t>::suLsuLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdRsdRgg=susy_QCD_base<model_t,value_t>::sdLsdLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scRscRgg=susy_QCD_base<model_t,value_t>::scLscLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssRssRgg=susy_QCD_base<model_t,value_t>::ssLssLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stRstRgg=susy_QCD_base<model_t,value_t>::stLstLgg;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbRsbRgg=susy_QCD_base<model_t,value_t>::sbLsbLgg;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suLsuLsuLsuL(0,-(value_t)(model_t::N_c-1)*susy_QCD_base<model_t,value_t>::g_s*susy_QCD_base<model_t,value_t>::g_s/value_t(2*model_t::N_c));
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdLsdLsdLsdL=susy_QCD_base<model_t,value_t>::suLsuLsuLsuL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scLscLscLscL=susy_QCD_base<model_t,value_t>::suLsuLsuLsuL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssLssLssLssL=susy_QCD_base<model_t,value_t>::suLsuLsuLsuL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stLstLstLstL=susy_QCD_base<model_t,value_t>::suLsuLsuLsuL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbLsbLsbLsbL=susy_QCD_base<model_t,value_t>::suLsuLsuLsuL;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suRsuRsuRsuR=susy_QCD_base<model_t,value_t>::suLsuLsuLsuL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdRsdRsdRsdR=susy_QCD_base<model_t,value_t>::sdLsdLsdLsdL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scRscRscRscR=susy_QCD_base<model_t,value_t>::scLscLscLscL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssRssRssRssR=susy_QCD_base<model_t,value_t>::ssLssLssLssL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stRstRstRstR=susy_QCD_base<model_t,value_t>::stLstLstLstL;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbRsbRsbRsbR=susy_QCD_base<model_t,value_t>::sbLsbLsbLsbL;

    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::suLsuLsuRsuR(0,susy_QCD_base<model_t,value_t>::g_s*susy_QCD_base<model_t,value_t>::g_s);
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sdLsdLsdRsdR=susy_QCD_base<model_t,value_t>::suLsuLsuRsuR;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::scLscLscRscR=susy_QCD_base<model_t,value_t>::suLsuLsuRsuR;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::ssLssLssRssR=susy_QCD_base<model_t,value_t>::suLsuLsuRsuR;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::stLstLstRstR=susy_QCD_base<model_t,value_t>::suLsuLsuRsuR;
    template<class model_t,class value_t>std::complex<value_t>susy_QCD_base<model_t,value_t>::sbLsbLsbRsbR=susy_QCD_base<model_t,value_t>::suLsuLsuRsuR;

    /* Gauge parameter initialisation: */

    template<class model_t,class value_t>int susy_QCD_base<model_t,value_t>::gauge=0;

    /* Four-gluon vertex initialisation: */

    template<class model_t,class value_t>bool susy_QCD_base<model_t,value_t>::four_gluon_vertex=false;
}

#endif /*CAMGEN_SUSY_QCD_BASE_H_*/

