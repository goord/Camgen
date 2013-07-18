//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/susy_QCD.h>
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

namespace Camgen
{
    /* Compile-time constants: */

    const std::size_t susy_QCD::dimension;
    const bool susy_QCD::coloured;
    const std::size_t susy_QCD::N_c;
    
    const bool susy_QCD::continuous_colours;
    const bool susy_QCD::continuous_helicities;
    
    const int susy_QCD::beam_direction;
    
    const susy_QCD::value_type susy_QCD::pi=std::acos(-(susy_QCD::value_type)1);
    const susy_QCD::value_type susy_QCD::alpha=-(susy_QCD::value_type)1;
    
    /* Initialisation of real input parameters: */

    susy_QCD::value_type susy_QCD::QCD_scale=(susy_QCD::value_type)SM_params::QCD_scale;
    susy_QCD::value_type susy_QCD::alpha_s=(susy_QCD::value_type)SM_params::alpha_s;
    susy_QCD::value_type susy_QCD::g_s=std::sqrt((susy_QCD::value_type)4*susy_QCD::pi*susy_QCD::alpha_s);

    susy_QCD::value_type susy_QCD::M_u=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_d=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_c=(susy_QCD::value_type)SM_params::M_c;
    susy_QCD::value_type susy_QCD::M_s=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_t=(susy_QCD::value_type)SM_params::M_t;
    susy_QCD::value_type susy_QCD::M_b=(susy_QCD::value_type)SM_params::M_b;
    susy_QCD::value_type susy_QCD::M_suL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_sdL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_scL=(susy_QCD::value_type)SM_params::M_c;
    susy_QCD::value_type susy_QCD::M_ssL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_stL=(susy_QCD::value_type)SM_params::M_t;
    susy_QCD::value_type susy_QCD::M_sbL=(susy_QCD::value_type)SM_params::M_b;
    susy_QCD::value_type susy_QCD::M_suR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_sdR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_scR=(susy_QCD::value_type)SM_params::M_c;
    susy_QCD::value_type susy_QCD::M_ssR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::M_stR=(susy_QCD::value_type)SM_params::M_t;
    susy_QCD::value_type susy_QCD::M_sbR=(susy_QCD::value_type)SM_params::M_b;
    susy_QCD::value_type susy_QCD::M_sg=(susy_QCD::value_type)0;

    susy_QCD::value_type susy_QCD::W_t=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_suL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_sdL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_scL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_ssL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_stL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_sbL=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_suR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_sdR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_scR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_ssR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_stR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_sbR=(susy_QCD::value_type)0;
    susy_QCD::value_type susy_QCD::W_sg=(susy_QCD::value_type)0;
    
    /* Initialisation of coupling constants: */

    std::complex<susy_QCD::value_type>susy_QCD::ggg(susy_QCD::g_s,0);
    std::complex<susy_QCD::value_type>susy_QCD::gggg(0,-susy_QCD::g_s*susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::Tgg(0,std::sqrt((susy_QCD::value_type)0.5)*susy_QCD::g_s);
    
    std::complex<susy_QCD::value_type>susy_QCD::guu(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gdd(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gcc(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gss(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gtt(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gbb(0,susy_QCD::g_s);
    
    std::complex<susy_QCD::value_type>susy_QCD::gsuLsuL(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gsdLsdL(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gscLscL(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gssLssL(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gstLstL(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gsbLsbL(0,susy_QCD::g_s);
    
    std::complex<susy_QCD::value_type>susy_QCD::gsuRsuR(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gsdRsdR(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gscRscR(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gssRssR(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gstRstR(0,susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::gsbRsbR(0,susy_QCD::g_s);

    std::complex<susy_QCD::value_type>susy_QCD::suLsgu(0,-std::sqrt((susy_QCD::value_type)2)*susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::sdLsgd=susy_QCD::suLsgu;
    std::complex<susy_QCD::value_type>susy_QCD::scLsgc=susy_QCD::suLsgu;
    std::complex<susy_QCD::value_type>susy_QCD::ssLsgs=susy_QCD::suLsgu;
    std::complex<susy_QCD::value_type>susy_QCD::stLsgt=susy_QCD::suLsgu;
    std::complex<susy_QCD::value_type>susy_QCD::sbLsgb=susy_QCD::suLsgu;

    std::complex<susy_QCD::value_type>susy_QCD::suRsgu=susy_QCD::suLsgu;
    std::complex<susy_QCD::value_type>susy_QCD::sdRsgd=susy_QCD::sdLsgd;
    std::complex<susy_QCD::value_type>susy_QCD::scRsgc=susy_QCD::scLsgc;
    std::complex<susy_QCD::value_type>susy_QCD::ssRsgs=susy_QCD::ssLsgs;
    std::complex<susy_QCD::value_type>susy_QCD::stRsgt=susy_QCD::stLsgt;
    std::complex<susy_QCD::value_type>susy_QCD::sbRsgb=susy_QCD::sbLsgb;

    std::complex<susy_QCD::value_type>susy_QCD::suLusg(0,std::sqrt((susy_QCD::value_type)2)*susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::sdLdsg=susy_QCD::suLusg;
    std::complex<susy_QCD::value_type>susy_QCD::scLcsg=susy_QCD::suLusg;
    std::complex<susy_QCD::value_type>susy_QCD::ssLssg=susy_QCD::suLusg;
    std::complex<susy_QCD::value_type>susy_QCD::stLtsg=susy_QCD::suLusg;
    std::complex<susy_QCD::value_type>susy_QCD::sbLbsg=susy_QCD::suLusg;

    std::complex<susy_QCD::value_type>susy_QCD::suRusg=susy_QCD::suLusg;
    std::complex<susy_QCD::value_type>susy_QCD::sdRdsg=susy_QCD::sdLdsg;
    std::complex<susy_QCD::value_type>susy_QCD::scRcsg=susy_QCD::scLcsg;
    std::complex<susy_QCD::value_type>susy_QCD::ssRssg=susy_QCD::ssLssg;
    std::complex<susy_QCD::value_type>susy_QCD::stRtsg=susy_QCD::stLtsg;
    std::complex<susy_QCD::value_type>susy_QCD::sbRbsg=susy_QCD::sbLbsg;

    std::complex<susy_QCD::value_type>susy_QCD::suLsuLgg(0,susy_QCD::g_s*susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::sdLsdLgg=susy_QCD::suLsuLgg;
    std::complex<susy_QCD::value_type>susy_QCD::scLscLgg=susy_QCD::suLsuLgg;
    std::complex<susy_QCD::value_type>susy_QCD::ssLssLgg=susy_QCD::suLsuLgg;
    std::complex<susy_QCD::value_type>susy_QCD::stLstLgg=susy_QCD::suLsuLgg;
    std::complex<susy_QCD::value_type>susy_QCD::sbLsbLgg=susy_QCD::suLsuLgg;

    std::complex<susy_QCD::value_type>susy_QCD::suRsuRgg=susy_QCD::suLsuLgg;
    std::complex<susy_QCD::value_type>susy_QCD::sdRsdRgg=susy_QCD::sdLsdLgg;
    std::complex<susy_QCD::value_type>susy_QCD::scRscRgg=susy_QCD::scLscLgg;
    std::complex<susy_QCD::value_type>susy_QCD::ssRssRgg=susy_QCD::ssLssLgg;
    std::complex<susy_QCD::value_type>susy_QCD::stRstRgg=susy_QCD::stLstLgg;
    std::complex<susy_QCD::value_type>susy_QCD::sbRsbRgg=susy_QCD::sbLsbLgg;

    std::complex<susy_QCD::value_type>susy_QCD::suLsuLsuLsuL(0,-(susy_QCD::N_c-1)*susy_QCD::g_s*susy_QCD::g_s/susy_QCD::value_type(2*susy_QCD::N_c));
    std::complex<susy_QCD::value_type>susy_QCD::sdLsdLsdLsdL=susy_QCD::suLsuLsuLsuL;
    std::complex<susy_QCD::value_type>susy_QCD::scLscLscLscL=susy_QCD::suLsuLsuLsuL;
    std::complex<susy_QCD::value_type>susy_QCD::ssLssLssLssL=susy_QCD::suLsuLsuLsuL;
    std::complex<susy_QCD::value_type>susy_QCD::stLstLstLstL=susy_QCD::suLsuLsuLsuL;
    std::complex<susy_QCD::value_type>susy_QCD::sbLsbLsbLsbL=susy_QCD::suLsuLsuLsuL;

    std::complex<susy_QCD::value_type>susy_QCD::suRsuRsuRsuR=susy_QCD::suLsuLsuLsuL;
    std::complex<susy_QCD::value_type>susy_QCD::sdRsdRsdRsdR=susy_QCD::sdLsdLsdLsdL;
    std::complex<susy_QCD::value_type>susy_QCD::scRscRscRscR=susy_QCD::scLscLscLscL;
    std::complex<susy_QCD::value_type>susy_QCD::ssRssRssRssR=susy_QCD::ssLssLssLssL;
    std::complex<susy_QCD::value_type>susy_QCD::stRstRstRstR=susy_QCD::stLstLstLstL;
    std::complex<susy_QCD::value_type>susy_QCD::sbRsbRsbRsbR=susy_QCD::sbLsbLsbLsbL;

    std::complex<susy_QCD::value_type>susy_QCD::suLsuLsuRsuR(0,susy_QCD::g_s*susy_QCD::g_s);
    std::complex<susy_QCD::value_type>susy_QCD::sdLsdLsdRsdR=susy_QCD::suLsuLsuRsuR;
    std::complex<susy_QCD::value_type>susy_QCD::scLscLscRscR=susy_QCD::suLsuLsuRsuR;
    std::complex<susy_QCD::value_type>susy_QCD::ssLssLssRssR=susy_QCD::suLsuLsuRsuR;
    std::complex<susy_QCD::value_type>susy_QCD::stLstLstRstR=susy_QCD::suLsuLsuRsuR;
    std::complex<susy_QCD::value_type>susy_QCD::sbLsbLsbRsbR=susy_QCD::suLsuLsuRsuR;

    /* Gauge parameter initialisation: */

    int susy_QCD::gauge=0;

    /* Four-gluon vertex initialisation: */

    bool susy_QCD::four_gluon_vertex=false;

    /* Constructor: */

    susy_QCD::susy_QCD()
    {
	/* Gluon definition: */

	switch(gauge)
	{
	    case 0:
		add_gluon<Feynman_gauge,SU<N_c> >();
		break;
	    case 1:
		add_gluon<unitary_gauge,SU<N_c> >();
		break;
	    case 2:
		add_gluon<R_vector_gauge,SU<N_c> >();
		break;
	    default:
		add_gluon<Feynman_gauge,SU<N_c> >();
	}

	/* Gluino definition: */

	(M_sg>(value_type)0)?add_fermion< adjoint_rep< SU<N_c> > >("~g",&M_sg,&W_sg,1000021):add_fermion< adjoint_rep< SU<N_c> > >("~g",1000021);
	
	/* Quark definitions: */

	(M_d>(value_type)0)?add_quarks< SU<N_c> >("d","dbar",&M_d,1):add_quarks< SU<N_c> >("d","dbar",1);
	(M_u>(value_type)0)?add_quarks< SU<N_c> >("u","ubar",&M_u,2):add_quarks< SU<N_c> >("u","ubar",2);
	(M_s>(value_type)0)?add_quarks< SU<N_c> >("s","sbar",&M_s,3):add_quarks< SU<N_c> >("s","sbar",3);
	(M_c>(value_type)0)?add_quarks< SU<N_c> >("c","cbar",&M_c,4):add_quarks< SU<N_c> >("c","cbar",4);
	(M_b>(value_type)0)?add_quarks< SU<N_c> >("b","bbar",&M_b,5):add_quarks< SU<N_c> >("b","bbar",5);
	(M_t>(value_type)0)?add_quarks< SU<N_c> >("t","tbar",&M_t,&W_t,6):add_quarks< SU<N_c> >("t","tbar",6);

	/* Squark definitions: */

	(M_sdL>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~d_L+","~d_L-",&M_sdL,&W_sdL,1000001):add_scalars< fundamental_rep< SU<N_c> > >("~d_L+","~d_L-",1000001);
	(M_sdR>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~d_R+","~d_R-",&M_sdR,&W_sdR,2000001):add_scalars< fundamental_rep< SU<N_c> > >("~d_R+","~d_R-",2000001);
	(M_suL>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~u_L+","~u_L-",&M_suL,&W_suL,1000002):add_scalars< fundamental_rep< SU<N_c> > >("~u_L+","~u_L-",1000002);
	(M_suR>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~u_R+","~u_R-",&M_suR,&W_suR,2000002):add_scalars< fundamental_rep< SU<N_c> > >("~u_R+","~u_R-",2000002);
	(M_ssL>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~s_L+","~s_L-",&M_ssL,&W_ssL,1000003):add_scalars< fundamental_rep< SU<N_c> > >("~s_L+","~s_L-",1000003);
	(M_ssR>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~s_R+","~s_R-",&M_ssR,&W_ssR,2000003):add_scalars< fundamental_rep< SU<N_c> > >("~s_R+","~s_R-",2000003);
	(M_scL>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~c_L+","~c_L-",&M_scL,&W_scL,1000004):add_scalars< fundamental_rep< SU<N_c> > >("~c_L+","~c_L-",1000004);
	(M_scR>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~c_R+","~c_R-",&M_scR,&W_scR,2000004):add_scalars< fundamental_rep< SU<N_c> > >("~c_R+","~c_R-",2000004);
	(M_sbL>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~b_L+","~b_L-",&M_sbL,&W_sbL,1000005):add_scalars< fundamental_rep< SU<N_c> > >("~b_L+","~b_L-",1000005);
	(M_sbR>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~b_R+","~b_R-",&M_sbR,&W_sbR,2000005):add_scalars< fundamental_rep< SU<N_c> > >("~b_R+","~b_R-",2000005);
	(M_sbL>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~t_L+","~t_L-",&M_stL,&W_stL,1000006):add_scalars< fundamental_rep< SU<N_c> > >("~t_L+","~t_L-",1000006);
	(M_sbR>(value_type)0)?add_scalars< fundamental_rep< SU<N_c> > >("~t_R+","~t_R-",&M_stR,&W_stR,2000006):add_scalars< fundamental_rep< SU<N_c> > >("~t_R+","~t_R-",2000006);
	
	/* Triple gluon vertex definition: */

	add_vertex<colour_tensor::f< SU<N_c> >,vvv>("g","g","g",&ggg);

	/* Gluon-squark-squark vertices: */

	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~d_L+","~d_L-",&gsdLsdL);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~d_R+","~d_R-",&gsdRsdR);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~u_L+","~u_L-",&gsuLsuL);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~u_R+","~u_R-",&gsuRsuR);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~s_L+","~s_L-",&gssLssL);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~s_R+","~s_R-",&gssRssR);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~c_L+","~c_L-",&gscLscL);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~c_R+","~c_R-",&gscRscR);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~b_L+","~b_L-",&gsbLsbL);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~b_R+","~b_R-",&gsbRsbR);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~t_L+","~t_L-",&gstLstL);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,0,2,1>,vss>("g","~t_R+","~t_R-",&gstRstR);

	/* Gluon-quark-quark vertices: */
	
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","dbar","d",&gdd);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","ubar","u",&guu);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","sbar","s",&gss);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","cbar","c",&gcc);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","bbar","b",&gbb);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","tbar","t",&gtt);

	/* Gluon-gluino-gluino vertex: */

	add_vertex<colour_tensor::f< SU<N_c> >,vff>("g","~g","~g",&ggg);

	/* Squark-gluino-quark vertices: */

	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffL>("~d_L+","~g","d",&sdLsgd);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffR>("~d_L-","dbar","~g",&sdLdsg);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffR>("~d_R+","~g","d",&sdRsgd);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffL>("~d_R-","dbar","~g",&sdRdsg);

	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffL>("~u_L-","~g","u",&suLsgu);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffR>("~u_L+","ubar","~g",&suLusg);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffR>("~u_R-","~g","u",&suRsgu);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffL>("~u_R+","ubar","~g",&suRusg);
	
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffL>("~s_L+","~g","s",&ssLsgs);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffR>("~s_L-","sbar","~g",&ssLssg);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffR>("~s_R+","~g","s",&ssRsgs);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffL>("~s_R-","sbar","~g",&ssRssg);
	
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffL>("~c_L-","~g","c",&scLsgc);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffR>("~c_L+","cbar","~g",&scLcsg);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffR>("~c_R-","~g","c",&scRsgc);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffL>("~c_R+","cbar","~g",&scRcsg);
	
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffL>("~b_L+","~g","b",&sbLsgb);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffR>("~b_L-","bbar","~g",&sbLbsg);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffR>("~b_R+","~g","b",&sbRsgb);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffL>("~b_R-","bbar","~g",&sbRbsg);
	
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffL>("~t_L-","~g","t",&stLsgt);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffR>("~t_L+","tbar","~g",&stLtsg);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,1,0,2>,sffR>("~t_R-","~g","t",&stRsgt);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> >,2,1,0>,sffL>("~t_R+","tbar","~g",&stRtsg);
	
	/* 4-gluon vertex definition: */
	
	if(four_gluon_vertex)
	{
	    add_vertex<colour_tensor::ff_contr< SU<N_c> >,vvvv>("g","g","g","g",&gggg);
	}
	else
	{
	    add_fast_4g_vertex< SU<N_c> >("g",&Tgg);
	}

	/* Squark-squark-gluon-gluon vertices: */

	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~d_L-","~d_L+","g","g",&sdLsdLgg);
	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~d_R-","~d_R+","g","g",&sdRsdRgg);

	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~u_L-","~u_L+","g","g",&suLsuLgg);
	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~u_R-","~u_R+","g","g",&suRsuRgg);

	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~s_L-","~s_L+","g","g",&ssLssLgg);
	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~s_R-","~s_R+","g","g",&ssRssRgg);

	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~c_L-","~c_L+","g","g",&scLscLgg);
	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~c_R-","~c_R+","g","g",&scRscRgg);

	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~b_L-","~b_L+","g","g",&sbLsbLgg);
	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~b_R-","~b_R+","g","g",&sbRsbRgg);

	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~t_L-","~t_L+","g","g",&stLstLgg);
	add_vertex<colour_tensor::TT_plus<fundamental_rep< SU<N_c> >,2,3,0,1>,ssvv>("~t_R-","~t_R+","g","g",&stRstRgg);

	/* 4-Squark vertices: */

	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~d_L+","~d_L-","~d_L+","~d_L-",&sdLsdLsdLsdL);
	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~d_R+","~d_R-","~d_R+","~d_R-",&sdRsdRsdRsdR);
	add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<N_c> > >,ssss>("~d_L+","~d_L-","~d_R+","~d_R-",&sdLsdLsdRsdR);

	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~u_L+","~u_L-","~u_L+","~u_L-",&suLsuLsuLsuL);
	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~u_R+","~u_R-","~u_R+","~u_R-",&suRsuRsuRsuR);
	add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<N_c> > >,ssss>("~u_L+","~u_L-","~u_R+","~u_R-",&suLsuLsuRsuR);

	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~s_L+","~s_L-","~s_L+","~s_L-",&ssLssLssLssL);
	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~s_R+","~s_R-","~s_R+","~s_R-",&ssRssRssRssR);
	add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<N_c> > >,ssss>("~s_L+","~s_L-","~s_R+","~s_R-",&ssLssLssRssR);

	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~c_L+","~c_L-","~c_L+","~c_L-",&scLscLscLscL);
	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~c_R+","~c_R-","~c_R+","~c_R-",&scRscRscRscR);
	add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<N_c> > >,ssss>("~c_L+","~c_L-","~c_R+","~c_R-",&scLscLscRscR);

	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~b_L+","~b_L-","~b_L+","~b_L-",&sbLsbLsbLsbL);
	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~b_R+","~b_R-","~b_R+","~b_R-",&sbRsbRsbRsbR);
	add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<N_c> > >,ssss>("~b_L+","~b_L-","~b_R+","~b_R-",&sbLsbLsbRsbR);

	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~t_L+","~t_L-","~t_L+","~t_L-",&stLstLstLstL);
	add_vertex<colour_tensor::dd_plus<fundamental_rep< SU<N_c> > >,ssss>("~t_R+","~t_R-","~t_R+","~t_R-",&stRstRstRstR);
	add_vertex<colour_tensor::TT_contr<fundamental_rep< SU<N_c> > >,ssss>("~t_L+","~t_L-","~t_R+","~t_R-",&stLstLstRstR);

	/* Quark family definition: */

	construct_family("q","u,d,c,s");

	/* Quark family definition including b's: */

	construct_family("Q","u,d,c,s,b");

	/* Quark family definition including b's and tops: */

	construct_family("Qt","u,d,c,s,b,t");

	/* Antiquark family definition: */

	construct_family("qbar","dbar,ubar,sbar,cbar");

	/* Antiquark family definition, including b's: */

	construct_family("Qbar","dbar,ubar,sbar,cbar,bbar");

	/* Antiquark family definition, including b's and tops: */

	construct_family("Qtbar","dbar,ubar,sbar,cbar,bbar,tbar");

	/* Up-type quark family definition: */

	construct_family("q_up","u,c");

	/* Up-type quark family definition, including b's: */

	construct_family("Q_up","u,c");

	/* Up-type quark family definition, including b's and tops: */

	construct_family("Qt_up","u,c,t");

	/* Up-type antiquark family definition: */

	construct_family("qbar_up","ubar,cbar");

	/* Up-type antiquark family definition, including b's: */

	construct_family("Qbar_up","ubar,cbar");

	/* Up-type antiquark family definition, including b's and tops: */

	construct_family("Qtbar_up","ubar,cbar,tbar");

	/* Down-type quark family definition: */

	construct_family("q_down","d,s");

	/* Down-type quark family definition, including b's: */

	construct_family("Q_down","d,s,b");

	/* Down-type quark family definition, including b's and tops: */

	construct_family("Qt_down","d,s,b");

	/* Down-type antiquark family definition: */

	construct_family("qbar_down","dbar,sbar");

	/* Down-type antiquark family definition, including b's: */

	construct_family("Qbar_down","dbar,sbar,bbar");

	/* Down-type antiquark family definition, including b's and tops: */

	construct_family("Qtbar_down","dbar,sbar,bbar");

	/* Positively-charged quark family definition: */

	construct_family("q+","u,dbar,c,sbar");

	/* Positively-charged quark family definition including b's: */

	construct_family("Q+","u,dbar,c,sbar,bbar");

	/* Positively-charged quark family definition including b's and
	 * tops: */

	construct_family("Qt+","u,dbar,c,sbar,bbar,t");

	/* Negatively-charged quark family definition: */

	construct_family("q-","d,ubar,s,cbar");

	/* Negatively-charged quark family definition, including b's: */

	construct_family("Q-","d,ubar,s,cbar,b");

	/* Negatively-charged quark family definition, including b's and
	 * tops: */

	construct_family("Qt-","d,ubar,s,cbar,b,tbar");

	/* Proton constituent parton family definition: */

	construct_family("p","u,ubar,d,dbar,c,cbar,s,sbar,g");

	/* Proton constituent parton family definition, including b's: */

	construct_family("P","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

	/* QCD jet parton family definition: */

	construct_family("j","u,ubar,d,dbar,c,cbar,s,sbar,g");

	/* QCD jet parton family definition, including b's: */

	construct_family("J","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar,g");

	/* QCD jet parton family definition, including b's and tops: */

	construct_family("Jt","u,ubar,d,dbar,c,cbar,s,sbar,t,tbar,b,bbar,g");

	/* Left-handed positive squark family definitions: */

	construct_family("~Q_L+","~d_L+,~u_L+,~s_L+,~c_L+,~b_L+,~t_L+");

	/* Left-handed negative squark family definitions: */

	construct_family("~Q_L-","~d_L-,~u_L-,~s_L-,~c_L-,~b_L-,~t_L-");

	/* Right-handed positive squark family definitions: */

	construct_family("~Q_R+","~d_R+,~u_R+,~s_R+,~c_R+,~b_R+,~t_R+");

	/* Right-handed negative squark family definitions: */

	construct_family("~Q_R-","~d_R-,~u_R-,~s_R-,~c_R-,~b_R-,~t_R-");

	/* Left-handed squark family definitions: */

	construct_family("~Q_L","~d_L+,~d_L-,~u_L+,~u_L-,~s_L+,~s_L-,~c_L+,~c_L-,~b_L+,~b_L-,~t_L+,~t_L-");

	/* Right-handed squark family definitions: */

	construct_family("~Q_R","~d_R+,~d_R-,~u_R+,~u_R-,~s_R+,~s_R-,~c_R+,~c_R-,~b_R+,~b_R-,~t_R+,~t_R-");

	/* Positive squark family definitions: */

	construct_family("~Q+","~d_L+,~d_R+,~u_L+,~u_R+,~s_L+,~s_R+,~c_L+,~c_R+,~b_L+,~b_R+,~t_L+,~t_R+");

	/* Negative squark family definitions: */

	construct_family("~Q-","~d_L-,~d_R-,~u_L-,~u_R-,~s_L-,~s_R-,~c_L-,~c_R-,~b_L-,~b_R-,~t_L-,~t_R-");
    }

    /* Function computing the couplings from the current alpha_s value:
     * */
	    
    void susy_QCD::refresh_couplings()
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

	suLsuLsuLsuL=std::complex<value_type>(0,-(value_type)(N_c-1)*g_s*g_s/value_type(2*N_c));
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

    void susy_QCD::refresh_masses()
    {
	if(initialised())
	{
	    (M_d>(value_type)0)?set_mass("d",&M_d):set_massless("d");
	    (M_u>(value_type)0)?set_mass("u",&M_u):set_massless("u");
	    (M_s>(value_type)0)?set_mass("s",&M_s):set_massless("s");
	    (M_c>(value_type)0)?set_mass("c",&M_c):set_massless("c");
	    (M_b>(value_type)0)?set_mass("b",&M_b):set_massless("b");
	    (M_t>(value_type)0)?set_mass("t",&M_t):set_massless("t");

	    (M_sdL>(value_type)0)?set_mass("~d_L",&M_sdL):set_massless("~d_L");
	    (M_suL>(value_type)0)?set_mass("~u_L",&M_suL):set_massless("~u_L");
	    (M_ssL>(value_type)0)?set_mass("~s_L",&M_ssL):set_massless("~s_L");
	    (M_scL>(value_type)0)?set_mass("~c_L",&M_scL):set_massless("~c_L");
	    (M_sbL>(value_type)0)?set_mass("~b_L",&M_sbL):set_massless("~b_L");
	    (M_stL>(value_type)0)?set_mass("~t_L",&M_stL):set_massless("~t_L");

	    (M_sdR>(value_type)0)?set_mass("~d_R",&M_sdR):set_massless("~d_R");
	    (M_suR>(value_type)0)?set_mass("~u_R",&M_suR):set_massless("~u_R");
	    (M_ssR>(value_type)0)?set_mass("~s_R",&M_ssR):set_massless("~s_R");
	    (M_scR>(value_type)0)?set_mass("~c_R",&M_scR):set_massless("~c_R");
	    (M_sbR>(value_type)0)?set_mass("~b_R",&M_sbR):set_massless("~b_R");
	    (M_stR>(value_type)0)?set_mass("~t_R",&M_stR):set_massless("~t_R");

	    (M_sg>(value_type)0)?set_mass("~g",&M_sg):set_massless("~g");
	}
    }

    /* Function refreshing the particle decay widths: */

    void susy_QCD::refresh_widths()
    {
	if(initialised())
	{
	    (W_t>(value_type)0)?set_width("t",&W_t):set_widthless("t");
	    (W_sdL>(value_type)0)?set_width("~d_L",&W_sdL):set_widthless("~d_L");
	    (W_suL>(value_type)0)?set_width("~u_L",&W_suL):set_widthless("~u_L");
	    (W_ssL>(value_type)0)?set_width("~s_L",&W_ssL):set_widthless("~s_L");
	    (W_scL>(value_type)0)?set_width("~c_L",&W_scL):set_widthless("~c_L");
	    (W_sbL>(value_type)0)?set_width("~b_L",&W_sbL):set_widthless("~b_L");
	    (W_stL>(value_type)0)?set_width("~t_L",&W_stL):set_widthless("~t_L");
	    (W_sdR>(value_type)0)?set_width("~d_R",&W_sdR):set_widthless("~d_R");
	    (W_suR>(value_type)0)?set_width("~u_R",&W_suR):set_widthless("~u_R");
	    (W_ssR>(value_type)0)?set_width("~s_R",&W_ssR):set_widthless("~s_R");
	    (W_scR>(value_type)0)?set_width("~c_R",&W_scR):set_widthless("~c_R");
	    (W_sbR>(value_type)0)?set_width("~b_R",&W_sbR):set_widthless("~b_R");
	    (W_stR>(value_type)0)?set_width("~t_R",&W_stR):set_widthless("~t_R");
	    (W_sg>(value_type)0)?set_width("~g",&W_sg):set_widthless("~g");
	}
    }

    /* Function returning the number of QCD colours: */

    std::size_t susy_QCD::QCD_colours()
    {
	return N_c;
    }

    /* Dummy function setting the EM fine-structure constant: */

    void susy_QCD::set_alpha(const value_type& a){}

    /* Function setting alpha strong and computing the couplings: */

    void susy_QCD::set_alpha_s(const value_type& a)
    {
	alpha_s=a;
	refresh_couplings();
    }

    /* Function setting the scale and computing alpha_s (at LO) and the
     * couplings: */

    void susy_QCD::set_QCD_scale(const value_type& mu)
    {
	QCD_scale=mu;
	value_type rho=mu/(value_type)SM_params::QCD_scale;
	alpha_s=SM_params::alpha_s/((value_type)1+(value_type)SM_params::alpha_s*value_type(11*N_c-12)*std::log(rho)/((value_type)6*pi));
	refresh_couplings();
    }

    /* Function setting the gluon-propagator to the Feynman gauge: */

    void susy_QCD::set_Feynman_gauge()
    {
	if(initialised() and gauge!=0)
	{
	    set_propagator<Feynman_gauge>("g");
	}
	gauge=0;
    }

    /* Function setting the gluon-propagator to the unitary gauge: */

    void susy_QCD::set_unitary_gauge()
    {
	if(initialised() and gauge!=1)
	{
	    set_propagator<unitary_gauge>("g");
	}
	gauge=1;
    }

    /* Function setting the gluon-propagator to the R-xi gauge: */

    void susy_QCD::set_R_xi_gauge()
    {
	if(initialised() and gauge!=2)
	{
	    set_propagator<R_vector_gauge>("g");
	}
	gauge=2;
    }
	    
    /* Function setting the xi-parameter: */
	    
    void susy_QCD::set_xi(const value_type& x)
    {
	R_gauge<susy_QCD>::xi=x;
    }

    /* Function switching to a 4-gluon vertex Lagrangian: */

    void susy_QCD::set_4_gluon_vertex()
    {
	if(initialised() and !four_gluon_vertex)
	{
	    erase_vertex("H_qcd","g","g");
	    add_vertex<colour_tensor::ff_contr< SU<N_c> >,vvvv>("g","g","g","g",&gggg);
	}
	four_gluon_vertex=true;
    }

    /* Function switching to an auxiliary tensor field description of
     * the 4-gluon vertex: */

    void susy_QCD::set_auxiliary_QCD_field()
    {
	if(initialised() and four_gluon_vertex)
	{
	    erase_vertex("g","g","g","g");
	    add_fast_4g_vertex< SU<N_c> >("g",&Tgg);
	}
	four_gluon_vertex=false;
    }
}

