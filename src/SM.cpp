//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/SM.h>
#include <Camgen/scalar_particle.h>
#include <Camgen/vector_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/vff.h>
#include <Camgen/sff.h>
#include <Camgen/vffVA.h>
#include <Camgen/vffR.h>
#include <Camgen/svv.h>
#include <Camgen/vss.h>
#include <Camgen/ssvv.h>
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

namespace Camgen
{
    /* Compile-time constant definitions: */

    const std::size_t SM::dimension;
    const bool SM::coloured;
    const bool SM::continuous_colours;
    const bool SM::continuous_helicities;
    const int SM::beam_direction;
    const std::size_t SM::N_c;

    const SM::value_type SM::pi=std::acos(-(SM::value_type)1);
    
    /* Real input parameters: */

    SM::value_type SM::M_e=(SM::value_type)0;
    SM::value_type SM::M_mu=(SM::value_type)0;
    SM::value_type SM::M_tau=(SM::value_type)SM_params::M_tau;
    SM::value_type SM::M_u=(SM::value_type)0;
    SM::value_type SM::M_d=(SM::value_type)0;
    SM::value_type SM::M_c=(SM::value_type)SM_params::M_c;
    SM::value_type SM::M_s=(SM::value_type)0;
    SM::value_type SM::M_t=(SM::value_type)SM_params::M_t;
    SM::value_type SM::M_b=(SM::value_type)SM_params::M_b;
    SM::value_type SM::M_Z=(SM::value_type)SM_params::M_Z;
    SM::value_type SM::M_h0=(SM::value_type)120;
    SM::value_type SM::alpha=(SM::value_type)SM_params::alpha;
    SM::value_type SM::G_F=(SM::value_type)SM_params::G_F;
    SM::value_type SM::QCD_scale=(SM::value_type)SM_params::QCD_scale;
    SM::value_type SM::alpha_s=(SM::value_type)SM_params::alpha_s;

    SM::value_type SM::s12=(SM::value_type)SM_params::s12;
    SM::value_type SM::s23=(SM::value_type)SM_params::s23;
    SM::value_type SM::s13=(SM::value_type)SM_params::s13;
    SM::value_type SM::delta=(SM::value_type)0;

    /* Computed parameters: */
    
    SM::value_type SM::Q_e=std::sqrt((SM::value_type)4*SM::pi*SM::alpha);
    SM::value_type SM::g_s=std::sqrt((SM::value_type)4*SM::pi*SM::alpha_s);
    SM::value_type SM::M_W=SM::M_Z*std::sqrt((SM::value_type)0.5+(SM::value_type)0.5*std::sqrt((SM::value_type)1-(SM::value_type)4*SM::pi*SM::alpha/(std::sqrt((SM::value_type)2)*SM::G_F*SM::M_Z*SM::M_Z)));
    SM::value_type SM::W_W=(3+2*N_c)*SM::G_F*std::pow(SM::M_W,3)/(SM::pi*std::sqrt((SM::value_type)2))/(SM::value_type)6;
    SM::value_type SM::W_Z=((SM::value_type)(162+49*N_c)-(SM::value_type)(324+92*N_c)*std::pow(SM::M_W/SM::M_Z,2)+(SM::value_type)(216+88*N_c)*std::pow(SM::M_W/SM::M_Z,4))*SM::G_F*std::pow(SM::M_Z,3)/((SM::value_type)108*SM::pi*std::sqrt((SM::value_type)2));
    SM::value_type SM::W_h0=N_c*SM::G_F*SM::M_b*SM::M_b*SM::M_h0*std::pow((SM::value_type)1-std::pow((SM::value_type)2*SM::M_b/SM::M_h0,2),1.5)/(std::sqrt((SM::value_type)2)*(SM::value_type)4*SM::pi);
    SM::value_type SM::W_t=SM::G_F*std::pow(SM::M_t,3)*((SM::value_type)1-std::pow(SM::M_W/SM::M_t,2))*((SM::value_type)1+std::pow(SM::M_W/SM::M_t,2)-2*std::pow(SM::M_W/SM::M_t,4))/(std::sqrt((SM::value_type)128)*SM::pi);

    std::complex<SM::value_type>SM::cM_Z=std::sqrt(std::complex<SM::value_type>(SM::M_Z*SM::M_Z,-SM::M_Z*SM::W_Z));
    std::complex<SM::value_type>SM::cM_W=std::sqrt(std::complex<SM::value_type>(SM::M_W*SM::M_W,-SM::M_W*SM::W_W));
    std::complex<SM::value_type>SM::cos_W=SM::cM_W/SM::cM_Z;
    std::complex<SM::value_type>SM::sin_W=std::sqrt(std::complex<SM::value_type>(1,0)-SM::cos_W*SM::cos_W);
    std::complex<SM::value_type>SM::tan_W=SM::sin_W/SM::cos_W;
    std::complex<SM::value_type>SM::cos2_W=SM::cos_W*SM::cos_W-SM::sin_W*SM::sin_W;
    std::complex<SM::value_type>SM::sin2_W=(SM::value_type)2*SM::sin_W*SM::cos_W;
    
    SM::value_type SM::c12=std::sqrt(((SM::value_type)1+SM::s12)*((SM::value_type)1-SM::s12));
    SM::value_type SM::c23=std::sqrt(((SM::value_type)1+SM::s23)*((SM::value_type)1-SM::s23));
    SM::value_type SM::c13=std::sqrt(((SM::value_type)1+SM::s13)*((SM::value_type)1-SM::s13));
    
    std::complex<SM::value_type>SM::V_CKM[3][3]={{SM::c12*SM::c13,			   SM::s12*SM::c13,			    SM::s13},
						 {-SM::s12*SM::c23-SM::c12*SM::s23*SM::s13,SM::c12*SM::c23-SM::s12*SM::s23*SM::s13, SM::s23*SM::c13},
						 {SM::s12*SM::s23-SM::c12*SM::c23*SM::s13, -SM::c12*SM::s23-SM::s12*SM::c23*SM::s13,SM::c23*SM::c13}};
    
    /* Gauge switch initialisation: */

    int SM::gauge=0;

    /* Four-gluon vertex switch initialisation: */

    bool SM::four_gluon_vertex=false;

    /* VVVV couplings: */

    std::complex<SM::value_type>SM::WWWW=std::complex<SM::value_type>(0,SM::Q_e*SM::Q_e)/(SM::sin_W*SM::sin_W);
    std::complex<SM::value_type>SM::WZWZ=-SM::WWWW*SM::cos_W*SM::cos_W;
    std::complex<SM::value_type>SM::WgammaWZ=SM::WWWW*SM::cos_W*SM::sin_W;
    std::complex<SM::value_type>SM::WgammaWgamma=-SM::WWWW*SM::sin_W*SM::sin_W;
    std::complex<SM::value_type>SM::gggg(0,-SM::g_s*SM::g_s);

    /* VVV couplings: */

    std::complex<SM::value_type>SM::gammaWW(0,SM::Q_e);
    std::complex<SM::value_type>SM::ZWW=-SM::gammaWW/SM::tan_W;
    std::complex<SM::value_type>SM::ggg(SM::g_s,0);

    /* SSSS couplings: */

    std::complex<SM::value_type>SM::hhhh=-(SM::value_type)0.75*std::complex<SM::value_type>(0,SM::Q_e*SM::Q_e)*std::pow(SM::M_h0/(SM::sin_W*SM::cM_W),2);
    std::complex<SM::value_type>SM::hhchichi=SM::hhhh/(SM::value_type)3;
    std::complex<SM::value_type>SM::hhphiphi=SM::hhchichi;
    std::complex<SM::value_type>SM::chichichichi=SM::hhhh;
    std::complex<SM::value_type>SM::chichiphiphi=SM::hhphiphi;
    std::complex<SM::value_type>SM::phiphiphiphi=(SM::value_type)2*SM::hhchichi;

    /* SSS couplings: */

    std::complex<SM::value_type>SM::hhh=-(SM::value_type)1.5*std::complex<SM::value_type>(0,SM::Q_e)*SM::M_h0*SM::M_h0/(SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::hchichi=SM::hhh/(SM::value_type)3;
    std::complex<SM::value_type>SM::hphiphi=SM::hchichi;
    
    /* VFF couplings: */

    std::complex<SM::value_type>SM::gammaee(0,SM::Q_e);
    std::complex<SM::value_type>SM::gammamumu=SM::gammaee;
    std::complex<SM::value_type>SM::gammatautau=SM::gammaee;
    std::complex<SM::value_type>SM::gammauu(0,-(SM::value_type)2*SM::Q_e/(SM::value_type)3);
    std::complex<SM::value_type>SM::gammadd(0,SM::Q_e/(SM::value_type)3);
    std::complex<SM::value_type>SM::gammacc=SM::gammauu;
    std::complex<SM::value_type>SM::gammass=SM::gammadd;
    std::complex<SM::value_type>SM::gammatt=SM::gammauu;
    std::complex<SM::value_type>SM::gammabb=SM::gammadd;
    
    std::complex<SM::value_type>SM::Znene=SM::gammaee/SM::sin2_W;
    std::complex<SM::value_type>SM::Zee_V=SM::Znene*((SM::value_type)2*SM::sin_W*SM::sin_W-(SM::value_type)0.5);
    std::complex<SM::value_type>SM::Zee_A=(SM::value_type)0.5*SM::Znene;
    std::complex<SM::value_type>SM::Znmunmu=SM::Znene;
    std::complex<SM::value_type>SM::Zmumu_V=SM::Zee_V;
    std::complex<SM::value_type>SM::Zmumu_A=SM::Zee_A;
    std::complex<SM::value_type>SM::Zntauntau=SM::Znene;
    std::complex<SM::value_type>SM::Ztautau_V=SM::Zee_V;
    std::complex<SM::value_type>SM::Ztautau_A=SM::Zee_A;
    std::complex<SM::value_type>SM::Zuu_V=SM::Znene*((SM::value_type)0.5-(SM::value_type)4*SM::sin_W*SM::sin_W/(SM::value_type)3);
    std::complex<SM::value_type>SM::Zuu_A=-SM::Zee_A;
    std::complex<SM::value_type>SM::Zdd_V=SM::Znene*(-(SM::value_type)0.5+(SM::value_type)2*SM::sin_W*SM::sin_W/(SM::value_type)3);
    std::complex<SM::value_type>SM::Zdd_A=SM::Zee_A;
    std::complex<SM::value_type>SM::Zcc_V=SM::Zuu_V;
    std::complex<SM::value_type>SM::Zcc_A=SM::Zuu_A;
    std::complex<SM::value_type>SM::Zss_V=SM::Zdd_V;
    std::complex<SM::value_type>SM::Zss_A=SM::Zdd_A;
    std::complex<SM::value_type>SM::Ztt_V=SM::Zuu_V;
    std::complex<SM::value_type>SM::Ztt_A=SM::Zuu_A;
    std::complex<SM::value_type>SM::Zbb_V=SM::Zdd_V;
    std::complex<SM::value_type>SM::Zbb_A=SM::Zdd_A;

    std::complex<SM::value_type>SM::Wnee=SM::gammaee/(std::sqrt((SM::value_type)2)*SM::sin_W);
    std::complex<SM::value_type>SM::Wnmumu=SM::Wnee;
    std::complex<SM::value_type>SM::Wntautau=SM::Wnee;
    std::complex<SM::value_type>SM::Wene=SM::Wnee;
    std::complex<SM::value_type>SM::Wmunmu=SM::Wnmumu;
    std::complex<SM::value_type>SM::Wtauntau=SM::Wntautau;
    std::complex<SM::value_type>SM::Wud=SM::Wnee*SM::V_CKM[0][0];
    std::complex<SM::value_type>SM::Wdu=SM::Wud;
    std::complex<SM::value_type>SM::Wus=SM::Wnee*SM::V_CKM[0][1];
    std::complex<SM::value_type>SM::Wsu=SM::Wus;
    std::complex<SM::value_type>SM::Wub=(SM::value_type)0;
    std::complex<SM::value_type>SM::Wbu=(SM::value_type)0;
    std::complex<SM::value_type>SM::Wcd=SM::Wnee*SM::V_CKM[1][0];
    std::complex<SM::value_type>SM::Wdc=SM::Wcd;
    std::complex<SM::value_type>SM::Wcs=SM::Wnee*SM::V_CKM[1][1];
    std::complex<SM::value_type>SM::Wsc=SM::Wcs;
    std::complex<SM::value_type>SM::Wcb=SM::Wnee*SM::V_CKM[1][2];
    std::complex<SM::value_type>SM::Wbc=SM::Wcb;
    std::complex<SM::value_type>SM::Wtd=(SM::value_type)0;
    std::complex<SM::value_type>SM::Wdt=(SM::value_type)0;
    std::complex<SM::value_type>SM::Wts=SM::Wnee*SM::V_CKM[2][1];
    std::complex<SM::value_type>SM::Wst=SM::Wts;
    std::complex<SM::value_type>SM::Wtb=SM::Wnee*SM::V_CKM[2][2];
    std::complex<SM::value_type>SM::Wbt=SM::Wtb;
    
    std::complex<SM::value_type>SM::guu(0,SM::g_s);
    std::complex<SM::value_type>SM::gdd(0,SM::g_s);
    std::complex<SM::value_type>SM::gcc(0,SM::g_s);
    std::complex<SM::value_type>SM::gss(0,SM::g_s);
    std::complex<SM::value_type>SM::gtt(0,SM::g_s);
    std::complex<SM::value_type>SM::gbb(0,SM::g_s);

    /* SFF couplings: */

    std::complex<SM::value_type>SM::hee=(SM::value_type)0;
    std::complex<SM::value_type>SM::hmumu=(SM::value_type)0;
    std::complex<SM::value_type>SM::htautau=-(SM::value_type)0.5*SM::gammaee*SM::M_tau/(SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::huu=(SM::value_type)0;
    std::complex<SM::value_type>SM::hdd=(SM::value_type)0;
    std::complex<SM::value_type>SM::hcc=-(SM::value_type)0.5*SM::gammaee*SM::M_c/(SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::hss=(SM::value_type)0;
    std::complex<SM::value_type>SM::htt=-(SM::value_type)0.5*SM::gammaee*SM::M_t/(SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::hbb=-(SM::value_type)0.5*SM::gammaee*SM::M_b/(SM::sin_W*SM::cM_W);

    std::complex<SM::value_type>SM::chiee=(SM::value_type)0;
    std::complex<SM::value_type>SM::chimumu=(SM::value_type)0;
    std::complex<SM::value_type>SM::chitautau=(SM::value_type)0.5*SM::Q_e*SM::M_tau/(SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::chiuu=(SM::value_type)0;
    std::complex<SM::value_type>SM::chidd=(SM::value_type)0;
    std::complex<SM::value_type>SM::chicc=-(SM::value_type)0.5*SM::Q_e*SM::M_t/(SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::chiss=(SM::value_type)0;
    std::complex<SM::value_type>SM::chitt=-(SM::value_type)0.5*SM::Q_e*SM::M_t/(SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::chibb=(SM::value_type)0.5*SM::Q_e*SM::M_b/(SM::sin_W*SM::cM_W);

    std::complex<SM::value_type>SM::phinee=(SM::value_type)0;
    std::complex<SM::value_type>SM::phiene=(SM::value_type)0;
    std::complex<SM::value_type>SM::phinmumu=(SM::value_type)0;
    std::complex<SM::value_type>SM::phimunmu=(SM::value_type)0;
    std::complex<SM::value_type>SM::phintautau=-SM::gammaee*SM::M_tau/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phitauntau=SM::phintautau;
    
    std::complex<SM::value_type>SM::phiud_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phiud_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phidu_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phidu_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phius_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phius_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phisu_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phisu_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phiub_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phiub_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phibu_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phibu_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phicd_S=(SM::value_type)0.5*SM::gammaee*SM::M_c*SM::V_CKM[1][0]/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phicd_A=-SM::phicd_S;
    std::complex<SM::value_type>SM::phidc_S=SM::phicd_S;
    std::complex<SM::value_type>SM::phidc_A=-SM::phicd_A;
    std::complex<SM::value_type>SM::phics_S=(SM::value_type)0.5*SM::gammaee*SM::M_c*SM::V_CKM[1][1]/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phics_A=-SM::phics_S;
    std::complex<SM::value_type>SM::phisc_S=SM::phics_S;
    std::complex<SM::value_type>SM::phisc_A=-SM::phics_A;
    std::complex<SM::value_type>SM::phicb_S=(SM::value_type)0.5*SM::gammaee*(SM::M_c-SM::M_b)*SM::V_CKM[1][2]/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phicb_A=-(SM::value_type)0.5*SM::gammaee*(SM::M_c+SM::M_b)*SM::V_CKM[1][2]/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phibc_S=SM::phicb_S;
    std::complex<SM::value_type>SM::phibc_A=-SM::phicb_A;
    std::complex<SM::value_type>SM::phitd_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phitd_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phidt_S=(SM::value_type)0;
    std::complex<SM::value_type>SM::phidt_A=(SM::value_type)0;
    std::complex<SM::value_type>SM::phits_S=(SM::value_type)0.5*SM::gammaee*SM::M_t*SM::V_CKM[2][1]/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phits_A=-SM::phits_S;
    std::complex<SM::value_type>SM::phist_S=SM::phits_S;
    std::complex<SM::value_type>SM::phist_A=-SM::phits_A;
    std::complex<SM::value_type>SM::phitb_S=(SM::value_type)0.5*SM::gammaee*(SM::M_t-SM::M_b)/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phitb_A=-(SM::value_type)0.5*SM::gammaee*(SM::M_t+SM::M_b)/(std::sqrt((SM::value_type)2)*SM::sin_W*SM::cM_W);
    std::complex<SM::value_type>SM::phibt_S=SM::phitb_S;
    std::complex<SM::value_type>SM::phibt_A=-SM::phitb_A;

    /* VSS couplings: */

    std::complex<SM::value_type>SM::Zchih=SM::Q_e/SM::sin2_W;
    std::complex<SM::value_type>SM::gammaphiphi=-SM::gammaee;
    std::complex<SM::value_type>SM::Zphiphi=SM::gammaee*SM::cos2_W/SM::sin2_W;
    std::complex<SM::value_type>SM::Wpphimh=-(SM::value_type)0.5*SM::gammaee/SM::sin_W;
    std::complex<SM::value_type>SM::Wmphiph=-SM::Wpphimh;
    std::complex<SM::value_type>SM::Wpphimchi=(SM::value_type)0.5*SM::Q_e/SM::sin_W;
    std::complex<SM::value_type>SM::Wmphipchi=SM::Wpphimchi;
    
    /* SVV couplings: */

    std::complex<SM::value_type>SM::hWW=SM::gammaee*SM::cM_W/SM::sin_W;
    std::complex<SM::value_type>SM::hZZ=(SM::value_type)2*SM::gammaee*SM::cM_Z/SM::sin2_W;
    std::complex<SM::value_type>SM::phipWmZ=-SM::gammaee*SM::cM_W*SM::tan_W;
    std::complex<SM::value_type>SM::phimWpZ=SM::phipWmZ;
    std::complex<SM::value_type>SM::phipWmgamma=-SM::gammaee*SM::cM_W;
    std::complex<SM::value_type>SM::phimWpgamma=SM::phipWmgamma;

    /* SSVV couplings: */

    std::complex<SM::value_type>SM::hhWW=std::complex<SM::value_type>((SM::value_type)0,(SM::value_type)0.5)*std::pow(SM::Q_e/SM::sin_W,2);
    std::complex<SM::value_type>SM::chichiWW=SM::hhWW;
    std::complex<SM::value_type>SM::phiphiWW=SM::hhWW;
    std::complex<SM::value_type>SM::hhZZ=SM::hhWW/(SM::cos_W*SM::cos_W);
    std::complex<SM::value_type>SM::chichiZZ=SM::hhZZ;
    std::complex<SM::value_type>SM::phiphiZZ=SM::hhZZ*SM::cos2_W*SM::cos2_W;
    std::complex<SM::value_type>SM::phiphigammaZ=-(SM::value_type)2*SM::hhWW*SM::cos2_W*SM::tan_W;
    std::complex<SM::value_type>SM::phiphigammagamma=(SM::value_type)2*SM::Q_e*SM::gammaee;
    std::complex<SM::value_type>SM::phiphWmZ=-(SM::value_type)0.5*SM::Q_e*SM::gammaee/SM::cos_W;
    std::complex<SM::value_type>SM::phimhWpZ=SM::phiphWmZ;
    std::complex<SM::value_type>SM::phiphWmgamma=SM::phimhWpZ/SM::tan_W;
    std::complex<SM::value_type>SM::phimhWpgamma=SM::phiphWmgamma;
    std::complex<SM::value_type>SM::phimchiWpZ=std::complex<SM::value_type>((SM::value_type)0,(SM::value_type)1)*SM::phiphWmZ;
    std::complex<SM::value_type>SM::phipchiWmZ=-SM::phimchiWpZ;
    std::complex<SM::value_type>SM::phimchiWpgamma=std::complex<SM::value_type>((SM::value_type)0,(SM::value_type)1)*SM::phiphWmgamma;
    std::complex<SM::value_type>SM::phipchiWmgamma=-SM::phimchiWpgamma;

    /* TVV couplings: */

    std::complex<SM::value_type>SM::Tgg(0,std::sqrt((SM::value_type)0.5)*SM::g_s);

    /* Standard model constructor: */

    SM::SM()
    {
	/* Lepton definitions: */
	
	(M_e>(value_type)0)?add_fermions("e-","e+",&M_e,11):add_fermions("e-","e+",11);
	(M_mu>(value_type)0)?add_fermions("mu-","mu+",&M_mu,13):add_fermions("mu-","mu+",13);
	(M_tau>(value_type)0)?add_fermions("tau-","tau+",&M_tau,15):add_fermions("tau-","tau+",15);
	add_fermions("nu_e","nu_ebar",12);
	add_fermions("nu_mu","nu_mubar",14);
	add_fermions("nu_tau","nu_taubar",16);
	
	/* Quark definitions: */
	
	(M_d>(value_type)0)?add_quarks< SU<N_c> >("d","dbar",&M_d,1):add_quarks< SU<N_c> >("d","dbar",1);
	(M_u>(value_type)0)?add_quarks< SU<N_c> >("u","ubar",&M_u,2):add_quarks< SU<N_c> >("u","ubar",2);
	(M_s>(value_type)0)?add_quarks< SU<N_c> >("s","sbar",&M_s,3):add_quarks< SU<N_c> >("s","sbar",3);
	(M_c>(value_type)0)?add_quarks< SU<N_c> >("c","cbar",&M_c,4):add_quarks< SU<N_c> >("c","cbar",4);
	(M_b>(value_type)0)?add_quarks< SU<N_c> >("b","bbar",&M_b,5):add_quarks< SU<N_c> >("b","bbar",5);
	(M_t>(value_type)0)?add_quarks< SU<N_c> >("t","tbar",&M_t,&W_t,6):add_quarks< SU<N_c> >("t","tbar",6);
	
	/* Gauge field insertions: */

	add_vector<Feynman_gauge>("gamma",22);
	add_vectors<Feynman_gauge>("W+","W-",&M_W,&W_W,24);
	add_vector<Feynman_gauge>("Z",&M_Z,&W_Z,23);
	add_gluon<Feynman_gauge,SU<N_c> >();
	
	/* Scalar particle insertions: */
	
	add_scalar("h0",&M_h0,&W_h0,25);
	add_scalar("chi",&M_Z,&W_Z);
	add_scalars("phi+","phi-",&M_W,&W_W);

	/* 4-vector vertices: */
	
	add_vertex<vvvv>("W+","W-","W+","W-",&WWWW);
	add_vertex<vvvv>("W+","Z","W-","Z",&WZWZ);
	add_vertex<vvvv>("W+","gamma","W-","Z",&WgammaWZ);
	add_vertex<vvvv>("W+","gamma","W-","gamma",&WgammaWgamma);
	
	if(four_gluon_vertex)
	{
	    add_vertex<colour_tensor::ff_contr< SU<N_c> >,vvvv>("g","g","g","g",&gggg);
	}
	else
	{
	    add_fast_4g_vertex< SU<N_c> >("g",&Tgg);
	}

	/* 3-vector vertices: */

	add_vertex<vvv>("gamma","W+","W-",&gammaee);
	add_vertex<vvv>("Z","W+","W-",&ZWW);
	add_vertex<colour_tensor::f< SU<N_c> >,vvv>("g","g","g",&ggg);

	/* 4-scalar vertices: */

	add_vertex<ssss>("h0","h0","h0","h0",&hhhh);
	add_vertex<ssss>("h0","h0","chi","chi",&hhchichi);
	add_vertex<ssss>("h0","h0","phi+","phi-",&hhphiphi);
	add_vertex<ssss>("chi","chi","chi","chi",&chichichichi);
	add_vertex<ssss>("chi","chi","phi+","phi-",&chichiphiphi);
	add_vertex<ssss>("phi+","phi-","phi+","phi-",&phiphiphiphi);

	/* 3-scalar vertices: */
	
	add_vertex<sss>("h0","h0","h0",&hhh);
	add_vertex<sss>("h0","chi","chi",&hchichi);
	add_vertex<sss>("h0","phi+","phi-",&hphiphi);

	/* Vector-fermion-fermion vertices: */

	add_vertex<vff>("gamma","e+","e-",&gammaee);
	add_vertex<vff>("gamma","mu+","mu-",&gammamumu);
	add_vertex<vff>("gamma","tau+","tau-",&gammatautau);

	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vff>("gamma","ubar","u",&gammauu);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vff>("gamma","dbar","d",&gammadd);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vff>("gamma","cbar","c",&gammacc);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vff>("gamma","sbar","s",&gammass);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vff>("gamma","tbar","t",&gammatt);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vff>("gamma","bbar","b",&gammabb);
	
	add_vertex<vffVA>("Z","e+","e-",&Zee_V,&Zee_A);
	add_vertex<vffR>("Z","nu_ebar","nu_e",&Znene);
	add_vertex<vffVA>("Z","mu+","mu-",&Zmumu_V,&Zmumu_A);
	add_vertex<vffR>("Z","nu_mubar","nu_mu",&Znmunmu);
	add_vertex<vffVA>("Z","tau+","tau-",&Ztautau_V,&Ztautau_A);
	add_vertex<vffR>("Z","nu_taubar","nu_tau",&Zntauntau);
	
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffVA>("Z","ubar","u",&Zuu_V,&Zuu_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffVA>("Z","dbar","d",&Zdd_V,&Zdd_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffVA>("Z","cbar","c",&Zcc_V,&Zcc_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffVA>("Z","sbar","s",&Zss_V,&Zss_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffVA>("Z","tbar","t",&Ztt_V,&Ztt_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffVA>("Z","bbar","b",&Zbb_V,&Zbb_A);
	
	add_vertex<vffR>("W+","nu_ebar","e-",&Wnee);
	add_vertex<vffR>("W-","e+","nu_e",&Wene);
	add_vertex<vffR>("W+","nu_mubar","mu-",&Wnmumu);
	add_vertex<vffR>("W-","mu+","nu_mu",&Wmunmu);
	add_vertex<vffR>("W+","nu_taubar","tau-",&Wntautau);
	add_vertex<vffR>("W-","tau+","nu_tau",&Wtauntau);

	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","ubar","d",&Wud);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","dbar","u",&Wdu);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","ubar","s",&Wus);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","sbar","u",&Wsu);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","ubar","b",&Wub);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","bbar","u",&Wbu);

	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","cbar","d",&Wcd);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","dbar","c",&Wdc);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","cbar","s",&Wcs);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","sbar","c",&Wsc);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","cbar","b",&Wcb);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","bbar","c",&Wbc);
	
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","tbar","d",&Wtd);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","dbar","t",&Wdt);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","tbar","s",&Wts);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","sbar","t",&Wst);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W+","tbar","b",&Wtb);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,vffR>("W-","bbar","t",&Wbt);

	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","ubar","u",&guu);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","dbar","d",&gdd);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","cbar","c",&gcc);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","sbar","s",&gss);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","tbar","t",&gtt);
	add_vertex<colour_tensor::T<fundamental_rep< SU<N_c> > >,vff>("g","bbar","b",&gbb);

	/* Scalar-fermion-fermion vertices: */
	
	add_vertex<sff>("h0","e+","e-",&hee);
	add_vertex<sff>("h0","mu+","mu-",&hmumu);
	add_vertex<sff>("h0","tau+","tau-",&htautau);
	
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff>("h0","ubar","u",&huu);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff>("h0","dbar","d",&hdd);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff>("h0","cbar","c",&hcc);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff>("h0","sbar","s",&hss);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff>("h0","tbar","t",&htt);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff>("h0","bbar","b",&hbb);
	
	add_vertex<sff5>("chi","e+","e-",&chiee);
	add_vertex<sff5>("chi","mu+","mu-",&chimumu);
	add_vertex<sff5>("chi","tau+","tau-",&chitautau);
	
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff5>("chi","ubar","u",&chiuu);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff5>("chi","dbar","d",&chidd);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff5>("chi","cbar","c",&chicc);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff5>("chi","sbar","s",&chiss);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff5>("chi","tbar","t",&chitt);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sff5>("chi","bbar","b",&chibb);
	
	add_vertex<sffL>("phi+","nu_ebar","e-",&phinee);
	add_vertex<sffR>("phi-","e+","nu_e",&phiene);
	add_vertex<sffL>("phi+","nu_mubar","mu-",&phinmumu);
	add_vertex<sffR>("phi-","mu+","nu_mu",&phimunmu);
	add_vertex<sffL>("phi+","nu_taubar","tau-",&phintautau);
	add_vertex<sffR>("phi-","tau+","nu_tau",&phitauntau);
	
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","ubar","d",&phiud_S,&phiud_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","dbar","u",&phidu_S,&phidu_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","ubar","s",&phius_S,&phius_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","sbar","u",&phisu_S,&phisu_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","ubar","b",&phiub_S,&phiub_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","bbar","u",&phibu_S,&phibu_A);

	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","cbar","d",&phicd_S,&phicd_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","dbar","c",&phidc_S,&phidc_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","cbar","s",&phics_S,&phics_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","sbar","c",&phisc_S,&phisc_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","cbar","b",&phicb_S,&phicb_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","bbar","c",&phibc_S,&phibc_A);
	
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","tbar","d",&phitd_S,&phitd_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","dbar","t",&phidt_S,&phidt_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","tbar","s",&phits_S,&phits_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","sbar","t",&phist_S,&phist_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi+","tbar","b",&phitb_S,&phitb_A);
	add_vertex<colour_tensor::d<fundamental_rep< SU<N_c> >,1,2>,sffVA>("phi-","bbar","t",&phibt_S,&phibt_A);

	/* Vector-scalar-scalar vertices: */

	add_vertex<vss>("Z","chi","h0",&Zchih);
	add_vertex<vss>("gamma","phi+","phi-",&gammaphiphi);
	add_vertex<vss>("Z","phi+","phi-",&Zphiphi);
	add_vertex<vss>("W+","phi-","h0",&Wpphimh);
	add_vertex<vss>("W-","phi+","h0",&Wmphiph);
	add_vertex<vss>("W+","phi-","chi",&Wpphimchi);
	add_vertex<vss>("W-","phi+","chi",&Wmphipchi);
	
	/* Scalar-vector-vector vertex insertions: */

	add_vertex<svv>("h0","W+","W-",&hWW);
	add_vertex<svv>("h0","Z","Z",&hZZ);
	add_vertex<svv>("phi+","W-","Z",&phipWmZ);
	add_vertex<svv>("phi-","W+","Z",&phimWpZ);
	add_vertex<svv>("phi+","W-","gamma",&phipWmgamma);
	add_vertex<svv>("phi-","W+","gamma",&phimWpgamma);

	/* Scalar-scalar-vector-vector vertex insertions: */

	add_vertex<ssvv>("h0","h0","W+","W-",&hhWW);
	add_vertex<ssvv>("h0","h0","Z","Z",&hhZZ);
	add_vertex<ssvv>("phi+","phi-","W+","W-",&phiphiWW);
	add_vertex<ssvv>("phi+","phi-","Z","Z",&phiphiZZ);
	add_vertex<ssvv>("chi","chi","W+","W-",&chichiWW);
	add_vertex<ssvv>("chi","chi","Z","Z",&chichiZZ);
	add_vertex<ssvv>("phi+","phi-","gamma","Z",&phiphigammaZ);
	add_vertex<ssvv>("phi+","phi-","gamma","gamma",&phiphigammagamma);
	add_vertex<ssvv>("phi+","h0","W-","Z",&phiphWmZ);
	add_vertex<ssvv>("phi-","h0","W+","Z",&phimhWpZ);
	add_vertex<ssvv>("phi+","h0","W-","gamma",&phiphWmgamma);
	add_vertex<ssvv>("phi-","h0","W+","gamma",&phimhWpgamma);
	add_vertex<ssvv>("phi-","chi","W+","Z",&phimchiWpZ);
	add_vertex<ssvv>("phi+","chi","W-","Z",&phipchiWmZ);
	add_vertex<ssvv>("phi-","chi","W+","gamma",&phimchiWpgamma);
	add_vertex<ssvv>("phi+","chi","W-","gamma",&phipchiWmgamma);

	decouple_zero_vertices();

	if(gauge==1)
	{
	    decouple_particle("phi+");
	    decouple_particle("chi");
	    set_propagator<unitary_gauge>("W+");
	    set_propagator<unitary_gauge>("W-");
	    set_propagator<unitary_gauge>("Z");
	    set_propagator<unitary_gauge>("gamma");
	    set_propagator<unitary_gauge>("g");
	}
	if(gauge==2)
	{
	    set_propagator<R_vector_gauge>("W+");
	    set_propagator<R_vector_gauge>("W-");
	    set_propagator<R_vector_gauge>("Z");
	    set_propagator<R_vector_gauge>("gamma");
	    set_propagator<R_vector_gauge>("g");
	    set_propagator<R_scalar_gauge>("phi+");
	    set_propagator<R_scalar_gauge>("phi-");
	    set_propagator<R_scalar_gauge>("chi");
	}

	/* Positively-charged lepton family definition: */

	construct_family("l+","e+,mu+");

	/* Positively-charged lepton family definition, including taus: */

	construct_family("L+","e+,mu+,tau+");

	/* Negatively-charged lepton family definition: */

	construct_family("l-","e-,mu-");

	/* Negatively-charged lepton family definition, including taus: */

	construct_family("L-","e-,mu-,tau-");

	/* Neutrino family definition: */

	construct_family("nu","nu_e,nu_mu,nu_tau");

	/* Antineutrino family definition: */

	construct_family("nubar","nu_ebar,nu_mubar,nu_taubar");

	/* Massive gauge boson definition: */

	construct_family("V","Z,W+,W-");

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
    }

    /* Recompute vertex couplings: */

    void SM::refresh_couplings()
    {
	refresh_weak_couplings();
	refresh_strong_couplings();
    }

    /* Recompute Yukawa vertex couplings: */

    void SM::refresh_Yukawa_couplings()
    {
	std::complex<value_type>prefactor=-(value_type)0.5*gammaee/(sin_W*cM_W);

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

    /* Recompute strong vertex couplings: */

    void SM::refresh_strong_couplings()
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

    /* Refresh CKM matrix: */

    void SM::refresh_CKM_matrix()
    {
	c12=std::sqrt(((value_type)1+s12)*((value_type)1-s12));
	c13=std::sqrt(((value_type)1+s13)*((value_type)1-s13));
	c23=std::sqrt(((value_type)1+s23)*((value_type)1-s23));
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

    /* Refresh particle masses (to use when along the run massless fermions are
     * to be made massive): */

    void SM::refresh_fermion_masses()
    {
	if(initialised())
	{
	    (M_e>(value_type)0)?set_mass("e-",&M_e):set_massless("e-");
	    (M_mu>(value_type)0)?set_mass("mu-",&M_mu):set_massless("mu-");
	    (M_tau>(value_type)0)?set_mass("tau-",&M_tau):set_massless("tau-");
	    (M_u>(value_type)0)?set_mass("u",&M_u):set_massless("u");
	    (M_d>(value_type)0)?set_mass("d",&M_d):set_massless("d");
	    (M_c>(value_type)0)?set_mass("c",&M_c):set_massless("c");
	    (M_s>(value_type)0)?set_mass("s",&M_s):set_massless("s");
	}
    }

    /* Recompute weak couplings: */

    void SM::refresh_weak_couplings()
    {
	Q_e=std::sqrt((value_type)4*pi*SM::alpha);
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

    /* Recompute particle decay widths: */

    void SM::refresh_widths()
    {
	refresh_W_width();
	refresh_Z_width();
	refresh_Higgs_width();
	refresh_top_width();
    }

    /* Recompute W-boson decay width: */

    void SM::refresh_W_width()
    {
	W_W=value_type(3+2*N_c)*G_F*std::pow(M_W,3)/((value_type)6*pi*std::sqrt((value_type)2));
    }

    /* Recompute Z-boson decay width: */

    void SM::refresh_Z_width()
    {
	value_type cw2=(M_W*M_W)/(M_Z*M_Z);
	value_type factor=G_F*std::pow(M_Z,3)/(108*std::sqrt((value_type)2)*pi);
	W_Z=(value_type(81+49*N_c)-(value_type)(162+92*N_c)*cw2+(value_type)(108+88*N_c)*cw2*cw2)*factor;
    }

    /* Recompute Higgs decay width: */

    void SM::refresh_Higgs_width()
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
	prefactor*=(N_c);
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

    /* Recompute top quark decay width: */

    void SM::refresh_top_width()
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

    /* Function resetting all decay widths to zero: */

    void SM::discard_widths()
    {
	W_Z=(value_type)0;
	W_W=(value_type)0;
	W_h0=(value_type)0;
	W_t=(value_type)0;
	refresh_weak_couplings();
    }

    /* Higgs mass value input: */

    void SM::set_Higgs_mass(const value_type& m)
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

    /* Fine structure input: */

    void SM::set_alpha(const value_type& a)
    {
	alpha=a;
	refresh_weak_couplings();
    }

    /* Fine structure input: */

    void SM::set_alpha_s(const value_type& a)
    {
	alpha_s=a;
	refresh_strong_couplings();
    }

    /* Strong scale value input: */

    void SM::set_QCD_scale(const value_type& mu)
    {
	QCD_scale=mu;
	value_type s=mu/(value_type)SM_params::QCD_scale;
	alpha_s=(value_type)SM_params::alpha_s/((value_type)1+(value_type)SM_params::alpha_s*(value_type)(11*N_c-12)*std::log(s)/((value_type)6*pi));
	refresh_strong_couplings();
    }

    /* Reset all parameters using the default input parameters in the
     * SM_params.h file: */

    void SM::set_default_params()
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
	delta=(SM::value_type)delta;
	
	refresh_couplings();
	set_Feynman_gauge();
	set_auxiliary_QCD_field();
    }

    /* Function setting all vectors to the unitary gauge: */

    void SM::set_unitary_gauge()
    {
	if(initialised() and gauge!=1)
	{
	    decouple_particle("phi+");
	    decouple_particle("chi");
	    set_propagator<unitary_gauge>("W+");
	    set_propagator<unitary_gauge>("W-");
	    set_propagator<unitary_gauge>("Z");
	    set_propagator<unitary_gauge>("gamma");
	    set_propagator<unitary_gauge>("g");
	}
	gauge=1;
    }

    /* Function setting all vectors to the Feynman gauge and add would-be
     * Goldstone bosons if necessary: */

    void SM::set_Feynman_gauge()
    {
	if(initialised() and gauge!=0)
	{
	    if(gauge==1)
	    {
		couple_particle("chi");
		couple_particle("phi+");
		decouple_zero_vertices();
	    }
	    set_propagator<Feynman_gauge>("W+");
	    set_propagator<Feynman_gauge>("W-");
	    set_propagator<Feynman_gauge>("Z");
	    set_propagator<Feynman_gauge>("gamma");
	    set_propagator<Feynman_gauge>("g");
	    set_propagator<scalar_propagator>("chi");
	    set_propagator<scalar_propagator>("phi+");
	    set_propagator<scalar_propagator>("phi-");
	}
	gauge=0;
    }
    
    /* Function setting all vectors to the R-xi gauge and add would-be
     * Goldstone bosons if necessary: */

    void SM::set_R_xi_gauge()
    {
	if(initialised() and gauge!=2)
	{
	    if(gauge==1)
	    {
		couple_particle("chi");
		couple_particle("phi+");
		decouple_zero_vertices();
	    }
	    set_propagator<R_vector_gauge>("W+");
	    set_propagator<R_vector_gauge>("W-");
	    set_propagator<R_vector_gauge>("Z");
	    set_propagator<R_vector_gauge>("gamma");
	    set_propagator<R_vector_gauge>("g");
	    set_propagator<R_scalar_gauge>("phi+");
	    set_propagator<R_scalar_gauge>("phi-");
	    set_propagator<R_scalar_gauge>("chi");
	}
	gauge=2;
    }

    /* Access the xi-value in the R-xi propagators: */

    void SM::set_xi(const value_type& x)
    {
	R_gauge<SM>::xi=x;
    }

    /* Function returning the number of QCD colours: */

    std::size_t SM::QCD_colours()
    {
	return N_c;
    }

    /* Function switching to a 4-gluon vertex Lagrangian: */

    void SM::set_4_gluon_vertex()
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

    void SM::set_auxiliary_QCD_field()
    {
	if(initialised() and four_gluon_vertex)
	{
	    erase_vertex("g","g","g","g");
	    add_fast_4g_vertex< SU<N_c> >("g",&Tgg);
	}
	four_gluon_vertex=false;
    }
}

