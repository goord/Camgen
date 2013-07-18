//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/EWSM.h>
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

namespace Camgen
{
    /* Compile-time constant definitions: */

    const std::size_t EWSM::dimension;
    const bool EWSM::coloured;
    const bool EWSM::continuous_helicities;
    const int EWSM::beam_direction;
    const std::size_t EWSM::N_c;

    const EWSM::value_type EWSM::pi=std::acos(-(EWSM::value_type)1);
    
    /* Real input parameters: */

    EWSM::value_type EWSM::M_e=(EWSM::value_type)0;
    EWSM::value_type EWSM::M_mu=(EWSM::value_type)0;
    EWSM::value_type EWSM::M_tau=(EWSM::value_type)SM_params::M_tau;
    EWSM::value_type EWSM::M_u=(EWSM::value_type)0;
    EWSM::value_type EWSM::M_d=(EWSM::value_type)0;
    EWSM::value_type EWSM::M_c=(EWSM::value_type)SM_params::M_c;
    EWSM::value_type EWSM::M_s=(EWSM::value_type)0;
    EWSM::value_type EWSM::M_t=(EWSM::value_type)SM_params::M_t;
    EWSM::value_type EWSM::M_b=(EWSM::value_type)SM_params::M_b;
    EWSM::value_type EWSM::M_Z=(EWSM::value_type)SM_params::M_Z;
    EWSM::value_type EWSM::M_h0=(EWSM::value_type)120;
    EWSM::value_type EWSM::alpha=(EWSM::value_type)SM_params::alpha;
    EWSM::value_type EWSM::G_F=(EWSM::value_type)SM_params::G_F;
    const EWSM::value_type EWSM::QCD_scale=-(EWSM::value_type)1;
    const EWSM::value_type EWSM::alpha_s=-(EWSM::value_type)1;

    EWSM::value_type EWSM::s12=(EWSM::value_type)SM_params::s12;
    EWSM::value_type EWSM::s23=(EWSM::value_type)SM_params::s23;
    EWSM::value_type EWSM::s13=(EWSM::value_type)SM_params::s13;
    EWSM::value_type EWSM::delta=(EWSM::value_type)0;

    /* Computed parameters: */
    
    EWSM::value_type EWSM::Q_e = std::sqrt((EWSM::value_type)4*EWSM::pi*EWSM::alpha);
    EWSM::value_type EWSM::M_W = EWSM::M_Z*std::sqrt((EWSM::value_type)0.5+(EWSM::value_type)0.5*std::sqrt((EWSM::value_type)1-(EWSM::value_type)4*EWSM::pi*EWSM::alpha/(std::sqrt((EWSM::value_type)2)*EWSM::G_F*EWSM::M_Z*EWSM::M_Z)));
    EWSM::value_type EWSM::W_W = (3+2*N_c)*EWSM::G_F*std::pow(EWSM::M_W,3)/(EWSM::pi*std::sqrt((EWSM::value_type)2))/(EWSM::value_type)6;
    EWSM::value_type EWSM::W_Z = ((EWSM::value_type)(162+49*N_c)-(EWSM::value_type)(324+92*N_c)*std::pow(EWSM::M_W/EWSM::M_Z,2)+(EWSM::value_type)(216+88*N_c)*std::pow(EWSM::M_W/EWSM::M_Z,4))*EWSM::G_F*std::pow(EWSM::M_Z,3)/((EWSM::value_type)108*EWSM::pi*std::sqrt((EWSM::value_type)2));
    EWSM::value_type EWSM::W_h0 = N_c*EWSM::G_F*EWSM::M_b*EWSM::M_b*EWSM::M_h0*std::pow((EWSM::value_type)1-std::pow((EWSM::value_type)2*EWSM::M_b/EWSM::M_h0,2),1.5)/(std::sqrt((EWSM::value_type)2)*(EWSM::value_type)4*EWSM::pi);
    EWSM::value_type EWSM::W_t = EWSM::G_F*std::pow(EWSM::M_t,3)*((EWSM::value_type)1-std::pow(EWSM::M_W/EWSM::M_t,2))*((EWSM::value_type)1+std::pow(EWSM::M_W/EWSM::M_t,2)-2*std::pow(EWSM::M_W/EWSM::M_t,4))/(std::sqrt((EWSM::value_type)128)*EWSM::pi);

    std::complex<EWSM::value_type>EWSM::cM_Z=std::sqrt(std::complex<EWSM::value_type>(EWSM::M_Z*EWSM::M_Z,-EWSM::M_Z*EWSM::W_Z));
    std::complex<EWSM::value_type>EWSM::cM_W=std::sqrt(std::complex<EWSM::value_type>(EWSM::M_W*EWSM::M_W,-EWSM::M_W*EWSM::W_W));
    std::complex<EWSM::value_type>EWSM::cos_W=EWSM::cM_W/EWSM::cM_Z;
    std::complex<EWSM::value_type>EWSM::sin_W=std::sqrt(std::complex<EWSM::value_type>(1,0)-EWSM::cos_W*EWSM::cos_W);
    std::complex<EWSM::value_type>EWSM::tan_W=EWSM::sin_W/EWSM::cos_W;
    std::complex<EWSM::value_type>EWSM::cos2_W=EWSM::cos_W*EWSM::cos_W-EWSM::sin_W*EWSM::sin_W;
    std::complex<EWSM::value_type>EWSM::sin2_W=(EWSM::value_type)2*EWSM::sin_W*EWSM::cos_W;
    
    EWSM::value_type EWSM::c12=std::sqrt(((EWSM::value_type)1+EWSM::s12)*((EWSM::value_type)1-EWSM::s12));
    EWSM::value_type EWSM::c23=std::sqrt(((EWSM::value_type)1+EWSM::s23)*((EWSM::value_type)1-EWSM::s23));
    EWSM::value_type EWSM::c13=std::sqrt(((EWSM::value_type)1+EWSM::s13)*((EWSM::value_type)1-EWSM::s13));
    
    std::complex<EWSM::value_type>EWSM::V_CKM[3][3]={{EWSM::c12*EWSM::c13,			   
						      EWSM::s12*EWSM::c13,			    
						      EWSM::s13},
						    {-EWSM::s12*EWSM::c23-EWSM::c12*EWSM::s23*EWSM::s13,
						      EWSM::c12*EWSM::c23-EWSM::s12*EWSM::s23*EWSM::s13,
						      EWSM::s23*EWSM::c13},
						     {EWSM::s12*EWSM::s23-EWSM::c12*EWSM::c23*EWSM::s13,
					             -EWSM::c12*EWSM::s23-EWSM::s12*EWSM::c23*EWSM::s13,
						      EWSM::c23*EWSM::c13}};
    
    /* Gauge switch initialisation: */

    int EWSM::gauge=0;

    /* VVVV couplings: */

    std::complex<EWSM::value_type>EWSM::WWWW=std::complex<EWSM::value_type>(0,EWSM::Q_e*EWSM::Q_e)/(EWSM::sin_W*EWSM::sin_W);
    std::complex<EWSM::value_type>EWSM::WZWZ=-EWSM::WWWW*EWSM::cos_W*EWSM::cos_W;
    std::complex<EWSM::value_type>EWSM::WgammaWZ=EWSM::WWWW*EWSM::cos_W*EWSM::sin_W;
    std::complex<EWSM::value_type>EWSM::WgammaWgamma=-EWSM::WWWW*EWSM::sin_W*sin_W;

    /* VVV couplings: */

    std::complex<EWSM::value_type>EWSM::gammaWW(0,EWSM::Q_e);
    std::complex<EWSM::value_type>EWSM::ZWW=-EWSM::gammaWW/EWSM::tan_W;

    /* SSSS couplings: */

    std::complex<EWSM::value_type>EWSM::hhhh=-(EWSM::value_type)0.75*std::complex<EWSM::value_type>(0,EWSM::Q_e*EWSM::Q_e)*std::pow(EWSM::M_h0/(EWSM::sin_W*EWSM::cM_W),2);
    std::complex<EWSM::value_type>EWSM::hhchichi=EWSM::hhhh/(EWSM::value_type)3;
    std::complex<EWSM::value_type>EWSM::hhphiphi=EWSM::hhchichi;
    std::complex<EWSM::value_type>EWSM::chichichichi=EWSM::hhhh;
    std::complex<EWSM::value_type>EWSM::chichiphiphi=EWSM::hhphiphi;
    std::complex<EWSM::value_type>EWSM::phiphiphiphi=(EWSM::value_type)2*EWSM::hhchichi;

    /* SSS couplings: */

    std::complex<EWSM::value_type>EWSM::hhh=-(EWSM::value_type)1.5*std::complex<EWSM::value_type>(0,EWSM::Q_e)*EWSM::M_h0*EWSM::M_h0/(EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::hchichi=EWSM::hhh/(EWSM::value_type)3;
    std::complex<EWSM::value_type>EWSM::hphiphi=EWSM::hchichi;
    
    /* VFF couplings: */

    std::complex<EWSM::value_type>EWSM::gammaee(0,EWSM::Q_e);
    std::complex<EWSM::value_type>EWSM::gammamumu=EWSM::gammaee;
    std::complex<EWSM::value_type>EWSM::gammatautau=EWSM::gammaee;
    std::complex<EWSM::value_type>EWSM::gammauu(0,-(EWSM::value_type)2*EWSM::Q_e/(EWSM::value_type)3);
    std::complex<EWSM::value_type>EWSM::gammadd(0,EWSM::Q_e/(EWSM::value_type)3);
    std::complex<EWSM::value_type>EWSM::gammacc=EWSM::gammauu;
    std::complex<EWSM::value_type>EWSM::gammass=EWSM::gammadd;
    std::complex<EWSM::value_type>EWSM::gammatt=EWSM::gammauu;
    std::complex<EWSM::value_type>EWSM::gammabb=EWSM::gammadd;
    
    std::complex<EWSM::value_type>EWSM::Znene=EWSM::gammaee/EWSM::sin2_W;
    std::complex<EWSM::value_type>EWSM::Zee_V=EWSM::Znene*((EWSM::value_type)2*EWSM::sin_W*EWSM::sin_W-(EWSM::value_type)0.5);
    std::complex<EWSM::value_type>EWSM::Zee_A=(EWSM::value_type)0.5*EWSM::Znene;
    std::complex<EWSM::value_type>EWSM::Znmunmu=EWSM::Znene;
    std::complex<EWSM::value_type>EWSM::Zmumu_V=EWSM::Zee_V;
    std::complex<EWSM::value_type>EWSM::Zmumu_A=EWSM::Zee_A;
    std::complex<EWSM::value_type>EWSM::Zntauntau=EWSM::Znene;
    std::complex<EWSM::value_type>EWSM::Ztautau_V=EWSM::Zee_V;
    std::complex<EWSM::value_type>EWSM::Ztautau_A=EWSM::Zee_A;
    std::complex<EWSM::value_type>EWSM::Zuu_V=EWSM::Znene*((EWSM::value_type)0.5-(EWSM::value_type)4*EWSM::sin_W*EWSM::sin_W/(EWSM::value_type)3);
    std::complex<EWSM::value_type>EWSM::Zuu_A=-EWSM::Zee_A;
    std::complex<EWSM::value_type>EWSM::Zdd_V=EWSM::Znene*(-(EWSM::value_type)0.5+(EWSM::value_type)2*EWSM::sin_W*EWSM::sin_W/(EWSM::value_type)3);
    std::complex<EWSM::value_type>EWSM::Zdd_A=EWSM::Zee_A;
    std::complex<EWSM::value_type>EWSM::Zcc_V=EWSM::Zuu_V;
    std::complex<EWSM::value_type>EWSM::Zcc_A=EWSM::Zuu_A;
    std::complex<EWSM::value_type>EWSM::Zss_V=EWSM::Zdd_V;
    std::complex<EWSM::value_type>EWSM::Zss_A=EWSM::Zdd_A;
    std::complex<EWSM::value_type>EWSM::Ztt_V=EWSM::Zuu_V;
    std::complex<EWSM::value_type>EWSM::Ztt_A=EWSM::Zuu_A;
    std::complex<EWSM::value_type>EWSM::Zbb_V=EWSM::Zdd_V;
    std::complex<EWSM::value_type>EWSM::Zbb_A=EWSM::Zdd_A;

    std::complex<EWSM::value_type>EWSM::Wnee=EWSM::gammaee/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W);
    std::complex<EWSM::value_type>EWSM::Wnmumu=EWSM::Wnee;
    std::complex<EWSM::value_type>EWSM::Wntautau=EWSM::Wnee;
    std::complex<EWSM::value_type>EWSM::Wene=EWSM::Wnee;
    std::complex<EWSM::value_type>EWSM::Wmunmu=EWSM::Wnmumu;
    std::complex<EWSM::value_type>EWSM::Wtauntau=EWSM::Wntautau;
    std::complex<EWSM::value_type>EWSM::Wud=EWSM::Wnee*EWSM::V_CKM[0][0];
    std::complex<EWSM::value_type>EWSM::Wdu=EWSM::Wud;
    std::complex<EWSM::value_type>EWSM::Wus=EWSM::Wnee*EWSM::V_CKM[0][1];
    std::complex<EWSM::value_type>EWSM::Wsu=EWSM::Wus;
    std::complex<EWSM::value_type>EWSM::Wub=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::Wbu=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::Wcd=EWSM::Wnee*EWSM::V_CKM[1][0];
    std::complex<EWSM::value_type>EWSM::Wdc=EWSM::Wcd;
    std::complex<EWSM::value_type>EWSM::Wcs=EWSM::Wnee*EWSM::V_CKM[1][1];
    std::complex<EWSM::value_type>EWSM::Wsc=EWSM::Wcs;
    std::complex<EWSM::value_type>EWSM::Wcb=EWSM::Wnee*EWSM::V_CKM[1][2];
    std::complex<EWSM::value_type>EWSM::Wbc=EWSM::Wcb;
    std::complex<EWSM::value_type>EWSM::Wtd=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::Wdt=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::Wts=EWSM::Wnee*EWSM::V_CKM[2][1];
    std::complex<EWSM::value_type>EWSM::Wst=EWSM::Wts;
    std::complex<EWSM::value_type>EWSM::Wtb=EWSM::Wnee*EWSM::V_CKM[2][2];
    std::complex<EWSM::value_type>EWSM::Wbt=EWSM::Wtb;

    /* SFF couplings: */

    std::complex<EWSM::value_type>EWSM::hee=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::hmumu=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::htautau=-(EWSM::value_type)0.5*EWSM::gammaee*EWSM::M_tau/(EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::huu=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::hdd=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::hcc=-(EWSM::value_type)0.5*EWSM::gammaee*EWSM::M_c/(EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::hss=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::htt=-(EWSM::value_type)0.5*EWSM::gammaee*EWSM::M_t/(EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::hbb=-(EWSM::value_type)0.5*EWSM::gammaee*EWSM::M_b/(EWSM::sin_W*EWSM::cM_W);

    std::complex<EWSM::value_type>EWSM::chiee=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::chimumu=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::chitautau=(EWSM::value_type)0.5*EWSM::Q_e*EWSM::M_tau/(EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::chiuu=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::chidd=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::chicc=-(EWSM::value_type)0.5*EWSM::Q_e*EWSM::M_t/(EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::chiss=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::chitt=-(EWSM::value_type)0.5*EWSM::Q_e*EWSM::M_t/(EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::chibb=(EWSM::value_type)0.5*EWSM::Q_e*EWSM::M_b/(EWSM::sin_W*EWSM::cM_W);

    std::complex<EWSM::value_type>EWSM::phinee=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phiene=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phinmumu=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phimunmu=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phintautau=-EWSM::gammaee*EWSM::M_tau/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phitauntau=EWSM::phintautau;
    
    std::complex<EWSM::value_type>EWSM::phiud_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phiud_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phidu_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phidu_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phius_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phius_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phisu_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phisu_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phiub_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phiub_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phibu_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phibu_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phicd_S=(EWSM::value_type)0.5*EWSM::gammaee*EWSM::M_c*EWSM::V_CKM[1][0]/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phicd_A=-EWSM::phicd_S;
    std::complex<EWSM::value_type>EWSM::phidc_S=EWSM::phicd_S;
    std::complex<EWSM::value_type>EWSM::phidc_A=-EWSM::phicd_A;
    std::complex<EWSM::value_type>EWSM::phics_S=(EWSM::value_type)0.5*EWSM::gammaee*EWSM::M_c*EWSM::V_CKM[1][1]/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phics_A=-EWSM::phics_S;
    std::complex<EWSM::value_type>EWSM::phisc_S=EWSM::phics_S;
    std::complex<EWSM::value_type>EWSM::phisc_A=-EWSM::phics_A;
    std::complex<EWSM::value_type>EWSM::phicb_S=(EWSM::value_type)0.5*EWSM::gammaee*(EWSM::M_c-EWSM::M_b)*EWSM::V_CKM[1][2]/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phicb_A=-(EWSM::value_type)0.5*EWSM::gammaee*(EWSM::M_c+EWSM::M_b)*EWSM::V_CKM[1][2]/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phibc_S=EWSM::phicb_S;
    std::complex<EWSM::value_type>EWSM::phibc_A=-EWSM::phicb_A;
    std::complex<EWSM::value_type>EWSM::phitd_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phitd_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phidt_S=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phidt_A=(EWSM::value_type)0;
    std::complex<EWSM::value_type>EWSM::phits_S=(EWSM::value_type)0.5*EWSM::gammaee*EWSM::M_t*EWSM::V_CKM[2][1]/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phits_A=-EWSM::phits_S;
    std::complex<EWSM::value_type>EWSM::phist_S=EWSM::phits_S;
    std::complex<EWSM::value_type>EWSM::phist_A=-EWSM::phits_A;
    std::complex<EWSM::value_type>EWSM::phitb_S=(EWSM::value_type)0.5*EWSM::gammaee*(EWSM::M_t-EWSM::M_b)/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phitb_A=-(EWSM::value_type)0.5*EWSM::gammaee*(EWSM::M_t+EWSM::M_b)/(std::sqrt((EWSM::value_type)2)*EWSM::sin_W*EWSM::cM_W);
    std::complex<EWSM::value_type>EWSM::phibt_S=EWSM::phitb_S;
    std::complex<EWSM::value_type>EWSM::phibt_A=-EWSM::phitb_A;

    /* VSS couplings: */

    std::complex<EWSM::value_type>EWSM::Zchih=EWSM::Q_e/EWSM::sin2_W;
    std::complex<EWSM::value_type>EWSM::gammaphiphi=-EWSM::gammaee;
    std::complex<EWSM::value_type>EWSM::Zphiphi=EWSM::gammaee*EWSM::cos2_W/EWSM::sin2_W;
    std::complex<EWSM::value_type>EWSM::Wpphimh=-(EWSM::value_type)0.5*EWSM::gammaee/EWSM::sin_W;
    std::complex<EWSM::value_type>EWSM::Wmphiph=-EWSM::Wpphimh;
    std::complex<EWSM::value_type>EWSM::Wpphimchi=(EWSM::value_type)0.5*EWSM::Q_e/EWSM::sin_W;
    std::complex<EWSM::value_type>EWSM::Wmphipchi=EWSM::Wpphimchi;
    
    /* SVV couplings: */

    std::complex<EWSM::value_type>EWSM::hWW=EWSM::gammaee*EWSM::cM_W/EWSM::sin_W;
    std::complex<EWSM::value_type>EWSM::hZZ=(EWSM::value_type)2*EWSM::gammaee*EWSM::cM_Z/EWSM::sin2_W;
    std::complex<EWSM::value_type>EWSM::phipWmZ=-EWSM::gammaee*EWSM::cM_W*EWSM::tan_W;
    std::complex<EWSM::value_type>EWSM::phimWpZ=EWSM::phipWmZ;
    std::complex<EWSM::value_type>EWSM::phipWmgamma=-EWSM::gammaee*EWSM::cM_W;
    std::complex<EWSM::value_type>EWSM::phimWpgamma=EWSM::phipWmgamma;

    /* SSVV couplings: */

    std::complex<EWSM::value_type>EWSM::hhWW=std::complex<EWSM::value_type>((EWSM::value_type)0,(EWSM::value_type)0.5)*std::pow(EWSM::Q_e/EWSM::sin_W,2);
    std::complex<EWSM::value_type>EWSM::chichiWW=EWSM::hhWW;
    std::complex<EWSM::value_type>EWSM::phiphiWW=EWSM::hhWW;
    std::complex<EWSM::value_type>EWSM::hhZZ=EWSM::hhWW/(EWSM::cos_W*EWSM::cos_W);
    std::complex<EWSM::value_type>EWSM::chichiZZ=EWSM::hhZZ;
    std::complex<EWSM::value_type>EWSM::phiphiZZ=EWSM::hhZZ*EWSM::cos2_W*EWSM::cos2_W;
    std::complex<EWSM::value_type>EWSM::phiphigammaZ=-(EWSM::value_type)2*EWSM::hhWW*EWSM::cos2_W*EWSM::tan_W;
    std::complex<EWSM::value_type>EWSM::phiphigammagamma=(EWSM::value_type)2*EWSM::Q_e*EWSM::gammaee;
    std::complex<EWSM::value_type>EWSM::phiphWmZ=-(EWSM::value_type)0.5*EWSM::Q_e*EWSM::gammaee/EWSM::cos_W;
    std::complex<EWSM::value_type>EWSM::phimhWpZ=EWSM::phiphWmZ;
    std::complex<EWSM::value_type>EWSM::phiphWmgamma=EWSM::phimhWpZ/EWSM::tan_W;
    std::complex<EWSM::value_type>EWSM::phimhWpgamma=EWSM::phiphWmgamma;
    std::complex<EWSM::value_type>EWSM::phimchiWpZ=std::complex<EWSM::value_type>((EWSM::value_type)0,(EWSM::value_type)1)*EWSM::phiphWmZ;
    std::complex<EWSM::value_type>EWSM::phipchiWmZ=-EWSM::phimchiWpZ;
    std::complex<EWSM::value_type>EWSM::phimchiWpgamma=std::complex<EWSM::value_type>((EWSM::value_type)0,(EWSM::value_type)1)*EWSM::phiphWmgamma;
    std::complex<EWSM::value_type>EWSM::phipchiWmgamma=-EWSM::phimchiWpgamma;

    /* Standard model constructor: */

    EWSM::EWSM()
    {
	/* Lepton definitions: */
	
	(M_e>(value_type)0)?add_fermions("e-","e+",&M_e,11):add_fermions("e-","e+",11);
	(M_mu>(value_type)0)?add_fermions("mu-","mu+",&M_mu,13):add_fermions("mu-","mu+",13);
	(M_tau>(value_type)0)?add_fermions("tau-","tau+",&M_tau,15):add_fermions("tau-","tau+",15);
	add_fermions("nu_e","nu_ebar",12);
	add_fermions("nu_mu","nu_mubar",14);
	add_fermions("nu_tau","nu_taubar",16);
	
	/* Quark definitions: */
	
	(M_d>(value_type)0)?add_fermions("d","dbar",&M_d,1):add_fermions("d","dbar",1);
	(M_u>(value_type)0)?add_fermions("u","ubar",&M_u,2):add_fermions("u","ubar",2);
	(M_s>(value_type)0)?add_fermions("s","sbar",&M_s,3):add_fermions("s","sbar",3);
	(M_c>(value_type)0)?add_fermions("c","cbar",&M_c,4):add_fermions("c","cbar",4);
	(M_b>(value_type)0)?add_fermions("b","bbar",&M_b,5):add_fermions("b","bbar",5);
	(M_t>(value_type)0)?add_fermions("t","tbar",&M_t,&W_t,6):add_fermions("t","tbar",6);
	
	/* Gauge field insertions: */

	add_vector<Feynman_gauge>("gamma",22);
	add_vectors<Feynman_gauge>("W+","W-",&M_W,&W_W,24);
	add_vector<Feynman_gauge>("Z",&M_Z,&W_Z,23);
	
	/* Scalar particle insertions: */
	
	add_scalar("h0",&M_h0,&W_h0,25);
	add_scalar("chi",&M_Z,&W_Z);
	add_scalars("phi+","phi-",&M_W,&W_W);

	/* 4-vector vertices: */
	
	add_vertex<vvvv>("W+","W-","W+","W-",&WWWW);
	add_vertex<vvvv>("W+","Z","W-","Z",&WZWZ);
	add_vertex<vvvv>("W+","gamma","W-","Z",&WgammaWZ);
	add_vertex<vvvv>("W+","gamma","W-","gamma",&WgammaWgamma);

	/* 3-vector vertices: */

	add_vertex<vvv>("gamma","W+","W-",&gammaee);
	add_vertex<vvv>("Z","W+","W-",&ZWW);

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

	add_vertex<vff>("gamma","ubar","u",&gammauu);
	add_vertex<vff>("gamma","dbar","d",&gammadd);
	add_vertex<vff>("gamma","cbar","c",&gammacc);
	add_vertex<vff>("gamma","sbar","s",&gammass);
	add_vertex<vff>("gamma","tbar","t",&gammatt);
	add_vertex<vff>("gamma","bbar","b",&gammabb);
	
	add_vertex<vffVA>("Z","e+","e-",&Zee_V,&Zee_A);
	add_vertex<vffR>("Z","nu_ebar","nu_e",&Znene);
	add_vertex<vffVA>("Z","mu+","mu-",&Zmumu_V,&Zmumu_A);
	add_vertex<vffR>("Z","nu_mubar","nu_mu",&Znmunmu);
	add_vertex<vffVA>("Z","tau+","tau-",&Ztautau_V,&Ztautau_A);
	add_vertex<vffR>("Z","nu_taubar","nu_tau",&Zntauntau);
	
	add_vertex<vffVA>("Z","ubar","u",&Zuu_V,&Zuu_A);
	add_vertex<vffVA>("Z","dbar","d",&Zdd_V,&Zdd_A);
	add_vertex<vffVA>("Z","cbar","c",&Zcc_V,&Zcc_A);
	add_vertex<vffVA>("Z","sbar","s",&Zss_V,&Zss_A);
	add_vertex<vffVA>("Z","tbar","t",&Ztt_V,&Ztt_A);
	add_vertex<vffVA>("Z","bbar","b",&Zbb_V,&Zbb_A);
	
	add_vertex<vffR>("W+","nu_ebar","e-",&Wnee);
	add_vertex<vffR>("W-","e+","nu_e",&Wene);
	add_vertex<vffR>("W+","nu_mubar","mu-",&Wnmumu);
	add_vertex<vffR>("W-","mu+","nu_mu",&Wmunmu);
	add_vertex<vffR>("W+","nu_taubar","tau-",&Wntautau);
	add_vertex<vffR>("W-","tau+","nu_tau",&Wtauntau);

	add_vertex<vffR>("W+","ubar","d",&Wud);
	add_vertex<vffR>("W-","dbar","u",&Wdu);
	add_vertex<vffR>("W+","ubar","s",&Wus);
	add_vertex<vffR>("W-","sbar","u",&Wsu);
	add_vertex<vffR>("W+","ubar","b",&Wub);
	add_vertex<vffR>("W-","bbar","u",&Wbu);

	add_vertex<vffR>("W+","cbar","d",&Wcd);
	add_vertex<vffR>("W-","dbar","c",&Wdc);
	add_vertex<vffR>("W+","cbar","s",&Wcs);
	add_vertex<vffR>("W-","sbar","c",&Wsc);
	add_vertex<vffR>("W+","cbar","b",&Wcb);
	add_vertex<vffR>("W-","bbar","c",&Wbc);
	
	add_vertex<vffR>("W+","tbar","d",&Wtd);
	add_vertex<vffR>("W-","dbar","t",&Wdt);
	add_vertex<vffR>("W+","tbar","s",&Wts);
	add_vertex<vffR>("W-","sbar","t",&Wst);
	add_vertex<vffR>("W+","tbar","b",&Wtb);
	add_vertex<vffR>("W-","bbar","t",&Wbt);

	/* Scalar-fermion-fermion vertices: */
	
	add_vertex<sff>("h0","e+","e-",&hee);
	add_vertex<sff>("h0","mu+","mu-",&hmumu);
	add_vertex<sff>("h0","tau+","tau-",&htautau);
	
	add_vertex<sff>("h0","ubar","u",&huu);
	add_vertex<sff>("h0","dbar","d",&hdd);
	add_vertex<sff>("h0","cbar","c",&hcc);
	add_vertex<sff>("h0","sbar","s",&hss);
	add_vertex<sff>("h0","tbar","t",&htt);
	add_vertex<sff>("h0","bbar","b",&hbb);
	
	add_vertex<sff5>("chi","e+","e-",&chiee);
	add_vertex<sff5>("chi","mu+","mu-",&chimumu);
	add_vertex<sff5>("chi","tau+","tau-",&chitautau);
	
	add_vertex<sff5>("chi","ubar","u",&chiuu);
	add_vertex<sff5>("chi","dbar","d",&chidd);
	add_vertex<sff5>("chi","cbar","c",&chicc);
	add_vertex<sff5>("chi","sbar","s",&chiss);
	add_vertex<sff5>("chi","tbar","t",&chitt);
	add_vertex<sff5>("chi","bbar","b",&chibb);
	
	add_vertex<sffL>("phi+","nu_ebar","e-",&phinee);
	add_vertex<sffR>("phi-","e+","nu_e",&phiene);
	add_vertex<sffL>("phi+","nu_mubar","mu-",&phinmumu);
	add_vertex<sffR>("phi-","mu+","nu_mu",&phimunmu);
	add_vertex<sffL>("phi+","nu_taubar","tau-",&phintautau);
	add_vertex<sffR>("phi-","tau+","nu_tau",&phitauntau);
	
	add_vertex<sffVA>("phi+","ubar","d",&phiud_S,&phiud_A);
	add_vertex<sffVA>("phi-","dbar","u",&phidu_S,&phidu_A);
	add_vertex<sffVA>("phi+","ubar","s",&phius_S,&phius_A);
	add_vertex<sffVA>("phi-","sbar","u",&phisu_S,&phisu_A);
	add_vertex<sffVA>("phi+","ubar","b",&phiub_S,&phiub_A);
	add_vertex<sffVA>("phi-","bbar","u",&phibu_S,&phibu_A);

	add_vertex<sffVA>("phi+","cbar","d",&phicd_S,&phicd_A);
	add_vertex<sffVA>("phi-","dbar","c",&phidc_S,&phidc_A);
	add_vertex<sffVA>("phi+","cbar","s",&phics_S,&phics_A);
	add_vertex<sffVA>("phi-","sbar","c",&phisc_S,&phisc_A);
	add_vertex<sffVA>("phi+","cbar","b",&phicb_S,&phicb_A);
	add_vertex<sffVA>("phi-","bbar","c",&phibc_S,&phibc_A);
	
	add_vertex<sffVA>("phi+","tbar","d",&phitd_S,&phitd_A);
	add_vertex<sffVA>("phi-","dbar","t",&phidt_S,&phidt_A);
	add_vertex<sffVA>("phi+","tbar","s",&phits_S,&phits_A);
	add_vertex<sffVA>("phi-","sbar","t",&phist_S,&phist_A);
	add_vertex<sffVA>("phi+","tbar","b",&phitb_S,&phitb_A);
	add_vertex<sffVA>("phi-","bbar","t",&phibt_S,&phibt_A);

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
	}
	if(gauge==2)
	{
	    set_propagator<R_vector_gauge>("W+");
	    set_propagator<R_vector_gauge>("W-");
	    set_propagator<R_vector_gauge>("Z");
	    set_propagator<R_vector_gauge>("gamma");
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

	construct_family("p","u,ubar,d,dbar,c,cbar,s,sbar");

	/* Proton constituent parton family definition, including b's: */

	construct_family("P","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar");

	/* QCD jet parton family definition: */

	construct_family("j","u,ubar,d,dbar,c,cbar,s,sbar");

	/* QCD jet parton family definition, including b's: */

	construct_family("J","u,ubar,d,dbar,c,cbar,s,sbar,b,bbar");

	/* QCD jet parton family definition, including b's and tops: */

	construct_family("Jt","u,ubar,d,dbar,c,cbar,s,sbar,t,tbar,b,bbar");
    }

    /* Recompute vertex couplings: */

    void EWSM::refresh_couplings()
    {
	refresh_weak_couplings();
    }

    /* Recompute Yukawa vertex couplings: */

    void EWSM::refresh_Yukawa_couplings()
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

    /* Refresh CKM matrix: */

    void EWSM::refresh_CKM_matrix()
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

    void EWSM::refresh_fermion_masses()
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

    void EWSM::refresh_weak_couplings()
    {
	Q_e=std::sqrt((value_type)4*pi*EWSM::alpha);
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

    void EWSM::refresh_widths()
    {
	refresh_W_width();
	refresh_Z_width();
	refresh_Higgs_width();
	refresh_top_width();
    }

    /* Recompute W-boson decay width: */

    void EWSM::refresh_W_width()
    {
	W_W=value_type(3+2*N_c)*G_F*std::pow(M_W,3)/((value_type)6*pi*std::sqrt((value_type)2));
    }

    /* Recompute Z-boson decay width: */

    void EWSM::refresh_Z_width()
    {
	value_type cw2=(M_W*M_W)/(M_Z*M_Z);
	value_type factor=G_F*std::pow(M_Z,3)/(108*std::sqrt((value_type)2)*pi);
	W_Z=(value_type(81+49*N_c)-(value_type)(162+92*N_c)*cw2+(value_type)(108+88*N_c)*cw2*cw2)*factor;
    }

    /* Recompute Higgs decay width: */

    void EWSM::refresh_Higgs_width()
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
	prefactor*=(SM_params::N_c);
	for(int i=0;i<6;++i)
	{
	    if(mq[i]>(value_type)0 and mq[i]<=(value_type)0.5*M_h0)
	    {
		W_h0+=prefactor*mq[i]*mq[i]*std::pow((value_type)1-std::pow(2*mq[i]/M_h0,(int)2),(value_type)1.5);
	    }
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

    void EWSM::refresh_top_width()
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

    void EWSM::discard_widths()
    {
	W_Z=(value_type)0;
	W_W=(value_type)0;
	W_h0=(value_type)0;
	W_t=(value_type)0;
	refresh_weak_couplings();
    }

    /* Higgs mass value input: */

    void EWSM::set_Higgs_mass(const value_type& m)
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

    void EWSM::set_alpha(const value_type& a)
    {
	alpha=a;
	refresh_weak_couplings();
    }

    /* Fine structure input: */

    void EWSM::set_alpha_s(const value_type& a){}

    /* Strong scale value input: */

    void EWSM::set_QCD_scale(const value_type& mu){}

    /* Reset all parameters using the default input parameters in the
     * SM_params.h file: */

    void EWSM::set_default_params()
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
	s12=(value_type)SM_params::s12;
	s23=(value_type)SM_params::s23;
	s13=(value_type)SM_params::s13;
	delta=(value_type)0;
	refresh_couplings();
	set_Feynman_gauge();
    }

    /* Function setting all vectors to the unitary gauge: */

    void EWSM::set_unitary_gauge()
    {
	if(initialised() and gauge!=1)
	{
	    decouple_particle("phi+");
	    decouple_particle("chi");
	    set_propagator<unitary_gauge>("W+");
	    set_propagator<unitary_gauge>("W-");
	    set_propagator<unitary_gauge>("Z");
	    set_propagator<unitary_gauge>("gamma");
	}
	gauge=1;
    }

    /* Function setting all vectors to the Feynman gauge and add would-be
     * Goldstone bosons if necessary: */

    void EWSM::set_Feynman_gauge()
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
	    set_propagator<scalar_propagator>("chi");
	    set_propagator<scalar_propagator>("phi+");
	    set_propagator<scalar_propagator>("phi-");
	}
	gauge=0;
    }
    
    /* Function setting all vectors to the R-xi gauge and add would-be
     * Goldstone bosons if necessary: */

    void EWSM::set_R_xi_gauge()
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

    void EWSM::set_xi(const value_type& x)
    {
	R_gauge<EWSM>::xi=x;
    }

    /* Function returning the number of QCD colours: */

    std::size_t EWSM::QCD_colours()
    {
	return N_c;
    }
}

