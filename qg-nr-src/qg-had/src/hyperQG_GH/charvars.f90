!$Id: charvars.f90,v 1.28 2008-01-21 14:43:04 carlos Exp $


#include "cctk.h"

#include "basic.h"
#include "Dus.h"
#include "GH_boundaries.h"

#if defined GH_MHD || defined GH_FLOWER
  module gh_charvars
#else
  module charvars
#endif

#if defined GH_MHD || defined GH_FLOWER
  use gh_rhs_mine
#else
  use rhs_mine
#endif

  use params
  implicit none
  private
#if defined GH_MHD || defined GH_FLOWER
  public gh_prim2char, gh_char2prim
#else
  public prim2char, char2prim
#endif
  
contains

  !----------------------------------------------------------------!
  !##################   PRIM2CHAR ##################################
  !----------------------------------------------------------------!


#if defined GH_MHD || defined GH_FLOWER
  subroutine gh_prim2char(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,&
                         uevolved_pt,uexact,dxu_pt,dyu_pt,dzu_pt,sources) 
#else
  subroutine prim2char(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,&
                       uevolved_pt,uexact,dxu_pt,dyu_pt,dzu_pt,sources) 
#endif
    implicit none

    CCTK_REAL, dimension(:), intent(inout) :: par
    CCTK_REAL, dimension(:), intent(in)    :: w_pt
    CCTK_REAL, dimension(:), intent(in)    :: v_pt
    character(len=*), intent(in)           :: direction   
    CCTK_REAL, dimension(:), intent(in)    :: u_pt,uevolved_pt,uexact,dxu_pt,dyu_pt,dzu_pt
    CCTK_REAL, dimension(:), intent(out)   :: sources    
    CCTK_REAL, dimension(:), intent(out)   :: w_in, w_out
    CCTK_REAL, dimension(:), intent(out)   :: w_0
    CCTK_REAL, dimension(3,3) :: T

    CCTK_INT                       :: a, b, c    
    CCTK_REAL, dimension(3)        :: nvec,pvec,qvec,nvecu           
    CCTK_REAL :: temp1, temp2, temp3, temp4                    
    CCTK_INT :: gauge_type, boundary_type, CPBC_type
  
!   parameters and constants
    CCTK_REAL :: alp, b1, b2, b3, bn, alp2, sqdetg, norm2
    CCTK_REAL :: sigma2, deth, facp, facm, sq2
    CCTK_REAL :: x, y, z, r    
    CCTK_REAL :: bg00, bg01, bg02, bg03
    CCTK_REAL :: dt1H0,dt2H0,dt3H0,dnH0,Gp,Gm
    CCTK_REAL :: dtphir1,dtphir2,dtphir3 
    CCTK_REAL :: dtphic1,dtphic2,dtphic3     
    CCTK_REAL :: dtphim1,dtphim2,dtphim3     
    CCTK_REAL :: dtphin1,dtphin2,dtphin3     
    CCTK_REAL :: Dn00,Dn01,Dn02,Dn03,Dn11,Dn12,Dn13,Dn22,Dn23,Dn33
    CCTK_REAL :: Dt100,Dt101,Dt102,Dt103,Dt111,Dt112,Dt113,Dt122,Dt123,Dt133    
    CCTK_REAL :: Dt200,Dt201,Dt202,Dt203,Dt211,Dt212,Dt213,Dt222,Dt223,Dt233
    CCTK_REAL :: Dt300,Dt301,Dt302,Dt303,Dt311,Dt312,Dt313,Dt322,Dt323,Dt333             
    CCTK_REAL :: Smm, Smp, Srm, Srp, Scm, Scp, Snp, Snm 
    CCTK_REAL :: dnphim, dnphin, dnphir, dnphic
    CCTK_REAL :: Lp00,Lp01,Lp02,Lp03,Lp11,Lp12,Lp13,Lp22,Lp23,Lp33
    CCTK_REAL :: Lm00,Lm01,Lm02,Lm03,Lm11,Lm12,Lm13,Lm22,Lm23,Lm33    

!   definitions for the CPBC    
    CCTK_REAL :: Kdd00,Kdd01,Kdd02,Kdd03,Kdd11,Kdd12,Kdd13,Kdd22,Kdd23,Kdd33
    CCTK_REAL :: Dddd100,Dddd101,Dddd102,Dddd103
    CCTK_REAL :: Dddd200,Dddd201,Dddd202,Dddd203
    CCTK_REAL :: Dddd300,Dddd301,Dddd302,Dddd303
    CCTK_REAL :: Dddd111,Dddd112,Dddd113,Dddd122,Dddd123,Dddd133
    CCTK_REAL :: Dddd211,Dddd212,Dddd213,Dddd222,Dddd223,Dddd233
    CCTK_REAL :: Dddd311,Dddd312,Dddd313,Dddd322,Dddd323,Dddd333
    CCTK_REAL :: Hd0,Hd1,Hd2,Hd3,Hu0,Hu1,Hu2,Hu3
    CCTK_REAL :: gdd00,gdd01,gdd02,gdd03,gdd11,gdd12,gdd13,gdd22,gdd23,gdd33
    CCTK_REAL :: guu00,guu01,guu02,guu03,guu11,guu12,guu13,guu22,guu23,guu33
    CCTK_REAL :: td0,td1,td2,td3,tu0,tu1,tu2,tu3
    CCTK_REAL :: nd0,nd1,nd2,nd3,nu0,nu1,nu2,nu3
    CCTK_REAL :: hdd00,hdd01,hdd02,hdd03,hdd11,hdd12,hdd13,hdd22,hdd23,hdd33    
    CCTK_REAL :: huu00,huu01,huu02,huu03,huu11,huu12,huu13,huu22,huu23,huu33
    CCTK_REAL :: hud10,hud11,hud12,hud13,hud20,hud21,hud22,hud23
    CCTK_REAL :: hud30,hud31,hud32,hud33,hud00,hud01,hud02,hud03
    CCTK_REAL :: kd0,kd1,kd2,kd3,ku0,ku1,ku2,ku3
    CCTK_REAL :: ld0,ld1,ld2,ld3,lu0,lu1,lu2,lu3    
    CCTK_REAL :: Pdd00,Pdd01,Pdd02,Pdd03,Pdd11,Pdd12,Pdd13,Pdd22,Pdd23,Pdd33
    CCTK_REAL :: Puu00,Puu01,Puu02,Puu03,Puu11,Puu12,Puu13,Puu22,Puu23,Puu33        
    CCTK_REAL :: Pud00,Pud01,Pud02,Pud03,Pud10,Pud11,Pud12,Pud13
    CCTK_REAL :: Pud20,Pud21,Pud22,Pud23,Pud30,Pud31,Pud32,Pud33                            
    CCTK_REAL :: Dddu100,Dddu101,Dddu102,Dddu103,Dddu110,Dddu111,Dddu112,Dddu113
    CCTK_REAL :: Dddu120,Dddu121,Dddu122,Dddu123,Dddu130,Dddu131,Dddu132,Dddu133
    CCTK_REAL :: Dddu200,Dddu201,Dddu202,Dddu203,Dddu210,Dddu211,Dddu212,Dddu213
    CCTK_REAL :: Dddu220,Dddu221,Dddu222,Dddu223,Dddu230,Dddu231,Dddu232,Dddu233
    CCTK_REAL :: Dddu300,Dddu301,Dddu302,Dddu303,Dddu310,Dddu311,Dddu312,Dddu313
    CCTK_REAL :: Dddu320,Dddu321,Dddu322,Dddu323,Dddu330,Dddu331,Dddu332,Dddu333
    
        
    CCTK_REAL :: D000,D001,D002,D003,D011,D012,D013,D022,D023,D033                  
    CCTK_REAL :: tD000,tD001,tD002,tD003,tD011,tD012,tD013,tD022,tD023,tD033    
    CCTK_REAL :: tD100,tD101,tD102,tD103,tD111,tD112,tD113,tD122,tD123,tD133    
    CCTK_REAL :: tD200,tD201,tD202,tD203,tD211,tD212,tD213,tD222,tD223,tD233    
    CCTK_REAL :: tD300,tD301,tD302,tD303,tD311,tD312,tD313,tD322,tD323,tD333    
    CCTK_REAL :: tK00,tK01,tK02,tK03,tK11,tK12,tK13,tK22,tK23,tK33            
        
    CCTK_REAL :: sDt100,sDt101,sDt102,sDt103,sDt111,sDt112,sDt113,sDt122,sDt123,sDt133    
    CCTK_REAL :: sDt200,sDt201,sDt202,sDt203,sDt211,sDt212,sDt213,sDt222,sDt223,sDt233
    CCTK_REAL :: sDt300,sDt301,sDt302,sDt303,sDt311,sDt312,sDt313,sDt322,sDt323,sDt333         
    CCTK_REAL :: sdtphir1,sdtphir2,sdtphir3 
    CCTK_REAL :: sdtphic1,sdtphic2,sdtphic3     
    CCTK_REAL :: sdtphim1,sdtphim2,sdtphim3     
    CCTK_REAL :: sdtphin1,sdtphin2,sdtphin3     
    CCTK_REAL :: sdt1H0,sdt2H0,sdt3H0
        
    CCTK_REAL :: CD100,CD101,CD102,CD103,CD111,CD112,CD113,CD122,CD123,CD133
    CCTK_REAL :: CD200,CD201,CD202,CD203,CD211,CD212,CD213,CD222,CD223,CD233
    CCTK_REAL :: CD300,CD301,CD302,CD303,CD311,CD312,CD313,CD322,CD323,CD333
    
    CCTK_REAL :: CDD1100,CDD1101,CDD1102,CDD1103
    CCTK_REAL :: CDD1200,CDD1201,CDD1202,CDD1203
    CCTK_REAL :: CDD1300,CDD1301,CDD1302,CDD1303
    CCTK_REAL :: CDD2100,CDD2101,CDD2102,CDD2103    
    CCTK_REAL :: CDD2200,CDD2201,CDD2202,CDD2203
    CCTK_REAL :: CDD2300,CDD2301,CDD2302,CDD2303
    CCTK_REAL :: CDD3100,CDD3101,CDD3102,CDD3103
    CCTK_REAL :: CDD3200,CDD3201,CDD3202,CDD3203        
    CCTK_REAL :: CDD3300,CDD3301,CDD3302,CDD3303
    CCTK_REAL :: CDD1111,CDD1112,CDD1113,CDD1122,CDD1123,CDD1133
    CCTK_REAL :: CDD1211,CDD1212,CDD1213,CDD1222,CDD1223,CDD1233
    CCTK_REAL :: CDD1311,CDD1312,CDD1313,CDD1322,CDD1323,CDD1333
    CCTK_REAL :: CDD2111,CDD2112,CDD2113,CDD2122,CDD2123,CDD2133    
    CCTK_REAL :: CDD2211,CDD2212,CDD2213,CDD2222,CDD2223,CDD2233
    CCTK_REAL :: CDD2311,CDD2312,CDD2313,CDD2322,CDD2323,CDD2333
    CCTK_REAL :: CDD3111,CDD3112,CDD3113,CDD3122,CDD3123,CDD3133
    CCTK_REAL :: CDD3211,CDD3212,CDD3213,CDD3222,CDD3223,CDD3233
    CCTK_REAL :: CDD3311,CDD3312,CDD3313,CDD3322,CDD3323,CDD3333        
    
    CCTK_REAL :: CDDn100,CDDn101,CDDn102,CDDn103
    CCTK_REAL :: CDDn200,CDDn201,CDDn202,CDDn203
    CCTK_REAL :: CDDn300,CDDn301,CDDn302,CDDn303    
    CCTK_REAL :: CDDn111,CDDn112,CDDn113,CDDn122,CDDn123,CDDn133
    CCTK_REAL :: CDDn211,CDDn212,CDDn213,CDDn222,CDDn223,CDDn233
    CCTK_REAL :: CDDn311,CDDn312,CDDn313,CDDn322,CDDn323,CDDn333

    CCTK_REAL :: CPPm11,CPPm12,CPPm13,CPPm21,CPPm22,CPPm23,CPPm31,CPPm32,CPPm33
    CCTK_REAL :: CPPn11,CPPn12,CPPn13,CPPn21,CPPn22,CPPn23,CPPn31,CPPn32,CPPn33
    CCTK_REAL :: CPPr11,CPPr12,CPPr13,CPPr21,CPPr22,CPPr23,CPPr31,CPPr32,CPPr33
    CCTK_REAL :: CPPc11,CPPc12,CPPc13,CPPc21,CPPc22,CPPc23,CPPc31,CPPc32,CPPc33
    CCTK_REAL :: CGGm11,CGGm12,CGGm13,CGGm21,CGGm22,CGGm23,CGGm31,CGGm32,CGGm33
            
    CCTK_REAL :: CPPmn1,CPPmn2,CPPmn3,CPPmnt1,CPPmnt2,CPPmnt3
    CCTK_REAL :: CPPnn1,CPPnn2,CPPnn3,CPPnnt1,CPPnnt2,CPPnnt3
    CCTK_REAL :: CPPrn1,CPPrn2,CPPrn3,CPPrnt1,CPPrnt2,CPPrnt3
    CCTK_REAL :: CPPcn1,CPPcn2,CPPcn3,CPPcnt1,CPPcnt2,CPPcnt3             
    CCTK_REAL :: CGGmn1,CGGmn2,CGGmn3,CGGmnt1,CGGmnt2,CGGmnt3    
                        
    CCTK_REAL :: trD1,trD2,trD3,ttD1,ttD2,ttD3
    CCTK_REAL :: tD10,tD11,tD12,tD13,tD20,tD21,tD22,tD23,tD30,tD31,tD32,tD33
    CCTK_REAL :: Kud00,Kud01,Kud02,Kud03,Kud10,Kud11,Kud12,Kud13
    CCTK_REAL :: Kud20,Kud21,Kud22,Kud23,Kud30,Kud31,Kud32,Kud33                        
    CCTK_REAL :: trK,tK0,tK1,tK2,tK3,trCD1,trCD2,trCD3
    
    CCTK_REAL :: Fd0,Fd1,Fd2,Fd3,Cm0,Cm1,Cm2,Cm3
    CCTK_REAL :: Cdd10,Cdd11,Cdd12,Cdd13,Cdd20,Cdd21,Cdd22,Cdd23
    CCTK_REAL :: Cdd30,Cdd31,Cdd32,Cdd33
    CCTK_REAL :: sL00,sL01,sL02,sL03,sL11,sL12,sL13,sL22,sL23,sL33    
    CCTK_REAL :: sSr,sSc,sSm,sSn,sG    
        
!-  reading the parameters ------------------       
        
    sigma2        = par(P_SIGMA2)
    gauge_type    = int(par(P_GAUGE_TYPE))    
    boundary_type = int(par(P_BOUNDARY_CONDITIONS))
    CPBC_type     = int(par(P_CPBC_TYPE))    
    sq2           = sqrt(2.0d0)

!-  reading the position ------------------           
    x = w_pt(H_XPHYS)
    y = w_pt(H_YPHYS)
    z = w_pt(H_ZPHYS)   
    r = sqrt(x**2+y**2+z**2)    
    
!-  computing the 3+1 quantities--------------
   
    bg00 = uevolved_pt(H_G00)       
    bg01 = uevolved_pt(H_G01)   
    bg02 = uevolved_pt(H_G02)   
    bg03 = uevolved_pt(H_G03)           
    hdd11 = uevolved_pt(H_G11)
    hdd12 = uevolved_pt(H_G12)
    hdd13 = uevolved_pt(H_G13)
    hdd22 = uevolved_pt(H_G22)
    hdd23 = uevolved_pt(H_G23)
    hdd33 = uevolved_pt(H_G33)
                                          
    huu11 = -hdd23**2 + hdd22*hdd33
    huu12 = hdd13*hdd23 - hdd12*hdd33
    huu13 = -(hdd13*hdd22) + hdd12*hdd23
    huu22 = -hdd13**2 + hdd11*hdd33
    huu23 = hdd12*hdd13 - hdd11*hdd23
    huu33 = -hdd12**2 + hdd11*hdd22
    deth = hdd11*huu11 + hdd12*huu12 + hdd13*huu13
    huu11 = huu11/deth
    huu12 = huu12/deth
    huu13 = huu13/deth
    huu22 = huu22/deth
    huu23 = huu23/deth
    huu33 = huu33/deth

    b1 = bg01*huu11 + bg02*huu12 + bg03*huu13
    b2 = bg01*huu12 + bg02*huu22 + bg03*huu23
    b3 = bg01*huu13 + bg02*huu23 + bg03*huu33
    alp2 = -bg00 + b1**2*hdd11 + b2**2*hdd22 + 2*(b1*(b2*hdd12 + b3*hdd13) &
   &     + b2*b3*hdd23) + b3**2*hdd33
    alp = sqrt(alp2)       

!-  computing the normal to the surface--------------

    call what_boundary(direction,nvec,pvec,qvec)        
           
    nd0 = 0.0d0; nd1 = nvec(1); nd2 = nvec(2); nd3 = nvec(3)     
        
    norm2 = huu11*nd1**2 + huu22*nd2**2 + huu33*nd3**2 + 2*(huu23*nd2*nd3&
     & + nd1*(huu12*nd2 + huu13*nd3))
    nd1 = nd1/sqrt(norm2)
    nd2 = nd2/sqrt(norm2)
    nd3 = nd3/sqrt(norm2)
    nu1 = huu11*nd1 + huu12*nd2 + huu13*nd3
    nu2 = huu12*nd1 + huu22*nd2 + huu23*nd3
    nu3 = huu13*nd1 + huu23*nd2 + huu33*nd3
    
! computing projections -------------------------------  

    bn = b1*nd1 + b2*nd2 + b3*nd3

    dnphir = dphir1*nu1 + dphir2*nu2 + dphir3*nu3
    dnphic = dphic1*nu1 + dphic2*nu2 + dphic3*nu3
    dnphim = dphim1*nu1 + dphim2*nu2 + dphim3*nu3    
    dnphin = dphin1*nu1 + dphin2*nu2 + dphin3*nu3    
    dnH0 = d1H0*nu1 + d2H0*nu2 + d3H0*nu3    
            
    dtphir1 = dphir1 - dnphir*nd1
    dtphir2 = dphir2 - dnphir*nd2
    dtphir3 = dphir3 - dnphir*nd3
    
    dtphic1 = dphic1 - dnphic*nd1
    dtphic2 = dphic2 - dnphic*nd2
    dtphic3 = dphic3 - dnphic*nd3

    dtphim1 = dphim1 - dnphim*nd1
    dtphim2 = dphim2 - dnphim*nd2
    dtphim3 = dphim3 - dnphim*nd3

    dtphin1 = dphin1 - dnphin*nd1
    dtphin2 = dphin2 - dnphin*nd2
    dtphin3 = dphin3 - dnphin*nd3

    dt1H0 = d1H0 - dnH0*nd1
    dt2H0 = d2H0 - dnH0*nd2
    dt3H0 = d3H0 - dnH0*nd3
        
!------------------------------------    
    
    Dn00 = nu1*D100 + nu2*D200 + nu3*D300    
    Dn01 = nu1*D101 + nu2*D201 + nu3*D301
    Dn02 = nu1*D102 + nu2*D202 + nu3*D302
    Dn03 = nu1*D103 + nu2*D203 + nu3*D303
    Dn11 = nu1*D111 + nu2*D211 + nu3*D311
    Dn12 = nu1*D112 + nu2*D212 + nu3*D312
    Dn13 = nu1*D113 + nu2*D213 + nu3*D313
    Dn22 = nu1*D122 + nu2*D222 + nu3*D322
    Dn23 = nu1*D123 + nu2*D223 + nu3*D323
    Dn33 = nu1*D133 + nu2*D233 + nu3*D333
  
    Dt100 = D100 - nd1*Dn00    
    Dt101 = D101 - nd1*Dn01
    Dt102 = D102 - nd1*Dn02
    Dt103 = D103 - nd1*Dn03
    Dt111 = D111 - nd1*Dn11
    Dt112 = D112 - nd1*Dn12
    Dt113 = D113 - nd1*Dn13
    Dt122 = D122 - nd1*Dn22
    Dt123 = D123 - nd1*Dn23
    Dt133 = D133 - nd1*Dn33

    Dt200 = D200 - nd2*Dn00    
    Dt201 = D201 - nd2*Dn01
    Dt202 = D202 - nd2*Dn02
    Dt203 = D203 - nd2*Dn03
    Dt211 = D211 - nd2*Dn11
    Dt212 = D212 - nd2*Dn12
    Dt213 = D213 - nd2*Dn13   
    Dt222 = D222 - nd2*Dn22
    Dt223 = D223 - nd2*Dn23
    Dt233 = D233 - nd2*Dn33

    Dt300 = D300 - nd3*Dn00    
    Dt301 = D301 - nd3*Dn01
    Dt302 = D302 - nd3*Dn02
    Dt303 = D303 - nd3*Dn03  
    Dt311 = D311 - nd3*Dn11
    Dt312 = D312 - nd3*Dn12
    Dt313 = D313 - nd3*Dn13
    Dt322 = D322 - nd3*Dn22
    Dt323 = D323 - nd3*Dn23
    Dt333 = D333 - nd3*Dn33

! the eigenvectors -------------------------------------  
    
    Lp00 = ( K00 - sigma2*g00 ) + ( Dn00 )
    Lp01 = ( K01 - sigma2*g01 ) + ( Dn01 )
    Lp02 = ( K02 - sigma2*g02 ) + ( Dn02 )
    Lp03 = ( K03 - sigma2*g03 ) + ( Dn03 )
    Lp11 = ( K11 - sigma2*g11 ) + ( Dn11 )
    Lp12 = ( K12 - sigma2*g12 ) + ( Dn12 )
    Lp13 = ( K13 - sigma2*g13 ) + ( Dn13 )
    Lp22 = ( K22 - sigma2*g22 ) + ( Dn22 )
    Lp23 = ( K23 - sigma2*g23 ) + ( Dn23 )
    Lp33 = ( K33 - sigma2*g33 ) + ( Dn33 )        
                            
    Lm00 = ( K00 - sigma2*g00 ) - ( Dn00 )
    Lm01 = ( K01 - sigma2*g01 ) - ( Dn01 )
    Lm02 = ( K02 - sigma2*g02 ) - ( Dn02 )
    Lm03 = ( K03 - sigma2*g03 ) - ( Dn03 )
    Lm11 = ( K11 - sigma2*g11 ) - ( Dn11 )
    Lm12 = ( K12 - sigma2*g12 ) - ( Dn12 )
    Lm13 = ( K13 - sigma2*g13 ) - ( Dn13 )
    Lm22 = ( K22 - sigma2*g22 ) - ( Dn22 )
    Lm23 = ( K23 - sigma2*g23 ) - ( Dn23 )
    Lm33 = ( K33 - sigma2*g33 ) - ( Dn33 )        

    Smp = ( pim - sigma2*phim ) + ( dnphim )
    Smm = ( pim - sigma2*phim ) - ( dnphim )

    Snp = ( pin - sigma2*phin ) + ( dnphin )
    Snm = ( pin - sigma2*phin ) - ( dnphin )
    
    Srp = ( pir - sigma2*phir ) + ( dnphir )
    Srm = ( pir - sigma2*phir ) - ( dnphir )

    Scp = ( pic - sigma2*phic ) + ( dnphic )
    Scm = ( pic - sigma2*phic ) - ( dnphic )
    
    Gp = ( G0 - sigma2*H0 ) + ( dnH0 )
    Gm = ( G0 - sigma2*H0 ) - ( dnH0 )
            
! copying the eigenvectors -------------------------  
!---------------------------------------------------

!--ingoing eigenvectors-----------------------------      
    w_in(1)  = Lm00
    w_in(2)  = Lm01
    w_in(3)  = Lm02
    w_in(4)  = Lm03
    w_in(5)  = Lm11
    w_in(6)  = Lm12
    w_in(7)  = Lm13
    w_in(8)  = Lm22      
    w_in(9)  = Lm23        
    w_in(10) = Lm33
    w_in(11) = Smm
    w_in(12) = Snm
    w_in(13) = Srm
    w_in(14) = Scm        
    w_in(15) = Gm        
    
!--outgoing eigenvectors-----------------------------        
    w_out(1)  = Lp00
    w_out(2)  = Lp01
    w_out(3)  = Lp02
    w_out(4)  = Lp03
    w_out(5)  = Lp11
    w_out(6)  = Lp12         
    w_out(7)  = Lp13        
    w_out(8)  = Lp22      
    w_out(9)  = Lp23        
    w_out(10) = Lp33
    w_out(11) = Smp
    w_out(12) = Snp
    w_out(13) = Srp
    w_out(14) = Scp
    w_out(15) = Gp
    
!--standing eigenvectors-----------------------------    
    w_0(1)  = g00
    w_0(2)  = g01
    w_0(3)  = g02
    w_0(4)  = g03
    w_0(5)  = g11
    w_0(6)  = g12        
    w_0(7)  = g13       
    w_0(8)  = g22    
    w_0(9)  = g23        
    w_0(10) = g33
    
    w_0(11)  = H0
    w_0(12)  = H1
    w_0(13)  = H2
    w_0(14)  = H3
    
!-------------------------------------------------------------------------- 
!--if bn > 0 or CPBC then the u^2_kab modes are ingoing -------------------
!-------------------------------------------------------------------------- 

    if (bn .GT. 0) then
      
      w_in(21) = Dt100
      w_in(22) = Dt101
      w_in(23) = Dt102
      w_in(24) = Dt103
      w_in(25) = Dt111
      w_in(26) = Dt112        
      w_in(27) = Dt113       
      w_in(28) = Dt122    
      w_in(29) = Dt123        
      w_in(30) = Dt133

      w_in(31) = Dt200
      w_in(32) = Dt201
      w_in(33) = Dt202
      w_in(34) = Dt203
      w_in(35) = Dt211
      w_in(36) = Dt212        
      w_in(37) = Dt213       
      w_in(38) = Dt222    
      w_in(39) = Dt223        
      w_in(40) = Dt233

      w_in(41) = Dt300
      w_in(42) = Dt301
      w_in(43) = Dt302
      w_in(44) = Dt303
      w_in(45) = Dt311
      w_in(46) = Dt312        
      w_in(47) = Dt313       
      w_in(48) = Dt322    
      w_in(49) = Dt323        
      w_in(50) = Dt333

      w_in(51) = dtphim1
      w_in(52) = dtphim2
      w_in(53) = dtphim3
      w_in(54) = dtphin1
      w_in(55) = dtphin2
      w_in(56) = dtphin3
      w_in(57) = dtphir1
      w_in(58) = dtphir2
      w_in(59) = dtphir3
      w_in(60) = dtphic1
      w_in(61) = dtphic2
      w_in(62) = dtphic3

      w_in(63) = dt1H0
      w_in(64) = dt2H0
      w_in(65) = dt3H0
      
                  
! otherwise they are outgoing, easier                  
    else 
    
      w_0(21) = Dt100
      w_0(22) = Dt101
      w_0(23) = Dt102
      w_0(24) = Dt103
      w_0(25) = Dt111
      w_0(26) = Dt112        
      w_0(27) = Dt113       
      w_0(28) = Dt122    
      w_0(29) = Dt123        
      w_0(30) = Dt133

      w_0(31) = Dt200
      w_0(32) = Dt201
      w_0(33) = Dt202
      w_0(34) = Dt203
      w_0(35) = Dt211
      w_0(36) = Dt212        
      w_0(37) = Dt213       
      w_0(38) = Dt222    
      w_0(39) = Dt223        
      w_0(40) = Dt233

      w_0(41) = Dt300
      w_0(42) = Dt301
      w_0(43) = Dt302
      w_0(44) = Dt303
      w_0(45) = Dt311
      w_0(46) = Dt312        
      w_0(47) = Dt313       
      w_0(48) = Dt322    
      w_0(49) = Dt323        
      w_0(50) = Dt333

      w_0(51) = dtphim1
      w_0(52) = dtphim2
      w_0(53) = dtphim3
      w_0(54) = dtphin1
      w_0(55) = dtphin2
      w_0(56) = dtphin3
      w_0(57) = dtphir1
      w_0(58) = dtphir2
      w_0(59) = dtphir3
      w_0(60) = dtphic1
      w_0(61) = dtphic2
      w_0(62) = dtphic3

      w_0(63) = dt1H0
      w_0(64) = dt2H0
      w_0(65) = dt3H0
      
                          
    end if      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  now the sources
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    sources = 0.0
    
    if (boundary_type.eq.7) then                        
    
#include "GH_boundaries_CPBC.inc"    
    
      sources(1)  = sL00
      sources(2)  = sL01
      sources(3)  = sL02
      sources(4)  = sL03
      sources(5)  = sL11
      sources(6)  = sL12
      sources(7)  = sL13
      sources(8)  = sL22
      sources(9)  = sL23
      sources(10) = sL33    
      sources(11) = sSm
      sources(12) = sSn
      sources(13) = sSr
      sources(14) = sSc
      sources(15) = sG
                                               
      sources(21) = sDt100
      sources(22) = sDt101
      sources(23) = sDt102
      sources(24) = sDt103
      sources(25) = sDt111
      sources(26) = sDt112
      sources(27) = sDt113
      sources(28) = sDt122
      sources(29) = sDt123
      sources(30) = sDt133

      sources(31) = sDt200
      sources(32) = sDt201
      sources(33) = sDt202
      sources(34) = sDt203
      sources(35) = sDt211
      sources(36) = sDt212
      sources(37) = sDt213
      sources(38) = sDt222
      sources(39) = sDt223
      sources(40) = sDt233

      sources(41) = sDt300
      sources(42) = sDt301
      sources(43) = sDt302
      sources(44) = sDt303
      sources(45) = sDt311
      sources(46) = sDt312
      sources(47) = sDt313
      sources(48) = sDt322
      sources(49) = sDt323
      sources(50) = sDt333
                        
      sources(51) = sdtphim1
      sources(52) = sdtphim2
      sources(53) = sdtphim3
      sources(54) = sdtphin1
      sources(55) = sdtphin2
      sources(56) = sdtphin3
      sources(57) = sdtphir1
      sources(58) = sdtphir2
      sources(59) = sdtphir3
      sources(60) = sdtphic1
      sources(61) = sdtphic2
      sources(62) = sdtphic3      

      sources(63) = sdt1H0
      sources(64) = sdt2H0
      sources(65) = sdt3H0
      
                
    end if
            
  end subroutine

  
  !----------------------------------------------------------------!
  !#################     CHAR2PRIM    #####################################!
  !----------------------------------------------------------------!  

#if defined GH_MHD || defined GH_FLOWER
  subroutine gh_char2prim(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,&
                       T,uevolved_pt,uexact) 
#else
  subroutine char2prim(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,&
                       T,uevolved_pt,uexact) 
#endif
    implicit none

    CCTK_REAL, dimension(:), intent(inout) :: par
    CCTK_REAL, dimension(:), intent(in)    :: w_pt, uevolved_pt,uexact,v_pt
    character(len=*), intent(in)           :: direction   
    CCTK_REAL, dimension(:), intent(out)   :: u_pt
    CCTK_REAL, dimension(:), intent(in)    :: w_in, w_out
    CCTK_REAL, dimension(:), intent(in)    :: w_0
    CCTK_REAL, dimension(3,3) :: T
 
    CCTK_REAL, dimension(3)    :: nvec,pvec,qvec,nvecu       
    CCTK_INT :: gauge_type, boundary_type

!   parameters and constants    
    CCTK_REAL :: alp, b1, b2, b3, bn, alp2, sqdetg, norm2
    CCTK_REAL :: sigma2, deth, facp, facm, sq2
    CCTK_REAL :: bg00,bg01,bg02,bg03
    CCTK_REAL :: dt1H0,dt2H0,dt3H0,dnH0,Gp,Gm    
    CCTK_REAL :: dtphir1,dtphir2,dtphir3 
    CCTK_REAL :: dtphic1,dtphic2,dtphic3     
    CCTK_REAL :: dtphim1,dtphim2,dtphim3     
    CCTK_REAL :: dtphin1,dtphin2,dtphin3     
    CCTK_REAL :: Dn00,Dn01,Dn02,Dn03,Dn11,Dn12,Dn13,Dn22,Dn23,Dn33
    CCTK_REAL :: Dt100,Dt101,Dt102,Dt103,Dt111,Dt112,Dt113,Dt122,Dt123,Dt133    
    CCTK_REAL :: Dt200,Dt201,Dt202,Dt203,Dt211,Dt212,Dt213,Dt222,Dt223,Dt233
    CCTK_REAL :: Dt300,Dt301,Dt302,Dt303,Dt311,Dt312,Dt313,Dt322,Dt323,Dt333         
    CCTK_REAL :: Smm, Snm, Smp, Snp, Srm, Srp, Scm, Scp 
    CCTK_REAL :: dnphim, dnphin, dnphir, dnphic
    CCTK_REAL :: Lp00,Lp01,Lp02,Lp03,Lp11,Lp12,Lp13,Lp22,Lp23,Lp33
    CCTK_REAL :: Lm00,Lm01,Lm02,Lm03,Lm11,Lm12,Lm13,Lm22,Lm23,Lm33    
    
    
!   definitions for the CPBC    
    CCTK_REAL :: Kdd00,Kdd01,Kdd02,Kdd03,Kdd11,Kdd12,Kdd13,Kdd22,Kdd23,Kdd33
    CCTK_REAL :: Dddd100,Dddd101,Dddd102,Dddd103
    CCTK_REAL :: Dddd200,Dddd201,Dddd202,Dddd203
    CCTK_REAL :: Dddd300,Dddd301,Dddd302,Dddd303
    CCTK_REAL :: Dddd111,Dddd112,Dddd113,Dddd122,Dddd123,Dddd133
    CCTK_REAL :: Dddd211,Dddd212,Dddd213,Dddd222,Dddd223,Dddd233
    CCTK_REAL :: Dddd311,Dddd312,Dddd313,Dddd322,Dddd323,Dddd333
    CCTK_REAL :: Hd0,Hd1,Hd2,Hd3,Hu0,Hu1,Hu2,Hu3
    CCTK_REAL :: gdd00,gdd01,gdd02,gdd03,gdd11,gdd12,gdd13,gdd22,gdd23,gdd33
    CCTK_REAL :: guu00,guu01,guu02,guu03,guu11,guu12,guu13,guu22,guu23,guu33
    CCTK_REAL :: td0,td1,td2,td3,tu0,tu1,tu2,tu3
    CCTK_REAL :: nd0,nd1,nd2,nd3,nu0,nu1,nu2,nu3
    CCTK_REAL :: hdd00,hdd01,hdd02,hdd03,hdd11,hdd12,hdd13,hdd22,hdd23,hdd33    
    CCTK_REAL :: huu00,huu01,huu02,huu03,huu11,huu12,huu13,huu22,huu23,huu33
    CCTK_REAL :: hud10,hud11,hud12,hud13,hud20,hud21,hud22,hud23
    CCTK_REAL :: hud30,hud31,hud32,hud33
    CCTK_REAL :: kd0,kd1,kd2,kd3,ku0,ku1,ku2,ku3
    CCTK_REAL :: ld0,ld1,ld2,ld3,lu0,lu1,lu2,lu3    
    CCTK_REAL :: Pdd00,Pdd01,Pdd02,Pdd03,Pdd11,Pdd12,Pdd13,Pdd22,Pdd23,Pdd33
    CCTK_REAL :: Puu00,Puu01,Puu02,Puu03,Puu11,Puu12,Puu13,Puu22,Puu23,Puu33        
    CCTK_REAL :: Pud00,Pud01,Pud02,Pud03,Pud10,Pud11,Pud12,Pud13
    CCTK_REAL :: Pud20,Pud21,Pud22,Pud23,Pud30,Pud31,Pud32,Pud33                        

    
!-  reading the parameters ------------------       
        
    sigma2        = par(P_SIGMA2)
    gauge_type    = int(par(P_GAUGE_TYPE))    
    boundary_type = int(par(P_BOUNDARY_CONDITIONS))
    sq2           = sqrt(2.0d0)
                   
!-  computing the 3+1 quantities--------------
   
    bg00 = uevolved_pt(H_G00)       
    bg01 = uevolved_pt(H_G01)   
    bg02 = uevolved_pt(H_G02)   
    bg03 = uevolved_pt(H_G03)           
    hdd11 = uevolved_pt(H_G11)
    hdd12 = uevolved_pt(H_G12)
    hdd13 = uevolved_pt(H_G13)
    hdd22 = uevolved_pt(H_G22)
    hdd23 = uevolved_pt(H_G23)
    hdd33 = uevolved_pt(H_G33)
                                          
    huu11 = -hdd23**2 + hdd22*hdd33
    huu12 = hdd13*hdd23 - hdd12*hdd33
    huu13 = -(hdd13*hdd22) + hdd12*hdd23
    huu22 = -hdd13**2 + hdd11*hdd33
    huu23 = hdd12*hdd13 - hdd11*hdd23
    huu33 = -hdd12**2 + hdd11*hdd22
    deth = hdd11*huu11 + hdd12*huu12 + hdd13*huu13
    huu11 = huu11/deth
    huu12 = huu12/deth
    huu13 = huu13/deth
    huu22 = huu22/deth
    huu23 = huu23/deth
    huu33 = huu33/deth

    b1 = bg01*huu11 + bg02*huu12 + bg03*huu13
    b2 = bg01*huu12 + bg02*huu22 + bg03*huu23
    b3 = bg01*huu13 + bg02*huu23 + bg03*huu33
    alp2 = -bg00 + b1**2*hdd11 + b2**2*hdd22 + 2*(b1*(b2*hdd12 + b3*hdd13) &
   &     + b2*b3*hdd23) + b3**2*hdd33
    alp = sqrt(alp2)       

!-  computing the normal to the surface--------------

    call what_boundary(direction,nvec,pvec,qvec)        
           
    nd0 = 0.0d0; nd1 = nvec(1); nd2 = nvec(2); nd3 = nvec(3) 
    
    norm2 = huu11*nd1**2 + huu22*nd2**2 + huu33*nd3**2 + 2*(huu23*nd2*nd3&
     & + nd1*(huu12*nd2 + huu13*nd3))
    nd1 = nd1/sqrt(norm2)
    nd2 = nd2/sqrt(norm2)
    nd3 = nd3/sqrt(norm2)
    nu1 = huu11*nd1 + huu12*nd2 + huu13*nd3
    nu2 = huu12*nd1 + huu22*nd2 + huu23*nd3
    nu3 = huu13*nd1 + huu23*nd2 + huu33*nd3
    
    bn = b1*nd1 + b2*nd2 + b3*nd3

!----------------------------------------------------------
! copying the eigenvectors -------------------------  
    
    Lm00 = w_in(1)
    Lm01 = w_in(2)
    Lm02 = w_in(3)
    Lm03 = w_in(4)
    Lm11 = w_in(5)
    Lm12 = w_in(6)
    Lm13 = w_in(7)
    Lm22 = w_in(8)
    Lm23 = w_in(9)
    Lm33 = w_in(10)
    Smm  = w_in(11)                                        
    Snm  = w_in(12) 
    Srm  = w_in(13)
    Scm  = w_in(14)        
    Gm   = w_in(15)            

    Lp00 = w_out(1)
    Lp01 = w_out(2)
    Lp02 = w_out(3)
    Lp03 = w_out(4)
    Lp11 = w_out(5)
    Lp12 = w_out(6)
    Lp13 = w_out(7)
    Lp22 = w_out(8)
    Lp23 = w_out(9)
    Lp33 = w_out(10)
    Smp  = w_out(11)                                        
    Snp  = w_out(12)                                        
    Srp  = w_out(13)
    Scp  = w_out(14)        
    Gp   = w_out(15)                
          
    if (bn .GT. 0)  then
      
      Dt100 = w_in(21)
      Dt101 = w_in(22)
      Dt102 = w_in(23)
      Dt103 = w_in(24)
      Dt111 = w_in(25)
      Dt112 = w_in(26)
      Dt113 = w_in(27)
      Dt122 = w_in(28)
      Dt123 = w_in(29)
      Dt133 = w_in(30)

      Dt200 = w_in(31)
      Dt201 = w_in(32)
      Dt202 = w_in(33)
      Dt203 = w_in(34)
      Dt211 = w_in(35)
      Dt212 = w_in(36)
      Dt213 = w_in(37)
      Dt222 = w_in(38)
      Dt223 = w_in(39)
      Dt233 = w_in(40)

      Dt300 = w_in(41)
      Dt301 = w_in(42)
      Dt302 = w_in(43)
      Dt303 = w_in(44)
      Dt311 = w_in(45)
      Dt312 = w_in(46)
      Dt313 = w_in(47)
      Dt322 = w_in(48)
      Dt323 = w_in(49)
      Dt333 = w_in(50)
      
      dtphim1 = w_in(51)
      dtphim2 = w_in(52)
      dtphim3 = w_in(53)            
      dtphin1 = w_in(54)
      dtphin2 = w_in(55)
      dtphin3 = w_in(56)            
      dtphir1 = w_in(57)
      dtphir2 = w_in(58)
      dtphir3 = w_in(59)            
      dtphic1 = w_in(60)
      dtphic2 = w_in(61)
      dtphic3 = w_in(62)            

      dt1H0 = w_in(63)
      dt2H0 = w_in(64)
      dt3H0 = w_in(65)            
      
                  
    else 
    
      Dt100 = w_0(21)
      Dt101 = w_0(22)
      Dt102 = w_0(23)
      Dt103 = w_0(24)
      Dt111 = w_0(25)
      Dt112 = w_0(26)
      Dt113 = w_0(27)
      Dt122 = w_0(28)
      Dt123 = w_0(29)
      Dt133 = w_0(30)

      Dt200 = w_0(31)
      Dt201 = w_0(32)
      Dt202 = w_0(33)
      Dt203 = w_0(34)
      Dt211 = w_0(35)
      Dt212 = w_0(36)
      Dt213 = w_0(37)
      Dt222 = w_0(38)
      Dt223 = w_0(39)
      Dt233 = w_0(40)

      Dt300 = w_0(41)
      Dt301 = w_0(42)
      Dt302 = w_0(43)
      Dt303 = w_0(44)
      Dt311 = w_0(45)
      Dt312 = w_0(46)
      Dt313 = w_0(47)
      Dt322 = w_0(48)
      Dt323 = w_0(49)
      Dt333 = w_0(50)

      dtphim1 = w_0(51)
      dtphim2 = w_0(52)
      dtphim3 = w_0(53)            
      dtphin1 = w_0(54)
      dtphin2 = w_0(55)
      dtphin3 = w_0(56)            
      dtphir1 = w_0(57)
      dtphir2 = w_0(58)
      dtphir3 = w_0(59)            
      dtphic1 = w_0(60)
      dtphic2 = w_0(61)
      dtphic3 = w_0(62)            

      dt1H0 = w_0(63)
      dt2H0 = w_0(64)
      dt3H0 = w_0(65)                  
                          
    end if                     
    
!-  recovering the primitive fields --------------
    
    K00 = 0.5d0*(Lp00 + Lm00) + sigma2*g00 
    K01 = 0.5d0*(Lp01 + Lm01) + sigma2*g01 
    K02 = 0.5d0*(Lp02 + Lm02) + sigma2*g02 
    K03 = 0.5d0*(Lp03 + Lm03) + sigma2*g03 
    K11 = 0.5d0*(Lp11 + Lm11) + sigma2*g11 
    K12 = 0.5d0*(Lp12 + Lm12) + sigma2*g12 
    K13 = 0.5d0*(Lp13 + Lm13) + sigma2*g13 
    K22 = 0.5d0*(Lp22 + Lm22) + sigma2*g22 
    K23 = 0.5d0*(Lp23 + Lm23) + sigma2*g23 
    K33 = 0.5d0*(Lp33 + Lm33) + sigma2*g33
    
    Dn00 = 0.5d0*(Lp00 - Lm00) 
    Dn01 = 0.5d0*(Lp01 - Lm01) 
    Dn02 = 0.5d0*(Lp02 - Lm02) 
    Dn03 = 0.5d0*(Lp03 - Lm03) 
    Dn11 = 0.5d0*(Lp11 - Lm11) 
    Dn12 = 0.5d0*(Lp12 - Lm12) 
    Dn13 = 0.5d0*(Lp13 - Lm13) 
    Dn22 = 0.5d0*(Lp22 - Lm22) 
    Dn23 = 0.5d0*(Lp23 - Lm23) 
    Dn33 = 0.5d0*(Lp33 - Lm33) 
    
    D100 = Dt100 + nd1*Dn00    
    D101 = Dt101 + nd1*Dn01
    D102 = Dt102 + nd1*Dn02
    D103 = Dt103 + nd1*Dn03
    D111 = Dt111 + nd1*Dn11
    D112 = Dt112 + nd1*Dn12
    D113 = Dt113 + nd1*Dn13
    D122 = Dt122 + nd1*Dn22
    D123 = Dt123 + nd1*Dn23
    D133 = Dt133 + nd1*Dn33

    D200 = Dt200 + nd2*Dn00    
    D201 = Dt201 + nd2*Dn01
    D202 = Dt202 + nd2*Dn02
    D203 = Dt203 + nd2*Dn03
    D211 = Dt211 + nd2*Dn11
    D212 = Dt212 + nd2*Dn12
    D213 = Dt213 + nd2*Dn13   
    D222 = Dt222 + nd2*Dn22
    D223 = Dt223 + nd2*Dn23
    D233 = Dt233 + nd2*Dn33

    D300 = Dt300 + nd3*Dn00    
    D301 = Dt301 + nd3*Dn01
    D302 = Dt302 + nd3*Dn02
    D303 = Dt303 + nd3*Dn03  
    D311 = Dt311 + nd3*Dn11
    D312 = Dt312 + nd3*Dn12
    D313 = Dt313 + nd3*Dn13
    D322 = Dt322 + nd3*Dn22
    D323 = Dt323 + nd3*Dn23
    D333 = Dt333 + nd3*Dn33

!-  the primitive scalar --------------
    
    pim    = 0.5d0*(Smp + Smm) + sigma2*phim
    pin    = 0.5d0*(Snp + Snm) + sigma2*phin
    pir    = 0.5d0*(Srp + Srm) + sigma2*phir
    pic    = 0.5d0*(Scp + Scm) + sigma2*phic
    dnphim = 0.5d0*(Smp - Smm)     
    dnphin = 0.5d0*(Snp - Snm)     
    dnphir = 0.5d0*(Srp - Srm)     
    dnphic = 0.5d0*(Scp - Scm)     

    dphim1 = dtphim1 + dnphim*nd1
    dphim2 = dtphim2 + dnphim*nd2
    dphim3 = dtphim3 + dnphim*nd3

    dphin1 = dtphin1 + dnphin*nd1
    dphin2 = dtphin2 + dnphin*nd2
    dphin3 = dtphin3 + dnphin*nd3
    
    dphir1 = dtphir1 + dnphir*nd1
    dphir2 = dtphir2 + dnphir*nd2
    dphir3 = dtphir3 + dnphir*nd3
    
    dphic1 = dtphic1 + dnphic*nd1
    dphic2 = dtphic2 + dnphic*nd2
    dphic3 = dtphic3 + dnphic*nd3

!-- the Frans gauge ---------------------

    G0   = 0.5d0*(Gp + Gm) + sigma2*H0
    dnH0 = 0.5d0*(Gp - Gm)     
    d1H0 = dt1H0 + dnH0*nd1
    d2H0 = dt2H0 + dnH0*nd2
    d3H0 = dt3H0 + dnH0*nd3
        
        
!--------------------------------------------------------------
!--------------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------


  end subroutine


!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
  subroutine what_boundary(direction,n,p,q)
    implicit none
    character(len=*), intent(in):: direction
    CCTK_REAL, dimension(3), intent(out):: n,p,q

    !----- construct a triad adapted to the boundary ---------!
    ! n is the normal, p and q are transverse to the boundary
    
    !^^^^^^^^^^^^^^^^^^^^^  FACES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    choose_boundary : if (direction=="xmax") then
       n=(/1.,0.,0./)
       p=(/0.,1.,0./) 
    else if (direction=="xmin") then
       n=(/-1.,0.,0./)
       p=(/0.,1.,0./) 
    else if (direction=="ymax") then
       n=(/0.,1.,0./)
       p=(/0.,0.,1./) 
    else if (direction=="ymin") then
       n=(/0.,-1.,0./)
       p=(/0.,0.,1./)
    else if (direction=="zmax") then
       n=(/0.,0.,1./) 
       p=(/1.,0.,0./)
    else if (direction=="zmin") then
       n=(/0.,0.,-1./) 
       p=(/1.,0.,0./)
       
       !^^^^^^^^^^^^^^^^^^^^^  EDGES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    else if (direction=="xmin_ymin") then
       n=(/-1./sqrt(2.),-1./sqrt(2.),0./)
       p=(/-1./sqrt(2.),1./sqrt(2.),0./)
    else if (direction=="xmin_ymax") then
       n=(/-1./sqrt(2.),1./sqrt(2.),0./)
       p=(/-1./sqrt(2.),-1./sqrt(2.),0./)
    else if (direction=="xmin_zmin") then
       n=(/-1./sqrt(2.),0.,-1./sqrt(2.)/)
       p=(/-1./sqrt(2.),0.,1./sqrt(2.)/)
    else if (direction=="xmin_zmax") then
       n=(/-1./sqrt(2.),0.,1./sqrt(2.)/)
       p=(/-1./sqrt(2.),0.,-1./sqrt(2.)/)

    else if (direction=="xmax_ymin") then
       n=(/1./sqrt(2.),-1./sqrt(2.),0./)
       p=(/1./sqrt(2.),1./sqrt(2.),0./)
    else if (direction=="xmax_ymax") then
       n=(/1./sqrt(2.),1./sqrt(2.),0./)
       p=(/1./sqrt(2.),-1./sqrt(2.),0./)
    else if (direction=="xmax_zmin") then
       n=(/1./sqrt(2.),0.,-1./sqrt(2.)/)
       p=(/1./sqrt(2.),0.,1./sqrt(2.)/)
    else if (direction=="xmax_zmax") then
       n=(/1./sqrt(2.),0.,1./sqrt(2.)/)
       p=(/1./sqrt(2.),0.,-1./sqrt(2.)/)
       
    else if (direction=="ymin_zmin") then
       n=(/0.,-1./sqrt(2.),-1./sqrt(2.)/)
       p=(/0.,-1./sqrt(2.),1./sqrt(2.)/)
    else if (direction=="ymin_zmax") then
       n=(/0.,-1./sqrt(2.),1./sqrt(2.)/)
       p=(/0.,-1./sqrt(2.),-1./sqrt(2.)/)
    else if (direction=="ymax_zmin") then
       n=(/0.,1./sqrt(2.),-1./sqrt(2.)/)
       p=(/0.,1./sqrt(2.),1./sqrt(2.)/)
    else if (direction=="ymax_zmax") then
       n=(/0.,1./sqrt(2.),1./sqrt(2.)/)
       p=(/0.,1./sqrt(2.),-1./sqrt(2.)/)

       !^^^^^^^^^^^^^^^^^^^^^  CORNERS ^^^^^^^^^^^^^^^^^^^^^^^^^!
     
    else if (direction=="xmin_ymin_zmin") then
       n=(/-1./sqrt(3.),-1./sqrt(3.),-1./sqrt(3.)/)
       p=(/-1./sqrt(2.),1./sqrt(2.),0./)
    else if (direction=="xmin_ymin_zmax") then
       n=(/-1./sqrt(3.),-1./sqrt(3.),1./sqrt(3.)/)
       p=(/-1./sqrt(2.),1./sqrt(2.),0./)
    else if (direction=="xmin_ymax_zmin") then
       n=(/-1./sqrt(3.),1./sqrt(3.),-1./sqrt(3.)/)
       p=(/-1./sqrt(2.),-1./sqrt(2.),0./)
    else if (direction=="xmin_ymax_zmax") then
       n=(/-1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)/)
       p=(/-1./sqrt(2.),-1./sqrt(2.),0./)
    else if (direction=="xmax_ymin_zmin") then
       n=(/1./sqrt(3.),-1./sqrt(3.),-1./sqrt(3.)/)
       p=(/1./sqrt(2.),1./sqrt(2.),0./)
    else if (direction=="xmax_ymin_zmax") then
       n=(/1./sqrt(3.),-1./sqrt(3.),1./sqrt(3.)/)
       p=(/1./sqrt(2.),1./sqrt(2.),0./)
    else if (direction=="xmax_ymax_zmin") then
       n=(/1./sqrt(3.),1./sqrt(3.),-1./sqrt(3.)/)
       p=(/1./sqrt(2.),-1./sqrt(2.),0./)
    else if (direction=="xmax_ymax_zmax") then
       n=(/1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)/)
       p=(/1./sqrt(2.),-1./sqrt(2.),0./)

    else 
       write(*,*) "******** wrong direction at the boundary *************"
       STOP
    end if choose_boundary    

    q(1) = n(2)*p(3) - n(3)*p(2)
    q(2) = n(3)*p(1) - n(1)*p(3)
    q(3) = n(1)*p(2) - n(2)*p(1)


  end subroutine what_boundary

  !#####################################################################################!
               

  end module
