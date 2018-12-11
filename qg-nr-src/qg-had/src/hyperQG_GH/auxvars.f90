!----------------------------------------------------------------------
!
!  $Id: auxvars.f90,v 1.18 2008-07-01 12:06:00 carlos Exp $
!
!----------------------------------------------------------------------

#include "cctk.h"
#include "basic.h"
#include "GH_rhs.h"


#if defined GH_MHD || defined GH_FLOWER
  module gh_auxvars
#else
  module auxvars
#endif


  use params

  implicit none
  CONTAINS

  !--------------------------------------------------------------------
  !
  ! Define the stress energy tensor for 1 complex and
  ! 1 real massless scalar field
  !
  !--------------------------------------------------------------------
  
#if defined GH_MHD || defined GH_FLOWER
  subroutine gh_define_aux_vars(v_pt, w_pt, u_pt, par)
#else
  subroutine define_aux_vars(v_pt, w_pt, u_pt, par)
#endif

    CCTK_REAL, dimension(NU_G), intent(IN)    :: u_pt   
    CCTK_REAL, dimension(NV_G), intent(INOUT) :: v_pt
    CCTK_REAL, dimension(NW)                :: w_pt
    CCTK_REAL, dimension(:)                 :: par
    CCTK_REAL :: x,y,z,r,t
    
    ! Function arrays and time derivatives
    
    CCTK_REAL                             :: alp, alp2, deth, trT
    CCTK_REAL, dimension(I_0:I_3)         :: beta, trG
    CCTK_REAL, dimension(I_0:I_3)         :: Duphim, Duphin, Duphir, Duphic    
    CCTK_REAL, dimension(IS3_11:IS3_33)   :: h, huu    
    CCTK_REAL, dimension(IS4_00:IS4_33)   :: guu
    CCTK_REAL, dimension(IS4_000:IS4_333) :: Gudd    
    CCTK_REAL :: dphir0, dphim0, dphin0, dphic0
    CCTK_REAL :: dphidphim, dphidphin, dphidphir, dphidphic, phi2
    CCTK_REAL :: chi, cf, cm, cp, sf_m, sf_lam, cetaP
    
    CCTK_REAL :: sg00,sg01,sg02,sg03,sg11,sg12,sg13,sg22,sg23,sg33    
    CCTK_REAL :: sK00,sK01,sK02,sK03,sK11,sK12,sK13,sK22,sK23,sK33        
    CCTK_REAL :: sD000,sD001,sD002,sD003,sD011,sD012,sD013,sD022,sD023,sD033    
    CCTK_REAL :: sD100,sD101,sD102,sD103,sD111,sD112,sD113,sD122,sD123,sD133    
    CCTK_REAL :: sD200,sD201,sD202,sD203,sD211,sD212,sD213,sD222,sD223,sD233    
    CCTK_REAL :: sD300,sD301,sD302,sD303,sD311,sD312,sD313,sD322,sD323,sD333                

    
    CCTK_REAL :: r2, dtr, dxr, dyr, dzr, Hs, nult, nulx, nuly, nulz, M, vx, vy
    CCTK_REAL :: dtH, dtnult, dtnulx, dtnuly, dtnulz
    CCTK_REAL :: dxH, dxnult, dxnulx, dxnuly, dxnulz
    CCTK_REAL :: dyH, dynult, dynulx, dynuly, dynulz
    CCTK_REAL :: dzH, dznult, dznulx, dznuly, dznulz            
    
    CCTK_REAL :: s1,s2,s3,s4,s5,cgamma
    
    !! local vars
    logical, parameter               :: ltrace2 = .false.
    CCTK_INT :: evolve_fluxes,gauge_type
    
    logical, parameter :: nancheck = .false.
    ! Save code from lapse getting too small
    !    (not wanted when doing a critical search):
    logical, parameter :: rescue   = .true.
    CCTK_INT           :: myid, proc_return_myid, mm
    logical            :: isanan, nanfound

    myid     = proc_return_myid()
    nanfound = .false.
    if (nancheck) then
       !
       ! Checking for nans on input
       !
       do mm = 1, NU_G
          if ( isanan(u_pt(mm)) ) then
             nanfound = .true.
          end if
       end do
       do mm = 1, NV_G
          if ( isanan(v_pt(mm)) ) then
             nanfound = .true.
          end if
       end do
       do mm = 1, NW
          if ( isanan(w_pt(mm)) ) then
             nanfound = .true.
          end if
       end do
       if (nanfound) then
          write(0,*)myid,']  auxvars: Nan found on  input!!!!!'
          do mm = 1, NU_G
             write(0,*)myid,']  auxvars: u_pt(mm) = ',mm,u_pt(mm)
          end do
          do mm = 1, NV_G
             write(0,*)myid,']  auxvars: v_pt(mm) = ',mm,v_pt(mm)
          end do
          do mm = 1, NW
             write(0,*)myid,']  auxvars: w_pt(mm) = ',mm,w_pt(mm)
          end do
          call my_exit('auxvars: Nan found on input')
       end if
    end if

    x = w_pt(H_XPHYS)
    y = w_pt(H_YPHYS)
    z = w_pt(H_ZPHYS)

    r = sqrt(x*x + y*y + z*z)
    if (r .LT. 1E-4) r=1E-4
    
    sf_m   = par(P_SF_MASS)
    sf_lam = par(P_SF_LAMBDA)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---- compute the 3+ metric components ----------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    h11 = g11
    h12 = g12
    h13 = g13
    h22 = g22
    h23 = g23
    h33 = g33
                                          
    huu11 = -h23**2 + h22*h33
    huu12 = h13*h23 - h12*h33
    huu13 = -(h13*h22) + h12*h23
    huu22 = -h13**2 + h11*h33
    huu23 = h12*h13 - h11*h23
    huu33 = -h12**2 + h11*h22
    deth = h11*huu11 + h12*huu12 + h13*huu13
    if(deth.lt.0.1) deth = 0.1
    huu11 = huu11/deth
    huu12 = huu12/deth
    huu13 = huu13/deth
    huu22 = huu22/deth
    huu23 = huu23/deth
    huu33 = huu33/deth
       
    b1 = g01*huu11 + g02*huu12 + g03*huu13
    b2 = g01*huu12 + g02*huu22 + g03*huu23
    b3 = g01*huu13 + g02*huu23 + g03*huu33
    alp2 = -g00 + b1**2*h11 + b2**2*h22 + 2*(b1*(b2*h12 + b3*h13) &
   &      + b2*b3*h23) + b3**2*h33
    if(rescue .and. alp2.lt.0.01) alp2 = 0.01
    alp = sqrt(alp2)
    
    if (alp2.le.0) then
!          write(0,*)myid,']  auxvars: alp2 <=0, Old alp2: ',alp2
!          write(0,*)myid,']  auxvars: x,y,z: ',x,y,z
!          write(0,*)myid,']  auxvars: g00      = ',g00
!          write(0,*)myid,']  auxvars: dx       = ',par(p_dx)
          !
          alp2 = 0.001d0*par(p_dx)
          alp  = sqrt(alp2)
          !
!          write(0,*)myid,']  auxvars: New alp2 = ',alp2
!          write(0,*)myid,']  auxvars: New alp  = ',alp 
!          write(0,*)myid,']  auxvars: hxyz0    =',par(p_hxyz0)
!          write(0,*)myid,']  auxvars: h        =',par(p_h)
    end if

    !------------ set the values of the lapse and the shift --------
    !-- it could be done in the include, but here is more evident
    !--- for the GH-MHD joining of the codes -----------------------
    v_pt(H_ALPHA)  = alp       
    v_pt(H_SHIFT1) = b1       
    v_pt(H_SHIFT2) = b2       
    v_pt(H_SHIFT3) = b3                           
            
    !------------- invert the metric --------------------------------           
    guu00 = -(1.0d0/alp2)
    guu01 = b1/alp2
    guu02 = b2/alp2
    guu03 = b3/alp2
    guu11 = -(b1**2/alp2) + huu11
    guu12 = -((b1*b2)/alp2) + huu12
    guu13 = -((b1*b3)/alp2) + huu13
    guu22 = -(b2**2/alp2) + huu22
    guu23 = -((b2*b3)/alp2) + huu23
    guu33 = -(b3**2/alp2) + huu33

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---- compute the stress-energy tensor ----------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    phi2 = phic**2 + phir**2 + phim**2 + phin**2
   
    dphim0 = b1*dphim1 + b2*dphim2 + b3*dphim3 - alp*pim
    dphin0 = b1*dphin1 + b2*dphin2 + b3*dphin3 - alp*pin
    dphir0 = b1*dphir1 + b2*dphir2 + b3*dphir3 - alp*pir
    dphic0 = b1*dphic1 + b2*dphic2 + b3*dphic3 - alp*pic
            
    duphim0 = dphim0*guu00 + dphim1*guu01 + dphim2*guu02 + dphim3*guu03
    duphim1 = dphim0*guu01 + dphim1*guu11 + dphim2*guu12 + dphim3*guu13
    duphim2 = dphim0*guu02 + dphim1*guu12 + dphim2*guu22 + dphim3*guu23
    duphim3 = dphim0*guu03 + dphim1*guu13 + dphim2*guu23 + dphim3*guu33

    duphin0 = dphin0*guu00 + dphin1*guu01 + dphin2*guu02 + dphin3*guu03
    duphin1 = dphin0*guu01 + dphin1*guu11 + dphin2*guu12 + dphin3*guu13
    duphin2 = dphin0*guu02 + dphin1*guu12 + dphin2*guu22 + dphin3*guu23
    duphin3 = dphin0*guu03 + dphin1*guu13 + dphin2*guu23 + dphin3*guu33

    duphir0 = dphir0*guu00 + dphir1*guu01 + dphir2*guu02 + dphir3*guu03
    duphir1 = dphir0*guu01 + dphir1*guu11 + dphir2*guu12 + dphir3*guu13
    duphir2 = dphir0*guu02 + dphir1*guu12 + dphir2*guu22 + dphir3*guu23
    duphir3 = dphir0*guu03 + dphir1*guu13 + dphir2*guu23 + dphir3*guu33

    duphic0 = dphic0*guu00 + dphic1*guu01 + dphic2*guu02 + dphic3*guu03
    duphic1 = dphic0*guu01 + dphic1*guu11 + dphic2*guu12 + dphic3*guu13
    duphic2 = dphic0*guu02 + dphic1*guu12 + dphic2*guu22 + dphic3*guu23
    duphic3 = dphic0*guu03 + dphic1*guu13 + dphic2*guu23 + dphic3*guu33

    dphidphim = dphim0*duphim0 + dphim1*duphim1 + dphim2*duphim2 +&
     & dphim3*duphim3+0.d0
    dphidphin = dphin0*duphin0 + dphin1*duphin1 + dphin2*duphin2 +&
     & dphin3*duphin3+0.d0
    dphidphir = dphir0*duphir0 + dphir1*duphir1 + dphir2*duphir2 +&
     & dphir3*duphir3+0.d0
    dphidphic = dphic0*duphic0 + dphic1*duphic1 + dphic2*duphic2 +&
     & dphic3*duphic3+0.d0
    

#ifndef GH_MHD
#ifndef GH_FLOWER 
    T00 = 0.0; T01 = 0.0; T02 = 0.0; T03 = 0.0; T11 = 0.0
    T12 = 0.0; T13 = 0.0; T22 = 0.0; T23 = 0.0; T33 = 0.0
#endif
#endif
        
    T00 = T00 + dphic0**2 + dphir0**2 + dphim0**2 + dphin0**2 -&
     & g00*((phi2**2*sf_lam)/4. + (dphidphic + dphidphim + dphidphir +&
     & dphidphin + phi2*sf_m**2)/2.)+0.d0
    T01 = T01 + dphic0*dphic1 + dphim0*dphim1 + dphin0*dphin1 + &
     & dphir0*dphir1 - g01*((phi2**2*sf_lam)/4. + (dphidphic + &
     & dphidphim + dphidphin + dphidphir + phi2*sf_m**2)/2.)+0.d0
    T02 = T02 + dphic0*dphic2 + dphim0*dphim2 + dphin0*dphin2 + &
     & dphir0*dphir2 - g02*((phi2**2*sf_lam)/4. + (dphidphic + &
     & dphidphim + dphidphin + dphidphir + phi2*sf_m**2)/2.)+0.d0
    T03 = T03 + dphic0*dphic3 + dphim0*dphim3 + dphin0*dphin3 + &
     & dphir0*dphir3 - g03*((phi2**2*sf_lam)/4. + (dphidphic + &
     & dphidphim + dphidphin + dphidphir + phi2*sf_m**2)/2.)+0.d0
    T11 = T11 + dphic1**2 + dphim1**2 + dphin1**2 + dphir1**2 -&
     & g11*((phi2**2*sf_lam)/4. + (dphidphic + dphidphim + &
     & dphidphin + dphidphir + phi2*sf_m**2)/2.)+0.d0
    T12 = T12 + dphic1*dphic2 + dphim1*dphim2 + &
     & dphin1*dphin2 + dphir1*dphir2 -&
     & g12*((phi2**2*sf_lam)/4. +&
     & (dphidphic + dphidphim + dphidphin + dphidphir +&
     & phi2*sf_m**2)/2.)+0.d0
    T13 = T13 + dphic1*dphic3 + dphim1*dphim3 + &
     & dphin1*dphin3 +&
     &  dphir1*dphir3 - g13*((phi2**2*sf_lam)/4. +&
     &  (dphidphic + dphidphim + dphidphin + dphidphir +&
     &  phi2*sf_m**2)/2.)+0.d0
    T22 = T22 + dphic2**2 + dphim2**2 + dphin2**2 + dphir2**2 -&
     & g22*((phi2**2*sf_lam)/4. + (dphidphic + dphidphim + &
     & dphidphin + dphidphir + phi2*sf_m**2)/2.)+0.d0
    T23 = T23 + dphic2*dphic3 + dphim2*dphim3 +&
     & dphin2*dphin3 + &
     & dphir2*dphir3 - g23*((phi2**2*sf_lam)/4. &
     &  + (dphidphic + dphidphim + dphidphin + dphidphir +&
     & phi2*sf_m**2)/2.)+0.d0
    T33 = T33 + dphic3**2 + dphim3**2 + dphin3**2 + dphir3**2 -&
     & g33*((phi2**2*sf_lam)/4. + (dphidphic + dphidphim + &
     & dphidphin + dphidphir + phi2*sf_m**2)/2.)+0.d0
                  

     if (r .LT. 1E-2) then
        print*, "in auxvars"
        print*, "x,y,z",x,y,z
        print*, "T00,T01,T02,T03=", T00,T01,T02,T03
        print*, "T11,T12,T13=", T11,T12,T13
        print*, "T22,T23,T33=", T22,T23,T33

        print*, "g00,g01,g02,g03=", guu00,guu01,guu02,guu03
        print*, "guu11,guu12,guu13=", guu11,guu12,guu13
        print*, "guu22,guu23,guu33=", guu22,guu23,guu33

!    dphic0 = b1*dphic1 + b2*dphic2 + b3*dphic3 - alp*pic

        print*, "alp=", alp
        print*, "pic=", pic/sqrt(2.0)
        print*, "phir=", phir/sqrt(2.0)
        print*, "dphir=", dphir0/sqrt(2.0),dphir1/sqrt(2.0),&
     &           dphir2/sqrt(2.0),dphir3/sqrt(2.0)
        print*, "dphic=", dphic0/sqrt(2.0),dphic1/sqrt(2.0),&
     &           dphic2/sqrt(2.0),dphic3/sqrt(2.0)
        print*, "duphir=", duphir0/sqrt(2.0),duphir1/sqrt(2.0),&
     &           duphir2/sqrt(2.0),duphir3/sqrt(2.0)
        print*, "duphic=", duphic0/sqrt(2.0),duphic1/sqrt(2.0),&
     &           duphic2/sqrt(2.0),duphic3/sqrt(2.0)
     end if


    if (ltrace2) write(0,*)'DEFINE_AUX_VARS: begin'
    
    if (nancheck) then
       !
       ! Checking for nans, but probably do not
       ! want to bother unless problems arise:
       !
       do mm = 1, NU_G
          if ( isanan(u_pt(mm)) ) then
             nanfound = .true.
          end if
       end do
       do mm = 1, NV_G
          if ( isanan(v_pt(mm)) ) then
             nanfound = .true.
          end if
       end do
       do mm = 1, NW
          if ( isanan(w_pt(mm)) ) then
             nanfound = .true.
          end if
       end do
       if (nanfound) then
          write(0,*)myid,']  auxvars: Nan found on output!!!!!'
          write(0,*)myid,']  auxvars: deth = ',deth
          write(0,*)myid,']  auxvars: alp  = ',alp
          write(0,*)myid,']  auxvars: alp2 = ',alp2
          do mm = 1, NU_G
             write(0,*)myid,']  auxvars: u_pt(mm) = ',mm,u_pt(mm)
          end do
          do mm = 1, NV_G
             write(0,*)myid,']  auxvars: v_pt(mm) = ',mm,v_pt(mm)
          end do
          do mm = 1, NW
             write(0,*)myid,']  auxvars: w_pt(mm) = ',mm,w_pt(mm)
          end do
          call my_exit('auxvars: Nan found on output')
       end if
    end if
    
               
    return
  end subroutine

  !--------------------------------------------------------------------
  !
  !
  !
  !--------------------------------------------------------------------
#ifndef GH_MHD
#ifndef GH_FLOWER
  subroutine SolvePrimVars(v, w, u, par)
    use GF

    type(gridfunction), dimension(NV_G) :: v
    type(gridfunction), dimension(NU_G) :: u
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par

  end subroutine SolvePrimVars
#endif
#endif


end module
