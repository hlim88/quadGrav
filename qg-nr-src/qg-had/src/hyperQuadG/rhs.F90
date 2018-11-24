#include "cctk.h"
!---------------------------------------------------------------------------
!
! $Id: rhs.F90,v 1.4 2010-08-03 18:04:48 dneilsen Exp $
!
!---------------------------------------------------------------------------

#ifdef GR_MHD
      module bssn_rhs
#else
      module rhs_mine
#endif
        use GF
        use params
        use OUTPUT
        use HYPER_DERIVS_42
        !   use OUTPUT
        CCTK_REAL, parameter :: r0 = 10.0d0 
        CCTK_REAL, parameter :: a1 = 5.0d0
        CCTK_REAL, parameter :: a2 = 6.0d0
        CCTK_REAL, parameter :: a3 = 7.0d0
        CCTK_REAL, parameter :: a4 = 8.0d0
        CCTK_REAL, parameter :: a5 = 4.0d0

        contains

!----------------------------------------------------------------------
!  This routine is empty
!----------------------------------------------------------------------
#ifdef GR_MHD
      subroutine empty_bssn_calcrhs (RHS, u, v, dxu, dyu, dzu, dxv, dyv, 
     &                         dzv, w, par)
#else
      subroutine calcrhs(RHS, u, v, dxu, dyu, dzu, dxv, dyv, 
     &                         dzv, w, par)
#endif
        implicit   none

        CCTK_REAL, dimension(NU_G)              :: RHS
        CCTK_REAL, dimension(NU_G)              :: u,dxu,dyu,dzu
        CCTK_REAL, dimension(NV_G)              :: dxv,dyv,dzv, v
        CCTK_REAL, dimension(NW)              :: w
        CCTK_REAL, dimension(NPAR)            :: par

!     Nothing goes here for now.
        return
      end subroutine

!----------------------------------------------------------------------
!
!  CALCRHS2
!
!----------------------------------------------------------------------
#ifdef GR_MHD
      subroutine bssn_calcrhs(dtu, u, v, w, imask, par)
#else
      subroutine calcrhs2(dtu, u, v, w, dxu, dyu, dzu, dxv, 
     &                       dyv, dzv, imask, par)
#endif
        implicit none

        include 'mpif.h'

        type(gridfunction), dimension(NU_G)       :: dtu
        type(gridfunction), dimension(NU_G)       :: u
        type(gridfunction), dimension(NV_G)       :: v
        type(gridfunction), dimension(NW)         :: w
        CCTK_INT, dimension(:,:,:)                :: imask
        CCTK_REAL, dimension(NPAR),intent(in)     :: par
#ifndef GR_MHD
        type(gridfunction), dimension(NU_G)       :: dxu, dyu, dzu
        type(gridfunction), dimension(NV_G)       :: dxv, dyv, dzv
#endif

        ! local vars
        CCTK_INT  :: nx, ny, nz, nd, shp(3)
        CCTK_INT  :: i

        CCTK_REAL :: dx, dy, dz
        CCTK_INT  :: d_type, dd_type

        logical, parameter :: ltrace  = .false.
        logical, parameter :: ltrace2 = .false.

        CCTK_REAL, dimension(:,:,:), pointer :: Alpha
        CCTK_REAL, dimension(:,:,:), pointer :: shiftx 
        CCTK_REAL, dimension(:,:,:), pointer :: shifty
        CCTK_REAL, dimension(:,:,:), pointer :: shiftz
        CCTK_REAL, dimension(:,:,:), pointer :: gtxx 
        CCTK_REAL, dimension(:,:,:), pointer :: gtxy 
        CCTK_REAL, dimension(:,:,:), pointer :: gtxz
        CCTK_REAL, dimension(:,:,:), pointer :: gtyy 
        CCTK_REAL, dimension(:,:,:), pointer :: gtyz 
        CCTK_REAL, dimension(:,:,:), pointer :: gtzz
        CCTK_REAL, dimension(:,:,:), pointer :: Atxx 
        CCTK_REAL, dimension(:,:,:), pointer :: Atxy 
        CCTK_REAL, dimension(:,:,:), pointer :: Atxz
        CCTK_REAL, dimension(:,:,:), pointer :: Atyy 
        CCTK_REAL, dimension(:,:,:), pointer :: Atyz 
        CCTK_REAL, dimension(:,:,:), pointer :: Atzz
        CCTK_REAL, dimension(:,:,:), pointer :: chi 
        CCTK_REAL, dimension(:,:,:), pointer :: trK
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtx 
        CCTK_REAL, dimension(:,:,:), pointer :: Gamty 
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtz
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn00 Tmn01 Tmn02
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn03 Tmn11 Tmn12
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn13 Tmn22 Tmn23
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn33
        CCTK_REAL, dimension(:,:,:), pointer :: gxx 
        CCTK_REAL, dimension(:,:,:), pointer :: gxy 
        CCTK_REAL, dimension(:,:,:), pointer :: gxz
        CCTK_REAL, dimension(:,:,:), pointer :: gyy 
        CCTK_REAL, dimension(:,:,:), pointer :: gyz 
        CCTK_REAL, dimension(:,:,:), pointer :: gzz
        CCTK_REAL, dimension(:,:,:), pointer :: gbx 
        CCTK_REAL, dimension(:,:,:), pointer :: gby 
        CCTK_REAL, dimension(:,:,:), pointer :: gbz

        CCTK_REAL, dimension(:,:,:), pointer :: Ex 
        CCTK_REAL, dimension(:,:,:), pointer :: Ey 
        CCTK_REAL, dimension(:,:,:), pointer :: Ez
        CCTK_REAL, dimension(:,:,:), pointer :: Bx 
        CCTK_REAL, dimension(:,:,:), pointer :: By 
        CCTK_REAL, dimension(:,:,:), pointer :: Bz
        CCTK_REAL, dimension(:,:,:), pointer :: Phi_em 
        CCTK_REAL, dimension(:,:,:), pointer :: Psi_em
        CCTK_REAL, dimension(:,:,:), pointer :: phiR
        CCTK_REAL, dimension(:,:,:), pointer :: phiI
        CCTK_REAL, dimension(:,:,:), pointer :: piR
        CCTK_REAL, dimension(:,:,:), pointer :: piI


        CCTK_REAL, dimension(:,:,:), pointer :: x 
        CCTK_REAL, dimension(:,:,:), pointer :: y 
        CCTK_REAL, dimension(:,:,:), pointer :: z
        CCTK_REAL, allocatable, dimension(:) :: x1d
        CCTK_REAL, allocatable, dimension(:) :: y1d
        CCTK_REAL, allocatable, dimension(:) :: z1d


        CCTK_REAL :: time
        CCTK_INT :: myproc, ierr
        CCTK_REAL, dimension(:), allocatable :: tempcoord
        CCTK_INT  :: rank, rc,  j, k
        CCTK_REAL :: nrm1, nrm2, nrm3

        ! some vars for debugging derivatives
        logical, parameter  :: debug_derivs = .false.
        CCTK_INT            :: deriv_error
        CCTK_INT            :: direc
        CCTK_INT            :: bssn_adv_derivs
        CCTK_REAL, save     :: deriv_time = 0.0d0
        CCTK_REAL           :: deriv_dt   = 1.0d-1
        character(len=64)   :: fname
        logical             :: apply_evil_bc = .false.

        


        
!       CCTK_INT, external :: gft_out_full 

        ! memory pointers
#include "MR_declare_ptrs.h"

        include 'mem.inc'
        CCTK_INT  ::  mem_alloc
        CCTK_REAL :: myl2norm3d
        external  :: mem_alloc, myl2norm3d

        if (debug_derivs) then
          ! This is a time coordinate for debugging derivs with SDFs
          deriv_time = deriv_time + deriv_dt
          deriv_error = 0
        end if

        if (ltrace) write(0,*)'*** BEGIN bssn_calcrhs'

        shp     = shape(u(1)%d)
        nx      = shp(1)
        ny      = shp(2)
        nz      = shp(3)
        nd      = nx*ny*nz

        if (nx .ne. nint(par(P_NX))) then
          write(0,*)'*** rhs.f90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'nx        = ',nx
          write(0,*) 'par(P_NX) = ',nint(par(P_NX))
        end if
        if (ny .ne. nint(par(P_NY))) then
          write(0,*)'*** rhs.f90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'ny        = ',ny
          write(0,*) 'par(P_NY) = ',nint(par(P_NY))
        end if
        if (nz .ne. nint(par(P_NZ))) then
          write(0,*)'*** rhs.f90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'nz        = ',nz
          write(0,*) 'par(P_NZ) = ',nint(par(P_NZ))
        end if


        bssn_adv_derivs = nint(par(P_BSSN_ADV_DERIVS))

!       write(0,*) 'AAAAAAAAAAAAAAAAAAAAAAAA'
        if (ltrace2) then
          write(0,*)'*** BSSN_RHS: writing norms. shp=',nx,ny,nz
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: A ||u||',i,
     &                            myl2norm3d(u(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NV_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: A ||v||',i,
     &                            myl2norm3d(v(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NW
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: A ||w||',i,
     &                            myl2norm3d(w(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: A ||dtu||',i,
     &                            myl2norm3d(dtu(i)%d, nx, ny, nz)
!            endif 
             
          end do
        end if

        x => w(H_XPHYS)%d
        y => w(H_YPHYS)%d
        z => w(H_ZPHYS)%d

        allocate(x1d(nx), y1d(ny), z1d(nz))
        x1d(:) = x(:,1,1)
        y1d(:) = y(1,:,1)
        z1d(:) = z(1,1,:)

        time = par(P_TIME)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myproc,ierr)

!       call test_set_funcs(u(G_PHI)%d, u(G_ALPHA)%d,
!    &                      u(G_SHIFT1)%d, u(G_SHIFT2)%d, u(G_SHIFT3)%d,
!    &                      u(G_TRK)%d, x1d, y1d, z1d, nx, ny, nz)

        dx      = par(P_DX)
        dy      = par(P_DY)
        dz      = par(P_DZ)
        d_type  = nint(par(P_DERIV_ORDER))
        dd_type = nint(par(P_DERIV_ORDER))

        Alpha       => u(G_ALPHA)%d
        shiftx      => u(G_SHIFT1)%d
        shifty      => u(G_SHIFT2)%d
        shiftz      => u(G_SHIFT3)%d
        chi         => u(G_CHI)%d
        trK         => u(G_TRK)%d
        gtxx        => u(G_GT11)%d
        gtxy        => u(G_GT12)%d
        gtxz        => u(G_GT13)%d
        gtyy        => u(G_GT22)%d
        gtyz        => u(G_GT23)%d
        gtzz        => u(G_GT33)%d
        Atxx        => u(G_A11)%d
        Atxy        => u(G_A12)%d
        Atxz        => u(G_A13)%d
        Atyy        => u(G_A22)%d
        Atyz        => u(G_A23)%d
        Atzz        => u(G_A33)%d
        Gamtx       => u(G_GAM1)%d
        Gamty       => u(G_GAM2)%d
        Gamtz       => u(G_GAM3)%d
        gbx         => u(G_GB1)%d
        gby         => u(G_GB2)%d
        gbz         => u(G_GB3)%d

        Ex          => u(G_EX)%d
        Ey          => u(G_EY)%d
        Ez          => u(G_EZ)%d
        Bx          => u(G_BX)%d
        By          => u(G_BY)%d
        Bz          => u(G_BZ)%d
        Phi_em      => u(G_PHI_EM)%d
        Psi_em      => u(G_PSI_EM)%d

        phiR        => u(G_PHIR)%d
        phiI        => u(G_PHII)%d
        piR         => u(G_PIR)%d
        piI         => u(G_PII)%d


        if (ltrace) write(0,*)'BSSN_RHS: allocating memory'
#include "MR_alloc.h"

        if (ltrace) write(0,*)'BSSN_RHS: calling derivatives'
        deriv_error = 0
#include "MR_derivs.h"
!       if (deriv_error .eq. 1) then
!         call my_exit('Die on deriv_error')
!       end if

#if 0

        if (ltrace2) then
          ! This is low level debugging stuff.  When the code is production
          ! ready, this should just be removed.

         call output3D(u(33)%d,"xphiR",x1d,y1d,z1d,time,myproc)

         call output1d_as_3D(q(dx_phiR),'xdxphiR',x1d,y1d, z1d,
     &                           nx, ny, nz, nd, time,myproc)
         call output1d_as_3D(q(dy_phiR),'xdyphiR',x1d,y1d, z1d,
     &                           nx, ny, nz, nd, time,myproc)
         call output1d_as_3D(q(dz_phiR),'xdzphiR',x1d,y1d, z1d,
     &                           nx, ny, nz, nd, time,myproc)



!         call output3D(u(18)%d,"xalpha",x1d,y1d,z1d,time,myproc)
!         call output3D(u(19)%d,"xbetax",x1d,y1d,z1d,time,myproc)
!         call output3D(u(20)%d,"xbetay",x1d,y1d,z1d,time,myproc)
!         call output3D(u(21)%d,"xbetaz",x1d,y1d,z1d,time,myproc)
  
!         call output1d_as_3D(q(dx_Atxx),'x_dx_Atxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)
!         call output1d_as_3D(q(dy_Atxx),'x_dy_Atxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)
!         call output1d_as_3D(q(dz_Atxx),'x_dz_Atxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)

!         call output1d_as_3D(q(dxx_gtxx),'x_dxx_gtxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)
!         call output1d_as_3D(q(dxy_gtxx),'x_dxy_gtxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)
!         call output1d_as_3D(q(dxz_gtxx),'x_dxz_gtxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)
!         call output1d_as_3D(q(dyy_gtxx),'x_dyy_gtxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)
!         call output1d_as_3D(q(dyz_gtxx),'x_dyz_gtxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)
!         call output1d_as_3D(q(dzz_gtxx),'x_dzz_gtxx',x1d,y1d, z1d,
!    &                           nx, ny, nz, nd, time,myproc)

!          call output3D(dtu(25)%d,"xxEx",x1d,y1d,z1d,deriv_time,myproc)
!          call output3D(dtu(26)%d,"xxEy",x1d,y1d,z1d,deriv_time,myproc)
!          call output3D(dtu(27)%d,"xxEz",x1d,y1d,z1d,deriv_time,myproc)
        end if
#endif


        if (ltrace) write(0,*)'BSSN_RHS: calling cal_bssn_rhs'

        call cal_bssn_rhs( dtu, u, v, w, imask,
#include "MR_call_bssn.h"
     &                     x1d, y1d, z1d, nx, ny, nz, par)

        
!       call output3D(dtu(13)%d,"xdotphi_b1",x1d,y1d,z1d,time,myproc)
        if (ltrace2) then
          write(0,*)'*** BSSN_RHS: writing norms. shp=',nx,ny,nz
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: B ||u||',i,
     &                            myl2norm3d(u(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NV_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: B ||v||',i,
     &                            myl2norm3d(v(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NW
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: B ||w||',i,
     &                            myl2norm3d(w(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: B ||dtu||',i,
     &                            myl2norm3d(dtu(i)%d, nx, ny, nz)
!            endif 
             
          end do
        end if


        if (ltrace) write(0,*)'BSSN_RHS: calling cal_bssn_bcs_rhs'
        call calc_bssn_bcs_rhs(dtu, u, v, w, imask,
#include "MR_call_bcs.h"
     &                     x1d, y1d, z1d, nx, ny, nz, par)

!       call output3D(dtu(1)%d,"xdotgtxx_b4",x1d,y1d,z1d,time,myproc)
!       call output3D(dtu(7)%d,"xdotaxx_b4",x1d,y1d,z1d,time,myproc)
!       call output3D(dtu(13)%d,"xdotphi_b4",x1d,y1d,z1d,time,myproc)
        if (ltrace2) then
          call output3D(dtu(25)%d,"yyEx",x1d,y1d,z1d,deriv_time,myproc)
          call output3D(dtu(26)%d,"yyEy",x1d,y1d,z1d,deriv_time,myproc)
          call output3D(dtu(27)%d,"yyEz",x1d,y1d,z1d,deriv_time,myproc)
          write(0,*)'*** BSSN_RHS: writing norms. shp=',nx,ny,nz
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: C ||u||',i,
     &                            myl2norm3d(u(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NV_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: C ||v||',i,
     &                            myl2norm3d(v(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NW
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: C ||w||',i,
     &                            myl2norm3d(w(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: C ||dtu||',i,
     &                            myl2norm3d(dtu(i)%d, nx, ny, nz)
!            endif 
             
          end do
        end if

        if (bssn_adv_derivs .eq. 0) then
          if (ltrace) write(0,*)'BSSN_RHS: calling bssn_rhs_lie'
            call cal_bssn_rhs_lie( dtu, u, v, w, imask,
#include "MR_call_bssn_adv_centered.h"
     &                             nx, ny, nz, par)
          if (ltrace) write(0,*)'BSSN_RHS: return bssn_rhs_lie'
        end if


        if (ltrace) write(0,*)'BSSN_RHS: deallocating I'
#include "MR_dealloc.h"

        if (bssn_adv_derivs .eq. 1) then
          if (ltrace) write(0,*)'BSSN_RHS: allocating advective derivs'
#include "MR_alloc_adv.h"

          if (ltrace) write(0,*)'BSSN_RHS: calling advective derivs'
#include "MR_derivs_adv.h"

          if (ltrace) write(0,*)'BSSN_RHS: calling bssn_rhs_lie'
          call cal_bssn_rhs_lie( dtu, u, v, w, imask,
#include "MR_call_bssn_adv.h"
     &                     nx, ny, nz, par)

          if (ltrace) write(0,*)'BSSN_RHS: deallocating last time'
#include "MR_dealloc_adv.h"
        end if

        if (ltrace2) then
          write(0,*)'*** BSSN_RHS: writing norms. shp=',nx,ny,nz
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: D ||u||',i,
     &                            myl2norm3d(u(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NV_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: D ||v||',i,
     &                            myl2norm3d(v(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NW
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: D ||w||',i,
     &                            myl2norm3d(w(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: D ||dtu||',i,
     &                            myl2norm3d(dtu(i)%d, nx, ny, nz)
!            endif 
             
          end do
        end if


!       call test_correct_rhs(dtu(G_PHI)%d, dtu(G_ALPHA)%d,
!    &                            x1d, y1d, z1d, nx, ny, nz)

        if (apply_evil_bc) then
          ! This boundary condition puts NaNs in places ghostzones and
          ! MPI decomposition zones. It is used to test for problems.
          print *,'calling evil BC'
          call bssn_evil_bcs( dtu, u, v, w, imask,
     &                        x1d, y1d, z1d, nx, ny, nz, par )
        end if

#if 1
        if (ltrace2) then
!         call output3D(dtu(13)%d,"xdotphi",x1d,y1d,z1d,time,myproc)
!         call output3D(dtu(18)%d,"xdotalpha",x1d,y1d,z1d,time,myproc)

!         call output3D(dtu(1)%d,"xdotgtxx",x1d,y1d,z1d,time,myproc)
!         call output3D(dtu(2)%d,"xdotgtxy",x1d,y1d,z1d,time,myproc)
!         call output3D(dtu(7)%d,"xdotaxx",x1d,y1d,z1d,time,myproc)

          call output3D(dtu(25)%d,"zzEx",x1d,y1d,z1d,deriv_time,myproc)
          call output3D(dtu(26)%d,"zzEy",x1d,y1d,z1d,deriv_time,myproc)
          call output3D(dtu(27)%d,"zzEz",x1d,y1d,z1d,deriv_time,myproc)
        end if
#endif

        if (ltrace2) then
          write(0,*)'*** BSSN_RHS: writing norms. shp=',nx,ny,nz
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: END ||u||',i,
     &                            myl2norm3d(u(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NV_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: END ||v||',i,
     &                            myl2norm3d(v(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NW
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: END ||w||',i,
     &                            myl2norm3d(w(i)%d, nx, ny, nz)
!            endif 
             
          end do
          do i = 1, NU_G
            
!            if (i.eq.29) then 
             write(0,*) 'BSSN_RHS: END ||dtu||',i,
     &                            myl2norm3d(dtu(i)%d, nx, ny, nz)
!            endif 
             
          end do
        end if



        if (ltrace) write(0,*)'*** END bssn_calcrhs'

        deallocate(x1d, y1d, z1d)
        
        return
      end subroutine

      !-----------------------------------------------------------------
      !
      !
      !
      !-----------------------------------------------------------------
      subroutine test_set_funcs(chi, alpha, betax, betay, betaz,
     &                          trK, x, y, z, nx, ny, nz)
        implicit none
        CCTK_INT :: nx, ny, nz
        CCTK_REAL, dimension(nx,ny,nz) :: chi, alpha, trK
        CCTK_REAL, dimension(nx,ny,nz) :: betax, betay, betaz
        CCTK_REAL :: x(nx), y(ny), z(nz)

        ! local vars
        CCTK_INT  :: i, j, k
        CCTK_REAL :: r

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx

          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
          chi(i,j,k)    = exp(-(r/r0)**2)*sin(r/a1)
          alpha(i,j,k)  = r**4*exp(-(r/r0)**2)*cos(r/a2)
          betax(i,j,k) = r**2*exp(-(r/r0)**2)*sin(x(i)/a3)
          betay(i,j,k) = r**4*exp(-(r/r0)**2)*cos(y(j)/a4)
          betaz(i,j,k) = r**2*exp(-(r/r0)**2)*sin(z(k)/a5)
          trK(i,j,k)   = 0.0d0
        end do
        end do
        end do

        return
      end subroutine

      !-----------------------------------------------------------------
      !
      !
      !
      !-----------------------------------------------------------------
      subroutine test_correct_rhs(dtchi, dtalpha,
     &                            x1d, y1d, z1d, nx, ny, nz)
        implicit none
        CCTK_INT :: nx, ny, nz
        CCTK_REAL, dimension(nx,ny,nz) :: dtchi, dtalpha
        CCTK_REAL :: x1d(nx), y1d(ny), z1d(nz)

        ! local vars
        CCTK_INT  :: i, j, k
        CCTK_REAL :: r, x, y, z
        CCTK_REAL :: alpharhs, chirhs

        real*8  t1, t2, t3, t4, t5
        real*8  t6, t8, t9, t10, t11
        real*8  t12, s, t14, t15, t16
        real*8  t17, t18, t21, t22, t24
        real*8  t30, t31, t32, t33, t34
        real*8  t43, t44, t45, t58, t62
        real*8  t73, t84, t13, t20, t23
        real*8  t27, t28, t29, t38, t52

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          x = x1d(i)
          y = y1d(j)
          z = z1d(k)
          r = sqrt(x**2 + y**2 + z**2)
          if (abs(r) .gt. 1.0d-10) then
      t1 = x**2
      t2 = y**2
      t3 = z**2
      t4 = t1+t2+t3
      t5 = r0**2
      t6 = 1/t5
      t8 = exp(-t4*t6)
      t9 = t4*t8
      t10 = 1/a3
      t11 = x*t10
      t12 = sin(t11)
      t14 = sqrt(t4)
      t15 = 1/a1
      t16 = t14*t15
      t17 = sin(t16)
      t18 = t8*t17
      t21 = cos(t16)
      t22 = t8*t21
      t24 = 1/t14*t15
      t30 = t4**2
      t31 = t30*t8
      t32 = 1/a4
      t33 = y*t32
      t34 = cos(t33)
      t43 = 1/a5
      t44 = z*t43
      t45 = sin(t44)
      t58 = t6*t8
      t62 = cos(t11)
      t73 = sin(t33)
      t84 = cos(t44)
      chirhs = t9*t12*(-2*x*t6*t18+t22*t24*x)+t31*t34*(-2*y*t6*t18+t22*t
     #24*y)+t9*t45*(-2*z*t6*t18+t22*t24*z)+x*t8*t12/3-t4*x*t58*t12/3+t9*
     #t62*t10/6+2.D0/3.D0*t9*t34*y-t30*y*t58*t34/3-t31*t73*t32/6+z*t8*t4
     #5/3-t4*z*t58*t45/3+t9*t84*t43/6
      t1 = x**2
      t2 = y**2
      t3 = z**2
      t4 = t1+t2+t3
      t5 = r0**2
      t6 = 1/t5
      t8 = exp(-t4*t6)
      t9 = t4*t8
      t12 = sin(x/a3)
      t13 = sqrt(t4)
      t14 = 1/a2
      t15 = t13*t14
      t16 = cos(t15)
      t20 = t4**2
      t23 = t6*t8*t16
      t27 = t13*t4*t8
      t28 = sin(t15)
      t29 = t28*t14
      t38 = cos(y/a4)
      t52 = sin(z/a5)
      alpharhs = t9*t12*(4*t9*t16*x-2*t20*x*t23-t27*t29*x)+t20*t8*t38*(4
     #*t9*t16*y-2*t20*y*t23-t27*t29*y)+t9*t52*(4*t9*t16*z-2*t20*z*t23-t2
     #7*t29*z)
          else 
            alpharhs = 0.0d0
            chirhs   = 0.0d0
          end if
          dtchi(i,j,k)   = dtchi(i,j,k)   - chirhs
          dtalpha(i,j,k) = dtalpha(i,j,k) - alpharhs
        end do
        end do
        end do

        do k = 1, nz
        do j = 1, ny
        do i = 1, 3
          dtalpha(i,j,k) = 0.0d0
          dtchi(i,j,k) = 0.0d0
          dtalpha(nx-i+1,j,k) = 0.0d0
          dtchi(nx-i+1,j,k) = 0.0d0
        end do
        end do
        end do

        do k = 1, nz
        do i = 1, nx
        do j = 1, 3
          dtalpha(i,j,k) = 0.0d0
          dtchi(i,j,k) = 0.0d0
          dtalpha(i,ny-j+1,k) = 0.0d0
          dtchi(i,ny-j+1,k) = 0.0d0
        end do
        end do
        end do

        do j = 1, ny
        do i = 1, nx
        do k = 1, 3
          dtalpha(i,j,k) = 0.0d0
          dtchi(i,j,k) = 0.0d0
          dtalpha(i,j,nz-k+1) = 0.0d0
          dtchi(i,j,nz-k+1) = 0.0d0
        end do
        end do
        end do




        return
      end subroutine

      !-----------------------------------------------------------------
      !
      !
      !
      !-----------------------------------------------------------------
      subroutine  output1d_as_3D(f,fname,x1d,y1d, z1d,
     &                           nx, ny, nz, nd, time,myproc)

        CCTK_INT                                 :: nx, ny, nz, nd
        CCTK_REAL, dimension(nx)                 :: x1d
        CCTK_REAL, dimension(ny)                 :: y1d
        CCTK_REAL, dimension(nz)                 :: z1d
        CCTK_REAL, dimension(nd)                 :: f
        CCTK_REAL                                :: time
        CCTK_INT                                 :: myproc
        character(len=*)                         :: fname

        ! local vars
        CCTK_REAL, dimension(:,:,:), allocatable :: tempf
        CCTK_INT  :: i, j, k


        allocate(tempf(nx,ny,nz))

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          tempf(i,j,k) = f(i + nx*(j-1) + nx*ny*(k-1))
        end do
        end do
        end do

        call output3D(tempf,fname,x1d,y1d,z1d,time,myproc)

        deallocate(tempf) 
        return
      end subroutine

      end module

