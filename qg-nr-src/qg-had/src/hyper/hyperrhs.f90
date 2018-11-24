!-----------------------------------------------------------------
!
!  $Id$
!
!-----------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_RHS
    use params
    use GF
    use rhs_mine
    use HYPER_DERIVS
    use HYPER_AUXVARS
  contains

!-----------------------------------------------------------------
!
! SUBROUTINE hyperrhs
!
!-----------------------------------------------------------------
    subroutine hyperrhs(u_dot, u, v, w, dxu, dyu, dzu, dxv, dyv, dzv, &
         imask, par)

      implicit none

      type(gridfunction), dimension(NU)   :: dxu, dyu, dzu
      type(gridfunction), dimension(NV)   :: dxv, dyv, dzv
      type(gridfunction), dimension(NU)   :: u, u_dot
      type(gridfunction), dimension(NW)   :: w
      type(gridfunction), dimension(NV)   :: v
      CCTK_REAL, dimension(:)             :: par
      CCTK_INT, dimension(:,:,:)          :: imask

      CCTK_REAL, dimension(NU)            :: u_dot_pt, u_pt, dxu_pt, dyu_pt, dzu_pt
      CCTK_REAL, dimension(NV)            :: dxv_pt, dyv_pt, dzv_pt, v_pt
      CCTK_REAL, dimension(NW)            :: w_pt
      CCTK_INT                            :: ii,jj,kk, nx, ny, nz, m
      CCTK_INT                            :: use_mask
      CCTK_REAL                           :: err, myl2norm3d
      CCTK_REAL, pointer, dimension(:,:,:):: mask, x, y, z

      logical, parameter                  :: ltrace  = .false.
      !logical                             :: ltrace
      logical, parameter                  :: ltrace2 = .false.
      CCTK_INT                            :: myid, proc_return_myid
      logical                             :: double_equal
      external                               double_equal
      
      mask => w(H_MASK)%d
      x    => w(H_XPHYS)%d
      y    => w(H_YPHYS)%d
      z    => w(H_ZPHYS)%d


      nx = nint(par(P_NX))
      ny = nint(par(P_NY))
      nz = nint(par(P_NZ))

      use_mask = nint(par(P_USE_MASK))
      myid     = proc_return_myid()

      !ltrace = .false.
      !if (myid.eq.1.and.nx.eq.39.and.ny.eq.41.and.nz.eq.39.and.double_equal(x(1,1,1),-5.0).and.double_equal(y(1,1,1),-1.66666666666667)) ltrace=.true.


      !---------------------------------------------------------------------
      ! We first calculate the auxiliary variables.
      !---------------------------------------------------------------------
      if (ltrace) then
         write(0,*)myid,'] >>> HYPERRHS: calling hyperauxvars'
      end if

      call hyperauxvars(v, w, u, par)      
      
      !---------------------------------------------------------------------
      ! Check to see if derivatives are stored in 3-D arrays.
      !---------------------------------------------------------------------
      
      if (P_POINT_WISE_DERIVATIVES .eq. 0 .and. P_HYPER_DERIVS .eq. 1) then
         if (ltrace) then
            write(0,*)myid,'] >>> HYPERRHS: calling hyperderivs'
         end if
      if (ltrace) then
         !write(0,*)myid, '] >>> HYPERRHS 0.45 u_dot(H_B1)=',u_dot(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45     u(H_B1)=',    u(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45   dxu(H_B1)=',  dxu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45   dyu(H_B1)=',  dyu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45   dzu(H_B1)=',  dzu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45 u_dot(H_E3)=',u_dot(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45     u(H_E3)=',    u(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45   dxu(H_E3)=',  dxu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45   dyu(H_E3)=',  dyu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.45   dzu(H_E3)=',  dzu(H_E3)%d(31,32,14)
         write(0,*)myid, '] >>> HYPERRHS 0.45   calling hyperderivs now'
      end if
         call hyperderivs(dxu,   dyu,  dzu,  dxv,  dyv,  dzv, &
              u,     v,    w,    imask,      par)
      end if
      
      if (ltrace) then
         !write(0,*)myid, '] >>> HYPERRHS 0.5 u_dot(H_B1)%d(31,32,14)=',u_dot(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5     u(H_B1)%d(31,32,14)=',    u(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5   dxu(H_B1)%d(31,32,14)=',  dxu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5   dyu(H_B1)%d(31,32,14)=',  dyu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5   dzu(H_B1)%d(31,32,14)=',  dzu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5 u_dot(H_E3)=',u_dot(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5     u(H_E3)=',    u(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5   dxu(H_E3)=',  dxu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5   dyu(H_E3)=',  dyu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 0.5   dzu(H_E3)=',  dzu(H_E3)%d(31,32,14)
      end if

       if (ltrace2) write(0,*)myid,'] >>> HYPERRHS: not sure what this does...'
       if (ltrace2) write(0,*)myid,'] >>> HYPERRHS: Apparently doing things pointwise'
 
       if (P_HYPER_RHS1 .eq. 1) then
         if (ltrace) write(0,*)myid,'] >>> HYPERRHS: doing things pointwise'
!$OMP PARALLEL DO default(shared) private(kk,jj,ii,m)
       do kk = 1, nz
          do jj = 1, ny
             do ii = 1, nx
 
 
                !-----------------------------------------------------------
                ! Check to see if we are in an excision region
                !-----------------------------------------------------------
                if (use_mask .eq. 1 .and. mask(ii,jj,kk) < 0.0) then
                   do m=1, NU
                      u_dot(m)%d(ii,jj,kk) = 0.0
                   end do
                else
 
                   do m = 1, NU
                      u_pt(m) = u(m)%d(ii,jj,kk)
                   end do
                   do m = 1, NU
                      u_dot_pt(m) = u_dot(m)%d(ii,jj,kk)
                   end do
                   do m = 1, P_NV_REAL_SIZE
                      v_pt(m) = v(m)%d(ii,jj,kk)
                   end do
                   do m = 1, NW
                      w_pt(m) = w(m)%d(ii,jj,kk)
                   end do
 
                   if (P_POINT_WISE_DERIVATIVES .eq. 1) then
                   if (ltrace .and. ii.eq.31.and.jj.eq.32.and.kk.eq.14) then
                      !write(0,*)myid, '] >>> HYPERRHS 0.61 u_dot(H_B1)=',u_dot_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61     u(H_B1)=',    u_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61   dxu(H_B1)=',  dxu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61   dyu(H_B1)=',  dyu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61   dzu(H_B1)=',  dzu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61 u_dot(H_E3)=',u_dot_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61     u(H_E3)=',    u_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61   dxu(H_E3)=',  dxu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61   dyu(H_E3)=',  dyu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.61   dzu(H_E3)=',  dzu_pt(H_E3)
                   end if
                      call hyperderivs_pt(dxu_pt,   dyu_pt,  dzu_pt,        &
                           dxv_pt,   dyv_pt,  dzv_pt,        &
                           u,   v,   w_pt,    imask(ii,jj,kk), & 
                           ii,  jj,  kk,      par)
                   if (ltrace .and. ii.eq.31.and.jj.eq.32.and.kk.eq.14) then
                      !write(0,*)myid, '] >>> HYPERRHS 0.63 u_dot(H_B1)=',u_dot_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63     u(H_B1)=',    u_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63   dxu(H_B1)=',  dxu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63   dyu(H_B1)=',  dyu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63   dzu(H_B1)=',  dzu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63 u_dot(H_E3)=',u_dot_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63     u(H_E3)=',    u_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63   dxu(H_E3)=',  dxu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63   dyu(H_E3)=',  dyu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.63   dzu(H_E3)=',  dzu_pt(H_E3)
                   end if
 
                   else
                      do m = 1, NU
                         u_pt(m)   = u(m)%d(ii,jj,kk)
 
                         if (u(m)%take_dx == 1) then
                            dxu_pt(m) = dxu(m)%d(ii,jj,kk)
                         else
                            dxu_pt(m) = 0.0
                         end if
 
                         if (u(m)%take_dy == 1) then
                            dyu_pt(m) = dyu(m)%d(ii,jj,kk)
                         else
                            dyu_pt(m) = 0.0
                         end if
 
                         if (u(m)%take_dz == 1) then
                            dzu_pt(m) = dzu(m)%d(ii,jj,kk)
                         else
                            dzu_pt(m) = 0.0
                         end if
                      end do

                      do m = 1, P_NV_REAL_SIZE
                          if (v(m)%take_dx == 1) then
                            dxv_pt(m) = dxv(m)%d(ii,jj,kk)
                         else
                            dxv_pt(m) = 0.0
                         end if
                         if (v(m)%take_dy == 1) then
                            dyv_pt(m) = dyv(m)%d(ii,jj,kk)
                         else
                            dyv_pt(m) = 0.0
                         end if
 
                         if (v(m)%take_dz == 1) then
                            dzv_pt(m) = dzv(m)%d(ii,jj,kk)
                         else
                            dzv_pt(m) = 0.0
                         end if
                      end do
                   end if
                   
                   if (ltrace2) write(0,*)myid,'] >>> HYPERRHS: calcrhs on pt:',ii,jj,kk
                   if (ltrace .and. ii.eq.31.and.jj.eq.32.and.kk.eq.14) then
                      !write(0,*)myid, '] >>> HYPERRHS 0.6 u_dot(H_B1)=',u_dot_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6     u(H_B1)=',    u_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6   dxu(H_B1)=',  dxu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6   dyu(H_B1)=',  dyu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6   dzu(H_B1)=',  dzu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6 u_dot(H_E3)=',u_dot_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6     u(H_E3)=',    u_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6   dxu(H_E3)=',  dxu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6   dyu(H_E3)=',  dyu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.6   dzu(H_E3)=',  dzu_pt(H_E3)
                   end if
 
                   call calcrhs(u_dot_pt, u_pt, v_pt, dxu_pt, dyu_pt, dzu_pt,  &
                        dxv_pt, dyv_pt, dzv_pt, w_pt, par)
                    
                   if (ltrace .and. ii.eq.31.and.jj.eq.32.and.kk.eq.14) then
                      !write(0,*)myid, '] >>> HYPERRHS 0.7 u_dot(H_B1)=',u_dot_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7     u(H_B1)=',    u_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7   dxu(H_B1)=',  dxu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7   dyu(H_B1)=',  dyu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7   dzu(H_B1)=',  dzu_pt(H_B1)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7 u_dot(H_E3)=',u_dot_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7     u(H_E3)=',    u_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7   dxu(H_E3)=',  dxu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7   dyu(H_E3)=',  dyu_pt(H_E3)
                      !write(0,*)myid, '] >>> HYPERRHS 0.7   dzu(H_E3)=',  dzu_pt(H_E3)
                   end if
                   do m=1, NU
                      u_dot(m)%d(ii,jj,kk) = u_dot_pt(m)
                   end do
 
                end if
 
             end do
          end do
       end do
!$OMP END PARALLEL DO

       else if (P_HYPER_RHS2 .eq. 1) then
          if (ltrace) write(0,*)myid,'] >>> HYPERRHS: calling calcrhs2'
          call calcrhs2(u_dot, u, v, w, dxu, dyu, dzu, dxv, dyv, dzv, &
               imask, par)
       end if

      if (ltrace) then
         !write(0,*)myid, '] >>> HYPERRHS 1.0 u_dot(H_B1)%d(31,32,14)=',u_dot(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 1.0     u(H_B1)%d(31,32,14)=',    u(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 1.0   dxu(H_B1)%d(31,32,14)=',  dxu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 1.0   dyu(H_B1)%d(31,32,14)=',  dyu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERRHS 1.0   dzu(H_B1)%d(31,32,14)=',  dzu(H_B1)%d(31,32,14)
      end if

    end subroutine hyperrhs
end MODULE HYPER_RHS
