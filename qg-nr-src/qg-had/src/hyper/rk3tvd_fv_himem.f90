!-------------------------------------------------------------------------
!
!  $Id$
!
!
!
!-------------------------------------------------------------------------

#include "cctk.h"


subroutine TVDRK3_FV_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
     urk2, urk3, par)
  use params
  use UTILEQS
  use GF
  use HYPER_RHS
  use HYPER_DISSIPATION
  use HYPER_BOUNDARY
  use MOD_MASK
  use hypercoords


  implicit none

  CCTK_INT     istat, rc, nt, nx, ny, nz, err, m, use_mask, dissipation
  CCTK_INT     update_scheme, amrbound_prepost
  CCTK_INT     boundary_conditions
  type(gridfunction),     dimension(NU)    :: u0,   u2,         &
                                              urk1, urk2, urk3, &
                                              dxu,  dyu,  dzu

  type(gridfunction),     dimension(NV)    :: v, dxv, dyv, dzv
  type(gridfunction),     dimension(NW)    :: w
  CCTK_REAL,              dimension(NPAR)  :: par


  CCTK_INT  :: iter
  CCTK_REAL                                :: dt, dx, sigma, sigma_diss, normv
  CCTK_INT,  allocatable, dimension(:,:,:) :: imask, imask2
  CCTK_REAL, allocatable, dimension(:,:,:) :: dxu_diss, dyu_diss, dzu_diss
  logical, parameter                       :: ltrace  = .false.
  logical, parameter                       :: ltrace2 = .false.
  logical, parameter                       :: ltraceB = .false. ! Trace boundary
  logical, parameter                       :: no_bc   = .false. ! For debugging ONLY, do NOTHING at boundary.
  CCTK_INT CCTK_Equals
  CCTK_INT                                 :: handle, err_red, i, j, k, kk, order
  CCTK_REAL                                :: norm_rk, max_g11, min_g11
  CCTK_REAL, pointer,     dimension(:,:,:) :: chr
  CCTK_REAL, pointer,     dimension(:,:,:) :: xx, yy, zz
  CCTK_INT                                 :: bc_type

  CCTK_REAL                                :: myl2norm3d
  CCTK_REAL, parameter                     :: one_third = 0.33333333333333333

!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================

      iter = nint(par(P_RK_ITER))
      chr => w(H_CHR)%d
      xx  => w(H_XPHYS)%d
      yy  => w(H_YPHYS)%d
      zz  => w(H_ZPHYS)%d
      sigma =0. ! move this later 

      if (ltrace) then
         write(0,*)'>>>  RK3 begins'
      end if
      nx                    = nint(par(P_NX))
      ny                    = nint(par(P_NY))
      nz                    = nint(par(P_NZ))
      update_scheme         = nint(par(P_update_scheme))
      amrbound_prepost      = nint(par(P_amrbound_prepost))
      use_mask              = nint(par(P_USE_MASK))
      order                 = nint(par(P_DERIV_ORDER))  
      dissipation           = nint(par(P_DISSIPATION))
      sigma_diss            = par(P_SIGMA_DISS)
      boundary_conditions   = nint(par(P_BOUNDARY_CONDITIONS))
      bc_type               = nint(par(P_BC_TYPE))

      allocate(dzu_diss(nx,ny,nz), STAT=istat)
      if (istat .ne. 0) then
        write(0,*)'>>> RK3TVD_FV_himem: Could not allocate memory for dzu_diss'
        write(0,*)'               size : ',nx,ny,nz
      end if
      allocate(dyu_diss(nx,ny,nz), STAT=istat)
      if (istat .ne. 0) then
        write(0,*)'>>> RK3TVD_FV_himem: Could not allocate memory for dyu_diss'
        write(0,*)'               size : ',nx,ny,nz
      end if
      allocate(dxu_diss(nx,ny,nz), STAT=istat)
      if (istat .ne. 0) then
        write(0,*)'>>> RK3TVD_FV_himem: Could not allocate memory for dxu_diss'
        write(0,*)'                size : ',nx,ny,nz
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Define coordinates if not using the default cartesian computational 
!   coords.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (P_ENABLE_ALT_COORDS==1) then
        if (ltrace) then
           write(0,*)'>>> RK3TVD_FV_himem: Calling def_alt_coords'
        end if
        call def_alt_coords(w, par)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This is not (always) the right value,
!   it is the global dx:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dx = par(P_DX)
      dt = par(P_DT)

      if (use_mask .eq. 1) then
        allocate(imask(nx,ny,nz),STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'RK3TVD_FV_himem: Could not allocate memory for imask array'
          write(0,*)'           size : ',nx,ny,nz
        end if

        allocate(imask2(nx,ny,nz),STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'RK3TVD_FV_himem: Could not allocate memory for imask2 array'
          write(0,*)'           size : ',nx,ny,nz
        end if
        call set_mask(w, imask,imask2, par)
      end if

      if ( ltrace2 ) then
          write(0,*)'>>> RK3TVD_FV_himem----params:'
          write(0,*)'       iter        = ',iter
          write(0,*)'       nx, ny, nz  = ',nx, ny, nz
          write(0,*)'       use_mask    = ',use_mask
          write(0,*)'       dissipation = ',dissipation
          write(0,*)'       boundary_conditions = ',boundary_conditions
          write(0,*)'       bc_type = ',bc_type     
          write(0,*)'       update_scheme = ', update_scheme
          write(0,*)'>>> RK3TVD_FV_himem----end params:'
       end if

!#######################################################################
!####################### 1st RK step ###################################
!#######################################################################

    if (iter .eq. 1) then
       if (ltrace) then
           write(0,*)'>>> RK3TVD_FV_himem: begin first step'
           write(0,*)'>>> RK3TVD_FV_himem: dt,  dx = ',dt, dx
       end if
       if ( ltrace2 ) then
           write(0,*)'>>>RK3TVD_FV_himem: initial data'
           do m=1,NU
              write(0,*) '||u||  ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
           end do
           do m=1,NV
              write(0,*) '||v||  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
           end do
       end if

       if ( boundary_conditions .eq. 1 .or.             &
            boundary_conditions .eq. 2 .and. .not.no_bc &
          ) then
            if ( ltraceB ) write(0,*)'>>> RK3TVD_FV_himem: hyperboundary()'
            call hyperboundary(u0, dxu, dyu, dzu, urk1, v, w, imask, par)
       end if

       if ( ltrace  ) write(0,*)'>>> RK3TVD_FV_himem: hyperrhs()'
       if ( ltrace2 ) then
          do m=1,NU
             write(0,*) ' u0   - 1 ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
          end do
          do m=1,NV
             write(0,*) ' v    - 1 ',m,myl2norm3d(v(m)%d, nx, ny, nz)
          end do
       end if
       if ( ltrace ) then
          write(0,*) 'RK3TVD_FV_himem: iter = ', iter 
       endif
       call SolvePrimVars(v, w, u0, par)
!      if ( ltrace ) then
!         write(0,*) 'RK3TVD_FV_himem: After 1st SolvePrimVars'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(1,1,1) = ',u0(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(1,1,1) = ',u0(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(1,1,1) = ',u0(H_Sz)%d(1,1,1)
!      endif
       call hyperrhs(   urk1,  u0,   v,    w,              & 
   &                           dxu,  dyu,  dzu,            &
   &                           dxv,  dyv,  dzv, imask2, par)
       
!      if ( ltrace ) then
!         write(0,*) 'RK3TVD_FV_himem: After Calculation of hyperrhs'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(1,1,1)   = ',u0(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(1,1,1)   = ',u0(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(1,1,1)   = ',u0(H_Sz)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: urk1(H_Sx)%d(1,1,1) = ',urk1(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: urk1(H_Sy)%d(1,1,1) = ',urk1(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: urk1(H_Sz)%d(1,1,1) = ',urk1(H_Sz)%d(1,1,1)
!      endif
       if ( ltrace2 ) then
          do m=1,NU
             write(0,*) ' u0   - A ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
             write(0,*) ' urk1 - A ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
          end do
          write(0,*) ''
       end if

       if (dissipation > 0) then
          if ( ltrace ) write(0,*)'>>> RK3TVD_FV_himem: hyperdissipation()', sigma_diss
          do m=1,NU
             if (u0(m)%dissipation .eq. 1) then
                call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,      &
           &                          u0(m)%d,   w,         imask,      par)
                urk1(m)%d = urk1(m)%d + sigma_diss*(dxu_diss +  dyu_diss +  dzu_diss)
             end if
          end do
       end if
       if ( ltrace2 ) then
           write(0,*)'>>> RK3TVD_FV_himem: after hyperdissipation()'
       endif
       if ( boundary_conditions .gt. 2 .and.  &
            bc_type             .ne. 3 .and.  .not.no_bc) then
!          if ( ltraceB ) then  
!               write(0,*)'>>> RK3TVD_FV_himem: hyperboundary'
!               write(0,*)'>>> RK3TVD_FV_himem: dt_gxx = ',urk1(H_G11)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_gxy = ',urk1(H_G12)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_gxz = ',urk1(H_G13)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_gyy = ',urk1(H_G22)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_gyz = ',urk1(H_G23)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_gzz = ',urk1(H_G33)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_Sx  = ',urk1(H_Sx)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_Sy  = ',urk1(H_Sy)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_D   = ',urk1(H_D)%d(1,1,1)
!               write(0,*)'>>> RK3TVD_FV_himem: dt_tau = ',urk1(H_TAU)%d(1,1,1)
!           endif
            call hyperboundary(u0, dxu, dyu, dzu, urk1, v, w, imask, par)
        end if
!      if ( ltrace ) then
!         write(0,*) 'RK3TVD_FV_himem: After hyperboundary'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(1,1,1) = ',u0(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(1,1,1) = ',u0(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(1,1,1) = ',u0(H_Sz)%d(1,1,1)
!       endif
        if ( ltrace ) write(0,*)'>>> RK3TVD_FV_himem: update u2()'
        do k = 1, nz
           do j = 1, ny
              do i = 1, nx
                 if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))   &
                      .and. (      amrbound_prepost.eq.1                             &
                      .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))  &
                    ) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     If setting amr boundaries from parents after the evolution step,
!                     then we can go ahead and write over the boundary points:
!                     Otherwise, as long as the point is not an amr or dd boundary point,
!                     then go ahead and write over:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        do m=1, NU
                           if      (update_scheme .eq. 1) then
                              u2(m)%d(i,j,k) = u0(m)%d(i,j,k) + dt *       urk1(m)%d(i,j,k)
                           else if (update_scheme .eq. 0) then
                              u2(m)%d(i,j,k) = u0(m)%d(i,j,k) + dt * 0.5 * urk1(m)%d(i,j,k)
                           else if (update_scheme .eq. 2) then
                              u2(m)%d(i,j,k) = urk1(m)%d(i,j,k)
                           else if (update_scheme .eq. 3) then
                              u2(m)%d(i,j,k) = u0(m)%d(i,j,k) + dt *       urk1(m)%d(i,j,k)
                           else
                              write(0,*) "Error, update_scheme not known"            
                           end if
                        end do

                      else if (nint(chr(i,j,k)) .eq. 1 .and. amrbound_prepost.eq.0) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     Need to interpolate on the boundaries so that
!                     when complete step taken, u2(i,j,k) has the values
!                     given from the parents:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        do m=1, NU
                           u2(m)%d(i,j,k) = 0.5d0*( u0(m)%d(i,j,k) + u2(m)%d(i,j,k) )
                        end do
                      end if
              end do
           end do
        end do
!      if ( ltrace ) then
!         write(0,*) 'RK3TVD_FV_himem: After update of u2'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(1,1,1) = ',u0(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(1,1,1) = ',u0(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(1,1,1) = ',u0(H_Sz)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sx)%d(1,1,1) = ',u2(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sy)%d(1,1,1) = ',u2(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sz)%d(1,1,1) = ',u2(H_Sz)%d(1,1,1)
!       endif
  
        if ( bc_type.eq.3 .and. .not.no_bc) then
            if (ltraceB) write(*,*) 'rk3: Using copy boundary condition:', iter
            do m=1, NU
               call copy_bc(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
            end do
        end if

        if (use_mask .eq. 1) then
           if (ltrace) write(0,*)'>>> RK3TVD_FV_himem: use_mask'
           call mask_gfuncs(u2, w, par)
        end if
       if ( ltrace ) then
          write(0,*) 'RK3TVD_FV_himem: iter = ', iter 
       endif
!      if ( ltrace ) then
!         write(0,*) 'RK3TVD_FV_himem: iter = ', iter
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(1,1,1) = ',u0(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(1,1,1) = ',u0(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(1,1,1) = ',u0(H_Sz)%d(1,1,1)
!         write(0,*) '- - -'
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sx)%d(1,1,1) = ',u2(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sy)%d(1,1,1) = ',u2(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sz)%d(1,1,1) = ',u2(H_Sz)%d(1,1,1)
!         write(0,*) '-----'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(2,2,2) = ',u0(H_Sx)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(2,2,2) = ',u0(H_Sy)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(2,2,2) = ',u0(H_Sz)%d(2,2,2)
!         write(0,*) '- - -'
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sx)%d(2,2,2) = ',u2(H_Sx)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sy)%d(2,2,2) = ',u2(H_Sy)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sz)%d(2,2,2) = ',u2(H_Sz)%d(2,2,2)
!         write(0,*) '------'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(3,3,3) = ',u0(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(3,3,3) = ',u0(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(3,3,3) = ',u0(H_Sx)%d(3,3,3)
!         write(0,*) '- - -'
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sx)%d(3,3,3) = ',u2(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sy)%d(3,3,3) = ',u2(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sz)%d(3,3,3) = ',u2(H_Sx)%d(3,3,3)
!      endif
        call SolvePrimVars(v, w, u2, par)
        if (ltrace) write(0,*)'>>> RK3TVD_FV_himem: Done step 1'
 
!#######################################################################
!####################### 2nd RK step ###################################
!#######################################################################
    else if (iter .eq. 2) then
          if (ltrace) then
              write(0,*)'>>> RK3TVD_FV_himem: begin second step'
          end if

        if ( boundary_conditions .eq. 1 .or.             &
             boundary_conditions .eq. 2 .and. .not.no_bc &
           ) then

            if (ltraceB) write(0,*)'>>> RK3TVD_FV_himem: hyperboundary()'
            call hyperboundary(u2, dxu, dyu, dzu, urk2, v, w, imask, par)
        end if

        call hyperrhs(  urk2,  u2,   v,    w,              &
   &                           dxu,  dyu,  dzu,            &
   &                           dxv,  dyv,  dzv, imask2, par)

        if (dissipation > 0) then
           if (ltrace) write(0,*)'>>> RK3TVD_FV_himem: hyperdissipation()', sigma_diss
           do m=1,NU
              if (u2(m)%dissipation .eq. 1) then
                 call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,     &
   &                                   u2(m)%d,   w,         imask,     par)  
                 urk2(m)%d = urk2(m)%d + sigma_diss*(dxu_diss +  dyu_diss +  dzu_diss)
              end if
           end do
        end if

        if ( boundary_conditions .gt. 2 .and.  &
   &         bc_type             .ne. 3 .and. .not.no_bc ) then
             if (ltraceB) write(0,*)'>>> RK3TVD_FV_himem: hyperboundary()'
             call hyperboundary(u2, dxu, dyu, dzu, urk2, v, w, imask, par)
        end if

        do k = 1, nz
           do j = 1, ny
              do i = 1, nx
                 if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))   &
                      .and. (      amrbound_prepost.eq.1                             &
                      .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))  &
                    ) then

                      do m=1, NU
                         if (update_scheme.eq.1) then
                            u2(m)%d(i,j,k) =u0(m)%d(i,j,k) + dt * 0.5 * ( urk2(m)%d(i,j,k) + &
                                                                          urk1(m)%d(i,j,k) )
                         else if (update_scheme .eq. 0) then
                            u2(m)%d(i,j,k) =u0(m)%d(i,j,k) + dt * 0.75 * urk2(m)%d(i,j,k)
                         else if (update_scheme .eq. 2) then
                            continue
                         else if (update_scheme .eq. 3) then
                            u2(m)%d(i,j,k) = 0.75* u0(m)%d(i,j,k) + 0.25 * u2(m)%d(i,j,k) &
                                                           + 0.25 * dt * urk2(m)%d(i,j,k)
                         else
                            write(0,*) "Error, update_scheme not known"            
                         end if
                      end do
                 else if (nint(chr(i,j,k)) .eq. 1 .and. amrbound_prepost.eq.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Need to interpolate on the boundaries so that
!                when complete step taken, u2(i,j,k) has the values
!                given from the parents:
!                Since we ve already overwritten the original u2,
!                we need to reconstruct it.
!                We need to stick into u2() the value:
!                       (3/4) u2tilde + (1/4) u0
!                where u2tilde() is the original value of u2.
!                Since we overwrote u2() with:
!                      u2 = (1/2) u2tilde + (1/2) u0
!                so
!                      u2tilde = 2 u2 - u0
!                so that
!                      u2 = (3/2) u2tilde - (1/2) u0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      do m=1, NU
                         u2(m)%d(i,j,k) = 1.5d0*u2(m)%d(i,j,k) - 0.5d0*u0(m)%d(i,j,k)
                      end do
                 end if
              end do
           end do
        end do

        if (bc_type.eq.3.and. .not.no_bc) then
           if (ltraceB) write(*,*) 'rk3: Using copy boundary condition:',iter
           do m=1, NU
              call copy_bc(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
           end do
        end if

        if (use_mask .eq. 1) then
           call mask_gfuncs(u2, w, par)
        end if
       if ( ltrace ) then
          write(0,*) 'RK3TVD_FV_himem: iter = ', iter 
       endif
        call SolvePrimVars(v, w, u2, par) 
        if (ltrace) write(0,*)'>>> RK3TVD_FV_himem: Done step 2' 

!#######################################################################
!####################### 3rd RK step ###################################
!#######################################################################
    else if (iter .eq. 3) then
        if (ltrace) then
            write(0,*)'>>> RK3TVD_FV_himem: begin third step'
        end if

        if (boundary_conditions .eq. 1 .or.    &
            boundary_conditions .eq. 2 .and. .not.no_bc) then
            if (ltraceB) write(0,*)'>>> RK3TVD_FV_himem: hyperboundary()'
            call hyperboundary(u2, dxu, dyu, dzu, urk3, v, w, imask, par)
        end if

        call hyperrhs(  urk3,  u2,   v,    w,              &
   &                           dxu,  dyu,  dzu,            &
   &                           dxv,  dyv,  dzv, imask2, par) 

        if (dissipation > 0) then
           if (ltrace) write(0,*)'>>> RK3TVD_FV_himem: hyperdissipation()', sigma_diss
           do m=1,NU
              if (u2(m)%dissipation .eq. 1) then
                 call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,     &
   &                                   u2(m)%d,   w,         imask,     par)
                 urk3(m)%d = urk3(m)%d + sigma_diss*(dxu_diss +  dyu_diss +  dzu_diss)
              end if
           end do
        end if

        if (boundary_conditions .gt. 2 .and.            &
   &        bc_type             .ne. 3 .and. .not.no_bc &
   &       ) then
             if (ltraceB) write(0,*)'>>> RK3TVD_FV_himem: hyperboundary()'
             call hyperboundary(u2, dxu, dyu, dzu, urk3, v, w, imask, par)
        end if

        do k = 1, nz
           do j = 1, ny
              do i = 1, nx
                 if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))    &
                      .and. (      amrbound_prepost.eq.1                              &
                      .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))   &
                    ) then
                      do m=1, NU
                         if (update_scheme.eq.1) then
                            u2(m)%d(i,j,k) =u0(m)%d(i,j,k) + dt * 0.5 * ( urk3(m)%d(i,j,k) + &
                                                                          urk1(m)%d(i,j,k) )
                         else if (update_scheme .eq. 0) then
                            u2(m)%d(i,j,k) = u0(m)%d(i,j,k)                                    &
                                           + dt/9.0 * (2.*urk1(m)%d(i,j,k)+3.*urk2(m)%d(i,j,k) &
                                           + 4.0 * urk3(m)%d(i,j,k))
                         else if (update_scheme .eq. 2) then
                            continue
                         else if (update_scheme .eq. 3) then
                            u2(m)%d(i,j,k) =      one_third * u0(m)%d(i,j,k)     &
                                             + 2.*one_third * u2(m)%d(i,j,k)     &
                                             + 2.*one_third * dt * urk3(m)%d(i,j,k)
                         else
                            write(0,*) "Error, update_scheme not known"            
                         end if
                      end do
                 else if (nint(chr(i,j,k)) .eq. 1 .and. amrbound_prepost.eq.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Need to interpolate on the boundaries so that
!                when complete step taken, u2(i,j,k) has the values
!                given from the parents:
!                Interpolation not straightforward, see comments above.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      do m=1, NU
                         u2(m)%d(i,j,k) = (4.d0*u2(m)%d(i,j,k) - u0(m)%d(i,j,k) ) /3.d0
                      end do
                 end if
              end do
           end do
        end do
        if (bc_type.eq.3.and. .not.no_bc) then
            if (ltraceB) write(*,*) 'rk3: Using copy boundary condition:',iter
            do m=1, NU
              call copy_bc(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
            end do
        end if

        if (use_mask .eq. 1) then
           call mask_gfuncs(u2, w, par)
        end if

        if (boundary_conditions .eq. 1 .or.            &    
            boundary_conditions .eq. 2 .and.           &
            bc_type             .ne. 3 .and. .not.no_bc) then
            if (ltraceB) write(0,*)'>>> RK3TVD_FV_himem: hyperboundary()'
            call hyperboundary(u2, dxu, dyu, dzu, urk3, v, w, imask, par)
        end if
!      if ( ltrace ) then
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(1,1,1) = ',u0(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(1,1,1) = ',u0(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(1,1,1) = ',u0(H_Sz)%d(1,1,1)
!         write(0,*) '- - -'
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sx)%d(1,1,1) = ',u2(H_Sx)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sy)%d(1,1,1) = ',u2(H_Sy)%d(1,1,1)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sz)%d(1,1,1) = ',u2(H_Sz)%d(1,1,1)
!         write(0,*) '-----'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(2,2,2) = ',u0(H_Sx)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(2,2,2) = ',u0(H_Sy)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(2,2,2) = ',u0(H_Sz)%d(2,2,2)
!         write(0,*) '- - -'
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sx)%d(2,2,2) = ',u2(H_Sx)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sy)%d(2,2,2) = ',u2(H_Sy)%d(2,2,2)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sz)%d(2,2,2) = ',u2(H_Sz)%d(2,2,2)
!         write(0,*) '------'
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sx)%d(3,3,3) = ',u0(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sy)%d(3,3,3) = ',u0(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u0(H_Sz)%d(3,3,3) = ',u0(H_Sx)%d(3,3,3)
!         write(0,*) '- - -'
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sx)%d(3,3,3) = ',u2(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sy)%d(3,3,3) = ',u2(H_Sx)%d(3,3,3)
!         write(0,*) 'RK3TVD_FV_himem: u2(H_Sz)%d(3,3,3) = ',u2(H_Sx)%d(3,3,3)
!      endif
       call SolvePrimVars(v, w, u2, par)
!      call HyperAnalysis(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, &
!                          urk2, urk3, par)


    end if !---iter

!#######################################################################
!####################### End of RK   ###################################
!#######################################################################
    
  
    if (allocated(dxu_diss) .and.  &
        allocated(dyu_diss) .and.  &
        allocated(dzu_diss) ) then
        if (ltrace2) write(0,*)'>>> RK3TVD_FV_himem: about to deallocate diss vars...'
        deallocate(dzu_diss, STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'>>> RK3TVD_FV_himem: Could not deallocate dzu_diss: istat=', istat
        end if
        deallocate(dxu_diss, STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'>>> RK3TVD_FV_himem: Could not deallocate dxu_diss: istat=', istat
        end if
        deallocate(dyu_diss, STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'>>> RK3TVD_FV_himem: Could not deallocate dyu_diss: istat=', istat
        end if
        if (ltrace2) write(0,*)'>>> RK3TVD_FV_himem: ...diss vars deallocated'
    else
        write(0,*)'>>> RK3TVD_FV_himem: diss vars not allocated as expected!'
    end if


    if (use_mask .eq. 1) then
      if (allocated(imask)) then
        deallocate(imask,imask2)
      else
         write(0,*)'>>> RK3TVD_FV_himem: imask not allocated as expected!'
      end if
    end if

    if (ltrace) then
       write(0,*)'>>>RK3TVD_FV_himem: return'
       if (ltrace2) then
          do m=1,NU
             write(0,*) ' u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
             write(0,*) ' u0   - ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
             write(0,*) ' urk1 - ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
             write(0,*) ' urk2 - ',m,myl2norm3d(urk2(m)%d, nx, ny, nz)
             write(0,*) ' urk3 - ',m,myl2norm3d(urk3(m)%d, nx, ny, nz)
             write(0,*) ' '
          end do
       end if
    end if
  return
end subroutine TVDRK3_FV_himem

!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################

