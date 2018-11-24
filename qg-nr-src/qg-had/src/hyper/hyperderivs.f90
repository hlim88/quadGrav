!-----------------------------------------------------------------
!
!  $Id$
!
!-----------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_DERIVS

  use params
  use GF
  use HYPER_DERIVS_21
  use HYPER_DERIVS_21_PT
  use HYPER_DERIVS_42  
  use HYPER_DERIVS_42_PT  
  use HYPER_DERIVS_43
  use HYPER_DERIVS_43_PT
 
contains

!-----------------------------------------------------------------
!
! SUBROUTINE HYPERDERIVS
!
!-----------------------------------------------------------------
  subroutine hyperderivs(dxu, dyu, dzu, dxv, dyv, dzv, &
                         u, v, w, imask, par)
    implicit none
 
    type(gridfunction), dimension(:)   :: dxu,dyu,dzu,u
    type(gridfunction), dimension(:)   :: dxv,dyv,dzv,v
    type(gridfunction), dimension(NW)   :: w
    CCTK_INT, dimension(:,:,:)          :: imask
    CCTK_REAL, dimension(NPAR)             :: par

    ! local vars
    CCTK_INT                            :: direc, m, use_mask, order
    CCTK_INT                            :: lnu, lnv
    CCTK_REAL, pointer, dimension(:,:,:):: mask, x, y, z
    target                              :: w
    CCTK_INT                            :: nx, ny, nz, myid
    logical, parameter                  :: ltrace  = .false.
    logical, parameter                  :: ltrace2 = .false.
    logical                             :: mytrace, double_equal
    CCTK_REAL                           :: myl2norm3d
    CCTK_REAL                           :: minx,miny,minz,maxx,maxy,maxz, time
    CCTK_INT                            :: proc_return_myid, shp(3)



    lnu = size(u)
    lnv = size(v)
     
    use_mask = nint(par(P_USE_MASK))
    order    = nint(par(P_DERIV_ORDER))

    mask => w(H_MASK)%d
    !x    => w(H_X)%d
    !y    => w(H_Y)%d
    !z    => w(H_Z)%d
    x    => w(H_XPHYS)%d
    y    => w(H_YPHYS)%d
    z    => w(H_ZPHYS)%d
    shp = shape(z)


    myid = proc_return_myid()
    nx   = nint(par(P_NX))
    ny   = nint(par(P_NY))
    nz   = nint(par(P_NZ))
    nx   = shp(1)
    ny   = shp(2)
    nz   = shp(3)
    minx = x( 1,1,1)
    maxx = x(nx,1,1)
    miny = y(1, 1,1)
    maxy = y(1,ny,1)
    minz = z(1,1, 1)
    maxz = z(1,1,nz)
    time = par(P_TIME)
    if (ltrace) then
      write(0,*)myid,'] ************** hyperderivs *******************'
      write(0,*)myid,'] hyperderivs: lnu, lnuv: ',lnu,lnv
      write(0,*)myid,'] hyperderivs: nx/y/z:    ',nx,ny,nz
      write(0,*)myid,'] hyperderivs: use_mask:  ',use_mask
      write(0,*)myid,'] hyperderivs: order:     ',order
      write(0,*)myid,'] hyperderivs: size(dxu/dxv):',size(dxu),size(dxv)
      write(0,*)myid,'] hyperderivs: minx/y/z:',minx,miny,minz
      write(0,*)myid,'] hyperderivs: miny:',miny,y(1,1,1)
      write(0,*)myid,'] hyperderivs: minz:',minz,z(1,1,1)
      write(0,*)myid,'] hyperderivs: time:',time
    end if

    mytrace = .false.
    !if (myid.eq.1.and.nx.eq.39.and.ny.eq.41.and.nz.eq.39.and.double_equal(minx,-5.0)) mytrace=.true.
    !if (myid.eq.1.and.nx.eq.39.and.ny.eq.41.and.nz.eq.39.and.double_equal(minx,-5.0).and.double_equal(miny,-1.66666666666667)) mytrace=.true.
    if (mytrace) write(0,*)myid,'] >> HYPERDERIVS mytrace on'


    if (use_mask .eq. 1) then
      if (ltrace) write(0,*)myid, '] >>> HYPERDERIVS: use_mask == 1'
      !---------------------------------------
      ! MASKED DERIVATIVES
      !---------------------------------------
      if (order .eq. 2) then
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          direc = 1
          if (u(m)%take_dx .eq. 1) then
            call D21(dxu(m)%d,u(m)%d,direc,x,y,z,mask,imask,par)
          end if

          direc = 2
          if (u(m)%take_dy .eq. 1) then
            call D21(dyu(m)%d,u(m)%d,direc,x,y,z,mask,imask,par)
          end if

          direc = 3
          if (u(m)%take_dz .eq. 1) then
            call D21(dzu(m)%d,u(m)%d,direc,x,y,z,mask,imask,par)
          end if
        end do
        !
        ! Derivatives of V
        !
        do m = 1, lnv
          direc = 1
          if (v(m)%take_dx .eq. 1) then
            call D21(dxv(m)%d,v(m)%d,direc,x,y,z,mask,imask,par)
          end if

          direc = 2
          if (v(m)%take_dy .eq. 1) then
            call D21(dyv(m)%d,v(m)%d,direc,x,y,z,mask,imask,par)
          end if

          direc = 3
          if (v(m)%take_dz .eq. 1) then
            call D21(dzv(m)%d,v(m)%d,direc,x,y,z,mask,imask,par)
          end if
        end do

      else if (order .eq. 42) then   
        if (ltrace) write(0,*)myid,'] >>> HYPERDERIVS: order == 42 masked, U derivs'
        !
        ! Derivatives of U
        !
        if (mytrace) write(0,*)myid,'] >> HYPERDERIVS: outputting pre'
        if (mytrace) call  output3Di(imask,      "HDimask",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(dyu(H_E3)%d,"HDdy_E3pre",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(  u(H_E3)%d,"HDE3",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(  u(H_B2)%d,"HDB2pre",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(  u(H_B3)%d,"HDB3pre",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(dxu(H_B2)%d,"HDdx_B2pre",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(dxu(H_B3)%d,"HDdx_B3pre",minx,maxx,miny,maxy,minz,maxz,time,myid)
        if (mytrace) then
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4     u(H_B1)=',    u(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dxu(H_B1)=',  dxu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dyu(H_B1)=',  dyu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dzu(H_B1)=',  dzu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4     u(H_E3)=',    u(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dxu(H_E3)=',  dxu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dyu(H_E3)=',  dyu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dzu(H_E3)=',  dzu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4     u(H_B2)=',    u(H_B2)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dxu(H_B2)=',  dxu(H_B2)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4     u(H_B3)=',    u(H_B3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.4   dxu(H_B3)=',  dxu(H_B3)%d(31,32,14)
        end if
        do m = 1, lnu
          direc = 1
          if (u(m)%take_dx .eq. 1) then
            call D42(dxu(m)%d,u(m)%d,direc,x,y,z,mask,imask,par)
          end if

          direc = 2
          if (u(m)%take_dy .eq. 1) then
            call D42(dyu(m)%d,u(m)%d,direc,x,y,z,mask,imask,par)
          end if

          direc = 3
          if (u(m)%take_dz .eq. 1) then
            call D42(dzu(m)%d,u(m)%d,direc,x,y,z,mask,imask,par)
          end if
        end do
        if (mytrace) write(0,*)myid,'] >> HYPERDERIVS: outputting pos'
        !if (mytrace) call  output3DB(dyu(H_E3)%d,"HDdy_E3pos",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(  u(H_B2)%d,"HDB2pos",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(  u(H_B3)%d,"HDB3pos",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(dxu(H_B2)%d,"HDdx_B2pos",minx,maxx,miny,maxy,minz,maxz,time,myid)
        !if (mytrace) call  output3DB(dxu(H_B3)%d,"HDdx_B3pos",minx,maxx,miny,maxy,minz,maxz,time,myid)
        if (mytrace) then
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6     u(H_B1)=',    u(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dxu(H_B1)=',  dxu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dyu(H_B1)=',  dyu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dzu(H_B1)=',  dzu(H_B1)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6     u(H_E3)=',    u(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dxu(H_E3)=',  dxu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dyu(H_E3)=',  dyu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dzu(H_E3)=',  dzu(H_E3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6     u(H_B2)=',    u(H_B2)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dxu(H_B2)=',  dxu(H_B2)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6     u(H_B3)=',    u(H_B3)%d(31,32,14)
         !write(0,*)myid, '] >>> HYPERDERIVS 0.6   dxu(H_B3)=',  dxu(H_B3)%d(31,32,14)
        end if
        !
        ! Derivatives of V
        !
        if (ltrace) write(0,*)myid,'] >>> HYPERDERIVS: order == 42 masked, V derivs'
        do m = 1, lnv
          direc = 1
          if (v(m)%take_dx .eq. 1) then
            call D42(dxv(m)%d,v(m)%d,direc,x,y,z,mask,imask,par)
          end if

          direc = 2
         if (v(m)%take_dy .eq. 1) then
            call D42(dyv(m)%d,v(m)%d,direc,x,y,z,mask,imask,par)
         end if

          direc = 3
         if (v(m)%take_dz .eq. 1) then
            call D42(dzv(m)%d,v(m)%d,direc,x,y,z,mask,imask,par)
          end if
        end do
        !

      else
        write(0,*)myid,'A: unmasked derivatives not implemented for order ',order
        call my_exit('hyperderivs: A: unmasked derivatives not implemented')
        !stop
      end if
    else
      if (ltrace) write(0,*)myid,'] >>> HYPERDERIVS: use_mask != 1'
      !---------------------------------------
      ! UNMASKED DERIVATIVES
      !---------------------------------------
      if (order .eq. 2) then
        if (ltrace) write(0,*)myid,'] >>> HYPERDERIVS: order == 2, derivs of U:',lnu
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          direc = 1
          if (u(m)%take_dx .eq. 1) then
            call D21NoMask(dxu(m)%d,u(m)%d,direc,par)
          end if

          direc = 2
          if (u(m)%take_dy .eq. 1) then
            call D21NoMask(dyu(m)%d,u(m)%d,direc,par)
          end if

          direc = 3
          if (u(m)%take_dz .eq. 1) then
            call D21NoMask(dzu(m)%d,u(m)%d,direc,par)
          end if
        end do
        !
        ! Derivatives of V
        !
        if (ltrace) write(0,*)myid,'] >>> HYPERDERIVS:              derivs of V:',lnv
        do m = 1, lnv
          direc = 1
          if (v(m)%take_dx .eq. 1) then
            !if (ltrace) write(0,*)'>>> HYPERDERIVS: x-deriv for v',m
            call D21NoMask(dxv(m)%d,v(m)%d,direc,par)
          end if

          direc = 2
          if (v(m)%take_dy .eq. 1) then
            !if (ltrace) write(0,*)'>>> HYPERDERIVS: y-deriv for v',m
            call D21NoMask(dyv(m)%d,v(m)%d,direc,par)
          end if

          direc = 3
          if (v(m)%take_dz .eq. 1) then
            !if (ltrace) write(0,*)'>>> HYPERDERIVS: z-deriv for v',m
            call D21NoMask(dzv(m)%d,v(m)%d,direc,par)
          end if
        end do

!4th - 3rd order
     else if (order .eq. 4) then
        if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: order == 4, derivs of U:',lnu
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          direc = 1
          if (u(m)%take_dx .eq. 1) then
            call D43NoMask(dxu(m)%d,u(m)%d,direc,par)
          end if

          direc = 2
          if (u(m)%take_dy .eq. 1) then
            call D43NoMask(dyu(m)%d,u(m)%d,direc,par)
          end if

          direc = 3
          if (u(m)%take_dz .eq. 1) then
            call D43NoMask(dzu(m)%d,u(m)%d,direc,par)
          end if
        end do
        !
        ! Derivatives of V
        !
        if (ltrace) write(0,*)myid,'>>> HYPERDERIVS:             derivs of V:',lnv
        do m = 1, lnv
          direc = 1
          if (v(m)%take_dx .eq. 1) then
            !if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: x-deriv for v',m
            call D43NoMask(dxv(m)%d,v(m)%d,direc,par)
          end if

          direc = 2
          if (v(m)%take_dy .eq. 1) then
            !if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: y-deriv for v',m
            call D43NoMask(dyv(m)%d,v(m)%d,direc,par)
          end if

          direc = 3
          if (v(m)%take_dz .eq. 1) then
            !if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: z-deriv for v',m
            call D43NoMask(dzv(m)%d,v(m)%d,direc,par)
          end if
        end do


!4th - 2rd order
     else if (order .eq. 42 .or.  &
  &           order .eq. 44 .or.  &
  &           order .eq. 642  .or.  &
  &           order .eq. 666) then
        ! The BSSN code uses orders 44, 642, 666, and others. In the
        ! BSSN code, however, these derivative routines are not called
        ! through hyper, but rather in the BSSN code itself.  Thus, the
        ! the distinction between these different options is only active
        ! in the BSSN code.  The derivative
        ! operators are set through the deriv_order parameter, and thus
        ! we need to include options for these other derivative operators
        ! here.  As the derivatives taken through hyper are not used in
        ! the calculation of the BSSN evolution equations, we simply call
        ! the 42 operators for now.  Maybe in the future, someone will
        ! be ambitious enough to add the derivative operators here.
        if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: order == 42, derivs of U:',lnu
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          direc = 1
          if (u(m)%take_dx .eq. 1) then
            call D42NoMask(dxu(m)%d,u(m)%d,direc,par)
          end if

          direc = 2
          if (u(m)%take_dy .eq. 1) then
            call D42NoMask(dyu(m)%d,u(m)%d,direc,par)
          end if

          direc = 3
          if (u(m)%take_dz .eq. 1) then
            call D42NoMask(dzu(m)%d,u(m)%d,direc,par)
          end if
        end do
        !
        ! Derivatives of V
        !
        if (ltrace) write(0,*)myid,'>>> HYPERDERIVS:             derivs of V:',lnv
        do m = 1, lnv
          direc = 1
          if (v(m)%take_dx .eq. 1) then
            !if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: x-deriv for v',m
            call D42NoMask(dxv(m)%d,v(m)%d,direc,par)
          end if

          direc = 2
          if (v(m)%take_dy .eq. 1) then
            !if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: y-deriv for v',m
            call D42NoMask(dyv(m)%d,v(m)%d,direc,par)
          end if

          direc = 3
          if (v(m)%take_dz .eq. 1) then
            !if (ltrace) write(0,*)myid,'>>> HYPERDERIVS: z-deriv for v',m
            call D42NoMask(dzv(m)%d,v(m)%d,direc,par)
          end if
        end do

     else if (order .eq. 44 .or. order .eq. 642 .or. order .eq. 666) then

     else 
        write(0,*)'B: unmasked derivatives not implemented for order ',order
        stop
       end if
     end if

  if (ltrace2) then
     write(0,*) myid,'] >>>Hyperderivs: Norms aftewards:'
     do m=1,lNU
        write(0,*) myid,'] Der:  ||u||   ',m,myl2norm3d(u(m)%d, nx, ny, nz)
        if (u(m)%take_dx .eq. 1) then
           write(0,*) myid,'] Der: ||u_x||  ',m,myl2norm3d(dxu(m)%d, nx, ny, nz)
        end if
        if (u(m)%take_dy .eq. 1) then
           write(0,*) myid,'] Der: ||u_y||  ',m,myl2norm3d(dyu(m)%d, nx, ny, nz)
        end if
        if (u(m)%take_dz .eq. 1) then
           write(0,*) myid,'] Der: ||u_z||  ',m,myl2norm3d(dzu(m)%d, nx, ny, nz)
        end if
     end do
     do m=1,lNV
        write(0,*) myid,'] Der:  ||v||   ',m,myl2norm3d(v(m)%d, nx, ny, nz)
        if (v(m)%take_dx .eq. 1) then
           write(0,*) myid,'] Der: ||v_x||  ',m,myl2norm3d(dxv(m)%d, nx, ny, nz)
        end if
        if (v(m)%take_dy .eq. 1) then
           write(0,*) myid,'] Der: ||v_y||  ',m,myl2norm3d(dyv(m)%d, nx, ny, nz)
        end if
        if (v(m)%take_dz .eq. 1) then
           write(0,*) myid,'] Der: ||v_z||  ',m,myl2norm3d(dzv(m)%d, nx, ny, nz)
        end if
     end do
     write(0,*)myid,'] >>> HYPERDERIVS: Done.'
  end if

    return
  end subroutine hyperderivs


!-----------------------------------------------------------------
!
! SUBROUTINE HYPERDERIVS_PT
!
! Point-wise derivatives.
!
!-----------------------------------------------------------------
  subroutine hyperderivs_pt(dxu, dyu, dzu, dxv, dyv, dzv, &
                            u, v, wpt, imask, i, j, k, par)
    implicit none
 
    CCTK_REAL, dimension(:)   :: dxu, dyu, dzu
    CCTK_REAL, dimension(:)   :: dxv, dyv, dzv
    type(gridfunction), dimension(:)   :: u
    type(gridfunction), dimension(:)   :: v
    CCTK_REAL, dimension(NW)            :: wpt
    CCTK_INT                            :: imask, point_wise_analysis
    CCTK_REAL, dimension(:)             :: par
    CCTK_INT                            :: i, j, k

    CCTK_INT                            :: nx, ny, nz
    CCTK_INT                            :: lnu, lnv
    CCTK_INT                            :: m, use_mask, order
    CCTK_REAL, parameter                :: ERROR_VALUE = 1.0d200
    logical, parameter                  :: ltrace  = .false.
    logical, parameter                  :: ltrace2 = .false.

    if (ltrace) then
      write(0,*)'Entering hyperderivs_pt'
    end if

    lnu = size(u)
    lnv = size(v)

    use_mask = nint(par(P_USE_MASK))
    order    = nint(par(P_DERIV_ORDER))
    point_wise_analysis = nint(par(P_POINT_WISE_ANALYSIS))

    nx = nint(par(P_NX))
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))

    if (i .lt. 1 .or. i .gt. nx) then
      write(0,*)'HYPERDERIVS_PT >> i out of range: i/nx = ',i,nx
    end if
    if (j .lt. 1 .or. j .gt. ny) then
      write(0,*)'HYPERDERIVS_PT >> j out of range: y/ny = ',j,ny
    end if
    if (k .lt. 1 .or. k .gt. nz) then
      write(0,*)'HYPERDERIVS_PT >> k out of range: k/nz = ',k,nz
    end if

    if (use_mask .eq. 1) then
      !---------------------------------------
      ! MASKED DERIVATIVES FOR 2-1 
      !---------------------------------------
      if (order .eq. 2) then
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          if (u(m)%take_dx .eq. 1 .or. u(m)%take_dy .eq. 1  &
                                  .or. u(m)%take_dz .eq. 1  &
	                          .or. point_wise_analysis .eq. 1) then
            call D21_pt(dxu(m), dyu(m), dzu(m), u(m)%d,   &
                        wpt,    imask,  i,      j,      k,      par)
          else
            dxu(m) = ERROR_VALUE
            dyu(m) = ERROR_VALUE
            dzu(m) = ERROR_VALUE
          end if
        end do
        !
        ! Derivatives of V
        !
        do m = 1, lnv
          if (v(m)%take_dx .eq. 1 .or. v(m)%take_dy .eq. 1  &
                                  .or. v(m)%take_dz .eq. 1) then
            call D21_pt(dxv(m),  dyv(m),  dzv(m),  v(m)%d,  &
                        wpt,     imask,   i,       j,       k,    par)
          else
            dxv(m) = ERROR_VALUE
            dyv(m) = ERROR_VALUE
            dzv(m) = ERROR_VALUE
          end if
        end do


      !---------------------------------------
      ! MASKED DERIVATIVES FOR 4-2
      !---------------------------------------
      else if (order .eq. 42) then
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          if (u(m)%take_dx .eq. 1 .or. u(m)%take_dy .eq. 1  &
                                  .or. u(m)%take_dz .eq. 1  &
	                          .or. point_wise_analysis .eq. 1) then
            call D42_pt(dxu(m), dyu(m), dzu(m), u(m)%d,   &
                        wpt,    imask,  i,      j,      k,      par)
          else
            dxu(m) = ERROR_VALUE
            dyu(m) = ERROR_VALUE
            dzu(m) = ERROR_VALUE
          end if
        end do
        !
        ! Derivatives of V
        !
        do m = 1, lnv
          if (v(m)%take_dx .eq. 1 .or. v(m)%take_dy .eq. 1  &
                                  .or. v(m)%take_dz .eq. 1) then
            call D42_pt(dxv(m),  dyv(m),  dzv(m),  v(m)%d,  &
                        wpt,     imask,   i,       j,       k,    par)
          else
            dxv(m) = ERROR_VALUE
            dyv(m) = ERROR_VALUE
            dzv(m) = ERROR_VALUE
          end if
        end do
	
	
      else
        write(0,*)'C: masked point-wise derivatives not implemented'
        write(0,*)'   for order ',order
        stop
      end if
    else
      !---------------------------------------
      ! UNMASKED DERIVATIVES FOR 2-1
      !---------------------------------------
      if (order .eq. 2) then
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          if (u(m)%take_dx .eq. 1 .or. u(m)%take_dy .eq. 1  &
                                  .or. u(m)%take_dz .eq. 1  &
                                  .or. point_wise_analysis .eq. 1) then
            call D21NoMask_pt(dxu(m),dyu(m),dzu(m),u(m)%d,wpt,i,j,k,par)
          else
            dxu(m) = ERROR_VALUE
            dyu(m) = ERROR_VALUE
            dzu(m) = ERROR_VALUE
          end if
        end do
        !
        ! Derivatives of V
        !
        do m = 1, lnv
          if (v(m)%take_dx .eq. 1 .or. v(m)%take_dy .eq. 1   &
                                  .or. v(m)%take_dz .eq. 1) then
            call D21NoMask_pt(dxv(m),dyv(m),dzv(m),v(m)%d,wpt,i,j,k,par)
          else
            dxv(m) = ERROR_VALUE
            dyv(m) = ERROR_VALUE
            dzv(m) = ERROR_VALUE
          end if
        end do

      !---------------------------------------
      ! UNMASKED DERIVATIVES FOR 4-3
      !---------------------------------------	
      else  if (order .eq. 4) then
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          if (u(m)%take_dx .eq. 1 .or. u(m)%take_dy .eq. 1  &
                                  .or. u(m)%take_dz .eq. 1  &
                                  .or. point_wise_analysis .eq. 1) then
            call D43NoMask_pt(dxu(m),dyu(m),dzu(m),u(m)%d,i,j,k,par)
          else
            dxu(m) = ERROR_VALUE
            dyu(m) = ERROR_VALUE
            dzu(m) = ERROR_VALUE
          end if
        end do
        !
        ! Derivatives of V
        !
        do m = 1, lnv
          if (v(m)%take_dx .eq. 1 .or. v(m)%take_dy .eq. 1   &
                                  .or. v(m)%take_dz .eq. 1) then
            call D43NoMask_pt(dxv(m),dyv(m),dzv(m),v(m)%d,i,j,k,par)
          else
            dxv(m) = ERROR_VALUE
            dyv(m) = ERROR_VALUE
            dzv(m) = ERROR_VALUE
          end if
        end do

      !---------------------------------------
      ! UNMASKED DERIVATIVES FOR 4-2
      !---------------------------------------	
      else  if (order .eq. 42) then
        !
        ! Derivatives of U
        !
        do m = 1, lnu
          if (u(m)%take_dx .eq. 1 .or. u(m)%take_dy .eq. 1  &
                                  .or. u(m)%take_dz .eq. 1  &
                                  .or. point_wise_analysis .eq. 1) then
            call D42NoMask_pt(dxu(m),dyu(m),dzu(m),u(m)%d,wpt,i,j,k,par)
          else
            dxu(m) = ERROR_VALUE
            dyu(m) = ERROR_VALUE
            dzu(m) = ERROR_VALUE
          end if
        end do
        !
        ! Derivatives of V
        !
        do m = 1, lnv
          if (v(m)%take_dx .eq. 1 .or. v(m)%take_dy .eq. 1   &
                                  .or. v(m)%take_dz .eq. 1) then
            call D42NoMask_pt(dxv(m),dyv(m),dzv(m),v(m)%d,wpt,i,j,k,par)
          else
            dxv(m) = ERROR_VALUE
            dyv(m) = ERROR_VALUE
            dzv(m) = ERROR_VALUE
          end if
        end do
      
      
      else 
        write(0,*)'D: unmasked point-wise derivatives not implemented '
        write(0,*)'   for order ',order
        stop
      end if
    end if

    return
  end subroutine hyperderivs_pt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output3DB(var,name,minx,maxx,miny,maxy,minz,maxz,time,myproc)
   implicit none
   CCTK_REAL, dimension(:,:,:) :: var
   CCTK_REAL                   :: minx,miny,minz,maxx,maxy,maxz
   CCTK_REAL :: time
   character(len=*) :: name
   
   CCTK_INT myproc, istat
      
   real*8, allocatable, dimension(:) :: tempcoord
   real*8                            :: bbox(6)
   integer :: nx, ny, nz, ret, gft_out_bbox, shp(3)
   character(20):: g11f
   character(32):: cnames
   
   nx=size(var(:,1,1))
   ny=size(var(1,:,1))
   nz=size(var(1,1,:))
   shp = shape(var)
   nx=shp(1)
   ny=shp(2)
   nz=shp(3)
   
   cnames = 'x|y|z'

   bbox(1) = minx
   bbox(2) = maxx
   bbox(3) = miny
   bbox(4) = maxy
   bbox(5) = minz
   bbox(6) = maxz
   !bbox(1) = 0.0
   !bbox(2) = 1.0
   !bbox(3) = 0.0
   !bbox(4) = 1.0
   !bbox(5) = 0.0
   !bbox(6) = 1.0
   
   write (g11f,'(A,1i3)') name,myproc
   if(myproc.lt.0) g11f = name
   
   
   write(0,*)'output3DB >>> shp: ',shp
   write(0,*)'output3DB >>> bbox: ',bbox
   write(0,*)'output3DB >>> nx/y/z: ',nx,ny,nz
   write(0,*)'output3DB >>> minx,maxx: ',minx,maxx
   write(0,*)'output3DB >>> miny,maxy: ',miny,maxy
   write(0,*)'output3DB >>> minz,maxz: ',minz,maxz
   write(0,*)'output3DB >>> time: ',time
   write(0,*)'output3DB >>> cnames: ',cnames
   write(0,*)'output3DB >>> outputting ',g11f,nx,ny,nz
   ret=gft_out_bbox(g11f,time,(/nx,ny,nz/),3,bbox,var)
   write(0,*)'output3DB >>> ret=',ret
            
            
   end subroutine output3DB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output3Di(var,name,minx,maxx,miny,maxy,minz,maxz,time,myproc)
   implicit none
   CCTK_INT, dimension(:,:,:) :: var
   CCTK_REAL                   :: minx,miny,minz,maxx,maxy,maxz
   CCTK_REAL :: time
   character(len=*) :: name
   
   CCTK_INT myproc, istat
      
   real*8, allocatable, dimension(:) :: tempdata
   real*8                            :: bbox(6)
   integer :: nx, ny, nz, ret, gft_out_bbox, shp(3), i,j,k,l
   character(20):: g11f
   character(32):: cnames
   
   nx=size(var(:,1,1))
   ny=size(var(1,:,1))
   nz=size(var(1,1,:))
   shp = shape(var)
   nx=shp(1)
   ny=shp(2)
   nz=shp(3)
   
   cnames = 'x|y|z'

   bbox(1) = minx
   bbox(2) = maxx
   bbox(3) = miny
   bbox(4) = maxy
   bbox(5) = minz
   bbox(6) = maxz
   !bbox(1) = 0.0
   !bbox(2) = 1.0
   !bbox(3) = 0.0
   !bbox(4) = 1.0
   !bbox(5) = 0.0
   !bbox(6) = 1.0
   
   allocate(tempdata(nx*ny*nz),STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'output3Di >>> can not allocate tempcoord'
          write(0,*)'             size: ',nx*ny*nz
        end if

   l = 0
   do i = 1, nx
   do j = 1, ny
   do k = 1, nz
     l = l + 1
     tempdata(l) = 1.0*var(i,j,k)
   end do
   end do
   end do
   

   write (g11f,'(A,1i3)') name,myproc
   if(myproc.lt.0) g11f = name
   
   
   write(0,*)'output3Di >>> shp: ',shp
   write(0,*)'output3Di >>> bbox: ',bbox
   write(0,*)'output3Di >>> nx/y/z: ',nx,ny,nz
   write(0,*)'output3Di >>> minx,maxx: ',minx,maxx
   write(0,*)'output3Di >>> miny,maxy: ',miny,maxy
   write(0,*)'output3Di >>> minz,maxz: ',minz,maxz
   write(0,*)'output3Di >>> time: ',time
   write(0,*)'output3Di >>> cnames: ',cnames
   write(0,*)'output3Di >>> outputting ',g11f,nx,ny,nz
   ret=gft_out_bbox(g11f,time,(/nx,ny,nz/),3,bbox,tempdata)
   write(0,*)'output3Di >>> ret=',ret
            
   if (allocated(tempdata)) then
     deallocate(tempdata)
   else
      write(0,*)'output3Di >>> tempdata not allocated as expected!'
   end if

            
   end subroutine output3Di
   

end MODULE HYPER_DERIVS

