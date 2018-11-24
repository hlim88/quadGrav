!------------------------------------------------------------------------
!
!   $Id$
!
!------------------------------------------------------------------------
#include "cctk.h"

#define IMASK_MODX(A,B) A = mod(A,10000) + B
#define IMASK_MODY(A,B) A = (A/10000)*10000 + mod(A,100) + B
#define IMASK_MODZ(A,B) A = (A/100)*100 + B

#define IMASK_MODXY(A,B,C) A = mod(A,100) + B + C
#define IMASK_MODXZ(A,B,C) A = mod(A/100,100)*100 + B + C
#define IMASK_MODYZ(A,B,C) A = (A/10000)*10000 + B + C

      MODULE MOD_MASK

      use params
      USE GF
      implicit none

      CCTK_REAL, parameter  :: P_INSIDE_MASK = -1.0
      CCTK_INT,  parameter  :: P_INT_INSIDE_MASK = -1
      CCTK_INT,  parameter  :: MAXDOMAINS  = 10
      logical,   parameter  :: fixup_masks = .true.

      CONTAINS
!------------------------------------------------------------------------
!
!     SUBROUTINE SET_MASK
!
!------------------------------------------------------------------------
      subroutine set_mask(w, imask, imask2, par)
      implicit   none

      type(gridfunction), dimension(NW)          :: w
      CCTK_INT,  dimension(:,:,:)                :: imask,imask2
      CCTK_REAL, dimension(NPAR)                    :: par

      ! local vars
      CCTK_REAL, dimension(:,:,:), pointer       :: x, y, z
      CCTK_INT                                   :: order, nholes, mask_type
      CCTK_INT                                   :: dissipation
      CCTK_INT                                   :: shp(3), nx, ny, nz
      CCTK_REAL, pointer, dimension(:,:,:)       :: mask
!     CCTK_REAL               :: posx, posy, posz, rad_exc
      logical                                    :: ltrace = .false.
      integer                                    :: myid, proc_return_myid

      target                                     :: w


      myid = proc_return_myid()
      if (ltrace) write(0,*) myid,'set_mask: Enter'
      mask => w(H_MASK)%d
      x    => w(H_X)%d
      y    => w(H_Y)%d
      z    => w(H_Z)%d

      shp = shape(w(H_MASK)%d)

      nx = shp(1)
      ny = shp(2)
      nz = shp(3)

      nholes = nint(par(P_NBHOLES))
      mask_type = nint(par(P_MASK_TYPE))


      ! fixup_mask goes through the mask to find faces, edges and verticies
      ! of the mask.  After calling fixup_mask, the mask values will be:
      !    1 - exterior of excision region
      !   -1 - interior of excision region
      !   -2 - face of excision region
      !   -3 - edge of excision region
      !   -4 - vertex of excision region
      !
      ! NOTE THERE WILL BE A PROBLEM IF THE MASK IS ONLY A SINGLE POINT
      ! THICK ALONG ANY LINE.  THIS COULD OCCUR WITH A SPHERICAL MASK,
      ! FOR EXAMPLE.  A single point will be counted as a boundary twice. 
      if (fixup_masks) then
         !
         ! These are in ../amr, the first re-normalizes to +/-1
         !
         call unfixup_mask(mask, nx,ny,nz)
         call myfixup_mask(mask, nx,ny,nz)
      end if

      ! Now create the imask (index mask) for the derivative routines.
      ! Different routines are needed depending on the stencil size,
      ! and dissipation stencil size.
      !
      order = nint(par(P_DERIV_ORDER))
      dissipation = nint(par(P_DISSIPATION))
      if (ltrace) write(0,*) myid,'set_mask: Setting imask'
      if (order .eq. 2) then
           call SetImask2_KO(w, imask, par)
           imask2 = imask 
!           call SetImask2_KO(w, imask2, par)
      else if(order .eq. 42) then
           call SetImask2_KO(w, imask, par)
           call SetImask42(w, imask2, par)
      else if(order.eq.43) then
           call SetImask2_KO(w, imask, par)
           call SetImask2(w, imask2, par)
      else 
         write(0,*)myid,'SET_MASK: Masked derivatives not implemented for order ',& 
                    order
      end if
      !if (ltrace) write(*,*) 'calling SetImask2...done'

      ! This routine changes the derivative stencils stored in imask
      ! for proper differencing around edges and verticies.
      if (fixup_masks) then
         if (ltrace) write(0,*) myid,'set_mask: calling fixup_imask'
         !write(*,*) 'calling fixup_imask'
         call fixup_imask(imask, mask, par)
         !write(*,*) 'calling fixup_imask Done'
      end if

      if (ltrace) write(0,*) myid,'set_mask: Done.'

      return
      end subroutine set_mask

!------------------------------------------------------------------------
!
!     SET_MASK_SPHERE
!
!------------------------------------------------------------------------
!     subroutine set_mask_sphere(mask, rad_exc, posx, posy, posz,  &
!                                x,    y,       z,    par)

!     implicit   none

!     CCTK_REAL, dimension(:,:,:)      :: mask
!     CCTK_REAL, dimension(:,:,:)      :: x, y, z
!     CCTK_REAL, dimension(:)          :: par
!     CCTK_REAL                        :: rad_exc, posx, posy, posz

      ! local vars
!     CCTK_INT, dimension(3)     :: shp
!     CCTK_INT                   :: i, j, k, nx, ny, nz
!     CCTK_REAL                  :: rex, dx, fuzz, radsq
!     CCTK_REAL, parameter       :: fuzz_scale = 0.5

      ! For now have this not do anything.
      ! Instead, let code in had/src/amr take care of this
      !    because I am trying to have the mask move, and this
      !    routine keeps writing over what I am trying to do.
      !    We can work out later where this should go once and
      !    for all (SLL 06/06/06).
!     return

!     shp = shape(mask)
!     nx = shp(1)
!     ny = shp(2)
!     nz = shp(3)

      ! Make sure that rx is positive.  I added some abs()
!     dx = par(P_DX)
!     fuzz = fuzz_scale * dx
!     rex = (rad_exc + fuzz)**2

!     do k = 1, nz
!     do j = 1, ny
!     do i = 1, nx
!        radsq = (x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
!        if (radsq < rex) then
!           mask(i,j,k) = P_INSIDE_MASK
!        end if
!     end do
!     end do
!     end do

!     return
!     end subroutine set_mask_sphere

!------------------------------------------------------------------------
!
!     SET_MASK_CUBE
!
!
!------------------------------------------------------------------------
      subroutine set_mask_cube(mask, rad_exc, posx, posy, posz,  &
                                 x,    y,       z,    par)
      implicit   none
      CCTK_REAL, dimension(:,:,:)    :: mask
      CCTK_REAL, dimension(:,:,:)    :: x, y, z
      CCTK_REAL, dimension(NPAR)        :: par
      CCTK_REAL                      :: rad_exc, posx, posy, posz
      
      ! local vars
      CCTK_INT, dimension(3)         :: shp
      CCTK_INT                       :: i, j, k
      CCTK_REAL                      :: exsq, dx, fuzz, xsq, ysq, zsq
      CCTK_REAL, parameter           :: fuzz_scale = 0.1d-4

      ! External function
!     logical                        :: double_equal

      ! For now have this not do anything.
      ! Instead, let code in had/src/amr take care of this
      !    because I am trying to have the mask move, and this
      !    routine keeps writing over what I am trying to do.
      !    We can work out later where this should go once and
      !    for all (SLL 12/2/05).
       return

         shp = shape(mask)
   
         dx = par(P_DX)
         fuzz = fuzz_scale * dx
         exsq = (rad_exc + fuzz)**2
   
         do k = 1, shp(3)
         do j = 1, shp(2)
         do i = 1, shp(1)
            xsq = (x(i,j,k) - posx)**2
            ysq = (y(i,j,k) - posy)**2
            zsq = (z(i,j,k) - posz)**2
            if (xsq < exsq .and. ysq < exsq .and. zsq < exsq) then
               mask(i,j,k) = P_INSIDE_MASK
            end if
         end do
         end do
         end do


      return
      end subroutine set_mask_cube

!------------------------------------------------------------------------
!
!     MASK_GFUNCS
!
!------------------------------------------------------------------------
      subroutine mask_gfuncs(u, w, par)
      use GF
      implicit   none

      type(gridfunction), dimension(NU)   :: u
      type(gridfunction), dimension(NW)   :: w
      CCTK_REAL,          dimension(NPAR) :: par

      ! local vars
      CCTK_INT, dimension(3)              :: shp
      CCTK_INT                            :: i, j, k, m
      CCTK_REAL, dimension(:,:,:), pointer:: mask
      
      mask => w(H_MASK)%d
      shp = shape(mask)

      !write(*,*) 'mask_gfuncs: Enter'
      do k = 1, shp(3)
      do j = 1, shp(2)
      do i = 1, shp(1)
         if (mask(i,j,k) .lt. 0.0) then
            do m = 1, NU
               u(m)%d(i,j,k) = u(m)%excise_val
            end do
         end if
      end do
      end do
      end do
      !write(*,*) 'mask_gfuncs: Done'


      return
      end subroutine mask_gfuncs

!------------------------------------------------------------------------
!
!     OUTPUT_MASK
!
!     Create an output mask from a standard mask.  The output mask
!     has only 1's and 0's.
!
!------------------------------------------------------------------------
!     subroutine output_mask(w, par)
!     use GF
!     implicit   none

!     type(gridfunction), dimension(NW)   :: w
!     CCTK_REAL,          dimension(NPAR) :: par

!     ! local vars
!     CCTK_INT, dimension(3)              :: shp
!     CCTK_INT                            :: i, j, k
!     CCTK_REAL, dimension(:,:,:), pointer:: mask


!     mask => w(H_MASK)%d
!     shp = shape(mask)

!     do k = 1, shp(3)
!     do j = 1, shp(2)
!     do i = 1, shp(1)
!        if (mask(i,j,k) .gt. 0.0) then
!           mask(i,j,k) = 1.0
!        else
!           mask(i,j,k) = 0.0
!        end if
!     end do
!     end do
!     end do


!     return
!     end subroutine output_mask


!------------------------------------------------------------------------
!
!     SetImask2
!
!------------------------------------------------------------------------
      subroutine SetImask2(w, imask, par)
      use GF
      implicit   none

      type(gridfunction), dimension(NW)    :: w
      CCTK_INT, dimension(:,:,:)           :: imask
      CCTK_REAL, dimension(NPAR)           :: par

      ! local vars
      CCTK_INT                             :: i, j, k
      CCTK_INT                             :: im, iml, imr, mask_stencil
      CCTK_INT                             :: nx, ny, nz
      CCTK_INT, dimension(3)               :: shp
      CCTK_REAL, pointer, dimension(:,:,:) :: mask
      logical, parameter                   :: ltrace2 = .false.
      logical, parameter                   :: ltrace3 = .false.
      logical, parameter                   :: paranoid = .false.

      if (ltrace3) then
         write(0,*)'Entering SetImask2'
      end if

      mask  => w(H_MASK)%d


      shp = shape(mask)
      nx = shp(1)
      ny = shp(2)
      nz = shp(3)

      if (ltrace3) then
         write(0,*)'SETIMASK2: minval/maxval(mask) ', &
                   minval(w(H_MASK)%d), maxval(w(H_MASK)%d)
      end if

      !print *,'SetImask2: nx, ny, nz ',nx, ny,nz
      !
      ! Initialize the index mask
      !
      imask = 0
      !print *,'mask: initial imask(1,1,1) ',imask(1,1,1)

      !--------------------------------------------------------
      ! March along x and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do j = 1, ny

         imask(1,j,k) = imask(1,j,k) + P_STENCIL_X_LEFT
         do i = 2, nx - 1

            im = nint(mask(i,j,k))
            iml = nint(mask(i-1,j,k))
            imr = nint(mask(i+1,j,k))
            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !

               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !

                  if (imr .eq. 1) then
                     !
                     ! Check that the point to the right is valid and on 
                     ! the domain
                     !

                     imask(i,j,k) = imask(i,j,k) + P_STENCIL_X_LEFT
                  else
                     imask(i,j,k) = imask(i,j,k) + P_STENCIL_X_UNDEF
                     if (ltrace2) then 
                        write(0,*)'Could not construct a valid stencil for i+1'
                     end if
                  end if

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !

                  if (i > 2 .and. nint(mask(i-2,j,k)) .eq. 1) then
                     !
                     ! Check that the point to the left (i-2) is valid
                     ! and on the domain
                     !

                     imask(i-1,j,k) = imask(i-1,j,k) + P_STENCIL_X_RIGHT
                  else
                     imask(i-1,j,k) = imask(i-1,j,k) + P_STENCIL_X_UNDEF
                     if (ltrace2) then 
                        write(0,*)'Could not construct a valid stencil for i-2'
                     end if
                  end if

               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in x sweep',im, iml,mask(i,j,k)
               end if
            end if
         end do
         imask(nx,j,k) = imask(nx,j,k) + P_STENCIL_X_RIGHT

      end do
      end do

      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'SETIMASK2: Check_imask found illegal value after x sweep'
            stop
         end if
      end if


      !--------------------------------------------------------
      ! March along y and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do i = 1, nx

         imask(i,1,k) = imask(i,1,k) + P_STENCIL_Y_LEFT

         do j = 2, ny - 1
            im = nint(mask(i,j,k))
            iml = nint(mask(i,j-1,k))
            imr = nint(mask(i,j+1,k))
            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !

               if (im .eq. 1 .and. iml .lt. 0) then 
                  !
                  ! The mask is to the left
                  !

                  if (imr .eq. 1) then
                     !
                     ! Check that the point to the right is valid and on
                     ! the domain
                     !

                     imask(i,j,k) = imask(i,j,k) + P_STENCIL_Y_LEFT
                  else
                     imask(i,j,k) = imask(i,j,k) + P_STENCIL_Y_UNDEF
                     if (ltrace2) then
                        write(0,*)'Could not construct a valid stencil for j+1'
                     end if
                  end if

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !

                  if (j > 2 .and. nint(mask(i,j-2,k)) .eq. 1) then
                     !
                     ! Check that the point to the left (j-2) is valid
                     ! and on the domain
                     !

                     imask(i,j-1,k) = imask(i,j-1,k) + P_STENCIL_Y_RIGHT
                  else
                     imask(i,j-1,k) = imask(i,j-1,k) + P_STENCIL_Y_UNDEF
                     if (ltrace2) then
                        write(0,*)'Could not construct a valid stencil for j-2'
                     end if
                  end if

               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in y sweep',im,iml
               end if
            end if
         end do
         imask(i,ny,k) = imask(i,ny,k) + P_STENCIL_Y_RIGHT

      end do
      end do

      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'SETIMASK2: Check_imask found illegal value after y sweep'
            stop
         end if
      end if

      
      !--------------------------------------------------------
      ! March along z and check for mask boundaries
      !--------------------------------------------------------
      do j = 1, ny
      do i = 1, nx

         imask(i,j,1) = imask(i,j,1) + P_STENCIL_Z_LEFT

         do k = 2, nz - 1
            im = nint(mask(i,j,k))
            iml = nint(mask(i,j,k-1))
            imr = nint(mask(i,j,k+1))
            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !

               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !

                  if (imr .eq. 1) then
                     !
                     ! Check that the point to the right is valid and on
                     ! the domain
                     !

                     imask(i,j,k) = imask(i,j,k) + P_STENCIL_Z_LEFT
                  else
                     imask(i,j,k) = imask(i,j,k) + P_STENCIL_Z_UNDEF
                     if (ltrace2) then
                        write(0,*)'Could not construct a valid stencil for k+1'
                     end if
                  end if

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !

                  if (k > 2 .and. nint(mask(i,j,k-2)) .eq. 1) then
                     !
                     ! Check that the point to the left (k-2) is valid
                     ! and on the domain
                     !

                     imask(i,j,k-1) = imask(i,j,k-1) + P_STENCIL_Z_RIGHT
                  else
                     imask(i,j,k-1) = imask(i,j,k-1) + P_STENCIL_Z_UNDEF
                     if (ltrace2) then
                        write(0,*)'Could not construct a valid stencil for k-2'
                     end if
                  end if

               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in z sweep',im,iml
               end if
            end if
         end do
         imask(i,j,nz) = imask(i,j,nz) + P_STENCIL_Z_RIGHT

      end do
      end do

      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'SETIMASK2: Check_imask found illegal value after z sweep'
            stop
         end if
      end if


      !print *,'mask: imask ',imask(1,4,23)
      !print *,'mask: imask ',imask(8,60,23) 

      !
      ! Set imask values inside the mask
      !
      mask_stencil = P_STENCIL_X_MASK + P_STENCIL_Y_MASK + P_STENCIL_Z_MASK
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (mask(i,j,k) .lt. 0.1) then
            imask(i,j,k) = mask_stencil
         end if
      end do
      end do
      end do

      !print *, 'return from SetImask2'
      !print *,'mask: imask ',imask(1,1,1)
      !print *,'mask: imask ',imask(8,60,23) 
      if (Check_imask(imask, par) == 0) then
         write(0,*)'SETIMASK2: Check_imask found illegal value.'
         stop
      end if



      return
      end subroutine SetImask2

!------------------------------------------------------------------------
!
!     SetImask2_KO
!
!     This routine generates the imask for a five point stencil
!     appropriate for Kreiss-Oliger dissipation.
!
!------------------------------------------------------------------------
      subroutine SetImask2_KO(w, imask, par)
      use GF
      implicit   none

      type(gridfunction), dimension(NW)    :: w
      CCTK_INT, dimension(:,:,:)           :: imask
      CCTK_REAL, dimension(NPAR)           :: par

      ! local vars
      CCTK_INT                             :: i, j, k, m, ii, q, istat
      !CTK_INT                             :: im, iml, imr, mask_stencil
      CCTK_INT                             :: im, iml, mask_stencil
      CCTK_INT                             :: nx, ny, nz, nsten
      CCTK_INT                             :: stenleft(2), stenright(2)
      CCTK_INT                             :: stenbg(2), stenend(2)
      CCTK_INT, dimension(3)               :: shp 
      CCTK_REAL, pointer, dimension(:,:,:) :: mask, xx, yy, zz
      CCTK_INT,  allocatable, dimension(:,:,:) :: ighost
      CCTK_INT                             :: bbox1, bbox2, bbox3, bbox4, &
                                              bbox5, bbox6, nxgz, nygz, nzgz
      logical                              :: status
      ! This better be set to true
      logical, parameter                   :: ltrace  = .true.
      logical, parameter                   :: ltrace3 = .false.
      logical, parameter                   :: paranoid = .false.
 


      if (ltrace3) then
         write(0,*)'Entering SetImask2_KO'
      end if

      mask  => w(H_MASK)%d
      xx    => w(H_X)%d
      yy    => w(H_Y)%d
      zz    => w(H_Z)%d

      shp = shape(mask)
      nx = shp(1)
      ny = shp(2)
      nz = shp(3)

      if (ltrace) then
        allocate(ighost(nx,ny,nz),STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'SetImask2_KO: can not allocate memory for ighost '
          write(0,*)'   size ',nx,ny,nz
        end if
        ighost = 0
        nxgz = nint(par(P_NGHOSTZONES_X))
        nygz = nint(par(P_NGHOSTZONES_Y))
        nzgz = nint(par(P_NGHOSTZONES_Z))
        bbox1 = nint(par(P_BBOX1))
        bbox2 = nint(par(P_BBOX2))
        bbox3 = nint(par(P_BBOX3))
        bbox4 = nint(par(P_BBOX4))
        bbox5 = nint(par(P_BBOX5))
        bbox6 = nint(par(P_BBOX6))
        if (bbox1 == 0) then
          do k = 1, nz
          do j = 1, ny
          do i = 1, nxgz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox2 == 0) then
          do k = 1, nz
          do j = 1, ny
          do i = nx+1-nxgz, nx
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox3 == 0) then
          do k = 1, nz
          do i = 1, nx
          do j = 1, nygz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox4 == 0) then
          do k = 1, nz
          do i = 1, nx
          do j = ny+1-nygz, ny
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox5 == 0) then
          do j = 1, ny
          do i = 1, nx
          do k = 1, nzgz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox6 == 0) then
          do j = 1, ny
          do i = 1, nx
          do k = nz+1-nzgz, nz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
      end if

      if (ltrace3) then
         write(0,*)'SETIMASK2_KO: minval/maxval(mask) ', &
                   minval(w(H_MASK)%d), maxval(w(H_MASK)%d)
      end if

      ! Initialize the index mask
      !
      imask = 0

      !--------------------------------------------------------
      ! Do the outer boundaries first
      !--------------------------------------------------------


      stenbg(1) = 0
      stenend(1) = 2
      stenbg(2) = -1
      stenend(2) = 2
     
      stenleft(1)  = P_STENCIL_X_LEFT
      stenleft(2)  = P_STENCIL_X_CENTER_LEFT
      stenright(1) = P_STENCIL_X_RIGHT
      stenright(2) = P_STENCIL_X_CENTER_RIGHT           
      
      do k = 1, nz
      do j = 1, ny
        do i=1,2
	  if(mask(i+stenend(i),j,k).gt.0) then
             imask(i,j,k)    = stenleft(i)
	  else
             imask(i,j,k)    = P_STENCIL_X_UNDEF
	  end if
	     
	 q = nx+1-i
	  if(mask(q-stenend(i),j,k).gt.0) then
             imask(q,j,k)    = stenright(i)
	  else
             imask(q,j,k)    =  P_STENCIL_X_UNDEF
	  end if	  
	end do
      end do
      end do
      if(ltrace3)write(*,*)'setimask2_ko: done w/ x bounds'
      
      stenleft(1)  = P_STENCIL_Y_LEFT
      stenleft(2)  = P_STENCIL_Y_CENTER_LEFT
      stenright(1) = P_STENCIL_Y_RIGHT
      stenright(2) = P_STENCIL_Y_CENTER_RIGHT 
      do k = 1, nz
      do i = 1, nx
        do j=1,2
	  if(mask(i,j+stenend(j),k).gt.0) then
	    imask(i,j,k) =      stenleft(j) + (imask(i,j,k)/10000)*10000 &
                              + mod(imask(i,j,k),100)
	  else
	    imask(i,j,k) = P_STENCIL_Y_UNDEF+ (imask(i,j,k)/10000)*10000 &
                              + mod(imask(i,j,k),100)
	  end if
	q = ny+1-j
	  if(mask(i,q-stenend(j),k).gt.0) then
	    imask(i,q,k) =      stenright(j)+ (imask(i,q,k)/10000)*10000 &
                              + mod(imask(i,q,k),100)
	  else
	    imask(i,q,k) = P_STENCIL_Y_UNDEF+ (imask(i,q,k)/10000)*10000 &
                              + mod(imask(i,q,k),100)
	  end if	  
	end do
 
      end do
      end do
      if(ltrace3)write(*,*)'setimask2_ko: done w/ y bounds'

      stenleft(1)  = P_STENCIL_Z_LEFT
      stenleft(2)  = P_STENCIL_Z_CENTER_LEFT
      stenright(1) = P_STENCIL_Z_RIGHT
      stenright(2) = P_STENCIL_Z_CENTER_RIGHT
      do j = 1, ny
      do i = 1, nx
       do k=1,2
 	 if(mask(i,j,k+stenend(k)).gt.0) then
	    imask(i,j,k) = stenleft(k) + (imask(i,j,k)/100)*100
	  else
	    imask(i,j,k) = P_STENCIL_Z_UNDEF + (imask(i,j,k)/100)*100
	  end if	       
        q = nz+1-k
 	 if(mask(i,j,q-stenend(k)).gt.0) then
	    imask(i,j,q) = stenright(k)+ (imask(i,j,q)/100)*100
	  else
	    imask(i,j,q) = P_STENCIL_Z_UNDEF+ (imask(i,j,q)/100)*100
	  end if	
       end do
       
      end do
      end do
      if(ltrace3)write(*,*)'setimask2_ko: done w/ z bounds'


      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)&
                 'Check_imask found illegal value after outer boundary sweep'
            stop
         end if
      end if


      stenleft(1)  = P_STENCIL_X_LEFT
      stenleft(2)  = P_STENCIL_X_CENTER_LEFT
      stenright(1) = P_STENCIL_X_RIGHT
      stenright(2) = P_STENCIL_X_CENTER_RIGHT
      stenbg(1) = 0
      stenend(1) = 2
      stenbg(2) = -1
      stenend(2) = 2
     

      !--------------------------------------------------------
      ! March along x and check for mask boundaries
      !--------------------------------------------------------
      if(ltrace3)write(*,*)'setimask2_ko: marching along x...'
      do k = 1, nz
      do j = 1, ny

         do i = 2, nx
            im  = nint(mask(i,j,k))
            iml = nint(mask(i-1,j,k))

            if (im *iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  if (i < nx) then
                     nsten = 2
                  else
                     nsten = 1
                  end if
                  do m = 1, nsten
                     ii = i + m - 1
                     status = .true.
                     do q = ii+stenbg(m), ii+stenend(m)
                       if (q < 1 .or. q > nx ) then
                          status = .false.
                       else if (nint(mask(q,j,k)).eq.0) then
                          status = .false.
                       end if
                       !if (q < 1 .or. q > nx .or. nint(mask(q,j,k)) .eq. 0) then
                       !   status = .false.
                       !end if
                     end do
                     if (status) then
                        imask(ii,j,k) = stenleft(m) + mod(imask(ii,j,k),10000)
                     !else
                     else if (ii.gt.0 .and. ii.le.nx) then
                        imask(ii,j,k) = P_STENCIL_X_UNDEF + mod(imask(ii,j,k),10000)
                     else 
                        if (ltrace) write(*,*) 'Set_Imask2_KO: Lft x would write: =',ii
                     end if
                  end do

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  if (i > 2) then
                     nsten = 2
                  else
                     nsten = 1
                  end if 
                  do m = 1, nsten
                  !do m = 1, 2
                     ii = i - m
                     status = .true.
                     do q = ii-stenbg(m), ii-stenend(m),-1
                       if (q < 1 .or. q > nx ) then
                          status = .false.
                       else if (nint(mask(q,j,k)).eq.0) then
                          status = .false.
                       end if
                       !if (q < 1 .or. q > nx .or. nint(mask(q,j,k)) .eq. 0) then
                       !   status = .false.
                       !end if
                     end do
                     if (status) then
                        imask(ii,j,k) = stenright(m) + mod(imask(ii,j,k),10000)
                     !else
                     else if (ii.gt.0 .and. ii.le.nx) then
                        imask(ii,j,k) = P_STENCIL_X_UNDEF + mod(imask(ii,j,k),10000)
                     else 
                        if (ltrace) write(*,*) 'Set_Imask2_KO: Rht x would write: =',ii
                     end if

                  end do


               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in x sweep',im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do


      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after x sweep'
            stop
         end if
      end if

      stenleft(1)  = P_STENCIL_Y_LEFT
      stenleft(2)  = P_STENCIL_Y_CENTER_LEFT
      stenright(1) = P_STENCIL_Y_RIGHT
      stenright(2) = P_STENCIL_Y_CENTER_RIGHT

      !--------------------------------------------------------
      ! March along y and check for mask boundaries
      !--------------------------------------------------------
      if(ltrace3)write(*,*)'setimask2_ko: marching along y...'
      do k = 1, nz
      do i = 1, nx

         do j = 2, ny
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j-1,k))

            if (im * iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  if (j < ny) then
                     nsten = 2
                  else
                     nsten = 1
                  end if
                  do m = 1, nsten
                     ii = j + m - 1
                     status = .true.
                     do q = ii+stenbg(m), ii+stenend(m)
                       if (q < 1 .or. q > ny ) then
                          status = .false.
                       else if (nint(mask(i,q,k)).eq.0) then
                          status = .false.
                       end if
                       !if (q < 1 .or. q > ny .or. nint(mask(i,q,k)) .eq. 0) then
                       !   status = .false.
                       !end if
                     end do
                     if (status) then
                        imask(i,ii,k) = stenleft(m) &
                              + (imask(i,ii,k)/10000)*10000 &
                              + mod(imask(i,ii,k),100)
                     !else
                     else if (ii.gt.0 .and. ii.le.ny) then
                        imask(i,ii,k) = P_STENCIL_Y_UNDEF &
                              + (imask(i,ii,k)/10000)*10000 &
                              +  mod(imask(i,ii,k),100)
                        if (ighost(i,ii,k) == 0) then
                          write(0,*)'Error on left: Undef derivKO at j=',ii,ny
                          write(0,*)'iml, im ',iml,im
                          write(0,*)'mask    ',mask(i,j,k),mask(i,j-1,k)
                          write(0,*)'bbox3,bbox4  ',par(P_BBOX3),par(P_BBOX4)
                          write(0,*)'x,y,z ',xx(i,j,k),yy(i,j,k),zz(i,j,k)
                          write(0,*)'xl,yl,zl',xx(i,j-1,k),yy(i,j-1,k),zz(i,j-1,k)
                        end if
                     else 
                        if (ltrace) write(*,*) 'Set_Imask2_KO: Lft y would write: =',ii
                     end if
                  end do

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  if (j > 2) then
                     nsten = 2
                  else
                     nsten = 1
                  end if 
                  do m = 1, nsten
                     ii = j - m
                     status = .true.
                     do q = ii-stenbg(m), ii-stenend(m),-1
                       if (q < 1 .or. q > ny ) then
                          status = .false.
                       else if (nint(mask(i,q,k)).eq.0) then
                          status = .false.
                       end if
                       !if (q < 1 .or. q > ny .or. nint(mask(i,q,k)) .eq. 0) then
                       !   status = .false.
                       !end if
                     end do
                     if (status) then
                        imask(i,ii,k) = stenright(m) &
                              + (imask(i,ii,k)/10000)*10000 &
                              + mod(imask(i,ii,k),100)
                     !else
                     else if (ii.gt.0 .and. ii.le.ny) then
                        imask(i,ii,k) = P_STENCIL_Y_UNDEF &
                              + (imask(i,ii,k)/10000)*10000 &
                              + mod(imask(i,ii,k),100)
                        if (ltrace) then
                           if (ighost(i,ii,k) == 0) then
                              write(0,*) &
                                 'Error on right: Undef derivKO at j=',ii
                              write(0,*)'iml, im ',iml,im
                              write(0,*)'mask    ',mask(i,j,k),mask(i,j-1,k)
                              write(0,*)'bbox3,bbox4  ',par(P_BBOX3),&
                                                        par(P_BBOX4)
                              write(0,*)'x,y,z ',xx(i,j,k),yy(i,j,k),zz(i,j,k)
                              write(0,*)'xl,yl,zl',xx(i,j-1,k),&
                                         yy(i,j-1,k),zz(i,j-1,k)
                           end if
                        end if
                     else 
                        if (ltrace) write(*,*) 'Set_Imask2_KO: Rht y would write: =',ii
                     end if
                  end do

               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in y sweep',im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do
 
      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after y sweep'
            stop
         end if
      end if

      stenleft(1)  = P_STENCIL_Z_LEFT
      stenleft(2)  = P_STENCIL_Z_CENTER_LEFT
      stenright(1) = P_STENCIL_Z_RIGHT
      stenright(2) = P_STENCIL_Z_CENTER_RIGHT

      !--------------------------------------------------------
      ! March along z and check for mask boundaries
      !--------------------------------------------------------
      if(ltrace3)write(*,*)'setimask2_ko: marching along z...'
      do j = 1, ny
      do i = 1, nx

         do k = 2, nz
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j,k-1))

            if (im * iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  if (k < nz) then
                     nsten = 2
                  else
                     nsten = 1
                  end if
                  do m = 1, nsten 
                     ii = k + m - 1
                     status = .true.
                     do q = ii+stenbg(m), ii+stenend(m)
                       if (q < 1 .or. q > nz ) then
                          status = .false.
                       else if (nint(mask(i,j,q)).eq.0) then
                          status = .false.
                       end if
                       !if (q < 1 .or. q > nz .or. nint(mask(i,j,q)) .eq. 0) then
                       !   status = .false.
                       !end if
                     end do
                     if (status) then
                        imask(i,j,ii) = stenleft(m) + (imask(i,j,ii)/100)*100
                     !else
                     else if (ii.gt.0 .and. ii.le.nz) then
                        imask(i,j,ii) = P_STENCIL_Z_UNDEF + (imask(i,j,ii)/100)*100
                     else 
                        if (ltrace) write(*,*) 'Set_Imask2_KO: Lft z would write: =',ii
                     end if
                  end do

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  if (k > 2) then
                     nsten = 2
                  else
                     nsten = 1
                  end if 
                  do m = 1, nsten
                     ii = k - m
                     status = .true.
                     do q = ii-stenbg(m), ii-stenend(m),-1
                       if (q < 1 .or. q > nz ) then
                          status = .false.
                       else if (nint(mask(i,j,q)).eq.0) then
                          status = .false.
                       end if
                       !if (q < 1 .or. q > nz .or. nint(mask(i,j,q)) .eq. 0) then
                       !   status = .false.
                       !end if
                     end do
                     if (status) then
                        imask(i,j,ii) = stenright(m) + (imask(i,j,ii)/100)*100
                     !else
                     else if (ii.gt.0 .and. ii.le.nz) then
                        imask(i,j,ii) = P_STENCIL_Z_UNDEF + (imask(i,j,ii)/100)*100
                     else 
                        if (ltrace) write(*,*) 'Set_Imask2_KO: Rht z would write: =',ii
                     end if
                  end do

               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in z sweep',im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after z sweep'
            stop
         end if
      end if

      !
      ! Set imask values inside the mask
      !
      mask_stencil = P_STENCIL_X_MASK + P_STENCIL_Y_MASK + P_STENCIL_Z_MASK
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (mask(i,j,k) .lt. 0.1) then
            imask(i,j,k) = mask_stencil
         end if
      end do
      end do
      end do

      
      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after mask sweep'
            stop
         end if
      end if

      if (Check_imask(imask, par) == 0) then
         write(0,*)'SETIMASK2_KO: Check_imask found illegal value in mask'
         stop
      end if

      if (ltrace) then
        if (allocated(ighost)) then
          deallocate(ighost)
        else
          write(0,*)'>>> SetImask2_KO: ighost not allocated as expected!'
        end if

      end if

      return
      end subroutine SetImask2_KO


!------------------------------------------------------------------------
!
!     SetImask42
!
!     This routine generates the imask for a 5 point stencil
!     appropriate for 4th order accurate derivatives.
!
!------------------------------------------------------------------------
      subroutine SetImask42(w, imask, par)
      use GF
      implicit   none

      type(gridfunction), dimension(NW)    :: w
      CCTK_INT, dimension(:,:,:)           :: imask
      CCTK_REAL, dimension(NPAR)           :: par

      ! local vars
      CCTK_INT                             :: i, j, k, m, ii, q, istat
      CCTK_INT                             :: im, iml, mask_stencil,iml1,iml2
      CCTK_INT                             :: nx, ny, nz, nsten
      CCTK_INT                             :: stenleft(4), stenright(4)
      CCTK_INT                             :: stenbg(4), stenend(4)
      CCTK_INT, dimension(3)               :: shp 
      CCTK_REAL, pointer, dimension(:,:,:) :: mask, xx, yy, zz
      CCTK_INT,  allocatable, dimension(:,:,:) :: ighost
      CCTK_INT                             :: bbox1, bbox2, bbox3, bbox4, &
                                              bbox5, bbox6, nxgz, nygz, nzgz
      logical                              :: status
      logical, parameter                   :: ltrace  = .true.
      logical, parameter                   :: ltrace3 = .false.
      logical, parameter                   :: paranoid = .false.
 


      if (ltrace3) then
         write(0,*)'Entering SetImask42'
      end if

      mask  => w(H_MASK)%d
      xx    => w(H_X)%d
      yy    => w(H_Y)%d
      zz    => w(H_Z)%d

      shp = shape(mask)
      nx = shp(1)
      ny = shp(2)
      nz = shp(3)

      if (ltrace) then
        allocate(ighost(nx,ny,nz),STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'SetImask42: can not allocate memory for ighost '
          write(0,*)'   size ',nx,ny,nz
        end if
        ighost = 0
        nxgz = nint(par(P_NGHOSTZONES_X))
        nygz = nint(par(P_NGHOSTZONES_Y))
        nzgz = nint(par(P_NGHOSTZONES_Z))
        bbox1 = nint(par(P_BBOX1))
        bbox2 = nint(par(P_BBOX2))
        bbox3 = nint(par(P_BBOX3))
        bbox4 = nint(par(P_BBOX4))
        bbox5 = nint(par(P_BBOX5))
        bbox6 = nint(par(P_BBOX6))
        if (bbox1 == 0) then
          do k = 1, nz
          do j = 1, ny
          do i = 1, nxgz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox2 == 0) then
          do k = 1, nz
          do j = 1, ny
          do i = nx+1-nxgz, nx
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox3 == 0) then
          do k = 1, nz
          do i = 1, nx
          do j = 1, nygz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox4 == 0) then
          do k = 1, nz
          do i = 1, nx
          do j = ny+1-nygz, ny
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox5 == 0) then
          do j = 1, ny
          do i = 1, nx
          do k = 1, nzgz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
        if (bbox6 == 0) then
          do j = 1, ny
          do i = 1, nx
          do k = nz+1-nzgz, nz
            ighost(i,j,k) = 1
          end do
          end do
          end do
        end if
      end if

      if (ltrace3) then
         write(0,*)'SETIMASK42: minval/maxval(mask) ', &
                   minval(w(H_MASK)%d), maxval(w(H_MASK)%d)
      end if

      ! Initialize the index mask
      !
      imask = 0

      stenbg(1) = 0
      stenend(1) = 3
      stenbg(2) = -1
      stenend(2) = 1
      stenbg(3) = -2
      stenend(3) = 2
      stenbg(4) = -3
      stenend(4) = 2 
      
      stenleft(1)  = P_STENCIL_X_LEFT
      stenleft(2)  = P_STENCIL_X_CENTER_LEFT
      stenleft(3)  = P_STENCIL_X_CENTER_LEFT1
      stenleft(4)  = P_STENCIL_X_CENTER_LEFT2
    
      stenright(1) = P_STENCIL_X_RIGHT
      stenright(2) = P_STENCIL_X_CENTER_RIGHT
      stenright(3) = P_STENCIL_X_CENTER_RIGHT1
      stenright(4) = P_STENCIL_X_CENTER_RIGHT2 
            
 
      do k = 1, nz
      do j = 1, ny
        do i=1,4
	  if(mask(i+stenend(i),j,k).gt.0) then
             imask(i,j,k)    = stenleft(i)
	  else
             imask(i,j,k)    = P_STENCIL_X_UNDEF
	  end if
	     
	 q = nx+1-i
	  if(mask(q-stenend(i),j,k).gt.0) then
             imask(q,j,k)    = stenright(i)
	  else
             imask(q,j,k)    =  P_STENCIL_X_UNDEF
	  end if	  
	end do
      end do
      end do
      if(ltrace3)write(*,*)'setimask42_ko: done w/ x bounds'
      
      stenleft(1)  = P_STENCIL_Y_LEFT
      stenleft(2)  = P_STENCIL_Y_CENTER_LEFT
      stenleft(3)  = P_STENCIL_Y_CENTER_LEFT1
      stenleft(4)  = P_STENCIL_Y_CENTER_LEFT2
      
      stenright(1) = P_STENCIL_Y_RIGHT
      stenright(2) = P_STENCIL_Y_CENTER_RIGHT
      stenright(3) = P_STENCIL_Y_CENTER_RIGHT1
      stenright(4) = P_STENCIL_Y_CENTER_RIGHT2  
      do k = 1, nz
      do i = 1, nx
        do j=1,4
	  if(mask(i,j+stenend(j),k).gt.0) then
	    imask(i,j,k) =      stenleft(j) + (imask(i,j,k)/10000)*10000 &
                              + mod(imask(i,j,k),100)
	  else
	    imask(i,j,k) = P_STENCIL_Y_UNDEF+ (imask(i,j,k)/10000)*10000 &
                              + mod(imask(i,j,k),100)
	  end if
	q = ny+1-j
	  if(mask(i,q-stenend(j),k).gt.0) then
	    imask(i,q,k) =      stenright(j)+ (imask(i,q,k)/10000)*10000 &
                              + mod(imask(i,q,k),100)
	  else
	    imask(i,q,k) = P_STENCIL_Y_UNDEF+ (imask(i,q,k)/10000)*10000 &
                              + mod(imask(i,q,k),100)
	  end if	  
	end do
 
      end do
      end do
      if(ltrace3)write(*,*)'setimask42_ko: done w/ y bounds'

      stenleft(1)  = P_STENCIL_Z_LEFT
      stenleft(2)  = P_STENCIL_Z_CENTER_LEFT
      stenleft(3)  = P_STENCIL_Z_CENTER_LEFT1
      stenleft(4)  = P_STENCIL_Z_CENTER_LEFT2
    
      stenright(1) = P_STENCIL_Z_RIGHT
      stenright(2) = P_STENCIL_Z_CENTER_RIGHT
      stenright(3) = P_STENCIL_Z_CENTER_RIGHT1
      stenright(4) = P_STENCIL_Z_CENTER_RIGHT2 

      do j = 1, ny
      do i = 1, nx
       do k=1,4
 	 if(mask(i,j,k+stenend(k)).gt.0) then
	    imask(i,j,k) = stenleft(k) + (imask(i,j,k)/100)*100
	  else
	    imask(i,j,k) = P_STENCIL_Z_UNDEF + (imask(i,j,k)/100)*100
	  end if	       
        q = nz+1-k
 	 if(mask(i,j,q-stenend(k)).gt.0) then
	    imask(i,j,q) = stenright(k)+ (imask(i,j,q)/100)*100
	  else
	    imask(i,j,q) = P_STENCIL_Z_UNDEF+ (imask(i,j,q)/100)*100
	  end if	
       end do
       
      end do
      end do
      if(ltrace3)write(*,*)'setimask42_ko: done w/ z bounds'
      
      
      
      
      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)&
                 'Check_imask found illegal value after outer boundary sweep'
            stop
         end if
      end if

      stenleft(1)  = P_STENCIL_X_LEFT
      stenleft(2)  = P_STENCIL_X_CENTER_LEFT
      stenleft(3)  = P_STENCIL_X_CENTER_LEFT1
      stenleft(4)  = P_STENCIL_X_CENTER_LEFT2
      
    
      stenright(1) = P_STENCIL_X_RIGHT
      stenright(2) = P_STENCIL_X_CENTER_RIGHT
      stenright(3) = P_STENCIL_X_CENTER_RIGHT1
      stenright(4) = P_STENCIL_X_CENTER_RIGHT2 
            
      stenbg(1) = 0
      stenend(1) = 3
      stenbg(2) = -1
      stenend(2) = 1
      stenbg(3) = -2
      stenend(3) = 2
      stenbg(4) = -3
      stenend(4) = 2            
     

      !--------------------------------------------------------
      ! March along x and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do j = 1, ny

         do i = 2, nx
            im  = nint(mask(i,j,k))
            iml = nint(mask(i-1,j,k))

            if (im *iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  if (i < nx-3) then
                     nsten = 4
                  else
                     nsten = 4
                  end if
                  do m = 1, nsten
                     ii = i + m - 1
                     status = .true.
                     do q = ii+stenbg(m), ii+stenend(m)
                       if (q < 1 .or. q > nx) then
                          status = .false.
                       else if(nint(mask(q,j,k)).eq.0) then
                          status = .false.
                       end if
                     end do
                     if (status) then
                        imask(ii,j,k) =       stenleft(m) + mod(imask(ii,j,k),10000) 
                     else if(ii.ge.1.and.ii.le.nx) then
                        imask(ii,j,k) = P_STENCIL_X_UNDEF + mod(imask(ii,j,k),10000)
                     end if
                  end do

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  if (i > 4) then
                     nsten = 4
                  else
                     nsten = 1
                  end if 
                  do m = 1, nsten
                     ii = i - m
                     status = .true.
                     do q = ii-stenbg(m), ii-stenend(m),-1
                       if (q < 1 .or. q > nx) then
                          status = .false.
                       else if(nint(mask(q,j,k)).eq.0) then
                          status = .false.
                       end if
                     end do
                     if (status) then
                        imask(ii,j,k) =      stenright(m) + mod(imask(ii,j,k),10000) 
                     else if (ii.ge.1.and.ii.le.nx) then
                        imask(ii,j,k) = P_STENCIL_X_UNDEF + mod(imask(ii,j,k),10000)
                     end if
                  end do


               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in x sweep',im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after x sweep'
            stop
         end if
      end if

      stenleft(1)  = P_STENCIL_Y_LEFT
      stenleft(2)  = P_STENCIL_Y_CENTER_LEFT
      stenleft(3)  = P_STENCIL_Y_CENTER_LEFT1
      stenleft(4)  = P_STENCIL_Y_CENTER_LEFT2      
      stenright(1) = P_STENCIL_Y_RIGHT
      stenright(2) = P_STENCIL_Y_CENTER_RIGHT
      stenright(3) = P_STENCIL_Y_CENTER_RIGHT1
      stenright(4) = P_STENCIL_Y_CENTER_RIGHT2

      !--------------------------------------------------------
      ! March along y and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do i = 1, nx

         do j = 2, ny
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j-1,k))
            iml1 = nint(mask(i,j-2,k))
            iml2 = nint(mask(i,j-3,k))

            if (im * iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  if (j < ny-3) then
                     nsten = 4
                  else
                     nsten = 4
                  end if
                  do m = 1, nsten
                     ii = j + m - 1
                     status = .true.
                     do q = ii+stenbg(m), ii+stenend(m)
                       if (q < 1 .or. q > ny ) then
                          status = .false.
                       else if(nint(mask(i,q,k)).eq.0) then
                          status = .false.
                       end if
                     end do
                     if (status) then
                        imask(i,ii,k) = stenleft(m) &
                                           + (imask(i,ii,k)/10000)*10000 &
			                   + mod(imask(i,ii,k),100) 
                     else if(ii.ge.1.and.ii.le.ny) then
                        imask(i,ii,k) = P_STENCIL_Y_UNDEF &
                                           + (imask(i,ii,k)/10000)*10000 &
                                           + mod(imask(i,ii,k),100)
                        if (ighost(i,ii,k) == 0) then
                          write(0,*)'Error on left: Undefined deriv at j=',ii,ny
                          write(0,*)'iml, im ',iml,im
                          write(0,*)'mask    ',mask(i,j,k),mask(i,j-1,k)
                          write(0,*)'bbox3,bbox4  ',par(P_BBOX3),par(P_BBOX4)
                          write(0,*)'x,y,z ',xx(i,j,k),yy(i,j,k),zz(i,j,k)
                          write(0,*)'xl,yl,zl',xx(i,j-1,k),yy(i,j-1,k),zz(i,j-1,k)
                       end if
                     end if
                  end do

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  if (j > 4) then
                     nsten = 4
                  else
                     nsten = 4
                  end if 
                  do m = 1, nsten
                     ii = j - m
                     status = .true.
                     do q = ii-stenbg(m), ii-stenend(m),-1
                       if (q < 1 .or. q > ny ) then
                          status = .false.
                       else if(nint(mask(i,q,k)).eq.0) then
                          status = .false.
                       end if
                     end do
                     if (status) then
                        imask(i,ii,k) = stenright(m)&
                                           + (imask(i,ii,k)/10000)*10000 &
                                           + mod(imask(i,ii,k),100)
                     else if(ii.ge.1.and.ii.le.ny) then 
                        imask(i,ii,k) = P_STENCIL_Y_UNDEF &
                                           + (imask(i,ii,k)/10000)*10000 &
                                           + mod(imask(i,ii,k),100)
                        if (ltrace) then
                           if (ighost(i,ii,k) == 0) then
                              write(0,*) &
                                 'Error on right: Undefined deriv at j=',ii
                              write(0,*)'iml, im ',iml,im
                              write(0,*)'mask    ',mask(i,j,k),mask(i,j-1,k)
                              write(0,*)'bbox3,bbox4  ',par(P_BBOX3),&
                                                        par(P_BBOX4)
                              write(0,*)'x,y,z ',xx(i,j,k),yy(i,j,k),zz(i,j,k)
                              write(0,*)'xl,yl,zl',xx(i,j-1,k),&
                                         yy(i,j-1,k),zz(i,j-1,k)
                           end if
                        end if

                     end if
                  end do

               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in y sweep',im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do
 
      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after y sweep'
            stop
         end if
      end if

      stenleft(1)  = P_STENCIL_Z_LEFT
      stenleft(2)  = P_STENCIL_Z_CENTER_LEFT
      stenleft(3)  = P_STENCIL_Z_CENTER_LEFT1
      stenleft(4)  = P_STENCIL_Z_CENTER_LEFT2

      stenright(1) = P_STENCIL_Z_RIGHT
      stenright(2) = P_STENCIL_Z_CENTER_RIGHT
      stenright(3) = P_STENCIL_Z_CENTER_RIGHT1
      stenright(4) = P_STENCIL_Z_CENTER_RIGHT2

      !--------------------------------------------------------
      ! March along z and check for mask boundaries
      !--------------------------------------------------------
      do j = 1, ny
      do i = 1, nx

         do k = 2, nz
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j,k-1))

            if (im * iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  if (k < nz-3) then
                     nsten = 4
                  else
                     nsten = 1
                  end if
                  do m = 1, nsten 
                     ii = k + m - 1
                     status = .true.
                     do q = ii+stenbg(m), ii+stenend(m)
                       if (q < 1 .or. q > nz ) then
                          status = .false.
                       else if (nint(mask(i,j,q)).eq.0) then
                          status = .false.
                       end if
                     end do
                     if (status) then
                        imask(i,j,ii) =       stenleft(m) + (imask(i,j,ii)/100)*100
                     else if(ii.ge.1.and.ii.le.nz) then
                        imask(i,j,ii) = P_STENCIL_Z_UNDEF + (imask(i,j,ii)/100)*100
                     end if
                  end do

               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  if (k > 4) then
                     nsten = 4
                  else
                     nsten = 1
                  end if 
                  do m = 1, nsten
                     ii = k - m
                     status = .true.
                     do q = ii-stenbg(m), ii-stenend(m),-1
                       if (q < 1 .or. q > nz ) then
                          status = .false.
                       else if (nint(mask(i,j,q)).eq.0) then
                          status = .false.
                       end if
                     end do
                     if (status) then
                        imask(i,j,ii) =      stenright(m) + (imask(i,j,ii)/100)*100
                     else if(ii.ge.1.and.ii.le.nz) then
                        imask(i,j,ii) = P_STENCIL_Z_UNDEF + (imask(i,j,ii)/100)*100
                     end if
                  end do

               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'mask problem in z sweep',im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after z sweep'
            stop
         end if
      end if

      !
      ! Set imask values inside the mask
      !
      mask_stencil = P_STENCIL_X_MASK + P_STENCIL_Y_MASK + P_STENCIL_Z_MASK
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (mask(i,j,k) .lt. 0.1) then
            imask(i,j,k) = mask_stencil
         end if
      end do
      end do
      end do

      
      if (paranoid) then
         if (Check_imask(imask, par) == 0) then
            write(0,*)'Check_imask found illegal value after mask sweep'
            stop
         end if
      end if

!      if (Check_imask(imask, par) == 0) then
!         write(0,*)'SETIMASK2_KO: Check_imask found illegal value in mask'
!         stop
!      end if

      if (ltrace) then
        if (allocated(ighost)) then
          deallocate(ighost)
        else
          write(0,*)'>>> SetImask2_KO: ighost not allocated as expected!'
        end if

      end if

      return
      end subroutine SetImask42




!------------------------------------------------------------------------
!
!     Check_imask
!
!     This routine is for debugging purposes.  It goes through imask
!     and checks for illegal values.
!
!------------------------------------------------------------------------
      integer function Check_imask(imask, par)
      implicit   none

      CCTK_INT, dimension(:,:,:)           :: imask
      CCTK_REAL, dimension(NPAR)           :: par

      ! local vars
      CCTK_INT                             :: i, j, k, m, nx, ny, nz
      CCTK_INT, dimension(3)               :: shp
      CCTK_INT                             :: dx, dy, dz
      CCTK_INT  :: gz_xmin, gz_xmax, gz_ymin, gz_ymax, gz_zmin, gz_zmax
      CCTK_INT  :: bbox1, bbox2, bbox3, bbox4, bbox5, bbox6

      logical                              :: fx, fy, fz
      logical, parameter                   :: ltrace = .false.

      Check_imask = 1

      shp = shape(imask)
      if (nint(par(P_NX)) .ne. shp(1) .or. nint(par(P_NY)) .ne. shp(2)&
          .or. nint(par(P_NZ)) .ne. shp(3)) then
         write(0,*)'Check_imask: imask has wrong shape'
         write(0,*)'Check_imask: shape of imask ', shp(1), shp(2), shp(3)
         write(0,*)'Check_imask: nx,ny,nz       ',nx, ny, nz
        
      end if


      nx = shp(1)
      ny = shp(2)
      nz = shp(3)

      gz_xmin = 1
      gz_xmax = nx
      gz_ymin = 1
      gz_ymax = ny
      gz_zmin = 1
      gz_zmax = nz

      bbox1 = nint(par(P_BBOX1))
      bbox2 = nint(par(P_BBOX2))
      bbox3 = nint(par(P_BBOX3))
      bbox4 = nint(par(P_BBOX4))
      bbox5 = nint(par(P_BBOX5))
      bbox6 = nint(par(P_BBOX6))

      if (bbox1 .eq. 0) then
        gz_xmin = gz_xmin + nint(par(P_NGHOSTZONES_X))
      end if
      if (bbox2 .eq. 0) then
        gz_xmax = gz_xmax - nint(par(P_NGHOSTZONES_X))
      end if
  
      if (bbox3 .eq. 0) then
        gz_ymin = gz_ymin + nint(par(P_NGHOSTZONES_Y))
      end if
      if (bbox4 .eq. 0) then
        gz_ymax = gz_ymax - nint(par(P_NGHOSTZONES_Y))
      end if
  
      if (bbox5 .eq. 0) then
        gz_zmin = gz_zmin + nint(par(P_NGHOSTZONES_Z))
      end if
      if (bbox6 .eq. 0) then
        gz_zmax = gz_zmax - nint(par(P_NGHOSTZONES_Z))
      end if

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         dx = imask(i,j,k) / 10000
         dy = mod(imask(i,j,k)/100,100)
         dz = mod(imask(i,j,k),100)
         fx = .false.
         fy = .false.
         fz = .false.
         do m = 0, P_STENCIL_UNDEF
           if (dx == m) then
             fx = .true.
           end if
         end do
         do m = 0, P_STENCIL_UNDEF
           if (dy == m) then
             fy = .true.
           end if
         end do
         do m = 0, P_STENCIL_UNDEF
           if (dz == m) then
             fz = .true.
           end if
         end do
         if (fx.eqv..false. .or. fy.eqv..false. .or. fz.eqv..false.) then
           write(0,*) 'Check_imask:: Illegal value found at ',i,j,k
           write(0,*) '              imask                  ',imask(i,j,k)
           write(0,*) '              dx                     ',dx
           write(0,*) '              dy                     ',dy
           write(0,*) '              dz                     ',dz
           Check_imask = 0
         end if
         ! Now check that undefined derivatives are in ghost zones
         if (dx == P_STENCIL_UNDEF) then
           if ((bbox1==0 .and. i < gz_xmin) .or.  &
               (bbox2==0 .and. i > gz_xmax)) then
              if (ltrace) then
              write(0,*)'Check_imask: OK---Undefined x derivative in ghostzone'
              end if
           else
              write(0,*)'Check_imask: ### Undefined x derivative outside of ghost zone!'
              write(0,*)'Check_imask: imask                          ',imask(i,j,k)
              write(0,*)'Check_imask: i, nx                          ',i, nx
              write(0,*)'Check_imask: gz_xmin                        ',gz_xmin
              write(0,*)'Check_imask: gz_xmax                        ',gz_xmax
              write(0,*)'Check_imask: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_X)
              write(0,*)'Check_imask: par(bbox1), par(bbox2)         ',&
                                                par(P_BBOX1), par(P_BBOX2)
              write(0,*)'Check_imask: bbox1, bbox2                   ',&
                                                bbox1,bbox2
              Check_imask = 0
           end if
         end if
         if (dy == P_STENCIL_UNDEF) then
            if ((bbox3==0 .and. j < gz_ymin) .or. &
                (bbox4==0 .and. j > gz_ymax)) then
              if (ltrace) then
              write(0,*)'Check_imask: OK---Undefined y derivative in ghostzone'
              end if
            else
              write(0,*)'Check_imask: ### Undefined y derivative outside of ghost zone!'
              write(0,*)'Check_imask: imask                          ',imask(i,j,k)
              write(0,*)'Check_imask: j, ny                          ',j, ny
              write(0,*)'Check_imask: gz_ymin                        ',gz_ymin
              write(0,*)'Check_imask: gz_ymax                        ',gz_ymax
              write(0,*)'Check_imask: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_Y)
              write(0,*)'Check_imask: par(bbox3), par(bbox4)         ',&
                                                par(P_BBOX3), par(P_BBOX4)
              write(0,*)'Check_imask: bbox3, bbox4                   ',&
                                                bbox3,bbox4
              Check_imask = 0
            end if

         end if
         if (dz == P_STENCIL_UNDEF) then
            if ((bbox5==0 .and. k < gz_zmin) .or. &
                (bbox6==0 .and. k > gz_zmax)) then
              if (ltrace) then
              write(0,*)'Check_imask: OK---Undefined z derivative in ghostzone'
              end if
            else
              write(0,*)'Check_imask: ### Undefined z derivative outside of ghost zone!'
              write(0,*)'Check_imask: imask                          ',imask(i,j,k)
              write(0,*)'Check_imask: k, nz                          ',k, nz
              write(0,*)'Check_imask: gz_zmin                        ',gz_zmin
              write(0,*)'Check_imask: gz_zmax                        ',gz_zmax
              write(0,*)'Check_imask: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_Z)
              write(0,*)'Check_imask: par(bbox5), par(bbox6)         ',&
                                                par(P_BBOX5), par(P_BBOX6)
              write(0,*)'Check_imask: bbox5, bbox6                   ',&
                                                bbox5,bbox6
              Check_imask = 0
            end if

         end if
      
      end do
      end do
      end do

100  CONTINUE

      return
      end function Check_imask


!------------------------------------------------------------------------
!
!     
!
!------------------------------------------------------------------------
      subroutine fixup_mask(mask, par)
      implicit   none

      CCTK_REAL, dimension(:,:,:)     :: mask
      CCTK_REAL, dimension(NPAR)      :: par

      !-------- local vars --------------
      CCTK_INT                        :: nx, ny, nz
      CCTK_INT                        :: i,  j,  k, im, iml
      logical, parameter              :: ltrace2=.false.

      ! For now have this not do anything.
      ! Instead, let code in had/src/amr take care of this
      !    because I am trying to have the mask move, and this
      !    routine keeps writing over what I am trying to do.
      !    We can work out later where this should go once and
      !    for all (SLL 12/2/05).
      return

      nx = nint(par(P_NX))
      ny = nint(par(P_NY))
      nz = nint(par(P_NZ))


      !--------------------------------------------------------
      ! March along x and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do j = 1, ny

         do i = 2, nx
            im  = nint(mask(i,j,k))
            iml = nint(mask(i-1,j,k))

            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  mask(i-1,j,k) = mask(i-1,j,k) - 1.0
                  if (ltrace2) then
                     write(0,*)'Found boundary in x ',mask(i-1,j,k)
                  end if
               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  mask(i,j,k) = mask(i,j,k) - 1.0
                  if (ltrace2) then
                     write(0,*)'Found boundary in x ',mask(i,j,k)
                  end if
               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'FIXUP_MASK >> mask problem in x sweep',&
                                          im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      !--------------------------------------------------------
      ! March along y and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do i = 1, nx

         do j = 2, ny
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j-1,k))

            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  mask(i,j-1,k) = mask(i,j-1,k) - 1.0
               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  mask(i,j,k) = mask(i,j,k) - 1.0
               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'FIXUP_MASK >> mask problem in y sweep',&
                                          im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      !--------------------------------------------------------
      ! March along y and check for mask boundaries
      !--------------------------------------------------------
      do j = 1, ny
      do i = 1, nx

         do k = 2, nz
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j,k-1))

            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  mask(i,j,k-1) = mask(i,j,k-1) - 1.0
               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  mask(i,j,k) = mask(i,j,k) - 1.0
               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'FIXUP_MASK >> mask problem in z sweep',&
                                          im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      return
      end subroutine fixup_mask



!------------------------------------------------------------------------
!
!     
!
!------------------------------------------------------------------------
      subroutine fixup_imask(imask, mask, par)
      implicit   none

      CCTK_INT,  dimension(:,:,:)       :: imask
      CCTK_REAL, dimension(:,:,:)       :: mask
      CCTK_REAL, dimension(:)           :: par

      !
      ! local vars
      !
      CCTK_INT                         :: i, j, k, nx, ny, nz
!     CCTK_INT                         :: ii, jj, kk, ib, ie, jb, je, kb, ke
      CCTK_INT                         :: ib, ie, jb, je, kb, ke
!     CCTK_INT                         :: bbox1, bbox2, bbox3
!     CCTK_INT                         :: bbox4, bbox5, bbox6
!     CCTK_INT                         :: xp, xm, yp, ym, zp, zm, iq, jq, kq
      CCTK_INT                         :: xp, xm, yp, ym, zp, zm
      CCTK_INT                         :: fail, dissipation
      logical                          :: modify_for_KO


      nx = nint(par(P_NX))
      ny = nint(par(P_NY))
      nz = nint(par(P_NZ))

      dissipation = nint(par(P_DISSIPATION))

      if (dissipation .eq. 2) then
         modify_for_KO = .true.
      else
         modify_for_KO = .false.
      end if

      if (modify_for_KO) then
         ib = 3
         ie = nx-2
         jb = 3
         je = ny-2
         kb = 3
         ke = nz-2
      else
         ib = 2
         ie = nx-1
         jb = 2
         je = ny-1
         kb = 2
         ke = nz-1
      end if

      !------------------------------------------------------
      ! make sure that index bounds are properly set.
      ! the stencil size depends on whether we are using 
      ! KO dissipation
      !------------------------------------------------------
      fail = 0
      if (modify_for_KO) then
         if (ib < 3 .or. ie > nx-2) then
            write(0,*)'FIXUP_IMASK >> problem here in edge x sweep'
            fail = 1
         end if
         if (jb < 3 .or. je > ny-2) then
            write(0,*)'FIXUP_IMASK >> problem here in edge y sweep'
            fail = 1
         end if
         if (kb < 3 .or. ke > nz-2) then
            write(0,*)'FIXUP_IMASK >> problem here in edge z sweep'
            fail = 1
         end if
      else
         if (ib < 2 .or. ie > nx-1) then
            write(0,*)'FIXUP_IMASK >> problem here in edge x sweep'
            fail = 1
         end if
         if (jb < 2 .or. je > ny-1) then
            write(0,*)'FIXUP_IMASK >> problem here in edge y sweep'
            fail = 1
         end if
         if (kb < 2 .or. ke > nz-1) then
            write(0,*)'FIXUP_IMASK >> problem here in edge z sweep'
            fail = 1
         end if
      end if
      if (fail==1) then
         !error message was printed earlier
         write(0,*) 'FIXUP_IMASK >> Mask is too close to a physical'
         write(0,*) '            >> boundary, or there are not enough'
         write(0,*) '            >> MPI ghost zones.'
         stop
      end if



      do k = kb, ke
      do j = jb, je
      do i = ib, ie
         ! 
         ! first check to see if mask(i,j,k) is less than -1.  This is done
         ! done for speed, as only a small number of points need the more
         ! detailed analysis.
         if (mask(i,j,k) < -1.0) then 

         ! Check to see if we are on an edge or vertex.
         if (nint(mask(i,j,k)) == P_INT_INSIDE_MASK-2) then
            !--------------------------------------------
            ! Edges
            !--------------------------------------------
            xm = 0
            xp = 0
            ym = 0
            yp = 0
            zm = 0
            zp = 0

            if (nint(mask(i+1,j,k)) > 0) xp = 1
            if (nint(mask(i-1,j,k)) > 0) xm = 1
            if (nint(mask(i,j+1,k)) > 0) yp = 1
            if (nint(mask(i,j-1,k)) > 0) ym = 1
            if (nint(mask(i,j,k+1)) > 0) zp = 1
            if (nint(mask(i,j,k-1)) > 0) zm = 1

            !----- Now figure out what edge we are on.  Painful but
            !      necessary for proper differencing with 
            !      Kreiss-Oliger dissipation.

            if (xp==1 .and. yp==1) then
               !
               !     o  o  l  o 
               !     o  R  C  L         H = here, the point i,j (within mask)
               !     x  H  r  o         C = the corner of the evolved points
               !     x  x  o  o         r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXY(imask(i+1,j+1,k), P_STENCIL_X_EDGE_P, P_STENCIL_Y_EDGE_P) 
               ! These are modified for KO
               if (modify_for_KO) then
                  IMASK_MODY(imask(i+1,j,k),   P_STENCIL_Y_CENTER_RIGHT)
                  IMASK_MODY(imask(i+1,j+2,k), P_STENCIL_Y_CENTER_LEFT)
                  IMASK_MODX(imask(i,j+1,k),   P_STENCIL_X_CENTER_RIGHT)
                  IMASK_MODX(imask(i+2,j+1,k), P_STENCIL_X_CENTER_LEFT)
               end if
            end if

            if (xp==1 .and. ym==1) then
               !
               !     x  x  o  o
               !     x  H  l  o       H = here, the point i,j (within mask)
               !     o  R  C  L       C = the corner of the evolved points
               !     o  o  r  o       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXY(imask(i+1,j-1,k), P_STENCIL_X_EDGE_P, P_STENCIL_Y_EDGE_M)
               if (modify_for_KO) then
                  IMASK_MODY(imask(i+1,j,k),   P_STENCIL_Y_CENTER_LEFT)
                  IMASK_MODY(imask(i+1,j-2,k), P_STENCIL_Y_CENTER_RIGHT)
                  IMASK_MODX(imask(i,j-1,k),   P_STENCIL_X_CENTER_RIGHT)
                  IMASK_MODX(imask(i+2,j-1,k), P_STENCIL_X_CENTER_LEFT)
               end if
            end if


            if (xm==1 .and. yp==1) then
               !
               !     o  l  o  o
               !     R  C  L  o       H = here, the point i,j (within mask)
               !     o  r  H  x       C = the corner of the evolved points
               !     o  o  x  x       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXY(imask(i-1,j+1,k), P_STENCIL_X_EDGE_M, P_STENCIL_Y_EDGE_P)
               if (modify_for_KO) then
                  IMASK_MODY(imask(i-1,j,k),   P_STENCIL_Y_CENTER_RIGHT)
                  IMASK_MODY(imask(i-1,j+2,k), P_STENCIL_Y_CENTER_LEFT)
                  IMASK_MODX(imask(i,j+1,k),   P_STENCIL_X_CENTER_LEFT)
                  IMASK_MODX(imask(i-2,j+1,k), P_STENCIL_X_CENTER_RIGHT)
               end if
            end if

            if (xm==1 .and. ym==1) then
               !
               !     o  o  x  x
               !     o  l  H  x       H = here, the point i,j (within mask)
               !     R  C  L  o       C = the corner of the evolved points
               !     o  r  o  o       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXY(imask(i-1,j-1,k), P_STENCIL_X_EDGE_M, P_STENCIL_Y_EDGE_M)
               if (modify_for_KO) then
                  IMASK_MODY(imask(i-1,j,k),   P_STENCIL_Y_CENTER_LEFT)
                  IMASK_MODY(imask(i-1,j-2,k), P_STENCIL_Y_CENTER_RIGHT)
                  IMASK_MODX(imask(i,j-1,k),   P_STENCIL_X_CENTER_LEFT)
                  IMASK_MODX(imask(i-2,j-1,k), P_STENCIL_X_CENTER_RIGHT)
               end if
            end if

            if (xp==1 .and. zp==1) then
               !
               !     o  o  l  o
               !     o  R  C  L         H = here, the point i,j (within mask)
               !     x  H  r  o         C = the corner of the evolved points
               !     x  x  o  o         r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXZ(imask(i+1,j,k+1), P_STENCIL_X_EDGE_P, P_STENCIL_Z_EDGE_P)
               if (modify_for_KO) then
                  ! These are modified for KO
                  IMASK_MODZ(imask(i+1,j,k),   P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODZ(imask(i+1,j,k+2), P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODX(imask(i,j,k+1),   P_STENCIL_X_CENTER_RIGHT)
                  IMASK_MODX(imask(i+2,j,k+1), P_STENCIL_X_CENTER_LEFT)
               end if
            end if

            if (xp==1 .and. zm==1) then
               !
               !     x  x  o  o
               !     x  H  l  o       H = here, the point i,j (within mask)
               !     o  R  C  L       C = the corner of the evolved points
               !     o  o  r  o       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXZ(imask(i+1,j,k-1), P_STENCIL_X_EDGE_P, P_STENCIL_Z_EDGE_M)
               if (modify_for_KO) then
                  IMASK_MODZ(imask(i+1,j,k),   P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODZ(imask(i+1,j,k-2), P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODX(imask(i,j,k-1),   P_STENCIL_X_CENTER_RIGHT)
                  IMASK_MODX(imask(i+2,j,k-1), P_STENCIL_X_CENTER_LEFT)
               end if
            end if

            if (xm==1 .and. zp==1) then
               !
               !     o  l  o  o
               !     R  C  L  o       H = here, the point i,j (within mask)
               !     o  r  H  x       C = the corner of the evolved points
               !     o  o  x  x       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXZ(imask(i-1,j,k+1), P_STENCIL_X_EDGE_M, P_STENCIL_Z_EDGE_P)
               if (modify_for_KO) then
                  IMASK_MODZ(imask(i-1,j,k),   P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODZ(imask(i-1,j,k+2), P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODX(imask(i,j,k+1),   P_STENCIL_X_CENTER_LEFT)
                  IMASK_MODX(imask(i-2,j,k+1), P_STENCIL_X_CENTER_RIGHT)
               end if
            end if

            if (xm==1 .and. zm==1) then
               !
               !     o  o  x  x
               !     o  l  H  x       H = here, the point i,j (within mask)
               !     R  C  L  o       C = the corner of the evolved points
               !     o  r  o  o       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODXZ(imask(i-1,j,k-1), P_STENCIL_X_EDGE_M, P_STENCIL_Z_EDGE_M)
               if (modify_for_KO) then
                  IMASK_MODZ(imask(i-1,j,k),   P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODZ(imask(i-1,j,k-2), P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODX(imask(i,j,k-1),   P_STENCIL_X_CENTER_LEFT)
                  IMASK_MODX(imask(i-2,j,k-1), P_STENCIL_X_CENTER_RIGHT)
               end if
            end if

            if (yp==1 .and. zp==1) then
               !
               !     o  o  l  o
               !     o  R  C  L         H = here, the point i,j (within mask)
               !     x  H  r  o         C = the corner of the evolved points
               !     x  x  o  o         r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODYZ(imask(i,j+1,k+1), P_STENCIL_Y_EDGE_P, P_STENCIL_Z_EDGE_P)
               if (modify_for_KO) then
                  ! These are modified for KO
                  IMASK_MODZ(imask(i,j+1,k),   P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODZ(imask(i,j+1,k+2), P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODY(imask(i,j,k+1),   P_STENCIL_Y_CENTER_RIGHT)
                  IMASK_MODY(imask(i,j+2,k+1), P_STENCIL_Y_CENTER_LEFT)
               end if
            end if

            if (yp==1 .and. zm==1) then
               !
               !     x  x  o  o
               !     x  H  l  o       H = here, the point i,j (within mask)
               !     o  R  C  L       C = the corner of the evolved points
               !     o  o  r  o       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODYZ(imask(i,j+1,k-1), P_STENCIL_Y_EDGE_P, P_STENCIL_Z_EDGE_M)
               if (modify_for_KO) then
                  IMASK_MODZ(imask(i,j+1,k),   P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODZ(imask(i,j+1,k-2), P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODY(imask(i,j,k-1),   P_STENCIL_Y_CENTER_RIGHT)
                  IMASK_MODY(imask(i,j+2,k-1), P_STENCIL_Y_CENTER_LEFT)
               end if
            end if

            if (ym==1 .and. zp==1) then
               !
               !     o  l  o  o
               !     R  C  L  o       H = here, the point i,j (within mask)
               !     o  r  H  x       C = the corner of the evolved points
               !     o  o  x  x       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODYZ(imask(i,j-1,k+1), P_STENCIL_Y_EDGE_M, P_STENCIL_Z_EDGE_P)
               if (modify_for_KO) then
                  IMASK_MODZ(imask(i,j-1,k),   P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODZ(imask(i,j-1,k+2), P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODY(imask(i,j,k+1),   P_STENCIL_Y_CENTER_LEFT)
                  IMASK_MODY(imask(i,j-2,k+1), P_STENCIL_Y_CENTER_RIGHT)
               end if
            end if

            if (ym==1 .and. zm==1) then
               !
               !     o  o  x  x
               !     o  l  H  x       H = here, the point i,j (within mask)
               !     R  C  L  o       C = the corner of the evolved points
               !     o  r  o  o       r,l,R,L = points to modify for KO diss.
               !
               IMASK_MODYZ(imask(i,j-1,k-1), P_STENCIL_Y_EDGE_M, P_STENCIL_Z_EDGE_M)
               if (modify_for_KO) then
                  IMASK_MODZ(imask(i,j-1,k),   P_STENCIL_Z_CENTER_LEFT)
                  IMASK_MODZ(imask(i,j-1,k-2), P_STENCIL_Z_CENTER_RIGHT)
                  IMASK_MODY(imask(i,j,k-1),   P_STENCIL_Y_CENTER_LEFT)
                  IMASK_MODY(imask(i,j-2,k-1), P_STENCIL_Y_CENTER_RIGHT)
               end if
            end if


         else if (nint(mask(i,j,k)) == P_INT_INSIDE_MASK-3) then
            !--------------------------------------------
            ! Verticies
            !--------------------------------------------
            xp = 0
            xm = 0
            yp = 0
            ym = 0
            zp = 0
            zm = 0
            if (nint(mask(i+1,j,k)) > 0) xp = 1
            if (nint(mask(i-1,j,k)) > 0) xm = 1
            if (nint(mask(i,j+1,k)) > 0) yp = 1
            if (nint(mask(i,j-1,k)) > 0) ym = 1
            if (nint(mask(i,j,k+1)) > 0) zp = 1
            if (nint(mask(i,j,k-1)) > 0) zm = 1

            !
            ! Vertex points
            !
            if (xp==1 .and. yp==1) then
               if (zp==1) then
                  imask(i+1,j+1,k+1) = P_STENCIL_X_VERTEX_P  &
                                     + P_STENCIL_Y_VERTEX_P  &
                                     + P_STENCIL_Z_VERTEX_P
                  if (modify_for_KO) then
                     ! vertex corrections
                     IMASK_MODX(imask(i+2,j+1,k+1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i+1,j+2,k+1), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i+1,j+1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k+1), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODX(imask(i+2,j,k+1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i,j+2,k+1), P_STENCIL_Y_CENTER_LEFT)

                     IMASK_MODXZ(imask(i,j+1,k), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODX(imask(i+2,j+1,k), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODZ(imask(i,j+1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     IMASK_MODYZ(imask(i+1,j,k), P_STENCIL_Y_CENTER_RIGHT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODY(imask(i+1,j+2,k), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i+1,j,k+2), P_STENCIL_Z_CENTER_LEFT)
                     imask(i+1,j,k+1) = P_STENCIL_Y_CENTER_RIGHT &
                                   + P_STENCIL_X_EDGE_P + P_STENCIL_Z_EDGE_P
                     imask(i,j+1,k+1) = P_STENCIL_X_CENTER_RIGHT &
                                   + P_STENCIL_Y_EDGE_P + P_STENCIL_Z_EDGE_P
                     imask(i+1,j+1,k) = P_STENCIL_Z_CENTER_RIGHT &
                                   + P_STENCIL_X_EDGE_P + P_STENCIL_Y_EDGE_P

                  else 
                     IMASK_MODXZ(imask(i+1,j,k+1), P_STENCIL_X_EDGE_P, P_STENCIL_Z_EDGE_P)
                     IMASK_MODYZ(imask(i,j+1,k+1), P_STENCIL_Y_EDGE_P, P_STENCIL_Z_EDGE_P)
                     IMASK_MODXY(imask(i+1,j+1,k), P_STENCIL_X_EDGE_P, P_STENCIL_Y_EDGE_P)
                  end if
               end if
               if (zm==1) then
                  imask(i+1,j+1,k-1) = P_STENCIL_X_VERTEX_P  &
                                     + P_STENCIL_Y_VERTEX_P  &
                                     + P_STENCIL_Z_VERTEX_M
                  if (modify_for_KO) then
                     ! Vertex corrections
                     IMASK_MODX(imask(i+2,j+1,k-1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i+1,j+2,k-1), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i+1,j+1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k-1), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODX(imask(i+2,j,k-1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i,j+2,k-1), P_STENCIL_Y_CENTER_LEFT)

                     IMASK_MODXZ(imask(i,j+1,k), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Z_CENTER_LEFT)
                     IMASK_MODX(imask(i+2,j+1,k), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODZ(imask(i,j+1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     IMASK_MODYZ(imask(i+1,j,k), P_STENCIL_Y_CENTER_RIGHT, P_STENCIL_Z_CENTER_LEFT)
                     IMASK_MODY(imask(i+1,j+2,k), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i+1,j,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     imask(i+1,j,k-1) = P_STENCIL_Y_CENTER_RIGHT &
                                   + P_STENCIL_X_EDGE_P + P_STENCIL_Z_EDGE_M
                     imask(i,j+1,k-1) = P_STENCIL_X_CENTER_RIGHT &
                                   + P_STENCIL_Y_EDGE_P + P_STENCIL_Z_EDGE_M
                     imask(i+1,j+1,k) = P_STENCIL_Z_CENTER_LEFT&
                                   + P_STENCIL_X_EDGE_P + P_STENCIL_Y_EDGE_P
                  else
                     IMASK_MODXZ(imask(i+1,j,k-1), P_STENCIL_X_EDGE_P, P_STENCIL_Z_EDGE_M)
                     IMASK_MODYZ(imask(i,j+1,k-1), P_STENCIL_Y_EDGE_P, P_STENCIL_Z_EDGE_M)
                     IMASK_MODXY(imask(i+1,j+1,k), P_STENCIL_X_EDGE_P, P_STENCIL_Y_EDGE_P)
                  end if
               end if
            end if

            if (xp==1 .and. ym==1) then
               if (zp==1) then
                  imask(i+1,j-1,k+1) = P_STENCIL_X_VERTEX_P  &
                                     + P_STENCIL_Y_VERTEX_M  &
                                     + P_STENCIL_Z_VERTEX_P
                  if (modify_for_KO) then
                     ! Vertex corrections
                     IMASK_MODX(imask(i+2,j-1,k+1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i+1,j-2,k+1), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i+1,j-1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k+1), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODX(imask(i+2,j,k+1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i,j-2,k+1), P_STENCIL_Y_CENTER_RIGHT)

                     IMASK_MODXZ(imask(i,j-1,k), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODX(imask(i+2,j-1,k), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODZ(imask(i,j-1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     IMASK_MODYZ(imask(i+1,j,k), P_STENCIL_Y_CENTER_LEFT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODY(imask(i+1,j-2,k), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i+1,j,k+2), P_STENCIL_Z_CENTER_LEFT)
                     imask(i+1,j,k+1) = P_STENCIL_Y_CENTER_LEFT &
                                 + P_STENCIL_X_EDGE_P + P_STENCIL_Z_EDGE_P
                     imask(i,j-1,k+1) = P_STENCIL_X_CENTER_RIGHT &
                                 + P_STENCIL_Y_EDGE_M + P_STENCIL_Z_EDGE_P
                     imask(i+1,j-1,k) = P_STENCIL_Z_CENTER_RIGHT &
                                   + P_STENCIL_X_EDGE_P + P_STENCIL_Y_EDGE_M
                  else
                     IMASK_MODXZ(imask(i+1,j,k+1), P_STENCIL_X_EDGE_P, P_STENCIL_Z_EDGE_P)
                     IMASK_MODYZ(imask(i,j-1,k+1), P_STENCIL_Y_EDGE_M, P_STENCIL_Z_EDGE_P)
                     IMASK_MODXY(imask(i+1,j-1,k), P_STENCIL_X_EDGE_P, P_STENCIL_Y_EDGE_M)
                  end if
               end if
               if (zm==1) then
                  imask(i+1,j-1,k-1) = P_STENCIL_X_VERTEX_P  &
                                     + P_STENCIL_Y_VERTEX_M  &
                                     + P_STENCIL_Z_VERTEX_M
                  if (modify_for_KO) then
                     ! Vertex corrections
                     IMASK_MODX(imask(i+2,j-1,k-1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i+1,j-2,k-1), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i+1,j-1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k-1), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODX(imask(i+2,j,k-1), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODY(imask(i,j-2,k-1), P_STENCIL_Y_CENTER_RIGHT)

                     IMASK_MODXZ(imask(i,j-1,k), P_STENCIL_X_CENTER_RIGHT, P_STENCIL_Z_CENTER_LEFT) 
                     IMASK_MODX(imask(i+2,j-1,k), P_STENCIL_X_CENTER_LEFT)
                     IMASK_MODZ(imask(i,j-1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     IMASK_MODYZ(imask(i+1,j,k), P_STENCIL_Y_CENTER_LEFT, P_STENCIL_Z_CENTER_LEFT) 
                     IMASK_MODY(imask(i+1,j-2,k), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i+1,j,k-2), P_STENCIL_Z_CENTER_RIGHT)
                     imask(i+1,j,k-1) = P_STENCIL_Y_CENTER_LEFT &
                                 + P_STENCIL_X_EDGE_P + P_STENCIL_Z_EDGE_M
                     imask(i,j-1,k-1) = P_STENCIL_X_CENTER_RIGHT &
                                 + P_STENCIL_Y_EDGE_M + P_STENCIL_Z_EDGE_M
                     imask(i+1,j-1,k) = P_STENCIL_Z_CENTER_LEFT &
                                   + P_STENCIL_X_EDGE_P + P_STENCIL_Y_EDGE_M
                  else
                     IMASK_MODXZ(imask(i+1,j,k-1), P_STENCIL_X_EDGE_P, P_STENCIL_Z_EDGE_M)
                     IMASK_MODYZ(imask(i,j-1,k-1), P_STENCIL_Y_EDGE_M, P_STENCIL_Z_EDGE_M)
                     IMASK_MODXY(imask(i+1,j-1,k), P_STENCIL_X_EDGE_P, P_STENCIL_Y_EDGE_M)
                  end if
               end if
            end if

            if (xm==1 .and. yp==1) then
               if (zp==1) then
                  imask(i-1,j+1,k+1) = P_STENCIL_X_VERTEX_M  &
                                     + P_STENCIL_Y_VERTEX_P  &
                                     + P_STENCIL_Z_VERTEX_P
                  if (modify_for_KO) then
                     ! Vertex corrections
                     IMASK_MODX(imask(i-2,j+1,k+1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i-1,j+2,k+1), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i-1,j+1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k+1), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODX(imask(i-2,j,k+1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i,j+2,k+1), P_STENCIL_Y_CENTER_LEFT)

                     IMASK_MODXZ(imask(i,j+1,k), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODX(imask(i-2,j+1,k), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODZ(imask(i,j+1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     IMASK_MODYZ(imask(i-1,j,k), P_STENCIL_Y_CENTER_RIGHT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODY(imask(i-1,j+2,k), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i-1,j,k+2), P_STENCIL_Z_CENTER_LEFT)

                     imask(i-1,j,k+1) = P_STENCIL_Y_CENTER_RIGHT &
                                 + P_STENCIL_X_EDGE_M + P_STENCIL_Z_EDGE_P
                     imask(i,j+1,k+1) = P_STENCIL_X_CENTER_LEFT &
                                 + P_STENCIL_Y_EDGE_P + P_STENCIL_Z_EDGE_P
                     imask(i-1,j+1,k) = P_STENCIL_Z_CENTER_RIGHT &
                                + P_STENCIL_X_EDGE_M + P_STENCIL_Y_EDGE_P
                  else
                     IMASK_MODXZ(imask(i-1,j,k+1), P_STENCIL_X_EDGE_M, P_STENCIL_Z_EDGE_P)
                     IMASK_MODYZ(imask(i,j+1,k+1), P_STENCIL_Y_EDGE_P, P_STENCIL_Z_EDGE_P)
                     IMASK_MODXY(imask(i-1,j+1,k), P_STENCIL_X_EDGE_M, P_STENCIL_Y_EDGE_P)
                  end if

               end if
               if (zm==1) then
                  imask(i-1,j+1,k-1) = P_STENCIL_X_VERTEX_M  &
                                     + P_STENCIL_Y_VERTEX_P  &
                                     + P_STENCIL_Z_VERTEX_M
                  if (modify_for_KO) then
                     ! Vertex corrections
                     IMASK_MODX(imask(i-2,j+1,k-1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i-1,j+2,k-1), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i-1,j+1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k-1), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODX(imask(i-2,j,k-1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i,j+2,k-1), P_STENCIL_Y_CENTER_LEFT)

                     IMASK_MODXZ(imask(i,j+1,k), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Z_CENTER_LEFT)
                     IMASK_MODX(imask(i-2,j+1,k), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODZ(imask(i,j+1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     IMASK_MODYZ(imask(i-1,j,k), P_STENCIL_Y_CENTER_RIGHT, P_STENCIL_Z_CENTER_LEFT)
                     IMASK_MODY(imask(i-1,j+2,k), P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODZ(imask(i-1,j,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     imask(i-1,j,k-1) = P_STENCIL_Y_CENTER_RIGHT &
                                 + P_STENCIL_X_EDGE_M + P_STENCIL_Z_EDGE_M
                     imask(i,j+1,k-1) = P_STENCIL_X_CENTER_LEFT &
                                 + P_STENCIL_Y_EDGE_P + P_STENCIL_Z_EDGE_M
                     imask(i-1,j+1,k) = P_STENCIL_Z_CENTER_LEFT &
                                + P_STENCIL_X_EDGE_M + P_STENCIL_Y_EDGE_P
                  else
                     IMASK_MODXZ(imask(i-1,j,k-1), P_STENCIL_X_EDGE_M, P_STENCIL_Z_EDGE_M)
                     IMASK_MODYZ(imask(i,j+1,k-1), P_STENCIL_Y_EDGE_P, P_STENCIL_Z_EDGE_M)
                     IMASK_MODXY(imask(i-1,j+1,k), P_STENCIL_X_EDGE_M, P_STENCIL_Y_EDGE_P)
                  end if
               end if
            end if

            if (xm==1 .and. ym==1) then
               if (zp==1) then
                  imask(i-1,j,k+1) = P_STENCIL_Y_CENTER_LEFT &
                                 + P_STENCIL_X_EDGE_M + P_STENCIL_Z_EDGE_P
                  imask(i,j-1,k+1) = P_STENCIL_X_CENTER_LEFT &
                                 + P_STENCIL_Y_EDGE_M + P_STENCIL_Z_EDGE_P
                  imask(i-1,j-1,k) = P_STENCIL_Z_CENTER_RIGHT &
                                   + P_STENCIL_X_EDGE_M + P_STENCIL_Y_EDGE_M
                  imask(i-1,j-1,k+1) = P_STENCIL_X_VERTEX_M  &
                                     + P_STENCIL_Y_VERTEX_M  &
                                     + P_STENCIL_Z_VERTEX_P
                  if (modify_for_KO) then
                     ! Vertex corrections
                     IMASK_MODX(imask(i-2,j-1,k+1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i-1,j-2,k+1), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i-1,j-1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k+1), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODX(imask(i-2,j,k+1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i,j-2,k+1), P_STENCIL_Y_CENTER_RIGHT)

                     IMASK_MODXZ(imask(i,j-1,k), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODX(imask(i-2,j-1,k), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODZ(imask(i,j-1,k+2), P_STENCIL_Z_CENTER_LEFT)

                     IMASK_MODYZ(imask(i-1,j,k), P_STENCIL_Y_CENTER_LEFT, P_STENCIL_Z_CENTER_RIGHT)
                     IMASK_MODY(imask(i-1,j-2,k), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i-1,j,k+2), P_STENCIL_Z_CENTER_LEFT)
                     imask(i-1,j,k+1) = P_STENCIL_Y_CENTER_LEFT &
                                 + P_STENCIL_X_EDGE_M + P_STENCIL_Z_EDGE_P
                     imask(i,j-1,k+1) = P_STENCIL_X_CENTER_LEFT &
                                 + P_STENCIL_Y_EDGE_M + P_STENCIL_Z_EDGE_P
                     imask(i-1,j-1,k) = P_STENCIL_Z_CENTER_RIGHT &
                                   + P_STENCIL_X_EDGE_M + P_STENCIL_Y_EDGE_M
                  else
                     IMASK_MODXZ(imask(i-1,j,k+1), P_STENCIL_X_EDGE_M, P_STENCIL_Z_EDGE_P)
                     IMASK_MODYZ(imask(i,j-1,k+1), P_STENCIL_Y_EDGE_M, P_STENCIL_Z_EDGE_P)
                     IMASK_MODXY(imask(i-1,j-1,k), P_STENCIL_X_EDGE_M, P_STENCIL_Y_EDGE_M)
                  end if

               end if
               if (zm==1) then
                  imask(i-1,j,k-1) = P_STENCIL_Y_CENTER_LEFT &
                                 + P_STENCIL_X_EDGE_M + P_STENCIL_Z_EDGE_M
                  imask(i,j-1,k-1) = P_STENCIL_X_CENTER_LEFT &
                                 + P_STENCIL_Y_EDGE_M + P_STENCIL_Z_EDGE_M
                  imask(i-1,j-1,k) = P_STENCIL_Z_CENTER_LEFT &
                                   + P_STENCIL_X_EDGE_M + P_STENCIL_Y_EDGE_M
                  imask(i-1,j-1,k-1) = P_STENCIL_X_VERTEX_M  &
                                     + P_STENCIL_Y_VERTEX_M  &
                                     + P_STENCIL_Z_VERTEX_M
                  if (modify_for_KO) then
                     ! Vertex corrections
                     IMASK_MODX(imask(i-2,j-1,k-1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i-1,j-2,k-1), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i-1,j-1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     ! Edge corrections
                     IMASK_MODXY(imask(i,j,k-1), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Y_CENTER_LEFT)
                     IMASK_MODX(imask(i-2,j,k-1), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODY(imask(i,j-2,k-1), P_STENCIL_Y_CENTER_RIGHT)

                     IMASK_MODXZ(imask(i,j-1,k), P_STENCIL_X_CENTER_LEFT, P_STENCIL_Z_CENTER_LEFT)
                     IMASK_MODX(imask(i-2,j-1,k), P_STENCIL_X_CENTER_RIGHT)
                     IMASK_MODZ(imask(i,j-1,k-2), P_STENCIL_Z_CENTER_RIGHT)

                     IMASK_MODYZ(imask(i-1,j,k), P_STENCIL_Y_CENTER_LEFT, P_STENCIL_Z_CENTER_LEFT)
                     IMASK_MODY(imask(i-1,j-2,k), P_STENCIL_Y_CENTER_RIGHT)
                     IMASK_MODZ(imask(i-1,j,k-2), P_STENCIL_Z_CENTER_RIGHT)
                     imask(i-1,j,k-1) = P_STENCIL_Y_CENTER_LEFT &
                                 + P_STENCIL_X_EDGE_M + P_STENCIL_Z_EDGE_M
                     imask(i,j-1,k-1) = P_STENCIL_X_CENTER_LEFT &
                                 + P_STENCIL_Y_EDGE_M + P_STENCIL_Z_EDGE_M
                     imask(i-1,j-1,k) = P_STENCIL_Z_CENTER_LEFT &
                                   + P_STENCIL_X_EDGE_M + P_STENCIL_Y_EDGE_M
                  else
                     IMASK_MODXZ(imask(i-1,j,k-1), P_STENCIL_X_EDGE_M, P_STENCIL_Z_EDGE_M)
                     IMASK_MODYZ(imask(i,j-1,k-1), P_STENCIL_Y_EDGE_M, P_STENCIL_Z_EDGE_M)
                     IMASK_MODXY(imask(i-1,j-1,k), P_STENCIL_X_EDGE_M, P_STENCIL_Y_EDGE_M)
                  end if
               end if
            end if
         end if
         end if
      end do
      end do
      end do


      return
      end subroutine fixup_imask

#if 0
!------------------------------------------------------------------------
!
!    This routine can be used to check that summation by parts
!    is numerically satisfied.  It checks this by calculating
!
!      residual = 2(u,du) - (u^2(x2) - u^2(x1)).  
!
!    along a line of y=const and z=const.  (j and k are input
!    parameters.)
!    If residual should be zero.  In the current condition the 
!    routine only partially works.  It works properly for no mask,
!    for lines that penetrate only faces of the mask.  It is
!    not clear if this function should be used again.
!
!------------------------------------------------------------------------
     CCTK_REAL function Check_line_integral_x(f, df, imask, j, k, &
               x, y, z, par)
     implicit  none
     
     CCTK_REAL, dimension(:,:,:)             :: f, df, x, y, z
     CCTK_INT                                :: j, k
     CCTK_INT, dimension(:,:,:)              :: imask
     CCTK_REAL, dimension(NPAR)                 :: par

     ! local vars
     CCTK_INT     :: i, nx, ny, nz, bpts, m
     CCTK_INT     :: center, face, medge, mvertex, undef, idt(3)
     CCTK_INT     :: mask_stencil, fpts, cpts, mvpts, mepts
     CCTK_REAL    :: sigma, dx, sum, sum2
     CCTK_REAL    :: int(MAXDOMAINS), endsum(MAXDOMAINS), endpt(MAXDOMAINS)
     CCTK_REAL    :: int2(MAXDOMAINS)

     dx = par(P_DX)

     nx = nint(par(P_NX))
     ny = nint(par(P_NY))
     nz = nint(par(P_NZ))

     if (j < 1 .or. j > ny) then
       write(0,*)'Check_line_integral_x:  j out of range ',j,ny
     end if
     if (k < 1 .or. k > nz) then
       write(0,*)'Check_line_integral_x:  k out of range ',k,nz
     end if

     mask_stencil = P_STENCIL_X_MASK + P_STENCIL_Y_MASK + P_STENCIL_Z_MASK

     do m = 1, MAXDOMAINS
       endsum(m) = 0.0
       endpt(m)  = 0.0
     end do
     bpts = 0

     cpts = 0
     fpts = 0
     mvpts = 0
     mepts = 0

     sum = 0.0
     do i = 1, nx
       idt(1) = 0
       idt(2) = 0
       idt(3) = 0
       if (imask(i,j,k) == 0) then
         ! interior point
         sigma = 1.0
         cpts = cpts + 1
       else if (imask(i,j,k) == mask_stencil) then
         ! point inside mask
         sigma = 0.0
       else
         center = 0
         face = 0
         medge = 0 
         mvertex = 0
         undef = 0
         idt(1) = imask(i,j,k) / 10000
         idt(2) = mod(imask(i,j,k)/100, 100)
         idt(3) = mod(imask(i,j,k), 100)
         do m = 1, 3
           if (idt(m) == P_STENCIL_CENTER_LEFT .or. &
               idt(m) == P_STENCIL_CENTER_RIGHT .or. &
               idt(m) == P_STENCIL_CENTER) then
             center = center + 1
           end if
           if (idt(m) == P_STENCIL_LEFT .or. &
               idt(m) == P_STENCIL_RIGHT) then
             face = face + 1
           end if
           if (idt(m) == P_STENCIL_EDGE_P .or. &
               idt(m) == P_STENCIL_EDGE_M) then
             medge = medge + 1
           end if
           if (idt(m) == P_STENCIL_VERTEX_P .or. &
               idt(m) == P_STENCIL_VERTEX_M) then
             mvertex = mvertex + 1
           end if
           if (idt(m) == P_STENCIL_UNDEF .or. &
               idt(m) == P_STENCIL_MASK) then
             undef = undef + 1
           end if
         end do
         if (center == 3) then
           ! interior point
           sigma = 1.0
           cpts = cpts + 1
         else if (undef > 0) then
           ! undefined point
           sigma = 0.0
         else if (mvertex > 0) then
           ! excision vertex (mask vertex)
           sigma = 7.0/8.0
           mvpts = mvpts + 1
         else if (medge > 0) then
           ! excision edge (mask edge)
           sigma = 0.75
           mepts = mepts + 1
         else if (face == 1) then
           ! boundary face
           sigma = 0.5
           fpts = fpts + 1
           !print *,'==> face point ',f(i,j,k),i,imask(i,j,k)
           !print *,'    ',x(i,j,k),y(i,j,k),z(i,j,k)
         else if (face == 2) then
           ! Outer edge (convex edge)
           sigma = 0.25
         else if (face == 3) then
           ! Outer vertex (convex vertex)
           sigma = 0.125
         else
           write(0,*)'Check_line_integral_x >>> unexpected point'
           stop
           sigma = 0.0
         end if
       end if

       sum = sum + sigma*f(i,j,k)*df(i,j,k)*dx

       if (idt(1) == P_STENCIL_LEFT .or. &
               idt(1) == P_STENCIL_RIGHT) then
         if (bpts == MAXDOMAINS) then
           write(0,*)'Check_line_integral_x>>> bpts is too large.'
           write(0,*)'increase MAXDOMAINS'
           stop
         end if
         bpts = bpts + 1
         !write(0,*)'adding domain ',idt(1), imask(i,j,k), bpts, i, j, k
         endpt(bpts) = f(i,j,k)
         endsum(bpts) = sum
         !if (bpts > 1 .and. mod(bpts,2) .ne. 0) sum = 0.0
         
       end if
     end do
     !print *,'XXXX sum ',sum

     if (mod(bpts,2) .ne. 0) then
       write(0,*) 'Check_line_integral_x>>> bpts is not even ',bpts
       write(0,*) '                         Error.'
       stop
     end if


!   sum = 0.5*f(1,j,k)*df(1,j,k)*dx
!   do i = 2, nx-1
!     sum = sum + f(i,j,k)*df(i,j,k)*dx
!   end do
!   sum = sum + 0.5*f(nx,j,k)*df(nx,j,k)*dx
!   int2(1) = 2.0*sum - (f(nx,j,k)**2 - f(1,j,k)**2)
!   print *,'==> int2=',int2(1),sum,f(1,j,k),f(nx,j,k)
!   print *,'==> fpts, cpts ',fpts,cpts

!    sum = 0.5d0*f(1,j,k)*(f(2,j,k)-f(1,j,k))
!    do i = 2, nx-1
!      sum = sum + f(i,j,k)*(f(i+1,j,k)-f(i-1,j,k))/2.d0
!    end do
!    sum = sum + 0.5d0*f(nx,j,k)*(f(nx,j,k) - f(nx-1,j,k))
!    int2(1) = 2.d0*sum - (f(nx,j,k)**2 - f(1,j,k)**2)

!    sum2 = 0.5*f(1,j,k)*df(1,j,k)*dx
!    do i = 2, 15
!      sum2 = sum2 + f(i,j,k)*df(i,j,k)*dx
!    end do
!    sum2 = sum2 + 0.5*f(16,j,k)*df(16,j,k)*dx
!    sum2 = sum2 + 0.5*f(26,j,k)*df(26,j,k)*dx
!    do i = 27, 40
!      sum2 = sum2 + f(i,j,k)*df(i,j,k)*dx
!    end do
!    sum2 = sum2 + 0.5*f(nx,j,k)*df(nx,j,k)*dx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Exact solution for very specific case
!
!    sum2 = 0.5*f(1,j,k)*(f(2,j,k)-f(1,j,k))
!    nx = 16
!    do i = 2, nx-1
!      sum2 = sum2 + f(i,j,k)*(f(i+1,j,k)-f(i-1,j,k))/2.d0
!    end do
!    sum2 = sum2 + 0.5d0*f(nx,j,k)*(f(nx,j,k) - f(nx-1,j,k))
!    int2(1) = 2.0*sum2 - (f(nx,j,k)**2 - f(1,j,k)**2)
!    print *,'==> 1st sum2 ',sum2, ' int2 ',int2(1)

!    nx = 41
!    sum2 = 0.5*f(26,j,k)*(f(27,j,k)-f(26,j,k))
!    do i = 27, nx-1
!      sum2 = sum2 + f(i,j,k)*(f(i+1,j,k)-f(i-1,j,k))/2.d0
!    end do
!    sum2 = sum2 + 0.5d0*f(nx,j,k)*(f(nx,j,k) - f(nx-1,j,k))
!    int2(2) = 2.0*sum2 - (f(nx,j,k)**2 - f(26,j,k)**2)
!    print *,'==> 2nd sum2 ',sum2, ' int2 ',int2(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Exact solution for very specific case

     sum2 = 0.5*f(1,j,k)*(f(2,j,k)-f(1,j,k))
     do i = 2, 15
       sum2 = sum2 + f(i,j,k)*(f(i+1,j,k)-f(i-1,j,k))/2.d0
     end do
     i = 16
     sum2 = sum2 + 0.75*f(i,j,k)*(-2.0*f(i-1,j,k)+f(i,j,k)+f(i+1,j,k))/3.d0
     do i = 17,25
       sum2 = sum2 + 0.5*f(i,j,k)*(f(i+1,j,k)-f(i-1,j,k))/2.d0
     end do
     i = 26
     sum2 = sum2 + 0.75*f(i,j,k)*(2.0*f(i+1,j,k)-f(i,j,k)-f(i-1,j,k))/3.d0
     do i = 27,40
       sum2 = sum2 + f(i,j,k)*(f(i+1,j,k)-f(i-1,j,k))/2.d0
     end do
     nx = 41
     sum2 = sum2 + 0.5d0*f(nx,j,k)*(f(nx,j,k) - f(nx-1,j,k))
     int2(1) = 2.0*sum2 - (f(nx,j,k)**2 - f(1,j,k)**2) - (f(16,j,k)**2 - f(26,j,k)**2)*0.5
     print *,'==> 1st sum2 ',sum2, ' int2 ',int2(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !write(*,*)'Integral:'
     !write(*,*)'...j,k,bpts = ',j,k, bpts
!    do i = 1, bpts-1, 2
!      int(i) = 2.0*endsum(i+1) - (endpt(i+1)**2 - endpt(i)**2)
!      if (int(i) .gt. 1.0) then
!        write(*,*)'Integral:   j,k,bpts = ',j,k,bpts
!        write(*,*)'...sum  ',endsum(i+1)
!        write(*,*)'...endpt2, endpt1, endpt2-endpt1 ',&
!                            endpt(i+1), endpt(i), endpt(i+1)-endpt(i)
!        write(*,*)'...int ',int(i)
!        !stop
!      end if
!      write(*,*)'...int ',int(i), int2(1),endsum(i+1)
!    end do

     int(1) = 2.0*sum
     !print *,'==> sum = ',sum
     do i = 1, bpts-1, 2
       int(1) = int(1) - (endpt(i+1)**2 - endpt(i)**2)
       !print *,'==> endpt2, endpt1 ',endpt(i+1), endpt(i)
     end do

     write(0,100)' ...fpts, cpts, mvpts, mepts, bpts, int ',&
                  fpts,cpts,mvpts,mepts,bpts,int(1)
 100 format(a,5i5,g14.3)

     Check_line_integral_x = 1

     return
     end function Check_line_integral_x
#endif

!------------------------------------------------------------------------
!
!     
!
!------------------------------------------------------------------------

      END MODULE MOD_MASK
