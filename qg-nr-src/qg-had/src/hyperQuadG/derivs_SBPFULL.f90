#include "cctk.h"

      !----------------------------------------------------------------
      !
      !  deriv_x  is the generic routine that outsources the actual 
      !           calculation of the x derivative of an array to 
      !           another routine.  
      !
      !----------------------------------------------------------------
      subroutine deriv_x(Dxu, u, dx, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT   nx, ny, nz
        CCTK_REAL  Dxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL  dx

        ! local vars
        CCTK_INT  d_type

        if ( d_type .eq. 42 ) then
          call deriv42_x(Dxu, u, dx, nx, ny, nz, d_type)

        else if ( d_type .eq. 44 ) then 
          call deriv44_x(Dxu, u, dx, nx, ny, nz, d_type)

        else if ( d_type .eq. 642 )  then
          call deriv642_x(Dxu, u, dx, nx, ny, nz, d_type)

        else if ( d_type .eq. 666 )  then
          call deriv666_x(Dxu, u, dx, nx, ny, nz, d_type)

        else if ( d_type .eq. 8642 )  then
          call deriv8642_x(Dxu, u, dx, nx, ny, nz, d_type)

        else if ( d_type .eq. 8888 )  then
          call deriv8888_x(Dxu, u, dx, nx, ny, nz, d_type)

        else
          write(0,*)'deriv_x: unknown d_type ',d_type

        end if

        return
      end subroutine

      !----------------------------------------------------------------
      !
      !  deriv_y  is the generic routine that outsources the actual 
      !           calculation of the y derivative of an array to 
      !           another routine.  
      !
      !----------------------------------------------------------------
      subroutine deriv_y(Dyu, u, dy, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy

        ! local vars
        CCTK_INT  d_type

             if ( d_type .eq. 42 ) then
          call deriv42_y(Dyu, u, dy, nx, ny, nz, d_type)

        else if ( d_type .eq. 44 )  then
          call deriv44_y(Dyu, u, dy, nx, ny, nz, d_type)

        else if ( d_type .eq. 642 )  then
          call deriv642_y(Dyu, u, dy, nx, ny, nz, d_type)

        else if ( d_type .eq. 666 )  then
          call deriv666_y(Dyu, u, dy, nx, ny, nz, d_type)

        else if ( d_type .eq. 8642 )  then
          call deriv8642_y(Dyu, u, dy, nx, ny, nz, d_type)

        else if ( d_type .eq. 8888 )  then
          call deriv8888_y(Dyu, u, dy, nx, ny, nz, d_type)

        else
          write(0,*)'deriv_y: unknown d_type ',d_type

        end if

        return
      end subroutine

      !----------------------------------------------------------------
      !
      !  deriv_z  is the generic routine that outsources the actual 
      !           calculation of the z derivative of an array to 
      !           another routine.  
      !
      !----------------------------------------------------------------
      subroutine deriv_z(Dzu, u, dz, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz

        ! local vars
        CCTK_INT  d_type

             if ( d_type .eq. 42 ) then
          call deriv42_z(Dzu, u, dz, nx, ny, nz, d_type)

        else if ( d_type .eq. 44 )  then
          call deriv44_z(Dzu, u, dz, nx, ny, nz, d_type)

        else if ( d_type .eq. 642 ) then
          call deriv642_z(Dzu, u, dz, nx, ny, nz, d_type)

        else if ( d_type .eq. 666 ) then
          call deriv666_z(Dzu, u, dz, nx, ny, nz, d_type)

        else if ( d_type .eq. 8642 ) then
          call deriv8642_z(Dzu, u, dz, nx, ny, nz, d_type)

        else if ( d_type .eq. 8888 ) then
          call deriv8888_z(Dzu, u, dz, nx, ny, nz, d_type)

        else
          write(0,*)'deriv_z: unknown d_type ',d_type

        end if

        return
      end subroutine

      !----------------------------------------------------------------
      !
      ! adv_deriv_x  is the generic routine that outsources the actual 
      !              calculation of the advective x derivative of an 
      !              array to another routine.  
      !
      !----------------------------------------------------------------
      subroutine adv_deriv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)
        implicit   none
        CCTK_INT   nx, ny, nz
        CCTK_REAL  Dxu(nx,ny,nz), u(nx,ny,nz), betax(nx,ny,nz)
        CCTK_REAL  dx

        ! local vars
        CCTK_INT  d_type

        if ( d_type .eq. 42 )  then
          !call deriv42adv_x(Dxu, u, dx, nx, ny, nz, par, betax)
          call deriv42adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)

        else if ( d_type .eq. 44 )  then
          call deriv44adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)

        else if ( d_type .eq. 642 )  then
          call deriv642adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)

        else if ( d_type .eq. 666 )  then
          call deriv666adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)

        else if ( d_type .eq. 8642 )  then
          call deriv8642adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)

        else if ( d_type .eq. 8888 )  then
          call deriv8888adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)

        else
          write(0,*)'adv_deriv_x: unknown d_type ',d_type

        end if

        return
      end subroutine

      !----------------------------------------------------------------
      !
      ! adv_deriv_y  is the generic routine that outsources the actual 
      !              calculation of the advective y derivative of an 
      !              array to another routine.  
      !
      !----------------------------------------------------------------
      subroutine adv_deriv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz), betay(nx,ny,nz)
        CCTK_REAL dy

        ! local vars
        CCTK_INT  d_type

        if ( d_type .eq. 42 )  then
          call deriv42adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)

        else if ( d_type .eq. 44 )  then
          call deriv44adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)

        else if ( d_type .eq. 642 )  then
          call deriv642adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)

        else if ( d_type .eq. 666 )  then
          call deriv666adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)

        else if ( d_type .eq. 8642 )  then
          call deriv8642adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)

        else if ( d_type .eq. 8888 )  then
          call deriv8888adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)

        else
          write(0,*)'adv_deriv_y: unknown d_type ',d_type

        end if

        return
      end subroutine

      !----------------------------------------------------------------
      !
      ! adv_deriv_z  is the generic routine that outsources the actual 
      !              calculation of the advective z derivative of an 
      !              array to another routine.  
      !
      !----------------------------------------------------------------
      subroutine adv_deriv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz), betaz(nx,ny,nz)
        CCTK_REAL dz

        ! local vars
        CCTK_INT  d_type

        if ( d_type .eq. 42 ) then
          call deriv42adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)

        else if ( d_type .eq. 44 ) then
          call deriv44adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)

        else if ( d_type .eq. 642 ) then
          call deriv642adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)

        else if ( d_type .eq. 666 ) then
          call deriv666adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)

        else if ( d_type .eq. 8642 )  then
          call deriv8642adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)

        else if ( d_type .eq. 8888 )  then
          call deriv8888adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)

        else
          write(0,*)'adv_deriv_z: unknown d_type ',d_type

        end if

        return
      end subroutine


      !----------------------------------------------------------------
      !
      !  deriv_xx  is the generic routine that outsources the actual 
      !            calculation of the second x derivative of an array 
      !            to another routine.  
      !
      !----------------------------------------------------------------
      subroutine deriv_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DxDxu(nx,ny,nz), u(nx,ny,nz) 
        CCTK_REAL dx

        ! local vars
        CCTK_INT  d_type, dd_type

	d_type=dd_type 

             if ( dd_type .eq. 42 ) then
          call deriv42_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
          !call deriv42_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          !call deriv42_x(DxDxu, Dxu, dx, nx, ny, nz, d_type)

        else if ( dd_type .eq. 44 ) then
          call deriv44_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
          !call deriv44_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          !call deriv44_x(DxDxu, Dxu, dx, nx, ny, nz, d_type)

        else if ( dd_type .eq. 642 ) then
          call deriv642_xx(DxDxu, u, dx, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 666 ) then
          call deriv666_xx(DxDxu, u, dx, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 8642 ) then
          call deriv8642_xx(DxDxu, u, dx, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 8888 ) then
          call deriv8888_xx(DxDxu, u, dx, nx, ny, nz, dd_type)

        end if

        return
      end subroutine

      !----------------------------------------------------------------
      !
      !  deriv_yy  is the generic routine that outsources the actual 
      !            calculation of the second y derivative of an array 
      !            to another routine.  
      !
      !----------------------------------------------------------------
      subroutine deriv_yy(DyDyu, u, dy, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy

        ! local vars
        CCTK_INT  dd_type

             if ( dd_type .eq. 42 ) then
          call deriv42_yy(DyDyu, u, dy, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 44 ) then
          call deriv44_yy(DyDyu, u, dy, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 642 ) then
          call deriv642_yy(DyDyu, u, dy, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 666 ) then
          call deriv666_yy(DyDyu, u, dy, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 8642 ) then
          call deriv8642_yy(DyDyu, u, dy, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 8888 ) then
          call deriv8888_yy(DyDyu, u, dy, nx, ny, nz, dd_type)

        end if

        return
      end subroutine

      !----------------------------------------------------------------
      !
      !  deriv_zz  is the generic routine that outsources the actual 
      !            calculation of the second z derivative of an array 
      !            to another routine.  
      !
      !----------------------------------------------------------------
      subroutine deriv_zz(DzDzu, u, dz, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DzDzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz

        ! local vars
        CCTK_INT  dd_type

             if ( dd_type .eq. 42 ) then
          call deriv42_zz(DzDzu, u, dz, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 44 ) then
          call deriv44_zz(DzDzu, u, dz, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 642 ) then
          call deriv642_zz(DzDzu, u, dz, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 666 ) then
          call deriv666_zz(DzDzu, u, dz, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 8642 ) then
          call deriv8642_zz(DzDzu, u, dz, nx, ny, nz, dd_type)

        else if ( dd_type .eq. 8888 ) then
          call deriv8888_zz(DzDzu, u, dz, nx, ny, nz, dd_type)

        end if

        return
      end subroutine


      !----------------------------------------------------------------
      !
      !  deriv_xy  is the generic routine that outsources the actual 
      !            calculation of the mixed partial derivatives with 
      !            respect to x and y of an array to another routine.  
      !
      !       Dxu is a work array
      !
      !----------------------------------------------------------------
      subroutine deriv_xy(DxDyu, u, Dxu, dx, dy, nx, ny, nz, dd_type)
!!      subroutine deriv_xy(DxDyu, u, dx, dy, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL dx, dy
        CCTK_REAL DxDyu(nx,ny,nz), u(nx,ny,nz), Dxu(nx,ny,nz)

        ! local vars
        CCTK_INT  d_type, dd_type

	d_type = dd_type 

             if ( dd_type .eq.  42 ) then
          call deriv42_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv42_y(DxDyu, Dxu, dy, nx, ny, nz, d_type)

        else if ( dd_type .eq. 44 ) then
          call deriv44_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv44_y(DxDyu, Dxu, dy, nx, ny, nz, d_type)

        else if ( dd_type .eq. 642 ) then
          call deriv642_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv642_y(DxDyu, Dxu, dy, nx, ny, nz, d_type)

        else if ( dd_type .eq. 666 ) then
          call deriv666_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv666_y(DxDyu, Dxu, dy, nx, ny, nz, d_type)

        else if ( dd_type .eq. 8642 ) then
          call deriv8642_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv8642_y(DxDyu, Dxu, dy, nx, ny, nz, d_type)

        else if ( dd_type .eq. 8888 ) then
          call deriv8888_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv8888_y(DxDyu, Dxu, dy, nx, ny, nz, d_type)

        end if

        return
      end subroutine


      !----------------------------------------------------------------
      !
      !  deriv_xz  is the generic routine that outsources the actual 
      !            calculation of the mixed partial derivatives with 
      !            respect to x and z of an array to another routine.  
      !
      !       Dxu is a work array
      !
      !----------------------------------------------------------------
      subroutine deriv_xz(DxDzu, u, Dxu, dx, dz, nx, ny, nz, dd_type)
!!      subroutine deriv_xz(DxDzu, u, dx, dz, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DxDzu(nx,ny,nz), u(nx,ny,nz), Dxu(nx,ny,nz)
        CCTK_REAL dx, dz

        ! local vars
        CCTK_INT  d_type, dd_type

	d_type = dd_type 

             if ( dd_type .eq. 42 ) then
          call deriv42_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv42_z(DxDzu, Dxu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 44 ) then
          call deriv44_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv44_z(DxDzu, Dxu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 642 ) then
          call deriv642_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv642_z(DxDzu, Dxu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 666 ) then
          call deriv666_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv666_z(DxDzu, Dxu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 8642 ) then
          call deriv8642_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv8642_z(DxDzu, Dxu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 8888 ) then
          call deriv8888_x(  Dxu,   u, dx, nx, ny, nz, d_type)
          call deriv8888_z(DxDzu, Dxu, dz, nx, ny, nz, d_type)

        end if

        return
      end subroutine


      !----------------------------------------------------------------
      !
      !  deriv_yz  is the generic routine that outsources the actual 
      !            calculation of the mixed partial derivatives with 
      !            respect to y and z of an array to another routine.  
      !
      !
      !       Dyu is a work array
      !
      !----------------------------------------------------------------
      subroutine deriv_yz(DyDzu, u, Dyu, dy, dz, nx, ny, nz, dd_type)
!!      subroutine deriv_yz(DyDzu, u, dy, dz, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDzu(nx,ny,nz), u(nx,ny,nz), Dyu(nx,ny,nz)
        CCTK_REAL dy, dz

        ! local vars
        CCTK_INT  d_type, dd_type

	d_type = dd_type 

             if ( dd_type .eq. 42 ) then
          call deriv42_y(  Dyu,   u, dy, nx, ny, nz, d_type)
          call deriv42_z(DyDzu, Dyu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 44 ) then
          call deriv44_y(  Dyu,   u, dy, nx, ny, nz, d_type)
          call deriv44_z(DyDzu, Dyu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 642 ) then
          call deriv642_y(  Dyu,   u, dy, nx, ny, nz, d_type)
          call deriv642_z(DyDzu, Dyu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 666 ) then
          call deriv666_y(  Dyu,   u, dy, nx, ny, nz, d_type)
          call deriv666_z(DyDzu, Dyu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 8642 ) then
          call deriv8642_y(  Dyu,   u, dy, nx, ny, nz, d_type)
          call deriv8642_z(DyDzu, Dyu, dz, nx, ny, nz, d_type)

        else if ( dd_type .eq. 8888 ) then
          call deriv8888_y(  Dyu,   u, dy, nx, ny, nz, d_type)
          call deriv8888_z(DyDzu, Dyu, dz, nx, ny, nz, d_type)

        end if

        return
      end subroutine







!======================================================================! 
!                                                                      ! 
!    FOURTH ORDER DERIVATIVE OPERATORS                                 ! 
!                                                                      ! 
!      these include:                                                  !  
!         deriv42_x, deriv42_y, deriv42_z,                             ! 
!         deriv44_x, deriv44_y, deriv44_z,                             ! 
!         deriv42adv_x, deriv42adv_y, deriv42adv_z,                    ! 
!         deriv44adv_x, deriv44adv_y, deriv44adv_z,                    ! 
!         deriv42_xx, deriv42_yy, deriv42_zz,                          ! 
!         deriv44_xx, deriv44_yy, deriv44_zz                           ! 
!                                                                      ! 
!======================================================================! 

      !----------------------------------------------------------------
      !
      !   deriv42_x  SBP, a la mattson, see:
      !  http://www.it.uu.se/research/publications/reports/2003-012/
      !
      !----------------------------------------------------------------
      subroutine deriv42_x(Dxu, u, dx, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idx, idx_by_2, idx_by_12 

        idx = 1.0d0/dx
        idx_by_2 = 0.5d0 * idx 
        idx_by_12 = idx / 12.d0  

        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  48.d0 * u(i  ,j,k)     &      
     &                   +  59.d0 * u(i+1,j,k)     &     
     &                   -   8.d0 * u(i+2,j,k)     &   
     &                   -   3.d0 * u(i+3,j,k)     &   
     &                 ) * idx/34.d0
          i =2
          Dxu(i,j,k) = ( - u(i-1,j,k)     &      
     &                   + u(i+1,j,k)     &      
     &                 ) * idx_by_2

          i=3
          Dxu(i,j,k) = (  8.d0*u(i-2,j,k)   &
     &                  -59.d0*u(i-1,j,k)   &
     &                  +59.d0*u(i+1,j,k)   &
     &                  - 8.d0*u(i+2,j,k)   &   
     &                 ) * idx/86.d0

	  i=4
          Dxu(i,j,k) = ( 3.d0  *u(i-3,j,k)   &
     &                  -59.d0* u(i-1,j,k)   &
     &                  +64.d0* u(i+1,j,k)   &
     &                  - 8.d0* u(i+2,j,k)   &
     &                 ) * idx/98.d0

          do i = 5, nx-4
            Dxu(i,j,k) = (          u(i-2,j,k)     &      
     &                     - 8.d0 * u(i-1,j,k)     &      
     &                     + 8.d0 * u(i+1,j,k)     &      
     &                     -        u(i+2,j,k)     &      
     &                   ) * idx_by_12  
          end do

	  i=nx-3
          Dxu(i,j,k) = -( 3.d0  *u(i+3,j,k)   &
     &                   -59.d0* u(i+1,j,k)   &
     &                   +64.d0* u(i-1,j,k)   &
     &                   - 8.d0* u(i-2,j,k)   &
     &                  ) * idx/98.d0


          i=nx-2
          Dxu(i,j,k) =-(  8.d0*u(i+2,j,k)   &
     &                  -59.d0*u(i+1,j,k)   &
     &                  +59.d0*u(i-1,j,k)   &
     &                  - 8.d0*u(i-2,j,k)   &   
     &                 ) * idx/86.d0


          i=nx-1 
          Dxu(i,j,k) = ( - u(i-1,j,k)     &      
     &                   + u(i+1,j,k)     &      
     &                 ) * idx_by_2 
          i=nx 
          Dxu(i,j,k) = - ( -  48.d0 * u(i  ,j,k)     &      
     &                     +  59.d0 * u(i-1,j,k)     &     
     &                     -   8.d0 * u(i-2,j,k)     &   
     &                     -   3.d0 * u(i-3,j,k)     &   
     &                   ) * idx/34.d0

        end do
        end do

        return
      end subroutine 
      
 


      !----------------------------------------------------------------
      !
      !   deriv42_y  SBP 
      !
      !----------------------------------------------------------------
      subroutine deriv42_y(Dyu, u, dy, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idy, idy_by_2, idy_by_12 

        idy = 1.0d0/dy
        idy_by_2 = 0.5d0 * idy 
        idy_by_12 = idy / 12.d0  

        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  48.d0 * u(i  ,j,k)     &      
     &                   +  59.d0 * u(i,j+1,k)     &     
     &                   -   8.d0 * u(i,j+2,k)     &   
     &                   -   3.d0 * u(i,j+3,k)     &   
     &                 ) * idy/34.d0
          j=2
          Dyu(i,j,k) = ( - u(i,j-1,k)     &      
     &                   + u(i,j+1,k)     &      
     &                 ) * idy_by_2

          j=3
          Dyu(i,j,k) = (  8.d0*u(i,j-2,k)   &
     &                  -59.d0*u(i,j-1,k)   &
     &                  +59.d0*u(i,j+1,k)   &
     &                  - 8.d0*u(i,j+2,k)   &   
     &                 ) * idy/86.d0

	  j=4
          Dyu(i,j,k) = ( 3.d0  *u(i,j-3,k)   &
     &                  -59.d0* u(i,j-1,k)   &
     &                  +64.d0* u(i,j+1,k)   &
     &                  - 8.d0* u(i,j+2,k)   &
     &                 ) * idy/98.d0

          do j = 5, ny-4
            Dyu(i,j,k) = (          u(i,j-2,k)     &      
     &                     - 8.d0 * u(i,j-1,k)     &      
     &                     + 8.d0 * u(i,j+1,k)     &      
     &                     -        u(i,j+2,k)     &      
     &                   ) * idy_by_12  
          end do

	  j=ny-3
          Dyu(i,j,k) = -( 3.  *u(i,j+3,k)   &
     &                   -59.* u(i,j+1,k)   &
     &                   +64.* u(i,j-1,k)   &
     &                   - 8.* u(i,j-2,k)   &
     &                  ) * idy/98.

          j=ny-2
          Dyu(i,j,k) =-(  8.d0*u(i,j+2,k)   &
     &                  -59.d0*u(i,j+1,k)   &
     &                  +59.d0*u(i,j-1,k)   &
     &                  - 8.d0*u(i,j-2,k)   &   
     &                 ) * idy/86.d0

          j=ny-1 
          Dyu(i,j,k) = ( - u(i,j-1,k)     &      
     &                   + u(i,j+1,k)     &      
     &                 ) * idy_by_2 
          j=ny 
          Dyu(i,j,k) = - ( -  48.d0 * u(i  ,j,k)     &      
     &                     +  59.d0 * u(i,j-1,k)     &     
     &                     -   8.d0 * u(i,j-2,k)     &   
     &                     -   3.d0 * u(i,j-3,k)     &   
     &                   ) * idy/34.d0


        end do
        end do

        return
      end subroutine 




      !----------------------------------------------------------------
      !
      !   deriv42_z  SBP
      !
      !----------------------------------------------------------------
      subroutine deriv42_z(Dzu, u, dz, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idz, idz_by_2, idz_by_12 

        idz = 1.0d0/dz
        idz_by_2 = 0.5d0 * idz 
        idz_by_12 = idz / 12.d0  

        do j = 1, ny
        do i = 1, nx

	 k=1 
          Dzu(i,j,k) = ( -  48.d0 * u(i  ,j,k)     &      
     &                   +  59.d0 * u(i,j,k+1)     &     
     &                   -   8.d0 * u(i,j,k+2)     &   
     &                   -   3.d0 * u(i,j,k+3)     &   
     &                 ) * idz/34.d0
          k=2
          Dzu(i,j,k) = ( - u(i,j,k-1)     &      
     &                   + u(i,j,k+1)     &      
     &                 ) * idz_by_2

          k=3
          Dzu(i,j,k) = (  8.d0*u(i,j,k-2)   &
     &                  -59.d0*u(i,j,k-1)   &
     &                  +59.d0*u(i,j,k+1)   &
     &                  - 8.d0*u(i,j,k+2)   &   
     &                 ) * idz/86.d0

	  k=4
          Dzu(i,j,k) = ( 3.d0  *u(i,j,k-3)   &
     &                  -59.d0* u(i,j,k-1)   &
     &                  +64.d0* u(i,j,k+1)   &
     &                  - 8.d0* u(i,j,k+2)   &
     &                 ) * idz/98.d0


          do k = 5, nz-4
            Dzu(i,j,k) = (          u(i,j,k-2)     &      
     &                     - 8.d0 * u(i,j,k-1)     &      
     &                     + 8.d0 * u(i,j,k+1)     &      
     &                     -        u(i,j,k+2)     &      
     &                   ) * idz_by_12  
          end do

	  k=nz-3
          Dzu(i,j,k) = -( 3.d0  *u(i,j,k+3)   &
     &                   -59.d0* u(i,j,k+1)   &
     &                   +64.d0* u(i,j,k-1)   &
     &                   - 8.d0* u(i,j,k-2)   &
     &                  ) * idz/98.d0

          k=nz-2
          Dzu(i,j,k) =-(  8.d0*u(i,j,k+2)   &
     &                  -59.d0*u(i,j,k+1)   &
     &                  +59.d0*u(i,j,k-1)   &
     &                  - 8.d0*u(i,j,k-2)   &   
     &                 ) * idz/86.d0

          k=nz-1 
          Dzu(i,j,k) = ( - u(i,j,k-1)     &      
     &                   + u(i,j,k+1)     &      
     &                 ) * idz_by_2 

          k=nz 
          Dzu(i,j,k) = - ( -  48.d0 * u(i  ,j,k)     &      
     &                     +  59.d0 * u(i,j,k-1)     &     
     &                     -   8.d0 * u(i,j,k-2)     &   
     &                     -   3.d0 * u(i,j,k-3)     &   
     &                   ) * idz/34.d0


        end do
        end do

        return
      end subroutine 




      !----------------------------------------------------------------
      !
      !   deriv44_x  calculates the first derivative with respect to 
      !              x using a centered fourth order finite difference 
      !              operator at all points except the boundaries (i=1
      !              and i=nx) and the penultimate boundary points 
      !              (i=2 and i=nx-1).  At these latter points, we use
      !              the appropriate up- or downwind fourth order 
      !              operators.  Note that each of these shifted 
      !              operators use five point stencils.  
      !
      !----------------------------------------------------------------      
      subroutine deriv44_x(Dxu, u, dx, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12 

        idx = 1.0d0/dx        
        idx_by_2 = 0.5d0 * idx
        idx_by_12 = idx / 12.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( - 25.d0 * u(i  ,j,k)     &      
     &                   + 48.d0 * u(i+1,j,k)     &               
     &                   - 36.d0 * u(i+2,j,k)     &               
     &                   + 16.d0 * u(i+3,j,k)     &               
     &                   -  3.d0 * u(i+4,j,k)     &      
     &                 ) * idx_by_12  
          i=2
          Dxu(i,j,k) = (  -  3.d0 * u(i-1,j,k)     &      
     &                    - 10.d0 * u(i  ,j,k)     &     
     &                    + 18.d0 * u(i+1,j,k)     &     
     &                    -  6.d0 * u(i+2,j,k)     &      
     &                    +         u(i+3,j,k)     &      
     &                 ) * idx_by_12
          do i = 3, nx-2
            Dxu(i,j,k) = (          u(i-2,j,k)     &      
     &                     - 8.d0 * u(i-1,j,k)     &     
     &                     + 8.d0 * u(i+1,j,k)     &     
     &                     -        u(i+2,j,k)     &     
     &                   ) * idx_by_12
          end do
          i=nx-1 
          Dxu(i,j,k) = (  -         u(i-3,j,k)     &      
     &                    +  6.d0 * u(i-2,j,k)     &     
     &                    - 18.d0 * u(i-1,j,k)     &     
     &                    + 10.d0 * u(i  ,j,k)     &      
     &                    +  3.d0 * u(i+1,j,k)     &      
     &                 ) * idx_by_12
          i=nx 
          Dxu(i,j,k) = (    3.d0 * u(i-4,j,k)     &      
     &                   - 16.d0 * u(i-3,j,k)     &               
     &                   + 36.d0 * u(i-2,j,k)     &               
     &                   - 48.d0 * u(i-1,j,k)     &               
     &                   + 25.d0 * u(i  ,j,k)     &      
     &                 ) * idx_by_12  
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      !   deriv44_y  calculates the first derivative with respect to 
      !              y using a centered fourth order finite difference 
      !              operator at all points except the boundaries (j=1
      !              and j=ny) and the penultimate boundary points 
      !              (j=2 and j=ny-1).  At these latter points, we use
      !              the appropriate up- or downwind fourth order 
      !              operators.  Note that each of these shifted 
      !              operators use five point stencils.  
      !
      !----------------------------------------------------------------      
      subroutine deriv44_y(Dyu, u, dy, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12 

        idy = 1.0d0/dy        
        idy_by_2 = 0.5d0 * idy
        idy_by_12 = idy / 12.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( - 25.d0 * u(i,j  ,k)     &      
     &                   + 48.d0 * u(i,j+1,k)     &               
     &                   - 36.d0 * u(i,j+2,k)     &               
     &                   + 16.d0 * u(i,j+3,k)     &               
     &                   -  3.d0 * u(i,j+4,k)     &      
     &                 ) * idy_by_12  
          j=2
          Dyu(i,j,k) = (  -  3.d0 * u(i,j-1,k)     &      
     &                    - 10.d0 * u(i,j  ,k)     &     
     &                    + 18.d0 * u(i,j+1,k)     &     
     &                    -  6.d0 * u(i,j+2,k)     &      
     &                    +         u(i,j+3,k)     &      
     &                 ) * idy_by_12
          do j = 3, ny-2
            Dyu(i,j,k) = (          u(i,j-2,k)     &      
     &                     - 8.d0 * u(i,j-1,k)     &     
     &                     + 8.d0 * u(i,j+1,k)     &     
     &                     -        u(i,j+2,k)     &     
     &                   ) * idy_by_12
          end do
          j=ny-1 
          Dyu(i,j,k) = (  -         u(i,j-3,k)     &      
     &                    +  6.d0 * u(i,j-2,k)     &     
     &                    - 18.d0 * u(i,j-1,k)     &     
     &                    + 10.d0 * u(i,j  ,k)     &      
     &                    +  3.d0 * u(i,j+1,k)     &      
     &                 ) * idy_by_12
          j=ny 
          Dyu(i,j,k) = (     3.d0 * u(i,j-4,k)     &      
     &                    - 16.d0 * u(i,j-3,k)     &               
     &                    + 36.d0 * u(i,j-2,k)     &               
     &                    - 48.d0 * u(i,j-1,k)     &               
     &                    + 25.d0 * u(i,j  ,k)     &      
     &                 ) * idy_by_12  
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      !   deriv44_z  calculates the first derivative with respect to 
      !              z using a centered fourth order finite difference 
      !              operator at all points except the boundaries (k=1
      !              and k=nz) and the penultimate boundary points 
      !              (k=2 and k=nz-1).  At these latter points, we use
      !              the appropriate up- or downwind fourth order 
      !              operators.  Note that each of these shifted 
      !              operators use five point stencils.  
      !
      !----------------------------------------------------------------      
      subroutine deriv44_z(Dzu, u, dz, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12 

        idz = 1.0d0/dz        
        idz_by_2 = 0.5d0 * idz
        idz_by_12 = idz / 12.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( - 25.d0 * u(i,j,k  )     &      
     &                   + 48.d0 * u(i,j,k+1)     &               
     &                   - 36.d0 * u(i,j,k+2)     &               
     &                   + 16.d0 * u(i,j,k+3)     &               
     &                   -  3.d0 * u(i,j,k+4)     &      
     &                 ) * idz_by_12  
          k=2
          Dzu(i,j,k) = (  -  3.d0 * u(i,j,k-1)     &      
     &                    - 10.d0 * u(i,j,k  )     &     
     &                    + 18.d0 * u(i,j,k+1)     &     
     &                    -  6.d0 * u(i,j,k+2)     &      
     &                    +         u(i,j,k+3)     &      
     &                 ) * idz_by_12
          do k = 3, nz-2
            Dzu(i,j,k) = (          u(i,j,k-2)     &      
     &                     - 8.d0 * u(i,j,k-1)     &     
     &                     + 8.d0 * u(i,j,k+1)     &     
     &                     -        u(i,j,k+2)     &     
     &                   ) * idz_by_12
          end do
          k=nz-1 
          Dzu(i,j,k) = (  -         u(i,j,k-3)     &      
     &                    +  6.d0 * u(i,j,k-2)     &     
     &                    - 18.d0 * u(i,j,k-1)     &     
     &                    + 10.d0 * u(i,j,k  )     &      
     &                    +  3.d0 * u(i,j,k+1)     &      
     &                 ) * idz_by_12
          k=nz 
          Dzu(i,j,k) = (     3.d0 * u(i,j,k-4)     &      
     &                    - 16.d0 * u(i,j,k-3)     &               
     &                    + 36.d0 * u(i,j,k-2)     &               
     &                    - 48.d0 * u(i,j,k-1)     &               
     &                    + 25.d0 * u(i,j,k  )     &      
     &                  ) * idz_by_12  
        end do
        end do
       
        return
      end subroutine 
     

      !----------------------------------------------------------------
      !
      !  deriv42adv_x  calculates the first derivative with respect to 
      !                x using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [4,nx-3].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^x \partial_x in the equations.  
      !                Further, it will depend on the sign of beta^x.  
      !                Namely, if beta^x is positive, we will use a 
      !                fourth order upwind operator while if beta^x 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At i=3, if beta^x is positive, we continue to  
      !                use the fourth order upwind operator but if  
      !                beta^x is negative, we use a second order downwind  
      !                operator.   
      !                  
      !                We reverse this at i=nx-2 such that if beta^x is  
      !                positive, we use a second order upwind operator 
      !                while if beta^x is negative we continue to use  
      !                the fourth order downwind operator.   
      !                 
      !                At i=2, if beta^x is positive, we use the   
      !                second order upwind operator and if beta^x is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At i=nx-1, if beta^x is positive, we use the   
      !                second order centered operator and if beta^x is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At i=1, we use the second order upwind operator 
      !                for both positive and negative values of beta^x.
      !                
      !                At i=nx, we use the second order downwind operator 
      !                for both positive and negative values of beta^x.
      !                 
      !----------------------------------------------------------------      
      subroutine deriv42adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz), betax(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12 

        idx = 1.0d0/dx        
        idx_by_2 = 0.5d0 * idx
        idx_by_12 = idx / 12.d0


        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  3.d0 * u(i  ,j,k)     &      
     &                   +  4.d0 * u(i+1,j,k)     &     
     &                   -         u(i+2,j,k)     &     
     &                 ) * idx_by_2

          i=2
          if (betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -  3.d0 * u(i  ,j,k)     &      
     &                     +  4.d0 * u(i+1,j,k)     &     
     &                     -         u(i+2,j,k)     &     
     &                   ) * idx_by_2 
          else 
            Dxu(i,j,k) = ( -         u(i-1,j,k)     &      
     &                     +         u(i+1,j,k)     &      
     &                   ) * idx_by_2 
          endif 

          i=3
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                     - 10.d0 * u(i  ,j,k)     &     
     &                     + 18.d0 * u(i+1,j,k)     &     
     &                     -  6.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_12 
          else 
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                     -  4.d0 * u(i-1,j,k)     &     
     &                     +  3.d0 * u(i  ,j,k)     &     
     &                   ) * idx_by_2 
          endif

          do i = 4, nx-3
            if ( betax(i,j,k) .ge. 0.d0 ) then
              Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                       - 10.d0 * u(i  ,j,k)     &     
     &                       + 18.d0 * u(i+1,j,k)     &     
     &                       -  6.d0 * u(i+2,j,k)     &     
     &                       +         u(i+3,j,k)     &     
     &                     ) * idx_by_12 
            else 
              Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                       +  6.d0 * u(i-2,j,k)     &     
     &                       - 18.d0 * u(i-1,j,k)     &     
     &                       + 10.d0 * u(i  ,j,k)     &     
     &                       +  3.d0 * u(i+1,j,k)     &     
     &                     ) * idx_by_12 
            endif 
          end do
 
          i=nx-2 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (  - 3.d0 * u(i  ,j,k)     &      
     &                      + 4.d0 * u(i+1,j,k)     &      
     &                      -        u(i+2,j,k)     &     
     &                   ) * idx_by_2 
          else 
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  6.d0 * u(i-2,j,k)     &     
     &                     - 18.d0 * u(i-1,j,k)     &     
     &                     + 10.d0 * u(i  ,j,k)     &     
     &                     +  3.d0 * u(i+1,j,k)     &     
     &                   ) * idx_by_12 
          endif

          i=nx-1 
          if ( betax(i,j,k) .ge. 0.d0 ) then 
            Dxu(i,j,k) = (  -         u(i-1,j,k)     &      
     &                      +         u(i+1,j,k)     &     
     &                    ) * idx_by_2
          else 
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                      - 4.d0 * u(i-1,j,k)     &      
     &                      + 3.d0 * u(i  ,j,k)     &      
     &                   ) * idx_by_2 
          endif

          i=nx 
          Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                    - 4.d0 * u(i-1,j,k)     &      
     &                    + 3.d0 * u(i  ,j,k)     &      
     &                 ) * idx_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      


      !----------------------------------------------------------------
      !
      !  deriv42adv_y  calculates the first derivative with respect to 
      !                y using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [4,ny-3].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^y \partial_y in the equations.  
      !                Further, it will depend on the sign of beta^y.  
      !                Namely, if beta^y is positive, we will use a 
      !                fourth order upwind operator while if beta^y 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At j=3, if beta^y is positive, we continue to  
      !                use the fourth order upwind operator but if  
      !                beta^y is negative, we use a second order downwind  
      !                operator.   
      !                  
      !                We reverse this at j=ny-2 such that if beta^y is  
      !                positive, we use a second order upwind operator 
      !                while if beta^y is negative we continue to use  
      !                the fourth order downwind operator.   
      !                 
      !                At j=2, if beta^y is positive, we use the   
      !                second order upwind operator and if beta^y is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At j=ny-1, if beta^y is positive, we use the   
      !                second order centered operator and if beta^y is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At j=1, we use the second order upwind operator 
      !                for both positive and negative values of beta^y.
      !                
      !                At j=ny, we use the second order downwind operator 
      !                for both positive and negative values of beta^y.
      !                 
      !----------------------------------------------------------------      
      subroutine deriv42adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz), betay(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12 

        idy = 1.0d0/dy        
        idy_by_2 = 0.5d0 * idy
        idy_by_12 = idy / 12.d0


        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  3.d0 * u(i,j  ,k)     &      
     &                   +  4.d0 * u(i,j+1,k)     &     
     &                   -         u(i,j+2,k)     &     
     &                 ) * idy_by_2

          j=2
          if (betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -  3.d0 * u(i,j  ,k)     &      
     &                     +  4.d0 * u(i,j+1,k)     &     
     &                     -         u(i,j+2,k)     &     
     &                   ) * idy_by_2 
          else 
            Dyu(i,j,k) = ( -         u(i,j-1,k)     &      
     &                     +         u(i,j+1,k)     &      
     &                   ) * idy_by_2 
          endif 

          j=3
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                     - 10.d0 * u(i,j  ,k)     &     
     &                     + 18.d0 * u(i,j+1,k)     &     
     &                     -  6.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_12 
          else 
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                     -  4.d0 * u(i,j-1,k)     &     
     &                     +  3.d0 * u(i,j  ,k)     &     
     &                   ) * idy_by_2 
          endif

          do j = 4, ny-3
            if ( betay(i,j,k) .ge. 0.d0 ) then
              Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                       - 10.d0 * u(i,j  ,k)     &     
     &                       + 18.d0 * u(i,j+1,k)     &     
     &                       -  6.d0 * u(i,j+2,k)     &     
     &                       +         u(i,j+3,k)     &     
     &                     ) * idy_by_12 
            else 
              Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                       +  6.d0 * u(i,j-2,k)     &     
     &                       - 18.d0 * u(i,j-1,k)     &     
     &                       + 10.d0 * u(i,j  ,k)     &     
     &                       +  3.d0 * u(i,j+1,k)     &     
     &                     ) * idy_by_12 
            endif 
          end do
 
          j=ny-2 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (  - 3.d0 * u(i,j  ,k)     &      
     &                      + 4.d0 * u(i,j+1,k)     &      
     &                      -        u(i,j+2,k)     &     
     &                   ) * idy_by_2 
          else 
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  6.d0 * u(i,j-2,k)     &     
     &                     - 18.d0 * u(i,j-1,k)     &     
     &                     + 10.d0 * u(i,j  ,k)     &     
     &                     +  3.d0 * u(i,j+1,k)     &     
     &                   ) * idy_by_12 
          endif

          j=ny-1 
          if ( betay(i,j,k) .ge. 0.d0 ) then 
            Dyu(i,j,k) = (  -         u(i,j-1,k)     &      
     &                      +         u(i,j+1,k)     &     
     &                   ) * idy_by_2
          else 
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                      - 4.d0 * u(i,j-1,k)     &      
     &                      + 3.d0 * u(i,j  ,k)     &      
     &                   ) * idy_by_2 
          endif

          j=ny 
          Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                    - 4.d0 * u(i,j-1,k)     &      
     &                    + 3.d0 * u(i,j  ,k)     &      
     &                 ) * idy_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      


      !----------------------------------------------------------------
      !
      !  deriv42adv_z  calculates the first derivative with respect to 
      !                z using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [4,nz-3].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^z \partial_z in the equations.  
      !                Further, it will depend on the sign of beta^z.  
      !                Namely, if beta^z is positive, we will use a 
      !                fourth order upwind operator while if beta^z 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At k=3, if beta^z is positive, we continue to  
      !                use the fourth order upwind operator but if  
      !                beta^z is negative, we use a second order downwind  
      !                operator.   
      !                  
      !                We reverse this at k=nz-2 such that if beta^z is  
      !                positive, we use a second order upwind operator 
      !                while if beta^z is negative we continue to use  
      !                the fourth order downwind operator.   
      !                 
      !                At k=2, if beta^z is positive, we use the   
      !                second order upwind operator and if beta^z is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At k=nz-1, if beta^z is positive, we use the   
      !                second order centered operator and if beta^z is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At k=1, we use the second order upwind operator 
      !                for both positive and negative values of beta^z.
      !                
      !                At k=nz, we use the second order downwind operator 
      !                for both positive and negative values of beta^z.
      !                 
      !----------------------------------------------------------------      
      subroutine deriv42adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz), betaz(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12 

        idz = 1.0d0/dz        
        idz_by_2 = 0.5d0 * idz
        idz_by_12 = idz / 12.d0


        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( -  3.d0 * u(i,j,k  )     &      
     &                   +  4.d0 * u(i,j,k+1)     &     
     &                   -         u(i,j,k+2)     &     
     &                 ) * idz_by_2

          k=2
          if (betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -  3.d0 * u(i,j,k  )     &      
     &                     +  4.d0 * u(i,j,k+1)     &     
     &                     -         u(i,j,k+2)     &     
     &                   ) * idz_by_2 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-1)     &      
     &                     +         u(i,j,k+1)     &      
     &                   ) * idz_by_2 
          endif 

          k=3
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                     - 10.d0 * u(i,j,k  )     &     
     &                     + 18.d0 * u(i,j,k+1)     &     
     &                     -  6.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_12 
          else 
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                     -  4.d0 * u(i,j,k-1)     &     
     &                     +  3.d0 * u(i,j,k  )     &     
     &                   ) * idz_by_2 
          endif

          do k = 4, nz-3
            if ( betaz(i,j,k) .ge. 0.d0 ) then
              Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                       - 10.d0 * u(i,j,k  )     &     
     &                       + 18.d0 * u(i,j,k+1)     &     
     &                       -  6.d0 * u(i,j,k+2)     &     
     &                       +         u(i,j,k+3)     &     
     &                     ) * idz_by_12 
            else 
              Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                       +  6.d0 * u(i,j,k-2)     &     
     &                       - 18.d0 * u(i,j,k-1)     &     
     &                       + 10.d0 * u(i,j,k  )     &     
     &                       +  3.d0 * u(i,j,k+1)     &     
     &                     ) * idz_by_12 
            endif 
          end do
 
          k=nz-2 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (  - 3.d0 * u(i,j,k  )     &      
     &                      + 4.d0 * u(i,j,k+1)     &      
     &                      -        u(i,j,k+2)     &     
     &                   ) * idz_by_2 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  6.d0 * u(i,j,k-2)     &     
     &                     - 18.d0 * u(i,j,k-1)     &     
     &                     + 10.d0 * u(i,j,k  )     &     
     &                     +  3.d0 * u(i,j,k+1)     &     
     &                   ) * idz_by_12 
          endif

          k=nz-1 
          if ( betaz(i,j,k) .ge. 0.d0 ) then 
            Dzu(i,j,k) = (  -         u(i,j,k-1)     &      
     &                      +         u(i,j,k+1)     &     
     &                   ) * idz_by_2
          else 
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                      - 4.d0 * u(i,j,k-1)     &      
     &                      + 3.d0 * u(i,j,k  )     &      
     &                   ) * idz_by_2 
          endif
          k=nz 
          
          Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                    - 4.d0 * u(i,j,k-1)     &      
     &                    + 3.d0 * u(i,j,k  )     &      
     &                  ) * idz_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      

      !----------------------------------------------------------------
      !
      !  deriv44adv_x  calculates the first derivative with respect to 
      !                x using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [4,nx-3].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^x \partial_x in the equations.  
      !                Further, it will depend on the sign of beta^x.  
      !                Namely, if beta^x is positive, we will use a 
      !                fourth order upwind operator while if beta^x 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At i=3, if beta^x is positive, we continue to  
      !                use the fourth order upwind operator but if  
      !                beta^x is negative, we use the fourth order 
      !                centered operator.   
      !                  
      !                We reverse this at i=nx-2 such that if beta^x is  
      !                positive, we use the fourth order centered operator 
      !                while if beta^x is negative we continue to use  
      !                the fourth order downwind operator.   
      !                 
      !                At i=2, we want a fourth order operator so it 
      !                does not matter what the sign of beta^x is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of  
      !                using the downwind operator for the case of 
      !                negative beta^x, but I do not see any way around 
      !                this if we want to have a fourth order operator.  
      !                 
      !                Likewise, at i=nx-1, the requirement of a fourth 
      !                order operator renders the sign of beta^x 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^x is positive.  
      !                 
      !                At i=1 and i=nx, we use the totally shifted fourth
      !                order operators irrespective of beta^xs sign.
      !
      !----------------------------------------------------------------      
      subroutine deriv44adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz), betax(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12 

        idx = 1.0d0/dx        
        idx_by_2 = 0.5d0 * idx
        idx_by_12 = idx / 12.d0


        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( - 25.d0 * u(i  ,j,k)     &      
     &                   + 48.d0 * u(i+1,j,k)     &               
     &                   - 36.d0 * u(i+2,j,k)     &               
     &                   + 16.d0 * u(i+3,j,k)     &               
     &                   -  3.d0 * u(i+4,j,k)     &      
     &                 ) * idx_by_12  

          i=2
          Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                   - 10.d0 * u(i  ,j,k)     &     
     &                   + 18.d0 * u(i+1,j,k)     &     
     &                   -  6.d0 * u(i+2,j,k)     &     
     &                   +         u(i+3,j,k)     &     
     &                 ) * idx_by_12 
        
          i=3
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                     - 10.d0 * u(i  ,j,k)     &     
     &                     + 18.d0 * u(i+1,j,k)     &     
     &                     -  6.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_12 
          else 
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                     -  8.d0 * u(i-1,j,k)     &     
     &                     +  8.d0 * u(i+1,j,k)     &     
     &                     -         u(i+2,j,k)     &     
     &                   ) * idx_by_12 
          endif
 
          do i = 4, nx-3
            if ( betax(i,j,k) .ge. 0.d0 ) then
              Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                       - 10.d0 * u(i  ,j,k)     &     
     &                       + 18.d0 * u(i+1,j,k)     &     
     &                       -  6.d0 * u(i+2,j,k)     &     
     &                       +         u(i+3,j,k)     &     
     &                     ) * idx_by_12 
            else 
              Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                       +  6.d0 * u(i-2,j,k)     &     
     &                       - 18.d0 * u(i-1,j,k)     &     
     &                       + 10.d0 * u(i  ,j,k)     &     
     &                       +  3.d0 * u(i+1,j,k)     &     
     &                     ) * idx_by_12 
            endif 
          end do

          i=nx-2 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                      - 8.d0 * u(i-1,j,k)     &      
     &                      + 8.d0 * u(i+1,j,k)     &     
     &                      -        u(i+2,j,k)     &     
     &                   ) * idx_by_12 
          else 
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  6.d0 * u(i-2,j,k)     &     
     &                     - 18.d0 * u(i-1,j,k)     &     
     &                     + 10.d0 * u(i  ,j,k)     &     
     &                     +  3.d0 * u(i+1,j,k)     &     
     &                   ) * idx_by_12 
          endif

          i=nx-1 
          Dxu(i,j,k) = (  -         u(i-3,j,k)     &      
     &                    +  6.d0 * u(i-2,j,k)     &     
     &                    - 18.d0 * u(i-1,j,k)     &     
     &                    + 10.d0 * u(i  ,j,k)     &      
     &                    +  3.d0 * u(i+1,j,k)     &      
     &                  ) * idx_by_12
          i=nx 
          Dxu(i,j,k) = (     3.d0 * u(i-4,j,k)     &      
     &                    - 16.d0 * u(i-3,j,k)     &               
     &                    + 36.d0 * u(i-2,j,k)     &               
     &                    - 48.d0 * u(i-1,j,k)     &               
     &                    + 25.d0 * u(i  ,j,k)     &      
     &                   ) * idx_by_12  
        end do
        end do
       
        return
      end subroutine 
      

      !----------------------------------------------------------------
      !
      !  deriv44adv_y  calculates the first derivative with respect to 
      !                y using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [4,ny-3].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^y \partial_y in the equations.  
      !                Further, it will depend on the sign of beta^y.  
      !                Namely, if beta^y is positive, we will use a 
      !                fourth order upwind operator while if beta^y 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At j=3, if beta^y is positive, we continue to  
      !                use the fourth order upwind operator but if  
      !                beta^y is negative, we use the fourth order 
      !                centered operator.   
      !                  
      !                We reverse this at j=ny-2 such that if beta^y is  
      !                positive, we use the fourth order centered operator 
      !                while if beta^y is negative we continue to use  
      !                the fourth order downwind operator.   
      !                 
      !                At j=2, we want a fourth order operator so it 
      !                does not matter what the sign of beta^y is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of  
      !                using the downwind operator for the case of 
      !                negative beta^y, but I do not see any way around 
      !                this if we want to have a fourth order operator.  
      !                 
      !                Likewise, at j=ny-1, the requirement of a fourth 
      !                order operator renders the sign of beta^y 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^y is positive.  
      !                 
      !                At j=1 and j=ny, we use the totally shifted fourth
      !                order operators irrespective of beta^ys sign.
      !
      !----------------------------------------------------------------      
      subroutine deriv44adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz), betay(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12 

        idy = 1.0d0/dy        
        idy_by_2 = 0.5d0 * idy
        idy_by_12 = idy / 12.d0


        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( - 25.d0 * u(i,j  ,k)     &      
     &                   + 48.d0 * u(i,j+1,k)     &               
     &                   - 36.d0 * u(i,j+2,k)     &               
     &                   + 16.d0 * u(i,j+3,k)     &               
     &                   -  3.d0 * u(i,j+4,k)     &      
     &                 ) * idy_by_12  

          j=2
          Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                   - 10.d0 * u(i,j  ,k)     &     
     &                   + 18.d0 * u(i,j+1,k)     &     
     &                   -  6.d0 * u(i,j+2,k)     &     
     &                   +         u(i,j+3,k)     &     
     &                 ) * idy_by_12 
        
          j=3
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                     - 10.d0 * u(i,j  ,k)     &     
     &                     + 18.d0 * u(i,j+1,k)     &     
     &                     -  6.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_12 
          else 
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                     -  8.d0 * u(i,j-1,k)     &     
     &                     +  8.d0 * u(i,j+1,k)     &     
     &                     -         u(i,j+2,k)     &     
     &                   ) * idy_by_12 
          endif
 
          do j = 4, ny-3
            if ( betay(i,j,k) .ge. 0.d0 ) then
              Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                       - 10.d0 * u(i,j  ,k)     &     
     &                       + 18.d0 * u(i,j+1,k)     &     
     &                       -  6.d0 * u(i,j+2,k)     &     
     &                       +         u(i,j+3,k)     &     
     &                     ) * idy_by_12 
            else 
              Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                       +  6.d0 * u(i,j-2,k)     &     
     &                       - 18.d0 * u(i,j-1,k)     &     
     &                       + 10.d0 * u(i,j  ,k)     &     
     &                       +  3.d0 * u(i,j+1,k)     &     
     &                     ) * idy_by_12 
            endif 
          end do

          j=ny-2 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                      - 8.d0 * u(i,j-1,k)     &      
     &                      + 8.d0 * u(i,j+1,k)     &     
     &                      -        u(i,j+2,k)     &     
     &                   ) * idy_by_12 
          else 
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  6.d0 * u(i,j-2,k)     &     
     &                     - 18.d0 * u(i,j-1,k)     &     
     &                     + 10.d0 * u(i,j  ,k)     &     
     &                     +  3.d0 * u(i,j+1,k)     &     
     &                   ) * idy_by_12 
          endif

          j=ny-1 
          Dyu(i,j,k) = (  -         u(i,j-3,k)     &      
     &                    +  6.d0 * u(i,j-2,k)     &     
     &                    - 18.d0 * u(i,j-1,k)     &     
     &                    + 10.d0 * u(i,j  ,k)     &      
     &                    +  3.d0 * u(i,j+1,k)     &      
     &                 ) * idy_by_12
          j=ny 
          Dyu(i,j,k) = (     3.d0 * u(i,j-4,k)     &      
     &                    - 16.d0 * u(i,j-3,k)     &               
     &                    + 36.d0 * u(i,j-2,k)     &               
     &                    - 48.d0 * u(i,j-1,k)     &               
     &                    + 25.d0 * u(i,j  ,k)     &      
     &                  ) * idy_by_12  
        end do
        end do
       
        return
      end subroutine 
      

      !----------------------------------------------------------------
      !
      !  deriv44adv_z  calculates the first derivative with respect to 
      !                z using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [4,nz-3].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^z \partial_z in the equations.  
      !                Further, it will depend on the sign of beta^z.  
      !                Namely, if beta^z is positive, we will use a 
      !                fourth order upwind operator while if beta^z 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At k=3, if beta^z is positive, we continue to  
      !                use the fourth order upwind operator but if  
      !                beta^z is negative, we use the fourth order 
      !                centered operator.   
      !                  
      !                We reverse this at k=nz-2 such that if beta^z is  
      !                positive, we use the fourth order centered operator 
      !                while if beta^z is negative we continue to use  
      !                the fourth order downwind operator.   
      !                 
      !                At k=2, we want a fourth order operator so it 
      !                does not matter what the sign of beta^z is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of  
      !                using the downwind operator for the case of 
      !                negative beta^z, but I do not see any way around 
      !                this if we want to have a fourth order operator.  
      !                 
      !                Likewise, at k=nz-1, the requirement of a fourth 
      !                order operator renders the sign of beta^z 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^z is positive.  
      !                 
      !                At k=1 and k=nz, we use the totally shifted fourth
      !                order operators irrespective of beta^zs sign.
      !
      !----------------------------------------------------------------      
      subroutine deriv44adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz), betaz(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12 

        idz = 1.0d0/dz        
        idz_by_2 = 0.5d0 * idz
        idz_by_12 = idz / 12.d0


        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( - 25.d0 * u(i,j,k  )     &      
     &                   + 48.d0 * u(i,j,k+1)     &               
     &                   - 36.d0 * u(i,j,k+2)     &               
     &                   + 16.d0 * u(i,j,k+3)     &               
     &                   -  3.d0 * u(i,j,k+4)     &      
     &                 ) * idz_by_12  

          k=2
          Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                   - 10.d0 * u(i,j,k  )     &     
     &                   + 18.d0 * u(i,j,k+1)     &     
     &                   -  6.d0 * u(i,j,k+2)     &     
     &                   +         u(i,j,k+3)     &     
     &                 ) * idz_by_12 
        
          k=3
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                     - 10.d0 * u(i,j,k  )     &     
     &                     + 18.d0 * u(i,j,k+1)     &     
     &                     -  6.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_12 
          else 
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                     -  8.d0 * u(i,j,k-1)     &     
     &                     +  8.d0 * u(i,j,k+1)     &     
     &                     -         u(i,j,k+2)     &     
     &                   ) * idz_by_12 
          endif
 
          do k = 4, nz-3
            if ( betaz(i,j,k) .ge. 0.d0 ) then
              Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                       - 10.d0 * u(i,j,k  )     &     
     &                       + 18.d0 * u(i,j,k+1)     &     
     &                       -  6.d0 * u(i,j,k+2)     &     
     &                       +         u(i,j,k+3)     &     
     &                     ) * idz_by_12 
            else 
              Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                       +  6.d0 * u(i,j,k-2)     &     
     &                       - 18.d0 * u(i,j,k-1)     &     
     &                       + 10.d0 * u(i,j,k  )     &     
     &                       +  3.d0 * u(i,j,k+1)     &     
     &                     ) * idz_by_12 
            endif 
          end do

          k=nz-2 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                      - 8.d0 * u(i,j,k-1)     &      
     &                      + 8.d0 * u(i,j,k+1)     &     
     &                      -        u(i,j,k+2)     &     
     &                   ) * idz_by_12 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  6.d0 * u(i,j,k-2)     &     
     &                     - 18.d0 * u(i,j,k-1)     &     
     &                     + 10.d0 * u(i,j,k  )     &     
     &                     +  3.d0 * u(i,j,k+1)     &     
     &                   ) * idz_by_12 
          endif

          k=nz-1 
          Dzu(i,j,k) = (  -         u(i,j,k-3)     &      
     &                    +  6.d0 * u(i,j,k-2)     &     
     &                    - 18.d0 * u(i,j,k-1)     &     
     &                    + 10.d0 * u(i,j,k  )     &      
     &                    +  3.d0 * u(i,j,k+1)     &      
     &                 ) * idz_by_12
          k=nz 
          Dzu(i,j,k) = (     3.d0 * u(i,j,k-4)     &      
     &                    - 16.d0 * u(i,j,k-3)     &               
     &                    + 36.d0 * u(i,j,k-2)     &               
     &                    - 48.d0 * u(i,j,k-1)     &               
     &                    + 25.d0 * u(i,j,k  )     &      
     &                  ) * idz_by_12  
        end do
        end do
       
        return
      end subroutine 
     


      !----------------------------------------------------------------
      !
      !  deriv42_xx  calculates the second derivative with respect to  
      !              x using a centered fourth order finite difference 
      !              operator at all points except 4 points near the boundary
      !              For those we impose SBP condition following.  
      !              Ken Mattsson and Jan Nordström, Technical Report 2003-012
      !              http://www.it.uu.se/research/publications/reports/2003-012/
      !              careful, they have a type in page 30, operator D2 that we use
      !              here. 4th row of D2 matrix, last relevan factor is -4/49 
      !              (not 43 as they have it)
      !
      !----------------------------------------------------------------
      subroutine deriv42_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DxDxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idx_sqrd, idx_sqrd_by_12

        idx_sqrd = 1.0d0/(dx*dx)
        idx_sqrd_by_12 = idx_sqrd / 12.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          DxDxu(i,j,k) = (   2.d0 * u(i  ,j,k)     &     
     &                     - 5.d0 * u(i+1,j,k)     &     
     &                     + 4.d0 * u(i+2,j,k)     &     
     &                     -        u(i+3,j,k)     &     
     &                   ) * idx_sqrd
!          DxDxu(i,j,k) = (    45.d0 * u(i  ,j,k)     &
!     &                     - 154.d0 * u(i+1,j,k)     &
!     &                     + 214.d0 * u(i+2,j,k)     &
!     &                     - 156.d0 * u(i+3,j,k)     &
!     &                     +  61.d0 * u(i+4,j,k)     &
!     &                     -  10.d0 * u(i+5,j,k)     &
!     &                   ) * idx_sqrd_by_12
          i=2
          DxDxu(i,j,k) = (          u(i-1,j,k)     &     
     &                     - 2.d0 * u(i  ,j,k)     &     
     &                     +        u(i+1,j,k)     &     
     &                   ) * idx_sqrd
!          DxDxu(i,j,k) = (   10.d0 * u(i-1,j,k)     &     
!     &                     - 15.d0 * u(i  ,j,k)     &     
!     &                     -  4.d0 * u(i+1,j,k)     &     
!     &                     + 14.d0 * u(i+2,j,k)     &     
!     &                     -  6.d0 * u(i+3,j,k)     &     
!     &                     +         u(i+4,j,k)     &     
!     &                   ) * idx_sqrd_by_12 
	
	   i=3
	   DxDxu(i,j,k) = (-4.d0  *u(i-2,j,k)      &
     &                     +59.d0 *u(i-1,j,k)      &
     &                     -110.d0*u(i  ,j,k)      &
     &                     +59.d0 *u(i+1,j,k)      &
     &                     -4.d0  *u(i+2,j,k)      &
     &                    ) * idx_sqrd/43.d0

	   i=4
	   DxDxu(i,j,k) = (-        u(i-3,j,k)         &
     &                     +  59.d0*u(i-1,j,k)         &
     &                     - 118.d0*u(i  ,j,k)         &
     &                     + 64.d0 *u(i+1,j,k)         &
     &                     -  4.d0 *u(i+2,j,k)         &
     &                    ) * idx_sqrd/49.d0

          do i = 5, nx-4
            DxDxu(i,j,k) = ( -         u(i-2,j,k)     &      
     &                       + 16.d0 * u(i-1,j,k)     &     
     &                       - 30.d0 * u(i  ,j,k)     &      
     &                       + 16.d0 * u(i+1,j,k)     &      
     &                       -         u(i+2,j,k)     &      
     &                     ) * idx_sqrd_by_12  
          end do

	  i=nx-3
	   DxDxu(i,j,k) = (-        u(i+3,j,k)         &
     &                     +  59.d0*u(i+1,j,k)         &
     &                     - 118.d0*u(i  ,j,k)         &
     &                     + 64.d0 *u(i-1,j,k)         &
     &                     -  4.d0 *u(i-2,j,k)         &
     &                    ) * idx_sqrd/49.d0

          i=nx-2
	   DxDxu(i,j,k) = (-4.d0  *u(i+2,j,k)      &
     &                     +59.d0 *u(i+1,j,k)      &
     &                     -110.d0*u(i  ,j,k)      &
     &                     +59.d0 *u(i-1,j,k)      &
     &                     -4.d0  *u(i-2,j,k)      &
     &                    ) * idx_sqrd/43.d0


          i=nx-1 
          DxDxu(i,j,k) = (          u(i-1,j,k)     &     
     &                     - 2.d0 * u(i  ,j,k)     &      
     &                     +        u(i+1,j,k)     &      
     &                   ) * idx_sqrd 
!          DxDxu(i,j,k) = (           u(i-4,j,k)     &     
!     &                     -  6.d0 * u(i-3,j,k)     &      
!     &                     + 14.d0 * u(i-2,j,k)     &      
!     &                     -  4.d0 * u(i-1,j,k)     &      
!     &                     - 15.d0 * u(i  ,j,k)     &      
!     &                     + 10.d0 * u(i+1,j,k)     &      
!     &                   ) * idx_sqrd_by_12 
          i=nx 
          DxDxu(i,j,k) = ( -        u(i-3,j,k)     &     
     &                     + 4.d0 * u(i-2,j,k)     &      
     &                     - 5.d0 * u(i-1,j,k)     &      
     &                     + 2.d0 * u(i  ,j,k)     &      
     &                   ) * idx_sqrd

!          DxDxu(i,j,k) = ( -  10.d0 * u(i-5,j,k)     &
!     &                     +  61.d0 * u(i-4,j,k)     &
!     &                     - 156.d0 * u(i-3,j,k)     &
!     &                     + 214.d0 * u(i-2,j,k)     &
!     &                     - 154.d0 * u(i-1,j,k)     &
!     &                     +  45.d0 * u(i  ,j,k)     &
!     &                   ) * idx_sqrd_by_12
        end do
        end do

        return
      end subroutine 
      
 

      !----------------------------------------------------------------
      !
      !  deriv42_yy SBP
      !
      !----------------------------------------------------------------
      subroutine deriv42_yy(DyDyu, u, dy, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idy_sqrd, idy_sqrd_by_12

        idy_sqrd = 1.0d0/(dy*dy)
        idy_sqrd_by_12 = idy_sqrd / 12.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          DyDyu(i,j,k) = (   2.d0 * u(i,j  ,k)     &     
     &                     - 5.d0 * u(i,j+1,k)     &     
     &                     + 4.d0 * u(i,j+2,k)     &     
     &                     -        u(i,j+3,k)     &     
     &                   ) * idy_sqrd
          j=2
          DyDyu(i,j,k) = (          u(i,j-1,k)     &     
     &                     - 2.d0 * u(i,j  ,k)     &     
     &                     +        u(i,j+1,k)     &     
     &                   ) * idy_sqrd

	   j=3
	   DyDyu(i,j,k) = (-4.d0  *u(i,j-2,k)      &
     &                     +59.d0 *u(i,j-1,k)      &
     &                     -110.d0*u(i  ,j,k)      &
     &                     +59.d0 *u(i,j+1,k)      &
     &                     -4.d0  *u(i,j+2,k)      &
     &                    ) * idy_sqrd/43.d0

	   j=4
	   DyDyu(i,j,k) = (-       u(i,j-3,k)         &
     &                     +  59.d0*u(i,j-1,k)         &
     &                     - 118.d0*u(i  ,j,k)         &
     &                     + 64.d0 *u(i,j+1,k)         &
     &                     -  4.d0 *u(i,j+2,k)         &
     &                    ) * idy_sqrd/49.d0


          do j = 3, ny-2
            DyDyu(i,j,k) = ( -         u(i,j-2,k)     &     
     &                       + 16.d0 * u(i,j-1,k)     &     
     &                       - 30.d0 * u(i,j  ,k)     &     
     &                       + 16.d0 * u(i,j+1,k)     &     
     &                       -         u(i,j+2,k)     &     
     &                     ) * idy_sqrd_by_12
          end do

	  j=ny-3
	   DyDyu(i,j,k) = (-        u(i,j+3,k)         &
     &                     +  59.d0*u(i,j+1,k)         &
     &                     - 118.d0*u(i  ,j,k)         &
     &                     + 64.d0 *u(i,j-1,k)         &
     &                     -  4.d0 *u(i,j-2,k)         &
     &                    ) * idy_sqrd/49.d0

          j=ny-2
	   DyDyu(i,j,k) = (-4.d0  *u(i,j+2,k)      &
     &                     +59.d0 *u(i,j+1,k)      &
     &                     -110.d0*u(i  ,j,k)      &
     &                     +59.d0 *u(i,j-1,k)      &
     &                     -4.d0  *u(i,j-2,k)      &
     &                    ) * idy_sqrd/43.d0

          j=ny-1
          DyDyu(i,j,k) = (          u(i,j-1,k)     &     
     &                     - 2.d0 * u(i,j  ,k)     &     
     &                     +        u(i,j+1,k)     &     
     &                   ) * idy_sqrd
          j=ny
          DyDyu(i,j,k) = ( -        u(i,j-3,k)     &     
     &                     + 4.d0 * u(i,j-2,k)     &      
     &                     - 5.d0 * u(i,j-1,k)     &      
     &                     + 2.d0 * u(i,j  ,k)     &      
     &                   ) * idy_sqrd
        end do
        end do

        return
      end subroutine



      !----------------------------------------------------------------
      !
      !  deriv42_zz  SBP
      !
      !----------------------------------------------------------------
      subroutine deriv42_zz(DzDzu, u, dz, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DzDzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idz_sqrd, idz_sqrd_by_12

        idz_sqrd = 1.0d0/(dz*dz)
        idz_sqrd_by_12 = idz_sqrd / 12.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          DzDzu(i,j,k) = (   2.d0 * u(i,j,k  )     &     
     &                     - 5.d0 * u(i,j,k+1)     &     
     &                     + 4.d0 * u(i,j,k+2)     &     
     &                     -        u(i,j,k+3)     &     
     &                   ) * idz_sqrd
          k=2
          DzDzu(i,j,k) = (          u(i,j,k-1)     &     
     &                     - 2.d0 * u(i,j,k  )     &     
     &                     +        u(i,j,k+1)     &     
     &                   ) * idz_sqrd


	   k=3
	   DzDzu(i,j,k) = (-4.d0  *u(i,j,k-2)      &
     &                     +59.d0 *u(i,j,k-1)      &
     &                     -110.d0*u(i  ,j,k)      &
     &                     +59.d0 *u(i,j,k+1)      &
     &                     -4.d0  *u(i,j,k+2)      &
     &                    ) * idz_sqrd/43.d0

	   k=4
	   DzDzu(i,j,k) = (-        u(i,j,k-3)         &
     &                     +  59.d0*u(i,j,k-1)         &
     &                     - 118.d0*u(i  ,j,k)         &
     &                     + 64.d0 *u(i,j,k+1)         &
     &                     -  4.d0 *u(i,j,k+2)         &
     &                    ) * idz_sqrd/49.d0


          do k = 3, nz-2
            DzDzu(i,j,k) = ( -         u(i,j,k-2)     &     
     &                       + 16.d0 * u(i,j,k-1)     &     
     &                       - 30.d0 * u(i,j,k  )     &     
     &                       + 16.d0 * u(i,j,k+1)     &     
     &                       -         u(i,j,k+2)     &     
     &                     ) * idz_sqrd_by_12
          end do

	  k=nz-3
	   DzDzu(i,j,k) = (-        u(i,j,k+3)         &
     &                     +  59.d0*u(i,j,k+1)         &
     &                     - 118.d0*u(i  ,j,k)         &
     &                     + 64.d0 *u(i,j,k-1)         &
     &                     -  4.d0 *u(i,j,k-2)         &
     &                    ) * idz_sqrd/49.d0 

          k=nz-2
	   DzDzu(i,j,k) = (-4.d0  *u(i,j,k+2)      &
     &                     +59.d0 *u(i,j,k+1)      &
     &                     -110.d0*u(i  ,j,k)      &
     &                     +59.d0 *u(i,j,k-1)      &
     &                     -4.d0  *u(i,j,k-2)      &
     &                    ) * idz_sqrd/43.d0


          k=nz-1 
          DzDzu(i,j,k) = (          u(i,j,k-1)     &     
     &                     - 2.d0 * u(i,j,k  )     &     
     &                     +        u(i,j,k+1)     &     
     &                   ) * idz_sqrd
          k=nz
          DzDzu(i,j,k) = ( -        u(i,j,k-3)     &     
     &                     + 4.d0 * u(i,j,k-2)     &      
     &                     - 5.d0 * u(i,j,k-1)     &      
     &                     + 2.d0 * u(i,j,k  )     &      
     &                   ) * idz_sqrd
        end do
        end do
     
        return          
      end subroutine  
          
      

      !----------------------------------------------------------------
      !                   
      !  deriv44_xx  calculates the second derivative with respect to  
      !              x using a centered fourth order finite difference 
      !              operator at all points except the boundaries (i=1
      !              and i=nx) and the penultimate boundary points 
      !              (i=2 and i=nx-1).  At the penultimate boundary 
      !              points, we use once-shifted fourth order operators.   
      !              Both the forward and backward versions have six
      !              point stencils.  At the boundary points themselves,  
      !              we use completely shifted fourth order operators, 
      !              again with six point stencils. 
      ! 
      !----------------------------------------------------------------      
      subroutine deriv44_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz     
        CCTK_REAL DxDxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx    
        CCTK_INT   d_type, dd_type 
     
        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx_sqrd, idx_sqrd_by_12
        
        idx_sqrd = 1.0d0/(dx*dx)
        idx_sqrd_by_12 = idx_sqrd / 12.d0
      
        do k = 1, nz
        do j = 1, ny
          i=1
          DxDxu(i,j,k) = (    45.d0 * u(i  ,j,k)     &     
     &                     - 154.d0 * u(i+1,j,k)     &     
     &                     + 214.d0 * u(i+2,j,k)     &     
     &                     - 156.d0 * u(i+3,j,k)     &               
     &                     +  61.d0 * u(i+4,j,k)     &      
     &                     -  10.d0 * u(i+5,j,k)     &     
     &                   ) * idx_sqrd_by_12
          i=2
          DxDxu(i,j,k) = (    10.d0 * u(i-1,j,k)     &     
     &                      - 15.d0 * u(i  ,j,k)     &     
     &                      -  4.d0 * u(i+1,j,k)     &     
     &                      + 14.d0 * u(i+2,j,k)     &     
     &                      -  6.d0 * u(i+3,j,k)     &     
     &                      +         u(i+4,j,k)     &     
     &                   ) * idx_sqrd_by_12
          do i = 3, nx-2 
            DxDxu(i,j,k) = ( -         u(i-2,j,k)     &     
     &                       + 16.d0 * u(i-1,j,k)     &     
     &                       - 30.d0 * u(i  ,j,k)     &     
     &                       + 16.d0 * u(i+1,j,k)     &     
     &                       -         u(i+2,j,k)     &     
     &                     ) * idx_sqrd_by_12
          end do
          i=nx-1        
          DxDxu(i,j,k) = (            u(i-4,j,k)     &     
     &                      -  6.d0 * u(i-3,j,k)     &     
     &                      + 14.d0 * u(i-2,j,k)     &     
     &                      -  4.d0 * u(i-1,j,k)     &     
     &                      - 15.d0 * u(i  ,j,k)     &     
     &                      + 10.d0 * u(i+1,j,k)     &     
     &                   ) * idx_sqrd_by_12
          i=nx
          DxDxu(i,j,k) = ( -  10.d0 * u(i-5,j,k)     &     
     &                     +  61.d0 * u(i-4,j,k)     &     
     &                     - 156.d0 * u(i-3,j,k)     &     
     &                     + 214.d0 * u(i-2,j,k)     &     
     &                     - 154.d0 * u(i-1,j,k)     &     
     &                     +  45.d0 * u(i  ,j,k)     &     
     &                   ) * idx_sqrd_by_12
        end do
        end do

        return
      end subroutine


      !----------------------------------------------------------------
      !
      !  deriv44_yy  calculates the second derivative with respect to  
      !              y using a centered fourth order finite difference 
      !              operator at all points except the boundaries (j=1
      !              and j=ny) and the penultimate boundary points 
      !              (j=2 and j=ny-1).  At the penultimate boundary 
      !              points, we use once-shifted fourth order operators.   
      !              Both the forward and backward versions have six
      !              point stencils.  At the boundary points themselves,  
      !              we use completely shifted fourth order operators, 
      !              again with six point stencils. 
      !
      !----------------------------------------------------------------      
      subroutine deriv44_yy(DyDyu, u, dy, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idy_sqrd, idy_sqrd_by_12

        idy_sqrd = 1.0d0/(dy*dy)
        idy_sqrd_by_12 = idy_sqrd / 12.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          DyDyu(i,j,k) = (    45.d0 * u(i,j  ,k)     &     
     &                     - 154.d0 * u(i,j+1,k)     &     
     &                     + 214.d0 * u(i,j+2,k)     &     
     &                     - 156.d0 * u(i,j+3,k)     &     
     &                     +  61.d0 * u(i,j+4,k)     &     
     &                     -  10.d0 * u(i,j+5,k)     &     
     &                   ) * idy_sqrd_by_12
          j=2
          DyDyu(i,j,k) = (    10.d0 * u(i,j-1,k)     &     
     &                      - 15.d0 * u(i,j  ,k)     &     
     &                      -  4.d0 * u(i,j+1,k)     &     
     &                      + 14.d0 * u(i,j+2,k)     &     
     &                      -  6.d0 * u(i,j+3,k)     &     
     &                      +         u(i,j+4,k)     &     
     &                   ) * idy_sqrd_by_12
          do j = 3, ny-2
            DyDyu(i,j,k) = ( -         u(i,j-2,k)     &     
     &                       + 16.d0 * u(i,j-1,k)     &     
     &                       - 30.d0 * u(i,j  ,k)     &     
     &                       + 16.d0 * u(i,j+1,k)     &     
     &                       -         u(i,j+2,k)     &     
     &                     ) * idy_sqrd_by_12
          end do
          j=ny-1
          DyDyu(i,j,k) = (            u(i,j-4,k)     &     
     &                      -  6.d0 * u(i,j-3,k)     &     
     &                      + 14.d0 * u(i,j-2,k)     &     
     &                      -  4.d0 * u(i,j-1,k)     &     
     &                      - 15.d0 * u(i,j  ,k)     &     
     &                      + 10.d0 * u(i,j+1,k)     &     
     &                   ) * idy_sqrd_by_12
          j=ny           
          DyDyu(i,j,k) = ( -  10.d0 * u(i,j-5,k)     &     
     &                     +  61.d0 * u(i,j-4,k)     &     
     &                     - 156.d0 * u(i,j-3,k)     &     
     &                     + 214.d0 * u(i,j-2,k)     &     
     &                     - 154.d0 * u(i,j-1,k)     &     
     &                     +  45.d0 * u(i,j  ,k)     &     
     &                   ) * idy_sqrd_by_12  
        end do          
        end do          
                      
        return
      end subroutine


      !----------------------------------------------------------------
      !
      !  deriv44_zz  calculates the second derivative with respect to  
      !              z using a centered fourth order finite difference 
      !              operator at all points except the boundaries (k=1
      !              and k=nz) and the penultimate boundary points 
      !              (k=2 and k=nz-1).  At the penultimate boundary 
      !              points, we use once-shifted fourth order operators.   
      !              Both the forward and backward versions have six
      !              point stencils.  At the boundary points themselves,  
      !              we use completely shifted fourth order operators, 
      !              again with six point stencils. 
      !
      !----------------------------------------------------------------      
      subroutine deriv44_zz(DzDzu, u, dz, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DzDzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 
        
        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idz_sqrd, idz_sqrd_by_12
        
        idz_sqrd = 1.0d0/(dz*dz)
        idz_sqrd_by_12 = idz_sqrd / 12.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          DzDzu(i,j,k) = (    45.d0 * u(i,j,k  )     &     
     &                     - 154.d0 * u(i,j,k+1)     &     
     &                     + 214.d0 * u(i,j,k+2)     &     
     &                     - 156.d0 * u(i,j,k+3)     &     
     &                     +  61.d0 * u(i,j,k+4)     &     
     &                     -  10.d0 * u(i,j,k+5)     &     
     &                   ) * idz_sqrd_by_12  
          k=2           
          DzDzu(i,j,k) = (    10.d0 * u(i,j,k-1)     &     
     &                      - 15.d0 * u(i,j,k  )     &     
     &                      -  4.d0 * u(i,j,k+1)     &     
     &                      + 14.d0 * u(i,j,k+2)     &     
     &                      -  6.d0 * u(i,j,k+3)     &     
     &                      +         u(i,j,k+4)     &     
     &                   ) * idz_sqrd_by_12
          do k = 3, nz-2 
            DzDzu(i,j,k) = ( -         u(i,j,k-2)     &     
     &                       + 16.d0 * u(i,j,k-1)     &     
     &                       - 30.d0 * u(i,j,k  )     &     
     &                       + 16.d0 * u(i,j,k+1)     &     
     &                       -         u(i,j,k+2)     &     
     &                     ) * idz_sqrd_by_12 
          end do          
          k=nz-1          
          DzDzu(i,j,k) = (            u(i,j,k-4)     &     
     &                      -  6.d0 * u(i,j,k-3)     &     
     &                      + 14.d0 * u(i,j,k-2)     &     
     &                      -  4.d0 * u(i,j,k-1)     &     
     &                      - 15.d0 * u(i,j,k  )     &     
     &                      + 10.d0 * u(i,j,k+1)     &     
     &                   ) * idz_sqrd_by_12
          k=nz
          DzDzu(i,j,k) = ( -  10.d0 * u(i,j,k-5)     &     
     &                     +  61.d0 * u(i,j,k-4)     &     
     &                     - 156.d0 * u(i,j,k-3)     &     
     &                     + 214.d0 * u(i,j,k-2)     &     
     &                     - 154.d0 * u(i,j,k-1)     &     
     &                     +  45.d0 * u(i,j,k  )     &     
     &                   ) * idz_sqrd_by_12
        end do
        end do

        return
      end subroutine





!      !----------------------------------------------------------------
!      !
!      !
!      !
!      !----------------------------------------------------------------
!      subroutine deriv43sbp_x(Dxu, u, dx, nx, ny, nz)
!        implicit   none
!        CCTK_INT  nx, ny, nz
!        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz)
!        CCTK_REAL dx
!
!        ! local vars
!        CCTK_INT  i, j, k
!        CCTK_REAL idx
!
!        idx = 1.0d0/dx
!        do k = 1, nz
!        do j = 1, ny
!          Dxu(1,j,k) = (d00*u(1,j,k)     &     + d01*u(2,j,k)     &     +d02*u(3,j,k)     &      &
!     &                  +d03*u(4,j,k)     &     )*idx
!          Dxu(2,j,k) = (d10*u(1,j,k)     &     +d11*u(2,j,k)     &     +d12*u(3,j,k)     &      &
!     &                  +d13*u(4,j,k)     &      +d14*u(5,j,k)     &     +d15*u(6,j,k)     &     )*idx
!          Dxu(3,j,k) = (d20*u(1,j,k)     &     +d21*u(2,j,k)     &     +d22*u(3,j,k)     &      &
!     &                  +d23*u(4,j,k)     &      +d24*u(5,j,k)     &     +d25*u(6,j,k)     &     )*idx
!          Dxu(4,j,k) = (d30*u(1,j,k)     &     +d31*u(2,j,k)     &     +d32*u(3,j,k)     &      &
!     &                  +d33*u(4,j,k)     &      +d34*u(5,j,k)     &     +d35*u(6,j,k)     &      &
!     &                  +d36*u(7,j,k)     &     )*idx
!          Dxu(5,j,k) = (d40*u(1,j,k)     &     +d41*u(2,j,k)     &     +d42*u(3,j,k)     &      &
!     &                  +d43*u(4,j,k)     &      +d44*u(5,j,k)     &     +d45*u(6,j,k)     &      &
!     &                  +d46*u(7,j,k)     &     )*idx
!    
!          do i = 6, nx-5 
!            Dxu(i,j,k) = (-d2*u(i-2,j,k)     &     -d1*u(i-1,j,k)     &     +d1*u(i+1,j,k)     &      &
!     &                  +d2*u(i+2,j,k)     &     )*idx
!          end do
!         
!          Dxu(nx-4,j,k) = -(d40*u(nx,j,k)     &     +d41*u(nx-1,j,k)     &      &
!     &                    +d42*u(nx-2,j,k)     &     +d43*u(nx-3,j,k)     &      &
!     &                    +d44*u(nx-4,j,k)     &     +d45*u(nx-5,j,k)     &      &
!     &                    +d46*u(nx-6,j,k)     &     )*idx
!          Dxu(nx-3,j,k) = -(d30*u(nx,j,k)     &     +d31*u(nx-1,j,k)     &      &
!     &                     +d32*u(nx-2,j,k)     &     +d33*u(nx-3,j,k)     &      &
!     &                     +d34*u(nx-4,j,k)     &     +d35*u(nx-5,j,k)     &      &
!     &                     +d36*u(nx-6,j,k)     &     )*idx
!          Dxu(nx-2,j,k)     &       = -(d20*u(nx,j,k)     &     +d21*u(nx-1,j,k)     &      &
!     &                     +d22*u(nx-2,j,k)     &     +d23*u(nx-3,j,k)     &      &
!     &                     +d24*u(nx-4,j,k)     &     +d25*u(nx-5,j,k)     &     )*idx
!          Dxu(nx-1,j,k) = -(d10*u(nx,j,k)     &     +d11*u(nx-1,j,k)     &      &
!     &                    +d12*u(nx-2,j,k)     &     +d13*u(nx-3,j,k)     &      &
!     &                    +d14*u(nx-4,j,k)     &     +d15*u(nx-5,j,k)     &     )*idx
!          Dxu(nx,j,k) = -(d00*u(nx,j,k)     &     +d01*u(nx-1,j,k)     &      &
!     &                    +d02*u(nx-2,j,k)     &     +d03*u(nx-3,j,k)     &     )*idx
!    
!        end do
!        end do
!    
!        return         
!      end subroutine   
!    
!    
!      !----------------------------------------------------------------
!      !              
!      !
!      !
!      !----------------------------------------------------------------
!      subroutine deriv43sbp_y(Dyu, u, dy, nx, ny, nz)
!        implicit   none
!        CCTK_INT  nx, ny, nz
!        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
!        CCTK_REAL dy
!
!        ! local vars
!        CCTK_INT  i, j, k
!        CCTK_REAL idy
!
!        idy = 1.0d0/dy
!
!        do k = 1, nz
!        do i = 1, nx
!          Dyu(i,1,k) = (d00*u(i,1,k)     &     + d01*u(i,2,k)     &     +d02*u(i,3,k)     &      &
!     &           +d03*u(i,4,k)     &     )*idy
!          Dyu(i,2,k) = (d10*u(i,1,k)     &     +d11*u(i,2,k)     &     +d12*u(i,3,k)     &      &
!     &           +d13*u(i,4,k)     &       +d14*u(i,5,k)     &     +d15*u(i,6,k)     &     )*idy
!          Dyu(i,3,k) = (d20*u(i,1,k)     &     +d21*u(i,2,k)     &     +d22*u(i,3,k)     &       &
!     &           +d23*u(i,4,k)     &      +d24*u(i,5,k)     &     +d25*u(i,6,k)     &     )*idy
!          Dyu(i,4,k) = (d30*u(i,1,k)     &     +d31*u(i,2,k)     &     +d32*u(i,3,k)     &       &
!     &           +d33*u(i,4,k)     &      +d34*u(i,5,k)     &     +d35*u(i,6,k)     &       &
!     &           +d36*u(i,7,k)     &     )*idy
!          Dyu(i,5,k) = (d40*u(i,1,k)     &     +d41*u(i,2,k)     &     +d42*u(i,3,k)     &       &
!     &           +d43*u(i,4,k)     &      +d44*u(i,5,k)     &     +d45*u(i,6,k)     &       &
!     &           +d46*u(i,7,k)     &     )*idy
!
!          do j = 6, ny-5
!            Dyu(i,j,k) = (-d2*u(i,j-2,k)     &     -d1*u(i,j-1,k)     &     +d1*u(i,j+1,k)     &       &
!     &           +d2*u(i,j+2,k)     &     )*idy
!          end do
!
!          Dyu(i,ny-4,k) = -(d40*u(i,ny,k)     &     +d41*u(i,ny-1,k)     &       &
!     &           +d42*u(i,ny-2,k)     &     +d43*u(i,ny-3,k)     &       &
!     &           +d44*u(i,ny-4,k)     &     +d45*u(i,ny-5,k)     &     +d46*u(i,ny-6,k)     &     )*idy
!          Dyu(i,ny-3,k) = -(d30*u(i,ny,k)     &     +d31*u(i,ny-1,k)     &       &
!     &           +d32*u(i,ny-2,k)     &     +d33*u(i,ny-3,k)     &       &
!     &           +d34*u(i,ny-4,k)     &     +d35*u(i,ny-5,k)     &     +d36*u(i,ny-6,k)     &     )*idy
!          Dyu(i,ny-2,k) = -(d20*u(i,ny,k)     &     +d21*u(i,ny-1,k)     &       &
!     &           +d22*u(i,ny-2,k)     &     +d23*u(i,ny-3,k)     &     +d24*u(i,ny-4,k)     &       &
!     &           +d25*u(i,ny-5,k)     &     )*idy
!          Dyu(i,ny-1,k) = -(d10*u(i,ny,k)     &     +d11*u(i,ny-1,k)     &       &
!     &            +d12*u(i,ny-2,k)     &     +d13*u(i,ny-3,k)     &     +d14*u(i,ny-4,k)     &       &
!     &            +d15*u(i,ny-5,k)     &     )*idy
!          Dyu(i,ny,k) = -(d00*u(i,ny,k)     &     +d01*u(i,ny-1,k)     &       &
!     &            +d02*u(i,ny-2,k)     &     +d03*u(i,ny-3,k)     &     )*idy
!        end do
!        end do
!
!        return
!      end subroutine

!      !----------------------------------------------------------------
!      !
!      !
!      !
!      !----------------------------------------------------------------
!      subroutine deriv43sbp_z(Dzu, u, dz, nx, ny, nz)
!        implicit   none
!        CCTK_INT  nx, ny, nz
!        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
!        CCTK_REAL dz
!
!        ! local vars
!        CCTK_INT  i, j, k
!        CCTK_REAL idz
!
!        idz = 1.0d0/dz
!
!        do j = 1, ny
!        do i = 1, nx
!          Dzu(i,j,1) = (d00*u(i,j,1)+ d01*u(i,j,2)+d02*u(i,j,3)  &
!     &          +d03*u(i,j,4))*idz
!          Dzu(i,j,2) = (d10*u(i,j,1)+d11*u(i,j,2)+d12*u(i,j,3)  &
!     &                 +d13*u(i,j,4) +d14*u(i,j,5)+d15*u(i,j,6))*idz
!          Dzu(i,j,3) = (d20*u(i,j,1)+d21*u(i,j,2)+d22*u(i,j,3) &
!     &                 +d23*u(i,j,4) +d24*u(i,j,5)+d25*u(i,j,6))*idz
!          Dzu(i,j,4) = (d30*u(i,j,1)+d31*u(i,j,2)+d32*u(i,j,3)
!     &                 +d33*u(i,j,4) +d34*u(i,j,5)+d35*u(i,j,6)
!     &                 +d36*u(i,j,7))*idz
!          Dzu(i,j,5) = (d40*u(i,j,1)+d41*u(i,j,2)+d42*u(i,j,3)
!     &                 +d43*u(i,j,4) +d44*u(i,j,5)+d45*u(i,j,6)
!     &                 +d46*u(i,j,7))*idz
!          do k = 6, nz-5
!            Dzu(i,j,k) = (-d2*u(i,j,k-2)     &     -d1*u(i,j,k-1)     &     +d1*u(i,j,k+1)     &     
!     &                    +d2*u(i,j,k+2)     &     )*idz
!          end do
!          Dzu(i,j,nz-4) = -(d40*u(i,j,nz)+d41*u(i,j,nz-1)
!     &                     +d42*u(i,j,nz-2)+d43*u(i,j,nz-3)
!     &                     +d44*u(i,j,nz-4)+d45*u(i,j,nz-5)
!     &                     +d46*u(i,j,nz-6))*idz
!          Dzu(i,j,nz-3) = -(d30*u(i,j,nz)+d31*u(i,j,nz-1)
!     &                     +d32*u(i,j,nz-2)+d33*u(i,j,nz-3)
!     &                     +d34*u(i,j,nz-4)+d35*u(i,j,nz-5)
!     &                     +d36*u(i,j,nz-6))*idz
!          Dzu(i,j,nz-2) = -(d20*u(i,j,nz)+d21*u(i,j,nz-1)
!     &                     +d22*u(i,j,nz-2)+d23*u(i,j,nz-3)
!     &                     +d24*u(i,j,nz-4)+d25*u(i,j,nz-5))*idz
!          Dzu(i,j,nz-1) = -(d10*u(i,j,nz)+d11*u(i,j,nz-1)
!     &                     +d12*u(i,j,nz-2)+d13*u(i,j,nz-3)
!     &                     +d14*u(i,j,nz-4)+d15*u(i,j,nz-5))*idz
!          Dzu(i,j,nz) = -(d00*u(i,j,nz)+d01*u(i,j,nz-1)
!     &                     +d02*u(i,j,nz-2)+d03*u(i,j,nz-3))*idz
!
!        end do
!        end do
!
!        return
!      end subroutine







!======================================================================! 
!                                                                      ! 
!                   SIXTH ORDER DERIVATIVE OPERATORS                   ! 
!                                                                      ! 
!======================================================================! 

      !----------------------------------------------------------------
      !
      !  deriv642_x  calculates the first derivative with respect to 
      !              x using a centered sixth order finite difference 
      !              operator on the interval [4,nx-3].  At the points 
      !              i=3 and i=nx-2, a centered fourth order operator
      !              is used.  At i=2 and i=nx-1, we use a centered 
      !              second order operator.  
      !     
      !              At the boundary points (i=1 and i=nx) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at i=1, we use a forward (upwind)  
      !              operator with a three point stencil.  At i=nx, we  
      !              a backward (downwind) operator, again with a three  
      !              point stencil. 
      !
      !----------------------------------------------------------------      
      subroutine deriv642_x(Dxu, u, dx, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60 

        idx = 1.0d0/dx        
        idx_by_2 = idx / 2.d0
        idx_by_12 = idx / 12.d0
        idx_by_60 = idx / 60.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  3.d0 * u(i  ,j,k)     &      
     &                   +  4.d0 * u(i+1,j,k)     &     
     &                   -         u(i+2,j,k)     &     
     &                 ) * idx_by_2
          i=2
          Dxu(i,j,k) = ( -         u(i-1,j,k)     &      
     &                   +         u(i+1,j,k)     &     
     &                 ) * idx_by_2
          i=3
          Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                   -  8.d0 * u(i-1,j,k)     &     
     &                   +  8.d0 * u(i+1,j,k)     &      
     &                   -         u(i+2,j,k)     &      
     &                 ) * idx_by_12 
          do i = 4, nx-3
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  9.d0 * u(i-2,j,k)     &     
     &                     - 45.d0 * u(i-1,j,k)     &     
     &                     + 45.d0 * u(i+1,j,k)     &     
     &                     -  9.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_60
          end do
          i=nx-2 
          Dxu(i,j,k) = (          u(i-2,j,k)     &      
     &                   - 8.d0 * u(i-1,j,k)     &     
     &                   + 8.d0 * u(i+1,j,k)     &      
     &                   -        u(i+2,j,k)     &      
     &                 ) * idx_by_12
          i=nx-1 
          Dxu(i,j,k) = ( -        u(i-1,j,k)     &      
     &                   +        u(i+1,j,k)     &      
     &                 ) * idx_by_2
          i=nx 
          Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                    - 4.d0 * u(i-1,j,k)     &      
     &                    + 3.d0 * u(i  ,j,k)     &      
     &                 ) * idx_by_2 
        end do
        end do
       
        return
      end subroutine 
 

      !----------------------------------------------------------------
      !
      !  deriv642_y  calculates the first derivative with respect to 
      !              y using a centered sixth order finite difference 
      !              operator on the interval [4,ny-3].  At the points 
      !              j=3 and j=ny-2, a centered fourth order operator
      !              is used.  At j=2 and j=ny-1, we use a centered 
      !              second order operator.  
      !     
      !              At the boundary points (j=1 and j=ny) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at j=1, we use a forward (upwind)  
      !              operator with a three point stencil.  At j=ny, we  
      !              a backward (downwind) operator, again with a three  
      !              point stencil. 
      !
      !----------------------------------------------------------------      
      subroutine deriv642_y(Dyu, u, dy, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12, idy_by_60 

        idy = 1.0d0/dy        
        idy_by_2 = idy / 2.d0
        idy_by_12 = idy / 12.d0
        idy_by_60 = idy / 60.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  3.d0 * u(i,j  ,k)     &      
     &                   +  4.d0 * u(i,j+1,k)     &     
     &                   -         u(i,j+2,k)     &     
     &                 ) * idy_by_2
          j=2
          Dyu(i,j,k) = ( -         u(i,j-1,k)     &      
     &                   +         u(i,j+1,k)     &     
     &                 ) * idy_by_2
          j=3
          Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                   -  8.d0 * u(i,j-1,k)     &     
     &                   +  8.d0 * u(i,j+1,k)     &      
     &                   -         u(i,j+2,k)     &      
     &                 ) * idy_by_12 
          do j = 4, ny-3
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  9.d0 * u(i,j-2,k)     &     
     &                     - 45.d0 * u(i,j-1,k)     &     
     &                     + 45.d0 * u(i,j+1,k)     &     
     &                     -  9.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_60
          end do
          j=ny-2 
          Dyu(i,j,k) = (          u(i,j-2,k)     &      
     &                   - 8.d0 * u(i,j-1,k)     &     
     &                   + 8.d0 * u(i,j+1,k)     &      
     &                   -        u(i,j+2,k)     &      
     &                 ) * idy_by_12
          j=ny-1 
          Dyu(i,j,k) = ( -        u(i,j-1,k)     &      
     &                   +        u(i,j+1,k)     &      
     &                 ) * idy_by_2
          j=ny 
          Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                    - 4.d0 * u(i,j-1,k)     &      
     &                    + 3.d0 * u(i,j  ,k)     &      
     &                 ) * idy_by_2 
        end do
        end do
       
        return
      end subroutine 
 

      !----------------------------------------------------------------
      !
      !  deriv642_z  calculates the first derivative with respect to 
      !              z using a centered sixth order finite difference 
      !              operator on the interval [4,nz-3].  At the points 
      !              k=3 and k=nz-2, a centered fourth order operator
      !              is used.  At k=2 and k=nz-1, we use a centered 
      !              second order operator.  
      !   
      !              At the boundary points (k=1 and k=nz) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at k=1, we use a forward (upwind)  
      !              operator with a three point stencil.  At k=nz, we  
      !              a backward (downwind) operator, again with a three  
      !              point stencil. 
      !
      !----------------------------------------------------------------      
      subroutine deriv642_z(Dzu, u, dz, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idz, idz_by_2, idz_by_12, idz_by_60 

        idz = 1.0d0/dz
        idz_by_2 = idz / 2.d0
        idz_by_12 = idz / 12.d0
        idz_by_60 = idz / 60.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( - 3.d0 * u(i,j,k  )     &     
     &                   + 4.d0 * u(i,j,k+1)     &     
     &                   -        u(i,j,k+2)     &     
     &                 ) * idz_by_2
          k=2
          Dzu(i,j,k) = ( -         u(i,j,k-1)     &     
     &                   +         u(i,j,k+1)     &     
     &                 ) * idz_by_2
          k=3
          Dzu(i,j,k) = (           u(i,j,k-2)     &     
     &                   -  8.d0 * u(i,j,k-1)     &     
     &                   +  8.d0 * u(i,j,k+1)     &     
     &                   -         u(i,j,k+2)     &     
     &                 ) * idz_by_12
          do k = 4, nz-3
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &     
     &                     +  9.d0 * u(i,j,k-2)     &     
     &                     - 45.d0 * u(i,j,k-1)     &     
     &                     + 45.d0 * u(i,j,k+1)     &     
     &                     -  9.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_60
          end do
          k=nz-2
          Dzu(i,j,k) = (          u(i,j,k-2)     &     
     &                   - 8.d0 * u(i,j,k-1)     &     
     &                   + 8.d0 * u(i,j,k+1)     &     
     &                   -        u(i,j,k+2)     &     
     &                 ) * idz_by_12
          k=nz-1
          Dzu(i,j,k) = ( -        u(i,j,k-1)     &     
     &                   +        u(i,j,k+1)     &     
     &                 ) * idz_by_2
          k=nz
          Dzu(i,j,k) = (          u(i,j,k-2)     &     
     &                   - 4.d0 * u(i,j,k-1)     &     
     &                   + 3.d0 * u(i,j,k  )     &     
     &                 ) * idz_by_2
        end do
        end do

        return
      end subroutine


        !!  a = sgn(beta(i,j,k)     &     ) 
    !!  
    !!  Du(i,j,k) = ( a*u(i+3*a,j,k)     &      - 6*a*u(i+2*a,j,k)     &      + 18*a*u(i+a,j,k)     &      - a*u)i,j,k)     &      - 3*a*u(i-a,j,k)     &       
     !!   return
      !!end subroutine




      !----------------------------------------------------------------
      !
      !  deriv666_x  calculates the first derivative of an array with  
      !              respect to x using a centered sixth order finite 
      !              difference operator on the interval [4,nx-3].  
      !              At the points i=3 and i=nx-2, a (once shifted) 
      !              upwind (respectively downwind) sixth order 
      !              operator is used.  Similarly, at the points i=2 
      !              and i=nx-1, (twice shifted) sixth order operators 
      !              are used.  And at the boundaries (i=1 and i=nx) 
      !              completely shifted sixth order operators are used.  
      !               
      !              Note that all of these sixth order operators use  
      !              seven point stencils.   
      !               
      !----------------------------------------------------------------      
      subroutine deriv666_x(Dxu, u, dx, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60 

        idx = 1.0d0/dx
        idx_by_60 = idx / 60.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( - 147.d0 * u(i  ,j,k)     &     
     &                   + 360.d0 * u(i+1,j,k)     &     
     &                   - 450.d0 * u(i+2,j,k)     &     
     &                   + 400.d0 * u(i+3,j,k)     &     
     &                   - 225.d0 * u(i+4,j,k)     &     
     &                   +  72.d0 * u(i+5,j,k)     &     
     &                   -  10.d0 * u(i+6,j,k)     &     
     &                 ) * idx_by_60
          i=2
          Dxu(2,j,k) = (  -  10.d0 * u(i-1,j,k)     &     
     &                    -  77.d0 * u(i  ,j,k)     &     
     &                    + 150.d0 * u(i+1,j,k)     &     
     &                    - 100.d0 * u(i+2,j,k)     &     
     &                    +  50.d0 * u(i+3,j,k)     &     
     &                    -  15.d0 * u(i+4,j,k)     &     
     &                    +   2.d0 * u(i+5,j,k)     &     
     &                 ) * idx_by_60
          i=3
          Dxu(3,j,k) = (      2.d0 * u(i-2,j,k)     &     
     &                    -  24.d0 * u(i-1,j,k)     &     
     &                    -  35.d0 * u(i  ,j,k)     &     
     &                    +  80.d0 * u(i+1,j,k)     &     
     &                    -  30.d0 * u(i+2,j,k)     &      
     &                    +   8.d0 * u(i+3,j,k)     &     
     &                    -          u(i+4,j,k)     &     
     &                 ) * idx_by_60
          do i = 4, nx-3
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &     
     &                     +  9.d0 * u(i-2,j,k)     &     
     &                     - 45.d0 * u(i-1,j,k)     &     
     &                     + 45.d0 * u(i+1,j,k)     &     
     &                     -  9.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_60
          end do
          i=nx-2 
          Dxu(i,j,k) = (            u(i-4,j,k)     &     
     &                    -  8.d0 * u(i-3,j,k)     &     
     &                    + 30.d0 * u(i-2,j,k)     &     
     &                    - 80.d0 * u(i-1,j,k)     &     
     &                    + 35.d0 * u(i  ,j,k)     &      
     &                    + 24.d0 * u(i+1,j,k)     &      
     &                    -  2.d0 * u(i+2,j,k)     &      
     &                 ) * idx_by_60
          i=nx-1 
          Dxu(i,j,k) = (  -   2.d0 * u(i-5,j,k)     &      
     &                    +  15.d0 * u(i-4,j,k)     &     
     &                    -  50.d0 * u(i-3,j,k)     &     
     &                    + 100.d0 * u(i-2,j,k)     &     
     &                    - 150.d0 * u(i-1,j,k)     &      
     &                    +  77.d0 * u(i  ,j,k)     &      
     &                    +  10.d0 * u(i+1,j,k)     &     
     &                 ) * idx_by_60
          i=nx 
          Dxu(i,j,k) = (     10.d0 * u(i-6,j,k)     &     
     &                    -  72.d0 * u(i-5,j,k)     &     
     &                    + 225.d0 * u(i-4,j,k)     &     
     &                    - 400.d0 * u(i-3,j,k)     &     
     &                    + 450.d0 * u(i-2,j,k)     &     
     &                    - 360.d0 * u(i-1,j,k)     &     
     &                    + 147.d0 * u(i  ,j,k)     &     
     &                  ) * idx_by_60
        end do
        end do

        return 
      end subroutine
          
          
      !----------------------------------------------------------------
      !                 
      !  deriv666_y  calculates the first derivative of an array with  
      !              respect to y using a centered sixth order finite 
      !              difference operator on the interval [4,ny-3].  
      !              At the points j=3 and j=ny-2, a (once shifted) 
      !              upwind (respectively downwind) sixth order 
      !              operator is used.  Similarly, at the points j=2 
      !              and j=ny-1, (twice shifted) sixth order operators 
      !              are used.  And at the boundaries (j=1 and j=ny) 
      !              completely shifted sixth order operators are used.  
      !                  
      !              Note that all of these sixth order operators use  
      !              seven point stencils.   
      !                  
      !----------------------------------------------------------------      
      subroutine deriv666_y(Dyu, u, dy, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idy, idy_by_60 

        idy = 1.0d0/dy
        idy_by_60 = idy / 60.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( - 147.d0 * u(i,j  ,k)     &     
     &                   + 360.d0 * u(i,j+1,k)     &     
     &                   - 450.d0 * u(i,j+2,k)     &     
     &                   + 400.d0 * u(i,j+3,k)     &     
     &                   - 225.d0 * u(i,j+4,k)     &     
     &                   +  72.d0 * u(i,j+5,k)     &     
     &                   -  10.d0 * u(i,j+6,k)     &     
     &                 ) * idy_by_60
          j=2
          Dyu(i,j,k) = (  -  10.d0 * u(i,j-1,k)     &     
     &                    -  77.d0 * u(i,j  ,k)     &     
     &                    + 150.d0 * u(i,j+1,k)     &     
     &                    - 100.d0 * u(i,j+2,k)     &     
     &                    +  50.d0 * u(i,j+3,k)     &     
     &                    -  15.d0 * u(i,j+4,k)     &     
     &                    +   2.d0 * u(i,j+5,k)     &     
     &                 ) * idy_by_60
          j=3
          Dyu(i,j,k) = (      2.d0 * u(i,j-2,k)     &     
     &                    -  24.d0 * u(i,j-1,k)     &     
     &                    -  35.d0 * u(i,j  ,k)     &     
     &                    +  80.d0 * u(i,j+1,k)     &     
     &                    -  30.d0 * u(i,j+2,k)     &     
     &                    +   8.d0 * u(i,j+3,k)     &     
     &                    -          u(i,j+4,k)     &     
     &                 ) * idy_by_60
          do j = 4, ny-3
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &     
     &                     +  9.d0 * u(i,j-2,k)     &     
     &                     - 45.d0 * u(i,j-1,k)     &     
     &                     + 45.d0 * u(i,j+1,k)     &     
     &                     -  9.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_60
          end do
          j=ny-2
          Dyu(i,j,k) = (            u(i,j-4,k)     &     
     &                    -  8.d0 * u(i,j-3,k)     &     
     &                    + 30.d0 * u(i,j-2,k)     &     
     &                    - 80.d0 * u(i,j-1,k)     &     
     &                    + 35.d0 * u(i,j  ,k)     &     
     &                    + 24.d0 * u(i,j+1,k)     &     
     &                    -  2.d0 * u(i,j+2,k)     &     
     &                 ) * idy_by_60
          j=ny-1
          Dyu(i,j,k) = (  -   2.d0 * u(i,j-5,k)     &     
     &                    +  15.d0 * u(i,j-4,k)     &     
     &                    -  50.d0 * u(i,j-3,k)     &     
     &                    + 100.d0 * u(i,j-2,k)     &     
     &                    - 150.d0 * u(i,j-1,k)     &     
     &                    +  77.d0 * u(i,j  ,k)     &     
     &                    +  10.d0 * u(i,j+1,k)     &     
     &                 ) * idy_by_60
          j=ny
          Dyu(i,j,k) = (    10.d0 * u(i,j-6,k)     &     
     &                   -  72.d0 * u(i,j-5,k)     &     
     &                   + 225.d0 * u(i,j-4,k)     &     
     &                   - 400.d0 * u(i,j-3,k)     &     
     &                   + 450.d0 * u(i,j-2,k)     &     
     &                   - 360.d0 * u(i,j-1,k)     &     
     &                   + 147.d0 * u(i,j  ,k)     &     
     &                 ) * idy_by_60
        end do
        end do
        
        return
      end subroutine
        
          
      !----------------------------------------------------------------
      !                 
      !  deriv666_z  calculates the first derivative of an array with  
      !              respect to z using a centered sixth order finite 
      !              difference operator on the interval [4,nz-3].  
      !              At the points k=3 and k=nz-2, a (once shifted) 
      !              upwind (respectively downwind) sixth order 
      !              operator is used.  Similarly, at the points k=2 
      !              and k=nz-1, (twice shifted) sixth order operators 
      !              are used.  And at the boundaries (k=1 and k=nz) 
      !              completely shifted sixth order operators are used.  
      !                  
      !              Note that all of these sixth order operators use  
      !              seven point stencils.   
      !                  
      !----------------------------------------------------------------      
      subroutine deriv666_z(Dzu, u, dz, nx, ny, nz, d_type)
        implicit   none
        CCTK_INT  nx, ny, nz 
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz     
        CCTK_INT   d_type, dd_type 
     
        ! local vars     
        CCTK_INT i, j, k        
        CCTK_REAL idz   , idz_by_60 
     
        idz = 1.0d0/dz   
        idz_by_60 = idz / 60.d0     
     
        do j = 1, ny      
        do i = 1, nx      
          k=1             
          Dzu(i,j,k) = ( - 147.d0 * u(i,j,k  )     &      
     &                   + 360.d0 * u(i,j,k+1)     &     
     &                   - 450.d0 * u(i,j,k+2)     &     
     &                   + 400.d0 * u(i,j,k+3)     &     
     &                   - 225.d0 * u(i,j,k+4)     &     
     &                   +  72.d0 * u(i,j,k+5)     &     
     &                   -  10.d0 * u(i,j,k+6)     &     
     &                 ) * idz_by_60  
          k=2            
          Dzu(i,j,k) = (  -  10.d0 * u(i,j,k-1)     &     
     &                    -  77.d0 * u(i,j,k  )     &     
     &                    + 150.d0 * u(i,j,k+1)     &     
     &                    - 100.d0 * u(i,j,k+2)     &     
     &                    +  50.d0 * u(i,j,k+3)     &     
     &                    -  15.d0 * u(i,j,k+4)     &     
     &                    +   2.d0 * u(i,j,k+5)     &     
     &                 ) * idz_by_60 
          k=3            
          Dzu(i,j,k) = (      2.d0 * u(i,j,k-2)     &     
     &                    -  24.d0 * u(i,j,k-1)     &     
     &                    -  35.d0 * u(i,j,k  )     &     
     &                    +  80.d0 * u(i,j,k+1)     &     
     &                    -  30.d0 * u(i,j,k+2)     &     
     &                    +   8.d0 * u(i,j,k+3)     &     
     &                    -          u(i,j,k+4)     &     
     &                 ) * idz_by_60
          do k = 4, nz-3
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &     
     &                     +  9.d0 * u(i,j,k-2)     &     
     &                     - 45.d0 * u(i,j,k-1)     &     
     &                     + 45.d0 * u(i,j,k+1)     &     
     &                     -  9.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_60
          end do
          k=nz-2
          Dzu(i,j,k) = (            u(i,j,k-4)     &     
     &                    -  8.d0 * u(i,j,k-3)     &     
     &                    + 30.d0 * u(i,j,k-2)     &     
     &                    - 80.d0 * u(i,j,k-1)     &     
     &                    + 35.d0 * u(i,j,k  )     &     
     &                    + 24.d0 * u(i,j,k+1)     &     
     &                    -  2.d0 * u(i,j,k+2)     &     
     &                 ) * idz_by_60
          k=nz-1
          Dzu(i,j,k) = (  -   2.d0 * u(i,j,k-5)     &     
     &                    +  15.d0 * u(i,j,k-4)     &     
     &                    -  50.d0 * u(i,j,k-3)     &     
     &                    + 100.d0 * u(i,j,k-2)     &     
     &                    - 150.d0 * u(i,j,k-1)     &     
     &                    +  77.d0 * u(i,j,k  )     &     
     &                    +  10.d0 * u(i,j,k+1)     &     
     &                 ) * idz_by_60
          k=nz
          Dzu(i,j,k) = (    10.d0 * u(i,j,k-6)     &     
     &                   -  72.d0 * u(i,j,k-5)     &     
     &                   + 225.d0 * u(i,j,k-4)     &     
     &                   - 400.d0 * u(i,j,k-3)     &     
     &                   + 450.d0 * u(i,j,k-2)     &     
     &                   - 360.d0 * u(i,j,k-1)     &     
     &                   + 147.d0 * u(i,j,k  )     &     
     &                 ) * idz_by_60
        end do
        end do

        return
      end subroutine




      !----------------------------------------------------------------
      !
      ! deriv642adv_x  calculates the first derivative with respect to 
      !                x using an "advective" sixth order finite 
      !                difference operator at all points on the 
      !                interval [5,nx-4].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^x \partial_x in the equations.  
      !                Further, it will depend on the sign of beta^x.  
      !                Namely, if beta^x is positive, we will use a 
      !                sixth order upwind operator while if beta^x 
      !                is negative, we use a sixth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At i=4, if beta^x is positive, we continue to  
      !                use the sixth order upwind operator but if  
      !                beta^x is negative, we use a fourth order downwind  
      !                operator.   
      !                  
      !                We reverse this at i=nx-3 such that if beta^x is  
      !                positive, we use a fourth order upwind operator 
      !                while if beta^x is negative we continue to use  
      !                the sixth order downwind operator.   
      !                 
      !                At i=3, if beta^x is positive, we use the   
      !                fourth order upwind operator and if beta^x is 
      !                negative, we use a second order downwind 
      !                operator.   
      !                  
      !                At i=nx-2 if beta^x is positive, we use a  
      !                fourth order centered operator while if beta^x 
      !                is negative we use the fourth order downwind 
      !                operator.   
      !                 
      !                At i=2, if beta^x is positive, we use the   
      !                second order upwind operator and if beta^x is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At i=nx-1, if beta^x is positive, we use the   
      !                second order centered operator and if beta^x is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At i=1 we use the second order upwind operator 
      !                for both positive and negative values of beta^x.
      !                
      !                At i=nx, we use the second order downwind operator 
      !                for both positive and negative values of beta^x.
      !                 
      !----------------------------------------------------------------      
      subroutine deriv642adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz), betax(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60 

        idx = 1.0d0/dx        
        idx_by_2 = 0.5d0 * idx
        idx_by_12 = idx / 12.d0
        idx_by_60 = idx / 60.d0


        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  3.d0 * u(i  ,j,k)     &      
     &                   +  4.d0 * u(i+1,j,k)     &     
     &                   -         u(i+2,j,k)     &     
     &                 ) * idx_by_2

          i=2
          if (betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( - 3.d0 * u(i  ,j,k)     &      
     &                     + 4.d0 * u(i+1,j,k)     &     
     &                     -        u(i+2,j,k)     &     
     &                   ) * idx_by_2 
          else 
            Dxu(i,j,k) = ( -         u(i-1,j,k)     &      
     &                     +         u(i+1,j,k)     &      
     &                   ) * idx_by_2 
          endif 

          i=3
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                     - 10.d0 * u(i  ,j,k)     &     
     &                     + 18.d0 * u(i+1,j,k)     &     
     &                     -  6.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_12 
          else 
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                     -  4.d0 * u(i-1,j,k)     &     
     &                     +  3.d0 * u(i  ,j,k)     &     
     &                   ) * idx_by_2 
          endif

          i=4
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (    2.d0 * u(i-2,j,k)     &      
     &                     - 24.d0 * u(i-1,j,k)     &     
     &                     - 35.d0 * u(i  ,j,k)     &     
     &                     + 80.d0 * u(i+1,j,k)     &     
     &                     - 30.d0 * u(i+2,j,k)     &     
     &                     +  8.d0 * u(i+3,j,k)     &     
     &                     -         u(i+4,j,k)     &     
     &                   ) * idx_by_60 
          else 
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  6.d0 * u(i-2,j,k)     &     
     &                     - 18.d0 * u(i-1,j,k)     &     
     &                     + 10.d0 * u(i  ,j,k)     &     
     &                     +  3.d0 * u(i+1,j,k)     &     
     &                   ) * idx_by_12 
          endif
          
          do i = 5, nx-4  ! see gr-qc/0505055v2.pdf 
            if ( betax(i,j,k) .ge. 0.d0 ) then
              Dxu(i,j,k) = (    2.d0 * u(i-2,j,k)     &      
     &                       - 24.d0 * u(i-1,j,k)     &     
     &                       - 35.d0 * u(i  ,j,k)     &     
     &                       + 80.d0 * u(i+1,j,k)     &     
     &                       - 30.d0 * u(i+2,j,k)     &     
     &                       +  8.d0 * u(i+3,j,k)     &     
     &                       -         u(i+4,j,k)     &     
     &                     ) * idx_by_60 
            else 
              Dxu(i,j,k) = (           u(i-4,j,k)     &      
     &                       -  8.d0 * u(i-3,j,k)     &     
     &                       + 30.d0 * u(i-2,j,k)     &     
     &                       - 80.d0 * u(i-1,j,k)     &     
     &                       + 35.d0 * u(i  ,j,k)     &     
     &                       + 24.d0 * u(i+1,j,k)     &     
     &                       -  2.d0 * u(i+2,j,k)     &     
     &                     ) * idx_by_60 
            endif 
          end do

          i=nx-3 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                     - 10.d0 * u(i  ,j,k)     &     
     &                     + 18.d0 * u(i+1,j,k)     &     
     &                     -  6.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_12 
          else 
            Dxu(i,j,k) = (           u(i-4,j,k)     &      
     &                     -  8.d0 * u(i-3,j,k)     &     
     &                     + 30.d0 * u(i-2,j,k)     &     
     &                     - 80.d0 * u(i-1,j,k)     &     
     &                     + 35.d0 * u(i  ,j,k)     &     
     &                     + 24.d0 * u(i+1,j,k)     &     
     &                     -  2.d0 * u(i+2,j,k)     &     
     &                   ) * idx_by_60 
          endif 

          i=nx-2 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (  - 3.d0 * u(i  ,j,k)     &      
     &                      + 4.d0 * u(i+1,j,k)     &      
     &                      -        u(i+2,j,k)     &     
     &                   ) * idx_by_2 
          else 
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  6.d0 * u(i-2,j,k)     &     
     &                     - 18.d0 * u(i-1,j,k)     &     
     &                     + 10.d0 * u(i  ,j,k)     &     
     &                     +  3.d0 * u(i+1,j,k)     &     
     &                   ) * idx_by_12 
          endif

          i=nx-1 
          if ( betax(i,j,k) .ge. 0.d0 ) then 
            Dxu(i,j,k) = (  -         u(i-1,j,k)     &      
     &                      +         u(i+1,j,k)     &     
     &                    ) * idx_by_2
          else 
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                      - 4.d0 * u(i-1,j,k)     &      
     &                      + 3.d0 * u(i  ,j,k)     &      
     &                   ) * idx_by_2 
          endif

          i=nx 
          Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                    - 4.d0 * u(i-1,j,k)     &      
     &                    + 3.d0 * u(i  ,j,k)     &      
     &                 ) * idx_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv642adv_y  calculates the first derivative with respect to 
      !                y using an "advective" sixth order finite 
      !                difference operator at all points on the 
      !                interval [5,ny-4].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^y \partial_y in the equations.  
      !                Further, it will depend on the sign of beta^y.  
      !                Namely, if beta^y is positive, we will use a 
      !                sixth order upwind operator while if beta^y 
      !                is negative, we use a sixth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At j=4, if beta^y is positive, we continue to  
      !                use the sixth order upwind operator but if  
      !                beta^y is negative, we use a fourth order downwind  
      !                operator.   
      !                  
      !                We reverse this at j=ny-3 such that if beta^y is  
      !                positive, we use a fourth order upwind operator 
      !                while if beta^y is negative we continue to use  
      !                the sixth order downwind operator.   
      !                 
      !                At j=3, if beta^y is positive, we use the   
      !                fourth order upwind operator and if beta^y is 
      !                negative, we use a second order downwind 
      !                operator.   
      !                  
      !                At j=ny-2 if beta^y is positive, we use a  
      !                fourth order centered operator while if beta^y 
      !                is negative we use the fourth order downwind 
      !                operator.   
      !                 
      !                At j=2, if beta^y is positive, we use the   
      !                second order upwind operator and if beta^y is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At j=ny-1, if beta^y is positive, we use the   
      !                second order centered operator and if beta^y is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At j=1 we use the second order upwind operator 
      !                for both positive and negative values of beta^y.
      !                
      !                At j=ny, we use the second order downwind operator 
      !                for both positive and negative values of beta^y.
      !                 
      !----------------------------------------------------------------
      subroutine deriv642adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz), betay(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12, idy_by_60 

        idy = 1.0d0/dy        
        idy_by_2 = 0.5d0 * idy
        idy_by_12 = idy / 12.d0
        idy_by_60 = idy / 60.d0


        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  3.d0 * u(i,j  ,k)     &      
     &                   +  4.d0 * u(i,j+1,k)     &     
     &                   -         u(i,j+2,k)     &     
     &                 ) * idy_by_2

          j=2
          if (betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( - 3.d0 * u(i,j  ,k)     &      
     &                     + 4.d0 * u(i,j+1,k)     &     
     &                     -        u(i,j+2,k)     &     
     &                   ) * idy_by_2 
          else 
            Dyu(i,j,k) = ( -         u(i,j-1,k)     &      
     &                     +         u(i,j+1,k)     &      
     &                   ) * idy_by_2 
          endif 

          j=3
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                     - 10.d0 * u(i,j  ,k)     &     
     &                     + 18.d0 * u(i,j+1,k)     &     
     &                     -  6.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_12 
          else 
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                     -  4.d0 * u(i,j-1,k)     &     
     &                     +  3.d0 * u(i,j  ,k)     &     
     &                   ) * idy_by_2 
          endif

          j=4
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (    2.d0 * u(i,j-2,k)     &      
     &                     - 24.d0 * u(i,j-1,k)     &     
     &                     - 35.d0 * u(i,j  ,k)     &     
     &                     + 80.d0 * u(i,j+1,k)     &     
     &                     - 30.d0 * u(i,j+2,k)     &     
     &                     +  8.d0 * u(i,j+3,k)     &     
     &                     -         u(i,j+4,k)     &     
     &                   ) * idy_by_60 
          else 
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  6.d0 * u(i,j-2,k)     &     
     &                     - 18.d0 * u(i,j-1,k)     &     
     &                     + 10.d0 * u(i,j  ,k)     &     
     &                     +  3.d0 * u(i,j+1,k)     &     
     &                   ) * idy_by_12 
          endif
          
          do j = 5, ny-4  ! see gr-qc/0505055v2.pdf 
            if ( betay(i,j,k) .ge. 0.d0 ) then
              Dyu(i,j,k) = (    2.d0 * u(i,j-2,k)     &      
     &                       - 24.d0 * u(i,j-1,k)     &     
     &                       - 35.d0 * u(i,j  ,k)     &     
     &                       + 80.d0 * u(i,j+1,k)     &     
     &                       - 30.d0 * u(i,j+2,k)     &     
     &                       +  8.d0 * u(i,j+3,k)     &     
     &                       -         u(i,j+4,k)     &     
     &                     ) * idy_by_60 
            else 
              Dyu(i,j,k) = (           u(i,j-4,k)     &      
     &                       -  8.d0 * u(i,j-3,k)     &     
     &                       + 30.d0 * u(i,j-2,k)     &     
     &                       - 80.d0 * u(i,j-1,k)     &     
     &                       + 35.d0 * u(i,j  ,k)     &     
     &                       + 24.d0 * u(i,j+1,k)     &     
     &                       -  2.d0 * u(i,j+2,k)     &     
     &                     ) * idy_by_60 
            endif 
          end do

          j=ny-3 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                     - 10.d0 * u(i,j  ,k)     &     
     &                     + 18.d0 * u(i,j+1,k)     &     
     &                     -  6.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_12 
          else 
            Dyu(i,j,k) = (           u(i,j-4,k)     &      
     &                     -  8.d0 * u(i,j-3,k)     &     
     &                     + 30.d0 * u(i,j-2,k)     &     
     &                     - 80.d0 * u(i,j-1,k)     &     
     &                     + 35.d0 * u(i,j  ,k)     &     
     &                     + 24.d0 * u(i,j+1,k)     &     
     &                     -  2.d0 * u(i,j+2,k)     &     
     &                   ) * idy_by_60 
          endif 

          j=ny-2 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (  - 3.d0 * u(i,j  ,k)     &      
     &                      + 4.d0 * u(i,j+1,k)     &      
     &                      -        u(i,j+2,k)     &     
     &                   ) * idy_by_2 
          else 
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  6.d0 * u(i,j-2,k)     &     
     &                     - 18.d0 * u(i,j-1,k)     &     
     &                     + 10.d0 * u(i,j  ,k)     &     
     &                     +  3.d0 * u(i,j+1,k)     &     
     &                   ) * idy_by_12 
          endif

          j=ny-1 
          if ( betay(i,j,k) .ge. 0.d0 ) then 
            Dyu(i,j,k) = (  -         u(i,j-1,k)     &      
     &                      +         u(i,j+1,k)     &     
     &                    ) * idy_by_2
          else 
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                      - 4.d0 * u(i,j-1,k)     &      
     &                      + 3.d0 * u(i,j  ,k)     &      
     &                   ) * idy_by_2 
          endif
          
          j=ny 
          Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                    - 4.d0 * u(i,j-1,k)     &      
     &                    + 3.d0 * u(i,j  ,k)     &      
     &                 ) * idy_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv642adv_z  calculates the first derivative with respect to 
      !                z using an "advective" sixth order finite 
      !                difference operator at all points on the 
      !                interval [5,nz-4].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^z \partial_z in the equations.  
      !                Further, it will depend on the sign of beta^z.  
      !                Namely, if beta^z is positive, we will use a 
      !                sixth order upwind operator while if beta^z 
      !                is negative, we use a sixth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At k=4, if beta^z is positive, we continue to  
      !                use the sixth order upwind operator but if  
      !                beta^z is negative, we use a fourth order 
      !                downwind operator.   
      !                  
      !                We reverse this at k=nz-3 such that if beta^z is  
      !                positive, we use a fourth order upwind operator 
      !                while if beta^z is negative we continue to use  
      !                the sixth order downwind operator.   
      !                 
      !                At k=3, if beta^z is positive, we use the   
      !                fourth order upwind operator and if beta^z is 
      !                negative, we use a second order downwind 
      !                operator.   
      !                  
      !                At k=nz-2 if beta^z is positive, we use a  
      !                fourth order centered operator while if beta^z 
      !                is negative we use the fourth order downwind 
      !                operator.   
      !                 
      !                At k=2, if beta^z is positive, we use the   
      !                second order upwind operator and if beta^z is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At k=nz-1, if beta^z is positive, we use the   
      !                second order centered operator and if beta^z is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At k=1 we use the second order upwind operator 
      !                for both positive and negative values of beta^z.
      !                
      !                At k=nz, we use the second order downwind operator 
      !                for both positive and negative values of beta^z.
      !                 
      !----------------------------------------------------------------
      subroutine deriv642adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz), betaz(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12, idz_by_60 

        idz = 1.0d0/dz        
        idz_by_2 = 0.5d0 * idz
        idz_by_12 = idz / 12.d0
        idz_by_60 = idz / 60.d0


        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( -  3.d0 * u(i,j,k  )     &      
     &                   +  4.d0 * u(i,j,k+1)     &     
     &                   -         u(i,j,k+2)     &     
     &                 ) * idz_by_2

          k=2
          if (betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( - 3.d0 * u(i,j,k  )     &      
     &                     + 4.d0 * u(i,j,k+1)     &     
     &                     -        u(i,j,k+2)     &     
     &                   ) * idz_by_2 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-1)     &      
     &                     +         u(i,j,k+1)     &      
     &                   ) * idz_by_2 
          endif 

          k=3
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                     - 10.d0 * u(i,j,k  )     &     
     &                     + 18.d0 * u(i,j,k+1)     &     
     &                     -  6.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_12 
          else 
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                     -  4.d0 * u(i,j,k-1)     &     
     &                     +  3.d0 * u(i,j,k  )     &     
     &                   ) * idz_by_2 
          endif

          k=4
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (    2.d0 * u(i,j,k-2)     &      
     &                     - 24.d0 * u(i,j,k-1)     &     
     &                     - 35.d0 * u(i,j,k  )     &     
     &                     + 80.d0 * u(i,j,k+1)     &     
     &                     - 30.d0 * u(i,j,k+2)     &     
     &                     +  8.d0 * u(i,j,k+3)     &     
     &                     -         u(i,j,k+4)     &     
     &                   ) * idz_by_60 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  6.d0 * u(i,j,k-2)     &     
     &                     - 18.d0 * u(i,j,k-1)     &     
     &                     + 10.d0 * u(i,j,k  )     &     
     &                     +  3.d0 * u(i,j,k+1)     &     
     &                   ) * idz_by_12 
          endif
          
          do k = 5, nz-4  ! see gr-qc/0505055v2.pdf 
            if ( betaz(i,j,k) .ge. 0.d0 ) then
              Dzu(i,j,k) = (    2.d0 * u(i,j,k-2)     &      
     &                       - 24.d0 * u(i,j,k-1)     &     
     &                       - 35.d0 * u(i,j,k  )     &     
     &                       + 80.d0 * u(i,j,k+1)     &     
     &                       - 30.d0 * u(i,j,k+2)     &     
     &                       +  8.d0 * u(i,j,k+3)     &     
     &                       -         u(i,j,k+4)     &     
     &                     ) * idz_by_60 
            else 
              Dzu(i,j,k) = (           u(i,j,k-4)     &      
     &                       -  8.d0 * u(i,j,k-3)     &     
     &                       + 30.d0 * u(i,j,k-2)     &     
     &                       - 80.d0 * u(i,j,k-1)     &     
     &                       + 35.d0 * u(i,j,k  )     &     
     &                       + 24.d0 * u(i,j,k+1)     &     
     &                       -  2.d0 * u(i,j,k+2)     &     
     &                     ) * idz_by_60 
            endif 
          end do

          k=nz-3 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                     - 10.d0 * u(i,j,k  )     &     
     &                     + 18.d0 * u(i,j,k+1)     &     
     &                     -  6.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_12 
          else 
            Dzu(i,j,k) = (           u(i,j,k-4)     &      
     &                     -  8.d0 * u(i,j,k-3)     &     
     &                     + 30.d0 * u(i,j,k-2)     &     
     &                     - 80.d0 * u(i,j,k-1)     &     
     &                     + 35.d0 * u(i,j,k  )     &     
     &                     + 24.d0 * u(i,j,k+1)     &     
     &                     -  2.d0 * u(i,j,k+2)     &     
     &                   ) * idz_by_60
          endif 

          k=nz-2 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (  - 3.d0 * u(i,j,k  )     &      
     &                      + 4.d0 * u(i,j,k+1)     &      
     &                      -        u(i,j,k+2)     &     
     &                   ) * idz_by_2 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  6.d0 * u(i,j,k-2)     &     
     &                     - 18.d0 * u(i,j,k-1)     &     
     &                     + 10.d0 * u(i,j,k  )     &     
     &                     +  3.d0 * u(i,j,k+1)     &     
     &                   ) * idz_by_12 
          endif

          k=nz-1 
          if ( betaz(i,j,k) .ge. 0.d0 ) then 
            Dzu(i,j,k) = (  -         u(i,j,k-1)     &      
     &                      +         u(i,j,k+1)     &     
     &                    ) * idz_by_2
          else 
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                      - 4.d0 * u(i,j,k-1)     &      
     &                      + 3.d0 * u(i,j,k  )     &      
     &                   ) * idz_by_2 
          endif

          k=nz 
          Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                    - 4.d0 * u(i,j,k-1)     &      
     &                    + 3.d0 * u(i,j,k  )     &      
     &                 ) * idz_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv666adv_x  calculates the first derivative with respect to 
      !                x using an "advective" sixth order finite 
      !                difference operator at all points on the 
      !                interval [5,nx-4].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^x \partial_x in the equations.  
      !                Further, it will depend on the sign of beta^x.  
      !                Namely, if beta^x is positive, we will use a 
      !                sixth order upwind operator while if beta^x 
      !                is negative, we use a sixth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At i=4, if beta^x is positive, we continue to  
      !                use the sixth order upwind operator but if  
      !                beta^x is negative, we use the sixth order 
      !                centered operator.   
      !                  
      !                We reverse this at i=nx-3 such that if beta^x is  
      !                positive, we use a sixth order centered operator 
      !                while if beta^x is negative we continue to use  
      !                the sixth order downwind operator.   
      !                 
      !                At i=3, we want a sixth order operator so it 
      !                does not matter what the sign of beta^x is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of  
      !                using the downwind operator for the case of 
      !                negative beta^x, but I do not see any way around 
      !                this if we want to have a sixth order operator.  
      !                 
      !                Likewise, at i=nx-2, the requirement of a sixth
      !                order operator renders the sign of beta^x 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^x is positive.  
      !                 
      !                At i=2 and i=nx-1, we use differently shifted 
      !                sixth order operators irrespective of beta^xs 
      !                sign.   
      !                 
      !                At i=1 and i=nx, we use the totally shifted sixth
      !                order operators irrespective of beta^xs sign.
      !
      !----------------------------------------------------------------      
      subroutine deriv666adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz), betax(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60 

        idx = 1.0d0/dx        
        idx_by_2 = 0.5d0 * idx
        idx_by_12 = idx / 12.d0
        idx_by_60 = idx / 60.d0


        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( - 147.d0 * u(i  ,j,k)     &      
     &                   + 360.d0 * u(i+1,j,k)     &     
     &                   - 450.d0 * u(i+2,j,k)     &     
     &                   + 400.d0 * u(i+3,j,k)     &     
     &                   - 225.d0 * u(i+4,j,k)     &     
     &                   +  72.d0 * u(i+5,j,k)     &     
     &                   -  10.d0 * u(i+6,j,k)     &     
     &                 ) * idx_by_60 

          i=2
          Dxu(i,j,k) = ( -  10.d0 * u(i-1,j,k)     &      
     &                   -  77.d0 * u(i  ,j,k)     &     
     &                   + 150.d0 * u(i+1,j,k)     &     
     &                   - 100.d0 * u(i+2,j,k)     &     
     &                   +  50.d0 * u(i+3,j,k)     &     
     &                   -  15.d0 * u(i+4,j,k)     &     
     &                   +   2.d0 * u(i+5,j,k)     &     
     &                 ) * idx_by_60 

          i=3
          Dxu(i,j,k) = (    2.d0 * u(i-2,j,k)     &      
     &                   - 24.d0 * u(i-1,j,k)     &     
     &                   - 35.d0 * u(i  ,j,k)     &     
     &                   + 80.d0 * u(i+1,j,k)     &     
     &                   - 30.d0 * u(i+2,j,k)     &     
     &                   +  8.d0 * u(i+3,j,k)     &     
     &                   -         u(i+4,j,k)     &     
     &                 ) * idx_by_60  

          i=4
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (    2.d0 * u(i-2,j,k)     &      
     &                     - 24.d0 * u(i-1,j,k)     &     
     &                     - 35.d0 * u(i  ,j,k)     &     
     &                     + 80.d0 * u(i+1,j,k)     &     
     &                     - 30.d0 * u(i+2,j,k)     &     
     &                     +  8.d0 * u(i+3,j,k)     &     
     &                     -         u(i+4,j,k)     &     
     &                   ) * idx_by_60 
          else 
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  9.d0 * u(i-2,j,k)     &     
     &                     - 45.d0 * u(i-1,j,k)     &     
     &                     + 45.d0 * u(i+1,j,k)     &     
     &                     -  9.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_60  
          endif
          
          do i = 5, nx-4  ! see gr-qc/0505055v2.pdf 
            if ( betax(i,j,k) .ge. 0.d0 ) then
              Dxu(i,j,k) = (    2.d0 * u(i-2,j,k)     &      
     &                       - 24.d0 * u(i-1,j,k)     &     
     &                       - 35.d0 * u(i  ,j,k)     &     
     &                       + 80.d0 * u(i+1,j,k)     &     
     &                       - 30.d0 * u(i+2,j,k)     &     
     &                       +  8.d0 * u(i+3,j,k)     &     
     &                       -         u(i+4,j,k)     &     
     &                     ) * idx_by_60 
            else 
              Dxu(i,j,k) = (           u(i-4,j,k)     &      
     &                       -  8.d0 * u(i-3,j,k)     &     
     &                       + 30.d0 * u(i-2,j,k)     &     
     &                       - 80.d0 * u(i-1,j,k)     &     
     &                       + 35.d0 * u(i  ,j,k)     &     
     &                       + 24.d0 * u(i+1,j,k)     &     
     &                       -  2.d0 * u(i+2,j,k)     &     
     &                     ) * idx_by_60 
            endif 
          end do

          i=nx-3 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  9.d0 * u(i-2,j,k)     &     
     &                     - 45.d0 * u(i-1,j,k)     &     
     &                     + 45.d0 * u(i+1,j,k)     &     
     &                     -  9.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_60  
          else 
            Dxu(i,j,k) = (           u(i-4,j,k)     &      
     &                     -  8.d0 * u(i-3,j,k)     &     
     &                     + 30.d0 * u(i-2,j,k)     &     
     &                     - 80.d0 * u(i-1,j,k)     &     
     &                     + 35.d0 * u(i  ,j,k)     &     
     &                     + 24.d0 * u(i+1,j,k)     &     
     &                     -  2.d0 * u(i+2,j,k)     &     
     &                   ) * idx_by_60 
          endif 

          i=nx-2 
          Dxu(i,j,k) = (           u(i-4,j,k)     &      
     &                   -  8.d0 * u(i-3,j,k)     &     
     &                   + 30.d0 * u(i-2,j,k)     &     
     &                   - 80.d0 * u(i-1,j,k)     &     
     &                   + 35.d0 * u(i  ,j,k)     &     
     &                   + 24.d0 * u(i+1,j,k)     &     
     &                   -  2.d0 * u(i+2,j,k)     &     
     &                 ) * idx_by_60 

          i=nx-1 
          Dxu(i,j,k) = ( -   2.d0 * u(i-5,j,k)     &      
     &                   +  15.d0 * u(i-4,j,k)     &     
     &                   -  50.d0 * u(i-3,j,k)     &     
     &                   + 100.d0 * u(i-2,j,k)     &     
     &                   - 150.d0 * u(i-1,j,k)     &     
     &                   +  77.d0 * u(i  ,j,k)     &     
     &                   +  10.d0 * u(i+1,j,k)     &     
     &                 ) * idx_by_60 

          i=nx 
          Dxu(i,j,k) = (    10.d0 * u(i-6,j,k)     &      
     &                   -  72.d0 * u(i-5,j,k)     &     
     &                   + 225.d0 * u(i-4,j,k)     &     
     &                   - 400.d0 * u(i-3,j,k)     &     
     &                   + 450.d0 * u(i-2,j,k)     &     
     &                   - 360.d0 * u(i-1,j,k)     &     
     &                   + 147.d0 * u(i  ,j,k)     &     
     &                 ) * idx_by_60 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv666adv_y  calculates the first derivative with respect to 
      !                y using an "advective" sixth order finite 
      !                difference operator at all points on the 
      !                interval [5,ny-4].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^y \partial_y in the equations.  
      !                Further, it will depend on the sign of beta^y.  
      !                Namely, if beta^y is positive, we will use a 
      !                sixth order upwind operator while if beta^y 
      !                is negative, we use a sixth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At j=4, if beta^y is positive, we continue to  
      !                use the sixth order upwind operator but if  
      !                beta^y is negative, we use the sixth order 
      !                centered operator.   
      !                  
      !                We reverse this at j=ny-3 such that if beta^y is  
      !                positive, we use a sixth order centered operator 
      !                while if beta^y is negative we continue to use  
      !                the sixth order downwind operator.   
      !                 
      !                At j=3, we want a sixth order operator so it 
      !                does not matter what the sign of beta^y is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of  
      !                using the downwind operator for the case of 
      !                negative beta^y, but I do not see any way around 
      !                this if we want to have a sixth order operator.  
      !                 
      !                Likewise, at j=ny-2, the requirement of a sixth
      !                order operator renders the sign of beta^y 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^y is positive.  
      !                 
      !                At j=2 and j=ny-1, we use differently shifted 
      !                sixth order operators irrespective of beta^ys 
      !                sign.   
      !                 
      !                At j=1 and j=ny, we use the totally shifted sixth
      !                order operators irrespective of beta^ys sign.
      !
      !----------------------------------------------------------------      
      subroutine deriv666adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz), betay(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12, idy_by_60 

        idy = 1.0d0/dy        
        idy_by_2 = 0.5d0 * idy
        idy_by_12 = idy / 12.d0
        idy_by_60 = idy / 60.d0


        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( - 147.d0 * u(i,j  ,k)     &      
     &                   + 360.d0 * u(i,j+1,k)     &     
     &                   - 450.d0 * u(i,j+2,k)     &     
     &                   + 400.d0 * u(i,j+3,k)     &     
     &                   - 225.d0 * u(i,j+4,k)     &     
     &                   +  72.d0 * u(i,j+5,k)     &     
     &                   -  10.d0 * u(i,j+6,k)     &     
     &                 ) * idy_by_60 

          j=2
          Dyu(i,j,k) = ( -  10.d0 * u(i,j-1,k)     &      
     &                   -  77.d0 * u(i,j  ,k)     &     
     &                   + 150.d0 * u(i,j+1,k)     &     
     &                   - 100.d0 * u(i,j+2,k)     &     
     &                   +  50.d0 * u(i,j+3,k)     &     
     &                   -  15.d0 * u(i,j+4,k)     &     
     &                   +   2.d0 * u(i,j+5,k)     &     
     &                 ) * idy_by_60 

          j=3
          Dyu(i,j,k) = (    2.d0 * u(i,j-2,k)     &      
     &                   - 24.d0 * u(i,j-1,k)     &     
     &                   - 35.d0 * u(i,j  ,k)     &     
     &                   + 80.d0 * u(i,j+1,k)     &     
     &                   - 30.d0 * u(i,j+2,k)     &     
     &                   +  8.d0 * u(i,j+3,k)     &     
     &                   -         u(i,j+4,k)     &     
     &                 ) * idy_by_60  

          j=4
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (    2.d0 * u(i,j-2,k)     &      
     &                     - 24.d0 * u(i,j-1,k)     &     
     &                     - 35.d0 * u(i,j  ,k)     &     
     &                     + 80.d0 * u(i,j+1,k)     &     
     &                     - 30.d0 * u(i,j+2,k)     &     
     &                     +  8.d0 * u(i,j+3,k)     &     
     &                     -         u(i,j+4,k)     &     
     &                   ) * idy_by_60 
          else 
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  9.d0 * u(i,j-2,k)     &     
     &                     - 45.d0 * u(i,j-1,k)     &     
     &                     + 45.d0 * u(i,j+1,k)     &     
     &                     -  9.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_60  
          endif
          
          do j = 5, ny-4  ! see gr-qc/0505055v2.pdf 
            if ( betay(i,j,k) .ge. 0.d0 ) then
              Dyu(i,j,k) = (    2.d0 * u(i,j-2,k)     &      
     &                       - 24.d0 * u(i,j-1,k)     &     
     &                       - 35.d0 * u(i,j  ,k)     &     
     &                       + 80.d0 * u(i,j+1,k)     &     
     &                       - 30.d0 * u(i,j+2,k)     &     
     &                       +  8.d0 * u(i,j+3,k)     &     
     &                       -         u(i,j+4,k)     &     
     &                     ) * idy_by_60 
            else 
              Dyu(i,j,k) = (           u(i,j-4,k)     &      
     &                       -  8.d0 * u(i,j-3,k)     &     
     &                       + 30.d0 * u(i,j-2,k)     &     
     &                       - 80.d0 * u(i,j-1,k)     &     
     &                       + 35.d0 * u(i,j  ,k)     &     
     &                       + 24.d0 * u(i,j+1,k)     &     
     &                       -  2.d0 * u(i,j+2,k)     &     
     &                     ) * idy_by_60 
            endif 
          end do

          j=ny-3 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  9.d0 * u(i,j-2,k)     &     
     &                     - 45.d0 * u(i,j-1,k)     &     
     &                     + 45.d0 * u(i,j+1,k)     &     
     &                     -  9.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_60  
          else 
            Dyu(i,j,k) = (           u(i,j-4,k)     &      
     &                     -  8.d0 * u(i,j-3,k)     &     
     &                     + 30.d0 * u(i,j-2,k)     &     
     &                     - 80.d0 * u(i,j-1,k)     &     
     &                     + 35.d0 * u(i,j  ,k)     &     
     &                     + 24.d0 * u(i,j+1,k)     &     
     &                     -  2.d0 * u(i,j+2,k)     &     
     &                   ) * idy_by_60 
          endif 

          j=ny-2 
          Dyu(i,j,k) = (           u(i,j-4,k)     &      
     &                   -  8.d0 * u(i,j-3,k)     &     
     &                   + 30.d0 * u(i,j-2,k)     &     
     &                   - 80.d0 * u(i,j-1,k)     &     
     &                   + 35.d0 * u(i,j  ,k)     &     
     &                   + 24.d0 * u(i,j+1,k)     &     
     &                   -  2.d0 * u(i,j+2,k)     &     
     &                 ) * idy_by_60 

          j=ny-1 
          Dyu(i,j,k) = ( -   2.d0 * u(i,j-5,k)     &      
     &                   +  15.d0 * u(i,j-4,k)     &     
     &                   -  50.d0 * u(i,j-3,k)     &     
     &                   + 100.d0 * u(i,j-2,k)     &     
     &                   - 150.d0 * u(i,j-1,k)     &     
     &                   +  77.d0 * u(i,j  ,k)     &     
     &                   +  10.d0 * u(i,j+1,k)     &     
     &                 ) * idy_by_60 

          j=ny 
          Dyu(i,j,k) = (    10.d0 * u(i,j-6,k)     &      
     &                   -  72.d0 * u(i,j-5,k)     &     
     &                   + 225.d0 * u(i,j-4,k)     &     
     &                   - 400.d0 * u(i,j-3,k)     &     
     &                   + 450.d0 * u(i,j-2,k)     &     
     &                   - 360.d0 * u(i,j-1,k)     &     
     &                   + 147.d0 * u(i,j  ,k)     &     
     &                 ) * idy_by_60 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv666adv_z  calculates the first derivative with respect to 
      !                z using an "advective" sixth order finite 
      !                difference operator at all points on the 
      !                interval [5,nz-4].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^z \partial_z in the equations.  
      !                Further, it will depend on the sign of beta^z.  
      !                Namely, if beta^z is positive, we will use a 
      !                sixth order upwind operator while if beta^z 
      !                is negative, we use a sixth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At k=4, if beta^z is positive, we continue to  
      !                use the sixth order upwind operator but if  
      !                beta^z is negative, we use the sixth order 
      !                centered operator.   
      !                  
      !                We reverse this at k=nz-3 such that if beta^z is  
      !                positive, we use a sixth order centered operator 
      !                while if beta^z is negative we continue to use  
      !                the sixth order downwind operator.   
      !                 
      !                At k=3, we want a sixth order operator so it 
      !                does not matter what the sign of beta^z is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of  
      !                using the downwind operator for the case of 
      !                negative beta^z, but I do not see any way around 
      !                this if we want to have a sixth order operator.  
      !                 
      !                Likewise, at k=nz-2, the requirement of a sixth
      !                order operator renders the sign of beta^z 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^z is positive.  
      !                 
      !                At k=2 and k=nz-1, we use differently shifted 
      !                sixth order operators irrespective of beta^zs 
      !                sign.   
      !                 
      !                At k=1 and k=nz, we use the totally shifted sixth
      !                order operators irrespective of beta^zs sign.
      !
      !----------------------------------------------------------------      
      subroutine deriv666adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz), betaz(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12, idz_by_60 

        idz = 1.0d0/dz        
        idz_by_2 = 0.5d0 * idz
        idz_by_12 = idz / 12.d0
        idz_by_60 = idz / 60.d0


        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( - 147.d0 * u(i,j,k  )     &      
     &                   + 360.d0 * u(i,j,k+1)     &     
     &                   - 450.d0 * u(i,j,k+2)     &     
     &                   + 400.d0 * u(i,j,k+3)     &     
     &                   - 225.d0 * u(i,j,k+4)     &     
     &                   +  72.d0 * u(i,j,k+5)     &     
     &                   -  10.d0 * u(i,j,k+6)     &     
     &                 ) * idz_by_60 

          k=2
          Dzu(i,j,k) = ( -  10.d0 * u(i,j,k-1)     &      
     &                   -  77.d0 * u(i,j,k  )     &     
     &                   + 150.d0 * u(i,j,k+1)     &     
     &                   - 100.d0 * u(i,j,k+2)     &     
     &                   +  50.d0 * u(i,j,k+3)     &     
     &                   -  15.d0 * u(i,j,k+4)     &     
     &                   +   2.d0 * u(i,j,k+5)     &     
     &                 ) * idz_by_60 

          k=3
          Dzu(i,j,k) = (    2.d0 * u(i,j,k-2)     &      
     &                   - 24.d0 * u(i,j,k-1)     &     
     &                   - 35.d0 * u(i,j,k  )     &     
     &                   + 80.d0 * u(i,j,k+1)     &     
     &                   - 30.d0 * u(i,j,k+2)     &     
     &                   +  8.d0 * u(i,j,k+3)     &     
     &                   -         u(i,j,k+4)     &     
     &                 ) * idz_by_60  

          k=4
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (    2.d0 * u(i,j,k-2)     &      
     &                     - 24.d0 * u(i,j,k-1)     &     
     &                     - 35.d0 * u(i,j,k  )     &     
     &                     + 80.d0 * u(i,j,k+1)     &     
     &                     - 30.d0 * u(i,j,k+2)     &     
     &                     +  8.d0 * u(i,j,k+3)     &     
     &                     -         u(i,j,k+4)     &     
     &                   ) * idz_by_60 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  9.d0 * u(i,j,k-2)     &     
     &                     - 45.d0 * u(i,j,k-1)     &     
     &                     + 45.d0 * u(i,j,k+1)     &     
     &                     -  9.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_60  
          endif
          
          do k = 5, nz-4  ! see gr-qc/0505055v2.pdf 
            if ( betaz(i,j,k) .ge. 0.d0 ) then
              Dzu(i,j,k) = (    2.d0 * u(i,j,k-2)     &      
     &                       - 24.d0 * u(i,j,k-1)     &     
     &                       - 35.d0 * u(i,j,k  )     &     
     &                       + 80.d0 * u(i,j,k+1)     &     
     &                       - 30.d0 * u(i,j,k+2)     &     
     &                       +  8.d0 * u(i,j,k+3)     &     
     &                       -         u(i,j,k+4)     &     
     &                     ) * idz_by_60 
            else 
              Dzu(i,j,k) = (           u(i,j,k-4)     &      
     &                       -  8.d0 * u(i,j,k-3)     &     
     &                       + 30.d0 * u(i,j,k-2)     &     
     &                       - 80.d0 * u(i,j,k-1)     &     
     &                       + 35.d0 * u(i,j,k  )     &     
     &                       + 24.d0 * u(i,j,k+1)     &     
     &                       -  2.d0 * u(i,j,k+2)     &     
     &                     ) * idz_by_60 
            endif 
          end do

          k=nz-3 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  9.d0 * u(i,j,k-2)     &     
     &                     - 45.d0 * u(i,j,k-1)     &     
     &                     + 45.d0 * u(i,j,k+1)     &     
     &                     -  9.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_60  
          else 
            Dzu(i,j,k) = (           u(i,j,k-4)     &      
     &                     -  8.d0 * u(i,j,k-3)     &     
     &                     + 30.d0 * u(i,j,k-2)     &     
     &                     - 80.d0 * u(i,j,k-1)     &     
     &                     + 35.d0 * u(i,j,k  )     &     
     &                     + 24.d0 * u(i,j,k+1)     &     
     &                     -  2.d0 * u(i,j,k+2)     &     
     &                   ) * idz_by_60 
          endif 

          k=nz-2 
          Dzu(i,j,k) = (           u(i,j,k-4)     &      
     &                   -  8.d0 * u(i,j,k-3)     &     
     &                   + 30.d0 * u(i,j,k-2)     &     
     &                   - 80.d0 * u(i,j,k-1)     &     
     &                   + 35.d0 * u(i,j,k  )     &     
     &                   + 24.d0 * u(i,j,k+1)     &     
     &                   -  2.d0 * u(i,j,k+2)     &     
     &                 ) * idz_by_60 

          k=nz-1 
          Dzu(i,j,k) = ( -   2.d0 * u(i,j,k-5)     &      
     &                   +  15.d0 * u(i,j,k-4)     &     
     &                   -  50.d0 * u(i,j,k-3)     &     
     &                   + 100.d0 * u(i,j,k-2)     &     
     &                   - 150.d0 * u(i,j,k-1)     &     
     &                   +  77.d0 * u(i,j,k  )     &     
     &                   +  10.d0 * u(i,j,k+1)     &     
     &                 ) * idz_by_60 

          k=nz 
          Dzu(i,j,k) = (    10.d0 * u(i,j,k-6)     &      
     &                   -  72.d0 * u(i,j,k-5)     &     
     &                   + 225.d0 * u(i,j,k-4)     &     
     &                   - 400.d0 * u(i,j,k-3)     &     
     &                   + 450.d0 * u(i,j,k-2)     &     
     &                   - 360.d0 * u(i,j,k-1)     &     
     &                   + 147.d0 * u(i,j,k  )     &     
     &                 ) * idz_by_60 
                       
        end do
        end do
       
        return
      end subroutine 
      


      !----------------------------------------------------------------
      !
      ! deriv642_xx  calculates the second derivative with respect to  
      !              x using a centered sixth order finite difference 
      !              operator on the interval [4,nx-3].  At other 
      !              points, we attempt to use centered difference 
      !              operators.  
      !   
      !              At i=3 and i=nx-2, we use centered fourth order 
      !              operators.  
      ! 
      !              At i=2 and i=nx-1, we use centered second order 
      !              operators.  
      ! 
      !              At the boundary points (i=1 and i=nx) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at i=1, we use a forward (upwind)  
      !              operator with a four point stencil.  At i=nx, we  
      !              use a backward (downwind) operator, again with 
      !              a four point stencil. 
      !
      !----------------------------------------------------------------
      subroutine deriv642_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DxDxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idx_sqrd, idx_sqrd_by_12, idx_sqrd_by_180

        idx_sqrd = 1.0d0/(dx*dx)
        idx_sqrd_by_12 = idx_sqrd / 12.d0
        idx_sqrd_by_180 = idx_sqrd / 180.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          DxDxu(i,j,k) = (   2.d0 * u(i  ,j,k)     &     
     &                     - 5.d0 * u(i+1,j,k)     &     
     &                     + 4.d0 * u(i+2,j,k)     &     
     &                     -        u(i+3,j,k)     &     
     &                   ) * idx_sqrd
          i=2
          DxDxu(i,j,k) = (          u(i-1,j,k)     &     
     &                     - 2.d0 * u(i  ,j,k)     &     
     &                     +        u(i+1,j,k)     &     
     &                   ) * idx_sqrd
          i=3
          DxDxu(i,j,k) = ( -         u(i-2,j,k)     &     
     &                     + 16.d0 * u(i-1,j,k)     &     
     &                     - 30.d0 * u(i  ,j,k)     &     
     &                     + 16.d0 * u(i+1,j,k)     &     
     &                     -         u(i+2,j,k)     &     
     &                   ) * idx_sqrd_by_12
          do i = 4, nx-3
            DxDxu(i,j,k) = (     2.d0 * u(i-3,j,k)     &     
     &                       -  27.d0 * u(i-2,j,k)     &     
     &                       + 270.d0 * u(i-1,j,k)     &     
     &                       - 490.d0 * u(i  ,j,k)     &     
     &                       + 270.d0 * u(i+1,j,k)     &     
     &                       -  27.d0 * u(i+2,j,k)     &     
     &                       +   2.d0 * u(i+3,j,k)     &     
                           ) * idx_sqrd_by_180
          end do
          i=nx-2
          DxDxu(i,j,k) = ( -         u(i-2,j,k)     &     
     &                     + 16.d0 * u(i-1,j,k)     &     
     &                     - 30.d0 * u(i  ,j,k)     &     
     &                     + 16.d0 * u(i+1,j,k)     &     
     &                     -         u(i+2,j,k)     &     
     &                   ) * idx_sqrd_by_12
          i=nx-1
          DxDxu(i,j,k) = (          u(i-1,j,k)     &     
     &                     - 2.d0 * u(i  ,j,k)     &     
     &                     +        u(i+1,j,k)     &     
     &                   ) * idx_sqrd 
          i=nx          
          DxDxu(i,j,k) = ( -        u(i-3,j,k)     &     
     &                     + 4.d0 * u(i-2,j,k)     &     
     &                     - 5.d0 * u(i-1,j,k)     &     
     &                     + 2.d0 * u(i  ,j,k)     &     
     &                   ) * idx_sqrd 
        end do          
        end do          
     
        return        
      end subroutine
        

      !----------------------------------------------------------------
      !
      ! deriv642_yy  calculates the second derivative with respect to  
      !              y using a centered sixth order finite difference 
      !              operator on the interval [4,ny-3].  At other 
      !              points, we attempt to use centered difference 
      !              operators.  
      !   
      !              At j=3 and j=ny-2, we use centered fourth order 
      !              operators.  
      ! 
      !              At j=2 and j=ny-1, we use centered second order 
      !              operators.  
      ! 
      !              At the boundary points (j=1 and j=ny) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at j=1, we use a forward (upwind)  
      !              operator with a four point stencil.  At j=ny, we  
      !              use a backward (downwind) operator, again with 
      !              a four point stencil. 
      !
      !----------------------------------------------------------------
      subroutine deriv642_yy(DyDyu, u, dy, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDyu(nx,ny,nz), u(nx,ny,nz) 
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 
        
        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idy_sqrd, idy_sqrd_by_12, idy_sqrd_by_180
        
        idy_sqrd = 1.0d0/(dy*dy)
        idy_sqrd_by_12 = idy_sqrd / 12.d0   
        idy_sqrd_by_180 = idy_sqrd / 180.d0
        
        do k = 1, nz
        do i = 1, nx
          j=1
          DyDyu(i,j,k) = (   2.d0 * u(i,j  ,k)     &     
     &                     - 5.d0 * u(i,j+1,k)     &     
     &                     + 4.d0 * u(i,j+2,k)     &     
     &                     -        u(i,j+3,k)     &     
     &                   ) * idy_sqrd
          j=2
          DyDyu(i,j,k) = (          u(i,j-1,k)     &     
     &                     - 2.d0 * u(i,j  ,k)     &      
     &                     +        u(i,j+1,k)     &     
     &                   ) * idy_sqrd
          j=3           
          DyDyu(i,j,k) = ( -         u(i,j-2,k)     &     
     &                     + 16.d0 * u(i,j-1,k)     &     
     &                     - 30.d0 * u(i,j  ,k)     &     
     &                     + 16.d0 * u(i,j+1,k)     &     
     &                     -         u(i,j+2,k)     &     
     &                   ) * idy_sqrd_by_12 
          do j = 4, ny-3
            DyDyu(i,j,k) = (     2.d0 * u(i,j-3,k)     &     
     &                       -  27.d0 * u(i,j-2,k)     &     
     &                       + 270.d0 * u(i,j-1,k)     &     
     &                       - 490.d0 * u(i,j  ,k)     &     
     &                       + 270.d0 * u(i,j+1,k)     &     
     &                       -  27.d0 * u(i,j+2,k)     &     
     &                       +   2.d0 * u(i,j+3,k)     &     
     &                     ) * idy_sqrd_by_180   
          end do          
          j=ny-2          
          DyDyu(i,j,k) = ( -         u(i,j-2,k)     &     
     &                     + 16.d0 * u(i,j-1,k)     &     
     &                     - 30.d0 * u(i,j  ,k)     &     
     &                     + 16.d0 * u(i,j+1,k)     &     
     &                     -         u(i,j+2,k)     &     
     &                   ) * idy_sqrd_by_12 
          j=ny-1        
          DyDyu(i,j,k) = (          u(i,j-1,k)     &      
     &                     - 2.d0 * u(i,j  ,k)     &     
     &                     +        u(i,j+1,k)     &     
     &                   ) * idy_sqrd 
          j=ny          
          DyDyu(i,j,k) = ( -        u(i,j-3,k)     &     
     &                     + 4.d0 * u(i,j-2,k)     &     
     &                     - 5.d0 * u(i,j-1,k)     &     
     &                     + 2.d0 * u(i,j  ,k)     &     
     &                   ) * idy_sqrd 
        end do
        end do

        return
      end subroutine


      !----------------------------------------------------------------
      !
      ! deriv642_zz  calculates the second derivative with respect to  
      !              z using a centered sixth order finite difference 
      !              operator on the interval [4,nz-3].  At other 
      !              points, we attempt to use centered difference 
      !              operators.  
      !   
      !              At k=3 and k=nz-2, we use centered fourth order 
      !              operators.  
      ! 
      !              At k=2 and k=nz-1, we use centered second order 
      !              operators.  
      ! 
      !              At the boundary points (k=1 and k=nz) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at k=1, we use a forward (upwind)  
      !              operator with a four point stencil.  At k=nz, we  
      !              use a backward (downwind) operator, again with 
      !              a four point stencil. 
      !
      !----------------------------------------------------------------
      subroutine deriv642_zz(DzDzu, u, dz, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DzDzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idz_sqrd, idz_sqrd_by_12, idz_sqrd_by_180

        idz_sqrd = 1.0d0/(dz*dz)
        idz_sqrd_by_12 = idz_sqrd / 12.d0
        idz_sqrd_by_180 = idz_sqrd / 180.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          DzDzu(i,j,k) = (   2.d0 * u(i,j,k  )     &     
     &                     - 5.d0 * u(i,j,k+1)     &     
     &                     + 4.d0 * u(i,j,k+2)     &     
     &                     -        u(i,j,k+3)     &     
     &                   ) * idz_sqrd
          k=2
          DzDzu(i,j,k) = (          u(i,j,k-1)     &     
     &                     - 2.d0 * u(i,j,k  )     &     
     &                     +        u(i,j,k+1)     &     
     &                   ) * idz_sqrd
          k=3
          DzDzu(i,j,k) = ( -         u(i,j,k-2)     &     
     &                     + 16.d0 * u(i,j,k-1)     &     
     &                     - 30.d0 * u(i,j,k  )     &     
     &                     + 16.d0 * u(i,j,k+1)     &     
     &                     -         u(i,j,k+2)     &     
     &                   ) * idz_sqrd_by_12
          do k = 4, nz-3
            DzDzu(i,j,k) = (     2.d0 * u(i,j,k-3)     &     
     &                       -  27.d0 * u(i,j,k-2)     &     
     &                       + 270.d0 * u(i,j,k-1)     &     
     &                       - 490.d0 * u(i,j,k  )     &     
     &                       + 270.d0 * u(i,j,k+1)     &     
     &                       -  27.d0 * u(i,j,k+2)     &     
     &                       +   2.d0 * u(i,j,k+3)     &     
     &                     ) * idz_sqrd_by_180
          end do
          k=nz-2
          DzDzu(i,j,k) = ( -         u(i,j,k-2)     &     
     &                     + 16.d0 * u(i,j,k-1)     &     
     &                     - 30.d0 * u(i,j,k  )     &     
     &                     + 16.d0 * u(i,j,k+1)     &     
     &                     -         u(i,j,k+2)     &     
     &                   ) * idz_sqrd_by_12
          k=nz-1
          DzDzu(i,j,k) = (          u(i,j,k-1)     &     
     &                     - 2.d0 * u(i,j,k  )     &     
     &                     +        u(i,j,k+1)     &     
     &                   ) * idz_sqrd
          k=nz
          DzDzu(i,j,k) = ( -        u(i,j,k-3)     &     
     &                     + 4.d0 * u(i,j,k-2)     &     
     &                     - 5.d0 * u(i,j,k-1)     &     
     &                     + 2.d0 * u(i,j,k  )     &     
     &                   ) * idz_sqrd 
        end do
        end do
          
        return
      end subroutine

        
      !----------------------------------------------------------------
      !
      ! deriv666_xx  calculates the second derivative with respect to  
      !              x using a centered sixth order finite difference 
      !              operator at all points on the interval [4,nx-3].  
      ! 
      !              At i=3 and i=nx-2, we use once-shifted sixth order 
      !              operators.  This has a eight point stencil.  
      ! 
      !              At i=2 and i=nx-1, we use twice-shifted sixth order
      !              operators.  This has a eight point stencil.  
      ! 
      !              At i=1 and i=nx, we use totally shifted sixth order
      !              operators.  This has a eight point stencil.  
      ! 
      !----------------------------------------------------------------
      subroutine deriv666_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DxDxu(nx,ny,nz), u(nx,ny,nz) 
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 
        
        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idx_sqrd, idx_sqrd_by_180
        
        idx_sqrd = 1.0d0/(dx*dx)
        idx_sqrd_by_180 = idx_sqrd / 180.d0 

        do k = 1, nz
        do j = 1, ny
          i=1
          DxDxu(i,j,k) = (    938.d0 * u(i  ,j,k)     &     
     &                     - 4014.d0 * u(i+1,j,k)     &     
     &                     + 7911.d0 * u(i+2,j,k)     &     
     &                     - 9490.d0 * u(i+3,j,k)     &     
     &                     + 7380.d0 * u(i+4,j,k)     &      
     &                     - 3618.d0 * u(i+5,j,k)     &     
     &                     + 1019.d0 * u(i+6,j,k)     &     
     &                     -  126.d0 * u(i+7,j,k)     &     
     &                   ) * idx_sqrd_by_180  
          i=2         
          DxDxu(i,j,k) = (    126.d0 * u(i-1,j,k)     &     
     &                      -  70.d0 * u(i  ,j,k)     &     
     &                      - 486.d0 * u(i+1,j,k)     &     
     &                      + 855.d0 * u(i+2,j,k)     &     
     &                      - 670.d0 * u(i+3,j,k)     &     
     &                      + 324.d0 * u(i+4,j,k)     &     
     &                      -  90.d0 * u(i+5,j,k)     &     
     &                      +  11.d0 * u(i+6,j,k)     &     
     &                   ) * idx_sqrd_by_180 
          i=3             
          DxDxu(i,j,k) = (  -  11.d0 * u(i-2,j,k)     &      
     &                      + 214.d0 * u(i-1,j,k)     &     
     &                      - 378.d0 * u(i  ,j,k)     &     
     &                      + 130.d0 * u(i+1,j,k)     &      
     &                      +  85.d0 * u(i+2,j,k)     &      
     &                      -  54.d0 * u(i+3,j,k)     &     
     &                      +  16.d0 * u(i+4,j,k)     &     
     &                      -   2.d0 * u(i+5,j,k)     &     
     &                   ) * idx_sqrd_by_180 
          do i = 4, nx-3
            DxDxu(i,j,k) = (     2.d0 * u(i-3,j,k)     &     
     &                       -  27.d0 * u(i-2,j,k)     &     
     &                       + 270.d0 * u(i-1,j,k)     &     
     &                       - 490.d0 * u(i  ,j,k)     &     
     &                       + 270.d0 * u(i+1,j,k)     &     
     &                       -  27.d0 * u(i+2,j,k)     &     
     &                       +   2.d0 * u(i+3,j,k)     &     
     &                     ) * idx_sqrd_by_180
          end do      
          i=nx-2
          DxDxu(i,j,k) = (  -   2.d0 * u(i-5,j,k)     &      
     &                      +  16.d0 * u(i-4,j,k)     &     
     &                      -  54.d0 * u(i-3,j,k)     &     
     &                      +  85.d0 * u(i-2,j,k)     &     
     &                      + 130.d0 * u(i-1,j,k)     &     
     &                      - 378.d0 * u(i  ,j,k)     &     
     &                      + 214.d0 * u(i+1,j,k)     &     
     &                      -  11.d0 * u(i+2,j,k)     &     
     &                   ) * idx_sqrd_by_180
          i=nx-1
          DxDxu(i,j,k) = (     11.d0 * u(i-6,j,k)     &     
     &                      -  90.d0 * u(i-5,j,k)     &     
     &                      + 324.d0 * u(i-4,j,k)     &     
     &                      - 670.d0 * u(i-3,j,k)     &     
     &                      + 855.d0 * u(i-2,j,k)     &     
     &                      - 486.d0 * u(i-1,j,k)     &     
     &                      -  70.d0 * u(i  ,j,k)     &     
     &                      + 126.d0 * u(i+1,j,k)     &     
     &                   ) * idx_sqrd_by_180
          i=nx
          DxDxu(i,j,k) = ( -  126.d0 * u(i-7,j,k)     &     
     &                     + 1019.d0 * u(i-6,j,k)     &     
     &                     - 3618.d0 * u(i-5,j,k)     &     
     &                     + 7380.d0 * u(i-4,j,k)     &     
     &                     - 9490.d0 * u(i-3,j,k)     &     
     &                     + 7911.d0 * u(i-2,j,k)     &     
     &                     - 4014.d0 * u(i-1,j,k)     &     
     &                     +  938.d0 * u(i  ,j,k)     &     
     &                      ) * idx_sqrd_by_180
        end do
        end do

        return
      end subroutine


      !----------------------------------------------------------------
      !
      ! deriv666_yy  calculates the second derivative with respect to  
      !              y using a centered sixth order finite difference 
      !              operator at all points on the interval [4,ny-3].  
      ! 
      !              At j=3 and j=ny-2, we use once-shifted sixth order 
      !              operators.  This has a eight point stencil.  
      ! 
      !              At j=2 and j=ny-1, we use twice-shifted sixth order
      !              operators.  This has a eight point stencil.  
      ! 
      !              At j=1 and j=ny, we use totally shifted sixth order
      !              operators.  This has a eight point stencil.  
      !
      !----------------------------------------------------------------      
      subroutine deriv666_yy(DyDyu, u, dy, nx, ny, nz, dd_type)
        implicit   none
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idy_sqrd, idy_sqrd_by_180

        idy_sqrd = 1.0d0/(dy*dy)
        idy_sqrd_by_180 = idy_sqrd / 180.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          DyDyu(i,j,k) = (    938.d0 * u(i,j  ,k)     &     
     &                     - 4014.d0 * u(i,j+1,k)     &     
     &                     + 7911.d0 * u(i,j+2,k)     &     
     &                     - 9490.d0 * u(i,j+3,k)     &     
     &                     + 7380.d0 * u(i,j+4,k)     &     
     &                     - 3618.d0 * u(i,j+5,k)     &     
     &                     + 1019.d0 * u(i,j+6,k)     &     
     &                     -  126.d0 * u(i,j+7,k)     &     
     &                   ) * idy_sqrd_by_180
          j=2
          DyDyu(i,j,k) = (    126.d0 * u(i,j-1,k)     &     
     &                      -  70.d0 * u(i,j  ,k)     &     
     &                      - 486.d0 * u(i,j+1,k)     &     
     &                      + 855.d0 * u(i,j+2,k)     &     
     &                      - 670.d0 * u(i,j+3,k)     &     
     &                      + 324.d0 * u(i,j+4,k)     &     
     &                      -  90.d0 * u(i,j+5,k)     &     
     &                      +  11.d0 * u(i,j+6,k)     &     
     &                   ) * idy_sqrd_by_180 
          j=3            
          DyDyu(i,j,k) = (  -  11.d0 * u(i,j-2,k)     &     
     &                      + 214.d0 * u(i,j-1,k)     &     
     &                      - 378.d0 * u(i,j  ,k)     &     
     &                      + 130.d0 * u(i,j+1,k)     &     
     &                      +  85.d0 * u(i,j+2,k)     &     
     &                      -  54.d0 * u(i,j+3,k)     &     
     &                      +  16.d0 * u(i,j+4,k)     &     
     &                      -   2.d0 * u(i,j+5,k)     &     
     &                   ) * idy_sqrd_by_180 
          do j = 4, ny-3 
            DyDyu(i,j,k) = (     2.d0 * u(i,j-3,k)     &     
     &                       -  27.d0 * u(i,j-2,k)     &     
     &                       + 270.d0 * u(i,j-1,k)     &     
     &                       - 490.d0 * u(i,j  ,k)     &     
     &                       + 270.d0 * u(i,j+1,k)     &     
     &                       -  27.d0 * u(i,j+2,k)     &     
     &                       +   2.d0 * u(i,j+3,k)     &     
     &                     ) * idy_sqrd_by_180
          end do        
          j=ny-2        
          DyDyu(i,j,k) = (  -   2.d0 * u(i,j-5,k)     &     
     &                      +  16.d0 * u(i,j-4,k)     &     
     &                      -  54.d0 * u(i,j-3,k)     &     
     &                      +  85.d0 * u(i,j-2,k)     &     
     &                      + 130.d0 * u(i,j-1,k)     &     
     &                      - 378.d0 * u(i,j  ,k)     &     
     &                      + 214.d0 * u(i,j+1,k)     &     
     &                      -  11.d0 * u(i,j+2,k)     &     
     &                   ) * idy_sqrd_by_180
          j=ny-1
          DyDyu(i,j,k) = (     11.d0 * u(i,j-6,k)     &     
     &                      -  90.d0 * u(i,j-5,k)     &     
     &                      + 324.d0 * u(i,j-4,k)     &     
     &                      - 670.d0 * u(i,j-3,k)     &      
     &                      + 855.d0 * u(i,j-2,k)     &     
     &                      - 486.d0 * u(i,j-1,k)     &     
     &                      -  70.d0 * u(i,j  ,k)     &     
     &                      + 126.d0 * u(i,j+1,k)     &     
     &                   ) * idy_sqrd_by_180
          j=ny 
          DyDyu(i,j,k) = ( -  126.d0 * u(i,j-7,k)     &     
     &                     + 1019.d0 * u(i,j-6,k)     &     
     &                     - 3618.d0 * u(i,j-5,k)     &     
     &                     + 7380.d0 * u(i,j-4,k)     &     
     &                     - 9490.d0 * u(i,j-3,k)     &     
     &                     + 7911.d0 * u(i,j-2,k)     &     
     &                     - 4014.d0 * u(i,j-1,k)     &     
     &                     +  938.d0 * u(i,j  ,k)     &     
     &                   ) * idy_sqrd_by_180
        end do
        end do          
                        
        return          
      end subroutine    
     
      !----------------------------------------------------------------
      !                 
      ! deriv666_zz  calculates the second derivative with respect to  
      !              z using a centered sixth order finite difference 
      !              operator at all points on the interval [4,nz-3].  
      ! 
      !              At k=3 and k=nz-2, we use once-shifted sixth order 
      !              operators.  This has a eight point stencil.  
      ! 
      !              At k=2 and k=nz-1, we use twice-shifted sixth order
      !              operators.  This has a eight point stencil.  
      ! 
      !              At k=1 and k=nz, we use totally shifted sixth order
      !              operators.  This has a eight point stencil.  
      !   
      !----------------------------------------------------------------      
      subroutine deriv666_zz(DzDzu, u, dz, nx, ny, nz, dd_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL DzDzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k
        CCTK_REAL idz_sqrd, idz_sqrd_by_180

        idz_sqrd = 1.0d0/(dz*dz)
        idz_sqrd_by_180 = idz_sqrd / 180.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          DzDzu(i,j,k) = (    938.d0 * u(i,j,k  )     &     
     &                     - 4014.d0 * u(i,j,k+1)     &     
     &                     + 7911.d0 * u(i,j,k+2)     &     
     &                     - 9490.d0 * u(i,j,k+3)     &     
     &                     + 7380.d0 * u(i,j,k+4)     &     
     &                     - 3618.d0 * u(i,j,k+5)     &     
     &                     + 1019.d0 * u(i,j,k+6)     &     
     &                     -  126.d0 * u(i,j,k+7)     &     
     &                   ) * idz_sqrd_by_180
          k=2
          DzDzu(i,j,k) = (    126.d0 * u(i,j,k-1)     &     
     &                      -  70.d0 * u(i,j,k  )     &     
     &                      - 486.d0 * u(i,j,k+1)     &     
     &                      + 855.d0 * u(i,j,k+2)     &     
     &                      - 670.d0 * u(i,j,k+3)     &     
     &                      + 324.d0 * u(i,j,k+4)     &     
     &                      -  90.d0 * u(i,j,k+5)     &     
     &                      +  11.d0 * u(i,j,k+6)     &     
     &                   ) * idz_sqrd_by_180
          k=3
          DzDzu(i,j,k) = (  -  11.d0 * u(i,j,k-2)     &     
     &                      + 214.d0 * u(i,j,k-1)     &     
     &                      - 378.d0 * u(i,j,k  )     &     
     &                      + 130.d0 * u(i,j,k+1)     &     
     &                      +  85.d0 * u(i,j,k+2)     &     
     &                      -  54.d0 * u(i,j,k+3)     &     
     &                      +  16.d0 * u(i,j,k+4)     &     
     &                      -   2.d0 * u(i,j,k+5)     &     
     &                   ) * idz_sqrd_by_180
          do k = 4, nz-3
            DzDzu(i,j,k) = (     2.d0 * u(i,j,k-3)     &     
     &                       -  27.d0 * u(i,j,k-2)     &     
     &                       + 270.d0 * u(i,j,k-1)     &     
     &                       - 490.d0 * u(i,j,k  )     &     
     &                       + 270.d0 * u(i,j,k+1)     &     
     &                       -  27.d0 * u(i,j,k+2)     &     
     &                       +   2.d0 * u(i,j,k+3)     &     
     &                     ) * idz_sqrd_by_180
          end do
          k=nz-2
          DzDzu(i,j,k) = (  -   2.d0 * u(i,j,k-5)     &     
     &                      +  16.d0 * u(i,j,k-4)     &     
     &                      -  54.d0 * u(i,j,k-3)     &     
     &                      +  85.d0 * u(i,j,k-2)     &     
     &                      + 130.d0 * u(i,j,k-1)     &     
     &                      - 378.d0 * u(i,j,k  )     &     
     &                      + 214.d0 * u(i,j,k+1)     &     
     &                      -  11.d0 * u(i,j,k+2)     &     
     &                   ) * idz_sqrd_by_180
          k=nz-1
          DzDzu(i,j,k) = (     11.d0 * u(i,j,k-6)     &     
     &                      -  90.d0 * u(i,j,k-5)     &     
     &                      + 324.d0 * u(i,j,k-4)     &     
     &                      - 670.d0 * u(i,j,k-3)     &     
     &                      + 855.d0 * u(i,j,k-2)     &     
     &                      - 486.d0 * u(i,j,k-1)     &     
     &                      -  70.d0 * u(i,j,k  )     &     
     &                      + 126.d0 * u(i,j,k+1)     &     
     &                   ) * idz_sqrd_by_180
          k=nz 
          DzDzu(i,j,k) = ( -  126.d0 * u(i,j,k-7)     &     
     &                     + 1019.d0 * u(i,j,k-6)     &     
     &                     - 3618.d0 * u(i,j,k-5)     &     
     &                     + 7380.d0 * u(i,j,k-4)     &     
     &                     - 9490.d0 * u(i,j,k-3)     &     
     &                     + 7911.d0 * u(i,j,k-2)     &     
     &                     - 4014.d0 * u(i,j,k-1)     &     
     &                     +  938.d0 * u(i,j,k  )     &     
     &                   ) * idz_sqrd_by_180
        end do
        end do          
                        
        return          
      end subroutine    
     
     








!======================================================================! 
!                                                                      ! 
!                  EIGHTH ORDER DERIVATIVE OPERATORS                   ! 
!                                                                      ! 
!======================================================================! 
      !----------------------------------------------------------------
      !
      ! deriv8642_x  calculates the first derivative with respect to 
      !              x using a centered eighth order finite difference 
      !              operator on the interval [5,nx-4].  At the points 
      !              i=4 and i=nx-3, a centered sixth order operator
      !              is used.  At i=3 and i=nx-2, we use a centered 
      !              fourth order operator.  At i=2 and i=nx-1, we use 
      !              a centered second order operator.  At the boundary 
      !              points, i=1 and i=nx, we spoil the spirit of this 
      !              effort to use only centered difference operators       
      !              and use instead the appropriate second order upwind  
      !              or downwind operator. 
      !
      !----------------------------------------------------------------      
      subroutine deriv8642_x(Dxu, u, dx, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60, idx_by_2520  

        idx = 1.0d0/dx        
        idx_by_2    = idx / 2.d0
        idx_by_12   = idx / 12.d0
        idx_by_60   = idx / 60.d0
        idx_by_2520 = idx / 2520.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  3.d0 * u(i  ,j,k)     &      
     &                   +  4.d0 * u(i+1,j,k)     &     
     &                   -         u(i+2,j,k)     &     
     &                 ) * idx_by_2
          i=2
          Dxu(i,j,k) = ( -         u(i-1,j,k)     &      
     &                   +         u(i+1,j,k)     &     
     &                 ) * idx_by_2
          i=3
          Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                   -  8.d0 * u(i-1,j,k)     &     
     &                   +  8.d0 * u(i+1,j,k)     &      
     &                   -         u(i+2,j,k)     &      
     &                 ) * idx_by_12 
          i=4
          Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                   +  9.d0 * u(i-2,j,k)     &     
     &                   - 45.d0 * u(i-1,j,k)     &     
     &                   + 45.d0 * u(i+1,j,k)     &     
     &                   -  9.d0 * u(i+2,j,k)     &     
     &                   +         u(i+3,j,k)     &     
     &                 ) * idx_by_60
          do i = 5, nx-4
            Dxu(i,j,k) = (      9.d0 * u(i-4,j,k)     &      
     &                     -   96.d0 * u(i-3,j,k)     &     
     &                     +  504.d0 * u(i-2,j,k)     &     
     &                     - 2016.d0 * u(i-1,j,k)     &     
     &                     + 2016.d0 * u(i+1,j,k)     &     
     &                     -  504.d0 * u(i+2,j,k)     &     
     &                     +   96.d0 * u(i+3,j,k)     &     
     &                     -    9.d0 * u(i+4,j,k)     &     
     &                   ) * idx_by_2520 
          end do
          i=nx-3
          Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                   +  9.d0 * u(i-2,j,k)     &     
     &                   - 45.d0 * u(i-1,j,k)     &     
     &                   + 45.d0 * u(i+1,j,k)     &     
     &                   -  9.d0 * u(i+2,j,k)     &     
     &                   +         u(i+3,j,k)     &     
     &                 ) * idx_by_60 
          i=nx-2 
          Dxu(i,j,k) = (          u(i-2,j,k)     &      
     &                   - 8.d0 * u(i-1,j,k)     &     
     &                   + 8.d0 * u(i+1,j,k)     &      
     &                   -        u(i+2,j,k)     &      
     &                 ) * idx_by_12
          i=nx-1 
          Dxu(i,j,k) = ( -        u(i-1,j,k)     &      
     &                   +        u(i+1,j,k)     &      
     &                 ) * idx_by_2
          i=nx 
          Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                   -  4.d0 * u(i-1,j,k)     &     
     &                   +  3.d0 * u(i  ,j,k)     &     
     &                 ) * idx_by_2
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      ! deriv8642_y  calculates the first derivative with respect to 
      !              y using a centered eighth order finite difference 
      !              operator on the interval [5,ny-4].  At the points 
      !              j=4 and j=ny-3, a centered sixth order operator
      !              is used.  At j=3 and j=ny-2, we use a centered 
      !              fourth order operator.  At j=2 and j=ny-1, we use 
      !              a centered second order operator.  At the boundary 
      !              points, j=1 and j=ny, we spoil the spirit of this 
      !              effort to use only centered difference operators       
      !              and use instead the appropriate second order upwind  
      !              or downwind operator. 
      !
      !----------------------------------------------------------------      
      subroutine deriv8642_y(Dyu, u, dy, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12, idy_by_60, idy_by_2520  

        idy = 1.0d0/dy        
        idy_by_2    = idy / 2.d0
        idy_by_12   = idy / 12.d0
        idy_by_60   = idy / 60.d0
        idy_by_2520 = idy / 2520.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  3.d0 * u(i,j  ,k)     &      
     &                   +  4.d0 * u(i,j+1,k)     &     
     &                   -         u(i,j+2,k)     &     
     &                 ) * idy_by_2
          j=2
          Dyu(i,j,k) = ( -         u(i,j-1,k)     &      
     &                   +         u(i,j+1,k)     &     
     &                 ) * idy_by_2
          j=3
          Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                   -  8.d0 * u(i,j-1,k)     &     
     &                   +  8.d0 * u(i,j+1,k)     &      
     &                   -         u(i,j+2,k)     &      
     &                 ) * idy_by_12 
          j=4
          Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                   +  9.d0 * u(i,j-2,k)     &     
     &                   - 45.d0 * u(i,j-1,k)     &     
     &                   + 45.d0 * u(i,j+1,k)     &     
     &                   -  9.d0 * u(i,j+2,k)     &     
     &                   +         u(i,j+3,k)     &     
     &                 ) * idy_by_60
          do j = 5, ny-4
            Dyu(i,j,k) = (      9.d0 * u(i,j-4,k)     &      
     &                     -   96.d0 * u(i,j-3,k)     &     
     &                     +  504.d0 * u(i,j-2,k)     &     
     &                     - 2016.d0 * u(i,j-1,k)     &     
     &                     + 2016.d0 * u(i,j+1,k)     &     
     &                     -  504.d0 * u(i,j+2,k)     &     
     &                     +   96.d0 * u(i,j+3,k)     &     
     &                     -    9.d0 * u(i,j+4,k)     &     
     &                   ) * idy_by_2520 
          end do
          j=ny-3
          Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                   +  9.d0 * u(i,j-2,k)     &     
     &                   - 45.d0 * u(i,j-1,k)     &     
     &                   + 45.d0 * u(i,j+1,k)     &     
     &                   -  9.d0 * u(i,j+2,k)     &     
     &                   +         u(i,j+3,k)     &     
     &                 ) * idy_by_60 
          j=ny-2 
          Dyu(i,j,k) = (          u(i,j-2,k)     &      
     &                   - 8.d0 * u(i,j-1,k)     &     
     &                   + 8.d0 * u(i,j+1,k)     &      
     &                   -        u(i,j+2,k)     &      
     &                 ) * idy_by_12
          j=ny-1 
          Dyu(i,j,k) = ( -        u(i,j-1,k)     &      
     &                   +        u(i,j+1,k)     &      
     &                 ) * idy_by_2
          j=ny 
          Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                   -  4.d0 * u(i,j-1,k)     &     
     &                   +  3.d0 * u(i,j  ,k)     &     
     &                 ) * idy_by_2
        end do
        end do
       
        return
      end subroutine 



      !----------------------------------------------------------------
      !
      ! deriv8642_z  calculates the first derivative with respect to 
      !              z using a centered eighth order finite difference 
      !              operator on the interval [5,nz-4].  At the points 
      !              k=4 and k=nz-3, a centered sixth order operator
      !              is used.  At k=3 and k=nz-2, we use a centered 
      !              fourth order operator.  At k=2 and k=nz-1, we use 
      !              a centered second order operator.  At the boundary 
      !              points, k=1 and k=nz, we spoil the spirit of this 
      !              effort to use only centered difference operators       
      !              and use instead the appropriate second order upwind  
      !              or downwind operator. 
      !
      !----------------------------------------------------------------      
      subroutine deriv8642_z(Dzu, u, dz, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12, idz_by_60, idz_by_2520 

        idz = 1.0d0/dz        
        idz_by_2    = idz / 2.d0
        idz_by_12   = idz / 12.d0
        idz_by_60   = idz / 60.d0
        idz_by_2520 = idz / 2520.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( -  3.d0 * u(i,j,k  )     &      
     &                   +  4.d0 * u(i,j,k+1)     &     
     &                   -         u(i,j,k+2)     &     
     &                 ) * idz_by_2
          k=2
          Dzu(i,j,k) = ( -         u(i,j,k-1)     &      
     &                   +         u(i,j,k+1)     &     
     &                 ) * idz_by_2
          k=3
          Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                   -  8.d0 * u(i,j,k-1)     &     
     &                   +  8.d0 * u(i,j,k+1)     &      
     &                   -         u(i,j,k+2)     &      
     &                 ) * idz_by_12 
          k=4
          Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                   +  9.d0 * u(i,j,k-2)     &     
     &                   - 45.d0 * u(i,j,k-1)     &     
     &                   + 45.d0 * u(i,j,k+1)     &     
     &                   -  9.d0 * u(i,j,k+2)     &     
     &                   +         u(i,j,k+3)     &     
     &                 ) * idz_by_60
          do k = 5, nz-4
            Dzu(i,j,k) = (      9.d0 * u(i,j,k-4)     &      
     &                     -   96.d0 * u(i,j,k-3)     &     
     &                     +  504.d0 * u(i,j,k-2)     &     
     &                     - 2016.d0 * u(i,j,k-1)     &     
     &                     + 2016.d0 * u(i,j,k+1)     &     
     &                     -  504.d0 * u(i,j,k+2)     &     
     &                     +   96.d0 * u(i,j,k+3)     &     
     &                     -    9.d0 * u(i,j,k+4)     &     
     &                   ) * idz_by_2520 
          end do
          k=nz-3
          Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                   +  9.d0 * u(i,j,k-2)     &     
     &                   - 45.d0 * u(i,j,k-1)     &     
     &                   + 45.d0 * u(i,j,k+1)     &     
     &                   -  9.d0 * u(i,j,k+2)     &     
     &                   +         u(i,j,k+3)     &     
     &                 ) * idz_by_60
          k=nz-2 
          Dzu(i,j,k) = (          u(i,j,k-2)     &      
     &                   - 8.d0 * u(i,j,k-1)     &     
     &                   + 8.d0 * u(i,j,k+1)     &      
     &                   -        u(i,j,k+2)     &      
     &                 ) * idz_by_12
          k=nz-1 
          Dzu(i,j,k) = ( -        u(i,j,k-1)     &      
     &                   +        u(i,j,k+1)     &      
     &                 ) * idz_by_2
          k=nz 
          Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                   -  4.d0 * u(i,j,k-1)     &     
     &                   +  3.d0 * u(i,j,k  )     &     
     &                 ) * idz_by_2
        end do
        end do
       
        return
      end subroutine 





      !----------------------------------------------------------------
      !
      ! deriv8888_x  calculates the first derivative with respect to 
      !              x using a centered eighth order finite difference 
      !              operator on the interval [5,nx-4].  At other points 
      !              appropriated shifted eighth order operators are used.  
      !
      !              Note that all of these eighth order operators use
      !              nine point stencils.
      !
      !----------------------------------------------------------------      
      subroutine deriv8888_x(Dxu, u, dx, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60, idx_by_2520  

        idx = 1.0d0/dx        
        idx_by_2520 = idx / 2520.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  6849.d0 * u(i  ,j,k)     &      
     &                   + 20160.d0 * u(i+1,j,k)     &     
     &                   - 35280.d0 * u(i+2,j,k)     &     
     &                   + 47040.d0 * u(i+3,j,k)     &     
     &                   - 44100.d0 * u(i+4,j,k)     &     
     &                   + 28224.d0 * u(i+5,j,k)     &     
     &                   - 11760.d0 * u(i+6,j,k)     &     
     &                   +  2880.d0 * u(i+7,j,k)     &     
     &                   -   315.d0 * u(i+8,j,k)     &     
     &                 ) * idx_by_2520
          i=2
          Dxu(i,j,k) = ( -  315.d0 * u(i-1,j,k)     &      
     &                   - 4014.d0 * u(i  ,j,k)     &     
     &                   + 8820.d0 * u(i+1,j,k)     &     
     &                   - 8820.d0 * u(i+2,j,k)     &     
     &                   + 7350.d0 * u(i+3,j,k)     &     
     &                   - 4410.d0 * u(i+4,j,k)     &     
     &                   + 1764.d0 * u(i+5,j,k)     &     
     &                   -  420.d0 * u(i+6,j,k)     &     
     &                   +   45.d0 * u(i+7,j,k)     &     
     &                 ) * idx_by_2520
          i=3
          Dxu(i,j,k) = (     45.d0 * u(i-2,j,k)     &      
     &                   -  720.d0 * u(i-1,j,k)     &     
     &                   - 2394.d0 * u(i  ,j,k)     &     
     &                   + 5040.d0 * u(i+1,j,k)     &     
     &                   - 3150.d0 * u(i+2,j,k)     &     
     &                   + 1680.d0 * u(i+3,j,k)     &     
     &                   -  630.d0 * u(i+4,j,k)     &     
     &                   +  144.d0 * u(i+5,j,k)     &     
     &                   -   15.d0 * u(i+6,j,k)     &     
     &                 ) * idx_by_2520
          i=4
          Dxu(i,j,k) = ( -   15.d0 * u(i-3,j,k)     &      
     &                   +  180.d0 * u(i-2,j,k)     &     
     &                   - 1260.d0 * u(i-1,j,k)     &     
     &                   - 1134.d0 * u(i  ,j,k)     &     
     &                   + 3150.d0 * u(i+1,j,k)     &     
     &                   - 1260.d0 * u(i+2,j,k)     &     
     &                   +  420.d0 * u(i+3,j,k)     &     
     &                   -   90.d0 * u(i+4,j,k)     &     
     &                   +    9.d0 * u(i+5,j,k)     &     
     &                 ) * idx_by_2520
          do i = 5, nx-4
            Dxu(i,j,k) = (      9.d0 * u(i-4,j,k)     &      
     &                     -   96.d0 * u(i-3,j,k)     &     
     &                     +  504.d0 * u(i-2,j,k)     &     
     &                     - 2016.d0 * u(i-1,j,k)     &     
     &                     + 2016.d0 * u(i+1,j,k)     &     
     &                     -  504.d0 * u(i+2,j,k)     &     
     &                     +   96.d0 * u(i+3,j,k)     &     
     &                     -    9.d0 * u(i+4,j,k)     &     
     &                   ) * idx_by_2520 
          end do
          i=nx-3
          Dxu(i,j,k) = ( -    9.d0 * u(i-5,j,k)     &      
     &                   +   90.d0 * u(i-4,j,k)     &     
     &                   -  420.d0 * u(i-3,j,k)     &     
     &                   + 1260.d0 * u(i-2,j,k)     &     
     &                   - 3150.d0 * u(i-1,j,k)     &     
     &                   + 1134.d0 * u(i  ,j,k)     &     
     &                   + 1260.d0 * u(i+1,j,k)     &     
     &                   -  180.d0 * u(i+2,j,k)     &     
     &                   +   15.d0 * u(i+3,j,k)     &     
     &                 ) * idx_by_2520 
          i=nx-2 
          Dxu(i,j,k) = (     15.d0 * u(i-6,j,k)     &      
     &                   -  144.d0 * u(i-5,j,k)     &     
     &                   +  630.d0 * u(i-4,j,k)     &     
     &                   - 1680.d0 * u(i-3,j,k)     &     
     &                   + 3150.d0 * u(i-2,j,k)     &     
     &                   - 5040.d0 * u(i-1,j,k)     &     
     &                   + 2394.d0 * u(i  ,j,k)     &     
     &                   +  720.d0 * u(i+1,j,k)     &     
     &                   -   45.d0 * u(i+2,j,k)     &     
     &                 ) * idx_by_2520
          i=nx-1 
          Dxu(i,j,k) = ( -   45.d0 * u(i-7,j,k)     &      
     &                   +  420.d0 * u(i-6,j,k)     &     
     &                   - 1764.d0 * u(i-5,j,k)     &     
     &                   + 4410.d0 * u(i-4,j,k)     &     
     &                   - 7350.d0 * u(i-3,j,k)     &     
     &                   + 8820.d0 * u(i-2,j,k)     &     
     &                   - 8820.d0 * u(i-1,j,k)     &     
     &                   + 4014.d0 * u(i  ,j,k)     &     
     &                   +  315.d0 * u(i+1,j,k)     &     
     &                 ) * idx_by_2520
          i=nx 
          Dxu(i,j,k) = (     315.d0 * u(i-8,j,k)     &      
     &                   -  2880.d0 * u(i-7,j,k)     &     
     &                   + 11760.d0 * u(i-6,j,k)     &     
     &                   - 28224.d0 * u(i-5,j,k)     &     
     &                   + 44100.d0 * u(i-4,j,k)     &     
     &                   - 47040.d0 * u(i-3,j,k)     &     
     &                   + 35280.d0 * u(i-2,j,k)     &     
     &                   - 20160.d0 * u(i-1,j,k)     &     
     &                   +  6849.d0 * u(i  ,j,k)     &     
     &                 ) * idx_by_2520
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      ! deriv8888_y  calculates the first derivative with respect to 
      !              y using a centered eighth order finite difference 
      !              operator on the interval [5,ny-4].  At other points 
      !              appropriated shifted eighth order operators are used.  
      !
      !----------------------------------------------------------------      
      subroutine deriv8888_y(Dyu, u, dy, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12, idy_by_60, idy_by_2520  

        idy = 1.0d0/dy        
        idy_by_2520 = idy / 2520.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  6849.d0 * u(i,j  ,k)     &      
     &                   + 20160.d0 * u(i,j+1,k)     &     
     &                   - 35280.d0 * u(i,j+2,k)     &     
     &                   + 47040.d0 * u(i,j+3,k)     &     
     &                   - 44100.d0 * u(i,j+4,k)     &     
     &                   + 28224.d0 * u(i,j+5,k)     &     
     &                   - 11760.d0 * u(i,j+6,k)     &     
     &                   +  2880.d0 * u(i,j+7,k)     &     
     &                   -   315.d0 * u(i,j+8,k)     &     
     &                 ) * idy_by_2520
          j=2
          Dyu(i,j,k) = ( -  315.d0 * u(i,j-1,k)     &      
     &                   - 4014.d0 * u(i,j  ,k)     &     
     &                   + 8820.d0 * u(i,j+1,k)     &     
     &                   - 8820.d0 * u(i,j+2,k)     &     
     &                   + 7350.d0 * u(i,j+3,k)     &     
     &                   - 4410.d0 * u(i,j+4,k)     &     
     &                   + 1764.d0 * u(i,j+5,k)     &     
     &                   -  420.d0 * u(i,j+6,k)     &     
     &                   +   45.d0 * u(i,j+7,k)     &     
     &                 ) * idy_by_2520
          j=3
          Dyu(i,j,k) = (     45.d0 * u(i,j-2,k)     &      
     &                   -  720.d0 * u(i,j-1,k)     &     
     &                   - 2394.d0 * u(i,j  ,k)     &     
     &                   + 5040.d0 * u(i,j+1,k)     &     
     &                   - 3150.d0 * u(i,j+2,k)     &     
     &                   + 1680.d0 * u(i,j+3,k)     &     
     &                   -  630.d0 * u(i,j+4,k)     &     
     &                   +  144.d0 * u(i,j+5,k)     &     
     &                   -   15.d0 * u(i,j+6,k)     &     
     &                 ) * idy_by_2520
          j=4
          Dyu(i,j,k) = ( -   15.d0 * u(i,j-3,k)     &      
     &                   +  180.d0 * u(i,j-2,k)     &     
     &                   - 1260.d0 * u(i,j-1,k)     &     
     &                   - 1134.d0 * u(i,j  ,k)     &     
     &                   + 3150.d0 * u(i,j+1,k)     &     
     &                   - 1260.d0 * u(i,j+2,k)     &     
     &                   +  420.d0 * u(i,j+3,k)     &     
     &                   -   90.d0 * u(i,j+4,k)     &     
     &                   +    9.d0 * u(i,j+5,k)     &     
     &                 ) * idy_by_2520
          do j = 5, ny-4
            Dyu(i,j,k) = (      9.d0 * u(i,j-4,k)     &      
     &                     -   96.d0 * u(i,j-3,k)     &     
     &                     +  504.d0 * u(i,j-2,k)     &     
     &                     - 2016.d0 * u(i,j-1,k)     &     
     &                     + 2016.d0 * u(i,j+1,k)     &     
     &                     -  504.d0 * u(i,j+2,k)     &     
     &                     +   96.d0 * u(i,j+3,k)     &     
     &                     -    9.d0 * u(i,j+4,k)     &     
     &                   ) * idy_by_2520 
          end do
          j=ny-3
          Dyu(i,j,k) = ( -    9.d0 * u(i,j-5,k)     &      
     &                   +   90.d0 * u(i,j-4,k)     &     
     &                   -  420.d0 * u(i,j-3,k)     &     
     &                   + 1260.d0 * u(i,j-2,k)     &     
     &                   - 3150.d0 * u(i,j-1,k)     &     
     &                   + 1134.d0 * u(i,j  ,k)     &     
     &                   + 1260.d0 * u(i,j+1,k)     &     
     &                   -  180.d0 * u(i,j+2,k)     &     
     &                   +   15.d0 * u(i,j+3,k)     &     
     &                 ) * idy_by_2520 
          j=ny-2 
          Dyu(i,j,k) = (     15.d0 * u(i,j-6,k)     &      
     &                   -  144.d0 * u(i,j-5,k)     &     
     &                   +  630.d0 * u(i,j-4,k)     &     
     &                   - 1680.d0 * u(i,j-3,k)     &     
     &                   + 3150.d0 * u(i,j-2,k)     &     
     &                   - 5040.d0 * u(i,j-1,k)     &     
     &                   + 2394.d0 * u(i,j  ,k)     &     
     &                   +  720.d0 * u(i,j+1,k)     &     
     &                   -   45.d0 * u(i,j+2,k)     &     
     &                 ) * idy_by_2520
          j=ny-1 
          Dyu(i,j,k) = ( -   45.d0 * u(i,j-7,k)     &      
     &                   +  420.d0 * u(i,j-6,k)     &     
     &                   - 1764.d0 * u(i,j-5,k)     &     
     &                   + 4410.d0 * u(i,j-4,k)     &     
     &                   - 7350.d0 * u(i,j-3,k)     &     
     &                   + 8820.d0 * u(i,j-2,k)     &     
     &                   - 8820.d0 * u(i,j-1,k)     &     
     &                   + 4014.d0 * u(i,j  ,k)     &     
     &                   +  315.d0 * u(i,j+1,k)     &     
     &                 ) * idy_by_2520
          j=ny 
          Dyu(i,j,k) = (     315.d0 * u(i,j-8,k)     &      
     &                   -  2880.d0 * u(i,j-7,k)     &     
     &                   + 11760.d0 * u(i,j-6,k)     &     
     &                   - 28224.d0 * u(i,j-5,k)     &     
     &                   + 44100.d0 * u(i,j-4,k)     &     
     &                   - 47040.d0 * u(i,j-3,k)     &     
     &                   + 35280.d0 * u(i,j-2,k)     &     
     &                   - 20160.d0 * u(i,j-1,k)     &     
     &                   +  6849.d0 * u(i,j  ,k)     &     
     &                 ) * idy_by_2520
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      ! deriv8888_z  calculates the first derivative with respect to 
      !              z using a centered eighth order finite difference 
      !              operator on the interval [5,nz-4].  At other points 
      !              appropriated shifted eighth order operators are used.  
      !
      !----------------------------------------------------------------      
      subroutine deriv8888_z(Dzu, u, dz, nx, ny, nz, d_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12, idz_by_60, idz_by_2520  

        idz = 1.0d0/dz        
        idz_by_2520 = idz / 2520.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( -  6849.d0 * u(i,j,k  )     &      
     &                   + 20160.d0 * u(i,j,k+1)     &     
     &                   - 35280.d0 * u(i,j,k+2)     &     
     &                   + 47040.d0 * u(i,j,k+3)     &     
     &                   - 44100.d0 * u(i,j,k+4)     &     
     &                   + 28224.d0 * u(i,j,k+5)     &     
     &                   - 11760.d0 * u(i,j,k+6)     &     
     &                   +  2880.d0 * u(i,j,k+7)     &     
     &                   -   315.d0 * u(i,j,k+8)     &     
     &                 ) * idz_by_2520
          k=2
          Dzu(i,j,k) = ( -  315.d0 * u(i,j,k-1)     &      
     &                   - 4014.d0 * u(i,j,k  )     &     
     &                   + 8820.d0 * u(i,j,k+1)     &     
     &                   - 8820.d0 * u(i,j,k+2)     &     
     &                   + 7350.d0 * u(i,j,k+3)     &     
     &                   - 4410.d0 * u(i,j,k+4)     &     
     &                   + 1764.d0 * u(i,j,k+5)     &     
     &                   -  420.d0 * u(i,j,k+6)     &     
     &                   +   45.d0 * u(i,j,k+7)     &     
     &                 ) * idz_by_2520
          k=3
          Dzu(i,j,k) = (     45.d0 * u(i,j,k-2)     &      
     &                   -  720.d0 * u(i,j,k-1)     &     
     &                   - 2394.d0 * u(i,j,k  )     &     
     &                   + 5040.d0 * u(i,j,k+1)     &     
     &                   - 3150.d0 * u(i,j,k+2)     &     
     &                   + 1680.d0 * u(i,j,k+3)     &     
     &                   -  630.d0 * u(i,j,k+4)     &     
     &                   +  144.d0 * u(i,j,k+5)     &     
     &                   -   15.d0 * u(i,j,k+6)     &     
     &                 ) * idz_by_2520
          k=4
          Dzu(i,j,k) = ( -   15.d0 * u(i,j,k-3)     &      
     &                   +  180.d0 * u(i,j,k-2)     &     
     &                   - 1260.d0 * u(i,j,k-1)     &     
     &                   - 1134.d0 * u(i,j,k  )     &     
     &                   + 3150.d0 * u(i,j,k+1)     &     
     &                   - 1260.d0 * u(i,j,k+2)     &     
     &                   +  420.d0 * u(i,j,k+3)     &     
     &                   -   90.d0 * u(i,j,k+4)     &     
     &                   +    9.d0 * u(i,j,k+5)     &     
     &                 ) * idz_by_2520
          do k = 5, nz-4
            Dzu(i,j,k) = (      9.d0 * u(i,j,k-4)     &      
     &                     -   96.d0 * u(i,j,k-3)     &     
     &                     +  504.d0 * u(i,j,k-2)     &     
     &                     - 2016.d0 * u(i,j,k-1)     &     
     &                     + 2016.d0 * u(i,j,k+1)     &     
     &                     -  504.d0 * u(i,j,k+2)     &     
     &                     +   96.d0 * u(i,j,k+3)     &     
     &                     -    9.d0 * u(i,j,k+4)     &     
     &                   ) * idz_by_2520 
          end do
          k=nz-3
          Dzu(i,j,k) = ( -    9.d0 * u(i,j,k-5)     &      
     &                   +   90.d0 * u(i,j,k-4)     &     
     &                   -  420.d0 * u(i,j,k-3)     &     
     &                   + 1260.d0 * u(i,j,k-2)     &     
     &                   - 3150.d0 * u(i,j,k-1)     &     
     &                   + 1134.d0 * u(i,j,k  )     &     
     &                   + 1260.d0 * u(i,j,k+1)     &     
     &                   -  180.d0 * u(i,j,k+2)     &     
     &                   +   15.d0 * u(i,j,k+3)     &     
     &                 ) * idz_by_2520 
          k=nz-2 
          Dzu(i,j,k) = (     15.d0 * u(i,j,k-6)     &      
     &                   -  144.d0 * u(i,j,k-5)     &     
     &                   +  630.d0 * u(i,j,k-4)     &     
     &                   - 1680.d0 * u(i,j,k-3)     &     
     &                   + 3150.d0 * u(i,j,k-2)     &     
     &                   - 5040.d0 * u(i,j,k-1)     &     
     &                   + 2394.d0 * u(i,j,k  )     &     
     &                   +  720.d0 * u(i,j,k+1)     &     
     &                   -   45.d0 * u(i,j,k+2)     &     
     &                 ) * idz_by_2520
          k=nz-1 
          Dzu(i,j,k) = ( -   45.d0 * u(i,j,k-7)     &      
     &                   +  420.d0 * u(i,j,k-6)     &     
     &                   - 1764.d0 * u(i,j,k-5)     &     
     &                   + 4410.d0 * u(i,j,k-4)     &     
     &                   - 7350.d0 * u(i,j,k-3)     &     
     &                   + 8820.d0 * u(i,j,k-2)     &     
     &                   - 8820.d0 * u(i,j,k-1)     &     
     &                   + 4014.d0 * u(i,j,k  )     &     
     &                   +  315.d0 * u(i,j,k+1)     &     
     &                 ) * idz_by_2520
          k=nz 
          Dzu(i,j,k) = (     315.d0 * u(i,j,k-8)     &      
     &                   -  2880.d0 * u(i,j,k-7)     &     
     &                   + 11760.d0 * u(i,j,k-6)     &     
     &                   - 28224.d0 * u(i,j,k-5)     &     
     &                   + 44100.d0 * u(i,j,k-4)     &     
     &                   - 47040.d0 * u(i,j,k-3)     &     
     &                   + 35280.d0 * u(i,j,k-2)     &     
     &                   - 20160.d0 * u(i,j,k-1)     &     
     &                   +  6849.d0 * u(i,j,k  )     &     
     &                 ) * idz_by_2520
        end do
        end do
       
        return
      end subroutine 



      !----------------------------------------------------------------
      !
      ! deriv8642adv_x calculates the first derivative with respect to 
      !                x using an "advective" eighth order finite 
      !                difference operator at all points on the 
      !                interval [6,nx-5].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^x \partial_x in the equations.  
      !                Further, it will depend on the sign of beta^x.  
      !                Namely, if beta^x is positive, we will use a 
      !                eighth order upwind operator while if beta^x 
      !                is negative, we use a eighth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At i=5, if beta^x is positive, we continue to  
      !                use the eighth order upwind operator but if  
      !                beta^x is negative, we use a sixth order downwind  
      !                operator.   
      !                  
      !                We reverse this at i=nx-4 such that if beta^x is  
      !                positive, we use a sixth order upwind operator 
      !                while if beta^x is negative we continue to use  
      !                the eighth order downwind operator.   
      !                 
      !                At i=4, if beta^x is positive, we use the   
      !                sixth order upwind operator and if beta^x is  
      !                negative, we use a fourth order downwind operator. 
      !                 
      !                At i=nx-3, if beta^x is positive, we use the   
      !                fourth order upwind operator and if beta^x is  
      !                negative, we use a sixth order downwind operator. 
      !                 
      !                At i=3, if beta^x is positive, we use the   
      !                fourth order upwind operator and if beta^x is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At i=nx-2, if beta^x is positive, we use the   
      !                second order upwind operator and if beta^x is  
      !                negative, we use a fourth order downwind operator. 
      !                 
      !                At i=2, if beta^x is positive, we use the   
      !                second order upwind operator and if beta^x is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At i=nx-1, if beta^x is positive, we use the   
      !                second order centered operator and if beta^x is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At i=1 we use the second order upwind operator 
      !                for both positive and negative values of beta^x.  
      !                 
      !                At i=nx, we use the second order downwind operator 
      !                for both positive and negative values of beta^x.  
      !                 
      !----------------------------------------------------------------      
      subroutine deriv8642adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz), betax(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60, idx_by_2520  

        idx = 1.0d0/dx        
        idx_by_2 = 0.5d0 * idx
        idx_by_12 = idx / 12.d0
        idx_by_60 = idx / 60.d0
        idx_by_2520 = idx / 2520.d0


        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  3.d0 * u(i  ,j,k)     &      
     &                   +  4.d0 * u(i+1,j,k)     &     
     &                   -         u(i+2,j,k)     &     
     &                 ) * idx_by_2

          i=2
          if (betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( - 3.d0 * u(i  ,j,k)     &      
     &                     + 4.d0 * u(i+1,j,k)     &     
     &                     -        u(i+2,j,k)     &     
     &                   ) * idx_by_2 
          else 
            Dxu(i,j,k) = ( -         u(i-1,j,k)     &      
     &                     +         u(i+1,j,k)     &      
     &                   ) * idx_by_2 
          endif 

          i=3
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                     - 10.d0 * u(i  ,j,k)     &     
     &                     + 18.d0 * u(i+1,j,k)     &     
     &                     -  6.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_12 
          else 
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                     -  4.d0 * u(i-1,j,k)     &     
     &                     +  3.d0 * u(i  ,j,k)     &     
     &                   ) * idx_by_2 
          endif

          i=4
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (    2.d0 * u(i-2,j,k)     &      
     &                     - 24.d0 * u(i-1,j,k)     &     
     &                     - 35.d0 * u(i  ,j,k)     &     
     &                     + 80.d0 * u(i+1,j,k)     &     
     &                     - 30.d0 * u(i+2,j,k)     &     
     &                     +  8.d0 * u(i+3,j,k)     &     
     &                     -         u(i+4,j,k)     &     
     &                   ) * idx_by_60 
          else 
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  6.d0 * u(i-2,j,k)     &     
     &                     - 18.d0 * u(i-1,j,k)     &     
     &                     + 10.d0 * u(i  ,j,k)     &     
     &                     +  3.d0 * u(i+1,j,k)     &     
     &                   ) * idx_by_12 
          endif
          
          i=5
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -   15.d0 * u(i-3,j,k)     &      
     &                     +  180.d0 * u(i-2,j,k)     &     
     &                     - 1260.d0 * u(i-1,j,k)     &     
     &                     - 1134.d0 * u(i  ,j,k)     &     
     &                     + 3150.d0 * u(i+1,j,k)     &     
     &                     - 1260.d0 * u(i+2,j,k)     &     
     &                     +  420.d0 * u(i+3,j,k)     &     
     &                     -   90.d0 * u(i+4,j,k)     &     
     &                     +    9.d0 * u(i+5,j,k)     &     
     &                   ) * idx_by_2520 
          else 
            Dxu(i,j,k) = (           u(i-4,j,k)     &      
     &                     -  8.d0 * u(i-3,j,k)     &     
     &                     + 30.d0 * u(i-2,j,k)     &     
     &                     - 80.d0 * u(i-1,j,k)     &     
     &                     + 35.d0 * u(i  ,j,k)     &     
     &                     + 24.d0 * u(i+1,j,k)     &     
     &                     -  2.d0 * u(i+2,j,k)     &     
     &                   ) * idx_by_60 
          endif 
          
          do i = 6, nx-5  ! see gr-qc/0505055v2.pdf 
            if ( betax(i,j,k) .ge. 0.d0 ) then
              Dxu(i,j,k) = ( -   15.d0 * u(i-3,j,k)     &      
     &                       +  180.d0 * u(i-2,j,k)     &     
     &                       - 1260.d0 * u(i-1,j,k)     &     
     &                       - 1134.d0 * u(i  ,j,k)     &     
     &                       + 3150.d0 * u(i+1,j,k)     &     
     &                       - 1260.d0 * u(i+2,j,k)     &     
     &                       +  420.d0 * u(i+3,j,k)     &     
     &                       -   90.d0 * u(i+4,j,k)     &     
     &                       +    9.d0 * u(i+5,j,k)     &     
     &                     ) * idx_by_2520 
            else 
              Dxu(i,j,k) = ( -    9.d0 * u(i-5,j,k)     &      
     &                       +   90.d0 * u(i-4,j,k)     &     
     &                       -  420.d0 * u(i-3,j,k)     &     
     &                       + 1260.d0 * u(i-2,j,k)     &     
     &                       - 3150.d0 * u(i-1,j,k)     &     
     &                       + 1134.d0 * u(i  ,j,k)     &     
     &                       + 1260.d0 * u(i+1,j,k)     &     
     &                       -  180.d0 * u(i+2,j,k)     &     
     &                       +   15.d0 * u(i+3,j,k)     &     
     &                     ) * idx_by_2520 
            endif 
          end do
          
          i=nx-4 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (    2.d0 * u(i-2,j,k)     &      
     &                     - 24.d0 * u(i-1,j,k)     &     
     &                     - 35.d0 * u(i  ,j,k)     &     
     &                     + 80.d0 * u(i+1,j,k)     &     
     &                     - 30.d0 * u(i+2,j,k)     &     
     &                     +  8.d0 * u(i+3,j,k)     &     
     &                     -         u(i+4,j,k)     &     
     &                   ) * idx_by_60 
          else 
            Dxu(i,j,k) = ( -    9.d0 * u(i-5,j,k)     &      
     &                     +   90.d0 * u(i-4,j,k)     &     
     &                     -  420.d0 * u(i-3,j,k)     &     
     &                     + 1260.d0 * u(i-2,j,k)     &     
     &                     - 3150.d0 * u(i-1,j,k)     &     
     &                     + 1134.d0 * u(i  ,j,k)     &     
     &                     + 1260.d0 * u(i+1,j,k)     &     
     &                     -  180.d0 * u(i+2,j,k)     &     
     &                     +   15.d0 * u(i+3,j,k)     &     
     &                   ) * idx_by_2520 
          endif 

          i=nx-3 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -  3.d0 * u(i-1,j,k)     &      
     &                     - 10.d0 * u(i  ,j,k)     &     
     &                     + 18.d0 * u(i+1,j,k)     &     
     &                     -  6.d0 * u(i+2,j,k)     &     
     &                     +         u(i+3,j,k)     &     
     &                   ) * idx_by_12 
          else 
            Dxu(i,j,k) = (           u(i-4,j,k)     &      
     &                     -  8.d0 * u(i-3,j,k)     &     
     &                     + 30.d0 * u(i-2,j,k)     &     
     &                     - 80.d0 * u(i-1,j,k)     &     
     &                     + 35.d0 * u(i  ,j,k)     &     
     &                     + 24.d0 * u(i+1,j,k)     &     
     &                     -  2.d0 * u(i+2,j,k)     &     
     &                   ) * idx_by_60 
          endif 

          i=nx-2 
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (  - 3.d0 * u(i  ,j,k)     &      
     &                      + 4.d0 * u(i+1,j,k)     &      
     &                      -        u(i+2,j,k)     &     
     &                   ) * idx_by_2 
          else 
            Dxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  6.d0 * u(i-2,j,k)     &     
     &                     - 18.d0 * u(i-1,j,k)     &     
     &                     + 10.d0 * u(i  ,j,k)     &     
     &                     +  3.d0 * u(i+1,j,k)     &     
     &                   ) * idx_by_12 
          endif

          i=nx-1 
          if ( betax(i,j,k) .ge. 0.d0 ) then 
            Dxu(i,j,k) = (  -         u(i-1,j,k)     &      
     &                      +         u(i+1,j,k)     &     
     &                    ) * idx_by_2
          else 
            Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                      - 4.d0 * u(i-1,j,k)     &      
     &                      + 3.d0 * u(i  ,j,k)     &      
     &                   ) * idx_by_2 
          endif

          i=nx 
          Dxu(i,j,k) = (           u(i-2,j,k)     &      
     &                    - 4.d0 * u(i-1,j,k)     &      
     &                    + 3.d0 * u(i  ,j,k)     &      
     &                 ) * idx_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      


      !----------------------------------------------------------------
      !
      ! deriv8642adv_y calculates the first derivative with respect to 
      !                y using an "advective" eighth order finite 
      !                difference operator at all points on the 
      !                interval [6,ny-5].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^y \partial_y in the equations.  
      !                Further, it will depend on the sign of beta^y.  
      !                Namely, if beta^y is positive, we will use a 
      !                eighth order upwind operator while if beta^y 
      !                is negative, we use a eighth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At j=5, if beta^y is positive, we continue to  
      !                use the eighth order upwind operator but if  
      !                beta^y is negative, we use a sixth order downwind  
      !                operator.   
      !                  
      !                We reverse this at j=ny-4 such that if beta^y is  
      !                positive, we use a sixth order upwind operator 
      !                while if beta^y is negative we continue to use  
      !                the eighth order downwind operator.   
      !                 
      !                At j=4, if beta^y is positive, we use the   
      !                sixth order upwind operator and if beta^y is  
      !                negative, we use a fourth order downwind operator. 
      !                 
      !                At j=ny-3, if beta^y is positive, we use the   
      !                fourth order upwind operator and if beta^y is  
      !                negative, we use a sixth order downwind operator. 
      !                 
      !                At j=3, if beta^y is positive, we use the   
      !                fourth order upwind operator and if beta^y is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At j=ny-2, if beta^y is positive, we use the   
      !                second order upwind operator and if beta^y is  
      !                negative, we use a fourth order downwind operator. 
      !                 
      !                At j=2, if beta^y is positive, we use the   
      !                second order upwind operator and if beta^y is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At j=ny-1, if beta^y is positive, we use the   
      !                second order centered operator and if beta^y is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At j=1 we use the second order upwind operator 
      !                for both positive and negative values of beta^y.  
      !                 
      !                At j=ny, we use the second order downwind operator 
      !                for both positive and negative values of beta^y.  
      !                 
      !----------------------------------------------------------------      
      subroutine deriv8642adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz), betay(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12, idy_by_60, idy_by_2520  

        idy = 1.0d0/dy        
        idy_by_2 = 0.5d0 * idy
        idy_by_12 = idy / 12.d0
        idy_by_60 = idy / 60.d0
        idy_by_2520 = idy / 2520.d0


        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  3.d0 * u(i,j  ,k)     &      
     &                   +  4.d0 * u(i,j+1,k)     &     
     &                   -         u(i,j+2,k)     &     
     &                 ) * idy_by_2

          j=2
          if (betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( - 3.d0 * u(i,j  ,k)     &      
     &                     + 4.d0 * u(i,j+1,k)     &     
     &                     -        u(i,j+2,k)     &     
     &                   ) * idy_by_2 
          else 
            Dyu(i,j,k) = ( -         u(i,j-1,k)     &      
     &                     +         u(i,j+1,k)     &      
     &                   ) * idy_by_2 
          endif 

          j=3
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                     - 10.d0 * u(i,j  ,k)     &     
     &                     + 18.d0 * u(i,j+1,k)     &     
     &                     -  6.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_12 
          else 
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                     -  4.d0 * u(i,j-1,k)     &     
     &                     +  3.d0 * u(i,j  ,k)     &     
     &                   ) * idy_by_2 
          endif

          j=4
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (    2.d0 * u(i,j-2,k)     &      
     &                     - 24.d0 * u(i,j-1,k)     &     
     &                     - 35.d0 * u(i,j  ,k)     &     
     &                     + 80.d0 * u(i,j+1,k)     &     
     &                     - 30.d0 * u(i,j+2,k)     &     
     &                     +  8.d0 * u(i,j+3,k)     &     
     &                     -         u(i,j+4,k)     &     
     &                   ) * idy_by_60 
          else 
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  6.d0 * u(i,j-2,k)     &     
     &                     - 18.d0 * u(i,j-1,k)     &     
     &                     + 10.d0 * u(i,j  ,k)     &     
     &                     +  3.d0 * u(i,j+1,k)     &     
     &                   ) * idy_by_12 
          endif
          
          j=5
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -   15.d0 * u(i,j-3,k)     &      
     &                     +  180.d0 * u(i,j-2,k)     &     
     &                     - 1260.d0 * u(i,j-1,k)     &     
     &                     - 1134.d0 * u(i,j  ,k)     &     
     &                     + 3150.d0 * u(i,j+1,k)     &     
     &                     - 1260.d0 * u(i,j+2,k)     &     
     &                     +  420.d0 * u(i,j+3,k)     &     
     &                     -   90.d0 * u(i,j+4,k)     &     
     &                     +    9.d0 * u(i,j+5,k)     &     
     &                   ) * idy_by_2520 
          else 
            Dyu(i,j,k) = (           u(i,j-4,k)     &      
     &                     -  8.d0 * u(i,j-3,k)     &     
     &                     + 30.d0 * u(i,j-2,k)     &     
     &                     - 80.d0 * u(i,j-1,k)     &     
     &                     + 35.d0 * u(i,j  ,k)     &     
     &                     + 24.d0 * u(i,j+1,k)     &     
     &                     -  2.d0 * u(i,j+2,k)     &     
     &                   ) * idy_by_60 
          endif 
          
          do j = 6, ny-5  ! see gr-qc/0505055v2.pdf 
            if ( betay(i,j,k) .ge. 0.d0 ) then
              Dyu(i,j,k) = ( -   15.d0 * u(i,j-3,k)     &      
     &                       +  180.d0 * u(i,j-2,k)     &     
     &                       - 1260.d0 * u(i,j-1,k)     &     
     &                       - 1134.d0 * u(i,j  ,k)     &     
     &                       + 3150.d0 * u(i,j+1,k)     &     
     &                       - 1260.d0 * u(i,j+2,k)     &     
     &                       +  420.d0 * u(i,j+3,k)     &     
     &                       -   90.d0 * u(i,j+4,k)     &     
     &                       +    9.d0 * u(i,j+5,k)     &     
     &                     ) * idy_by_2520 
            else 
              Dyu(i,j,k) = ( -    9.d0 * u(i,j-5,k)     &      
     &                       +   90.d0 * u(i,j-4,k)     &     
     &                       -  420.d0 * u(i,j-3,k)     &     
     &                       + 1260.d0 * u(i,j-2,k)     &     
     &                       - 3150.d0 * u(i,j-1,k)     &     
     &                       + 1134.d0 * u(i,j  ,k)     &     
     &                       + 1260.d0 * u(i,j+1,k)     &     
     &                       -  180.d0 * u(i,j+2,k)     &     
     &                       +   15.d0 * u(i,j+3,k)     &     
     &                     ) * idy_by_2520 
            endif 
          end do
          
          j=ny-4 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (    2.d0 * u(i,j-2,k)     &      
     &                     - 24.d0 * u(i,j-1,k)     &     
     &                     - 35.d0 * u(i,j  ,k)     &     
     &                     + 80.d0 * u(i,j+1,k)     &     
     &                     - 30.d0 * u(i,j+2,k)     &     
     &                     +  8.d0 * u(i,j+3,k)     &     
     &                     -         u(i,j+4,k)     &     
     &                   ) * idy_by_60 
          else 
            Dyu(i,j,k) = ( -    9.d0 * u(i,j-5,k)     &      
     &                     +   90.d0 * u(i,j-4,k)     &     
     &                     -  420.d0 * u(i,j-3,k)     &     
     &                     + 1260.d0 * u(i,j-2,k)     &     
     &                     - 3150.d0 * u(i,j-1,k)     &     
     &                     + 1134.d0 * u(i,j  ,k)     &     
     &                     + 1260.d0 * u(i,j+1,k)     &     
     &                     -  180.d0 * u(i,j+2,k)     &     
     &                     +   15.d0 * u(i,j+3,k)     &     
     &                   ) * idy_by_2520 
          endif 

          j=ny-3 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -  3.d0 * u(i,j-1,k)     &      
     &                     - 10.d0 * u(i,j  ,k)     &     
     &                     + 18.d0 * u(i,j+1,k)     &     
     &                     -  6.d0 * u(i,j+2,k)     &     
     &                     +         u(i,j+3,k)     &     
     &                   ) * idy_by_12 
          else 
            Dyu(i,j,k) = (           u(i,j-4,k)     &      
     &                     -  8.d0 * u(i,j-3,k)     &     
     &                     + 30.d0 * u(i,j-2,k)     &     
     &                     - 80.d0 * u(i,j-1,k)     &     
     &                     + 35.d0 * u(i,j  ,k)     &     
     &                     + 24.d0 * u(i,j+1,k)     &     
     &                     -  2.d0 * u(i,j+2,k)     &     
     &                   ) * idy_by_60 
          endif 

          j=ny-2 
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (  - 3.d0 * u(i,j  ,k)     &      
     &                      + 4.d0 * u(i,j+1,k)     &      
     &                      -        u(i,j+2,k)     &     
     &                   ) * idy_by_2 
          else 
            Dyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  6.d0 * u(i,j-2,k)     &     
     &                     - 18.d0 * u(i,j-1,k)     &     
     &                     + 10.d0 * u(i,j  ,k)     &     
     &                     +  3.d0 * u(i,j+1,k)     &     
     &                   ) * idy_by_12 
          endif

          j=ny-1 
          if ( betay(i,j,k) .ge. 0.d0 ) then 
            Dyu(i,j,k) = (  -         u(i,j-1,k)     &      
     &                      +         u(i,j+1,k)     &     
     &                    ) * idy_by_2
          else 
            Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                      - 4.d0 * u(i,j-1,k)     &      
     &                      + 3.d0 * u(i,j  ,k)     &      
     &                   ) * idy_by_2 
          endif

          j=ny 
          Dyu(i,j,k) = (           u(i,j-2,k)     &      
     &                    - 4.d0 * u(i,j-1,k)     &      
     &                    + 3.d0 * u(i,j  ,k)     &      
     &                 ) * idy_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      


      !----------------------------------------------------------------
      !
      ! deriv8642adv_z calculates the first derivative with respect to 
      !                z using an "advective" eighth order finite 
      !                difference operator at all points on the 
      !                interval [6,nz-5].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^z \partial_z in the equations.  
      !                Further, it will depend on the sign of beta^z.  
      !                Namely, if beta^z is positive, we will use a 
      !                eighth order upwind operator while if beta^z 
      !                is negative, we use a eighth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At k=5, if beta^z is positive, we continue to  
      !                use the eighth order upwind operator but if  
      !                beta^z is negative, we use a sixth order downwind  
      !                operator.   
      !                  
      !                We reverse this at k=nz-4 such that if beta^z is  
      !                positive, we use a sixth order upwind operator 
      !                while if beta^z is negative we continue to use  
      !                the eighth order downwind operator.   
      !                 
      !                At k=4, if beta^z is positive, we use the   
      !                sixth order upwind operator and if beta^z is  
      !                negative, we use a fourth order downwind operator. 
      !                 
      !                At k=nz-3, if beta^z is positive, we use the   
      !                fourth order upwind operator and if beta^z is  
      !                negative, we use a sixth order downwind operator. 
      !                 
      !                At k=3, if beta^z is positive, we use the   
      !                fourth order upwind operator and if beta^z is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At k=nz-2, if beta^z is positive, we use the   
      !                second order upwind operator and if beta^z is  
      !                negative, we use a fourth order downwind operator. 
      !                 
      !                At k=2, if beta^z is positive, we use the   
      !                second order upwind operator and if beta^z is  
      !                negative, we use a second order centered operator. 
      !                 
      !                At k=nz-1, if beta^z is positive, we use the   
      !                second order centered operator and if beta^z is  
      !                negative, we use a second order downwind operator. 
      !                 
      !                At k=1 we use the second order upwind operator 
      !                for both positive and negative values of beta^z.  
      !                 
      !                At k=nz, we use the second order downwind operator 
      !                for both positive and negative values of beta^z.  
      !                 
      !----------------------------------------------------------------      
      subroutine deriv8642adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz), betaz(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12, idz_by_60, idz_by_2520  

        idz = 1.0d0/dz        
        idz_by_2 = 0.5d0 * idz
        idz_by_12 = idz / 12.d0
        idz_by_60 = idz / 60.d0
        idz_by_2520 = idz / 2520.d0


        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( -  3.d0 * u(i,j,k  )     &      
     &                   +  4.d0 * u(i,j,k+1)     &     
     &                   -         u(i,j,k+2)     &     
     &                 ) * idz_by_2

          k=2
          if (betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( - 3.d0 * u(i,j,k  )     &      
     &                     + 4.d0 * u(i,j,k+1)     &     
     &                     -        u(i,j,k+2)     &     
     &                   ) * idz_by_2 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-1)     &      
     &                     +         u(i,j,k+1)     &      
     &                   ) * idz_by_2 
          endif 

          k=3
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                     - 10.d0 * u(i,j,k  )     &     
     &                     + 18.d0 * u(i,j,k+1)     &     
     &                     -  6.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_12 
          else 
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                     -  4.d0 * u(i,j,k-1)     &     
     &                     +  3.d0 * u(i,j,k  )     &     
     &                   ) * idz_by_2 
          endif

          k=4
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (    2.d0 * u(i,j,k-2)     &      
     &                     - 24.d0 * u(i,j,k-1)     &     
     &                     - 35.d0 * u(i,j,k  )     &     
     &                     + 80.d0 * u(i,j,k+1)     &     
     &                     - 30.d0 * u(i,j,k+2)     &     
     &                     +  8.d0 * u(i,j,k+3)     &     
     &                     -         u(i,j,k+4)     &     
     &                   ) * idz_by_60 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  6.d0 * u(i,j,k-2)     &     
     &                     - 18.d0 * u(i,j,k-1)     &     
     &                     + 10.d0 * u(i,j,k  )     &     
     &                     +  3.d0 * u(i,j,k+1)     &     
     &                   ) * idz_by_12 
          endif
          
          k=5
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -   15.d0 * u(i,j,k-3)     &      
     &                     +  180.d0 * u(i,j,k-2)     &     
     &                     - 1260.d0 * u(i,j,k-1)     &     
     &                     - 1134.d0 * u(i,j,k  )     &     
     &                     + 3150.d0 * u(i,j,k+1)     &     
     &                     - 1260.d0 * u(i,j,k+2)     &     
     &                     +  420.d0 * u(i,j,k+3)     &     
     &                     -   90.d0 * u(i,j,k+4)     &     
     &                     +    9.d0 * u(i,j,k+5)     &     
     &                   ) * idz_by_2520 
          else 
            Dzu(i,j,k) = (           u(i,j,k-4)     &      
     &                     -  8.d0 * u(i,j,k-3)     &     
     &                     + 30.d0 * u(i,j,k-2)     &     
     &                     - 80.d0 * u(i,j,k-1)     &     
     &                     + 35.d0 * u(i,j,k  )     &     
     &                     + 24.d0 * u(i,j,k+1)     &     
     &                     -  2.d0 * u(i,j,k+2)     &     
     &                   ) * idz_by_60 
          endif 
          
          do k = 6, nz-5  ! see gr-qc/0505055v2.pdf 
            if ( betaz(i,j,k) .ge. 0.d0 ) then
              Dzu(i,j,k) = ( -   15.d0 * u(i,j,k-3)     &      
     &                       +  180.d0 * u(i,j,k-2)     &     
     &                       - 1260.d0 * u(i,j,k-1)     &     
     &                       - 1134.d0 * u(i,j,k  )     &     
     &                       + 3150.d0 * u(i,j,k+1)     &     
     &                       - 1260.d0 * u(i,j,k+2)     &     
     &                       +  420.d0 * u(i,j,k+3)     &     
     &                       -   90.d0 * u(i,j,k+4)     &     
     &                       +    9.d0 * u(i,j,k+5)     &     
     &                     ) * idz_by_2520 
            else 
              Dzu(i,j,k) = ( -    9.d0 * u(i,j,k-5)     &      
     &                       +   90.d0 * u(i,j,k-4)     &     
     &                       -  420.d0 * u(i,j,k-3)     &     
     &                       + 1260.d0 * u(i,j,k-2)     &     
     &                       - 3150.d0 * u(i,j,k-1)     &     
     &                       + 1134.d0 * u(i,j,k  )     &     
     &                       + 1260.d0 * u(i,j,k+1)     &     
     &                       -  180.d0 * u(i,j,k+2)     &     
     &                       +   15.d0 * u(i,j,k+3)     &     
     &                     ) * idz_by_2520 
            endif 
          end do
          
          k=nz-4 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (    2.d0 * u(i,j,k-2)     &      
     &                     - 24.d0 * u(i,j,k-1)     &     
     &                     - 35.d0 * u(i,j,k  )     &     
     &                     + 80.d0 * u(i,j,k+1)     &     
     &                     - 30.d0 * u(i,j,k+2)     &     
     &                     +  8.d0 * u(i,j,k+3)     &     
     &                     -         u(i,j,k+4)     &     
     &                   ) * idz_by_60 
          else 
            Dzu(i,j,k) = ( -    9.d0 * u(i,j,k-5)     &      
     &                     +   90.d0 * u(i,j,k-4)     &     
     &                     -  420.d0 * u(i,j,k-3)     &     
     &                     + 1260.d0 * u(i,j,k-2)     &     
     &                     - 3150.d0 * u(i,j,k-1)     &     
     &                     + 1134.d0 * u(i,j,k  )     &     
     &                     + 1260.d0 * u(i,j,k+1)     &     
     &                     -  180.d0 * u(i,j,k+2)     &     
     &                     +   15.d0 * u(i,j,k+3)     &     
     &                   ) * idz_by_2520 
          endif 

          k=nz-3 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -  3.d0 * u(i,j,k-1)     &      
     &                     - 10.d0 * u(i,j,k  )     &     
     &                     + 18.d0 * u(i,j,k+1)     &     
     &                     -  6.d0 * u(i,j,k+2)     &     
     &                     +         u(i,j,k+3)     &     
     &                   ) * idz_by_12 
          else 
            Dzu(i,j,k) = (           u(i,j,k-4)     &      
     &                     -  8.d0 * u(i,j,k-3)     &     
     &                     + 30.d0 * u(i,j,k-2)     &     
     &                     - 80.d0 * u(i,j,k-1)     &     
     &                     + 35.d0 * u(i,j,k  )     &     
     &                     + 24.d0 * u(i,j,k+1)     &     
     &                     -  2.d0 * u(i,j,k+2)     &     
     &                   ) * idz_by_60 
          endif 

          k=nz-2 
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (  - 3.d0 * u(i,j,k  )     &      
     &                      + 4.d0 * u(i,j,k+1)     &      
     &                      -        u(i,j,k+2)     &     
     &                   ) * idz_by_2 
          else 
            Dzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  6.d0 * u(i,j,k-2)     &     
     &                     - 18.d0 * u(i,j,k-1)     &     
     &                     + 10.d0 * u(i,j,k  )     &     
     &                     +  3.d0 * u(i,j,k+1)     &     
     &                   ) * idz_by_12 
          endif

          k=nz-1 
          if ( betaz(i,j,k) .ge. 0.d0 ) then 
            Dzu(i,j,k) = (  -         u(i,j,k-1)     &      
     &                      +         u(i,j,k+1)     &     
     &                    ) * idz_by_2
          else 
            Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                      - 4.d0 * u(i,j,k-1)     &      
     &                      + 3.d0 * u(i,j,k  )     &      
     &                   ) * idz_by_2 
          endif

          k=nz 
          Dzu(i,j,k) = (           u(i,j,k-2)     &      
     &                    - 4.d0 * u(i,j,k-1)     &      
     &                    + 3.d0 * u(i,j,k  )     &      
     &                 ) * idz_by_2 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv8888adv_x calculates the first derivative with respect to 
      !                x using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [6,nx-5].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^x \partial_x in the equations.  
      !                Further, it will depend on the sign of beta^x.  
      !                Namely, if beta^x is positive, we will use a 
      !                fourth order upwind operator while if beta^x 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At i=5, if beta^x is positive, we continue to  
      !                use the eighth order upwind operator but if  
      !                beta^x is negative, we use the eighth order downwind  
      !                operator.   
      !                  
      !                We reverse this at i=nx-4 such that if beta^x is  
      !                positive, we use an eighth order centered operator 
      !                while if beta^x is negative we continue to use  
      !                the eighth order downwind operator.   
      !                 
      !                At i=4, we want an eighth order operator so it
      !                does not matter what the sign of beta^x is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of 
      !                using the downwind operator for the case of 
      !                negative beta^x, but I do not see any way around 
      !                this if we want to have an eighth order operator.  
      !                 
      !                Likewise, at i=nx-3, the requirement of an eighth 
      !                order operator renders the sign of beta^x 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^x is positive. 
      !                 
      !                At i=2,3 and i=nx-2,nx-1, we use differently shifted 
      !                eighth order operators irrespective of beta^xs
      !                sign.  
      !                
      !                At i=1 and i=nx, we use the totally shifted eighth
      !                order operators irrespective of beta^xs sign.  
      !            
      !----------------------------------------------------------------      
      subroutine deriv8888adv_x(Dxu, u, dx, nx, ny, nz, d_type, betax)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dxu(nx,ny,nz), u(nx,ny,nz), betax(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_by_2, idx_by_12, idx_by_60, idx_by_2520  

        idx = 1.0d0/dx        
        idx_by_2 = 0.5d0 * idx
        idx_by_12 = idx / 12.d0
        idx_by_60 = idx / 60.d0
        idx_by_2520 = idx / 2520.d0


        do k = 1, nz
        do j = 1, ny
          i=1
          Dxu(i,j,k) = ( -  6849.d0 * u(i  ,j,k)     &      
     &                   + 20160.d0 * u(i+1,j,k)     &     
     &                   - 35280.d0 * u(i+2,j,k)     &     
     &                   + 47040.d0 * u(i+3,j,k)     &     
     &                   - 44100.d0 * u(i+4,j,k)     &     
     &                   + 28224.d0 * u(i+5,j,k)     &     
     &                   - 11760.d0 * u(i+6,j,k)     &     
     &                   +  2880.d0 * u(i+7,j,k)     &     
     &                   -   315.d0 * u(i+8,j,k)     &     
     &                 ) * idx_by_2520 

          i=2
          Dxu(i,j,k) = ( -  315.d0 * u(i-1,j,k)     &      
     &                   - 4014.d0 * u(i  ,j,k)     &     
     &                   + 8820.d0 * u(i+1,j,k)     &     
     &                   - 8820.d0 * u(i+2,j,k)     &     
     &                   + 7350.d0 * u(i+3,j,k)     &     
     &                   - 4410.d0 * u(i+4,j,k)     &     
     &                   + 1764.d0 * u(i+5,j,k)     &     
     &                   -  420.d0 * u(i+6,j,k)     &     
     &                   +   45.d0 * u(i+7,j,k)     &     
     &                 ) * idx_by_2520 

          i=3
          Dxu(i,j,k) = (     45.d0 * u(i-2,j,k)     &      
     &                   -  720.d0 * u(i-1,j,k)     &     
     &                   - 2394.d0 * u(i  ,j,k)     &     
     &                   + 5040.d0 * u(i+1,j,k)     &     
     &                   - 3150.d0 * u(i+2,j,k)     &     
     &                   + 1680.d0 * u(i+3,j,k)     &     
     &                   -  630.d0 * u(i+4,j,k)     &     
     &                   +  144.d0 * u(i+5,j,k)     &     
     &                   -   15.d0 * u(i+6,j,k)     &     
     &                 ) * idx_by_2520 

          i=4
          Dxu(i,j,k) = ( -   15.d0 * u(i-3,j,k)     &      
     &                   +  180.d0 * u(i-2,j,k)     &     
     &                   - 1260.d0 * u(i-1,j,k)     &     
     &                   - 1134.d0 * u(i  ,j,k)     &     
     &                   + 3150.d0 * u(i+1,j,k)     &     
     &                   - 1260.d0 * u(i+2,j,k)     &     
     &                   +  420.d0 * u(i+3,j,k)     &     
     &                   -   90.d0 * u(i+4,j,k)     &     
     &                   +    9.d0 * u(i+5,j,k)     &     
     &                 ) * idx_by_2520 
          
          i=5  
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = ( -   15.d0 * u(i-3,j,k)     &      
     &                     +  180.d0 * u(i-2,j,k)     &     
     &                     - 1260.d0 * u(i-1,j,k)     &     
     &                     - 1134.d0 * u(i  ,j,k)     &     
     &                     + 3150.d0 * u(i+1,j,k)     &     
     &                     - 1260.d0 * u(i+2,j,k)     &     
     &                     +  420.d0 * u(i+3,j,k)     &     
     &                     -   90.d0 * u(i+4,j,k)     &     
     &                     +    9.d0 * u(i+5,j,k)     &     
     &                   ) * idx_by_2520 
          else 
            Dxu(i,j,k) = (      9.d0 * u(i-4,j,k)     &      
     &                     -   96.d0 * u(i-3,j,k)     &     
     &                     +  504.d0 * u(i-2,j,k)     &     
     &                     - 2016.d0 * u(i-1,j,k)     &     
     &                     + 2016.d0 * u(i+1,j,k)     &     
     &                     -  504.d0 * u(i+2,j,k)     &     
     &                     +   96.d0 * u(i+3,j,k)     &     
     &                     -    9.d0 * u(i+4,j,k)     &     
     &                   ) * idx_by_2520 
          endif 

          do i = 6, nx-5  ! see gr-qc/0505055v2.pdf 
            if ( betax(i,j,k) .ge. 0.d0 ) then
              Dxu(i,j,k) = ( -   15.d0 * u(i-3,j,k)     &      
     &                       +  180.d0 * u(i-2,j,k)     &     
     &                       - 1260.d0 * u(i-1,j,k)     &     
     &                       - 1134.d0 * u(i  ,j,k)     &     
     &                       + 3150.d0 * u(i+1,j,k)     &     
     &                       - 1260.d0 * u(i+2,j,k)     &     
     &                       +  420.d0 * u(i+3,j,k)     &     
     &                       -   90.d0 * u(i+4,j,k)     &     
     &                       +    9.d0 * u(i+5,j,k)     &     
     &                     ) * idx_by_2520 
            else 
              Dxu(i,j,k) = ( -    9.d0 * u(i-5,j,k)     &      
     &                       +   90.d0 * u(i-4,j,k)     &     
     &                       -  420.d0 * u(i-3,j,k)     &     
     &                       + 1260.d0 * u(i-2,j,k)     &     
     &                       - 3150.d0 * u(i-1,j,k)     &     
     &                       + 1134.d0 * u(i  ,j,k)     &     
     &                       + 1260.d0 * u(i+1,j,k)     &     
     &                       -  180.d0 * u(i+2,j,k)     &     
     &                       +   15.d0 * u(i+3,j,k)     &     
     &                     ) * idx_by_2520 
            endif 
          end do

          i=nx-4  
          if ( betax(i,j,k) .ge. 0.d0 ) then
            Dxu(i,j,k) = (      9.d0 * u(i-4,j,k)     &      
     &                     -   96.d0 * u(i-3,j,k)     &     
     &                     +  504.d0 * u(i-2,j,k)     &     
     &                     - 2016.d0 * u(i-1,j,k)     &     
     &                     + 2016.d0 * u(i+1,j,k)     &     
     &                     -  504.d0 * u(i+2,j,k)     &     
     &                     +   96.d0 * u(i+3,j,k)     &     
     &                     -    9.d0 * u(i+4,j,k)     &     
     &                   ) * idx_by_2520 
          else 
            Dxu(i,j,k) = ( -    9.d0 * u(i-5,j,k)     &      
     &                     +   90.d0 * u(i-4,j,k)     &     
     &                     -  420.d0 * u(i-3,j,k)     &     
     &                     + 1260.d0 * u(i-2,j,k)     &     
     &                     - 3150.d0 * u(i-1,j,k)     &     
     &                     + 1134.d0 * u(i  ,j,k)     &     
     &                     + 1260.d0 * u(i+1,j,k)     &     
     &                     -  180.d0 * u(i+2,j,k)     &     
     &                     +   15.d0 * u(i+3,j,k)     &     
     &                   ) * idx_by_2520 
          endif 

          i=nx-3 
          Dxu(i,j,k) = ( -    9.d0 * u(i-5,j,k)     &      
     &                   +   90.d0 * u(i-4,j,k)     &     
     &                   -  420.d0 * u(i-3,j,k)     &     
     &                   + 1260.d0 * u(i-2,j,k)     &     
     &                   - 3150.d0 * u(i-1,j,k)     &     
     &                   + 1134.d0 * u(i  ,j,k)     &     
     &                   + 1260.d0 * u(i+1,j,k)     &     
     &                   -  180.d0 * u(i+2,j,k)     &     
     &                   +   15.d0 * u(i+3,j,k)     &     
     &                 ) * idx_by_2520  

          i=nx-2 
          Dxu(i,j,k) = (     15.d0 * u(i-6,j,k)     &      
     &                   -  144.d0 * u(i-5,j,k)     &     
     &                   +  630.d0 * u(i-4,j,k)     &     
     &                   - 1680.d0 * u(i-3,j,k)     &     
     &                   + 3150.d0 * u(i-2,j,k)     &     
     &                   - 5040.d0 * u(i-1,j,k)     &     
     &                   + 2394.d0 * u(i  ,j,k)     &     
     &                   +  720.d0 * u(i+1,j,k)     &     
     &                   -   45.d0 * u(i+2,j,k)     &     
     &                 ) * idx_by_2520 

          i=nx-1 
          Dxu(i,j,k) = ( -   45.d0 * u(i-7,j,k)     &      
     &                   +  420.d0 * u(i-6,j,k)     &     
     &                   - 1764.d0 * u(i-5,j,k)     &     
     &                   + 4410.d0 * u(i-4,j,k)     &     
     &                   - 7350.d0 * u(i-3,j,k)     &     
     &                   + 8820.d0 * u(i-2,j,k)     &     
     &                   - 8820.d0 * u(i-1,j,k)     &     
     &                   + 4014.d0 * u(i  ,j,k)     &     
     &                   +  315.d0 * u(i+1,j,k)     &     
     &                 ) * idx_by_2520 

          i=nx 
          Dxu(i,j,k) = (     315.d0 * u(i-8,j,k)     &      
     &                   -  2880.d0 * u(i-7,j,k)     &     
     &                   + 11760.d0 * u(i-6,j,k)     &     
     &                   - 28224.d0 * u(i-5,j,k)     &     
     &                   + 44100.d0 * u(i-4,j,k)     &     
     &                   - 47040.d0 * u(i-3,j,k)     &     
     &                   + 35280.d0 * u(i-2,j,k)     &     
     &                   - 20160.d0 * u(i-1,j,k)     &     
     &                   +  6849.d0 * u(i  ,j,k)     &     
     &                 ) * idx_by_2520 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv8888adv_y calculates the first derivative with respect to 
      !                y using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [6,ny-5].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^y \partial_y in the equations.  
      !                Further, it will depend on the sign of beta^y.  
      !                Namely, if beta^y is positive, we will use a 
      !                fourth order upwind operator while if beta^y 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At j=5, if beta^y is positive, we continue to  
      !                use the eighth order upwind operator but if  
      !                beta^y is negative, we use the eighth order downwind  
      !                operator.   
      !                  
      !                We reverse this at j=ny-4 such that if beta^y is  
      !                positive, we use an eighth order centered operator 
      !                while if beta^y is negative we continue to use  
      !                the eighth order downwind operator.   
      !                 
      !                At j=4, we want an eighth order operator so it
      !                does not matter what the sign of beta^y is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of 
      !                using the downwind operator for the case of 
      !                negative beta^y, but I do not see any way around 
      !                this if we want to have an eighth order operator.  
      !                 
      !                Likewise, at j=ny-3, the requirement of an eighth 
      !                order operator renders the sign of beta^y 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^y is positive. 
      !                 
      !                At j=2,3 and j=ny-2,ny-1, we use differently shifted 
      !                eighth order operators irrespective of beta^ys
      !                sign.  
      !                
      !                At j=1 and j=ny, we use the totally shifted eighth
      !                order operators irrespective of beta^ys sign.  
      !            
      !----------------------------------------------------------------      
      subroutine deriv8888adv_y(Dyu, u, dy, nx, ny, nz, d_type, betay)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dyu(nx,ny,nz), u(nx,ny,nz), betay(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_by_2, idy_by_12, idy_by_60, idy_by_2520  

        idy = 1.0d0/dy        
        idy_by_2 = 0.5d0 * idy
        idy_by_12 = idy / 12.d0
        idy_by_60 = idy / 60.d0
        idy_by_2520 = idy / 2520.d0


        do k = 1, nz
        do i = 1, nx
          j=1
          Dyu(i,j,k) = ( -  6849.d0 * u(i,j  ,k)     &      
     &                   + 20160.d0 * u(i,j+1,k)     &     
     &                   - 35280.d0 * u(i,j+2,k)     &     
     &                   + 47040.d0 * u(i,j+3,k)     &     
     &                   - 44100.d0 * u(i,j+4,k)     &     
     &                   + 28224.d0 * u(i,j+5,k)     &     
     &                   - 11760.d0 * u(i,j+6,k)     &     
     &                   +  2880.d0 * u(i,j+7,k)     &     
     &                   -   315.d0 * u(i,j+8,k)     &     
     &                 ) * idy_by_2520 

          j=2
          Dyu(i,j,k) = ( -  315.d0 * u(i,j-1,k)     &      
     &                   - 4014.d0 * u(i,j  ,k)     &     
     &                   + 8820.d0 * u(i,j+1,k)     &     
     &                   - 8820.d0 * u(i,j+2,k)     &     
     &                   + 7350.d0 * u(i,j+3,k)     &     
     &                   - 4410.d0 * u(i,j+4,k)     &     
     &                   + 1764.d0 * u(i,j+5,k)     &     
     &                   -  420.d0 * u(i,j+6,k)     &     
     &                   +   45.d0 * u(i,j+7,k)     &     
     &                 ) * idy_by_2520 

          j=3
          Dyu(i,j,k) = (     45.d0 * u(i,j-2,k)     &      
     &                   -  720.d0 * u(i,j-1,k)     &     
     &                   - 2394.d0 * u(i,j  ,k)     &     
     &                   + 5040.d0 * u(i,j+1,k)     &     
     &                   - 3150.d0 * u(i,j+2,k)     &     
     &                   + 1680.d0 * u(i,j+3,k)     &     
     &                   -  630.d0 * u(i,j+4,k)     &     
     &                   +  144.d0 * u(i,j+5,k)     &     
     &                   -   15.d0 * u(i,j+6,k)     &     
     &                 ) * idy_by_2520 

          j=4
          Dyu(i,j,k) = ( -   15.d0 * u(i,j-3,k)     &      
     &                   +  180.d0 * u(i,j-2,k)     &     
     &                   - 1260.d0 * u(i,j-1,k)     &     
     &                   - 1134.d0 * u(i,j  ,k)     &     
     &                   + 3150.d0 * u(i,j+1,k)     &     
     &                   - 1260.d0 * u(i,j+2,k)     &     
     &                   +  420.d0 * u(i,j+3,k)     &     
     &                   -   90.d0 * u(i,j+4,k)     &     
     &                   +    9.d0 * u(i,j+5,k)     &     
     &                 ) * idy_by_2520 
          
          j=5  
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = ( -   15.d0 * u(i,j-3,k)     &      
     &                     +  180.d0 * u(i,j-2,k)     &     
     &                     - 1260.d0 * u(i,j-1,k)     &     
     &                     - 1134.d0 * u(i,j  ,k)     &     
     &                     + 3150.d0 * u(i,j+1,k)     &     
     &                     - 1260.d0 * u(i,j+2,k)     &     
     &                     +  420.d0 * u(i,j+3,k)     &     
     &                     -   90.d0 * u(i,j+4,k)     &     
     &                     +    9.d0 * u(i,j+5,k)     &     
     &                   ) * idy_by_2520 
          else 
            Dyu(i,j,k) = (      9.d0 * u(i,j-4,k)     &      
     &                     -   96.d0 * u(i,j-3,k)     &     
     &                     +  504.d0 * u(i,j-2,k)     &     
     &                     - 2016.d0 * u(i,j-1,k)     &     
     &                     + 2016.d0 * u(i,j+1,k)     &     
     &                     -  504.d0 * u(i,j+2,k)     &     
     &                     +   96.d0 * u(i,j+3,k)     &     
     &                     -    9.d0 * u(i,j+4,k)     &     
     &                   ) * idy_by_2520 
          endif 

          do j = 6, ny-5  ! see gr-qc/0505055v2.pdf 
            if ( betay(i,j,k) .ge. 0.d0 ) then
              Dyu(i,j,k) = ( -   15.d0 * u(i,j-3,k)     &      
     &                       +  180.d0 * u(i,j-2,k)     &     
     &                       - 1260.d0 * u(i,j-1,k)     &     
     &                       - 1134.d0 * u(i,j  ,k)     &     
     &                       + 3150.d0 * u(i,j+1,k)     &     
     &                       - 1260.d0 * u(i,j+2,k)     &     
     &                       +  420.d0 * u(i,j+3,k)     &     
     &                       -   90.d0 * u(i,j+4,k)     &     
     &                       +    9.d0 * u(i,j+5,k)     &     
     &                     ) * idy_by_2520 
            else 
              Dyu(i,j,k) = ( -    9.d0 * u(i,j-5,k)     &      
     &                       +   90.d0 * u(i,j-4,k)     &     
     &                       -  420.d0 * u(i,j-3,k)     &     
     &                       + 1260.d0 * u(i,j-2,k)     &     
     &                       - 3150.d0 * u(i,j-1,k)     &     
     &                       + 1134.d0 * u(i,j  ,k)     &     
     &                       + 1260.d0 * u(i,j+1,k)     &     
     &                       -  180.d0 * u(i,j+2,k)     &     
     &                       +   15.d0 * u(i,j+3,k)     &     
     &                     ) * idy_by_2520 
            endif 
          end do

          j=ny-4  
          if ( betay(i,j,k) .ge. 0.d0 ) then
            Dyu(i,j,k) = (      9.d0 * u(i,j-4,k)     &      
     &                     -   96.d0 * u(i,j-3,k)     &     
     &                     +  504.d0 * u(i,j-2,k)     &     
     &                     - 2016.d0 * u(i,j-1,k)     &     
     &                     + 2016.d0 * u(i,j+1,k)     &     
     &                     -  504.d0 * u(i,j+2,k)     &     
     &                     +   96.d0 * u(i,j+3,k)     &     
     &                     -    9.d0 * u(i,j+4,k)     &     
     &                   ) * idy_by_2520 
          else 
            Dyu(i,j,k) = ( -    9.d0 * u(i,j-5,k)     &      
     &                     +   90.d0 * u(i,j-4,k)     &     
     &                     -  420.d0 * u(i,j-3,k)     &     
     &                     + 1260.d0 * u(i,j-2,k)     &     
     &                     - 3150.d0 * u(i,j-1,k)     &     
     &                     + 1134.d0 * u(i,j  ,k)     &     
     &                     + 1260.d0 * u(i,j+1,k)     &     
     &                     -  180.d0 * u(i,j+2,k)     &     
     &                     +   15.d0 * u(i,j+3,k)     &     
     &                   ) * idy_by_2520 
          endif 

          j=ny-3 
          Dyu(i,j,k) = ( -    9.d0 * u(i,j-5,k)     &      
     &                   +   90.d0 * u(i,j-4,k)     &     
     &                   -  420.d0 * u(i,j-3,k)     &     
     &                   + 1260.d0 * u(i,j-2,k)     &     
     &                   - 3150.d0 * u(i,j-1,k)     &     
     &                   + 1134.d0 * u(i,j  ,k)     &     
     &                   + 1260.d0 * u(i,j+1,k)     &     
     &                   -  180.d0 * u(i,j+2,k)     &     
     &                   +   15.d0 * u(i,j+3,k)     &     
     &                 ) * idy_by_2520  

          j=ny-2 
          Dyu(i,j,k) = (     15.d0 * u(i,j-6,k)     &      
     &                   -  144.d0 * u(i,j-5,k)     &     
     &                   +  630.d0 * u(i,j-4,k)     &     
     &                   - 1680.d0 * u(i,j-3,k)     &     
     &                   + 3150.d0 * u(i,j-2,k)     &     
     &                   - 5040.d0 * u(i,j-1,k)     &     
     &                   + 2394.d0 * u(i,j  ,k)     &     
     &                   +  720.d0 * u(i,j+1,k)     &     
     &                   -   45.d0 * u(i,j+2,k)     &     
     &                 ) * idy_by_2520 

          j=ny-1 
          Dyu(i,j,k) = ( -   45.d0 * u(i,j-7,k)     &      
     &                   +  420.d0 * u(i,j-6,k)     &     
     &                   - 1764.d0 * u(i,j-5,k)     &     
     &                   + 4410.d0 * u(i,j-4,k)     &     
     &                   - 7350.d0 * u(i,j-3,k)     &     
     &                   + 8820.d0 * u(i,j-2,k)     &     
     &                   - 8820.d0 * u(i,j-1,k)     &     
     &                   + 4014.d0 * u(i,j  ,k)     &     
     &                   +  315.d0 * u(i,j+1,k)     &     
     &                 ) * idy_by_2520 

          j=ny 
          Dyu(i,j,k) = (     315.d0 * u(i,j-8,k)     &      
     &                   -  2880.d0 * u(i,j-7,k)     &     
     &                   + 11760.d0 * u(i,j-6,k)     &     
     &                   - 28224.d0 * u(i,j-5,k)     &     
     &                   + 44100.d0 * u(i,j-4,k)     &     
     &                   - 47040.d0 * u(i,j-3,k)     &     
     &                   + 35280.d0 * u(i,j-2,k)     &     
     &                   - 20160.d0 * u(i,j-1,k)     &     
     &                   +  6849.d0 * u(i,j  ,k)     &     
     &                 ) * idy_by_2520 
                       
        end do
        end do
       
        return
      end subroutine 
      



      !----------------------------------------------------------------
      !
      ! deriv8888adv_z calculates the first derivative with respect to 
      !                z using an "advective" fourth order finite 
      !                difference operator at all points on the 
      !                interval [6,nz-5].  This derivative, if used, 
      !                will only be applied on "Lie" terms in the equations 
      !                such as beta^z \partial_z in the equations.  
      !                Further, it will depend on the sign of beta^z.  
      !                Namely, if beta^z is positive, we will use a 
      !                fourth order upwind operator while if beta^z 
      !                is negative, we use a fourth order downwind 
      !                operator.  At other points, things get more 
      !                complicated.   
      !                 
      !                At k=5, if beta^z is positive, we continue to  
      !                use the eighth order upwind operator but if  
      !                beta^z is negative, we use the eighth order downwind  
      !                operator.   
      !                  
      !                We reverse this at k=nz-4 such that if beta^z is  
      !                positive, we use an eighth order centered operator 
      !                while if beta^z is negative we continue to use  
      !                the eighth order downwind operator.   
      !                 
      !                At k=4, we want an eighth order operator so it
      !                does not matter what the sign of beta^z is as 
      !                both signs will require an upwind operator.  
      !                This does seem to violate the spirit of 
      !                using the downwind operator for the case of 
      !                negative beta^z, but I do not see any way around 
      !                this if we want to have an eighth order operator.  
      !                 
      !                Likewise, at k=nz-3, the requirement of an eighth 
      !                order operator renders the sign of beta^z 
      !                irrelevant and we are forced to use a downwind 
      !                operator even for the case that beta^z is positive. 
      !                 
      !                At k=2,3 and k=nz-2,nz-1, we use differently shifted 
      !                eighth order operators irrespective of beta^zs
      !                sign.  
      !                
      !                At k=1 and k=nz, we use the totally shifted eighth
      !                order operators irrespective of beta^zs sign.  
      !            
      !----------------------------------------------------------------      
      subroutine deriv8888adv_z(Dzu, u, dz, nx, ny, nz, d_type, betaz)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL Dzu(nx,ny,nz), u(nx,ny,nz), betaz(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_by_2, idz_by_12, idz_by_60, idz_by_2520  

        idz = 1.0d0/dz        
        idz_by_2 = 0.5d0 * idz
        idz_by_12 = idz / 12.d0
        idz_by_60 = idz / 60.d0
        idz_by_2520 = idz / 2520.d0


        do j = 1, ny
        do i = 1, nx
          k=1
          Dzu(i,j,k) = ( -  6849.d0 * u(i,j,k  )     &      
     &                   + 20160.d0 * u(i,j,k+1)     &     
     &                   - 35280.d0 * u(i,j,k+2)     &     
     &                   + 47040.d0 * u(i,j,k+3)     &     
     &                   - 44100.d0 * u(i,j,k+4)     &     
     &                   + 28224.d0 * u(i,j,k+5)     &     
     &                   - 11760.d0 * u(i,j,k+6)     &     
     &                   +  2880.d0 * u(i,j,k+7)     &     
     &                   -   315.d0 * u(i,j,k+8)     &     
     &                 ) * idz_by_2520 

          k=2
          Dzu(i,j,k) = ( -  315.d0 * u(i,j,k-1)     &      
     &                   - 4014.d0 * u(i,j,k  )     &     
     &                   + 8820.d0 * u(i,j,k+1)     &     
     &                   - 8820.d0 * u(i,j,k+2)     &     
     &                   + 7350.d0 * u(i,j,k+3)     &     
     &                   - 4410.d0 * u(i,j,k+4)     &     
     &                   + 1764.d0 * u(i,j,k+5)     &     
     &                   -  420.d0 * u(i,j,k+6)     &     
     &                   +   45.d0 * u(i,j,k+7)     &     
     &                 ) * idz_by_2520 

          k=3
          Dzu(i,j,k) = (     45.d0 * u(i,j,k-2)     &      
     &                   -  720.d0 * u(i,j,k-1)     &     
     &                   - 2394.d0 * u(i,j,k  )     &     
     &                   + 5040.d0 * u(i,j,k+1)     &     
     &                   - 3150.d0 * u(i,j,k+2)     &     
     &                   + 1680.d0 * u(i,j,k+3)     &     
     &                   -  630.d0 * u(i,j,k+4)     &     
     &                   +  144.d0 * u(i,j,k+5)     &     
     &                   -   15.d0 * u(i,j,k+6)     &     
     &                 ) * idz_by_2520 

          k=4
          Dzu(i,j,k) = ( -   15.d0 * u(i,j,k-3)     &      
     &                   +  180.d0 * u(i,j,k-2)     &     
     &                   - 1260.d0 * u(i,j,k-1)     &     
     &                   - 1134.d0 * u(i,j,k  )     &     
     &                   + 3150.d0 * u(i,j,k+1)     &     
     &                   - 1260.d0 * u(i,j,k+2)     &     
     &                   +  420.d0 * u(i,j,k+3)     &     
     &                   -   90.d0 * u(i,j,k+4)     &     
     &                   +    9.d0 * u(i,j,k+5)     &     
     &                 ) * idz_by_2520 
          
          k=5  
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = ( -   15.d0 * u(i,j,k-3)     &      
     &                     +  180.d0 * u(i,j,k-2)     &     
     &                     - 1260.d0 * u(i,j,k-1)     &     
     &                     - 1134.d0 * u(i,j,k  )     &     
     &                     + 3150.d0 * u(i,j,k+1)     &     
     &                     - 1260.d0 * u(i,j,k+2)     &     
     &                     +  420.d0 * u(i,j,k+3)     &     
     &                     -   90.d0 * u(i,j,k+4)     &     
     &                     +    9.d0 * u(i,j,k+5)     &     
     &                   ) * idz_by_2520 
          else 
            Dzu(i,j,k) = (      9.d0 * u(i,j,k-4)     &      
     &                     -   96.d0 * u(i,j,k-3)     &     
     &                     +  504.d0 * u(i,j,k-2)     &     
     &                     - 2016.d0 * u(i,j,k-1)     &     
     &                     + 2016.d0 * u(i,j,k+1)     &     
     &                     -  504.d0 * u(i,j,k+2)     &     
     &                     +   96.d0 * u(i,j,k+3)     &     
     &                     -    9.d0 * u(i,j,k+4)     &     
     &                   ) * idz_by_2520 
          endif 

          do k = 6, nz-5  ! see gr-qc/0505055v2.pdf 
            if ( betaz(i,j,k) .ge. 0.d0 ) then
              Dzu(i,j,k) = ( -   15.d0 * u(i,j,k-3)     &      
     &                       +  180.d0 * u(i,j,k-2)     &     
     &                       - 1260.d0 * u(i,j,k-1)     &     
     &                       - 1134.d0 * u(i,j,k  )     &     
     &                       + 3150.d0 * u(i,j,k+1)     &     
     &                       - 1260.d0 * u(i,j,k+2)     &     
     &                       +  420.d0 * u(i,j,k+3)     &     
     &                       -   90.d0 * u(i,j,k+4)     &     
     &                       +    9.d0 * u(i,j,k+5)     &     
     &                     ) * idz_by_2520 
            else 
              Dzu(i,j,k) = ( -    9.d0 * u(i,j,k-5)     &      
     &                       +   90.d0 * u(i,j,k-4)     &     
     &                       -  420.d0 * u(i,j,k-3)     &     
     &                       + 1260.d0 * u(i,j,k-2)     &     
     &                       - 3150.d0 * u(i,j,k-1)     &     
     &                       + 1134.d0 * u(i,j,k  )     &     
     &                       + 1260.d0 * u(i,j,k+1)     &     
     &                       -  180.d0 * u(i,j,k+2)     &     
     &                       +   15.d0 * u(i,j,k+3)     &     
     &                     ) * idz_by_2520 
            endif 
          end do

          k=nz-4  
          if ( betaz(i,j,k) .ge. 0.d0 ) then
            Dzu(i,j,k) = (      9.d0 * u(i,j,k-4)     &      
     &                     -   96.d0 * u(i,j,k-3)     &     
     &                     +  504.d0 * u(i,j,k-2)     &     
     &                     - 2016.d0 * u(i,j,k-1)     &     
     &                     + 2016.d0 * u(i,j,k+1)     &     
     &                     -  504.d0 * u(i,j,k+2)     &     
     &                     +   96.d0 * u(i,j,k+3)     &     
     &                     -    9.d0 * u(i,j,k+4)     &     
     &                   ) * idz_by_2520 
          else 
            Dzu(i,j,k) = ( -    9.d0 * u(i,j,k-5)     &      
     &                     +   90.d0 * u(i,j,k-4)     &     
     &                     -  420.d0 * u(i,j,k-3)     &     
     &                     + 1260.d0 * u(i,j,k-2)     &     
     &                     - 3150.d0 * u(i,j,k-1)     &     
     &                     + 1134.d0 * u(i,j,k  )     &     
     &                     + 1260.d0 * u(i,j,k+1)     &     
     &                     -  180.d0 * u(i,j,k+2)     &     
     &                     +   15.d0 * u(i,j,k+3)     &     
     &                   ) * idz_by_2520 
          endif 

          k=nz-3 
          Dzu(i,j,k) = ( -    9.d0 * u(i,j,k-5)     &      
     &                   +   90.d0 * u(i,j,k-4)     &     
     &                   -  420.d0 * u(i,j,k-3)     &     
     &                   + 1260.d0 * u(i,j,k-2)     &     
     &                   - 3150.d0 * u(i,j,k-1)     &     
     &                   + 1134.d0 * u(i,j,k  )     &     
     &                   + 1260.d0 * u(i,j,k+1)     &     
     &                   -  180.d0 * u(i,j,k+2)     &     
     &                   +   15.d0 * u(i,j,k+3)     &     
     &                 ) * idz_by_2520  

          k=nz-2 
          Dzu(i,j,k) = (     15.d0 * u(i,j,k-6)     &      
     &                   -  144.d0 * u(i,j,k-5)     &     
     &                   +  630.d0 * u(i,j,k-4)     &     
     &                   - 1680.d0 * u(i,j,k-3)     &     
     &                   + 3150.d0 * u(i,j,k-2)     &     
     &                   - 5040.d0 * u(i,j,k-1)     &     
     &                   + 2394.d0 * u(i,j,k  )     &     
     &                   +  720.d0 * u(i,j,k+1)     &     
     &                   -   45.d0 * u(i,j,k+2)     &     
     &                 ) * idz_by_2520 

          k=nz-1 
          Dzu(i,j,k) = ( -   45.d0 * u(i,j,k-7)     &      
     &                   +  420.d0 * u(i,j,k-6)     &     
     &                   - 1764.d0 * u(i,j,k-5)     &     
     &                   + 4410.d0 * u(i,j,k-4)     &     
     &                   - 7350.d0 * u(i,j,k-3)     &     
     &                   + 8820.d0 * u(i,j,k-2)     &     
     &                   - 8820.d0 * u(i,j,k-1)     &     
     &                   + 4014.d0 * u(i,j,k  )     &     
     &                   +  315.d0 * u(i,j,k+1)     &     
     &                 ) * idz_by_2520 

          k=nz 
          Dzu(i,j,k) = (     315.d0 * u(i,j,k-8)     &      
     &                   -  2880.d0 * u(i,j,k-7)     &     
     &                   + 11760.d0 * u(i,j,k-6)     &     
     &                   - 28224.d0 * u(i,j,k-5)     &     
     &                   + 44100.d0 * u(i,j,k-4)     &     
     &                   - 47040.d0 * u(i,j,k-3)     &     
     &                   + 35280.d0 * u(i,j,k-2)     &     
     &                   - 20160.d0 * u(i,j,k-1)     &     
     &                   +  6849.d0 * u(i,j,k  )     &     
     &                 ) * idz_by_2520 
                       
        end do
        end do
       
        return
      end subroutine 
      





      !----------------------------------------------------------------
      !
      ! deriv8642_xx calculates the second derivative with respect to  
      !              x using a centered eighth order finite difference 
      !              operator on the interval [5,nx-4].  At other 
      !              points, we attempt to use centered difference 
      !              operators.  
      !   
      !              At i=4 and i=nx-3, we use centered sixth order 
      !              operators.  
      ! 
      !              At i=3 and i=nx-2, we use centered fourth order 
      !              operators.  
      ! 
      !              At i=2 and i=nx-1, we use centered second order 
      !              operators.  
      ! 
      !              At the boundary points (i=1 and i=nx) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at i=1, we use a forward (upwind)  
      !              operator with a four point stencil.  At i=nx, we  
      !              use a backward (downwind) operator, again with 
      !              a four point stencil. 
      !
      !----------------------------------------------------------------
      subroutine deriv8642_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL DxDxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_sqrd, idx_sqrd_by_12 
        CCTK_REAL idx_sqrd_by_180, idx_sqrd_by_5040  

        idx_sqrd         = 1.0d0 /(dx*dx)         
        idx_sqrd_by_12   = idx_sqrd / 12.d0
        idx_sqrd_by_180  = idx_sqrd / 180.d0
        idx_sqrd_by_5040 = idx_sqrd / 5040.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          DxDxu(i,j,k) = (    2.d0 * u(i  ,j,k)     &      
     &                     -  5.d0 * u(i+1,j,k)     &     
     &                     +  4.d0 * u(i+2,j,k)     &     
     &                     -         u(i+3,j,k)     &     
     &                   ) * idx_sqrd 
          i=2
          DxDxu(i,j,k) = (           u(i-1,j,k)     &      
     &                     -  2.d0 * u(i  ,j,k)     &     
     &                     +         u(i+1,j,k)     &     
     &                   ) * idx_sqrd
          i=3
          DxDxu(i,j,k) = ( -         u(i-2,j,k)     &      
     &                     + 16.d0 * u(i-1,j,k)     &     
     &                     - 30.d0 * u(i  ,j,k)     &     
     &                     + 16.d0 * u(i+1,j,k)     &      
     &                     -         u(i+2,j,k)     &      
     &                   ) * idx_sqrd_by_12 
          i=4
          DxDxu(i,j,k) = (     2.d0 * u(i-3,j,k)     &      
     &                     -  27.d0 * u(i-2,j,k)     &     
     &                     + 270.d0 * u(i-1,j,k)     &     
     &                     - 490.d0 * u(i  ,j,k)     &     
     &                     + 270.d0 * u(i+1,j,k)     &     
     &                     -  27.d0 * u(i+2,j,k)     &     
     &                     +   2.d0 * u(i+3,j,k)     &     
     &                   ) * idx_sqrd_by_180
          do i = 5, nx-4
            DxDxu(i,j,k) = ( -     9.d0 * u(i-4,j,k)     &      
     &                       +   128.d0 * u(i-3,j,k)     &     
     &                       -  1008.d0 * u(i-2,j,k)     &     
     &                       +  8064.d0 * u(i-1,j,k)     &     
     &                       - 14350.d0 * u(i  ,j,k)     &     
     &                       +  8064.d0 * u(i+1,j,k)     &     
     &                       -  1008.d0 * u(i+2,j,k)     &     
     &                       +   128.d0 * u(i+3,j,k)     &     
     &                       -     9.d0 * u(i+4,j,k)     &     
     &                     ) * idx_sqrd_by_5040 
          end do
          i=nx-3
          DxDxu(i,j,k) = (     2.d0 * u(i-3,j,k)     &      
     &                     -  27.d0 * u(i-2,j,k)     &     
     &                     + 270.d0 * u(i-1,j,k)     &     
     &                     - 490.d0 * u(i  ,j,k)     &     
     &                     + 270.d0 * u(i+1,j,k)     &     
     &                     -  27.d0 * u(i+2,j,k)     &     
     &                     +   2.d0 * u(i+3,j,k)     &     
     &                   ) * idx_sqrd_by_180 
          i=nx-2 
          DxDxu(i,j,k) = ( -         u(i-2,j,k)     &      
     &                     + 16.d0 * u(i-1,j,k)     &     
     &                     - 30.d0 * u(i  ,j,k)     &     
     &                     + 16.d0 * u(i+1,j,k)     &      
     &                     -         u(i+2,j,k)     &      
     &                   ) * idx_sqrd_by_12
          i=nx-1 
          DxDxu(i,j,k) = ( +        u(i-1,j,k)     &      
     &                     - 2.d0 * u(i  ,j,k)     &      
     &                     +        u(i+1,j,k)     &      
     &                   ) * idx_sqrd 
          i=nx 
          DxDxu(i,j,k) = ( -         u(i-3,j,k)     &      
     &                     +  4.d0 * u(i-2,j,k)     &     
     &                     -  5.d0 * u(i-1,j,k)     &     
     &                     +  2.d0 * u(i  ,j,k)     &     
     &                   ) * idx_sqrd
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      ! deriv8642_yy calculates the second derivative with respect to  
      !              y using a centered eighth order finite difference 
      !              operator on the interval [5,ny-4].  At other 
      !              points, we attempt to use centered difference 
      !              operators.  
      !   
      !              At j=4 and j=ny-3, we use centered sixth order 
      !              operators.  
      ! 
      !              At j=3 and j=ny-2, we use centered fourth order 
      !              operators.  
      ! 
      !              At j=2 and j=ny-1, we use centered second order 
      !              operators.  
      ! 
      !              At the boundary points (j=1 and j=ny) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at j=1, we use a forward (upwind)  
      !              operator with a four point stencil.  At j=ny, we  
      !              use a backward (downwind) operator, again with 
      !              a four point stencil. 
      !
      !----------------------------------------------------------------
      subroutine deriv8642_yy(DyDyu, u, dy, nx, ny, nz, dd_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_sqrd, idy_sqrd_by_12
        CCTK_REAL idy_sqrd_by_180, idy_sqrd_by_5040  

        idy_sqrd         = 1.0d0 /(dy*dy)         
        idy_sqrd_by_12   = idy_sqrd / 12.d0
        idy_sqrd_by_180  = idy_sqrd / 180.d0
        idy_sqrd_by_5040 = idy_sqrd / 5040.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          DyDyu(i,j,k) = (    2.d0 * u(i,j  ,k)     &      
     &                     -  5.d0 * u(i,j+1,k)     &     
     &                     +  4.d0 * u(i,j+2,k)     &     
     &                     -         u(i,j+3,k)     &     
     &                   ) * idy_sqrd 
          j=2
          DyDyu(i,j,k) = (           u(i,j-1,k)     &      
     &                     -  2.d0 * u(i,j  ,k)     &     
     &                     +         u(i,j+1,k)     &     
     &                   ) * idy_sqrd
          j=3
          DyDyu(i,j,k) = ( -         u(i,j-2,k)     &      
     &                     + 16.d0 * u(i,j-1,k)     &     
     &                     - 30.d0 * u(i,j  ,k)     &     
     &                     + 16.d0 * u(i,j+1,k)     &      
     &                     -         u(i,j+2,k)     &      
     &                   ) * idy_sqrd_by_12 
          j=4
          DyDyu(i,j,k) = (     2.d0 * u(i,j-3,k)     &      
     &                     -  27.d0 * u(i,j-2,k)     &     
     &                     + 270.d0 * u(i,j-1,k)     &     
     &                     - 490.d0 * u(i,j  ,k)     &     
     &                     + 270.d0 * u(i,j+1,k)     &     
     &                     -  27.d0 * u(i,j+2,k)     &     
     &                     +   2.d0 * u(i,j+3,k)     &     
     &                   ) * idy_sqrd_by_180
          do j = 5, ny-4
            DyDyu(i,j,k) = ( -     9.d0 * u(i,j-4,k)     &      
     &                       +   128.d0 * u(i,j-3,k)     &     
     &                       -  1008.d0 * u(i,j-2,k)     &     
     &                       +  8064.d0 * u(i,j-1,k)     &     
     &                       - 14350.d0 * u(i,j  ,k)     &     
     &                       +  8064.d0 * u(i,j+1,k)     &     
     &                       -  1008.d0 * u(i,j+2,k)     &     
     &                       +   128.d0 * u(i,j+3,k)     &     
     &                       -     9.d0 * u(i,j+4,k)     &     
     &                     ) * idy_sqrd_by_5040 
          end do
          j=ny-3
          DyDyu(i,j,k) = (     2.d0 * u(i,j-3,k)     &      
     &                     -  27.d0 * u(i,j-2,k)     &     
     &                     + 270.d0 * u(i,j-1,k)     &     
     &                     - 490.d0 * u(i,j  ,k)     &     
     &                     + 270.d0 * u(i,j+1,k)     &     
     &                     -  27.d0 * u(i,j+2,k)     &     
     &                     +   2.d0 * u(i,j+3,k)     &     
     &                   ) * idy_sqrd_by_180 
          j=ny-2 
          DyDyu(i,j,k) = ( -         u(i,j-2,k)     &      
     &                     + 16.d0 * u(i,j-1,k)     &     
     &                     - 30.d0 * u(i,j  ,k)     &     
     &                     + 16.d0 * u(i,j+1,k)     &      
     &                     -         u(i,j+2,k)     &      
     &                   ) * idy_sqrd_by_12
          j=ny-1 
          DyDyu(i,j,k) = ( +        u(i,j-1,k)     &      
     &                     - 2.d0 * u(i,j  ,k)     &      
     &                     +        u(i,j+1,k)     &      
     &                   ) * idy_sqrd 
          j=ny 
          DyDyu(i,j,k) = ( -         u(i,j-3,k)     &      
     &                     +  4.d0 * u(i,j-2,k)     &     
     &                     -  5.d0 * u(i,j-1,k)     &     
     &                     +  2.d0 * u(i,j  ,k)     &     
     &                   ) * idy_sqrd
        end do
        end do
       
        return
      end subroutine 



      !----------------------------------------------------------------
      !
      ! deriv8642_zz calculates the second derivative with respect to  
      !              z using a centered eighth order finite difference 
      !              operator on the interval [5,nz-4].  At other 
      !              points, we attempt to use centered difference 
      !              operators.  
      !   
      !              At k=4 and k=nz-3, we use centered sixth order 
      !              operators.  
      ! 
      !              At k=3 and k=nz-2, we use centered fourth order 
      !              operators.  
      ! 
      !              At k=2 and k=nz-1, we use centered second order 
      !              operators.  
      ! 
      !              At the boundary points (k=1 and k=nz) we can not,  
      !              of course, define centered difference operators.  
      !              As a result, we have to spoil our effort in this  
      !              routine to use only centered difference operators.  
      !              Instead, at the boundary points, we choose to use 
      !              appropriately shifted second order operators.  In
      !              particular, at k=1, we use a forward (upwind)  
      !              operator with a four point stencil.  At k=nz, we  
      !              use a backward (downwind) operator, again with 
      !              a four point stencil. 
      !
      !----------------------------------------------------------------
      subroutine deriv8642_zz(DzDzu, u, dz, nx, ny, nz, dd_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL DzDzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_sqrd, idz_sqrd_by_12
        CCTK_REAL idz_sqrd_by_180, idz_sqrd_by_5040  

        idz_sqrd         = 1.0d0 /(dz*dz)         
        idz_sqrd_by_12   = idz_sqrd / 12.d0
        idz_sqrd_by_180  = idz_sqrd / 180.d0
        idz_sqrd_by_5040 = idz_sqrd / 5040.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          DzDzu(i,j,k) = (    2.d0 * u(i,j,k  )     &      
     &                     -  5.d0 * u(i,j,k+1)     &     
     &                     +  4.d0 * u(i,j,k+2)     &     
     &                     -         u(i,j,k+3)     &     
     &                   ) * idz_sqrd 
          k=2
          DzDzu(i,j,k) = (           u(i,j,k-1)     &      
     &                     -  2.d0 * u(i,j,k  )     &     
     &                     +         u(i,j,k+1)     &     
     &                   ) * idz_sqrd
          k=3
          DzDzu(i,j,k) = ( -         u(i,j,k-2)     &      
     &                     + 16.d0 * u(i,j,k-1)     &     
     &                     - 30.d0 * u(i,j,k  )     &     
     &                     + 16.d0 * u(i,j,k+1)     &      
     &                     -         u(i,j,k+2)     &      
     &                   ) * idz_sqrd_by_12 
          k=4
          DzDzu(i,j,k) = (     2.d0 * u(i,j,k-3)     &      
     &                     -  27.d0 * u(i,j,k-2)     &     
     &                     + 270.d0 * u(i,j,k-1)     &     
     &                     - 490.d0 * u(i,j,k  )     &     
     &                     + 270.d0 * u(i,j,k+1)     &     
     &                     -  27.d0 * u(i,j,k+2)     &     
     &                     +   2.d0 * u(i,j,k+3)     &     
     &                   ) * idz_sqrd_by_180
          do k = 5, nz-4
            DzDzu(i,j,k) = ( -     9.d0 * u(i,j,k-4)     &      
     &                       +   128.d0 * u(i,j,k-3)     &     
     &                       -  1008.d0 * u(i,j,k-2)     &     
     &                       +  8064.d0 * u(i,j,k-1)     &     
     &                       - 14350.d0 * u(i,j,k  )     &     
     &                       +  8064.d0 * u(i,j,k+1)     &     
     &                       -  1008.d0 * u(i,j,k+2)     &     
     &                       +   128.d0 * u(i,j,k+3)     &     
     &                       -     9.d0 * u(i,j,k+4)     &     
     &                     ) * idz_sqrd_by_5040 
          end do
          k=nz-3
          DzDzu(i,j,k) = (     2.d0 * u(i,j,k-3)     &      
     &                     -  27.d0 * u(i,j,k-2)     &     
     &                     + 270.d0 * u(i,j,k-1)     &     
     &                     - 490.d0 * u(i,j,k  )     &     
     &                     + 270.d0 * u(i,j,k+1)     &     
     &                     -  27.d0 * u(i,j,k+2)     &     
     &                     +   2.d0 * u(i,j,k+3)     &     
     &                   ) * idz_sqrd_by_180 
          k=nz-2 
          DzDzu(i,j,k) = ( -         u(i,j,k-2)     &      
     &                     + 16.d0 * u(i,j,k-1)     &     
     &                     - 30.d0 * u(i,j,k  )     &     
     &                     + 16.d0 * u(i,j,k+1)     &      
     &                     -         u(i,j,k+2)     &      
     &                   ) * idz_sqrd_by_12
          k=nz-1 
          DzDzu(i,j,k) = ( +        u(i,j,k-1)     &      
     &                     - 2.d0 * u(i,j,k  )     &      
     &                     +        u(i,j,k+1)     &      
     &                   ) * idz_sqrd 
          k=nz 
          DzDzu(i,j,k) = ( -         u(i,j,k-3)     &      
     &                     +  4.d0 * u(i,j,k-2)     &     
     &                     -  5.d0 * u(i,j,k-1)     &     
     &                     +  2.d0 * u(i,j,k  )     &     
     &                   ) * idz_sqrd
        end do
        end do
       
        return
      end subroutine 






      !----------------------------------------------------------------
      !
      ! deriv8888_xx calculates the second derivative with respect to  
      !              x using a centered eighth order finite difference 
      !              operator at all points on the interval [5,nx-4].  
      ! 
      !              At i=4 and i=nx-3, we use once-shifted eighth order 
      !              operators.  This has a ten point stencil.  
      ! 
      !              At i=3 and i=nx-2, we use twice-shifted eighth order 
      !              operators.  This has a ten point stencil.  
      ! 
      !              At i=2 and i=nx-1, we use thrice-shifted eighth order
      !              operators.  This has a ten point stencil.  
      ! 
      !              At i=1 and i=nx, we use totally shifted eighth order
      !              operators.  This has a ten point stencil.  
      ! 
      !----------------------------------------------------------------
      subroutine deriv8888_xx(DxDxu, u, dx, nx, ny, nz, dd_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL DxDxu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dx
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idx, idx_sqrd, idx_sqrd_by_12 
        CCTK_REAL idx_sqrd_by_180, idx_sqrd_by_5040  

        idx_sqrd         = 1.0d0 /(dx*dx)         
        idx_sqrd_by_12   = idx_sqrd / 12.d0
        idx_sqrd_by_180  = idx_sqrd / 180.d0
        idx_sqrd_by_5040 = idx_sqrd / 5040.d0

        do k = 1, nz
        do j = 1, ny
          i=1
          DxDxu(i,j,k) = (    32575.d0 * u(i  ,j,k)     &      
     &                     - 165924.d0 * u(i+1,j,k)     &     
     &                     + 422568.d0 * u(i+2,j,k)     &     
     &                     - 704368.d0 * u(i+3,j,k)     &     
     &                     + 818874.d0 * u(i+4,j,k)     &     
     &                     - 667800.d0 * u(i+5,j,k)     &     
     &                     + 375704.d0 * u(i+6,j,k)     &     
     &                     - 139248.d0 * u(i+7,j,k)     &     
     &                     +  30663.d0 * u(i+8,j,k)     &     
     &                     -   3044.d0 * u(i+9,j,k)     &     
     &                   ) * idx_sqrd_by_5040  
          i=2
          DxDxu(i,j,k) = (    3044.d0 * u(i-1,j,k)     &      
     &                     +  2135.d0 * u(i  ,j,k)     &     
     &                     - 28944.d0 * u(i+1,j,k)     &     
     &                     + 57288.d0 * u(i+2,j,k)     &     
     &                     - 65128.d0 * u(i+3,j,k)     &     
     &                     + 51786.d0 * u(i+4,j,k)     &     
     &                     - 28560.d0 * u(i+5,j,k)     &     
     &                     + 10424.d0 * u(i+6,j,k)     &     
     &                     -  2268.d0 * u(i+7,j,k)     &     
     &                     +   223.d0 * u(i+8,j,k)     &     
     &                   ) * idx_sqrd_by_5040 
          i=3
          DxDxu(i,j,k) = ( -   223.d0 * u(i-2,j,k)     &      
     &                     +  5274.d0 * u(i-1,j,k)     &     
     &                     -  7900.d0 * u(i  ,j,k)     &     
     &                     -  2184.d0 * u(i+1,j,k)     &      
     &                     + 10458.d0 * u(i+2,j,k)     &      
     &                     -  8932.d0 * u(i+3,j,k)     &      
     &                     +  4956.d0 * u(i+4,j,k)     &      
     &                     -  1800.d0 * u(i+5,j,k)     &      
     &                     +   389.d0 * u(i+6,j,k)     &      
     &                     -    38.d0 * u(i+7,j,k)     &      
     &                   ) * idx_sqrd_by_5040  
          i=4
          DxDxu(i,j,k) = (      38.d0 * u(i-3,j,k)     &      
     &                     -   603.d0 * u(i-2,j,k)     &     
     &                     +  6984.d0 * u(i-1,j,k)     &     
     &                     - 12460.d0 * u(i  ,j,k)     &     
     &                     +  5796.d0 * u(i+1,j,k)     &     
     &                     +   882.d0 * u(i+2,j,k)     &     
     &                     -   952.d0 * u(i+3,j,k)     &     
     &                     +   396.d0 * u(i+4,j,k)     &     
     &                     -    90.d0 * u(i+5,j,k)     &     
     &                     +     9.d0 * u(i+6,j,k)     &     
     &                   ) * idx_sqrd_by_5040
          do i = 5, nx-4
            DxDxu(i,j,k) = ( -     9.d0 * u(i-4,j,k)     &      
     &                       +   128.d0 * u(i-3,j,k)     &     
     &                       -  1008.d0 * u(i-2,j,k)     &     
     &                       +  8064.d0 * u(i-1,j,k)     &     
     &                       - 14350.d0 * u(i  ,j,k)     &     
     &                       +  8064.d0 * u(i+1,j,k)     &     
     &                       -  1008.d0 * u(i+2,j,k)     &     
     &                       +   128.d0 * u(i+3,j,k)     &     
     &                       -     9.d0 * u(i+4,j,k)     &     
     &                     ) * idx_sqrd_by_5040 
          end do
          i=nx-3
          DxDxu(i,j,k) = (       9.d0 * u(i-6,j,k)     &      
     &                     -    90.d0 * u(i-5,j,k)     &     
     &                     +   396.d0 * u(i-4,j,k)     &     
     &                     -   952.d0 * u(i-3,j,k)     &     
     &                     +   882.d0 * u(i-2,j,k)     &     
     &                     +  5796.d0 * u(i-1,j,k)     &     
     &                     - 12460.d0 * u(i  ,j,k)     &     
     &                     +  6984.d0 * u(i+1,j,k)     &     
     &                     -   603.d0 * u(i+2,j,k)     &     
     &                     +    38.d0 * u(i+3,j,k)     &     
     &                   ) * idx_sqrd_by_5040 
          i=nx-2 
          DxDxu(i,j,k) = ( -    38.d0 * u(i-7,j,k)     &      
     &                     +   389.d0 * u(i-6,j,k)     &     
     &                     -  1800.d0 * u(i-5,j,k)     &     
     &                     +  4956.d0 * u(i-4,j,k)     &     
     &                     -  8932.d0 * u(i-3,j,k)     &     
     &                     + 10458.d0 * u(i-2,j,k)     &     
     &                     -  2184.d0 * u(i-1,j,k)     &     
     &                     -  7900.d0 * u(i  ,j,k)     &     
     &                     +  5274.d0 * u(i+1,j,k)     &      
     &                     -   223.d0 * u(i+2,j,k)     &      
     &                   ) * idx_sqrd_by_5040 
          i=nx-1 
          DxDxu(i,j,k) = (     223.d0 * u(i-8,j,k)     &      
     &                     -  2268.d0 * u(i-7,j,k)     &      
     &                     + 10424.d0 * u(i-6,j,k)     &      
     &                     - 28560.d0 * u(i-5,j,k)     &      
     &                     + 51786.d0 * u(i-4,j,k)     &      
     &                     - 65128.d0 * u(i-3,j,k)     &      
     &                     + 57288.d0 * u(i-2,j,k)     &      
     &                     - 28944.d0 * u(i-1,j,k)     &      
     &                     +  2135.d0 * u(i  ,j,k)     &      
     &                     +  3044.d0 * u(i+1,j,k)     &      
     &                   ) * idx_sqrd_by_5040  
          i=nx 
          DxDxu(i,j,k) = ( -   3044.d0 * u(i-9,j,k)     &      
     &                     +  30663.d0 * u(i-8,j,k)     &     
     &                     - 139248.d0 * u(i-7,j,k)     &     
     &                     + 375704.d0 * u(i-6,j,k)     &     
     &                     - 667800.d0 * u(i-5,j,k)     &     
     &                     + 818874.d0 * u(i-4,j,k)     &     
     &                     - 704368.d0 * u(i-3,j,k)     &     
     &                     + 422568.d0 * u(i-2,j,k)     &     
     &                     - 165924.d0 * u(i-1,j,k)     &     
     &                     +  32575.d0 * u(i  ,j,k)     &     
     &                   ) * idx_sqrd_by_5040 
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      ! deriv8888_yy calculates the second derivative with respect to  
      !              y using a centered eighth order finite difference 
      !              operator at all points on the interval [5,ny-4].  
      ! 
      !              At j=4 and j=ny-3, we use once-shifted eighth order 
      !              operators.  This has a ten point stencil.  
      ! 
      !              At j=3 and j=ny-2, we use twice-shifted eighth order 
      !              operators.  This has a ten point stencil.  
      ! 
      !              At j=2 and j=ny-1, we use thrice-shifted eighth order
      !              operators.  This has a ten point stencil.  
      ! 
      !              At j=1 and j=ny, we use totally shifted eighth order
      !              operators.  This has a ten point stencil.  
      ! 
      !----------------------------------------------------------------
      subroutine deriv8888_yy(DyDyu, u, dy, nx, ny, nz, dd_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL DyDyu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dy
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idy, idy_sqrd, idy_sqrd_by_12 
        CCTK_REAL idy_sqrd_by_180, idy_sqrd_by_5040  

        idy_sqrd         = 1.0d0 /(dy*dy)         
        idy_sqrd_by_12   = idy_sqrd / 12.d0
        idy_sqrd_by_180  = idy_sqrd / 180.d0
        idy_sqrd_by_5040 = idy_sqrd / 5040.d0

        do k = 1, nz
        do i = 1, nx
          j=1
          DyDyu(i,j,k) = (    32575.d0 * u(i,j  ,k)     &      
     &                     - 165924.d0 * u(i,j+1,k)     &     
     &                     + 422568.d0 * u(i,j+2,k)     &     
     &                     - 704368.d0 * u(i,j+3,k)     &     
     &                     + 818874.d0 * u(i,j+4,k)     &     
     &                     - 667800.d0 * u(i,j+5,k)     &     
     &                     + 375704.d0 * u(i,j+6,k)     &     
     &                     - 139248.d0 * u(i,j+7,k)     &     
     &                     +  30663.d0 * u(i,j+8,k)     &     
     &                     -   3044.d0 * u(i,j+9,k)     &     
     &                   ) * idy_sqrd_by_5040  
          j=2
          DyDyu(i,j,k) = (    3044.d0 * u(i,j-1,k)     &      
     &                     +  2135.d0 * u(i,j  ,k)     &     
     &                     - 28944.d0 * u(i,j+1,k)     &     
     &                     + 57288.d0 * u(i,j+2,k)     &     
     &                     - 65128.d0 * u(i,j+3,k)     &     
     &                     + 51786.d0 * u(i,j+4,k)     &     
     &                     - 28560.d0 * u(i,j+5,k)     &     
     &                     + 10424.d0 * u(i,j+6,k)     &     
     &                     -  2268.d0 * u(i,j+7,k)     &     
     &                     +   223.d0 * u(i,j+8,k)     &     
     &                   ) * idy_sqrd_by_5040 
          j=3
          DyDyu(i,j,k) = ( -   223.d0 * u(i,j-2,k)     &      
     &                     +  5274.d0 * u(i,j-1,k)     &     
     &                     -  7900.d0 * u(i,j  ,k)     &     
     &                     -  2184.d0 * u(i,j+1,k)     &      
     &                     + 10458.d0 * u(i,j+2,k)     &      
     &                     -  8932.d0 * u(i,j+3,k)     &      
     &                     +  4956.d0 * u(i,j+4,k)     &      
     &                     -  1800.d0 * u(i,j+5,k)     &      
     &                     +   389.d0 * u(i,j+6,k)     &      
     &                     -    38.d0 * u(i,j+7,k)     &      
     &                   ) * idy_sqrd_by_5040  
          j=4
          DyDyu(i,j,k) = (      38.d0 * u(i,j-3,k)     &      
     &                     -   603.d0 * u(i,j-2,k)     &     
     &                     +  6984.d0 * u(i,j-1,k)     &     
     &                     - 12460.d0 * u(i,j  ,k)     &     
     &                     +  5796.d0 * u(i,j+1,k)     &     
     &                     +   882.d0 * u(i,j+2,k)     &     
     &                     -   952.d0 * u(i,j+3,k)     &     
     &                     +   396.d0 * u(i,j+4,k)     &     
     &                     -    90.d0 * u(i,j+5,k)     &     
     &                     +     9.d0 * u(i,j+6,k)     &     
     &                   ) * idy_sqrd_by_5040
          do j = 5, ny-4
            DyDyu(i,j,k) = ( -     9.d0 * u(i,j-4,k)     &      
     &                       +   128.d0 * u(i,j-3,k)     &     
     &                       -  1008.d0 * u(i,j-2,k)     &     
     &                       +  8064.d0 * u(i,j-1,k)     &     
     &                       - 14350.d0 * u(i,j  ,k)     &     
     &                       +  8064.d0 * u(i,j+1,k)     &     
     &                       -  1008.d0 * u(i,j+2,k)     &     
     &                       +   128.d0 * u(i,j+3,k)     &     
     &                       -     9.d0 * u(i,j+4,k)     &     
     &                     ) * idy_sqrd_by_5040 
          end do
          j=ny-3
          DyDyu(i,j,k) = (       9.d0 * u(i,j-6,k)     &      
     &                     -    90.d0 * u(i,j-5,k)     &     
     &                     +   396.d0 * u(i,j-4,k)     &     
     &                     -   952.d0 * u(i,j-3,k)     &     
     &                     +   882.d0 * u(i,j-2,k)     &     
     &                     +  5796.d0 * u(i,j-1,k)     &     
     &                     - 12460.d0 * u(i,j  ,k)     &     
     &                     +  6984.d0 * u(i,j+1,k)     &     
     &                     -   603.d0 * u(i,j+2,k)     &     
     &                     +    38.d0 * u(i,j+3,k)     &     
     &                   ) * idy_sqrd_by_5040 
          j=ny-2 
          DyDyu(i,j,k) = ( -    38.d0 * u(i,j-7,k)     &      
     &                     +   389.d0 * u(i,j-6,k)     &     
     &                     -  1800.d0 * u(i,j-5,k)     &     
     &                     +  4956.d0 * u(i,j-4,k)     &     
     &                     -  8932.d0 * u(i,j-3,k)     &     
     &                     + 10458.d0 * u(i,j-2,k)     &     
     &                     -  2184.d0 * u(i,j-1,k)     &     
     &                     -  7900.d0 * u(i,j  ,k)     &     
     &                     +  5274.d0 * u(i,j+1,k)     &      
     &                     -   223.d0 * u(i,j+2,k)     &      
     &                   ) * idy_sqrd_by_5040 
          j=ny-1 
          DyDyu(i,j,k) = (     223.d0 * u(i,j-8,k)     &      
     &                     -  2268.d0 * u(i,j-7,k)     &      
     &                     + 10424.d0 * u(i,j-6,k)     &      
     &                     - 28560.d0 * u(i,j-5,k)     &      
     &                     + 51786.d0 * u(i,j-4,k)     &      
     &                     - 65128.d0 * u(i,j-3,k)     &      
     &                     + 57288.d0 * u(i,j-2,k)     &      
     &                     - 28944.d0 * u(i,j-1,k)     &      
     &                     +  2135.d0 * u(i,j  ,k)     &      
     &                     +  3044.d0 * u(i,j+1,k)     &      
     &                   ) * idy_sqrd_by_5040  
          j=ny 
          DyDyu(i,j,k) = ( -   3044.d0 * u(i,j-9,k)     &      
     &                     +  30663.d0 * u(i,j-8,k)     &     
     &                     - 139248.d0 * u(i,j-7,k)     &     
     &                     + 375704.d0 * u(i,j-6,k)     &     
     &                     - 667800.d0 * u(i,j-5,k)     &     
     &                     + 818874.d0 * u(i,j-4,k)     &     
     &                     - 704368.d0 * u(i,j-3,k)     &     
     &                     + 422568.d0 * u(i,j-2,k)     &     
     &                     - 165924.d0 * u(i,j-1,k)     &     
     &                     +  32575.d0 * u(i,j  ,k)     &     
     &                   ) * idy_sqrd_by_5040 
        end do
        end do
       
        return
      end subroutine 




      !----------------------------------------------------------------
      !
      ! deriv8888_zz calculates the second derivative with respect to  
      !              z using a centered eighth order finite difference 
      !              operator at all points on the interval [5,nz-4].  
      ! 
      !              At k=4 and k=nz-3, we use once-shifted eighth order 
      !              operators.  This has a ten point stencil.  
      ! 
      !              At k=3 and k=nz-2, we use twice-shifted eighth order 
      !              operators.  This has a ten point stencil.  
      ! 
      !              At k=2 and k=nz-1, we use thrice-shifted eighth order
      !              operators.  This has a ten point stencil.  
      ! 
      !              At k=1 and k=nz, we use totally shifted eighth order
      !              operators.  This has a ten point stencil.  
      ! 
      !----------------------------------------------------------------
      subroutine deriv8888_zz(DzDzu, u, dz, nx, ny, nz, dd_type)
        implicit   none        
        CCTK_INT  nx, ny, nz
        CCTK_REAL DzDzu(nx,ny,nz), u(nx,ny,nz)
        CCTK_REAL dz
        CCTK_INT   d_type, dd_type 

        ! local vars
        CCTK_INT i, j, k        
        CCTK_REAL idz, idz_sqrd, idz_sqrd_by_12 
        CCTK_REAL idz_sqrd_by_180, idz_sqrd_by_5040  

        idz_sqrd         = 1.0d0 /(dz*dz)         
        idz_sqrd_by_12   = idz_sqrd / 12.d0
        idz_sqrd_by_180  = idz_sqrd / 180.d0
        idz_sqrd_by_5040 = idz_sqrd / 5040.d0

        do j = 1, ny
        do i = 1, nx
          k=1
          DzDzu(i,j,k) = (    32575.d0 * u(i,j,k  )     &      
     &                     - 165924.d0 * u(i,j,k+1)     &     
     &                     + 422568.d0 * u(i,j,k+2)     &     
     &                     - 704368.d0 * u(i,j,k+3)     &     
     &                     + 818874.d0 * u(i,j,k+4)     &     
     &                     - 667800.d0 * u(i,j,k+5)     &     
     &                     + 375704.d0 * u(i,j,k+6)     &     
     &                     - 139248.d0 * u(i,j,k+7)     &     
     &                     +  30663.d0 * u(i,j,k+8)     &     
     &                     -   3044.d0 * u(i,j,k+9)     &     
     &                   ) * idz_sqrd_by_5040  
          k=2
          DzDzu(i,j,k) = (    3044.d0 * u(i,j,k-1)     &      
     &                     +  2135.d0 * u(i,j,k  )     &     
     &                     - 28944.d0 * u(i,j,k+1)     &     
     &                     + 57288.d0 * u(i,j,k+2)     &     
     &                     - 65128.d0 * u(i,j,k+3)     &     
     &                     + 51786.d0 * u(i,j,k+4)     &     
     &                     - 28560.d0 * u(i,j,k+5)     &     
     &                     + 10424.d0 * u(i,j,k+6)     &     
     &                     -  2268.d0 * u(i,j,k+7)     &     
     &                     +   223.d0 * u(i,j,k+8)     &     
     &                   ) * idz_sqrd_by_5040 
          k=3
          DzDzu(i,j,k) = ( -   223.d0 * u(i,j,k-2)     &      
     &                     +  5274.d0 * u(i,j,k-1)     &     
     &                     -  7900.d0 * u(i,j,k  )     &     
     &                     -  2184.d0 * u(i,j,k+1)     &      
     &                     + 10458.d0 * u(i,j,k+2)     &      
     &                     -  8932.d0 * u(i,j,k+3)     &      
     &                     +  4956.d0 * u(i,j,k+4)     &      
     &                     -  1800.d0 * u(i,j,k+5)     &      
     &                     +   389.d0 * u(i,j,k+6)     &      
     &                     -    38.d0 * u(i,j,k+7)     &      
     &                   ) * idz_sqrd_by_5040  
          k=4
          DzDzu(i,j,k) = (      38.d0 * u(i,j,k-3)     &      
     &                     -   603.d0 * u(i,j,k-2)     &     
     &                     +  6984.d0 * u(i,j,k-1)     &     
     &                     - 12460.d0 * u(i,j,k  )     &     
     &                     +  5796.d0 * u(i,j,k+1)     &     
     &                     +   882.d0 * u(i,j,k+2)     &     
     &                     -   952.d0 * u(i,j,k+3)     &     
     &                     +   396.d0 * u(i,j,k+4)     &     
     &                     -    90.d0 * u(i,j,k+5)     &     
     &                     +     9.d0 * u(i,j,k+6)     &     
     &                   ) * idz_sqrd_by_5040
          do k = 5, nz-4
            DzDzu(i,j,k) = ( -     9.d0 * u(i,j,k-4)     &      
     &                       +   128.d0 * u(i,j,k-3)     &     
     &                       -  1008.d0 * u(i,j,k-2)     &     
     &                       +  8064.d0 * u(i,j,k-1)     &     
     &                       - 14350.d0 * u(i,j,k  )     &     
     &                       +  8064.d0 * u(i,j,k+1)     &     
     &                       -  1008.d0 * u(i,j,k+2)     &     
     &                       +   128.d0 * u(i,j,k+3)     &     
     &                       -     9.d0 * u(i,j,k+4)     &     
     &                     ) * idz_sqrd_by_5040 
          end do
          k=nz-3
          DzDzu(i,j,k) = (       9.d0 * u(i,j,k-6)     &      
     &                     -    90.d0 * u(i,j,k-5)     &     
     &                     +   396.d0 * u(i,j,k-4)     &     
     &                     -   952.d0 * u(i,j,k-3)     &     
     &                     +   882.d0 * u(i,j,k-2)     &     
     &                     +  5796.d0 * u(i,j,k-1)     &     
     &                     - 12460.d0 * u(i,j,k  )     &     
     &                     +  6984.d0 * u(i,j,k+1)     &     
     &                     -   603.d0 * u(i,j,k+2)     &     
     &                     +    38.d0 * u(i,j,k+3)     &     
     &                   ) * idz_sqrd_by_5040 
          k=nz-2 
          DzDzu(i,j,k) = ( -    38.d0 * u(i,j,k-7)     &      
     &                     +   389.d0 * u(i,j,k-6)     &     
     &                     -  1800.d0 * u(i,j,k-5)     &     
     &                     +  4956.d0 * u(i,j,k-4)     &     
     &                     -  8932.d0 * u(i,j,k-3)     &     
     &                     + 10458.d0 * u(i,j,k-2)     &     
     &                     -  2184.d0 * u(i,j,k-1)     &     
     &                     -  7900.d0 * u(i,j,k  )     &     
     &                     +  5274.d0 * u(i,j,k+1)     &      
     &                     -   223.d0 * u(i,j,k+2)     &      
     &                   ) * idz_sqrd_by_5040 
          k=nz-1 
          DzDzu(i,j,k) = (     223.d0 * u(i,j,k-8)     &      
     &                     -  2268.d0 * u(i,j,k-7)     &      
     &                     + 10424.d0 * u(i,j,k-6)     &      
     &                     - 28560.d0 * u(i,j,k-5)     &      
     &                     + 51786.d0 * u(i,j,k-4)     &      
     &                     - 65128.d0 * u(i,j,k-3)     &      
     &                     + 57288.d0 * u(i,j,k-2)     &      
     &                     - 28944.d0 * u(i,j,k-1)     &      
     &                     +  2135.d0 * u(i,j,k  )     &      
     &                     +  3044.d0 * u(i,j,k+1)     &      
     &                   ) * idz_sqrd_by_5040  
          k=nz 
          DzDzu(i,j,k) = ( -   3044.d0 * u(i,j,k-9)     &      
     &                     +  30663.d0 * u(i,j,k-8)     &     
     &                     - 139248.d0 * u(i,j,k-7)     &     
     &                     + 375704.d0 * u(i,j,k-6)     &     
     &                     - 667800.d0 * u(i,j,k-5)     &     
     &                     + 818874.d0 * u(i,j,k-4)     &     
     &                     - 704368.d0 * u(i,j,k-3)     &     
     &                     + 422568.d0 * u(i,j,k-2)     &     
     &                     - 165924.d0 * u(i,j,k-1)     &     
     &                     +  32575.d0 * u(i,j,k  )     &     
     &                   ) * idz_sqrd_by_5040 
        end do
        end do
       
        return
      end subroutine 

      !-------------------------------------------------------------------
      !
      !
      !
      !-------------------------------------------------------------------
      subroutine D42NoMask_wrap(du, u, direction, par)
        use HYPER_DERIVS_42
        implicit none
        CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
        CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
        CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
        CCTK_INT :: direction

!       interface
!         subroutine D42(Du, u, direction, x, y, z, mask, imask, par)
!           implicit none
!           CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
!           CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
!           CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
!           CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: x, y, z
!           CCTK_INT :: direction
!           CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: mask
!           CCTK_INT, dimension(:,:,:), INTENT(IN)   :: imask
!         end subroutine
!         subroutine D42NoMask(Du, u, direction, par)
!           implicit none
!           CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
!           CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
!           CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
!           CCTK_INT :: direction
!         end subroutine
!       end interface

        call D42NoMask(Du, u, direction, par)

        return
      end subroutine



