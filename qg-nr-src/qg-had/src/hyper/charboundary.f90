!-------------------------------------------------------------------------
!
! $Id$
!
!-------------------------------------------------------------------------

#include "cctk.h"

MODULE M_CHAR_BOUNDARY
  use exact_u
  use charvars
  use params
  use GF
  implicit none

CONTAINS
!---------------------------------------------------------------------------
!
!  Some comments: what is called u in this subroutine is actually 
!  u or u_dot, depending on whether we are giving boundary conditions 
!  to the rhs or the vars.
!
!____________________________________________________________________________
  
  subroutine char_boundary(u, dxu, dyu, dzu, udot, v, w, imask, par)
    type(gridfunction), dimension(:), intent(inout) :: u, dxu,dyu,dzu,udot
    type(gridfunction), dimension(NV) :: v
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask

    ! local vars
    CCTK_INT :: lnu
    CCTK_REAL, dimension(:,:,:), pointer :: xx, yy, zz
    CCTK_REAL, DIMENSION(3) :: n  ! outward unit normal
    CCTK_INT :: imin,imax,jmin,jmax,kmin,kmax, IBBOX1, IBBOX2, IBBOX3,&
         IBBOX4, IBBOX5, IBBOX6, nx, ny, nz, ng,k, i, j, l, NB    
    CCTK_INT :: idx, idy, idz, mask_stencil
    logical bitantz
    logical double_equal

    logical, parameter :: ltrace = .false.
    

    if (ltrace) write(*,*)' @@@@@  charboundary '

    lnu = size(u)

    IBBOX1 = nint(par(P_BBOX1))
    IBBOX2 = nint(par(P_BBOX2))
    IBBOX3 = nint(par(P_BBOX3))
    IBBOX4 = nint(par(P_BBOX4))
    IBBOX5 = nint(par(P_BBOX5))
    IBBOX6 = nint(par(P_BBOX6)) 


    if (ltrace) then
       write(*,*)' charboundary: bbox=',ibbox1, ibbox2, ibbox3, &
                ibbox4, ibbox5, ibbox6
    end if

    ! check if we are using xy-plane bitant symmetry:  modify IBBOX
    ! variables accordingly
!    bitantz = ( (CCTK_Equals(domain,"bitant").eq.1) .and. &
!                (CCTK_Equals(bitant_plane,"xy").eq.1) )
!    if(bitantz) IBBOX5 = 0
   
    nx = int(par(P_NX))  ! local size of arrays
    ny = int(par(P_NY))
    nz = int(par(P_NZ))
    ng = nint(par(P_BOUND_WIDTH))

    xx => w(H_XPHYS)%d
    yy => w(H_YPHYS)%d
    zz => w(H_ZPHYS)%d

    !####### do through the boundary and give data ############
       if (ltrace) write(0,*)'### START FACES ************************'

       if ( par(P_BOUNDARY_CONDITIONS) .eq. 10 ) then
         do l=1,lnu
           call fix_outer(udot(l)%d,nx,ny,nz,ng,par)
         end do
       else
    !---------- faces -----------!

    x_min: if (IBBOX1 == 1) then
       if (ltrace.and. .not.double_equal(xx(1,1,1),-5.d0)) then
          write(*,*) 'charboundary: Xmin boundary do all the mess: ',nx,ny,nz
          write(*,*) 'charboundary: span x: ',xx(1,1,1),xx(nx,ny,nz)
	STOP
       end if
       imin=1
       imax=1
       jmin=2
       jmax=ny-1
       kmin=2
       kmax=nz-1
       call do_all_the_mess("xmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min
    
    x_max: if (IBBOX2 == 1) then
       if (ltrace.and. .not.double_equal(xx(nx,1,1),5.d0)) then
          write(*,*) 'charboundary: Xmax boundary do all the mess: ',nx,ny,nz
          write(*,*) 'charboundary: span x: ',xx(1,1,1),xx(nx,ny,nz)
	STOP
       end if
       imin=nx
       imax=nx
       jmin=2
       jmax=ny-1
       kmin=2
       kmax=nz-1
       call do_all_the_mess("xmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max

    y_min: if (IBBOX3 == 1) then
       if (ltrace.and. .not.double_equal(yy(nx,1,1),-5.d0)) then
          write(*,*) 'charboundary: Ymin boundary do all the mess: ',nx,ny,nz
          write(*,*) 'charboundary: span y: ',yy(1,1,1),yy(nx,ny,nz)
	STOP
       end if
       imin=2
       imax=nx-1
       jmin=1
       jmax=1
       kmin=2
       kmax=nz-1
       call do_all_the_mess("ymin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if y_min
    
    y_max: if (IBBOX4 == 1) then
       if (ltrace.and. .not.double_equal(yy(nx,ny,1),5.d0)) then
          write(*,*) 'charboundary: Ymax boundary do all the mess: ',nx,ny,nz
          write(*,*) 'charboundary: span y: ',yy(1,1,1),yy(nx,ny,nz)
	STOP
       end if
       imin=2
       imax=nx-1
       jmin=ny
       jmax=ny
       kmin=2
       kmax=nz-1
       call do_all_the_mess("ymax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if y_max

    z_min: if (IBBOX5 == 1) then 
       if (ltrace.and. .not.double_equal(zz(nx,ny,1),-5.d0)) then
          write(*,*) 'charboundary: Zmin boundary do all the mess: ',nx,ny,nz
          write(*,*) 'charboundary: span z: ',zz(1,1,1),zz(nx,ny,nz)
	STOP
       end if
       imin=2
       imax=nx-1
       jmin=2
       jmax=ny-1
       kmin=1
       kmax=1
       call do_all_the_mess("zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if z_min

    z_max: if (IBBOX6 == 1) then
       if (ltrace.and. .not.double_equal(zz(nx,ny,nz),5.d0)) then
          write(*,*) 'charboundary: Zmax boundary do all the mess: ',nx,ny,nz
          write(*,*) 'charboundary: span z: ',zz(1,1,1),zz(nx,ny,nz)
	STOP
       end if
       imin=2
       imax=nx-1
       jmin=2
       jmax=ny-1
       kmin=nz
       kmax=nz
       call do_all_the_mess("zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if z_max

       if (ltrace) write(0,*)'### END FACES, START EDGES *********************'

    !------------ edges ----------!
    x_min_y_min: if (IBBOX1*IBBOX3 == 1) then
       imin=1
       imax=1
       jmin=1
       jmax=1
       kmin=2
       kmax=nz-1
       call do_all_the_mess("xmin_ymin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_y_min

    x_max_y_min: if (IBBOX2*IBBOX3 == 1) then
       imin=nx
       imax=nx
       jmin=1
       jmax=1
       kmin=2
       kmax=nz-1
       call do_all_the_mess("xmax_ymin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_y_min

    x_min_y_max: if (IBBOX1*IBBOX4 == 1) then
       imin=1
       imax=1
       jmin=ny
       jmax=ny
       kmin=2
       kmax=nz-1
       call do_all_the_mess("xmin_ymax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_y_max

    x_max_y_max : if (IBBOX2*IBBOX4 == 1) then
       imin=nx
       imax=nx
       jmin=ny
       jmax=ny
       kmin=2
       kmax=nz-1
       call do_all_the_mess("xmax_ymax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_y_max

    x_min_z_min : if (IBBOX1*IBBOX5 == 1) then
       imin=1
       imax=1
       jmin=2
       jmax=ny-1
       kmin=1
       kmax=1
       call do_all_the_mess("xmin_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_z_min

    x_max_z_min : if (IBBOX2*IBBOX5 == 1) then
       imin=nx
       imax=nx
       jmin=2
       jmax=ny-1
       kmin=1
       kmax=1
       call do_all_the_mess("xmax_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_z_min

    x_min_z_max : if (IBBOX1*IBBOX6 == 1) then
       imin=1
       imax=1
       jmin=2
       jmax=ny-1
       kmin=nz
       kmax=nz
       call do_all_the_mess("xmin_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_z_max

    x_max_z_max : if (IBBOX2*IBBOX6 == 1) then
       imin=nx
       imax=nx
       jmin=2
       jmax=ny-1
       kmin=nz
       kmax=nz
       call do_all_the_mess("xmax_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_z_max

    y_min_z_min : if (IBBOX3*IBBOX5 == 1) then
       imin=2
       imax=nx-1
       jmin=1
       jmax=1
       kmin=1
       kmax=1
       call do_all_the_mess("ymin_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if y_min_z_min

    y_max_z_min : if (IBBOX4*IBBOX5 == 1) then
       imin=2
       imax=nx-1
       jmin=ny
       jmax=ny
       kmin=1
       kmax=1
       call do_all_the_mess("ymax_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if y_max_z_min

    y_min_z_max : if (IBBOX3*IBBOX6 == 1) then
       imin=2
       imax=nx-1
       jmin=1
       jmax=1
       kmin=nz
       kmax=nz
       call do_all_the_mess("ymin_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if y_min_z_max

    y_max_z_max : if (IBBOX4*IBBOX6 == 1) then
       imin=2
       imax=nx-1
       jmin=ny
       jmax=ny
       kmin=nz
       kmax=nz
       call do_all_the_mess("ymax_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if y_max_z_max

       if (ltrace) write(0,*)'### END EDGES, DO CORNERS'

    !--------- corners ---------!

    x_min_y_min_z_min : if (IBBOX1*IBBOX3*IBBOX5 == 1) then
       imin=1
       imax=1
       jmin=1
       jmax=1
       kmin=1
       kmax=1
       call do_all_the_mess("xmin_ymin_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_y_min_z_min

    x_max_y_min_z_min : if (IBBOX2*IBBOX3*IBBOX5 == 1) then
       imin=nx
       imax=nx
       jmin=1
       jmax=1
       kmin=1
       kmax=1
       call do_all_the_mess("xmax_ymin_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_y_min_z_min

    x_min_y_max_z_min : if (IBBOX1*IBBOX4*IBBOX5 == 1) then
       imin=1
       imax=1
       jmin=ny
       jmax=ny
       kmin=1
       kmax=1
       call do_all_the_mess("xmin_ymax_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_y_max_z_min

    x_min_y_min_z_max : if (IBBOX1*IBBOX3*IBBOX6 == 1) then
       imin=1
       imax=1
       jmin=1
       jmax=1
       kmin=nz
       kmax=nz
       call do_all_the_mess("xmin_ymin_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_y_min_z_max

    x_min_y_max_z_max : if (IBBOX1*IBBOX4*IBBOX6 == 1) then
       imin=1
       imax=1
       jmin=ny
       jmax=ny
       kmin=nz
       kmax=nz
       call do_all_the_mess("xmin_ymax_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_min_y_max_z_max

    x_max_y_min_z_max : if (IBBOX2*IBBOX3*IBBOX6 == 1) then
       imin=nx
       imax=nx
       jmin=1
       jmax=1
       kmin=nz
       kmax=nz
       call do_all_the_mess("xmax_ymin_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_y_min_z_max

    x_max_y_max_z_min : if (IBBOX2*IBBOX4*IBBOX5 == 1) then
       imin=nx
       imax=nx
       jmin=ny
       jmax=ny
       kmin=1
       kmax=1
       call do_all_the_mess("xmax_ymax_zmin",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_y_max_z_min

    x_max_y_max_z_max : if (IBBOX2*IBBOX4*IBBOX6 == 1) then
       imin=nx
       imax=nx
       jmin=ny
       jmax=ny
       kmin=nz
       kmax=nz
       call do_all_the_mess("xmax_ymax_zmax",imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
    end if x_max_y_max_z_max

!  End if statement to fix outer boundary
    end if
       if (ltrace) write(0,*)'### END CORNERS'


    mask_stencil = P_STENCIL_X_MASK + P_STENCIL_Y_MASK + P_STENCIL_Z_MASK

    if (nint(par(P_USE_MASK)) == 1 .and. &
        nint(par(P_INNER_BOUND_DATA)) == 1) then

       if (ltrace) write(0,*)'### INNER BOUNDARY****************************'

      ! loop over the interior of the box looking for mask boundaries.
      do k = 2, nz-1
      do j = 2, ny-1
      do i = 2, nx-1
        if (imask(i,j,k) .ne. 0) then
          if (imask(i,j,k) .ne. mask_stencil) then
            idx = imask(i,j,k) / 10000
            idy = mod(imask(i,j,k)/100, 100)
            idz = mod(imask(i,j,k), 100)

            !------------------------------------------------------
            ! Check faces
            !------------------------------------------------------

            ! x faces
            if (idx .eq. P_STENCIL_LEFT) then
              call do_all_the_mess("xmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idx .eq. P_STENCIL_RIGHT) then
              call do_all_the_mess("xmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)

            ! y faces
            else if (idy .eq. P_STENCIL_LEFT) then
              call do_all_the_mess("ymin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idy .eq. P_STENCIL_RIGHT) then
              call do_all_the_mess("ymax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)

            ! z faces
            else if (idz .eq. P_STENCIL_LEFT) then
              call do_all_the_mess("zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idz .eq. P_STENCIL_RIGHT) then
              call do_all_the_mess("zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)

            !------------------------------------------------------
            ! Check edges
            !------------------------------------------------------
 
            ! x-y edges
            else if (idx==P_STENCIL_EDGE_P .and. idy==P_STENCIL_EDGE_P) then
              call do_all_the_mess("xmin_ymin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idx==P_STENCIL_EDGE_M .and. idy==P_STENCIL_EDGE_P) then
              call do_all_the_mess("xmax_ymin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idx==P_STENCIL_EDGE_P .and. idy==P_STENCIL_EDGE_M) then
              call do_all_the_mess("xmin_ymax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idx==P_STENCIL_EDGE_M .and. idy==P_STENCIL_EDGE_M) then
              call do_all_the_mess("xmax_ymax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)

            ! x-z edges
            else if (idx==P_STENCIL_EDGE_P .and. idz==P_STENCIL_EDGE_P) then
              call do_all_the_mess("xmin_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idx==P_STENCIL_EDGE_M .and. idz==P_STENCIL_EDGE_P) then
              call do_all_the_mess("xmax_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idx==P_STENCIL_EDGE_P .and. idz==P_STENCIL_EDGE_M) then
              call do_all_the_mess("xmin_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idx==P_STENCIL_EDGE_M .and. idz==P_STENCIL_EDGE_M) then
              call do_all_the_mess("xmax_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)

            ! y-z edges
            else if (idy==P_STENCIL_EDGE_P .and. idz==P_STENCIL_EDGE_P) then
              call do_all_the_mess("ymin_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idy==P_STENCIL_EDGE_M .and. idz==P_STENCIL_EDGE_P) then
              call do_all_the_mess("ymax_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idy==P_STENCIL_EDGE_P .and. idz==P_STENCIL_EDGE_M) then
              call do_all_the_mess("ymin_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            else if (idy==P_STENCIL_EDGE_M .and. idz==P_STENCIL_EDGE_M) then
              call do_all_the_mess("ymax_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,v,w)
            
            !------------------------------------------------------
            ! Check vertices
            !------------------------------------------------------

            ! x+, y+, z+
            else if (idx==P_STENCIL_VERTEX_P .and. idy==P_STENCIL_VERTEX_P &
                .and. idz==P_STENCIL_VERTEX_P) then
              call do_all_the_mess("xmin_ymin_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            ! x+, y+, z-
            else if (idx==P_STENCIL_VERTEX_P .and. idy==P_STENCIL_VERTEX_P &
                .and. idz==P_STENCIL_VERTEX_M) then
              call do_all_the_mess("xmin_ymin_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            ! x+, y-, z+
            else if (idx==P_STENCIL_VERTEX_P .and. idy==P_STENCIL_VERTEX_M &
                .and. idz==P_STENCIL_VERTEX_P) then
              call do_all_the_mess("xmin_ymax_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            ! x+, y-, z-
            else if (idx==P_STENCIL_VERTEX_P .and. idy==P_STENCIL_VERTEX_M &
                .and. idz==P_STENCIL_VERTEX_M) then
              call do_all_the_mess("xmin_ymax_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            ! x-, y+, z+
            else if (idx==P_STENCIL_VERTEX_M .and. idy==P_STENCIL_VERTEX_P &
                .and. idz==P_STENCIL_VERTEX_P) then
              call do_all_the_mess("xmax_ymin_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            ! x-, y+, z-
            else if (idx==P_STENCIL_VERTEX_M .and. idy==P_STENCIL_VERTEX_P &
                .and. idz==P_STENCIL_VERTEX_M) then
              call do_all_the_mess("xmax_ymin_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            ! x-, y-, z+
            else if (idx==P_STENCIL_VERTEX_M .and. idy==P_STENCIL_VERTEX_M &
                .and. idz==P_STENCIL_VERTEX_P) then
              call do_all_the_mess("xmax_ymax_zmin",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            ! x-, y-, z-
            else if (idx==P_STENCIL_VERTEX_M .and. idy==P_STENCIL_VERTEX_M &
                .and. idz==P_STENCIL_VERTEX_M) then
              call do_all_the_mess("xmax_ymax_zmax",i,i,j,j,k,k,u,dxu,dyu,dzu,udot,par,&
                                    v,w)

            end if
 

          end if
        end if
      end do
      end do
      end do
    end if

  if (ltrace) write(*,*) 'charboundary: End.'

  end subroutine char_boundary


! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


  subroutine do_all_the_mess(direction,imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,v,w)
  !subroutine do_all_the_mess(direction,imin,imax,jmin,jmax,kmin,kmax,u,dxu,dyu,dzu,udot,par,w)
    use params
    use GF
    character(len=*), intent(in):: direction
    CCTK_INT, intent(in):: imin,imax,jmin,jmax,kmin,kmax
    type(gridfunction), dimension(:)  ::  u,dxu,dyu,dzu,udot
    type(gridfunction), dimension(NV) ::  v
    type(gridfunction), dimension(NW) ::  w
    CCTK_REAL, dimension(:)     :: par

    ! we are generous, we allocate the maximum posible number
    ! plus some extra
    CCTK_REAL, dimension(2*NU) :: w_in,w_out,w_0, wexact_in,wexact_out,wexact_0, sources
    CCTK_REAL, dimension(2*NU) :: wdot_in,wdot_out,wdot_0, wdotexact_in,wdotexact_out,wdotexact_0 

    CCTK_REAL, dimension(NU) :: u_pt,dxu_pt,dyu_pt,dzu_pt,uexact_pt,udot_pt,udotexact_pt,uexact_pt_orig,u_pt2
    CCTK_REAL, dimension(NW) :: w_pt
    CCTK_REAL, dimension(NV) :: v_pt
    CCTK_INT:: i,j,k,l,ii,jj,dir1,dir2,lnu
    CCTK_REAL :: x,y,z, random,dx,penalty,amp_random_bc, time, t0, sigma_t, amp_boundary, r, r0, sigma_rho,x1,x2, rho
    CCTK_REAL, dimension(3,3) :: T,Texact
    logical, parameter :: debug=.false.
    logical, parameter :: ltrace=.false.
    logical :: face

    !-----------------------------------------------

    lnu = size(u)
    dx = par(P_DX)
    penalty = par(P_PENALTY)
    amp_random_bc = par(P_AMP_RANDOM_BC)
    time = par(P_TIME)
    t0 = par(P_T0)
    sigma_t = par(P_SIGMA_T)
    sigma_rho = par(P_SIGMA_RHO)

    ! pick up the two tangential coords at each face (used to give pulses in time as bc)
    face=.false.
    ! These sometimes do not get set such as when direction is 'xmax_ymax_zmin' and then
    ! further down incorrect memory is accessed. So set to something reasonable and assume safe:
    dir1 = H_YPHYS
    dir2 = H_YPHYS
    if (direction=="xmin".or.direction=="xmax") then
       face=.true.
       dir1 = H_YPHYS
       dir2 = H_ZPHYS
    else if (direction=="ymin".or.direction=="ymax") then
       face=.true.
       dir1 = H_XPHYS
       dir2 = H_ZPHYS
    else if (direction=="zmin".or.direction=="zmax") then
       face=.true.
       dir1 = H_XPHYS
       dir2 = H_YPHYS
    end if

    if (face) then
       amp_boundary = par(P_AMP_BOUNDARY)
    else
       amp_boundary = 0.0
    end if

    if (ltrace) then
       write(*,*) 'do_all_the_mess: Begin'
       write(*,*) 'do_all_the_mess: dx        = ',dx
       write(*,*) 'do_all_the_mess: penalty   = ',penalty
       write(*,*) 'do_all_the_mess: sigma_rho = ',sigma_rho
       write(*,*) 'do_all_the_mess: NV        = ',NV
       write(*,*) 'do_all_the_mess: amp_boundary',amp_boundary
       write(*,*) 'do_all_the_mess: direction = ',direction
       write(*,*) 'do_all_the_mess: dir1      = ',dir1
       write(*,*) 'do_all_the_mess: dir2      = ',dir2
    end if

    ! loop over boundary
    do i=imin,imax
       do j=jmin,jmax
          do k=kmin,kmax
             
             if (debug) then
                do l=1,lnu
                   call random_number(random)
                   ! I do not redefine the metric as a random vars, since it would not be 
                   ! positive definite and charvars would complain when normalizing the 
                   ! vectors
                   u(l)%d(i,j,k) = u(l)%d(i,j,k)*(1.0 + 0.01*random)
                end do
             end if
             
             !write(*,*) 'do_all_the_mess: point A'
             do l=1,lnu
                u_pt(l)=u(l)%d(i,j,k)
             !write(*,*) 'do_all_the_mess: point B'
	        if(u(l)%take_dx.eq.1) then
                  dxu_pt(l) = dxu(l)%d(i,j,k)
	        else
	          dxu_pt(l) = 0.0
	        end if
             !write(*,*) 'do_all_the_mess: point c'
                if(u(l)%take_dy.eq.1) then
                  dyu_pt(l) = dyu(l)%d(i,j,k)
                else
                  dyu_pt(l) = 0.0
	        end if
             !write(*,*) 'do_all_the_mess: point d'
                if(u(l)%take_dz.eq.1) then
                  dzu_pt(l) = dzu(l)%d(i,j,k)
                else
                  dzu_pt(l) = 0.0
		end if
             !write(*,*) 'do_all_the_mess: point e'
                udot_pt(l)=udot(l)%d(i,j,k)
             end do
             if (ltrace) write(*,*) 'do_all_the_mess: point f'
             do l=1,NW
                w_pt(l)=w(l)%d(i,j,k)
             end do
             if (ltrace) write(*,*) 'do_all_the_mess: point g'
             do l=1,NV
                v_pt(l)=v(l)%d(i,j,k)
             end do
             if (ltrace) write(*,*) 'do_all_the_mess: point h',v_pt(1)

             x = w_pt(H_XPHYS)
             y = w_pt(H_YPHYS)
             z = w_pt(H_ZPHYS)   
             r = sqrt(x**2+y**2+z**2)

             if (ltrace) write(*,*) 'do_all_the_mess: point i'

             x1 = w_pt(DIR1)
             x2 = w_pt(DIR2)  
             if (ltrace) write(*,*) 'do_all_the_mess: ',DIR1, DIR2, NW
             if (ltrace) write(*,*) 'do_all_the_mess: ',x1,x2, direction
             rho = sqrt(x1**2 + x2**2)
             
             if (ltrace) write(*,*) 'do_all_the_mess: ',x,y,z,r,x1,x2,rho

             if (debug) then
                call prim2char(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt,uexact_pt,dxu_pt,dyu_pt,dzu_pt,sources)
             end if


             do_not_debug: if (.not.debug) then
                
                !######################  boundary conditions to the rhs ##################################
                
                if (ltrace) write(*,*) 'do_all_the_mess: exactu'
                call exactu(par,x,y,z,uexact_pt,w_pt) ! we might or might need this, compute it to simplify logic
                if (ltrace) write(*,*) 'do_all_the_mess: '
                
                uexact_pt_orig = uexact_pt ! a safety copy, for a check below
	        u_pt2 = u_pt

                boundary_cond_type: if (par(P_BOUNDARY_CONDITIONS).gt.2) then
                   !------- do not impose boundary conditions on the nonlinear terms of the principal part ------!

                   if (par(P_BOUNDARY_CONDITIONS).eq.3.or.par(P_BOUNDARY_CONDITIONS).eq.5) then
                      if (ltrace) write(*,*) 'do_all_the_mess: BC==5'
                      
                      ! characteristic decomposition for the rhs
                      call prim2char(par,v_pt,w_pt,direction,udot_pt,&
                                     wdot_in,wdot_out,wdot_0,T,u_pt2,&
                                     uexact_pt,dxu_pt,dyu_pt,dzu_pt,sources)
                      ! check the prim2char and char2prim are the inverse of each other
                      call char2prim(par,v_pt,w_pt,direction,udot_pt,wdot_in,wdot_out,wdot_0,T,u_pt2,uexact_pt)
                      if (ltrace) write(*,*) 'do_all_the_mess: BC==5'
                      do l=1,lnu
                         if (abs(udot(l)%d(i,j,k)-udot_pt(l)).gt.1.0D-02*(abs(udot_pt(l))+1.)) then 
                            write(*,*) "prim2char followed by char2prim does not give the identity dot",& 
                                 udot(l)%d(i,j,k), udot_pt(l), abs(udot(l)%d(i,j,k)-udot_pt(l))
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) abs(udot(l)%d(i,j,k)-udot_pt(l))/(abs(udot_pt(l))+1.d0)
                            write(*,*) ' '
                            call my_exit('prim2char on char2prim not identity')
                         end if
                      end do
                      if (ltrace) write(*,*) 'do_all_the_mess: BC==5'

                      ! characteristic decomposition for the vars
                      call prim2char(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt2,uexact_pt,dxu_pt,dyu_pt,dzu_pt,sources)
                      ! check the prim2char and char2prim are the inverse of each other
                      call char2prim(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt2,uexact_pt)
                      do l=1,lnu
                         !if (abs(u(l)%d(i,j,k)-u_pt(l)).gt.1.0E-05*(1.+abs(udot_pt(l)))) then 
                         if (abs(u(l)%d(i,j,k)-u_pt(l)).gt.1.0E-02*(1.+abs(u_pt(l)))) then 
                            write(*,*) "prim2char followed by char2prim does not give the identity func",& 
                                 u(l)%d(i,j,k), u_pt(l), abs(u(l)%d(i,j,k)-u_pt(l))
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) abs(u(l)%d(i,j,k)-u_pt(l))/(abs(u_pt(l))+1.d0)
                            write(*,*) ' '
                            call my_exit('prim2char on char2prim not identity')
                         end if
                      end do
                      if (ltrace) write(*,*) 'do_all_the_mess: BC==5'
                      
                      !*** A la Carpenter ****!
                      if (par(P_BOUNDARY_CONDITIONS).eq.3) then
                         do jj=1,lnu
                            wdot_in(jj) = wdot_in(jj) - penalty/dx*w_in(jj)
                         end do
!			    wdot_0(22) = wdot_0(22) - penalty/dx*(w_0(22)-1.)
                         !*** A la Olsson ***!
                      else if (par(P_BOUNDARY_CONDITIONS).eq.5) then
                         wdot_in = 0.0
                      end if
                      

                      !-------  Impose boundary conditions on the nonlinear terms of the principal part -------!
                      
                   else if (par(P_BOUNDARY_CONDITIONS).eq.4.or.par(P_BOUNDARY_CONDITIONS).eq.6) then    
                      
                      ! characteristic decomposition for the exact vars
                      call prim2char(par,v_pt,w_pt,direction,uexact_pt,wexact_in,wexact_out,wexact_0,Texact,u_pt,uexact_pt,&
                           dxu_pt, dyu_pt, dzu_pt, sources)

                      ! check the prim2char and char2prim are the inverse of each other
                      call char2prim(par,v_pt,w_pt,direction,uexact_pt,wexact_in,wexact_out,wexact_0,Texact,u_pt,uexact_pt)
                      do l=1,lnu
                         if (abs(uexact_pt_orig(l)-uexact_pt(l)).gt.1.0E-02*(1.+abs(uexact_pt(l)))) then 
                            write(*,*) "prim2char followed by char2prim does not give the identity, call 1",& 
                                 uexact_pt_orig(l), uexact_pt(l), abs(uexact_pt_orig(l)-uexact_pt(l)), l
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) ' '
                            call my_exit('prim2char on char2prim not identity')
                         end if
                      end do

                      ! characteristic decomposition for the rhs
                      call prim2char(par,v_pt,w_pt,direction,udot_pt,wdot_in,wdot_out,wdot_0,T,u_pt,uexact_pt,&
                           dxu_pt, dyu_pt, dzu_pt, sources)
                      ! check the prim2char and char2prim are the inverse of each other
                      call char2prim(par,v_pt,w_pt,direction,udot_pt,wdot_in,wdot_out,wdot_0,T,u_pt,uexact_pt)
                      do l=1,lnu
                         if (abs(udot(l)%d(i,j,k)-udot_pt(l)).gt.1.0E-02*(1.+abs(udot_pt(l)))) then 
                            write(*,*) "prim2char followed by char2prim does not give the identity, call 2",& 
                                 udot(l)%d(i,j,k), udot_pt(l), abs(udot(l)%d(i,j,k)-udot_pt(l)), l
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) ' '
                            call my_exit('prim2char on char2prim not identity')
                         end if
                      end do
                      
                      ! characteristic decomposition for the vars
                      call prim2char(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt,uexact_pt,dxu_pt,dyu_pt,dzu_pt,sources)
                      !--- make sure everything is Ok before giving boundary data--!
                      ! check the prim2char and char2prim are the inverse of each other
                      call char2prim(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt,uexact_pt)
                      do l=1,lnu
                         if (abs(u(l)%d(i,j,k)-u_pt(l)).gt.1.0E-02*(1.+abs(u_pt(l)))) then 
                            write(*,*) "prim2char followed by char2prim does not give the identity, call 3",& 
                                 u(l)%d(i,j,k), u_pt(l), abs(u(l)%d(i,j,k)-u_pt(l)), l
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) ' '
                            call my_exit('prim2char on char2prim not identity')
                         end if
                      end do

                      !*** A la Carpenter ****!
                      if (par(P_BOUNDARY_CONDITIONS).eq.4) then 
                         do jj=1,lnu
                            wdot_in(jj) = wdot_in(jj) -  penalty/dx*(w_in(jj)-wexact_in(jj))
                         end do
                         !*** A la Olsson ***!
                      else if (par(P_BOUNDARY_CONDITIONS).eq.6) then
                         wdot_in = 0.0
                      end if

		      
                   else if (par(P_BOUNDARY_CONDITIONS).eq.7.or.par(P_BOUNDARY_CONDITIONS).eq.8) then
                      
                      ! ****** Constraint preserving boundary conditions (only worked out for hyperGH) *********!
                      ! the "sources" variable when calling prim2char is used to give boundary conditions 
                      !in this case. 
                      ! Therefore, "sources" is not a priori given, as its name may imply (except for the freely 
                      ! specifiable sources in the cbpc, of course).
                      
                      call prim2char(par,v_pt,w_pt,direction,udot_pt,wdot_in,wdot_out,wdot_0,T,u_pt,uexact_pt,&
                           dxu_pt, dyu_pt, dzu_pt, sources)                      
		      if (debug) then   
                         ! check the prim2char and char2prim are the inverse of each other
                         call char2prim(par,v_pt,w_pt,direction,udot_pt,wdot_in,wdot_out,wdot_0,T,u_pt,uexact_pt)
                         do l=1,lnu
                            if (abs(udot(l)%d(i,j,k)-udot_pt(l)).gt.1.0E-02*(1.+abs(udot_pt(l)))) then 
                               write(*,*) "prim2char followed by char2prim does not give the identity",& 
                                    udot(l)%d(i,j,k), udot_pt(l), abs(udot(l)%d(i,j,k)-udot_pt(l))
                               write(*,*) "at point i/j/k: ",i,j,k
                               write(*,*) "for field l: ",l
                               write(*,*) ' '
                               call my_exit('prim2char on char2prim not identity')
                            end if
                         end do
		      end if	 

                      wdot_in = sources                      
                      
                   else 
                      write(*,*) "********** Unknown boundary conditions type ***********", par(P_BOUNDARY_CONDITIONS)
                      stop
                      
                   end if
                   
                   call char2prim(par,v_pt,w_pt,direction,udot_pt,wdot_in,wdot_out,wdot_0,T,u_pt,uexact_pt) 

                   !############# boundary conditions to the incoming characteristic variables ################
                else if (par(P_BOUNDARY_CONDITIONS).eq.2.or.par(P_BOUNDARY_CONDITIONS).eq.1.) then
	           u_pt2 = u_pt 
                   call prim2char(par,v_pt,w_pt,direction,uexact_pt,wexact_in,wexact_out,wexact_0,Texact,u_pt,uexact_pt,&
                        dxu_pt,dyu_pt,dzu_pt,sources)
                   ! check the prim2char and char2prim are the inverse of each other


                   call char2prim(par,v_pt,w_pt,direction,u_pt,wexact_in,wexact_out,wexact_0,Texact,u_pt,uexact_pt)
                   do l=1,lnu
                      if (abs(uexact_pt_orig(l)-uexact_pt(l)).gt.1.0E-02*(1.+abs(uexact_pt(l)))) then 
                         write(*,*) "prim2char followed by char2prim does not give the identity func",& 
                                 uexact_pt_orig(l), uexact_pt(l), abs(uexact_pt_orig(l)-uexact_pt(l))
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) ' '
                         stop
                      end if
                   end do
 
                   !------- do not impose boundary conditions on the nonlinear terms of the principal part ------!
                   if (par(P_BOUNDARY_CONDITIONS).eq.1) then 


                      call prim2char(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt2,uexact_pt,dxu_pt,dyu_pt,dzu_pt,sources)
                      ! check the prim2char and char2prim are the inverse of each other
                      call char2prim(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt2,uexact_pt)
                      do l=1,lnu
                         if (abs(u(l)%d(i,j,k)-u_pt(l)).gt.1.0E-02*(1.+abs(u_pt(l)))) then 
                            write(*,*) "prim2char followed by char2prim does not give the identity func",& 
                                 u(l)%d(i,j,k), u_pt(l), abs(u(l)%d(i,j,k)-u_pt(l))
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) ' '
                            stop
                         end if
                      end do
                      
                      !------- Impose boundary conditions on the nonlinear terms of the principal part ------!
                   else if (par(P_BOUNDARY_CONDITIONS).eq.2) then 
                      
                      call prim2char(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt2,uexact_pt,dxu_pt,dyu_pt,dzu_pt,sources)
                      ! check the prim2char and char2prim are the inverse of each other
                      call char2prim(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt2,uexact_pt)
                      do l=1,lnu
                         if (abs(u(l)%d(i,j,k)-u_pt(l)).gt.1.0E-02*(1.+abs(u_pt(l)))) then 
                            write(*,*) "prim2char followed by char2prim does not give the identity func",& 
                                 u(l)%d(i,j,k), u_pt(l), abs(u(l)%d(i,j,k)-u_pt(l))
                            write(*,*) "at point i/j/k: ",i,j,k
                            write(*,*) "for field l: ",l
                            write(*,*) ' '
                            stop
                         end if
                      end do
                   end if
                   
                   ! give boundary data
                   w_in = wexact_in
                   call char2prim(par,v_pt,w_pt,direction,u_pt,w_in,w_out,w_0,T,u_pt,uexact_pt)            
                   
                   else 
                      write(*,*) "********* Boundary condition seems to be of unknown type *********",&
                           par(P_BOUNDARY_CONDITIONS)
                      stop

                end if boundary_cond_type
             end if do_not_debug
             
             ! go back from pointwise quantities to the 3d arrays
             
             do l=1,lnu
                if (.not.debug) then
                   if (par(P_BOUNDARY_CONDITIONS).eq.2.or.par(P_BOUNDARY_CONDITIONS).eq.1.) then
                      u(l)%d(i,j,k)= u_pt(l)
                   else if (par(P_BOUNDARY_CONDITIONS).gt.2) then
                      udot(l)%d(i,j,k)= udot_pt(l)
                   else 
                      write(*,*) "********* Boundary condition seems to be of unknown type *********",&
                           par(P_BOUNDARY_CONDITIONS)
                      stop
                   end if
                else ! check that chars2prim and prim2chars are the inverse of each other, with random data 
                   if (abs(u(l)%d(i,j,k)-u_pt(l)).gt.1.0E-12) then 
                      write(*,*) "oops2", u(l)%d(i,j,k), u_pt(l), abs(u(l)%d(i,j,k)-u_pt(l))
                      stop
                   end if
                end if

                
             end do
          end do
       end do
    end do
    
    if (ltrace) write(*,*) 'do_all_the_mess: End.'

  end subroutine do_all_the_mess

  !--------------------------------------------------------------------
  ! 
  ! 
  !   
  !--------------------------------------------------------------------
  subroutine apply_outflow(f, chr, nx, ny, nz, ng,par)
    implicit   none
    CCTK_INT                       :: nx, ny, nz, ng
    CCTK_REAL, dimension(:)           :: par
    CCTK_REAL, dimension(nx,ny,nz) :: f, chr
!   CCTK_INT :: IBBOX1, IBBOX2, IBBOX3,&
!               IBBOX4, IBBOX5, IBBOX6
    ! local vars
    CCTK_INT                       :: i, j, k, mychr, myng
    CCTK_REAL                      :: val_nan, x
    logical, parameter :: ltrace = .false.
    logical, parameter :: evil_bc = .false.

    include 'chr.inc'

    ! These are antiquated/obsolete and should not be used anymore:
    !IBBOX1 = nint(par(P_BBOX1))
    !IBBOX2 = nint(par(P_BBOX2))
    !IBBOX3 = nint(par(P_BBOX3))
    !IBBOX4 = nint(par(P_BBOX4))
    !IBBOX5 = nint(par(P_BBOX5))
    !IBBOX6 = nint(par(P_BBOX6))

!    if(ltrace) write(*,*) 'apply_outflow: Enter:'
!    if(ltrace) write(*,*) 'apply_outflow: ng: ',ng
    ! physics should not depend on parameters related to domain decomposition
    ! Not sure how best to pick this, but for now use 5
    myng = 3
    if (myng .ge. nx) myng = nx-1
    if (myng .ge. ny) myng = ny-1
    if (myng .ge. nz) myng = nz-1
!    if(ltrace) write(*,*) 'apply_outflow: Using myng: ',myng

    x = -1.0d0
    val_nan = sqrt(x)

    ! xmin
    do k = 1, nz
      do j = 1, ny
        mychr = NINT(chr(1,j,k))
        if (          mychr.ne.NINT(CHR_interior)  &
                .and. mychr.ne.NINT(CHR_amr_bdy)   &
                .and. mychr.ne.NINT(CHR_deco_bdy)  &
                .and. mychr.ne.NINT(CHR_Refl_bdy)) then
          do i = 1, myng
            f(i,j,k) = f(myng+1,j,k)
          end do
        else
          if (evil_bc) then
            do i = 1, myng
              f(i,j,k) = val_nan
            end do
          end if
        end if

      end do
    end do

    ! xmax
    do k = 1, nz
      do j = 1, ny
        mychr = NINT(chr(nx,j,k))
        if (        mychr.ne.NINT(CHR_interior) .and. &
                    mychr.ne.NINT(CHR_amr_bdy) .and. &
                    mychr.ne.NINT(CHR_deco_bdy) .and. &
                    mychr.ne.NINT(CHR_Refl_bdy)) then
 
          do i = nx-myng+1, nx
            f(i,j,k) = f(nx-myng,j,k)
          end do
        else
          if (evil_bc) then
            do i = nx-myng+1, nx
              f(i,j,k) = val_nan
            end do
          end if
        end if
      end do
    end do

    ! ymin
    do k = 1, nz
      do i = 1, nx
        mychr = NINT(chr(i,1,k))
        if (        mychr.ne.NINT(CHR_interior) .and.  &
                    mychr.ne.NINT(CHR_amr_bdy) .and.  &
                    mychr.ne.NINT(CHR_deco_bdy) .and.  &
                    mychr.ne.NINT(CHR_Refl_bdy)) then
          do j = 1, myng
            f(i,j,k) = f(i,myng+1,k)
          end do
        else
          if (evil_bc) then
            do j = 1, myng
              f(i,j,k) = val_nan
            end do
          end if
        end if
      end do
    end do

    ! ymax
    do k = 1, nz
      do i = 1, nx
        mychr = NINT(chr(i,ny,k))
        if (      mychr.ne.NINT(CHR_interior) .and. &
                  mychr.ne.NINT(CHR_amr_bdy) .and. &
                  mychr.ne.NINT(CHR_deco_bdy) .and. &
                  mychr.ne.NINT(CHR_Refl_bdy)) then
          do j = ny-myng+1, ny
            f(i,j,k) = f(i,ny-myng,k)
          end do
        else
          if (evil_bc) then
            do j = ny-myng+1, ny
              f(i,j,k) = val_nan
            end do
          end if
        end if
      end do
    end do

    ! zmin
    do j = 1, ny
      do i = 1, nx
        mychr = NINT(chr(i,j,1))
        if (        mychr.ne.NINT(CHR_interior) .and. &
                    mychr.ne.NINT(CHR_amr_bdy) .and. &
                    mychr.ne.NINT(CHR_deco_bdy) .and. &
                    mychr.ne.NINT(CHR_Refl_bdy)) then
 
          do k = 1, myng
            f(i,j,k) = f(i,j,myng+1)
          end do
        else
          if (evil_bc) then
            do k = 1, myng
              f(i,j,k) = val_nan
            end do
          end if
        end if
      end do
    end do

    ! zmax
    do j = 1, ny
      do i = 1, nx
        mychr = NINT(chr(i,j,nz))
        if (       mychr.ne.NINT(CHR_interior) .and. &
                   mychr.ne.NINT(CHR_amr_bdy) .and. &
                   mychr.ne.NINT(CHR_deco_bdy) .and. &
                   mychr.ne.NINT(CHR_Refl_bdy)) then
          do k = nz-myng+1, nz
            f(i,j,k) = f(i,j,nz-myng)
          end do
        else
          if (evil_bc) then
            do k = nz-myng+1, nz
              f(i,j,k) = val_nan
            end do
          end if
        end if
      end do
    end do

!    if(ltrace) write(*,*) 'apply_outflow: Done:'
    return
  end subroutine apply_outflow

  !--------------------------------------------------------------------
  ! fix_outer : routine used to zero out time 
  ! derivatives for fixed outer boundaries.
  !   
  !--------------------------------------------------------------------
  subroutine fix_outer(f, nx, ny, nz, ng,par)
    implicit   none
    CCTK_INT                       :: nx, ny, nz, ng
    CCTK_REAL, dimension(:)           :: par
    CCTK_REAL, dimension(nx,ny,nz) :: f
    CCTK_INT :: IBBOX1, IBBOX2, IBBOX3,&
                IBBOX4, IBBOX5, IBBOX6

    ! local vars
    CCTK_INT                       :: i, j, k

    IBBOX1 = nint(par(P_BBOX1))
    IBBOX2 = nint(par(P_BBOX2))
    IBBOX3 = nint(par(P_BBOX3))
    IBBOX4 = nint(par(P_BBOX4))
    IBBOX5 = nint(par(P_BBOX5))
    IBBOX6 = nint(par(P_BBOX6))

    ! xmin
    if (IBBOX1 == 1) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, ng
            f(i,j,k) = 0.0d0
          end do
      end do
    end do
    end if

    ! xmax
    if (IBBOX2 == 1) then
      do k = 1, nz
        do j = 1, ny
          do i = nx-ng+1, nx
            f(i,j,k) = 0.0d0
          end do
        end do
      end do
    end if

    ! ymin
    if (IBBOX3 == 1) then
      do k = 1, nz
        do i = 1, nx
          do j = 1, ng
            f(i,j,k) = 0.0d0
          end do
        end do
      end do
    end if

    ! ymax
    if (IBBOX4 == 1) then
      do k = 1, nz
        do i = 1, nx
          do j = ny-ng+1, ny
            f(i,j,k) = 0.0d0
          end do
        end do
      end do
    end if

    ! zmin
    if (IBBOX5 == 1) then
      do j = 1, ny
        do i = 1, nx
          do k = 1, ng
            f(i,j,k) = 0.0d0
          end do
        end do
      end do
    end if

    ! zmax
    if (IBBOX6 == 1) then
      do j = 1, ny
        do i = 1, nx
          do k = nz-ng+1, nz
            f(i,j,k) = 0.0d0
          end do
        end do
      end do
    end if

    return
  end subroutine fix_outer


end module M_CHAR_BOUNDARY



