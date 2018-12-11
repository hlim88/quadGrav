!----------------------------------------------------------------
!
! $Id: idfile.f90,v 1.38 2012-01-28 22:46:01 carlos Exp $
!
!_________________________________________________________________

#include "cctk.h"
#include "cctk_DefineThorn.h"

module m_gr_id_file
  use GF
  use params
  implicit none
  contains

    !----------------------------------------------------------------
    !
    !
    !
    !_________________________________________________________________
    subroutine file_data_gr(u, u2, w, xx, yy, zz, nx, ny, nz, par)

      implicit none
      integer                                 :: nx, ny, nz
      type(gridfunction), dimension(NU_G)     :: u, u2
      type(gridfunction), dimension(NW)       :: w
      CCTK_REAL, dimension(:)                 :: par
      CCTK_REAL xx(nx,ny,nz), yy(nx,ny,nz), zz(nx,ny,nz)

      !! local vars
      CCTK_INT                                :: lb(3), gs(3)
      CCTK_INT                                :: lev
      CCTK_INT                                :: interp
      CCTK_REAL                               :: minx0,miny0,minz0, hi
      CCTK_REAL                               :: minx, miny, minz
      CCTK_REAL                               :: xcenter, ycenter, zcenter
      CCTK_REAL                               :: phase, omega1, omega2
      integer, parameter                      :: nfiles = 28
      character*64                            :: filename(nfiles), infile(2,3)
      CCTK_INT                                :: fileindex(nfiles)
      CCTK_INT                                :: readfromfile(nfiles)
      CCTK_REAL                               :: funcdefault(nfiles)
      CCTK_INT                                :: shp(3), dnx, dny, dnz, dnr
      CCTK_INT                                :: ifc, istat
      CCTK_INT                                :: read_g, read_k, read_phir, read_phic
      CCTK_INT                                :: read_alp, read_b, read_psi, read_phim
      CCTK_INT                                :: read_rho

      logical, parameter                      :: ltrace = .false.
      CCTK_REAL, allocatable, dimension(:)    :: x1d, y1d, z1d, coords
      CCTK_REAL, allocatable, dimension(:,:,:)  :: fdata
      CCTK_REAL, allocatable, dimension(:)      :: coords1d
      CCTK_REAL, allocatable, dimension(:)      :: fdata1d      
      
      CCTK_INT i,j,k,idtype
      CCTK_REAL myl2norm3d
      
     !added for the boosted BS      
      CCTK_INT ::  li, lj, lk, lm, ll, lp, boost_type
      CCTK_REAL :: vx, vy, vz, v2, cgamma, temp1, temp2, g00_bg, g01_bg, g02_bg, g03_bg
      CCTK_REAL :: gold(0:3,0:3), Dgold(0:3,0:3,0:3), gnew(0:3,0:3), Dgnew(0:3,0:3,0:3)
      CCTK_REAL :: phirold(0:3), phicold(0:3), Dphirold(0:3), Dphicold(0:3)
      CCTK_REAL :: phirnew(0:3), phicnew(0:3), Dphirnew(0:3), Dphicnew(0:3)
      CCTK_REAL :: phimold(0:3), phimnew(0:3), Dphimold(0:3), Dphimnew(0:3)
      CCTK_REAL :: Jud(0:3,0:3), DJud(0:3,0:3,0:3)

      CCTK_REAL :: u_d0,u_d1,u_d2,u_d3,u_up1,u_up2,u_up3,u_up0
      CCTK_REAL :: g00,g01,g02,g03,h11,h12,h13,h22,h23,h33
      CCTK_REAL :: huu11,huu12,huu13,huu22,huu23,huu33
      CCTK_REAL :: deth,b1,b2,b3,alp2,alp
      CCTK_REAL :: guu11,guu12,guu13,guu22,guu23,guu33
      CCTK_REAL :: guu00,guu01,guu02,guu03
      CCTK_REAL :: vd1,vd2,vd3,vu1,vu2,vu3,Wl            

      interface
        subroutine read_sdf_file(filename, f, xf, yf, zf, nx, ny, nz, &
                            level, fdata, coords, dnx, dny, dnz,    &
                            interp, cctk_llb, cctk_gshp)
          implicit none
          CCTK_INT        nx, ny, nz, level
          CCTK_INT        dnx, dny, dnz
          CCTK_INT        interp, cctk_llb(3), cctk_gshp(3)
          CCTK_REAL       f(nx,ny,nz)
          CCTK_REAL       xf(nx), yf(ny), zf(nz)
          CCTK_REAL       fdata(dnx,dny,dnz)
          CCTK_REAL       coords(dnx+dny+dnz)
          character*(*)   filename
        end subroutine read_sdf_file
	
        subroutine read_sdf_file1d(filename, f, xf, yf, zf, nx, ny, nz,    &
                            level, fdata1d, coords1d, dnr, interp, xc, yc, zc)
          implicit none
          CCTK_INT        nx, ny, nz, level
          CCTK_INT        dnr
          CCTK_INT        interp
          CCTK_REAL       f(nx,ny,nz)
          CCTK_REAL       xf(nx), yf(ny), zf(nz)
          CCTK_REAL       fdata1d(dnr)
          CCTK_REAL       coords1d(dnr)
          character*(*)   filename
          CCTK_REAL        xc, yc, zc
        end subroutine read_sdf_file1d		
		
        subroutine sdf_file_mem(filename, level, shp)
          implicit none
          character*(*)   filename
          CCTK_INT         ::  shp(3), level
        end subroutine sdf_file_mem
      end interface

      if (ltrace) then
         write(*,*)'-------------------file_data_gr--------------------'
      end if

      idtype     = nint(par(P_IDTYPE))      
      boost_type = nint(par(P_BOOST_TYPE))            
      
      ! create the list of functions that may be read from files
      ! Note, derivatives of the metric and beta are not included
      ! here.  This is in line with the Brill data runs.
      filename(1)  = 'g11.sdf'
      filename(2)  = 'g12.sdf'
      filename(3)  = 'g13.sdf'
      filename(4)  = 'g22.sdf'
      filename(5)  = 'g23.sdf'
      filename(6)  = 'g33.sdf'
      filename(7)  = 'k11.sdf'
      filename(8)  = 'k12.sdf'
      filename(9)  = 'k13.sdf'
      filename(10) = 'k22.sdf'
      filename(11) = 'k23.sdf'
      filename(12) = 'k33.sdf'
      filename(13) = 'alpha.sdf'
      filename(14) = 'beta1.sdf'
      filename(15) = 'beta2.sdf'
      filename(16) = 'beta3.sdf'
      filename(17) = 'psi.sdf'
      filename(18) = 'phir.sdf'      
      filename(19) = 'pir.sdf'            
      filename(20) = 'phic.sdf'      
      filename(21) = 'pic.sdf'            
      filename(22) = 'phim.sdf'      
      filename(23) = 'pim.sdf'            
      filename(24) = 'rho.sdf'            
      filename(25) = 'press.sdf'            
      filename(26) = 'vx.sdf'            
      filename(27) = 'vy.sdf'            
      filename(28) = 'vz.sdf'            
      
      fileindex(1)  = H_G11
      fileindex(2)  = H_G12
      fileindex(3)  = H_G13
      fileindex(4)  = H_G22
      fileindex(5)  = H_G23
      fileindex(6)  = H_G33
      fileindex(7)  = H_K11
      fileindex(8)  = H_K12
      fileindex(9)  = H_K13
      fileindex(10) = H_K22
      fileindex(11) = H_K23
      fileindex(12) = H_K33
      fileindex(13) = H_G00
      fileindex(14) = H_G01
      fileindex(15) = H_G02
      fileindex(16) = H_G03
      fileindex(17) = H_G00 ! save psi in alp initially      
      fileindex(18) = H_PHIR      
      fileindex(19) = H_PIR        
      fileindex(20) = H_PHIC      
      fileindex(21) = H_PIC        
      fileindex(22) = H_PHIM
      fileindex(23) = H_PIM
      !copy the fluid in the phim        
!      fileindex(24) = F_PHIM
!      fileindex(25) = F_PIM        
      fileindex(24) = H_PHIM
      fileindex(25) = H_PIM        

      ! Determine which functions should be read from sdf files
      ! These are set by parameters, which should be 0 or 1.
      !     P_READ_FILE_G = 1 -> read metric from files
      !     P_READ_FILE_K = 1 -> read extrinsic curvature from files
      !     P_READ_FILE_A = 1 -> read lapse from files
      !     P_READ_FILE_DA= 1 -> read derivatives of lapse from files
      read_g     = nint(par(P_READ_FILE_G))
      read_k     = nint(par(P_READ_FILE_K))
      read_alp   = nint(par(P_READ_FILE_ALP))
      read_b     = nint(par(P_READ_FILE_B))      
      read_psi   = nint(par(P_READ_FILE_PSI))
      read_phir  = nint(par(P_READ_FILE_PHIR))      
      read_phic  = nint(par(P_READ_FILE_PHIC))            
      read_phim  = nint(par(P_READ_FILE_PHIM))      
      read_rho   = nint(par(P_READ_FILE_RHO))      
      
      do i = 1, 6
        readfromfile(i) = read_g
      end do
      do i = 7, 12
        readfromfile(i) = read_k
      end do
      readfromfile(13) = read_alp
      readfromfile(14) = read_b
      readfromfile(15) = read_b      
      readfromfile(16) = read_b
      readfromfile(17) = read_psi
      readfromfile(18) = read_phir
      readfromfile(19) = read_phir      
      readfromfile(20) = read_phic
      readfromfile(21) = read_phic
      readfromfile(22) = read_phim
      readfromfile(23) = read_phim
      readfromfile(24) = read_rho
      
      if (ltrace) then
         write(*,*)' read_g     = ',read_g
         write(*,*)' read_k     = ',read_k
         write(*,*)' read_alp   = ',read_alp
         write(*,*)' read_b     = ',read_b
         write(*,*)' read_psi   = ',read_psi
         write(*,*)' read_phir  = ',read_phir
         write(*,*)' read_phic  = ',read_phic 
         write(*,*)' read_phim  = ',read_phim 
         write(*,*)' read_rho   = ',read_rho 
      end if
   
      ! This array has (constant) default values for each function.  
      ! The defaults are flat space values
      funcdefault(1)  = 1.0d0  ! gxx
      funcdefault(2)  = 0.0d0  ! gxy
      funcdefault(3)  = 0.0d0  ! gxz
      funcdefault(4)  = 1.0d0  ! gyy
      funcdefault(5)  = 0.0d0  ! gyz
      funcdefault(6)  = 1.0d0  ! gzz
      funcdefault(7)  = 0.0d0  ! kxx
      funcdefault(8)  = 0.0d0  ! kxy
      funcdefault(9)  = 0.0d0  ! kxz
      funcdefault(10) = 0.0d0  ! kyy
      funcdefault(11) = 0.0d0  ! kyz
      funcdefault(12) = 0.0d0  ! kzz
      funcdefault(13) = 1.0d0  ! alpha
      funcdefault(14) = 0.0d0  ! beta1
      funcdefault(15) = 0.0d0  ! beta2
      funcdefault(16) = 0.0d0  ! beta3
      funcdefault(17) = 1.0d0  ! psi
      funcdefault(18) = 0.0d0  ! phir
      funcdefault(19) = 0.0d0  ! pir
      funcdefault(20) = 0.0d0  ! phic
      funcdefault(21) = 0.0d0  ! pic
      funcdefault(22) = 0.0d0  ! phim
      funcdefault(23) = 0.0d0  ! pim
      funcdefault(24) = 0.0d0  ! rho
      funcdefault(25) = 0.0d0  ! press
      funcdefault(26) = 0.0d0  ! vx
      funcdefault(27) = 0.0d0  ! vy
      funcdefault(28) = 0.0d0  ! vz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     only the scalar fields of a binary boson star 
!     superposed with a Lorene binary NS later on
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      if (idtype .EQ. -5) then  

!       the center of the first boson star
!       look at the position of NS and set it by hand
!       in any case, they will be symmetric
        xcenter = par(P_BS1_X0)
	ycenter = par(P_BS1_Y0)	
	zcenter = par(P_BS1_Z0)
        omega1 = par(P_SF_OMEGA1)
        omega2 = par(P_SF_OMEGA2)
	
        lev = nint(par(P_READ_DATA_LEVEL))
        interp = nint(par(P_INTERP_ID))
      
        call sdf_file_mem("phir1d.sdf", lev, shp)

        allocate(fdata1d(shp(1)),coords1d(shp(1)),STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for arrays to read sdf files'
          stop
        end if

        ! it is most useful to have 1-d coordinate arrays.
        allocate (x1d(nx), y1d(ny), z1d(nz), STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for 1d coordinate arrays'
          stop
        end if
        do i = 1, nx
          x1d(i) = xx(i,1,1)
        end do
        do i = 1, ny
          y1d(i) = yy(1,i,1)
        end do
        do i = 1, nz
          z1d(i) = zz(1,1,i)
        end do

       minx0 = par(P_MINX0)
       miny0 = par(P_MINY0)
       minz0 = par(P_MINZ0)

       hi   = x1d(2) - x1d(1)
       !
       ! Need minimums for this grid alone:
       !
       minx = x1d(1)
       miny = y1d(1)
       minz = z1d(1)

       dnr = shp(1)

       ! for the lapses
       call read_sdf_file1d("alp1d.sdf", u(H_G01)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       call read_sdf_file1d("alp1d.sdf", u(H_G02)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, -xcenter, -ycenter, -zcenter)

       ! for the phis
       call read_sdf_file1d("phir1d.sdf", u(H_PHIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       call read_sdf_file1d("phim1d.sdf", u(H_PHIM)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, -xcenter, -ycenter, -zcenter)

       u(H_PIC)%d = omega1*u(H_PHIR)%d/u(H_G01)%d
       u(H_PIN)%d = omega2*u(H_PHIM)%d/u(H_G02)%d       
       u(H_G01)%d = 0.0d0
       u(H_G02)%d = 0.0d0

      ! setting the the first order derivates of the metric
       call default_data_derivs(u, w, par)      

       if (ltrace) then
	 do ifc=1,23
	   write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
	 end do 
       end if        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!    only a single boson star from 1D initial data
!       idtype=-6,   without anything else
!       idtype=-11,  with a BH (in Kerr Schild coordinates)
!       idtype=-12,  with a BH (in isotropic coordinates)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      else if ((idtype .EQ. -6) .OR. (idtype .EQ. -11) .OR. (idtype .EQ. -12))  then  

        xcenter = par(P_BS1_X0)
        ycenter = par(P_BS1_Y0)
        zcenter = par(P_BS1_Z0)
        omega1  = par(P_SF_OMEGA1)       
        vx      = par(P_BS1_VX)
        vy      = par(P_BS1_VY)
        vz      = par(P_BS1_VZ)
                  
        lev = nint(par(P_READ_DATA_LEVEL))
        interp = nint(par(P_INTERP_ID))
     
        call sdf_file_mem("psi1d.sdf", lev, shp)

        allocate(fdata1d(shp(1)),coords1d(shp(1)),STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for arrays to read sdf files'
          stop
        end if

        ! it is most useful to have 1-d coordinate arrays.
        allocate (x1d(nx), y1d(ny), z1d(nz), STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for 1d coordinate arrays'
          stop
        end if
        do i = 1, nx
          x1d(i) = xx(i,1,1)
        end do
        do i = 1, ny
          y1d(i) = yy(1,i,1)
        end do
        do i = 1, nz
          z1d(i) = zz(1,1,i)
        end do

       minx0 = par(P_MINX0)
       miny0 = par(P_MINY0)
       minz0 = par(P_MINZ0)

       hi   = x1d(2) - x1d(1)
       !
       ! Need minimums for this grid alone:
       !
       minx = x1d(1)
       miny = y1d(1)
       minz = z1d(1)

       dnr = shp(1)
      
       call read_sdf_file1d("psi1d.sdf", u(H_G00)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
       u(H_G11)%d = (u(H_G00)%d)**4
       u(H_G22)%d = u(H_G11)%d
       u(H_G33)%d = u(H_G11)%d
       u(H_G12)%d = 0.0d0
       u(H_G13)%d = 0.0d0
       u(H_G23)%d = 0.0d0      
       u(H_K11)%d = 0.0d0
       u(H_K22)%d = 0.0d0
       u(H_K33)%d = 0.0d0
       u(H_K12)%d = 0.0d0
       u(H_K13)%d = 0.0d0
       u(H_K23)%d = 0.0d0      
      
       call read_sdf_file1d("alpha1d.sdf", u(H_G01)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
       u(H_G00)%d =-(u(H_G01)%d)**2
       u(H_G01)%d = 0.0d0
       u(H_G02)%d = 0.0d0
       u(H_G03)%d = 0.0d0

       call read_sdf_file1d("phir1d.sdf", u(H_PIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

!     with the new stress-energy tensor which is twice the old one
       u(H_PHIR)%d = u(H_PIR)%d*sqrt(2.0)


!       u(H_PHIR)%d = u(H_PIR)%d

       u(H_PIC)%d = -omega1*u(H_PHIR)%d       
       u(H_PIR)%d = 0.0d0
       u(H_PIM)%d = 0.0d0       
       u(H_PHIC)%d = 0.0d0
       u(H_PHIM)%d = 0.0d0

       if (read_rho .EQ. 1) then
         call read_sdf_file1d("rho1d.sdf", u(H_PHIM)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)
         call read_sdf_file1d("press1d.sdf", u(H_PIM)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)
       end if

      ! setting the the first order derivates of the metric
      call default_data_derivs(u, w, par)              
       
       if (ltrace) then
         do ifc=1,23
           write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
         end do 
       end if        

      ! make the boost
       if (boost_type .EQ. 0) then
         print*, "STOP!!!!: boost_type has to be 1=lorentzian for the BH+BS superposition"
	 STOP
       end if
       	 
       do i = 1, nx; do j = 1, ny; do k = 1, nz
#include "boost.inc"
       end do; end do; end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     binary boson star 
!        idtype = -7 , from the same 1D ID, boson/boson pair or boson/antiboson pair
!        idtype = -8 , from different files, unequal mass pair
!        idtype = -9 , special case for boson/boson in phase opposition pair
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      else if ((idtype .EQ. -7) .OR. (idtype .EQ. -8) .OR. (idtype .EQ. -9)) then  

!       the center of the boson stars, symmetric         
        xcenter = par(P_BS1_X0)
	ycenter = par(P_BS1_Y0)	
	zcenter = par(P_BS1_Z0)
        omega1 = par(P_SF_OMEGA1)
        omega2 = par(P_SF_OMEGA2)
	
	infile(1,1) = "psi1d.sdf"
	infile(1,2) = "alpha1d.sdf"
	infile(1,3) = "phir1d.sdf"	

	if (idtype .EQ. -8) then
	  infile(2,1) = "spsi1d.sdf"
	  infile(2,2) = "salpha1d.sdf"
	  infile(2,3) = "sphir1d.sdf"
	else
	  infile(2,1) = "psi1d.sdf"
	  infile(2,2) = "alpha1d.sdf"
	  infile(2,3) = "phir1d.sdf"
        end if

        lev = nint(par(P_READ_DATA_LEVEL))
        interp = nint(par(P_INTERP_ID))
      
        call sdf_file_mem("psi1d.sdf", lev, shp)

        allocate(fdata1d(shp(1)),coords1d(shp(1)),STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for arrays to read sdf files'
          stop
        end if

        ! it is most useful to have 1-d coordinate arrays.
        allocate (x1d(nx), y1d(ny), z1d(nz), STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for 1d coordinate arrays'
          stop
        end if
        do i = 1, nx
          x1d(i) = xx(i,1,1)
        end do
        do i = 1, ny
          y1d(i) = yy(1,i,1)
        end do
        do i = 1, nz
          z1d(i) = zz(1,1,i)
        end do

       minx0 = par(P_MINX0)
       miny0 = par(P_MINY0)
       minz0 = par(P_MINZ0)

       hi   = x1d(2) - x1d(1)
       !
       ! Need minimums for this grid alone:
       !
       minx = x1d(1)
       miny = y1d(1)
       minz = z1d(1)

       dnr = shp(1)

!      psiT = psi1 + psi2 - 1      
       call read_sdf_file1d(infile(1,1), u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       call read_sdf_file1d(infile(2,1), u(H_K12)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, -xcenter, -ycenter, -zcenter)

       u(H_G11)%d = (u(H_K11)%d + u(H_K12)%d - 1.0d0)**4
       u(H_G22)%d = u(H_G11)%d
       u(H_G33)%d = u(H_G11)%d
       u(H_G12)%d = 0.0d0
       u(H_G13)%d = 0.0d0
       u(H_G23)%d = 0.0d0      
       u(H_K11)%d = 0.0d0

!      alpT = alp1 + alp2 - 1
!      g00T = g00_1 + g00_2 + 1
       call read_sdf_file1d(infile(1,2), u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       call read_sdf_file1d(infile(2,2), u(H_G00)%d,  &
                      x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		      coords1d, dnr, interp, -xcenter, -ycenter, -zcenter)
		      
       u(H_G00)%d = -(u(H_G00)%d)**2 -(u(H_K11)%d)**2 + 1.0d0      
       u(H_G01)%d = 0.0d0
       u(H_G02)%d = 0.0d0
       u(H_G03)%d = 0.0d0
       u(H_K11)%d = 0.0d0

       ! now the density and the pressure
       if (read_rho .EQ. 1) then
         call read_sdf_file1d("rho1d.sdf", u(H_PHIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)
         call read_sdf_file1d("press1d.sdf", u(H_PIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

         call read_sdf_file1d("rho1d.sdf", u(H_PHIC)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, -xcenter, -ycenter, -zcenter)
         call read_sdf_file1d("press1d.sdf", u(H_PIC)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, -xcenter, -ycenter, -zcenter)

         u(H_PHIM)%d = u(H_PHIR)%d + u(H_PHIC)%d
         u(H_PIM)%d = u(H_PIR)%d + u(H_PIC)%d
       end if


       ! for the phis
       call read_sdf_file1d(infile(1,3), u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       call read_sdf_file1d(infile(2,3), u(H_PHIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, -xcenter, -ycenter, -zcenter)

       ! special case for the oposition of phase       
       if (idtype .EQ. -9) then       
         u(H_PIC)%d = -omega1*(u(H_K11)%d - u(H_PHIR)%d) 
         u(H_PIR)%d = 0.0d0
          
         u(H_PHIR)%d = u(H_K11)%d - u(H_PHIR)%d 
         u(H_PHIC)%d = 0.0d0
       else 
         u(H_PIC)%d = -(omega1*u(H_K11)%d + omega2*u(H_PHIR)%d) 
         u(H_PIR)%d = 0.0d0	 

         u(H_PHIR)%d = u(H_K11)%d + u(H_PHIR)%d 
         u(H_PHIC)%d = 0.0d0       
       end if	               
       
       u(H_K11)%d = 0.0d0
       u(H_K22)%d = 0.0d0
       u(H_K33)%d = 0.0d0
       u(H_K12)%d = 0.0d0
       u(H_K13)%d = 0.0d0
       u(H_K23)%d = 0.0d0      
 
      ! setting the the first order derivates of the metric
       call default_data_derivs(u, w, par)      

       if (ltrace) then
	 do ifc=1,23
	   write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
	 end do 
       end if        

!      the particular case of a general phase
!       u(H_PIC)%d = u(H_K11)%d + u(H_PHIR)%d*dcos(phase) 
!       u(H_PIR)%d = u(H_PHIR)%d*dsin(phase) 
!       u(H_PIM)%d = 0.0d0    
          
!       u(H_PHIR)%d = u(H_K11)%d + u(H_PHIR)%d*dcos(phase) 
!       u(H_PHIC)%d = u(H_PHIR)%d*dsin(phase)
!       u(H_PHIM)%d = 0.0d0
                                
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     Galilean boosted boson star 
!        idtype = -21 , single one
!        idtype = -22 , any equal-mass boson pair from the same file
!        idtype = -23 , any equal-mass boson trinary from the same file
!        idtype = -24 , any boson trinary from different files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      else if ((idtype .EQ. -21) .OR. (idtype .EQ. -22) .OR. (idtype .EQ. -23) .OR. (idtype .EQ. -24)) then  
      
!       the center of the boson stars, symmetric         
        xcenter = par(P_BS1_X0)
	ycenter = par(P_BS1_Y0)	
	zcenter = par(P_BS1_Z0)
        omega1  = par(P_SF_OMEGA1)
        vx      = par(P_BS1_VX)
        vy      = par(P_BS1_VY)
        vz      = par(P_BS1_VZ)
        
        infile(1,1) = "psi1d.sdf"
        infile(1,2) = "alpha1d.sdf"
        infile(1,3) = "phir1d.sdf"

        infile(2,1) = "psi1d.sdf"
        infile(2,2) = "alpha1d.sdf"
        infile(2,3) = "phir1d.sdf"
  
        lev = nint(par(P_READ_DATA_LEVEL))
        interp = nint(par(P_INTERP_ID))
      	
        call sdf_file_mem("psi1d.sdf", lev, shp)

        allocate(fdata1d(shp(1)),coords1d(shp(1)),STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for arrays to read sdf files'
          stop
        end if

        ! it is most useful to have 1-d coordinate arrays.
        allocate (x1d(nx), y1d(ny), z1d(nz), STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for 1d coordinate arrays'
          stop
        end if
        do i = 1, nx
          x1d(i) = xx(i,1,1)
        end do
        do i = 1, ny
          y1d(i) = yy(1,i,1)
        end do
        do i = 1, nz
          z1d(i) = zz(1,1,i)
        end do

       minx0 = par(P_MINX0)
       miny0 = par(P_MINY0)
       minz0 = par(P_MINZ0)

       hi   = x1d(2) - x1d(1)
       !
       ! Need minimums for this grid alone:
       !
       minx = x1d(1)
       miny = y1d(1)
       minz = z1d(1)

       dnr = shp(1)

!-----LOAD THE FIRST BOSON--------------------------------------      
!----- G, K, PHIR, PHIC ----------------------------------------
       !  the intrinsic metric g_{ij}
       call read_sdf_file1d("psi1d.sdf", u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       u(H_G11)%d =(u(H_K11)%d)**4
       u(H_G22)%d = u(H_G11)%d
       u(H_G33)%d = u(H_G11)%d
       u(H_G12)%d = 0.0d0
       u(H_G13)%d = 0.0d0
       u(H_G23)%d = 0.0d0      
       u(H_K11)%d = 0.0d0

       !  the gauge g_{0a}
       call read_sdf_file1d("alpha1d.sdf", u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
		      
       u(H_G00)%d =-(u(H_K11)%d)**2
       u(H_G01)%d = 0.0d0
       u(H_G02)%d = 0.0d0
       u(H_G03)%d = 0.0d0
       u(H_K11)%d = 0.0d0

       !  the complex scalar field
       call read_sdf_file1d("phir1d.sdf", u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       u(H_PIC)%d = -omega1*u(H_K11)%d 
       u(H_PIR)%d = 0.0d0	 

       u(H_PHIR)%d = u(H_K11)%d
       u(H_PHIC)%d = 0.0d0       

       ! now the density ; the pressure almost at the end        
       if (read_rho .EQ. 1) then
          call read_sdf_file1d("rho1d.sdf", u(H_PHIM)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

          ! set the four-velocity u_a = W*(v_a + n_a)
          ! v_a = 0 in the initial data; save u_0 = -alp*W in pir 
          do i = 1, nx
          do j = 1, ny
          do k = 1, nz
            if (u(H_PHIM)%d(i,j,k) .GT. 1e-10) then
              alp = sqrt(-u(H_G00)%d(i,j,k))
              u(H_PIM)%d(i,j,k)    =-alp
              u(H_D1PHIM)%d(i,j,k) = 0.0d0
              u(H_D2PHIM)%d(i,j,k) = 0.0d0
              u(H_D3PHIM)%d(i,j,k) = 0.0d0   
            end if         
          end do
          end do
          end do

       end if

       ! the time derivative of the metric g_ij
       ! if the K_ij=beta^i=0 then it is 0
       u(H_K11)%d = 0.0d0
       u(H_K22)%d = 0.0d0
       u(H_K33)%d = 0.0d0
       u(H_K12)%d = 0.0d0
       u(H_K13)%d = 0.0d0
       u(H_K23)%d = 0.0d0      

       ! the first order derivates of the metric
       call default_data_derivs(u, w, par)      

       ! make the boost
       do i = 1, nx; do j = 1, ny; do k = 1, nz
#include "boost.inc"               
       end do; end do; end do      

       ! set the three velocity and the pressure
       if (read_rho .EQ. 1) then

         ! get the three-velocity v_i from u_a
         ! the MHD code reads v^i where rho>0
         do i = 1, nx
         do j = 1, ny
         do k = 1, nz
           if (u(H_PHIM)%d(i,j,k) .GT. 1e-10) then

             u_d0 =  u(H_PIM)%d(i,j,k) 
             u_d1 =  u(H_D1PHIM)%d(i,j,k) 
             u_d2 =  u(H_D2PHIM)%d(i,j,k) 
             u_d3 =  u(H_D3PHIM)%d(i,j,k) 

             g00 = u(H_G00)%d(i,j,k)
             g01 = u(H_G01)%d(i,j,k)
             g02 = u(H_G01)%d(i,j,k)
             g03 = u(H_G01)%d(i,j,k)

             h11 = u(H_G11)%d(i,j,k)
             h12 = u(H_G12)%d(i,j,k)
             h13 = u(H_G13)%d(i,j,k)
             h22 = u(H_G22)%d(i,j,k)
             h23 = u(H_G23)%d(i,j,k)
             h33 = u(H_G33)%d(i,j,k)

             huu11 = -h23**2 + h22*h33
             huu12 = h13*h23 - h12*h33
             huu13 = -(h13*h22) + h12*h23
             huu22 = -h13**2 + h11*h33
             huu23 = h12*h13 - h11*h23
             huu33 = -h12**2 + h11*h22
             deth = h11*huu11 + h12*huu12 + h13*huu13
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
     &        + b2*b3*h23) + b3**2*h33
             alp = sqrt(alp2)

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

             u_up1 = huu11*u_d1 + huu12*u_d2 + huu13*u_d3
             u_up2 = huu12*u_d1 + huu22*u_d2 + huu23*u_d3
             u_up3 = huu13*u_d1 + huu23*u_d2 + huu33*u_d3
             Wl = sqrt(1.0d0 + u_up1*u_d1 + u_up2*u_d2 + u_up3*u_d3)
     
             vd1 = u_d1/Wl
             vd2 = u_d2/Wl
             vd3 = u_d3/Wl     
     
             vu1 = huu11*vd1 + huu12*vd2 + huu13*vd3
             vu2 = huu12*vd1 + huu22*vd2 + huu23*vd3
             vu3 = huu13*vd1 + huu23*vd2 + huu33*vd3

             u(H_D1PHIM)%d(i,j,k) = vu1
             u(H_D2PHIM)%d(i,j,k) = vu2
             u(H_D3PHIM)%d(i,j,k) = vu3
             u(H_PIM)%d(i,j,k)    = 0.0d0

           else 

             u(H_D1PHIM)%d(i,j,k) = 0.0d0
             u(H_D2PHIM)%d(i,j,k) = 0.0d0
             u(H_D3PHIM)%d(i,j,k) = 0.0d0
             u(H_PIM)%d(i,j,k)    = 0.0d0

           end if

         end do
         end do
         end do

         ! save the the pressure in pir
         call read_sdf_file1d("press1d.sdf", u(H_PIM)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       end if

!---------------------------------------------------------------              
!-----ADD THE SECOND BOSON--------------------------------------             
!---------------------------------------------------------------       
       if ((idtype .EQ. -22) .OR. (idtype .EQ. -23) .OR. (idtype .EQ. -24)) then  

!--------MOVE THE FIRST BOSON TO A TEMP VAR-----------------------               
         u2(H_G00)%d = u(H_G00)%d
         u2(H_G01)%d = u(H_G01)%d         
         u2(H_G02)%d = u(H_G02)%d	 
         u2(H_G03)%d = u(H_G03)%d	 
         u2(H_G11)%d = u(H_G11)%d	 
         u2(H_G12)%d = u(H_G12)%d	 
         u2(H_G13)%d = u(H_G13)%d
         u2(H_G22)%d = u(H_G22)%d
         u2(H_G23)%d = u(H_G23)%d         
         u2(H_G33)%d = u(H_G33)%d	 

         u2(H_D1G00)%d = u(H_D1G00)%d
         u2(H_D1G01)%d = u(H_D1G01)%d         
         u2(H_D1G02)%d = u(H_D1G02)%d	 
         u2(H_D1G03)%d = u(H_D1G03)%d	 
         u2(H_D1G11)%d = u(H_D1G11)%d	 
         u2(H_D1G12)%d = u(H_D1G12)%d	 
         u2(H_D1G13)%d = u(H_D1G13)%d
         u2(H_D1G22)%d = u(H_D1G22)%d
         u2(H_D1G23)%d = u(H_D1G23)%d         
         u2(H_D1G33)%d = u(H_D1G33)%d	 

         u2(H_D2G00)%d = u(H_D2G00)%d
         u2(H_D2G01)%d = u(H_D2G01)%d         
         u2(H_D2G02)%d = u(H_D2G02)%d	 
         u2(H_D2G03)%d = u(H_D2G03)%d	 
         u2(H_D2G11)%d = u(H_D2G11)%d	 
         u2(H_D2G12)%d = u(H_D2G12)%d	 
         u2(H_D2G13)%d = u(H_D2G13)%d
         u2(H_D2G22)%d = u(H_D2G22)%d
         u2(H_D2G23)%d = u(H_D2G23)%d         
         u2(H_D2G33)%d = u(H_D2G33)%d	 
	 
         u2(H_D3G00)%d = u(H_D3G00)%d
         u2(H_D3G01)%d = u(H_D3G01)%d         
         u2(H_D3G02)%d = u(H_D3G02)%d	 
         u2(H_D3G03)%d = u(H_D3G03)%d	 
         u2(H_D3G11)%d = u(H_D3G11)%d	 
         u2(H_D3G12)%d = u(H_D3G12)%d	 
         u2(H_D3G13)%d = u(H_D3G13)%d
         u2(H_D3G22)%d = u(H_D3G22)%d
         u2(H_D3G23)%d = u(H_D3G23)%d         
         u2(H_D3G33)%d = u(H_D3G33)%d	 

         u2(H_K00)%d = u(H_K00)%d 
         u2(H_K01)%d = u(H_K01)%d 
         u2(H_K02)%d = u(H_K02)%d 
         u2(H_K03)%d = u(H_K03)%d
         u2(H_K11)%d = u(H_K11)%d 
         u2(H_K12)%d = u(H_K12)%d 
         u2(H_K13)%d = u(H_K13)%d
         u2(H_K22)%d = u(H_K22)%d
         u2(H_K23)%d = u(H_K23)%d         
         u2(H_K33)%d = u(H_K33)%d 
	 
         u2(H_PHIR)%d = u(H_PHIR)%d 
         u2(H_PHIC)%d = u(H_PHIC)%d 
         u2(H_PHIM)%d = u(H_PHIM)%d

         u2(H_D1PHIR)%d = u(H_D1PHIR)%d 
         u2(H_D2PHIR)%d = u(H_D2PHIR)%d 
         u2(H_D3PHIR)%d = u(H_D3PHIR)%d

         u2(H_D1PHIC)%d = u(H_D1PHIC)%d 
         u2(H_D2PHIC)%d = u(H_D2PHIC)%d 
         u2(H_D3PHIC)%d = u(H_D3PHIC)%d

         u2(H_D1PHIM)%d = u(H_D1PHIM)%d 
         u2(H_D2PHIM)%d = u(H_D2PHIM)%d 
         u2(H_D3PHIM)%d = u(H_D3PHIM)%d
	 	 
         u2(H_PIR)%d = u(H_PIR)%d	 
         u2(H_PIC)%d = u(H_PIC)%d	 
         u2(H_PIM)%d = u(H_PIM)%d

!-----LOAD THE SECOND BOSON--------------------------------------      
!       the center of the boson stars, symmetric         
         xcenter = par(P_BS2_X0)
  	 ycenter = par(P_BS2_Y0)	
	 zcenter = par(P_BS2_Z0)
         omega1  = par(P_SF_OMEGA2)
         vx      = par(P_BS2_VX)
         vy      = par(P_BS2_VY)               
         vz      = par(P_BS2_VZ)
	 phase   = par(P_SF_PHASE)

         g00_bg = -1.0d0
         g01_bg = 0.0d0; g02_bg = 0.0d0; g03_bg = 0.0d0
	
		    		 
!----- G, K, PHIR, PHIC ----------------------------------------
         !  the intrinsic metric g_{ij}
         call read_sdf_file1d("psi1d.sdf", u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

         u(H_G11)%d =(u(H_K11)%d)**4
         u(H_G22)%d = u(H_G11)%d
         u(H_G33)%d = u(H_G11)%d
         u(H_G12)%d = 0.0d0
         u(H_G13)%d = 0.0d0
         u(H_G23)%d = 0.0d0      
         u(H_K11)%d = 0.0d0

         !  the gauge g_{0a}
         call read_sdf_file1d("alpha1d.sdf", u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
		      
         u(H_G00)%d =-(u(H_K11)%d)**2
         u(H_G01)%d = 0.0d0
         u(H_G02)%d = 0.0d0
         u(H_G03)%d = 0.0d0
         u(H_K11)%d = 0.0d0

         !  the complex scalar field
         call read_sdf_file1d("phir1d.sdf", u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
		   
         u(H_PIC)%d = -omega1 * u(H_K11)%d * cos(phase)
         u(H_PIR)%d = -omega1 * u(H_K11)%d * sin(phase)
         u(H_PIM)%d = 0.0d0       

         u(H_PHIR)%d = u(H_K11)%d * cos(phase)
         u(H_PHIC)%d =-u(H_K11)%d * sin(phase)       
         u(H_PHIM)%d = 0.0d0

         ! now the density ; the pressure almost at the end        
       if (read_rho .EQ. 1) then
            call read_sdf_file1d("rho1d.sdf", u(H_PHIM)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

            ! set the four-velocity u_a = W*(v_a + n_a)
            ! v_a = 0 in the initial data; save u_0 = -alp*W in pir 
            do i = 1, nx
            do j = 1, ny
            do k = 1, nz
              if (u(H_PHIM)%d(i,j,k) .GT. 1e-10) then
                alp = sqrt(-u(H_G00)%d(i,j,k))
                u(H_PIM)%d(i,j,k)    =-alp
                u(H_D1PHIM)%d(i,j,k) = 0.0d0
                u(H_D2PHIM)%d(i,j,k) = 0.0d0
                u(H_D3PHIM)%d(i,j,k) = 0.0d0   
              end if         
            end do
            end do
            end do

         end if

!	 print* ,"phase=",phase
!	 print*,"cos=",dcos(phase),"sin=",dsin(phase)
	 
         !  the time derivative of the metric       
         u(H_K11)%d = 0.0d0
         u(H_K22)%d = 0.0d0
         u(H_K33)%d = 0.0d0
         u(H_K12)%d = 0.0d0
         u(H_K13)%d = 0.0d0
         u(H_K23)%d = 0.0d0      

         ! the first order derivates of the metric
         call default_data_derivs(u, w, par)      

         ! make the boost
         do i = 1, nx; do j = 1, ny; do k = 1, nz
#include "boost.inc"               
         end do; end do; end do      	 

         ! set the three velocity and the pressure
         if (read_rho .EQ. 1) then
           ! get the three-velocity v_i from u_a
           ! the MHD code reads v^i where rho>0
           do i = 1, nx
           do j = 1, ny
           do k = 1, nz
             if (u(H_PHIM)%d(i,j,k) .GT. 1e-9) then

               u_d0 =  u(H_PIM)%d(i,j,k) 
               u_d1 =  u(H_D1PHIM)%d(i,j,k) 
               u_d2 =  u(H_D2PHIM)%d(i,j,k) 
               u_d3 =  u(H_D3PHIM)%d(i,j,k) 

               g00 = u(H_G00)%d(i,j,k)
               g01 = u(H_G01)%d(i,j,k)
               g02 = u(H_G01)%d(i,j,k)
               g03 = u(H_G01)%d(i,j,k)

               h11 = u(H_G11)%d(i,j,k)
               h12 = u(H_G12)%d(i,j,k)
               h13 = u(H_G13)%d(i,j,k)
               h22 = u(H_G22)%d(i,j,k)
               h23 = u(H_G23)%d(i,j,k)
               h33 = u(H_G33)%d(i,j,k)

               huu11 = -h23**2 + h22*h33
               huu12 = h13*h23 - h12*h33
               huu13 = -(h13*h22) + h12*h23
               huu22 = -h13**2 + h11*h33
               huu23 = h12*h13 - h11*h23
               huu33 = -h12**2 + h11*h22
               deth = h11*huu11 + h12*huu12 + h13*huu13
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
     &          + b2*b3*h23) + b3**2*h33
               alp = sqrt(alp2)

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

               u_up0 = guu00*u_d0 + guu01*u_d1 + guu02*u_d2 &
     &               + guu03*u_d3
               Wl    = alp*u_up0

               vd1 =  u_d1/Wl
               vd2 =  u_d2/Wl
               vd3 =  u_d3/Wl

               vu1 = huu11*vd1 + huu12*vd2 + huu13*vd3
               vu2 = huu12*vd1 + huu22*vd2 + huu23*vd3
               vu3 = huu13*vd1 + huu23*vd2 + huu33*vd3

               u(H_D1PHIM)%d(i,j,k) = vu1
               u(H_D2PHIM)%d(i,j,k) = vu2
               u(H_D3PHIM)%d(i,j,k) = vu3
               u(H_PIM)%d(i,j,k)    = 0.0d0

             else 

               u(H_D1PHIM)%d(i,j,k) = 0.0d0
               u(H_D2PHIM)%d(i,j,k) = 0.0d0
               u(H_D3PHIM)%d(i,j,k) = 0.0d0
               u(H_PIM)%d(i,j,k)    = 0.0d0

             end if

           end do
           end do
           end do

           ! save the the pressure in pir
           call read_sdf_file1d("press1d.sdf", u(H_PIM)%d,  &
                  x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                  coords1d, dnr, interp, xcenter, ycenter, zcenter)

         end if

	 	 	 
!--------COMPUTING THE SUPERPOSITION-----------------------         
         u(H_G00)%d = u2(H_G00)%d + u(H_G00)%d - g00_bg
         u(H_G01)%d = u2(H_G01)%d + u(H_G01)%d - g01_bg        
         u(H_G02)%d = u2(H_G02)%d + u(H_G02)%d - g02_bg
         u(H_G03)%d = u2(H_G03)%d + u(H_G03)%d - g03_bg	 
         u(H_G11)%d = u2(H_G11)%d + u(H_G11)%d - 1.0d0
         u(H_G12)%d = u2(H_G12)%d + u(H_G12)%d	 
         u(H_G13)%d = u2(H_G13)%d + u(H_G13)%d
         u(H_G22)%d = u2(H_G22)%d + u(H_G22)%d - 1.0d0
         u(H_G23)%d = u2(H_G23)%d + u(H_G23)%d         
         u(H_G33)%d = u2(H_G33)%d + u(H_G33)%d - 1.0d0	 
	 
         u(H_D1G00)%d = u2(H_D1G00)%d + u(H_D1G00)%d
         u(H_D1G01)%d = u2(H_D1G01)%d + u(H_D1G01)%d         
         u(H_D1G02)%d = u2(H_D1G02)%d + u(H_D1G02)%d	 
         u(H_D1G03)%d = u2(H_D1G03)%d + u(H_D1G03)%d	 
         u(H_D1G11)%d = u2(H_D1G11)%d + u(H_D1G11)%d	 
         u(H_D1G12)%d = u2(H_D1G12)%d + u(H_D1G12)%d	 
         u(H_D1G13)%d = u2(H_D1G13)%d + u(H_D1G13)%d
         u(H_D1G22)%d = u2(H_D1G22)%d + u(H_D1G22)%d
         u(H_D1G23)%d = u2(H_D1G23)%d + u(H_D1G23)%d         
         u(H_D1G33)%d = u2(H_D1G33)%d + u(H_D1G33)%d	 

         u(H_D2G00)%d = u2(H_D2G00)%d + u(H_D2G00)%d
         u(H_D2G01)%d = u2(H_D2G01)%d + u(H_D2G01)%d         
         u(H_D2G02)%d = u2(H_D2G02)%d + u(H_D2G02)%d	 
         u(H_D2G03)%d = u2(H_D2G03)%d + u(H_D2G03)%d	 
         u(H_D2G11)%d = u2(H_D2G11)%d + u(H_D2G11)%d	 
         u(H_D2G12)%d = u2(H_D2G12)%d + u(H_D2G12)%d	 
         u(H_D2G13)%d = u2(H_D2G13)%d + u(H_D2G13)%d
         u(H_D2G22)%d = u2(H_D2G22)%d + u(H_D2G22)%d
         u(H_D2G23)%d = u2(H_D2G23)%d + u(H_D2G23)%d         
         u(H_D2G33)%d = u2(H_D2G33)%d + u(H_D2G33)%d	 
	 
         u(H_D3G00)%d = u2(H_D3G00)%d + u(H_D3G00)%d
         u(H_D3G01)%d = u2(H_D3G01)%d + u(H_D3G01)%d         
         u(H_D3G02)%d = u2(H_D3G02)%d + u(H_D3G02)%d	 
         u(H_D3G03)%d = u2(H_D3G03)%d + u(H_D3G03)%d	 
         u(H_D3G11)%d = u2(H_D3G11)%d + u(H_D3G11)%d	 
         u(H_D3G12)%d = u2(H_D3G12)%d + u(H_D3G12)%d	 
         u(H_D3G13)%d = u2(H_D3G13)%d + u(H_D3G13)%d
         u(H_D3G22)%d = u2(H_D3G22)%d + u(H_D3G22)%d
         u(H_D3G23)%d = u2(H_D3G23)%d + u(H_D3G23)%d         
         u(H_D3G33)%d = u2(H_D3G33)%d + u(H_D3G33)%d	 

         u(H_K00)%d = u2(H_K00)%d + u(H_K00)%d
         u(H_K01)%d = u2(H_K01)%d + u(H_K01)%d	 
         u(H_K02)%d = u2(H_K02)%d + u(H_K02)%d	 
         u(H_K03)%d = u2(H_K03)%d + u(H_K03)%d
         u(H_K11)%d = u2(H_K11)%d + u(H_K11)%d	 
         u(H_K12)%d = u2(H_K12)%d + u(H_K12)%d	 
         u(H_K13)%d = u2(H_K13)%d + u(H_K13)%d
         u(H_K22)%d = u2(H_K22)%d + u(H_K22)%d
         u(H_K23)%d = u2(H_K23)%d + u(H_K23)%d         
         u(H_K33)%d = u2(H_K33)%d + u(H_K33)%d	 
	 
         u(H_PHIR)%d = u2(H_PHIR)%d + u(H_PHIR)%d
         u(H_PHIC)%d = u2(H_PHIC)%d + u(H_PHIC)%d	 
         u(H_PHIM)%d = u2(H_PHIM)%d + u(H_PHIM)%d

         u(H_D1PHIR)%d = u2(H_D1PHIR)%d	+ u(H_D1PHIR)%d 
         u(H_D2PHIR)%d = u2(H_D2PHIR)%d	+ u(H_D2PHIR)%d 	 
         u(H_D3PHIR)%d = u2(H_D3PHIR)%d	+ u(H_D3PHIR)%d 

         u(H_D1PHIC)%d = u2(H_D1PHIC)%d	+ u(H_D1PHIC)%d 
         u(H_D2PHIC)%d = u2(H_D2PHIC)%d	+ u(H_D2PHIC)%d 	 
         u(H_D3PHIC)%d = u2(H_D3PHIC)%d	+ u(H_D3PHIC)%d 

         u(H_D1PHIM)%d = u2(H_D1PHIM)%d	+ u(H_D1PHIM)%d 
         u(H_D2PHIM)%d = u2(H_D2PHIM)%d	+ u(H_D2PHIM)%d 	 
         u(H_D3PHIM)%d = u2(H_D3PHIM)%d	+ u(H_D3PHIM)%d 
	 	 
         u(H_PIR)%d = u2(H_PIR)%d + u(H_PIR)%d
         u(H_PIC)%d = u2(H_PIC)%d + u(H_PIC)%d	 
         u(H_PIM)%d = u2(H_PIM)%d + u(H_PIM)%d
	 
       end if
              
!---------------------------------------------------------------              
!-----ADD THE THIRD BOSON--------------------------------------             
!---------------------------------------------------------------       
       if ((idtype .EQ. -23) .OR. (idtype .EQ. -24)) then  

           write(*,*) '@@@ initial data not working idtype=',idtype
           stop

       end if        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!    the oscillaton star from 1D initial data
!       idtype=-31,  single one
!       idtype=-32,  two of them
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
      else if ((idtype .EQ. -31) .OR. (idtype .EQ. -32))  then  

        xcenter = par(P_BS1_X0)
        ycenter = par(P_BS1_Y0)
        zcenter = par(P_BS1_Z0)
        omega1  = par(P_SF_OMEGA1)       
        vx      = par(P_BS1_VX)
        vy      = par(P_BS1_VY)
        vz      = par(P_BS1_VZ)

	infile(1,1) = "psi1d.sdf"
	infile(1,2) = "alpha1d.sdf"
	infile(1,3) = "phir1d.sdf"	

!	if (idtype .EQ. -32) then
!	  infile(2,1) = "spsi1d.sdf"
!	  infile(2,2) = "salpha1d.sdf"
!	  infile(2,3) = "sphir1d.sdf"
!	else
!	  infile(2,1) = "psi1d.sdf"
!	  infile(2,2) = "alpha1d.sdf"
!	  infile(2,3) = "phir1d.sdf"
!        end if	    

                  
        lev = nint(par(P_READ_DATA_LEVEL))
        interp = nint(par(P_INTERP_ID))
     
        call sdf_file_mem("psi1d.sdf", lev, shp)

        allocate(fdata1d(shp(1)),coords1d(shp(1)),STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for arrays to read sdf files'
          stop
        end if

        ! it is most useful to have 1-d coordinate arrays.
        allocate (x1d(nx), y1d(ny), z1d(nz), STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for 1d coordinate arrays'
          stop
        end if
        do i = 1, nx
          x1d(i) = xx(i,1,1)
        end do
        do i = 1, ny
          y1d(i) = yy(1,i,1)
        end do
        do i = 1, nz
          z1d(i) = zz(1,1,i)
        end do

       minx0 = par(P_MINX0)
       miny0 = par(P_MINY0)
       minz0 = par(P_MINZ0)

       hi   = x1d(2) - x1d(1)
       !
       ! Need minimums for this grid alone:
       !
       minx = x1d(1)
       miny = y1d(1)
       minz = z1d(1)

       dnr = shp(1)
      
       call read_sdf_file1d("psi1d.sdf", u(H_G00)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
       u(H_G11)%d = (u(H_G00)%d)**4
       u(H_G22)%d = u(H_G11)%d
       u(H_G33)%d = u(H_G11)%d
       u(H_G12)%d = 0.0d0
       u(H_G13)%d = 0.0d0
       u(H_G23)%d = 0.0d0      
       u(H_K11)%d = 0.0d0
       u(H_K22)%d = 0.0d0
       u(H_K33)%d = 0.0d0
       u(H_K12)%d = 0.0d0
       u(H_K13)%d = 0.0d0
       u(H_K23)%d = 0.0d0      
      
       call read_sdf_file1d("alpha1d.sdf", u(H_G01)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
       u(H_G00)%d =-(u(H_G01)%d)**2
       u(H_G01)%d = 0.0d0
       u(H_G02)%d = 0.0d0
       u(H_G03)%d = 0.0d0

       call read_sdf_file1d("phir1d.sdf", u(H_PHIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
		   
       u(H_PIC)%d = 0.0d0    
       u(H_PIR)%d = 0.0d0
       u(H_PIM)%d = 0.0d0       
       
       u(H_PHIC)%d = 0.0d0
       u(H_PHIM)%d = 0.0d0


      ! setting the the first order derivates of the metric
      call default_data_derivs(u, w, par)              
       
       if (ltrace) then
         do ifc=1,23
           write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
         end do 
       end if        

      ! make the boost
       if (boost_type .EQ. 0) then
         print*, "STOP!!!!: boost_type has to be 1=lorentzian for the BH+BS superposition"
	 STOP
       end if
       	 
       do i = 1, nx; do j = 1, ny; do k = 1, nz
#include "boost.inc"
       end do; end do; end do

!---------------------------------------------------------------              
!-----ADD THE SECOND OSCILLATON---------------------------------             
!---------------------------------------------------------------       
       if (idtype .EQ. -32) then  
!         print*, "second oscillaton"

!--------MOVE THE FIRST OSCILON TO A TEMP VAR-----------------------               
         u2(H_G00)%d = u(H_G00)%d
         u2(H_G01)%d = u(H_G01)%d         
         u2(H_G02)%d = u(H_G02)%d	 
         u2(H_G03)%d = u(H_G03)%d	 
         u2(H_G11)%d = u(H_G11)%d	 
         u2(H_G12)%d = u(H_G12)%d	 
         u2(H_G13)%d = u(H_G13)%d
         u2(H_G22)%d = u(H_G22)%d
         u2(H_G23)%d = u(H_G23)%d         
         u2(H_G33)%d = u(H_G33)%d	 

         u2(H_D1G00)%d = u(H_D1G00)%d
         u2(H_D1G01)%d = u(H_D1G01)%d         
         u2(H_D1G02)%d = u(H_D1G02)%d	 
         u2(H_D1G03)%d = u(H_D1G03)%d	 
         u2(H_D1G11)%d = u(H_D1G11)%d	 
         u2(H_D1G12)%d = u(H_D1G12)%d	 
         u2(H_D1G13)%d = u(H_D1G13)%d
         u2(H_D1G22)%d = u(H_D1G22)%d
         u2(H_D1G23)%d = u(H_D1G23)%d         
         u2(H_D1G33)%d = u(H_D1G33)%d	 

         u2(H_D2G00)%d = u(H_D2G00)%d
         u2(H_D2G01)%d = u(H_D2G01)%d         
         u2(H_D2G02)%d = u(H_D2G02)%d	 
         u2(H_D2G03)%d = u(H_D2G03)%d	 
         u2(H_D2G11)%d = u(H_D2G11)%d	 
         u2(H_D2G12)%d = u(H_D2G12)%d	 
         u2(H_D2G13)%d = u(H_D2G13)%d
         u2(H_D2G22)%d = u(H_D2G22)%d
         u2(H_D2G23)%d = u(H_D2G23)%d         
         u2(H_D2G33)%d = u(H_D2G33)%d	 
	 
         u2(H_D3G00)%d = u(H_D3G00)%d
         u2(H_D3G01)%d = u(H_D3G01)%d         
         u2(H_D3G02)%d = u(H_D3G02)%d	 
         u2(H_D3G03)%d = u(H_D3G03)%d	 
         u2(H_D3G11)%d = u(H_D3G11)%d	 
         u2(H_D3G12)%d = u(H_D3G12)%d	 
         u2(H_D3G13)%d = u(H_D3G13)%d
         u2(H_D3G22)%d = u(H_D3G22)%d
         u2(H_D3G23)%d = u(H_D3G23)%d         
         u2(H_D3G33)%d = u(H_D3G33)%d	 

         u2(H_K00)%d = u(H_K00)%d	 
         u2(H_K01)%d = u(H_K01)%d	 
         u2(H_K02)%d = u(H_K02)%d	 
         u2(H_K03)%d = u(H_K03)%d
         u2(H_K11)%d = u(H_K11)%d	 
         u2(H_K12)%d = u(H_K12)%d	 
         u2(H_K13)%d = u(H_K13)%d
         u2(H_K22)%d = u(H_K22)%d
         u2(H_K23)%d = u(H_K23)%d         
         u2(H_K33)%d = u(H_K33)%d	 
	 
         u2(H_PHIR)%d = u(H_PHIR)%d	 
         u2(H_PHIC)%d = u(H_PHIC)%d	 
         u2(H_PHIM)%d = u(H_PHIM)%d

         u2(H_D1PHIR)%d = u(H_D1PHIR)%d	 
         u2(H_D2PHIR)%d = u(H_D2PHIR)%d	 
         u2(H_D3PHIR)%d = u(H_D3PHIR)%d

         u2(H_D1PHIC)%d = u(H_D1PHIC)%d	 
         u2(H_D2PHIC)%d = u(H_D2PHIC)%d	 
         u2(H_D3PHIC)%d = u(H_D3PHIC)%d

         u2(H_D1PHIM)%d = u(H_D1PHIM)%d	 
         u2(H_D2PHIM)%d = u(H_D2PHIM)%d	 
         u2(H_D3PHIM)%d = u(H_D3PHIM)%d
	 	 
         u2(H_PIR)%d = u(H_PIR)%d	 
         u2(H_PIC)%d = u(H_PIC)%d	 
         u2(H_PIM)%d = u(H_PIM)%d

!-----LOAD THE SECOND BOSON--------------------------------------      
!       the center of the boson stars, symmetric         
         xcenter = par(P_BS2_X0)
  	 ycenter = par(P_BS2_Y0)	
	 zcenter = par(P_BS2_Z0)
         omega1  = par(P_SF_OMEGA2)
         vx      = par(P_BS2_VX)
         vy      = par(P_BS2_VY)               
         vz      = par(P_BS2_VZ)

         if (boost_type .EQ. 0 ) then
           g00_bg  = -1.0d0 + par(P_BS1_VX)*par(P_BS1_VX) + par(P_BS1_VY)*par(P_BS1_VY) + &
                      par(P_BS1_VZ)*par(P_BS1_VZ) + par(P_BS2_VX)*par(P_BS2_VX) + &
                      par(P_BS2_VY)*par(P_BS2_VY) + par(P_BS2_VZ)*par(P_BS2_VZ)
           g01_bg  = - par(P_BS1_VX) - par(P_BS2_VX)     
           g02_bg  = - par(P_BS1_VY) - par(P_BS2_VY)
           g03_bg  = - par(P_BS1_VZ) - par(P_BS2_VZ)
         else
           g00_bg = -1.0d0; g01_bg = 0.0d0; g02_bg = 0.0d0; g03_bg = 0.0d0
         end if  
	
		    		 
!----- G, K, PHIR, PHIC ----------------------------------------
         !  the intrinsic metric g_{ij}
         call read_sdf_file1d(infile(1,1), u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

         u(H_G11)%d =(u(H_K11)%d)**4
         u(H_G22)%d = u(H_G11)%d
         u(H_G33)%d = u(H_G11)%d
         u(H_G12)%d = 0.0d0
         u(H_G13)%d = 0.0d0
         u(H_G23)%d = 0.0d0      
         u(H_K11)%d = 0.0d0

         !  the gauge g_{0a}
         call read_sdf_file1d(infile(1,2), u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
		      
         u(H_G00)%d =-(u(H_K11)%d)**2
         u(H_G01)%d = 0.0d0
         u(H_G02)%d = 0.0d0
         u(H_G03)%d = 0.0d0
         u(H_K11)%d = 0.0d0

         !  the complex scalar field
         call read_sdf_file1d(infile(1,3), u(H_K11)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)
		   
         u(H_PIC)%d = 0.0d0
         u(H_PIR)%d = 0.0d0
         u(H_PIM)%d = 0.0d0       

         u(H_PHIR)%d = u(H_K11)%d
         u(H_PHIC)%d = 0.0d0
         u(H_PHIM)%d = 0.0d0
	 
         !  the time derivative of the metric       
         u(H_K11)%d = 0.0d0
         u(H_K22)%d = 0.0d0
         u(H_K33)%d = 0.0d0
         u(H_K12)%d = 0.0d0
         u(H_K13)%d = 0.0d0
         u(H_K23)%d = 0.0d0      

         ! the first order derivates of the metric
         call default_data_derivs(u, w, par)      

       ! make the boost
       do i = 1, nx; do j = 1, ny; do k = 1, nz
#include "boost.inc"               
       end do; end do; end do      	 
	 	 	 
!--------COMPUTING THE SUPERPOSITION-----------------------         
         u(H_G00)%d = u2(H_G00)%d + u(H_G00)%d - g00_bg
         u(H_G01)%d = u2(H_G01)%d + u(H_G01)%d - g01_bg        
         u(H_G02)%d = u2(H_G02)%d + u(H_G02)%d - g02_bg
         u(H_G03)%d = u2(H_G03)%d + u(H_G03)%d - g03_bg	 
         u(H_G11)%d = u2(H_G11)%d + u(H_G11)%d - 1.0d0
         u(H_G12)%d = u2(H_G12)%d + u(H_G12)%d	 
         u(H_G13)%d = u2(H_G13)%d + u(H_G13)%d
         u(H_G22)%d = u2(H_G22)%d + u(H_G22)%d - 1.0d0
         u(H_G23)%d = u2(H_G23)%d + u(H_G23)%d         
         u(H_G33)%d = u2(H_G33)%d + u(H_G33)%d - 1.0d0	 
	 
         u(H_D1G00)%d = u2(H_D1G00)%d + u(H_D1G00)%d
         u(H_D1G01)%d = u2(H_D1G01)%d + u(H_D1G01)%d         
         u(H_D1G02)%d = u2(H_D1G02)%d + u(H_D1G02)%d	 
         u(H_D1G03)%d = u2(H_D1G03)%d + u(H_D1G03)%d	 
         u(H_D1G11)%d = u2(H_D1G11)%d + u(H_D1G11)%d	 
         u(H_D1G12)%d = u2(H_D1G12)%d + u(H_D1G12)%d	 
         u(H_D1G13)%d = u2(H_D1G13)%d + u(H_D1G13)%d
         u(H_D1G22)%d = u2(H_D1G22)%d + u(H_D1G22)%d
         u(H_D1G23)%d = u2(H_D1G23)%d + u(H_D1G23)%d         
         u(H_D1G33)%d = u2(H_D1G33)%d + u(H_D1G33)%d	 

         u(H_D2G00)%d = u2(H_D2G00)%d + u(H_D2G00)%d
         u(H_D2G01)%d = u2(H_D2G01)%d + u(H_D2G01)%d         
         u(H_D2G02)%d = u2(H_D2G02)%d + u(H_D2G02)%d	 
         u(H_D2G03)%d = u2(H_D2G03)%d + u(H_D2G03)%d	 
         u(H_D2G11)%d = u2(H_D2G11)%d + u(H_D2G11)%d	 
         u(H_D2G12)%d = u2(H_D2G12)%d + u(H_D2G12)%d	 
         u(H_D2G13)%d = u2(H_D2G13)%d + u(H_D2G13)%d
         u(H_D2G22)%d = u2(H_D2G22)%d + u(H_D2G22)%d
         u(H_D2G23)%d = u2(H_D2G23)%d + u(H_D2G23)%d         
         u(H_D2G33)%d = u2(H_D2G33)%d + u(H_D2G33)%d	 
	 
         u(H_D3G00)%d = u2(H_D3G00)%d + u(H_D3G00)%d
         u(H_D3G01)%d = u2(H_D3G01)%d + u(H_D3G01)%d         
         u(H_D3G02)%d = u2(H_D3G02)%d + u(H_D3G02)%d	 
         u(H_D3G03)%d = u2(H_D3G03)%d + u(H_D3G03)%d	 
         u(H_D3G11)%d = u2(H_D3G11)%d + u(H_D3G11)%d	 
         u(H_D3G12)%d = u2(H_D3G12)%d + u(H_D3G12)%d	 
         u(H_D3G13)%d = u2(H_D3G13)%d + u(H_D3G13)%d
         u(H_D3G22)%d = u2(H_D3G22)%d + u(H_D3G22)%d
         u(H_D3G23)%d = u2(H_D3G23)%d + u(H_D3G23)%d         
         u(H_D3G33)%d = u2(H_D3G33)%d + u(H_D3G33)%d	 

         u(H_K00)%d = u2(H_K00)%d + u(H_K00)%d
         u(H_K01)%d = u2(H_K01)%d + u(H_K01)%d	 
         u(H_K02)%d = u2(H_K02)%d + u(H_K02)%d	 
         u(H_K03)%d = u2(H_K03)%d + u(H_K03)%d
         u(H_K11)%d = u2(H_K11)%d + u(H_K11)%d	 
         u(H_K12)%d = u2(H_K12)%d + u(H_K12)%d	 
         u(H_K13)%d = u2(H_K13)%d + u(H_K13)%d
         u(H_K22)%d = u2(H_K22)%d + u(H_K22)%d
         u(H_K23)%d = u2(H_K23)%d + u(H_K23)%d         
         u(H_K33)%d = u2(H_K33)%d + u(H_K33)%d	 
	 
         u(H_PHIR)%d = u2(H_PHIR)%d + u(H_PHIR)%d
         u(H_PHIC)%d = u2(H_PHIC)%d + u(H_PHIC)%d	 
         u(H_PHIM)%d = u2(H_PHIM)%d + u(H_PHIM)%d

         u(H_D1PHIR)%d = u2(H_D1PHIR)%d	+ u(H_D1PHIR)%d 
         u(H_D2PHIR)%d = u2(H_D2PHIR)%d	+ u(H_D2PHIR)%d 	 
         u(H_D3PHIR)%d = u2(H_D3PHIR)%d	+ u(H_D3PHIR)%d 

         u(H_D1PHIC)%d = u2(H_D1PHIC)%d	+ u(H_D1PHIC)%d 
         u(H_D2PHIC)%d = u2(H_D2PHIC)%d	+ u(H_D2PHIC)%d 	 
         u(H_D3PHIC)%d = u2(H_D3PHIC)%d	+ u(H_D3PHIC)%d 

         u(H_D1PHIM)%d = u2(H_D1PHIM)%d	+ u(H_D1PHIM)%d 
         u(H_D2PHIM)%d = u2(H_D2PHIM)%d	+ u(H_D2PHIM)%d 	 
         u(H_D3PHIM)%d = u2(H_D3PHIM)%d	+ u(H_D3PHIM)%d 
	 	 
         u(H_PIR)%d = u2(H_PIR)%d + u(H_PIR)%d
         u(H_PIC)%d = u2(H_PIC)%d + u(H_PIC)%d	 
         u(H_PIM)%d = u2(H_PIM)%d + u(H_PIM)%d
	 
       end if
      
      
      else
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!   otherwise, the standard shit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     
      lev = nint(par(P_READ_DATA_LEVEL))
      interp = nint(par(P_INTERP_ID))

      if (readfromfile(17) .eq. 1) then
        call sdf_file_mem(filename(17), lev, shp)
      else	
        call sdf_file_mem(filename(1), lev, shp)
      end if	
      
      
      allocate(fdata(shp(1),shp(2),shp(3)),&
               coords(shp(1)+shp(2)+shp(3)),&
               STAT=istat)
      if (istat .ne. 0) then
        write(*,*)'DEF_INITIAL:  Could not allocate memory'
        write(*,*)'for arrays to read sdf files'
        stop
      end if

       ! it is most useful to have 1-d coordinate arrays.
       allocate (x1d(nx), y1d(ny), z1d(nz), STAT=istat)
       if (istat .ne. 0) then
         write(*,*)'DEF_INITIAL:  Could not allocate memory'
         write(*,*)'for 1d coordinate arrays'
         stop
       end if
       do i = 1, nx
          x1d(i) = xx(i,1,1)
       end do
       do i = 1, ny
          y1d(i) = yy(1,i,1)
       end do
       do i = 1, nz
          z1d(i) = zz(1,1,i)
       end do

      minx0 = par(P_MINX0)
      miny0 = par(P_MINY0)
      minz0 = par(P_MINZ0)

      hi   = x1d(2) - x1d(1)
      !
      ! Need minimums for this grid alone:
      !
      minx = x1d(1)
      miny = y1d(1)
      minz = z1d(1)

      !
      ! Added by sll:
      !
      lb(1) = NINT( (minx - minx0) / hi )
      lb(2) = NINT( (miny - miny0) / hi )
      lb(3) = NINT( (minz - minz0) / hi )

      gs(1) = nint(par(P_GLOBAL_NX))
      gs(2) = nint(par(P_GLOBAL_NY))
      gs(3) = nint(par(P_GLOBAL_NZ))

      if (ltrace) then
         write(*,*)'@@@ FILE_DATA_G:  reading metric from files'
         write(*,*)'@@@ :  local array size ',nx,ny,nz
         write(*,*)'@@@ :  global array size',gs(1),gs(2),gs(3)
         write(*,*)'@@@ :  lower bounds     ',lb(1),lb(2),lb(3)
         write(*,*)'@@@ :  shapes           ',shp(1),shp(2),shp(3)	 
      end if

      dnx = shp(1)
      dny = shp(2)
      dnz = shp(3)
      
!     first let us set the metric and the extrinsic curvature either from a 
!     conformal factor or for a general metric      
      ifc =17  ! the conformal factor psi
      if (readfromfile(ifc) .eq. 1) then
          if (ltrace) write(*,*)'reading from a conformal factor!!!!'
          call read_sdf_file(filename(ifc), u(fileindex(ifc))%d,   &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata, coords,    &
                   dnx, dny, dnz, interp, lb, gs)
          u(H_G11)%d = (u(fileindex(ifc))%d)**4
	  u(H_G22)%d = u(H_G11)%d
          u(H_G33)%d = u(H_G11)%d
          u(H_G12)%d = 0.0d0
          u(H_G13)%d = 0.0d0
          u(H_G23)%d = 0.0d0      
          u(H_K11)%d = 0.0d0
          u(H_K22)%d = 0.0d0
          u(H_K33)%d = 0.0d0
          u(H_K12)%d = 0.0d0
          u(H_K13)%d = 0.0d0
          u(H_K23)%d = 0.0d0      
          if (ltrace) then
	    do ifc=1,12
	      write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
	    end do 
	  end if        
      else 
        if (ltrace) write(*,*)'reading the physical metric gij and Kij !!!!'      
	do ifc = 1, 12
          if (readfromfile(ifc) .eq. 1) then
            if (ltrace) write(*,*)'readfromfile = ',filename(ifc)
            call read_sdf_file(filename(ifc), u(fileindex(ifc))%d,   &
                     x1d, y1d, z1d, nx, ny, nz, lev, fdata, coords,    &
                     dnx, dny, dnz, interp, lb, gs)
          else
              u(fileindex(ifc))%d = funcdefault(ifc)
          end if
          if (ltrace) then
	    write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
	  end if  
        end do
      end if
                  
!     then set the gauge (lapse and shift) in an independent way      
      do ifc = 13, 16
        if (readfromfile(ifc) .eq. 1) then
          if (ltrace) write(*,*)'readfromfile = ',filename(ifc)  
          call read_sdf_file(filename(ifc), u(fileindex(ifc))%d,   &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata, coords,    &
                   dnx, dny, dnz, interp, lb, gs)
        else
            u(fileindex(ifc))%d = funcdefault(ifc)
        end if
        if (ltrace) then
          write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
	end if  
      end do

!     get the lapse for the boson star, that will be used for the pi      
      if (idtype .EQ. -2) then        
          call read_sdf_file("alpha_boson.sdf", u(H_K00)%d,   &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata, coords,    &
                   dnx, dny, dnz, interp, lb, gs)      
      end if
      
      
!     then set the scalar fields
      do ifc = 18, 23
        if (readfromfile(ifc) .eq. 1) then
          if (ltrace) write(*,*)'readfromfile = ',filename(ifc)  
          call read_sdf_file(filename(ifc), u(fileindex(ifc))%d,   &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata, coords,    &
                   dnx, dny, dnz, interp, lb, gs)
        else
            u(fileindex(ifc))%d = funcdefault(ifc)
        end if
        if (ltrace) then
          write(*,*) '@@@ file:  ||f|| = ',ifc, myl2norm3d(u(fileindex(ifc))%d, nx, ny, nz)
        end if  
      end do

      ! setting the the first order derivates of the metric
       call default_data_derivs(u, w, par)            

      end if      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!   more of the temporal crap   !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      if ((idtype .LE. -6)) then              
        !write(*,*) 'deallocating coords1d'
        if (allocated(coords1d)) deallocate(coords1d)
        if (allocated(fdata1d) ) deallocate(fdata1d)
        if (allocated(x1d)     ) deallocate(x1d)
        if (allocated(y1d)     ) deallocate(y1d)
        if (allocated(z1d)     ) deallocate(z1d)
        !deallocate(coords1d, fdata1d, x1d, y1d, z1d)
      else
        !write(*,*) 'deallocating coords'
        if (allocated(coords)  ) deallocate(coords)
        if (allocated(fdata)   ) deallocate(fdata)
        if (allocated(x1d)     ) deallocate(x1d)
        if (allocated(y1d)     ) deallocate(y1d)
        if (allocated(z1d)     ) deallocate(z1d)
        !deallocate(coords, fdata, x1d, y1d, z1d)
      end if

      
      return
    end subroutine file_data_gr

    !----------------------------------------------------------------
    !
    !  Calculate the ds as numerical spatial derivatives of the metric
    !
    !_________________________________________________________________
    subroutine default_data_derivs(u, w, par)

      implicit none
      type(gridfunction), dimension(NU_G)     :: u
      type(gridfunction), dimension(NW)       :: w
      CCTK_REAL, dimension(:)                 :: par

      interface
        subroutine D2_NoMask(Du, u, direction, par)
          use params
          implicit none
          CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
          CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
          CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
          CCTK_INT :: direction
        end subroutine D2_NoMask
	
        subroutine D42NoMask_simple(Du, u, direction, par)
          use params
          implicit none
          CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
          CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
          CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
          CCTK_INT :: direction
        end subroutine D42NoMask_simple
      end interface

      if ((nint(par(P_DERIV_ORDER)) .eq. 4) .OR. &
          (nint(par(P_DERIV_ORDER)) .eq. 42)) then

! 	print*, "4th order derivative in the initial data"  
	  
 	call D42NoMask_simple(u(H_D1G00)%d, u(H_G00)%d, 1, par)
        call D42NoMask_simple(u(H_D2G00)%d, u(H_G00)%d, 2, par)
        call D42NoMask_simple(u(H_D3G00)%d, u(H_G00)%d, 3, par)
 	call D42NoMask_simple(u(H_D1G01)%d, u(H_G01)%d, 1, par)
        call D42NoMask_simple(u(H_D2G01)%d, u(H_G01)%d, 2, par)
        call D42NoMask_simple(u(H_D3G01)%d, u(H_G01)%d, 3, par)
 	call D42NoMask_simple(u(H_D1G02)%d, u(H_G02)%d, 1, par)
        call D42NoMask_simple(u(H_D2G02)%d, u(H_G02)%d, 2, par)
        call D42NoMask_simple(u(H_D3G02)%d, u(H_G02)%d, 3, par)
 	call D42NoMask_simple(u(H_D1G03)%d, u(H_G03)%d, 1, par)
        call D42NoMask_simple(u(H_D2G03)%d, u(H_G03)%d, 2, par)
        call D42NoMask_simple(u(H_D3G03)%d, u(H_G03)%d, 3, par)

	call D42NoMask_simple(u(H_D1G11)%d, u(H_G11)%d, 1, par)
        call D42NoMask_simple(u(H_D2G11)%d, u(H_G11)%d, 2, par)
        call D42NoMask_simple(u(H_D3G11)%d, u(H_G11)%d, 3, par)
        call D42NoMask_simple(u(H_D1G12)%d, u(H_G12)%d, 1, par)
        call D42NoMask_simple(u(H_D2G12)%d, u(H_G12)%d, 2, par)
        call D42NoMask_simple(u(H_D3G12)%d, u(H_G12)%d, 3, par)
        call D42NoMask_simple(u(H_D1G13)%d, u(H_G13)%d, 1, par)
        call D42NoMask_simple(u(H_D2G13)%d, u(H_G13)%d, 2, par)
        call D42NoMask_simple(u(H_D3G13)%d, u(H_G13)%d, 3, par)
        call D42NoMask_simple(u(H_D1G22)%d, u(H_G22)%d, 1, par)
        call D42NoMask_simple(u(H_D2G22)%d, u(H_G22)%d, 2, par)
        call D42NoMask_simple(u(H_D3G22)%d, u(H_G22)%d, 3, par)
        call D42NoMask_simple(u(H_D1G23)%d, u(H_G23)%d, 1, par)
        call D42NoMask_simple(u(H_D2G23)%d, u(H_G23)%d, 2, par)
        call D42NoMask_simple(u(H_D3G23)%d, u(H_G23)%d, 3, par)
        call D42NoMask_simple(u(H_D1G33)%d, u(H_G33)%d, 1, par)
        call D42NoMask_simple(u(H_D2G33)%d, u(H_G33)%d, 2, par)
        call D42NoMask_simple(u(H_D3G33)%d, u(H_G33)%d, 3, par)

        call D42NoMask_simple(u(H_D1PHIR)%d, u(H_PHIR)%d, 1, par)
        call D42NoMask_simple(u(H_D2PHIR)%d, u(H_PHIR)%d, 2, par)
        call D42NoMask_simple(u(H_D3PHIR)%d, u(H_PHIR)%d, 3, par)
        call D42NoMask_simple(u(H_D1PHIC)%d, u(H_PHIC)%d, 1, par)
        call D42NoMask_simple(u(H_D2PHIC)%d, u(H_PHIC)%d, 2, par)
        call D42NoMask_simple(u(H_D3PHIC)%d, u(H_PHIC)%d, 3, par)

!        call D42NoMask_simple(u(H_D1PHIM)%d, u(H_PHIM)%d, 1, par)
!        call D42NoMask_simple(u(H_D2PHIM)%d, u(H_PHIM)%d, 2, par)
!        call D42NoMask_simple(u(H_D3PHIM)%d, u(H_PHIM)%d, 3, par)
!        call D42NoMask_simple(u(H_D1PHIN)%d, u(H_PHIN)%d, 1, par)
!        call D42NoMask_simple(u(H_D2PHIN)%d, u(H_PHIN)%d, 2, par)
!        call D42NoMask_simple(u(H_D3PHIN)%d, u(H_PHIN)%d, 3, par)

      else

!	print*, "2nd order derivative in the initial data"  
	 
	call D2_NoMask(u(H_D1G00)%d, u(H_G00)%d, 1, par)
        call D2_NoMask(u(H_D2G00)%d, u(H_G00)%d, 2, par)
        call D2_NoMask(u(H_D3G00)%d, u(H_G00)%d, 3, par)
 	call D2_NoMask(u(H_D1G01)%d, u(H_G01)%d, 1, par)
        call D2_NoMask(u(H_D2G01)%d, u(H_G01)%d, 2, par)
        call D2_NoMask(u(H_D3G01)%d, u(H_G01)%d, 3, par)
 	call D2_NoMask(u(H_D1G02)%d, u(H_G02)%d, 1, par)
        call D2_NoMask(u(H_D2G02)%d, u(H_G02)%d, 2, par)
        call D2_NoMask(u(H_D3G02)%d, u(H_G02)%d, 3, par)
 	call D2_NoMask(u(H_D1G03)%d, u(H_G03)%d, 1, par)
        call D2_NoMask(u(H_D2G03)%d, u(H_G03)%d, 2, par)
        call D2_NoMask(u(H_D3G03)%d, u(H_G03)%d, 3, par)
       	
	call D2_NoMask(u(H_D1G11)%d, u(H_G11)%d, 1, par)
        call D2_NoMask(u(H_D2G11)%d, u(H_G11)%d, 2, par)
        call D2_NoMask(u(H_D3G11)%d, u(H_G11)%d, 3, par)
        call D2_NoMask(u(H_D1G12)%d, u(H_G12)%d, 1, par)
        call D2_NoMask(u(H_D2G12)%d, u(H_G12)%d, 2, par)
        call D2_NoMask(u(H_D3G12)%d, u(H_G12)%d, 3, par)
        call D2_NoMask(u(H_D1G13)%d, u(H_G13)%d, 1, par)
        call D2_NoMask(u(H_D2G13)%d, u(H_G13)%d, 2, par)
        call D2_NoMask(u(H_D3G13)%d, u(H_G13)%d, 3, par)
        call D2_NoMask(u(H_D1G22)%d, u(H_G22)%d, 1, par)
        call D2_NoMask(u(H_D2G22)%d, u(H_G22)%d, 2, par)
        call D2_NoMask(u(H_D3G22)%d, u(H_G22)%d, 3, par)
        call D2_NoMask(u(H_D1G23)%d, u(H_G23)%d, 1, par)
        call D2_NoMask(u(H_D2G23)%d, u(H_G23)%d, 2, par)
        call D2_NoMask(u(H_D3G23)%d, u(H_G23)%d, 3, par)
        call D2_NoMask(u(H_D1G33)%d, u(H_G33)%d, 1, par)
        call D2_NoMask(u(H_D2G33)%d, u(H_G33)%d, 2, par)
        call D2_NoMask(u(H_D3G33)%d, u(H_G33)%d, 3, par)

        call D2_NoMask(u(H_D1PHIR)%d, u(H_PHIR)%d, 1, par)
        call D2_NoMask(u(H_D2PHIR)%d, u(H_PHIR)%d, 2, par)
        call D2_NoMask(u(H_D3PHIR)%d, u(H_PHIR)%d, 3, par)	        
        call D2_NoMask(u(H_D1PHIC)%d, u(H_PHIC)%d, 1, par)
        call D2_NoMask(u(H_D2PHIC)%d, u(H_PHIC)%d, 2, par)
        call D2_NoMask(u(H_D3PHIC)%d, u(H_PHIC)%d, 3, par)	        	

!        call D2_NoMask(u(H_D1PHIM)%d, u(H_PHIM)%d, 1, par)
!        call D2_NoMask(u(H_D2PHIM)%d, u(H_PHIM)%d, 2, par)
!        call D2_NoMask(u(H_D3PHIM)%d, u(H_PHIM)%d, 3, par)	        
!        call D2_NoMask(u(H_D1PHIN)%d, u(H_PHIN)%d, 1, par)
!        call D2_NoMask(u(H_D2PHIN)%d, u(H_PHIN)%d, 2, par)
!        call D2_NoMask(u(H_D3PHIN)%d, u(H_PHIN)%d, 3, par)	        	
	
      end if

    end subroutine default_data_derivs

    !----------------------------------------------------------------
    !
    !  Calculate the ds as numerical spatial derivatives of the metric
    !
    !_________________________________________________________________
    subroutine boost_vec(vnew, vold, Jud)
      implicit none
      real*8, dimension(0:3)     :: vnew, vold
      real*8, dimension(0:3,0:3) :: Jud
      ! local vars
      integer                    :: i, k

      do i = 0, 3
         vnew(i) = 0.0d0
         do k = 0, 3
            vnew(i) = vnew(i) + vold(k)*Jud(k,i)
         end do
      end do

      return
    end subroutine boost_vec

    !----------------------------------------------------------------
    !
    !  Calculate the ds as numerical spatial derivatives of the metric
    !
    !_________________________________________________________________
    subroutine boost_metric(gnew, gold, Jud)
      implicit none
      real*8, dimension(0:3,0:3) :: gnew, gold
      real*8, dimension(0:3,0:3) :: Jud
      ! local vars
      integer                    :: i, j, k, m

      do i = 0, 3
      do j = 0, 3
         gnew(i,j) = 0.0d0
         do k = 0, 3
         do m = 0, 3
            gnew(i,j) = gnew(i,j) + gold(k,m)*Jud(k,i)*Jud(m,j)
         end do
         end do
      end do
      end do

      return
    end subroutine boost_metric

    !----------------------------------------------------------------
    !
    !  Calculate the ds as numerical spatial derivatives of the metric
    !
    !_________________________________________________________________
    subroutine boost_metric_derivs(Dgnew, Dgold, gold, Jud, DJud)
      implicit none
      real*8, dimension(0:3,0:3,0:3) :: Dgnew, Dgold, DJud
      real*8, dimension(0:3,0:3) :: Jud, gold
      ! local vars
      integer                    :: li, lj, ll, lk, lm, lp
      real*8                     :: Temp1, Temp2

      do li = 0, 3
      do lj = 0, 3
      do ll = 0, 3
         Temp1 = 0.0d0
         Temp2 = 0.0d0
         do lk = 0, 3
         do lm = 0, 3
         do lp = 0, 3
            Temp1 =  Temp1   &
             &   + Dgold(lp,lk,lm) * Jud(lp,ll) * Jud(lk,li) * Jud(lm,lj)
         end do
         end do
         end do
         do lk = 0, 3
         do lm = 0, 3
            Temp2 =  Temp2  &
             &    + gold(lk,lm) * DJud(lk,li,ll) * Jud(lm,lj) &
             &    + gold(lk,lm) * Jud(lk,li) * DJud(lm,lj,ll)
         end do
         end do
         Dgnew(ll,li,lj) = Temp1 + Temp2
      end do
      end do
      end do

      return
    end subroutine boost_metric_derivs



 
end
