c-----------------------------------------------------------------------
c  This routine is called by level_lib to interface with Mijan's
c  horizon finder
c-----------------------------------------------------------------------
        subroutine mijans_horizon(level)
        implicit none
        integer  gi
        include     'grid_methods.inc'
        include     'grid.inc'
        include     'glob.inc'
        include     'mask.inc'
        include     'mpif.h'
        include     'mpi_stuff.inc'
        include     'param.inc'
        integer level,i,j,k,flag
        real*8, pointer, dimension(:,:,:) :: eg11
        real*8, pointer, dimension(:,:,:) :: eg12
        real*8, pointer, dimension(:,:,:) :: eg13
        real*8, pointer, dimension(:,:,:) :: eg22
        real*8, pointer, dimension(:,:,:) :: eg23
        real*8, pointer, dimension(:,:,:) :: eg33
        real*8, pointer, dimension(:,:,:) :: ek11
        real*8, pointer, dimension(:,:,:) :: ek12
        real*8, pointer, dimension(:,:,:) :: ek13
        real*8, pointer, dimension(:,:,:) :: ek22
        real*8, pointer, dimension(:,:,:) :: ek23
        real*8, pointer, dimension(:,:,:) :: ek33
        real*8, pointer, dimension(:,:,:) :: msk
        real*8, pointer, dimension(:) :: xx
        real*8, pointer, dimension(:) :: yy
        real*8, pointer, dimension(:) :: zz
        real(kind=8) hl
        real(kind=8) minx, maxx, miny, maxy, minz, maxz
        integer      nxl,nyl,nzl
        integer      readflag
        real(kind=8) bh_center(3)
        real(kind=8) horizon_area,horizon_ct,horizon_cp
        real(kind=8) Pi,polar_radius,equatorial_radius
        !
        character(5) fieldname
        logical      ltrace
        parameter (  ltrace = .false. )
        real(kind=8) SHRINK         ! How much smaller is region to excise?
        parameter  ( SHRINK = 0.75d0 )


 99   format('[',I3,'] ',A,3I5)

        Pi = 4.0d0*atan(1.0d0)

c find the bounds of the level being used for horizon finding
        call level_find_bounds(level,minx,maxx,miny,maxy,minz,maxz)

        gi = level_return_start(level)
        hl    = grid_return_h(gi)
        nxl   = NINT( (maxx-minx)/hl ) + 1
        nyl   = NINT( (maxy-miny)/hl ) + 1
        nzl   = NINT( (maxz-minz)/hl ) + 1

      if (ltrace) then
         write(*,99) myid, 'mijans_horizon: Called on level=',level
         write(*,99) myid, 'mijans_horizon: nx/y/zl:',nxl,nyl,nzl
      end if

c allocate memory big enough to hold what we need
        if (myid.eq.master) then
          if(ltrace)write(*,99)myid,'mijans_horizon: Allocating memoroy'
          allocate(eg11(nxl,nyl,nzl))
          allocate(eg12(nxl,nyl,nzl))
          allocate(eg13(nxl,nyl,nzl))
          allocate(eg22(nxl,nyl,nzl))
          allocate(eg23(nxl,nyl,nzl))
          allocate(eg33(nxl,nyl,nzl))
          allocate(ek11(nxl,nyl,nzl))
          allocate(ek12(nxl,nyl,nzl))
          allocate(ek13(nxl,nyl,nzl))
          allocate(ek22(nxl,nyl,nzl))
          allocate(ek23(nxl,nyl,nzl))
          allocate(ek33(nxl,nyl,nzl))
          allocate(msk(nxl,nyl,nzl))
          allocate(xx(nxl))
          allocate(yy(nyl))
          allocate(zz(nzl))
        end if

c   compute the ADM extrinsic curvature on the finest mesh
        gi = level_return_start(level)
 10     if (grid_return_existence(gi)) then     
          if (grid_is_local(gi)) then
            if(ltrace)write(*,*)myid,'    Computing on grid: ',gi
            call compute_ADM(gi)
          end if 
          gi = grid_return_sibling(gi)
          goto 10
        end if
        
c       gather the adm variables from the finest mesh onto 
c       the master processor
        if (ltrace) write(*,*) myid, 'mijans_horizon: Gathering'
c the fieldnumber is found in had/include/fields.inc
        call level_gather_field(level,9 ,eg11,nxl,nyl,nzl) 
        call level_gather_field(level,11,eg12,nxl,nyl,nzl) 
        call level_gather_field(level,13,eg13,nxl,nyl,nzl) 
        call level_gather_field(level,15,eg22,nxl,nyl,nzl) 
        call level_gather_field(level,17,eg23,nxl,nyl,nzl) 
        call level_gather_field(level,19,eg33,nxl,nyl,nzl) 
        call level_gather_field(level,526,ek11,nxl,nyl,nzl) 
        call level_gather_field(level,527,ek12,nxl,nyl,nzl) 
        call level_gather_field(level,528,ek13,nxl,nyl,nzl) 
        call level_gather_field(level,529,ek22,nxl,nyl,nzl) 
        call level_gather_field(level,530,ek23,nxl,nyl,nzl) 
        call level_gather_field(level,531,ek33,nxl,nyl,nzl) 
        call level_gather_field(level,540,msk,nxl,nyl,nzl) 
        if (ltrace) write(*,*) myid, 'mijans_horizon: Done Gathering'

        if (myid.eq.master) then
c  initialize the mask grid and coordinate arrays
          do k=1,nzl
          do j=1,nyl
          do i=1,nxl
c  Mijan's finder users different definitions for the mask
            if ( msk(i,j,k) .le. 0 ) then
              msk(i,j,k) = 0.0
            else if ( msk(i,j,k) .gt. 0 ) then 
              msk(i,j,k) = 1.0
            else
              msk(i,j,k) = 0.0
            end if
          end do
          end do
          end do
          do i=1,nxl
            xx(i) = minx + hl*(i-1)
          end do
          do j=1,nyl
            yy(j) = miny + hl*(j-1)
          end do
          do k=1,nzl
            zz(k) = minz + hl*(k-1)
          end do
          ! estimate for the location for one hole
          if ( local_time < 1.d-7 ) then
            bh_center(1) = bh1_x0
            bh_center(2) = bh1_y0
            bh_center(3) = bh1_z0
          else
            bh_center(1) = mask_center(1,1)
            bh_center(2) = mask_center(2,1)
            bh_center(3) = mask_center(3,1)
          end if
c  call the horizon finder
          if ( bh_true(1) ) then
c  read in the horizon from the previous successful search
            readflag = 0
          else 
            readflag = 0
          end if
          call mijan_horizon_interface(eg11,eg12,eg13,
     &                                 eg22,eg23,eg33,
     &                                 ek11,ek12,ek13,
     &                                 ek22,ek23,ek33,
     &                                 xx,yy,zz, 
     &                                 msk,nxl,nyl,nzl,
     &                                 horizon_ntheta,horizon_nphi,flag,
     &                                 bh_center,horizon_area,
     &                                 horizon_ct,horizon_cp,
     &                                 bh1_mass,2.0d0*bh1_mass,readflag)
c write it out to an sdf file from the master processor --debug
          if (ltrace) call field_out3d(eg11,0.0,"eg11",
     &             0.0,1.0,0.0,1.0,0.0,1.0,
     &            nxl,nyl,nzl,0)

          if ( flag .eq. 1 ) then
            bh_true(1) = .true.
c            bh_radius = sqrt(horizon_area/(4*Pi))
            polar_radius = horizon_ct/(2*Pi)
            equatorial_radius = horizon_cp/(2*Pi)
            mask_bbox(1,1)= bh_center(1) 
     &               - SHRINK*equatorial_radius/sqrt(2.0d0)
            mask_bbox(3,1)= bh_center(2)
     &               - SHRINK*equatorial_radius/sqrt(2.0d0)
            mask_bbox(5,1)= bh_center(3) 
     &               - SHRINK*polar_radius/sqrt(2.0d0)
            !
            mask_bbox(2,1)= bh_center(1)
     &               + SHRINK*equatorial_radius/sqrt(2.0d0)
            mask_bbox(4,1)= bh_center(2) 
     &               + SHRINK*equatorial_radius/sqrt(2.0d0)
            mask_bbox(6,1)= bh_center(3) 
     &               + SHRINK*polar_radius/sqrt(2.0d0)

            mask_center(1,1) = bh_center(1)
            mask_center(2,1) = bh_center(2)
            mask_center(3,1) = bh_center(3)

            mask_radius(1,1) = equatorial_radius
            mask_radius(2,1) = equatorial_radius
            mask_radius(3,1) = polar_radius
!          else
!            bh_true(1) = .false.
          end if

          deallocate(eg11)
          deallocate(eg12)
          deallocate(eg13)
          deallocate(eg22)
          deallocate(eg23)
          deallocate(eg33)
          deallocate(ek11)
          deallocate(ek12)
          deallocate(ek13)
          deallocate(ek22)
          deallocate(ek23)
          deallocate(ek33)
          deallocate(msk)
          deallocate(xx)
          deallocate(yy)
          deallocate(zz)
        end if

        call MPI_BCAST(bh_true, max_num_masks, MPI_LOGICAL,
     &              master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mask_bbox(1,1),6*max_num_masks,
     &              MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mask_center(1,1),3*max_num_masks,
     &              MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mask_radius(1,1),3*max_num_masks,
     &              MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)

        end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  computes the ADM extrinsic curvature
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_ADM(gi)
        implicit     none
        integer      gi
        include     'grid.inc'
        include     'mpif.h'
        include     'mpi_stuff.inc'
        include     'output.inc'

        integer      nx,ny,nz
        integer      i,j,k, index
        real(kind=8) hg
        real(kind=8) g11, g12, g13, g22, g23, g33, g00, g01, g02,
     *               g03, mask, d1g11, d1g12, d1g13,
     *               D100,D101,D102,D103,D111,D112,D113,D122,D123,D133,
     *               D200,D201,D202,D203,D211,D212,D213,D222,D223,D233,
     *               D300,D301,D302,D303,D311,D312,D313,D322,D323,D333,
     *               K11, K12, K13, K22, K23, K33

        real(kind=8) h11,h12,h13,h22,h23,h33,
     *               huu11,huu12,huu13,huu22,huu23,huu33,
     *               deth, b1, b2, b3, alp, alp2

        logical      ltrace
        parameter (  ltrace  = .false. )

        call load_pointers(gi)
        nx   = gr_nx(gi)
        ny   = gr_ny(gi)
        nz   = gr_nz(gi)
        hg   = gr_h(gi)
#ifdef GH
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          index = (k-1)*ny*nx + (j-1)*nx + (i-1)
          g00 = q(gr_g00+index)
          g01 = q(gr_g01+index)
          g02 = q(gr_g02+index)
          g03 = q(gr_g03+index)
          g11 = q(gr_g11+index)
          g12 = q(gr_g12+index)
          g13 = q(gr_g13+index)
          g22 = q(gr_g22+index)
          g23 = q(gr_g23+index)
          g33 = q(gr_g33+index)
          K11 = q(gr_K11+index)
          K12 = q(gr_K12+index)
          K13 = q(gr_K13+index)
          K22 = q(gr_K22+index)
          K23 = q(gr_K23+index)
          K33 = q(gr_K33+index)
          D100 = q(gr_d1g00+index)
          D101 = q(gr_d1g01+index)
          D102 = q(gr_d1g02+index)
          D103 = q(gr_d1g03+index)
          D111 = q(gr_d1g11+index)
          D112 = q(gr_d1g12+index)
          D113 = q(gr_d1g13+index)
          D122 = q(gr_d1g22+index)
          D123 = q(gr_d1g23+index)
          D133 = q(gr_d1g33+index)
          D200 = q(gr_d2g00+index)
          D201 = q(gr_d2g01+index)
          D202 = q(gr_d2g02+index)
          D203 = q(gr_d2g03+index)
          D211 = q(gr_d2g11+index)
          D212 = q(gr_d2g12+index)
          D213 = q(gr_d2g13+index)
          D222 = q(gr_d2g22+index)
          D223 = q(gr_d2g23+index)
          D233 = q(gr_d2g33+index)
          D300 = q(gr_d3g00+index)
          D301 = q(gr_d3g01+index)
          D302 = q(gr_d3g02+index)
          D303 = q(gr_d3g03+index)
          D311 = q(gr_d3g11+index)
          D312 = q(gr_d3g12+index)
          D313 = q(gr_d3g13+index)
          D322 = q(gr_d3g22+index)
          D323 = q(gr_d3g23+index)
          D333 = q(gr_d3g33+index)


          !
          !       compute the 3+1 quantities
          !
          h11 = g11
          h12 = g12
          h13 = g13
          h22 = g22
          h23 = g23
          h33 = g33
          !
          huu11 = -h23**2 + h22*h33
          huu12 = h13*h23 - h12*h33
          huu13 = -(h13*h22) + h12*h23
          huu22 = -h13**2 + h11*h33
          huu23 = h12*h13 - h11*h23
          huu33 = -h12**2 + h11*h22
          deth  = h11*huu11 + h12*huu12 + h13*huu13

          if(deth.lt.1e-6) then
            deth = 1.e-6
          end if
         
          huu11 = huu11/deth
          huu12 = huu12/deth
          huu13 = huu13/deth
          huu22 = huu22/deth
          huu23 = huu23/deth
          huu33 = huu33/deth

          !
          b1   =   g01*huu11 + g02*huu12 + g03*huu13
          b2   =   g01*huu12 + g02*huu22 + g03*huu23
          b3   =   g01*huu13 + g02*huu23 + g03*huu33
          alp2 = - g00 + b1**2 *h11 + b2**2*h22
     *     + 2*(b1*(b2*h12 + b3*h13)
     *     + b2*b3*h23) + b3**2*h33
          if(alp2.lt.1e-4) then
            alp2 = 1.e-4
          end if
          alp  = dsqrt(alp2)

          q(gr_ADM_k11+index) = -(-alp*K11 - D101 - D101
     *                        + b1*(D111 + D111) + b2*(D112 + D112)
     *                        + b3*(D113 + D113))/(2.0*alp)
          q(gr_ADM_k12+index) = -(-alp*K12 - D102 - D201
     *                        + b1*(D112 + D211) + b2*(D122 + D212)
     *                        + b3*(D123 + D213))/(2.0*alp)
          q(gr_ADM_k13+index) = -(-alp*K13 - D103 - D301
     *                        + b1*(D113 + D311) + b2*(D123 + D312)
     *                        + b3*(D133 + D313))/(2.0*alp)
          q(gr_ADM_k22+index) = -(-alp*K22 - D202 - D202
     *                        + b1*(D212 + D212) + b2*(D222 + D222)
     *                        + b3*(D223 + D223))/(2.0*alp)
          q(gr_ADM_k23+index) = -(-alp*K23 - D203 - D302
     *                        + b1*(D213 + D312) + b2*(D223 + D322)
     *                        + b3*(D233 + D323))/(2.0*alp)
          q(gr_ADM_k33+index) = -(-alp*K33 - D303 - D303
     *                        + b1*(D313 + D313) + b2*(D323 + D323)
     *                        + b3*(D333 + D333))/(2.0*alp)
        end do
        end do
        end do
#endif
      end subroutine

