cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  comp_amr_error:                                                           cc
cc                   Needs to set the field gr_error to have                  cc
cc                   the appropriate error function....                       cc
cc                   therefore needs to be non-negative function              cc
cc                   of the fields.                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_amr_error( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'grid_methods.inc'
      include    'glob.inc'
      include    'param.inc'
      integer     nx, ny, nz, numfields
      real*8      hg,  tmp, tmp_maxchi
      real*8      distance_wanted, distance_actual
      real*8      mask_distance_movedx
      external    mask_distance_movedx
      integer     tmp_maxchi_i, tmp_maxchi_j,tmp_maxchi_k
      integer     level, levelmax, myid, mylen
      integer     li,gi,own
      real*8      xmin, ymin, zmin, x, y, z, rr
      integer     proc_return_myid
      external    proc_return_myid
      character*12 fname, procname
      logical     ltrace
      parameter ( ltrace = .false. )
      integer    i,j,k,index, l


      call load_pointers(gridnum)

      if (ltrace) write(*,*) 'comp_amr_error: On grid: ',gridnum

      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      hg   = gr_h(gridnum)
      levelmax = 6
       
      level = grid_return_level(gridnum)
      xmin = q(gr_x(gridnum))
      ymin = q(gr_y(gridnum))
      zmin = q(gr_z(gridnum))

      if (.true.) then
         x = 0.d0
         y = 0.d0
         z = 0.d0
         !
         ! Find first grid on finest level which contains
         ! this point:
         !
         call level_finest_containing(x,y,z,li,gi,own)
         if (gi.eq.gridnum) then
            ! Find where/if r=0 is on grid:
            i = NINT( ( 0.d0 - xmin ) / hg ) + 1
            j = NINT( ( 0.d0 - ymin ) / hg ) + 1
            k = NINT( ( 0.d0 - zmin ) / hg ) + 1
            if (i.ge.1 .and. j .ge. 1 .and. k .ge. 1 .and.
     *          i.le.nx.and. j .le.ny .and. k .le.nz ) then
               rr = sqrt( q(gr_x(gridnum)+i-1)**2
     *                   +q(gr_y(gridnum)+j-1)**2
     *                   +q(gr_z(gridnum)+k-1)**2 )
               if(rr.le.5.e-4) then
                  index = (k-1)*ny*nx + (j-1)*nx + i-1
                  call int2str(proc_return_myid(),procname)
                  fname = 'fort.10.'
                  mylen = len_trim(fname)
                  fname(mylen+1:) = procname
                  open(unit=10,file=fname,STATUS='UNKNOWN')
                  write(10,"(2F11.5,4E15.7)") gr_t(gridnum), hg,
     *                   q(gr_alpha+index),q(gr_curvature+index),
     *                   q(gr_phim+index),q(gr_energy_const+index)
               end if
            end if
         end if
      end if

      if (shadow .eq. 0) then

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
	      x = xmin + hg*(i-1)
              y = ymin + hg*(j-1)
              z = zmin + hg*(k-1)
              index = (k-1)*ny*nx + (j-1)*nx + i

              q(gr_error + index - 1) =  (
!
!              q(gr_error + index - 1) = hg**2 * (
     *              + q(gr_K00+index-1)**2
     *              + q(gr_K01+index-1)**2
     *              + q(gr_K02+index-1)**2
     *              + q(gr_K03+index-1)**2
     *              + q(gr_K11+index-1)**2
     *              + q(gr_K12+index-1)**2
     *              + q(gr_K13+index-1)**2
     *              + q(gr_K22+index-1)**2
     *              + q(gr_K23+index-1)**2
     *              + q(gr_K33+index-1)**2
!
     *              + q(gr_g00+index-1)**2
     *              + q(gr_g01+index-1)**2
     *              + q(gr_g02+index-1)**2
     *              + q(gr_g03+index-1)**2     
     *              + q(gr_g11+index-1)**2
     *              + q(gr_g12+index-1)**2
     *              + q(gr_g13+index-1)**2
     *              + q(gr_g22+index-1)**2
     *              + q(gr_g23+index-1)**2
     *              + q(gr_g33+index-1)**2
!     
!     *              + q(gr_phir+index-1)**2
!     *              + q(gr_phic+index-1)**2
!     *              + q(gr_phim+index-1)**2          
!
!     *              + q(gr_pir+index-1)**2
!     *              + q(gr_pic+index-1)**2
!     *              + q(gr_pim+index-1)**2          
!
     *                        )
     
             !output the central value
	      if ((level .eq. levelmax) .AND. 
     *	           (abs(x)+abs(y)+abs(z).lt. 0.001)) then
                write(*,*) 'phir_0 = ', local_time, q(gr_phir+index-1)     
                write(*,*) 'phic_0 = ', local_time, q(gr_phic+index-1)
                write(*,*) 'alp_0  = ', local_time, q(gr_alpha+index-1)     
                write(*,*) 'g11_0  = ', local_time, q(gr_g11+index-1)              
	      end if     
            end do
          end do
        end do
	
      end if
         ! Keep track of the minimum of alpha....
         ! but since code is setup for maximums,
         ! store 1-min(alpha):
         !  
!         call find_absmin3d( q(gr_alpha),  tmp_maxchi,
          call find_absmin3d( q(gr_g11),  tmp_maxchi,
     *                    tmp_maxchi_i,
     *                    tmp_maxchi_j,
     *                    tmp_maxchi_k,
     *                    nx,ny,nz           )
!         tmp_maxchi = 1.d0 - tmp_maxchi
         ! 
         ! Check if we've found a new maximum:
         !        
         if ( tmp_maxchi .ge. maxchi ) then
           maxchi   = tmp_maxchi
            maxchi_i = tmp_maxchi_i
            maxchi_j = tmp_maxchi_j
            maxchi_k = tmp_maxchi_k
            maxchi_l = gr_level(gridnum)
            maxchi_t = gr_t(gridnum)
         end if


      return
      end    ! END: comp_amr_error
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_ahfind:                                                              cc
cc                                                                            cc
cc         Project specific routine to find where to excise and hence         cc
cc         should be located in comp_amr_error.f.                             cc
cc         This routines requires centers for up to masks to be specified.    cc
cc         All this is subject to change, but for now this routine is         cc
cc         responsible for determining on this grid whether the quantities    cc
cc             local_left(1..2,1..3)                                          cc
cc             local_right(1..2,1..3)                                         cc
cc         need to be adjusted. These are local to this processor, and        cc
cc         global bounds will be found in level_bhmask_local().               cc
cc                                                                            cc
cc         To see the data structures for this routine, look in:              cc
cc                  had/include/mask.inc                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_ahfind(gi)
      implicit none
      integer  gi
      include 'mpif.h'
      include 'mpi_stuff.inc'
      include 'grid.inc'
      include 'mask.inc'
      include 'param.inc'      
      integer  nx,ny,nz, i,j,k, kk
      real(kind=8) my_left(3,max_num_masks),
     *             my_right(3,max_num_masks),
     *             bh_mass, time
      logical     bh_found
      character*3  tmpstr
      character*12 myname

      logical     ltrace
      parameter ( ltrace  = .false.)
      logical     ltrace2
      parameter ( ltrace2 = .false.)
      logical     ltrace3
      parameter ( ltrace3 = .false.)
      logical     ltrace4
      parameter ( ltrace4 = .false.)

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

c      time = par(P_TIME)

      num_masks = 0

      !
      ! Do the real computation for this grid:
      !
      do kk = 1, max_num_masks
        !
        if (ltrace3) then
           write(*,99) myid, 'grid_ahfind: Center:',kk
           do i = 1, 3
              write(*,98) myid, '    center coord',mask_center(i,kk)
              write(*,98) myid, '    left/right: ',mask_left(i,kk),
     *                                            mask_right(i,kk)
           end do
           write(*,98) myid, '    radius:     ',mask_radius(1,kk)
        end if
        !
        bh_found    = .false.
	
	if (kk .EQ. 1) bh_mass = bh1_mass
	if (kk .EQ. 2) bh_mass = bh2_mass	
        if(forcedmerge) bh_mass = bh1_mass + bh2_mass
        if (bh1_mass.le.0.0) bh_mass = abs(bh1_mass)
		
c        if (time .GT. 39.5) bh_mass = bh1_mass + bh2_mass

	if (local_time .LT. 1E-10) then
          call  compute_boundary_bh(forcedmerge,q(gr_mask),
     &           q(gr_g00),   q(gr_g01),   q(gr_g02),
     &           q(gr_g03),
     &           q(gr_g11),   q(gr_g12),   q(gr_g13),
     &           q(gr_g22),   q(gr_g23),   q(gr_g33),
     &           q(gr_K11),   q(gr_K12),   q(gr_K13),
     &           q(gr_K22),   q(gr_K23),   q(gr_K33),
     &           q(gr_d1g00), q(gr_d1g01), q(gr_d1g02),
     &           q(gr_d1g03),
     &           q(gr_d1g11), q(gr_d1g12), q(gr_d1g13),
     &           q(gr_d1g22), q(gr_d1g23), q(gr_d1g33),
     &           q(gr_d2g00), q(gr_d2g01), q(gr_d2g02),
     &           q(gr_d2g03),
     &           q(gr_d2g11), q(gr_d2g12), q(gr_d2g13),
     &           q(gr_d2g22), q(gr_d2g23), q(gr_d2g33),
     &           q(gr_d3g00), q(gr_d3g01), q(gr_d3g02),
     &           q(gr_d3g03),
     &           q(gr_d3g11), q(gr_d3g12), q(gr_d3g13),
     &           q(gr_d3g22), q(gr_d3g23), q(gr_d3g33), 
     &           q(gr_rad_speed), q(gr_rad_exp),       
     &           q(gr_x(gi)), q(gr_y(gi)), q(gr_z(gi)),
     &           mask_center(:,kk), mask_radius(1,kk),
     &           mask_left, mask_right,
     &           my_left(:,kk), my_right(:,kk),
     &           bh_found,    nx, ny, nz, mask_coords, kk,
     &           bh_mass, 0)
        else
          call  compute_boundary_bh(forcedmerge,q(gr_mask),
     &           q(gr_g00_np1),   q(gr_g01_np1),   q(gr_g02_np1),
     &           q(gr_g03_np1),
     &           q(gr_g11_np1),   q(gr_g12_np1),   q(gr_g13_np1),
     &           q(gr_g22_np1),   q(gr_g23_np1),   q(gr_g33_np1),
     &           q(gr_K11_np1),   q(gr_K12_np1),   q(gr_K13_np1),
     &           q(gr_K22_np1),   q(gr_K23_np1),   q(gr_K33_np1),
     &           q(gr_d1g00_np1), q(gr_d1g01_np1), q(gr_d1g02_np1),
     &           q(gr_d1g03_np1),
     &           q(gr_d1g11_np1), q(gr_d1g12_np1), q(gr_d1g13_np1),
     &           q(gr_d1g22_np1), q(gr_d1g23_np1), q(gr_d1g33_np1),
     &           q(gr_d2g00_np1), q(gr_d2g01_np1), q(gr_d2g02_np1),
     &           q(gr_d2g03_np1),
     &           q(gr_d2g11_np1), q(gr_d2g12_np1), q(gr_d2g13_np1),
     &           q(gr_d2g22_np1), q(gr_d2g23_np1), q(gr_d2g33_np1),
     &           q(gr_d3g00_np1), q(gr_d3g01_np1), q(gr_d3g02_np1),
     &           q(gr_d3g03_np1),
     &           q(gr_d3g11_np1), q(gr_d3g12_np1), q(gr_d3g13_np1),
     &           q(gr_d3g22_np1), q(gr_d3g23_np1), q(gr_d3g33_np1), 
     &           q(gr_rad_speed), q(gr_rad_exp),       
     &           q(gr_x(gi)), q(gr_y(gi)), q(gr_z(gi)),
     &           mask_center(:,kk), mask_radius(1,kk),
     &		 mask_left, mask_right,
     &           my_left(:,kk), my_right(:,kk),
     &           bh_found,    nx, ny, nz, mask_coords,kk,
     &           bh_mass, 1)
        end if
        
	if ( bh_found ) then
            bh_true(kk) = .true.
            num_masks = num_masks + 1
            if(ltrace)write(*,99)myid, 'grid_ahfind: BH found: ',kk
            do i = 1, 3
               if (ltrace2)write(*,98)myid, ' l/right:',
     *                              my_left(i,kk),my_right(i,kk)
               if (my_left(i,kk)  .lt. local_left(i,kk)  )
     *             local_left(i,kk)  = my_left(i,kk)
               if (my_right(i,kk) .gt. local_right(i,kk) )
     *             local_right(i,kk) = my_right(i,kk)
               if(ltrace)write(*,98)myid,'      local_left/right(i,k)=',
     *            local_left(i,kk), local_right(i,kk)
            end do
        end if
         !
         !
        if (ltrace4) then
          call int2str(kk,tmpstr)
          myname = 'myrad_exp2D'//tmpstr
          call field_out2d(
     *              q(gr_rad_exp+((nz+1)/2-1)*nx*ny ),
     *              gr_t(gi), myname,
!    *             gr_t(gi), 'myrad_exp2D'//tmpstr,
     *              gr_minx(gi),gr_maxx(gi),
     *              gr_miny(gi),gr_maxy(gi),
     *              nx,ny,myid)
        end if 
         !
      end do

      if (ltrace2) then
         write(*,99) myid, 'grid_ahfind: gi     = ',gi
         write(*,99) myid, 'grid_ahfind: nx/y/z = ',nx,ny,nz
         write(*,99) myid, 'grid_ahfind: Done.'
      end if

      if (local_time .GT. 1E-10) then
        if(forcedmerge) then
          bh1_x0 = 0.5*(mask_coords(1,1)+mask_coords(2,1))
          bh1_y0 = 0.5*(mask_coords(3,1)+mask_coords(4,1))
          bh1_z0 = 0.5*(mask_coords(5,1)+mask_coords(6,1))
          bh2_x0 = bh1_x0
          bh2_y0 = bh1_y0
          bh2_z0 = bh1_z0
        else
          bh1_x0 = 0.5*(mask_coords(1,1)+mask_coords(2,1))
          bh1_y0 = 0.5*(mask_coords(3,1)+mask_coords(4,1))
          bh1_z0 = 0.5*(mask_coords(5,1)+mask_coords(6,1))
          bh2_x0 = 0.5*(mask_coords(1,2)+mask_coords(2,2))
          bh2_y0 = 0.5*(mask_coords(3,2)+mask_coords(4,2))
          bh2_z0 = 0.5*(mask_coords(5,2)+mask_coords(6,2))
        end if
       
        bh1_exc_rad = mask_radius(1,1)
        bh2_exc_rad = mask_radius(1,2)
      end if



  98  format('[',I3,'] ',A,3F18.6)
  99  format('[',I3,'] ',A,3I5)
      return
      end       ! END: grid_ahfind

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  compute_boundary_bh:                                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_boundary_bh( forcedmerge,
     *           mask, g00, g01, g02, g03,
     *           g11, g12, g13, g22, g23, g33,
     *           K11, K12, K13, K22, K23, K33,
     *           d1g00, d1g01, d1g02, d1g03,
     *           d1g11, d1g12, d1g13, d1g22, d1g23, d1g33,
     *           d2g00, d2g01, d2g02, d2g03,
     *           d2g11, d2g12, d2g13, d2g22, d2g23, d2g33,
     *           d3g00, d3g01, d3g02, d3g03,
     *           d3g11, d3g12, d3g13, d3g22, d3g23, d3g33, 
     *           rad_speed, rad_exp,       
     *           xx,   yy,   zz,
     *           bh_center, bh_radius, m_left, m_right, 
     *           bh_min, bh_max, bh_true,
     *           nx, ny, nz, mask_coords, khole,
     *           bhmass, stage)
      implicit   none
      integer        nx,ny,nz
      real(kind=8)   g11(nx,ny,nz), g12(nx,ny,nz), g13(nx,ny,nz),
     *               g22(nx,ny,nz), g23(nx,ny,nz), g33(nx,ny,nz),
     *               g00(nx,ny,nz), g01(nx,ny,nz), g02(nx,ny,nz),
     *               g03(nx,ny,nz), mask(nx,ny,nz), 
     *               d1g11(nx,ny,nz), d1g12(nx,ny,nz), d1g13(nx,ny,nz),
     *               d1g22(nx,ny,nz), d1g23(nx,ny,nz), d1g33(nx,ny,nz),
     *               d1g00(nx,ny,nz), d1g01(nx,ny,nz), d1g02(nx,ny,nz),
     *               d1g03(nx,ny,nz), 
     *               d2g11(nx,ny,nz), d2g12(nx,ny,nz), d2g13(nx,ny,nz),
     *               d2g22(nx,ny,nz), d2g23(nx,ny,nz), d2g33(nx,ny,nz),
     *               d2g00(nx,ny,nz), d2g01(nx,ny,nz), d2g02(nx,ny,nz),
     *               d2g03(nx,ny,nz), 
     *               d3g11(nx,ny,nz), d3g12(nx,ny,nz), d3g13(nx,ny,nz),
     *               d3g22(nx,ny,nz), d3g23(nx,ny,nz), d3g33(nx,ny,nz),
     *               d3g00(nx,ny,nz), d3g01(nx,ny,nz), d3g02(nx,ny,nz),
     *               d3g03(nx,ny,nz), 
     *               K11(nx,ny,nz), K12(nx,ny,nz), K13(nx,ny,nz),
     *               K22(nx,ny,nz), K23(nx,ny,nz), K33(nx,ny,nz),
     *               rad_speed(nx,ny,nz), rad_exp(nx,ny,nz),
     *               xx(nx),        yy(ny),        zz(nz)
      real(kind=8)   bh_center(3), bh_radius, mask_coords(6,2)
      real(kind=8)   m_left(3,2), m_right(3,2)
      real(kind=8)   bh_min(3),bh_max(3),fl,fr,tempval,sa,sb,sc
      integer	     stage
      logical        bh_true, forcedmerge
      ! Local vars:
      integer        i,j,k,ib,jb,kb,ISWEEP,IRUN,i_min,i_max,khole,
     *               j_min,j_max,k_min,k_max, in_il(3),in_ir(3),
     *               in_jl(3),in_jr(3),in_kl(3),in_kr(3), ISINGLE
      real(kind=8)   x,y,z,r,h11,h12,h13,h22,h23,h33,
     *               huu11,huu12,huu13,huu22,huu23,huu33,
     *               deth, b1, b2, b3, alp, alp2,
     *               eK11,eK12,eK13,eK22,eK23,eK33,etrK,     
     *               nr1, nr2, nr3, nru1, nru2, nru3, nrnorm,
     *               D100,D101,D102,D103,D111,D112,D113,D122,D123,D133,
     *               D200,D201,D202,D203,D211,D212,D213,D222,D223,D233,
     *               D300,D301,D302,D303,D311,D312,D313,D322,D323,D333,     
     *               Dddu111,Dddu112,Dddu113,Dddu121,Dddu122,Dddu123,
     *               Dddu131,Dddu132,Dddu133,    
     *               Dddu211,Dddu212,Dddu213,Dddu221,Dddu222,Dddu223,
     *               Dddu231,Dddu232,Dddu233,    
     *               Dddu311,Dddu312,Dddu313,Dddu321,Dddu322,Dddu323,
     *               Dddu331,Dddu332,Dddu333,    
     *               Dudd111,Dudd112,Dudd113,Dudd122,Dudd123,Dudd133,
     *               Dudd211,Dudd212,Dudd213,Dudd222,Dudd223,Dudd233,
     *               Dudd311,Dudd312,Dudd313,Dudd322,Dudd323,Dudd333,     
     *               G111,G112,G113,G122,G123,G133,
     *               G211,G212,G213,G222,G223,G233,
     *               G311,G312,G313,G322,G323,G333,
     *               Gndd11,Gndd12,Gndd13,Gndd22,Gndd23,Gndd33,
     *               Dndd11,Dndd12,Dndd13,Dndd22,Dndd23,Dndd33,
     *               d1nr1,d1nr2,d1nr3,d2nr1,d2nr2,
     *               d2nr3,d3nr1,d3nr2,d3nr3,
     *               bhmass, besafe

      logical     ltrace, interpolation
      parameter ( ltrace = .false.)
       
      real(kind=8)   bhcen(3), h, bhl(3), bhr(3), fboth, fsing 

      bh_true = .false.
      !
      !
      ! Initialize bounds:
      !
      bh_min(1) = xx(nx)
      bh_min(2) = yy(ny)
      bh_min(3) = zz(nz)
      !
      bh_max(1) = xx(1)
      bh_max(2) = yy(1)
      bh_max(3) = zz(1)

      h = xx(2) - xx(1)
      bhl = 0.0
      bhr = 0.0
      rad_speed = 0.0
      
      fboth = 1.
      fsing = 1.
      if(forcedmerge) fboth = -1.0
      if(forcedmerge) fsing = 1.0
      
      if(khole.eq.1) rad_exp = 0.0
      if(khole.eq.2) then
       rad_speed = rad_exp
       rad_exp = 0.0
      end if
    

      !
      ! Store indices as well:
      !
      in_il(1) = nx
      in_ir(1) = 1
      in_jl(2) = ny
      in_jr(2) = 1
      in_kl(3) = nz
      in_kr(3) = 1

c get the sizes to try an elliptical surface
c try an average to get approx sizes
      
      DO ISINGLE = 1, 1

      if(ISINGLE.eq.2) then
      sa = 0.65+.35*abs(mask_coords(2,2)-mask_coords(1,1))
      sb = 0.65+.35*abs(mask_coords(4,2)-mask_coords(3,1))
      sc = 0.65+.35*abs(mask_coords(6,2)-mask_coords(5,1))

c      bhcen(1) = 0.5*(mask_coords(2,2)+mask_coords(1,1))
c      bhcen(2) = 0.5*(mask_coords(4,2)+mask_coords(3,1))
c      bhcen(3) = 0.5*(mask_coords(6,2)+mask_coords(5,1))

      bhcen(1) = 0.5*(m_left(1,2)+m_right(1,1))
      bhcen(2) = 0.5*(m_left(2,2)+m_right(2,1))
      bhcen(3) = 0.5*(m_left(3,2)+m_right(3,1))


      else
      sa = 0.65+.35*abs(mask_coords(2,khole)-mask_coords(1,khole))
      sb = 0.65+.35*abs(mask_coords(4,khole)-mask_coords(3,khole))
      sc = 0.65+.35*abs(mask_coords(6,khole)-mask_coords(5,khole))

c      bhcen(1) = 0.5*(mask_coords(2,khole)+mask_coords(1,khole))
c      bhcen(2) = 0.5*(mask_coords(4,khole)+mask_coords(3,khole))
c      bhcen(3) = 0.5*(mask_coords(6,khole)+mask_coords(5,khole))

      bhcen(1) = 0.5*(m_left(1,khole)+m_right(1,khole))
      bhcen(2) = 0.5*(m_left(2,khole)+m_right(2,khole))
      bhcen(3) = 0.5*(m_left(3,khole)+m_right(3,khole))

      end if
      
c      if(stage.lt.1) bhcen = bh_center
      ! the default simplest values     
      bhcen = bh_center   ! the center is computed in the standard way
      interpolation = .true.   ! no interpolation to find the horizon
c      sa = 1.0; sb = 1.0; sc = 1.0   ! spherical surfaces for the AH
      
      do k = 1, nz
         z = zz(k)
         !
         do j = 1, ny
            y = yy(j)
            !
            do i = 1, nx
               x = xx(i)
                 besafe = 1.
               !
               !       compute the 3+1 quantities
               !
               h11 = g11(i,j,k)
               h12 = g12(i,j,k)
               h13 = g13(i,j,k)
               h22 = g22(i,j,k)
               h23 = g23(i,j,k)
               h33 = g33(i,j,k)
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
                 besafe = -1.
               end if

               huu11 = huu11/deth
               huu12 = huu12/deth
               huu13 = huu13/deth
               huu22 = huu22/deth
               huu23 = huu23/deth
               huu33 = huu33/deth
               !
               b1   =   g01(i,j,k)*huu11 + g02(i,j,k) *huu12
     *                + g03(i,j,k)*huu13
               b2   =   g01(i,j,k)*huu12 + g02(i,j,k) *huu22
     *                + g03(i,j,k)*huu23
               b3   =   g01(i,j,k)*huu13 + g02(i,j,k) *huu23
     *                + g03(i,j,k)*huu33
               alp2 = - g00(i,j,k)       + b1**2 *h11
     *                + b2**2*h22 + 2*(b1*(b2*h12 + b3*h13)
     *                + b2*b3*h23) + b3**2*h33
	       if(alp2.lt.1e-4) then
	        alp2 = 1.e-4
	        besafe = -1.
	       end if
 
               alp  = sqrt(alp2)
               !
               ! compute the speed in the radial direction
               !

               nr1    = x - bhcen(1)  
               nr2    = y - bhcen(2)  
               nr3    = z - bhcen(3) 

               if (ISINGLE.eq.2) then
                  nr1 = x 
	          nr2 = y 
	          nr3 = z 
	       end if
	       
               r      = sqrt(nr1**2 + nr2**2 + nr3**2)	     
               
	       ! careful cause the norm has x/a**2 , etc 
               nr1 = nr1 /sa**2
               nr2 = nr2 /sb**2
               nr3 = nr3 /sc**2
 
	       if ((r .lt. 2.5*bhmass*ISINGLE**2) .AND. (r.gt.1e-3)) then  
	       
	       nrnorm = sqrt(   huu11*nr1*nr1 + huu22*nr2*nr2
     *                         + huu33*nr3*nr3
     *                    +2.0d0*(huu12*nr1*nr2 + huu13*nr1*nr3
     *                            + huu23*nr2*nr3))
	       if(nrnorm.lt.1.e-8) then
	         nrnorm = 1.e-8
                 besafe = -1.
	       end if
 
                 nr1 = nr1/nrnorm
                 nr2 = nr2/nrnorm
                 nr3 = nr3/nrnorm      
	       
               !
               rad_speed(i,j,k) = 0.0d0 ! - (nr1*b1 + nr2*b2 + nr3*b3) + alp
               !
               !
               ! compute the expansion in the radial direction
               !
	       
	       nru1 = huu11*nr1 + huu12*nr2 + huu13*nr3
               nru2 = huu12*nr1 + huu22*nr2 + huu23*nr3
               nru3 = huu13*nr1 + huu23*nr2 + huu33*nr3

	       D111 = d1g11(i,j,k)
               D112 = d1g12(i,j,k)
               D113 = d1g13(i,j,k)
	       D122 = d1g22(i,j,k)
               D123 = d1g23(i,j,k)
               D133 = d1g33(i,j,k)
	       D100 = d1g00(i,j,k)
               D101 = d1g01(i,j,k)
               D102 = d1g02(i,j,k)
	       D103 = d1g03(i,j,k)
	       D211 = d2g11(i,j,k)
               D212 = d2g12(i,j,k)
               D213 = d2g13(i,j,k)
	       D222 = d2g22(i,j,k)
               D223 = d2g23(i,j,k)
               D233 = d2g33(i,j,k)
	       D200 = d2g00(i,j,k)
               D201 = d2g01(i,j,k)
               D202 = d2g02(i,j,k)
	       D203 = d2g03(i,j,k)
	       D311 = d3g11(i,j,k) 
               D312 = d3g12(i,j,k)
               D313 = d3g13(i,j,k)
	       D322 = d3g22(i,j,k)
               D323 = d3g23(i,j,k)
               D333 = d3g33(i,j,k)
	       D300 = d3g00(i,j,k)
               D301 = d3g01(i,j,k) 
               D302 = d3g02(i,j,k)
	       D303 = d3g03(i,j,k)

	       	       
               eK11 = -(-alp*K11(i,j,k) - D101 - D101 
     *              + b1*(D111 + D111) + b2*(D112 + D112) 
     *              + b3*(D113 + D113))/(2.0*alp)
               eK12 = -(-alp*K12(i,j,k) - D102 - D201 
     *              + b1*(D112 + D211) + b2*(D122 + D212) 
     *              + b3*(D123 + D213))/(2.0*alp)
               eK13 = -(-alp*K13(i,j,k) - D103 - D301 
     *              + b1*(D113 + D311) + b2*(D123 + D312) 
     *              + b3*(D133 + D313))/(2.0*alp)
               eK22 = -(-alp*K22(i,j,k) - D202 - D202 
     *              + b1*(D212 + D212) + b2*(D222 + D222) 
     *              + b3*(D223 + D223))/(2.0*alp)
               eK23 = -(-alp*K23(i,j,k) - D203 - D302 
     *              + b1*(D213 + D312) + b2*(D223 + D322)
     *              + b3*(D233 + D323))/(2.0*alp)
               eK33 = -(-alp*K33(i,j,k) - D303 - D303 
     *              + b1*(D313 + D313) + b2*(D323 + D323) 
     *              + b3*(D333 + D333))/(2.0*alp)
	       
               etrK = huu11*eK11 + huu12*eK12 + huu13*eK13
     *              + huu12*eK12 + huu22*eK22 + huu23*eK23
     *              + huu13*eK13 + huu23*eK23 + huu33*eK33
	       
               Dddu111 = huu11*D111 + huu12*D112 + huu13*D113
               Dddu112 = huu12*D111 + huu22*D112 + huu23*D113
               Dddu113 = huu13*D111 + huu23*D112 + huu33*D113
               Dddu121 = huu11*D112 + huu12*D122 + huu13*D123
               Dddu122 = huu12*D112 + huu22*D122 + huu23*D123
               Dddu123 = huu13*D112 + huu23*D122 + huu33*D123
               Dddu131 = huu11*D113 + huu12*D123 + huu13*D133
               Dddu132 = huu12*D113 + huu22*D123 + huu23*D133
               Dddu133 = huu13*D113 + huu23*D123 + huu33*D133
               Dddu211 = huu11*D211 + huu12*D212 + huu13*D213
               Dddu212 = huu12*D211 + huu22*D212 + huu23*D213
               Dddu213 = huu13*D211 + huu23*D212 + huu33*D213
               Dddu221 = huu11*D212 + huu12*D222 + huu13*D223
               Dddu222 = huu12*D212 + huu22*D222 + huu23*D223
               Dddu223 = huu13*D212 + huu23*D222 + huu33*D223
               Dddu231 = huu11*D213 + huu12*D223 + huu13*D233
               Dddu232 = huu12*D213 + huu22*D223 + huu23*D233
               Dddu233 = huu13*D213 + huu23*D223 + huu33*D233
               Dddu311 = huu11*D311 + huu12*D312 + huu13*D313
               Dddu312 = huu12*D311 + huu22*D312 + huu23*D313
               Dddu313 = huu13*D311 + huu23*D312 + huu33*D313
               Dddu321 = huu11*D312 + huu12*D322 + huu13*D323
               Dddu322 = huu12*D312 + huu22*D322 + huu23*D323
               Dddu323 = huu13*D312 + huu23*D322 + huu33*D323
               Dddu331 = huu11*D313 + huu12*D323 + huu13*D333
               Dddu332 = huu12*D313 + huu22*D323 + huu23*D333
               Dddu333 = huu13*D313 + huu23*D323 + huu33*D333

               Dudd111 = huu11*D111 + huu12*D211 + huu13*D311
               Dudd112 = huu11*D112 + huu12*D212 + huu13*D312
               Dudd113 = huu11*D113 + huu12*D213 + huu13*D313
               Dudd122 = huu11*D122 + huu12*D222 + huu13*D322
               Dudd123 = huu11*D123 + huu12*D223 + huu13*D323
               Dudd133 = huu11*D133 + huu12*D233 + huu13*D333
               Dudd211 = huu12*D111 + huu22*D211 + huu23*D311
               Dudd212 = huu12*D112 + huu22*D212 + huu23*D312
               Dudd213 = huu12*D113 + huu22*D213 + huu23*D313
               Dudd222 = huu12*D122 + huu22*D222 + huu23*D322
               Dudd223 = huu12*D123 + huu22*D223 + huu23*D323
               Dudd233 = huu12*D133 + huu22*D233 + huu23*D333
               Dudd311 = huu13*D111 + huu23*D211 + huu33*D311
               Dudd312 = huu13*D112 + huu23*D212 + huu33*D312
               Dudd313 = huu13*D113 + huu23*D213 + huu33*D313
               Dudd322 = huu13*D122 + huu23*D222 + huu33*D322
               Dudd323 = huu13*D123 + huu23*D223 + huu33*D323
               Dudd333 = huu13*D133 + huu23*D233 + huu33*D333
	       	       
              !with the first index up                                             
               G111 = Dddu111 - Dudd111/2.
               G112 = (Dddu121 + Dddu211 - Dudd112)/2.
               G113 = (Dddu131 + Dddu311 - Dudd113)/2.
               G122 = Dddu221 - Dudd122/2.
               G123 = (Dddu231 + Dddu321 - Dudd123)/2.
               G133 = Dddu331 - Dudd133/2.
               G211 = Dddu112 - Dudd211/2.
               G212 = (Dddu122 + Dddu212 - Dudd212)/2.
               G213 = (Dddu132 + Dddu312 - Dudd213)/2.
               G222 = Dddu222 - Dudd222/2.
               G223 = (Dddu232 + Dddu322 - Dudd223)/2.
               G233 = Dddu332 - Dudd233/2.
               G311 = Dddu113 - Dudd311/2.
               G312 = (Dddu123 + Dddu213 - Dudd312)/2.
               G313 = (Dddu133 + Dddu313 - Dudd313)/2.
               G322 = Dddu223 - Dudd322/2.
               G323 = (Dddu233 + Dddu323 - Dudd323)/2.
               G333 = Dddu333 - Dudd333/2.

               Gndd11 = nr1*G111 + nr2*G211 + nr3*G311
               Gndd12 = nr1*G112 + nr2*G212 + nr3*G312
               Gndd13 = nr1*G113 + nr2*G213 + nr3*G313
               Gndd22 = nr1*G122 + nr2*G222 + nr3*G322
               Gndd23 = nr1*G123 + nr2*G223 + nr3*G323
               Gndd33 = nr1*G133 + nr2*G233 + nr3*G333

               Dndd11 = nru1*D111 + nru2*D211 + nru3*D311
               Dndd12 = nru1*D112 + nru2*D212 + nru3*D312
               Dndd13 = nru1*D113 + nru2*D213 + nru3*D313
               Dndd22 = nru1*D122 + nru2*D222 + nru3*D322
               Dndd23 = nru1*D123 + nru2*D223 + nru3*D323
               Dndd33 = nru1*D133 + nru2*D233 + nru3*D333
	       
c	       IF(ISINGLE.eq.2) then
               d1nr1=1./sa**2
	       d2nr2=1./sb**2
	       d3nr3 =1./sc**2
               d1nr2 = 0.0d0
	       d1nr3 = 0.0d0
               d2nr1 = 0.0d0
	       d2nr3 = 0.0d0
               d3nr1 = 0.0d0
	       d3nr2 = 0.0d0
	       
c	       ELSE
c               d1nr1 = 1.0d0; d2nr2 = 1.0d0; d3nr3 = 1.0d0
c               d1nr2 = 0.0d0; d1nr3 = 0.0d0
c               d2nr1 = 0.0d0; d2nr3 = 0.0d0
c               d3nr1 = 0.0d0; d3nr2 = 0.0d0
c	       END IF
              
               tempval = - etrK +  
     *         ( (huu11 - nru1*nru1)*d1nr1 + (huu12 - nru1*nru2)*d1nr2
     *         + (huu13 - nru1*nru3)*d1nr3 + (huu12 - nru2*nru1)*d2nr1
     *         + (huu22 - nru2*nru2)*d2nr2 + (huu23 - nru2*nru3)*d2nr3
     *         + (huu13 - nru3*nru1)*d3nr1 + (huu23 - nru3*nru2)*d3nr2
     *         + (huu33 - nru3*nru3)*d3nr3 )/nrnorm
     *         - ( huu11*Gndd11 + huu22*Gndd22 + huu33*Gndd33
     *         + 2.0*(huu12*Gndd12 + huu13*Gndd13 + huu23*Gndd23) )
     *         + nru1*nru1*(eK11  + 0.5*Dndd11)
     *         + nru2*nru2*(eK22  + 0.5*Dndd22)
     *         + nru3*nru3*(eK33  + 0.5*Dndd33)
     *         + 2.0*nru1*nru2*(eK12  + 0.5*Dndd12)
     *         + 2.0*nru1*nru3*(eK13  + 0.5*Dndd13)
     *         + 2.0*nru2*nru3*(eK23  + 0.5*Dndd23)

           if(besafe.lt.0) tempval = 0.

	      IF(ISINGLE.EQ.1) rad_exp(i,j,k) = tempval
              IF(ISINGLE.EQ.2) rad_speed(i,j,k) = tempval
     
	       
!	      if ((rad_speed(i,j,k) .lt. 0) .AND. 
!     *            (rad_exp(i,j,k) .lt. 0) .AND.
!     *            (r .lt. 2.5*bh_radius))  then
	      if (rad_exp(i,j,k) .le. 0..or.alp.le.0.5) then                  
		  bh_true = .true.
c                  if (x .lt. bh_min(1)) bh_min(1) = x
c                  if (y .lt. bh_min(2)) bh_min(2) = y
c                  if (z .lt. bh_min(3)) bh_min(3) = z	    
c                  if (x .gt. bh_max(1)) bh_max(1) = x
c                  if (y .gt. bh_max(2)) bh_max(2) = y
c                  if (z .gt. bh_max(3)) bh_max(3) = z	

          if (x.lt.bh_min(1)) then
                bh_min(1) = x
                in_il(1)=i
                in_il(2)=j
                in_il(3)=k
          end if
          if (y.lt.bh_min(2)) then
                bh_min(2) = y
                in_jl(1)=i
                in_jl(2)=j
                in_jl(3)=k
          end if
          if (z.lt.bh_min(3)) then
                bh_min(3) = z
                in_kl(1)=i
                in_kl(2)=j
                in_kl(3)=k
          end if
          if (x.gt.bh_max(1)) then
                bh_max(1) = x
                in_ir(1)=i
                in_ir(2)=j
                in_ir(3)=k
          end if
          if (y.gt.bh_max(2)) then
                bh_max(2) = y
                in_jr(1)=i
                in_jr(2)=j
                in_jr(3)=k
          end if
          if (z.gt.bh_max(3)) then
                bh_max(3) = z
                in_kr(1)=i
                in_kr(2)=j
                in_kr(3)=k
          end if
	  
!this is the new stuff... push holes if they missbehave
             if(alp.lt.0.225.and.mask(i,j,k).gt.0.and.
     &          mask_coords(1,khole).ge.x) then
                bhl(1)=-2.*h*fsing
                bhr(1)=-0.*h*fboth
             end if
             if(alp.lt.0.225.and.mask(i,j,k).gt.0.and.
     &          mask_coords(3,khole).ge.y) then
                bhl(2)= - 2.*h*fsing
                bhr(2)= - 0.*h*fboth
             end if
             if(alp.lt.0.225.and.mask(i,j,k).gt.0.and.
     &          mask_coords(5,khole).ge.z) then
                bhl(3)=- 2.*h*fsing
                bhr(3)=- 0.*h*fboth
              end if
              if(alp.lt.0.225.and.mask(i,j,k).gt.0.and.
     &          mask_coords(2,khole).le.x) then
                bhr(1) =  2.*h*fsing
                bhl(1) =  0.*h*fboth
              end if
              if(alp.lt.0.225.and.mask(i,j,k).gt.0.and.
     &          mask_coords(4,khole).le.y) then
                bhr(2) =  2.*h*fsing
                bhl(2) =  0.*h*fboth
              end if
              if(alp.lt.0.225.and.mask(i,j,k).gt.0.and.
     &          mask_coords(6,khole).le.z) then
                bhr(3) = 2.*h*fsing
                bhl(3) = 0.*h*fboth
              end if	  
	  
	  
              end if  
	       
	      end if !  if (r .lt. 2.5*bhmass) then  
	      
            end do  
         end do  
      end do  

      !now get the right location from a linear interpolation
      if (interpolation .AND. (ISINGLE.EQ.1)) then
        if(in_il(1).lt.nx.and.in_il(1).gt.1) then
        i = in_il(1)-1
        j = in_il(2)
        k = in_il(3)
        fl = rad_exp(i,j,k)
        fr = rad_exp(i+1,j,k)
        if(fl*fr.lt.-1.e-7.and.abs(fl-fr).gt.1.e-6) then
       bh_min(1)=(fl*xx(i+1)-fr*xx(i))/(fl-fr)
	 end if
        end if

        if(in_ir(1).gt.1.and.in_ir(1).lt.nx) then
        i = in_ir(1)
        j = in_ir(2)
        k = in_ir(3)
        fl = rad_exp(i,j,k)
        fr = rad_exp(i+1,j,k)
        if(fl*fr.lt.-1.e-7.and.abs(fl-fr).gt.1.e-6)
     * bh_max(1)=(fl*xx(i+1)-fr*xx(i))/(fl-fr)
        end if

        if(in_jl(2).lt.ny.and.in_jl(2).gt.1) then
        i = in_jl(1)
        j = in_jl(2)-1
        k = in_jl(3)
        fl = rad_exp(i,j,k)
        fr = rad_exp(i,j+1,k)
        if(fl*fr.lt.-1.e-7.and.abs(fl-fr).gt.1.e-6)
     * bh_min(2)=(fl*yy(j+1)-fr*yy(j))/(fl-fr)
        end if

        if(in_jr(2).gt.1.and.in_jr(2).lt.ny) then
        i = in_jr(1)
        j = in_jr(2)
        k = in_jr(3)
        fl = rad_exp(i,j,k)
        fr = rad_exp(i,j+1,k)
        if(fl*fr.lt.-1.e-7.and.abs(fl-fr).gt.1.e-6)
     * bh_max(2)=(fl*yy(j+1)-fr*yy(j))/(fl-fr)
        end if

        if(in_kl(3).lt.nz.and.in_kl(3).gt.1) then
        i = in_kl(1)
        j = in_kl(2)
        k = in_kl(3)-1
        fl = rad_exp(i,j,k)
        fr = rad_exp(i,j,k+1)
        if(fl*fr.lt.-1.e-7.and.abs(fl-fr).gt.1.e-6)
     * bh_min(3)=(fl*zz(k+1)-fr*zz(k))/(fl-fr)
        end if

        if(in_kr(3).gt.1.and.in_kr(3).lt.nz) then
        i = in_kr(1)
        j = in_kr(2)
        k = in_kr(3)
        fl = rad_exp(i,j,k)
        fr = rad_exp(i,j,k+1)
        if(fl*fr.lt.-1.e-7.and.abs(fl-fr).gt.1.e-6)
     * bh_max(3)=(fl*zz(k+1)-fr*zz(k))/(fl-fr)
        end if


!-----do the push man-------------
!        bh_min = bh_min + bhl
!        bh_max = bh_max + bhr
!--------------------------------

        if (ltrace) then
         write(*,*) 'compute_boundary_bh: bh_true = ',bh_true
         if (bh_true) then
            do i = 1, 3
               write(*,*)'compute_boundary_bh: bh_center: ',bh_center(i)
            end do
            do i = 1, 3
               write(*,*)'compute_boundary_bh: min/max: ',
     *                                     bh_min(i),bh_max(i)
            end do
         end if
      end if

      END IF
    
      END DO


      if(khole.eq.2) rad_exp = rad_exp + rad_speed

      return
      end subroutine compute_boundary_bh


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_asurface_comp:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It needs to compute the various             cc
cc                   surface functions at that point.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_asurface_comp(gi, x,y,z, itheta, jphi)
      implicit     none
      integer      gi, itheta, jphi
      real(kind=8) x,y,z
      include     'grid.inc'
      include     'surfaces.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, rpsi4, ipsi4, fracx, fracy, fracz, R2scalar
      real(kind=8) interp, massADM, lin_momx, lin_momy, ang_momz
      real(kind=8) gthethe, gthephi, gphiphi, gur, tau
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ilow .lt. 1) then
         write(*,*) 'grid_surface_comp: x too small ',ilow
         call my_exit('grid_surface_comp: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_surface_comp: fraction cannot be negative'
         write(*,*) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_surface_comp: fraction cannot be negative')
      end if

      rpsi4 = interp(  q(gr_psi4R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ipsi4  = interp(  q(gr_psi4I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      massADM  = interp(  q(gr_massADM),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      lin_momx  = interp(  q(gr_lin_momx),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)     
      lin_momy  = interp(  q(gr_lin_momy),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ang_momz  = interp(  q(gr_ang_momz),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      R2scalar  = interp(  q(gr_R2scalar),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gthethe   = interp(  q(gr_gthethe),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gthephi   = interp(  q(gr_gthephi),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gphiphi   = interp(  q(gr_gphiphi),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gur       = interp(  q(gr_gur),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)     
      tau       = interp(  q(gr_energy_dens),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
          

      index = (jphi-1)*asf_shape(1) + (itheta-1)
      !index = (jphi-1)*asf_ntheta + (itheta-1)
      q(asf_rap4+index) = rpsi4
      q(asf_iap4+index) = ipsi4
      q(asf_aADM+index) = massADM 
      q(asf_aPx+index)  = lin_momx
      q(asf_aPy+index)  = lin_momy
      q(asf_aJz+index)  = ang_momz      
      q(asf_aR2+index)  = R2scalar            
      q(asf_agtt+index) = gthethe      
      q(asf_agtp+index) = gthephi
      q(asf_agpp+index) = gphiphi                  
      q(asf_agur+index) = gur 
      q(asf_atau+index) = tau                  
      
      if (ltrace) then
         write(*,*) 'grid_surface_comp: ********'
         write(*,*) 'grid_surface_comp: gi,hg     =',gi,hg
         write(*,9) 'grid_surface_comp: x,y,z     =',x,y,z
         write(*,8) 'grid_surface_comp: nx/y/z    =',nx,ny,nz
         write(*,8) 'grid_surface_comp: i/j/klow  =',ilow,jlow,klow
         write(*,9) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         write(*,9) 'grid_surface_comp: xmin/ymin/zmin =',
     *                gr_minx(gi),gr_miny(gi),gr_minz(gi)
         write(*,9) 'grid_surface_comp: rps4       =',rpsi4

       write(*,9) 'testA: ',q(gr_psi4r+(17-1)*nx*ny+(17-1)*nx+(17-1) ),
     *      interp(q(gr_psi4r), 1.d0,1.d0,1.d0,
     *        16,16,16, nx,ny,nz)
       write(*,9) 'testb: ',q(gr_psi4r+(16-1)*nx*ny+(16-1)*nx+(16-1) ),
     *      interp(q(gr_psi4r), 0.d0,0.d0,0.d0,
     *        16,16,16, nx,ny,nz)
      end if

 8    format(A,5I10)
 9    format(A,5F10.4)

      return
      end    ! END: grid_surface_comp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_bsurface_comp:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It needs to compute the various             cc
cc                   surface functions at that point.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_bsurface_comp(gi, x,y,z, itheta, jphi)
      implicit     none
      integer      gi, itheta, jphi
      real(kind=8) x,y,z
      include     'grid.inc'
      include     'surfaces.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, rpsi4, ipsi4, fracx, fracy, fracz, R2scalar
      real(kind=8) interp, massADM, lin_momx, lin_momy, ang_momz
      real(kind=8) gthethe, gthephi, gphiphi, gur, tau
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ilow .lt. 1) then
         write(*,*) 'grid_surface_comp: x too small ',ilow
         call my_exit('grid_surface_comp: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_surface_comp: fraction cannot be negative'
         write(*,*) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_surface_comp: fraction cannot be negative')
      end if

      rpsi4 = interp(  q(gr_psi4R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ipsi4  = interp(  q(gr_psi4I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      massADM  = interp(  q(gr_massADM),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      lin_momx  = interp(  q(gr_lin_momx),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)     
      lin_momy  = interp(  q(gr_lin_momy),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ang_momz  = interp(  q(gr_ang_momz),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)     
      R2scalar  = interp(  q(gr_R2scalar),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gthethe   = interp(  q(gr_gthethe),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gthephi   = interp(  q(gr_gthephi),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gphiphi   = interp(  q(gr_gphiphi),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gur       = interp(  q(gr_gur),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)          
      tau       = interp(  q(gr_energy_dens),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
     
      index = (jphi-1)*bsf_shape(1) + (itheta-1)
      !index = (jphi-1)*asf_ntheta + (itheta-1)
      q(bsf_rbp4+index) = rpsi4
      q(bsf_ibp4+index) = ipsi4
      q(bsf_bADM+index) = massADM
      q(bsf_bPx+index)  = lin_momx
      q(bsf_bPy+index)  = lin_momy
      q(bsf_bJz+index)  = ang_momz      
      q(bsf_bR2+index)  = R2scalar            
      q(bsf_bgtt+index) = gthethe      
      q(bsf_bgtp+index) = gthephi
      q(bsf_bgpp+index) = gphiphi  
      q(bsf_bgur+index) = gur                           
      q(bsf_btau+index) = tau                  
      
      if (ltrace) then
         write(*,*) 'grid_surface_comp: ********'
         write(*,*) 'grid_surface_comp: gi,hg     =',gi,hg
         write(*,9) 'grid_surface_comp: x,y,z     =',x,y,z
         write(*,8) 'grid_surface_comp: nx/y/z    =',nx,ny,nz
         write(*,8) 'grid_surface_comp: i/j/klow  =',ilow,jlow,klow
         write(*,9) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         write(*,9) 'grid_surface_comp: xmin/ymin/zmin =',
     *                gr_minx(gi),gr_miny(gi),gr_minz(gi)
         write(*,9) 'grid_surface_comp: rps4       =',rpsi4

       write(*,9) 'testA: ',q(gr_psi4r+(17-1)*nx*ny+(17-1)*nx+(17-1) ),
     *      interp(q(gr_psi4r), 1.d0,1.d0,1.d0,
     *        16,16,16, nx,ny,nz)
       write(*,9) 'testb: ',q(gr_psi4r+(16-1)*nx*ny+(16-1)*nx+(16-1) ),
     *      interp(q(gr_psi4r), 0.d0,0.d0,0.d0,
     *        16,16,16, nx,ny,nz)
      end if

 8    format(A,5I10)
 9    format(A,5F10.4)

      return
      end    ! END: grid_surface_comp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_csurface_comp:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It needs to compute the various             cc
cc                   surface functions at that point.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_csurface_comp(gi, x,y,z, itheta, jphi)
      implicit     none
      integer      gi, itheta, jphi
      real(kind=8) x,y,z
      include     'grid.inc'
      include     'surfaces.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, rpsi4, ipsi4, fracx, fracy, fracz, R2scalar
      real(kind=8) interp, massADM, lin_momx, lin_momy, ang_momz
      real(kind=8) gthethe, gthephi, gphiphi, gur, tau
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ilow .lt. 1) then
         write(*,*) 'grid_surface_comp: x too small ',ilow
         call my_exit('grid_surface_comp: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_surface_comp: fraction cannot be negative'
         write(*,*) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_surface_comp: fraction cannot be negative')
      end if

      rpsi4 = interp(  q(gr_psi4R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ipsi4  = interp(  q(gr_psi4I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      massADM  = interp(  q(gr_massADM),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      lin_momx  = interp(  q(gr_lin_momx),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)     
      lin_momy  = interp(  q(gr_lin_momy),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ang_momz  = interp(  q(gr_ang_momz),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      R2scalar  = interp(  q(gr_R2scalar),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gthethe   = interp(  q(gr_gthethe),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gthephi   = interp(  q(gr_gthephi),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gphiphi   = interp(  q(gr_gphiphi),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      gur       = interp(  q(gr_gur),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)          
      tau       = interp(  q(gr_energy_dens),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)

     
     
      index = (jphi-1)*csf_shape(1) + (itheta-1)
      !index = (jphi-1)*asf_ntheta + (itheta-1)
      q(csf_rcp4+index) = rpsi4
      q(csf_icp4+index) = ipsi4
      q(csf_cADM+index) = massADM
      q(csf_cPx+index)  = lin_momx
      q(csf_cPy+index)  = lin_momy
      q(csf_cJz+index)  = ang_momz      
      q(csf_cR2+index)  = R2scalar                  
      q(csf_cgtt+index) = gthethe      
      q(csf_cgtp+index) = gthephi
      q(csf_cgpp+index) = gphiphi                
      q(csf_cgur+index) = gur      
      q(csf_ctau+index) = tau                  

      if (ltrace) then
         write(*,*) 'grid_surface_comp: ********'
         write(*,*) 'grid_surface_comp: gi,hg     =',gi,hg
         write(*,9) 'grid_surface_comp: x,y,z     =',x,y,z
         write(*,8) 'grid_surface_comp: nx/y/z    =',nx,ny,nz
         write(*,8) 'grid_surface_comp: i/j/klow  =',ilow,jlow,klow
         write(*,9) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         write(*,9) 'grid_surface_comp: xmin/ymin/zmin =',
     *                gr_minx(gi),gr_miny(gi),gr_minz(gi)
         write(*,9) 'grid_surface_comp: rps4       =',rpsi4

       write(*,9) 'testA: ',q(gr_psi4r+(17-1)*nx*ny+(17-1)*nx+(17-1) ),
     *      interp(q(gr_psi4r), 1.d0,1.d0,1.d0,
     *        16,16,16, nx,ny,nz)
       write(*,9) 'testb: ',q(gr_psi4r+(16-1)*nx*ny+(16-1)*nx+(16-1) ),
     *      interp(q(gr_psi4r), 0.d0,0.d0,0.d0,
     *        16,16,16, nx,ny,nz)
      end if

 8    format(A,5I10)
 9    format(A,5F10.4)

      return
      end    ! END: grid_surface_comp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_modify_data:                                                         cc
cc                   Here one need not do anything. However, this             cc
cc                   gets called after reading in a checkpoint file           cc
cc                   and also after global analysis quantities are called     cc
cc                   and so one can use those to set data.                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_modify_data( gridnum )
      implicit    none
      integer     gridnum
      include    'mpif.h'
      include    'mpi_stuff.inc'
      include    'grid.inc'
      include    'grid_methods.inc'
      include    'glob.inc'
      include    'param.inc'
      integer     nx, ny, nz, numfields
      real*8      hg,  tmp, tmp_maxchi
      real*8      distance_wanted, distance_actual
      real*8      mask_distance_movedx
      external    mask_distance_movedx
      integer     tmp_maxchi_i, tmp_maxchi_j,tmp_maxchi_k
      integer     level
      real*8      xmin, ymin, zmin, x, y, z, xc, yc, zc
      logical     ltrace
      parameter ( ltrace = .true. )
      integer    i,j,k,index, l
      integer    dx_alpha, dy_alpha, dz_alpha, nd, d_type
      real*8     local_dx, local_dy, local_dz
      real*8     xs1, ys1, zs1, xs2, ys2, zs2, dist, t1, t2
      real*8     myP, myrho
      real*8     dxAx,dyAx,dzAx
      real*8     dxAy,dyAy,dzAy
      real*8     dxAz,dyAz,dzAz
      real*8     sdetg_pt, r1sq, r2sq, invsdetg_pt
      real*8     myl2norm3D
      external   myl2norm3D
      !
      ! Factor by which rho must be greater than vacuum to have nonzero Bfield:
      !
      real*8     RHOTHRESH
      parameter (RHOTHRESH = 1d3)
      logical    DONOTHING 
      parameter (DONOTHING = .true.) 
      logical    reset_magnetic_field
      parameter (reset_magnetic_field = .false.) 

      integer     mem_alloc
      external    mem_alloc


      if(ltrace)write(*,99)myid,'grid_modify_data: Done grid: ',gridnum

  98  format('[',I3,'] ',A,3F18.6)
  99  format('[',I3,'] ',A,3I5)
      return
      end    ! END: grid_modify_data

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_test:                                                                cc
cc                   Called at the end of grid_initialize() which             cc
cc                   sets the data on a newly created amr grid.               cc
cc                   This routine is just for testing to see if anything      cc
cc                   is wrong with how it is set (such as a bad metric).      cc
cc                                                                            cc
cc   DWN: This was stub was copied from hyperBSSN_XMHD. It needs to be
cc        modified specifically for XMHD.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_test( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'grid_methods.inc'
      include    'glob.inc'
      include    'param.inc'



      return
      end    ! END: grid_test


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_tracer_update:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It updatesthis point via convecting         cc
cc                   the particle with a "fluid".                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_tracer_update(gi,particlenum,dt)
      implicit     none
      integer      gi,particlenum
      real(kind=8) dt
      include     'grid.inc'
      include     'tracers.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, x, y, z
      real(kind=8) velx, vely, velz
      real(kind=8) fracx, fracy, fracz
      real(kind=8) interp
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Coordinates of particle:
      x = tracers_x(particlenum)
      y = tracers_y(particlenum)
      z = tracers_z(particlenum)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ltrace) then
         write(*,*) 'grid_tracer_update: Enter...dt=',dt
      end if
      if (ilow .lt. 1) then
         write(*,*) 'grid_tracer_update: x too small ',ilow
         call my_exit('grid_tracer_update: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_tracer_update: fraction cannot be negative'
         write(*,*) 'grid_tracer_update: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_tracer_update: fraction cannot be negative')
      end if

      ! Project-dependent: Choose what your particle velocity is:
      velx = interp(  q(gr_shift1),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      vely = interp(  q(gr_shift2),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      velz = interp(  q(gr_shift3),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)

      tracers_x(particlenum) = tracers_x(particlenum) + dt*velx
      tracers_y(particlenum) = tracers_y(particlenum) + dt*vely
      tracers_z(particlenum) = tracers_z(particlenum) + dt*velz

      return
      end    ! END: grid_tracer_update

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  tracers_init_sub                                                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine tracers_init_sub()
       implicit     none
       include    'grid.inc'
       include    'grid_methods.inc'
       include    'glob.inc'
       include     'param.inc'
       include     'tracers.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       include     'chr.inc'
       !
       integer      i, j, k, l, mynx,myny,mynz, index
       integer      coarselevel, gi, numcells, offset,part
       integer      numparticles(0:numprocs-1), mynumparticles
       integer      myarray(0:numprocs-1)
       !integer      numparticles(0:10), mynumparticles
       !integer      myarray(0:10)
       integer      SEED
       !parameter (  SEED = 78237 )
       real*8       myx,myy,myz
       real*8       mydx, mydy, allsumrho, sumrho
       real*8       myrho
       !  The range of densities in which  to put particles:
       real*8       MINRHO,        MAXRHO
       parameter (  MINRHO = 5d-8, MAXRHO = 3d-1)
       ! values that will get overwritten with MPI call:
       real*8       UNDEFINED
       parameter (  UNDEFINED = 1d99)
       logical     ltrace
       parameter ( ltrace  = .true. )
       logical     ltrace2
       parameter ( ltrace2  = .false. )

       return
       end     ! END: tracers_init_sub      
