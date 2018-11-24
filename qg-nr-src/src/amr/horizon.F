cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  horizon_evolve:                                                           cc
cc                 Assumes the horizon is already initialized with some       cc
cc                 reasonable guess outisde any true horizon.                 cc
cc                 Follows steps:                                             cc
cc                    (1) Compute the expansion along null geodesics          cc
cc                        normal to putative horizon surface                  cc
cc                    (2) Deforms the surface as a parabolic equation         cc
cc                        in some pseudo-time so that expansion approaches    cc
cc                        zero along the surface.                             cc
cc                    (3) Repeat until expansion is below some threshold.     cc
cc                                                                            cc
cc                        If hole is found to exist:                          cc
cc                                                                            cc
cc                    (4) Find new center for the horizon.                    cc
cc                    (5) Set excision region within horizon if hole exists.  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine horizon_evolve()
      implicit     none
      include     'glob.inc'
      include     'mem.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      include     'largesmall.inc'
      real(kind=8) time, r, phi, theta, x, y, z, avg_x, avg_y, avg_z
      real(kind=8) dphi, dtheta, old_center(3)
      real(kind=8) r_n, r_np1
      ! For finding excision region:
      real(kind=8) avgsq_x, avgsq_y, avgsq_z
      real(kind=8) xdistsq, ydistsq, zdistsq
      real(kind=8) theta_l2norm, r_l2norm, dlambda,dx, oldtheta_l2norm
      integer      kk, i, j, k, tmpi1, tmpi2
      ! Recentering variables:
      integer      i_p1, i_t1, i_p1t1, nk,nj, ni
      real(kind=8) p,t, p1, t1, r_t1, r_p1, r_p1t1
      real(kind=8) x_p1t1, y_p1t1, z_p1t1, x_p1, y_p1, z_p1,
     *             x_t1, y_t1, z_t1, nx,ny,nz, nx_t1,ny_t1,nz_t1,
     *             nx_p1,ny_p1,nz_p1, nx_p1t1,ny_p1t1,nz_p1t1,
     *             nr,nt,np, nr_p1,nr_t1,nr_p1t1, np1,nt1
      ! Results of this search
      logical      mybh_true(max_num_masks)
      !
      integer      ntheta, nphi, numpoints, index, iters
      integer      MAXITERS            ! Max number of evo ITERS
      !parameter  ( MAXITERS = 1 )
      parameter  ( MAXITERS = 2000 )
      real(kind=8) MINRNORM            ! Norm of R at which to terminate
      parameter  ( MINRNORM = 0.50d0 )
      !parameter  ( MINRNORM = 0.10d0 )
      !real(kind=8) THETATHRESH         ! Norm of Theta for a BH
      !parameter  ( THETATHRESH = 0.08d0 )
      real(kind=8) SHRINK         ! How much smaller is region to excise?
      parameter  ( SHRINK = 0.95d0 )
      ! Skip recentering the surface (mostly for debugging)
      !logical      SKIPRECENTER
      !parameter  ( SKIPRECENTER = .true. )
      ! Threshold for how far off-center results in recentering:
      !    (in units of finest level grid spacing)
      real(kind=8) offcenter
      real(kind=8) RECENTERTHRESH
      parameter  ( RECENTERTHRESH = 0.15d0 )
      !parameter  ( RECENTERTHRESH = 0.07d0 ) ! in absolute units
      !parameter  ( RECENTERTHRESH = 0.38d0 ) ! in units of gridpts
      !parameter  ( RECENTERTHRESH = 0.25d0 ) ! in units of gridpts
      ! How often to check for off-center:
      !integer      RECENTERPERIOD
      !parameter  ( RECENTERPERIOD = 15 )
      real(kind=8) interp2dNU
      external     interp2dNU
      real(kind=8) cpi
      parameter  ( cpi  =  3.14159 26535 89793 23846 26433 83279 d0)
      !parameter  ( cpi  =  3.14159265d0 )
      logical      ltrace
      parameter  ( ltrace  = .true. )
      logical      ltrace2            ! output on all procs
      parameter  ( ltrace2 = .false. )
      logical      vtrace             ! SDF output if hole found
      parameter  ( vtrace  = .true. )
      logical      vtrace2            ! SDF output for each step
      parameter  ( vtrace2 = .true. )

      time = grid_return_time(grid_return_finest())
      if (ltrace2 .or. (myid.eq.master.and.ltrace)) then
         write(*,98) myid, 'horizon_evolve: Enter at t=',time
         if (HORIZON_USEICN) write(*,98) myid,'horizon_evolve:Using ICN'
      end if

      ntheta    = horizon_shape(1)
      nphi      = horizon_shape(2)
      numpoints = ntheta * nphi

      ! Using conventions of Pretorius: gr-qc/0407110
      !       dtheta    = 1.d0*cPi / (ntheta-1)
      !       dphi      = 2.d0*cPi / (nphi  -1)
      dtheta    = 1.d0*cPi / (ntheta-1)
      dphi      = 2.d0*cPi / (nphi  -1)
      dx        = grid_return_h(level_return_start(
     *                      level_return_finest(Levelp, maxlev)))
      if (HORIZON_USEICN) then
         dlambda   = 0.5d0 * dx**2
         dlambda   = 0.3 * min(dphi,dtheta)*min(dphi,dtheta)
      else
         dlambda   = 0.5d0 * dx**2
         dlambda   = 0.3 * min(dphi,dtheta)*min(dphi,dtheta)
      end if
      if(ltrace2.or.(myid.eq.master.and.ltrace))then
         write(*,98) myid, 'horizon_evolve: dx      = ',dx
         write(*,98) myid, 'horizon_evolve: dphi    = ',dphi
         write(*,98) myid, 'horizon_evolve: dtheta  = ',dtheta
         write(*,98) myid, 'horizon_evolve: dlambda = ',dlambda
      end if

      do kk = 1, horizon_nholes
         mybh_true(kk) = .false.
         if(ltrace2.or.(myid.eq.master.and.ltrace))then
            write(*,99) myid, 'horizon_evolve: kk: ',kk
            write(*,98) myid, 'horizon_evolve: Center: ',
     *    horizon_center(1,kk),horizon_center(2,kk),horizon_center(3,kk)
         end if
         oldtheta_l2norm = LARGENUMBER
         iters = 1
      !..............................
      ! (1) Compute expansion
      !..............................
 10      continue
         if (HORIZON_USEICN .and. iters.gt.1) then
            ! Swap _np1 and _n points:
            !if(ltrace2)write(*,98)myid,' Swapping _n and _np1 ptrs'
            tmpi1                 = horizon_R(kk)
            tmpi2                 = horizon_Theta(kk)
            horizon_R(kk)         = horizon_R_np1(kk)
            horizon_theta(kk)     = horizon_theta_np1(kk)
            horizon_R_np1(kk)     = tmpi1
            horizon_theta_np1(kk) = tmpi2
         end if
         !
         call horizon_compute(theta_l2norm,r_l2norm,kk)
         !
         if (HORIZON_USEICN .and. iters.eq.1) then
            ! copy data to _np1 vectors:
            !if(ltrace2)write(*,98)myid,' Copying  _n to _np1 for ICN'
            call  mat_copy1d( q(horizon_R(kk)),
     *                        q(horizon_R_np1(kk)),     numpoints)
            call  mat_copy1d( q(horizon_Theta(kk)),
     *                        q(horizon_Theta_np1(kk)), numpoints)
         end if
         !
         if (vtrace2 .and. myid.eq.master) then
             call horizon_out(time+iters*1.d0,kk)
         end if
         !
         if(ltrace2.or.(myid.eq.master.and.ltrace))
     *         write(*,97)myid,'horizon_evolve: Step, Theta/R Norms:',
     *                                  iters,theta_l2norm,r_l2norm
      !..............................
      ! (2) Evolve
      !..............................
         if (theta_l2norm .lt. horizon_thresh) then
            mybh_true(kk) =  .true.
            !bh_true(kk) =  .true.
            if(ltrace2.or.(myid.eq.master.and.ltrace))
     *           write(*,97)myid,'horizon_evolve: Black hole Found:',kk
         else if ( iters .lt. MAXITERS .and.
     *             r_l2norm .gt. MINRNORM ) then
            iters = iters + 1
            do k = 1, nphi
            do j = 1, ntheta
               !
               index = (k-1)*ntheta+(j-1)
               if (HORIZON_USEICN) then
                  !
                  !  {}^(i+1) R^(n+1) = R^(n) - dLambda (1/2)
                  !                             (Theta^n + {}^(i) Theta^(n+1))
                  !
                  r_n   = q(horizon_R_np1(kk)+index)
                  r_np1 = r_n - dlambda * 0.5d0 *(
     *                                q(horizon_theta(kk)    +index)
     *                              + q(horizon_theta_np1(kk)+index) )
                  !
               else
                  r_n = q(horizon_R(kk)+index)
                  !
                  !  R^(n+1) = R^(n) - dLambda Theta
                  !
                  r_np1 = r_n - dlambda * q(horizon_theta(kk)+index)
                  !
               end if
               !
               q(horizon_R(kk)+index)     = r_np1
               if ( q(horizon_R(kk)+index).le.0 ) then
                  if(ltrace2.or.(myid.eq.master.and.ltrace))
     *           write(*,97)myid,'horizon_: R went negative:',kk
                 goto 20
               end if
               !
            end do
            end do
            if (HORIZON_USEICN) then
               ! For iCN, need more substeps:
               call horizon_compute(theta_l2norm,r_l2norm,kk)
               do k = 1, nphi
               do j = 1, ntheta
                  index = (k-1)*ntheta+(j-1)
                  r_n   = q(horizon_R_np1(kk)+index)
                  r_np1 = r_n - dlambda * 0.5d0 *(
     *                                q(horizon_theta(kk)    +index)
     *                              + q(horizon_theta_np1(kk)+index) )
                  q(horizon_R(kk)+index)     = r_np1
               end do
               end do
               ! For iCN, need more substeps:
               call horizon_compute(theta_l2norm,r_l2norm,kk)
               do k = 1, nphi
               do j = 1, ntheta
                  index = (k-1)*ntheta+(j-1)
                  r_n   = q(horizon_R_np1(kk)+index)
                  r_np1 = r_n - dlambda * 0.5d0 *(
     *                                q(horizon_theta(kk)    +index)
     *                              + q(horizon_theta_np1(kk)+index) )
                  q(horizon_R(kk)+index)     = r_np1
               end do
               end do
            end if
            !
            ! If Theta norm is increasing, decrease time step
            !
            if (theta_l2norm .gt. 1.5d0*oldtheta_l2norm) then
               dlambda = 0.5d0 * dlambda
               if(ltrace2.or.(myid.eq.master.and.ltrace))
     *           write(*,98)myid,'horizon_: Theta norm not converging:',
     *           dlambda
                 if (dlambda .lt. 0.0005) then
                    if (myid.eq.master) write(*,98) myid,
     *                    'horizon_evolve :Giving up, dlambda too small'
                    goto 20
                 End if
            end if
            oldtheta_l2norm = theta_l2norm
            !
            ! Monitor out-of-centeredness:
            !
            if (horizon_recenterp.gt.0 .and.
     *          mod(iters,horizon_recenterp).eq.0) then
               avg_x = 0.d0
               avg_y = 0.d0
               avg_z = 0.d0
               !
               do k = 1, nphi
               do j = 1, ntheta
                  !
                  index = (k-1)*ntheta+(j-1)
                  r     = q(horizon_R(kk)+index)
                  phi   = (k-1)*dphi
                  theta = (j-1)*dtheta
                  ! Using conventions of Pretorius: gr-qc/0407110
                  x     = horizon_center(1,kk) + r*cos(phi)*sin(theta)
                  y     = horizon_center(2,kk) + r*sin(phi)*sin(theta)
                  z     = horizon_center(3,kk) + r*         cos(theta)
                  avg_x = avg_x + x
                  avg_y = avg_y + y
                  avg_z = avg_z + z
                  !
               end do
               end do
               !
               avg_x = avg_x / numpoints
               avg_y = avg_y / numpoints
               avg_z = avg_z / numpoints
               offcenter = max(abs(avg_x - horizon_center(1,kk)),
     *                         abs(avg_y - horizon_center(2,kk)),
     *                         abs(avg_z - horizon_center(3,kk))  )
               if(ltrace2.or.(myid.eq.master.and.ltrace))
     *            write(*,98)myid,'horizon_: avg_x/y/z:',avg_x,avg_y,
     *                          avg_z,offcenter
               ! 
               if (offcenter .gt. RECENTERTHRESH .and.
     *             horizon_recenterp .gt. 0 .and.
     *             mod(iters,horizon_recenterp) .eq. 0 ) then
!    *                   .not. SKIPRECENTER) then
               !if (offcenter .gt. RECENTERTHRESH*dx) then
                  old_center(1)        = horizon_center(1,kk)
                  horizon_center(1,kk) = avg_x
                  old_center(2)        = horizon_center(2,kk)
                  horizon_center(2,kk) = avg_y
                  old_center(3)        = horizon_center(3,kk)
                  horizon_center(3,kk) = avg_z
                  !
                  call horizon_recenter( q(horizon_R_np1(kk)),
     *                                   q(horizon_R(kk)),
     *                                   horizon_center(1,kk),
     *                                   old_center(1),
     *                                   dtheta,dphi, ntheta, nphi )
                  !
                  call  mat_copy1d( q(horizon_R_np1(kk)),
     *                        q(horizon_R(kk)),     numpoints)
                  oldtheta_l2norm = LARGENUMBER
                  if(ltrace2.or.(myid.eq.master.and.ltrace))
     *             write(*,98)myid,'horizon_: Recentered'
               else
                  if(ltrace2.or.(myid.eq.master.and.ltrace))
     *             write(*,98)myid,'horizon_: Recentering not necessary'
               end if
            end if
            !
            ! Re-compute expansion Theta:
            goto 10
         else if ( iters .lt. MAXITERS) then
            if(ltrace2.or.(myid.eq.master.and.ltrace))
     *           write(*,97)myid,'horizon_: Too small, stopping',kk
         else 
            if(ltrace2.or.(myid.eq.master.and.ltrace))
     *           write(*,97)myid,'horizon_: Too many iters,stop',kk
         end if
         if (mybh_true(kk)) then
         !if (bh_true(kk)) then
            if (vtrace .and. myid.eq.master) then
                call horizon_out(time,kk)
            end if
      !..............................
      ! (4a) Find new center
      !..............................
            !
            avg_x = 0.d0
            avg_y = 0.d0
            avg_z = 0.d0
            !
            do k = 1, nphi
            do j = 1, ntheta
               !
               index = (k-1)*ntheta+(j-1)
               r     = q(horizon_R(kk)+index)
               phi   = (k-1)*dphi
               theta = (j-1)*dtheta
               ! Using conventions of Pretorius: gr-qc/0407110
               x     = horizon_center(1,kk) + r*cos(phi)*sin(theta)
               y     = horizon_center(2,kk) + r*sin(phi)*sin(theta)
               z     = horizon_center(3,kk) + r*         cos(theta)
               avg_x = avg_x + x
               avg_y = avg_y + y
               avg_z = avg_z + z
               !
            end do
            end do
            !
            avg_x = avg_x / numpoints
            avg_y = avg_y / numpoints
            avg_z = avg_z / numpoints
            !
            old_center(1) = horizon_center(1,kk)
            old_center(2) = horizon_center(2,kk)
            old_center(3) = horizon_center(3,kk)
            !
            !if (.not. SKIPRECENTER) then
            if ( horizon_recenterp .gt. 0 ) then
            !..............................
            ! (4b) Find new radii w/r/t new center
            !..............................
               horizon_center(1,kk) = avg_x
               horizon_center(2,kk) = avg_y
               horizon_center(3,kk) = avg_z
               !
               call horizon_recenter(    q(horizon_R_np1(kk)),
     *                                   q(horizon_R(kk)),
     *                                   horizon_center(1,kk),
     *                                   old_center(1),
     *                                   dtheta,dphi, ntheta, nphi )
               if (myid.eq.master) call horizon_out(time,kk)
               ! Copy _np1 data to _n:
               call  mat_copy1d( q(horizon_R_np1(kk)),
     *                        q(horizon_R(kk)),     numpoints)
               ! For debugging, compute the expansion before outputing:
               call horizon_compute(theta_l2norm,r_l2norm,kk)
               if (myid.eq.master) call horizon_out(time+0.0001d0,kk)
               !
               if (ltrace2.or.(myid.eq.master.and.ltrace)) then
                  write(*,99)myid,'horizon_evolve: New center',kk
                  do i = 1, 3
                     write(*,98)myid,'  Old--> New: ',old_center(i),
     *                                      horizon_center(i,kk)
                  end do
               end if
            else 
               if(ltrace2)write(*,99)myid,'horizon_evolve: keep center'
            end if
            !
      !..............................
      ! (5) Set excision region:
      !..............................
            avgsq_x = 0.d0
            avgsq_y = 0.d0
            avgsq_z = 0.d0
            do k = 1, nphi-1
            do j = 1, ntheta-1
               !
               i     = (k-1)*ntheta+(j-1)
               p     = (k-1)*dphi
               t     = (j-1)*dtheta
               r     = q(horizon_R(kk)+i)
               ! Using conventions of Pretorius: gr-qc/0407110
               x     = old_center(1) + r*cos(p)*sin(t)
               y     = old_center(2) + r*sin(p)*sin(t)
               z     = old_center(3) + r*       cos(t)
               ! for part (5) compute
               xdistsq = abs(x-horizon_center(1,kk))
               ydistsq = abs(y-horizon_center(2,kk))
               zdistsq = abs(z-horizon_center(3,kk))
               avgsq_x = avgsq_x + xdistsq
               avgsq_y = avgsq_y + ydistsq
               avgsq_z = avgsq_z + zdistsq
            end do
            end do
            avgsq_x = (avgsq_x/numpoints)
            avgsq_y = (avgsq_y/numpoints)
            avgsq_z = (avgsq_z/numpoints)
            if(ltrace2)write(*,98)myid,'horizon_evolve: avgsq:',
     *         avgsq_x,avgsq_y,avgsq_z
            !
            mask_bbox(1,kk) =  horizon_center(1,kk) - SHRINK*avgsq_x
            mask_bbox(2,kk) =  horizon_center(1,kk) + SHRINK*avgsq_x
            mask_bbox(3,kk) =  horizon_center(2,kk) - SHRINK*avgsq_y
            mask_bbox(4,kk) =  horizon_center(2,kk) + SHRINK*avgsq_y
            mask_bbox(5,kk) =  horizon_center(3,kk) - SHRINK*avgsq_z
            mask_bbox(6,kk) =  horizon_center(3,kk) + SHRINK*avgsq_z
            ! Following characterization is obsolete now:
            do i = 1, 3
               mask_center(i,kk) = horizon_center(i,kk)
            end do
            mask_radius(1,kk) = min(avgsq_x,avgsq_y,avgsq_z)
            !
            if (ltrace2.or.(myid.eq.master.and.ltrace)) then
               write(*,99)myid,'horizon_evolve: Mask bbox:',kk
               do i = 1, 3
               write(*,98)myid,'b',mask_bbox(2*i-1,kk),mask_bbox(2*i,kk)
               end do
            end if
            !
            ! Let the rest of the code know that a black hole was found
            !    (NB: sometimes a black hole was found in a previous
            !         search, but not this time, so do not set to false 
            !         if none found here)
            bh_true(kk) = .true.
            !
         end if   ! only if mybh_true
 20      continue
         !
      end do


      if (ltrace2.or.(myid.eq.master.and.ltrace)) then
         write(*,99) myid, 'horizon_evolve: Done'
      end if

 97   format('[',I3,'] ',A,I5,3F13.5)
 98   format('[',I3,'] ',A,4F13.5)
 99   format('[',I3,'] ',A,2I5)
      return
      end         ! END: horizon_evolve


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  horizon_out:                                                              cc
cc               output the horizon and its expansion.                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine horizon_out(t,kk)
      implicit     none
      real(kind=8) t
      integer      kk
      include     'mem.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gft_rc, mylen
      integer      gft_out_bbox
      external     gft_out_bbox
      character(9) tmps, fname
      real(kind=8) h_bbox(4)
      real(kind=8) cpi
      parameter  ( cpi  =  3.14159 26535 89793 23846 26433 83279 d0)
      !parameter  ( cpi  =  3.14159265d0 )
      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) then
         write(*,98) myid, 'horizon_out: Enter',t
         write(*,99) myid, 'horizon_out: kk:',kk
         write(*,99) myid, 'horizon_out: shape(1):',horizon_shape(1)
         write(*,99) myid, 'horizon_out: shape(2):',horizon_shape(2)
      end if

      ! Using conventions of Pretorius: gr-qc/0407110
      !       dtheta    = 1.d0*cPi / (ntheta-1)
      !       dphi      = 2.d0*cPi / (nphi  -1)
      h_bbox(1) = 0.d0
      h_bbox(2) = 1.d0 * cPi
      h_bbox(3) = 0.d0
      h_bbox(4) = 2.d0 * cPi

         call int2str(kk,tmps)
         !
         fname           = 'h_theta'
         mylen           = len_trim(fname)
         fname(mylen+1:) = tmps
         if(ltrace)write(*,98)myid, 'horizon_out: Output theta'
#ifdef AIX
         gft_rc = gft_out_bbox( trim(fname)/ /CHAR(0),
     *                               %VAL(t),
     *                               horizon_shape,
     *                               %VAL(2), h_bbox, 
     *                               q(horizon_theta(kk)) )
#else
         gft_rc = gft_out_bbox( fname, t,   horizon_shape,
     *                               2, h_bbox, q(horizon_theta(kk)) )
#endif  
         !
         fname           = 'h_R'
         mylen           = len_trim(fname)
         fname(mylen+1:) = tmps
         if(ltrace)write(*,98)myid, 'horizon_out: Output R'
#ifdef AIX
         gft_rc = gft_out_bbox( trim(fname)/ /CHAR(0),
     *                               %VAL(t),
     *                               horizon_shape,
     *                               %VAL(2), 
     *                               h_bbox, q(horizon_R(kk)) )
#else
         gft_rc = gft_out_bbox( fname,t,horizon_shape,
     *                               2, h_bbox, q(horizon_R(kk)) )
#endif

      if (ltrace) then
         write(*,99) myid, 'horizon_out: Done'
      end if

 98   format('[',I3,'] ',A,F13.7)
 99   format('[',I3,'] ',A,2I5)
      return
      end         ! END: horizon_out


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  horizon_compute:                                                          cc
cc               Compute expansion in "num" surface and return norms.         cc
cc     (a) initialize expansion on this process to a large number             cc
cc     (b) loop over local grids and compute what we can where finest         cc
cc     (c) do a global mpi_allreduce w/ the minimum so all procs have same datac
cc     (d) compute norms                                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine horizon_compute(thetanorm,rnorm,num)
      implicit     none
      real(kind=8) thetanorm, rnorm
      integer      num
      include     'mem.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      include     'largesmall.inc'
      integer      j,k, numpoints, index, ntheta,nphi
      integer      tmpptr
      integer      lev, gi
      real(kind=8) theta
      integer      mem_alloc
      external     mem_alloc
      !    Amount of disspiation to apply
      real(kind=8) HORIZON_EPSDIS
      !parameter  ( HORIZON_EPSDIS = 0.75d0 )  ! too high
      !parameter  ( HORIZON_EPSDIS = 0.5d0 )
       !!parameter  ( HORIZON_EPSDIS = 0.1d0 )
      parameter  ( HORIZON_EPSDIS = 0.05d0 )
       !!parameter  ( HORIZON_EPSDIS = 0.01d0 )
      !parameter  ( HORIZON_EPSDIS = 0.00d0 )
      logical      ltrace
      parameter (  ltrace  = .false. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )

      if (ltrace) then
         write(*,99) myid, 'horizon_compute: Enter, num=',num
      end if

      ntheta    = horizon_shape(1)
      nphi      = horizon_shape(2)
      numpoints = ntheta * nphi

      !..............
      ! (a) Initialize expansion to large number
      !..............
      if(ltrace)write(*,98)myid,'horizon_compute: Initializing.....'
      call load_scal1D(q(horizon_theta(num)),LARGENUMBER,numpoints)

      !..............
      ! (b) Compute expansion at points on the surface using
      !     the finest data available...only do so for those
      !     points located on locally owned grids
      !..............
      if(ltrace)write(*,98)myid,'horizon_compute: Computing........'
      lev = 0
  5   if (level_return_existence(lev))then
         if(ltrace)write(*,99)myid,'    Computing on level: ',lev
         gi = level_return_start(lev)
 10      if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               if(ltrace)write(*,99)myid,'    Computing on grid: ',gi
               call grid_horizon_comp(gi,num)
            end if
            gi = grid_return_sibling(gi)
            goto 10
         end if
         lev = lev + 1
         goto 5
      end if

      !..............
      ! (c) Gather data *from* all procs *to* all procs
      !..............
      if(ltrace)write(*,98)myid,'horizon_compute: Gathering data...'
      tmpptr = mem_alloc(numpoints)
      if(ltrace2)write(*,98)myid,'         calling allreduce'
      call MPI_AllReduce(q(horizon_theta(num)),q(tmpptr),numpoints,
     *              MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      if(ltrace2)write(*,98)myid,'         copying data'
      call mat_copy1d( q(tmpptr), q(horizon_theta(num)), numpoints )
      ! Apply dissipation on surface:
      call mat_copy1d( q(tmpptr), q(horizon_theta(num)), numpoints )

      !..............
      !  Smooth out residual; especially at poles
      !..............
      !call horizon_out(1000.d0,num)
         k = nphi
         do j = 1, theta
            index = (k-1)*ntheta+(j-1)
            if (.true.) then
            ! identify azimuthal sides
            q(horizon_theta(num)+index) =
     *           q(horizon_theta(num)+(1-1)*ntheta+(j-1))
            else
            ! smooth/average
            q(horizon_theta(num)+index) =
     *       0.5d0*q(horizon_theta(num)+(k-1-1)*ntheta+(j-1))
     *     + 0.5d0*q(horizon_theta(num)+(1+1-1)*ntheta+(j-1))
            q(horizon_theta(num)+(1-1)*ntheta+(j-1)) =
     *       0.5d0*q(horizon_theta(num)+(k-1-1)*ntheta+(j-1))
     *     + 0.5d0*q(horizon_theta(num)+(1+1-1)*ntheta+(j-1))
            end if
         end do
         ! Identify all points along poles as the same point
         do k = 2, nphi
            ! top:
            j = ntheta
            q(horizon_theta(num)+(k-1)*ntheta+(j-1) ) = 
     *      q(horizon_theta(num)+(1-1)*ntheta+(j-1))
            ! top:-1
!           j = ntheta-1
!           q(horizon_theta(num)+(k-1)*ntheta+(j-1) ) = 
!    *      q(horizon_theta(num)+(1-1)*ntheta+(j-1))
            ! bottom:
            j = 1
            q(horizon_theta(num)+(k-1)*ntheta+(j-1) ) = 
     *      q(horizon_theta(num)+(1-1)*ntheta+(j-1))
            ! bottom-1:
!           j = 2
!           q(horizon_theta(num)+(k-1)*ntheta+(j-1) ) = 
!    *      q(horizon_theta(num)+(1-1)*ntheta+(j-1))
         end do
!        ! Smooth via averaging near poles:
!        do k = 1, nphi
!           ! top:
!           q(horizon_theta(num)+(k-1)*ntheta+(ntheta-1-1) ) = 
!    *      0.5d0*q(horizon_theta(num)+(k-1)*ntheta+(ntheta  -1) )
!    *     +0.5d0*q(horizon_theta(num)+(k-1)*ntheta+(ntheta-2-1) )
!           ! bottom:
!           q(horizon_theta(num)+(k-1)*ntheta+(2-1) ) = 
!    *      0.5d0*q(horizon_theta(num)+(k-1)*ntheta+(1-1) )
!    *     +0.5d0*q(horizon_theta(num)+(k-1)*ntheta+(3-1) )
!        end do
      call apply_diss2D(q(horizon_theta(num)),
     *                  q(tmpptr),
     *                  HORIZON_EPSDIS, ntheta, nphi)
      !call horizon_out(1000.5d0,num)
      if(ltrace2)write(*,98)myid,'         freeing memory'
      call mem_dealloc(tmpptr,numpoints)

      !..............
      ! (d) Compute norms
      !..............
      if(ltrace)write(*,98)myid,'horizon_compute: Computing norms...'
      thetanorm = 0.d0
      rnorm     = 0.d0
      do k = 1, nphi
      do j = 1, ntheta
         !
         index     = (k-1)*ntheta+(j-1)
         theta     = q(horizon_theta(num)+index)
         if (theta .ge. LARGE_M_SMALL) then
            write(*,99)myid,'***Expansion not computed everywhere *',j,k
            thetanorm = LARGENUMBER
            rnorm     = LARGENUMBER
            goto 20
         end if
         thetanorm = thetanorm + theta**2
         rnorm     = rnorm + q(horizon_R(num)    +index)**2
         !
      end do
      end do
      thetanorm = sqrt(thetanorm/numpoints)
      rnorm     = sqrt(rnorm    /numpoints)

 20   if (ltrace) then
         write(*,98)myid,'horizon_compute: Done:',thetanorm,rnorm
      end if

 98   format('[',I3,'] ',A,3F13.7)
 99   format('[',I3,'] ',A,2I5)
      return
      end         ! END: horizon_compute

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_horizon_comp:                                                        cc
cc                      Compute the expansion on horizon surface kk           cc
cc                      for all points that ( i) lie on this grid             cc
cc                                     and  (ii) do not exist on a finer levelcc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_horizon_comp(gi,kk)
      implicit     none
      integer      gi, kk
      include     'grid.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'output.inc'
      integer      nx,ny,nz
      integer      ntheta, nphi, numpoints, i,j,k, index
      integer      ilow, jlow, klow
      integer      index_xp,index_yp,index_zp
      integer      index_xm,index_ym,index_zm

      integer      index_xpyp,index_xpzp,index_ypzp
      integer      index_xpym,index_xpzm,index_ypzm
      integer      index_xmyp,index_xmzp,index_ymzp
      integer      index_xmym,index_xmzm,index_ymzm

      real(kind=8) theta, phi, dphi, dtheta, r, x,y,z, hg
      real(kind=8) dx,dy,dz,n1,n2,n3,d1n1,d1n2,d1n3,d2n2,d2n3,d3n3
      real(kind=8) alpha, g00, F, rtilde, Rinterp, expansion
      real(kind=8) mybbox(6)
      ! For interpolations:
      real(kind=8) fracx, fracy, fracz, fracj, frack, distp,distt
      real(kind=8) cpi
      parameter  ( cpi  =  3.14159 26535 89793 23846 26433 83279 d0)
      !parameter  ( cpi  =  3.14159265d0 )
      real(kind=8) interp, interp2d
      external     interp, interp2d
      !integer      gft_out_bbox
      !external     gft_out_bbox
      integer      iteration
      save         iteration

      logical      ltrace
      parameter (  ltrace  = .false. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )
      logical      vtrace              ! output SDFs
      parameter (  vtrace  = .false. )

      if (ltrace) then
         write(*,99) myid, 'grid_horizon_comp: gi/kk=',gi,kk
      end if

      call load_pointers(gi)
      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ntheta    = horizon_shape(1)
      nphi      = horizon_shape(2)
      numpoints = ntheta * nphi
      ! Using conventions of Pretorius: gr-qc/0407110
      dphi      = 2.d0*cPi / (nphi  -1)
      dtheta    = 1.d0*cPi / (ntheta-1)

      ! Compute F(x,y,z) = rtilde - R(theta,phi)
      !    (storing it in gr_tmp)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         ! Index into 3D array
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         ! Get coordinates of each point:
         x     = q(gr_x(gi)+i-1) - horizon_center(1,kk)
         y     = q(gr_y(gi)+j-1) - horizon_center(2,kk)
         z     = q(gr_z(gi)+k-1) - horizon_center(3,kk)
         rtilde= sqrt( x**2 + y**2 + z**2 )
         !
         ! Ensure this is the inverse of
         ! the (x,y,z) equations down below:
         if (rtilde .gt. 0) then
            theta = acos( z / rtilde )
         else
            theta = 0.d0
         end if
         if (x.gt.0) then
            phi = atan( y / x )
         else if (x.lt.0) then
            phi = atan( y / x ) + cPi
         else if (y.lt.0) then
            phi = 1.5d0*cPi
         else
            phi = 0.5d0*cPi
         end if
         if (phi .lt. 0) phi   = phi + 2.*cPi
         !
         ! Now that we have all the coordinates for this
         ! point, need to find what points in the surface
         ! to use to interpolate an R value:
         distt = theta / dtheta
         distp = phi   / dphi
         jlow  = INT( distt ) + 1
         klow  = INT( distp ) + 1
         fracj = distt-jlow+1
         frack = distp-klow+1
         !
         if (jlow .lt.1 .or.klow.lt.1 .or.jlow.gt.ntheta
     *       .or. klow.gt.nphi .or. fracj.gt.1
     *       .or. frack.gt.1 .or. fracj.lt.0
     *       .or. fracj.lt.0
     *       .or. theta .lt.0 .or. phi.lt.0) then
            write(*,*) ' * * * * * * * * * *'
            write(*,*) ' rtilde: ',rtilde
            !write(*,*) ' x,y,z: ',x,y,z
            write(*,*) ' theta/phi: ',theta,phi
            write(*,*) 'dtheta/phi: ',dtheta,dphi
            write(*,*) 'j/klow: ',jlow,klow
            write(*,*) 'ntheta/phi: ',ntheta,nphi
            write(*,*) 'fracj/k: ',fracj, frack
            write(*,*) 'distt/p: ',distt,distp
            write(*,*) ' * * * * * * * * * *'
         end if
         !
         ! Linearly interoplate:
         Rinterp = interp2d( q(horizon_R(kk)),fracj,frack,
     *                       jlow,klow,ntheta,nphi)
         !
         q(gr_tmp+index)   = rtilde - Rinterp
         q(gr_flag +index) = rtilde
         q(gr_error+index) = Rinterp
         !
      end do
      end do
      end do
      if (vtrace) then
         ! Output F(x,y,z) for debugging purposes
         gft_shape(1) = nx
         gft_shape(2) = ny
         gft_shape(3) = nz
         mybbox(1) = gr_minx(gi)
         mybbox(3) = gr_miny(gi)
         mybbox(5) = gr_minz(gi)
         mybbox(2) = gr_maxx(gi)
         mybbox(4) = gr_maxy(gi)
         mybbox(6) = gr_maxz(gi)
!        gft_rc = gft_out_bbox( 'f_xyz', iteration*1.d0,    gft_shape,
!    *                               gft_rank, mybbox, q(gr_tmp) )
!        gft_rc = gft_out_bbox( 'theta_xyz', iteration*1.d0,gft_shape,
!    *                               gft_rank, mybbox, q(gr_flag) )
!        gft_rc = gft_out_bbox( 'phi_xyz', iteration*1.d0,  gft_shape,
!    *                               gft_rank, mybbox, q(gr_error) )
      end if

      ! Compute spatial derivative of F(x,y,z)
      ! and the expansion Theta(x,y,z)
      do k = 2, nz-1
      do j = 2, ny-1
      do i = 2, nx-1
         ! Index into 3D array
         index    = (k-1)*ny*nx + (j-1)*nx + (i-1)
         !
         index_xp = (k-1)*ny*nx + (j-1)*nx + (i-1+1)
         index_xm = (k-1)*ny*nx + (j-1)*nx + (i-1-1)
         index_yp = (k-1)*ny*nx + (j-1+1)*nx + (i-1)
         index_ym = (k-1)*ny*nx + (j-1-1)*nx + (i-1)
         index_zp = (k-1+1)*ny*nx + (j-1)*nx + (i-1)
         index_zm = (k-1-1)*ny*nx + (j-1)*nx + (i-1)
         !
         index_xpyp = (k-1)*ny*nx + (j-1+1)*nx + (i-1+1)
         index_xpym = (k-1)*ny*nx + (j-1-1)*nx + (i-1+1)
         index_xmyp = (k-1)*ny*nx + (j-1+1)*nx + (i-1-1)
         index_xmym = (k-1)*ny*nx + (j-1-1)*nx + (i-1-1)
         !
         index_xpzp = (k-1+1)*ny*nx + (j-1)*nx + (i-1+1)
         index_xpzm = (k-1-1)*ny*nx + (j-1)*nx + (i-1+1)
         index_xmzp = (k-1+1)*ny*nx + (j-1)*nx + (i-1-1)
         index_xmzm = (k-1-1)*ny*nx + (j-1)*nx + (i-1-1)
         !
         index_ypzp = (k-1+1)*ny*nx + (j-1+1)*nx + (i-1)
         index_ypzm = (k-1-1)*ny*nx + (j-1+1)*nx + (i-1)
         index_ymzp = (k-1+1)*ny*nx + (j-1-1)*nx + (i-1)
         index_ymzm = (k-1-1)*ny*nx + (j-1-1)*nx + (i-1)
         ! Get the dx,dy,dz at each point
         dx = hg
         dy = hg
         dz = hg
         ! the derivative of F is the normal vector n_i
         !
         n1 = (q(gr_tmp+index_xp) - q(gr_tmp+index_xm))/(2.0*dx) 
         n2 = (q(gr_tmp+index_yp) - q(gr_tmp+index_ym))/(2.0*dy) 
         n3 = (q(gr_tmp+index_zp) - q(gr_tmp+index_zm))/(2.0*dz) 
         ! the second derivative of F is the derivative
         ! of normal vector n_i, which are symmetric
         ! 
         d1n1 = (q(gr_tmp+index_xp) + q(gr_tmp+index_xm)
     *        -  2.0*q(gr_tmp+index) )/(dx*dx) 
         d2n2 = (q(gr_tmp+index_yp) + q(gr_tmp+index_ym)
     *        -  2.0*q(gr_tmp+index) )/(dy*dy) 
         d3n3 = (q(gr_tmp+index_zp) + q(gr_tmp+index_zm)
     *        -  2.0*q(gr_tmp+index) )/(dz*dz) 
         !
         d1n2 = (q(gr_tmp+index_xpyp) + q(gr_tmp+index_xmym)
     *        -  q(gr_tmp+index_xpym) - q(gr_tmp+index_xmyp) 
     *         )/(4.0*dx*dy) 
         d1n3 = (q(gr_tmp+index_xpzp) + q(gr_tmp+index_xmzm)
     *        -  q(gr_tmp+index_xpzm) - q(gr_tmp+index_xmzp) 
     *         )/(4.0*dx*dz) 
         d2n3 = (q(gr_tmp+index_ypzp) + q(gr_tmp+index_ymzm)
     *        -  q(gr_tmp+index_ypzm) - q(gr_tmp+index_ymzp)
     *        )/(4.0*dy*dz)               
         !
         ! Get coordinates of each point:
         x     = q(gr_x(gi)+i-1)
         y     = q(gr_y(gi)+j-1)
         z     = q(gr_z(gi)+k-1)
         if (.false.) then
             n1 = x
             n2 = y
             n3 = z
             d1n1 = 1.0
             d1n2 = 0.0
             d1n3 = 0.0
             d2n2 = 1.0
             d2n3 = 0.0
             d3n3 = 1.0
         end if
         !
#ifdef GH
         call compute_exp(n1,n2,n3,d1n1,d2n2,d3n3,d1n2,d1n3,d2n3,
     *           q(gr_mask+index),q(gr_g00+index),q(gr_g01+index),
     *           q(gr_g02+index),q(gr_g03+index),
     &           q(gr_g11+index),q(gr_g12+index),q(gr_g13+index),
     &           q(gr_g22+index),q(gr_g23+index),q(gr_g33+index),
     &           q(gr_K11+index),q(gr_K12+index),q(gr_K13+index),
     &           q(gr_K22+index),q(gr_K23+index),q(gr_K33+index),
     &           q(gr_d1g00+index),q(gr_d1g01+index),q(gr_d1g02+index),
     &           q(gr_d1g03+index),
     &           q(gr_d1g11+index),q(gr_d1g12+index),q(gr_d1g13+index),
     &           q(gr_d1g22+index),q(gr_d1g23+index),q(gr_d1g33+index),
     &           q(gr_d2g00+index),q(gr_d2g01+index),q(gr_d2g02+index),
     &           q(gr_d2g03+index),
     &           q(gr_d2g11+index),q(gr_d2g12+index),q(gr_d2g13+index),
     &           q(gr_d2g22+index),q(gr_d2g23+index),q(gr_d2g33+index),
     &           q(gr_d3g00+index),q(gr_d3g01+index),q(gr_d3g02+index),
     &           q(gr_d3g03+index),
     &           q(gr_d3g11+index),q(gr_d3g12+index),q(gr_d3g13+index),
     &           q(gr_d3g22+index),q(gr_d3g23+index),q(gr_d3g33+index), 
     &           q(gr_rad_speed+index), q(gr_rad_exp+index), kk)
#endif
         !
      end do
      end do
      end do

#ifdef GH
      if (vtrace) then
         gft_rc = gft_out_bbox( 'exp_xyz', iteration*1.d0,   gft_shape,
     *                               gft_rank, mybbox, q(gr_rad_exp) )
      end if
#endif
                
      ! Interpolate Theta(x,y,z) onto surface r=R
      do k = 1, nphi
         do j = 1, ntheta
            !
            index = (k-1)*ntheta+(j-1)
            ! Using conventions of Pretorius: gr-qc/0407110
            r     = q(horizon_R(kk)+index)
            phi   = (k-1)*dphi
            theta = (j-1)*dtheta
            ! Using conventions of Pretorius: gr-qc/0407110
            x     = horizon_center(1,kk) + r*cos(phi)*sin(theta)
            y     = horizon_center(2,kk) + r*sin(phi)*sin(theta)
            z     = horizon_center(3,kk) + r*         cos(theta)
            if(ltrace2)write(*,98)myid,'       point:',x,y,z
            !
            ! Index of lower left point of interpolation cell:
            ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
            jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
            klow  = INT( (z - gr_minz(gi)) / hg ) + 1
            if (ilow .lt. 1 .or. ilow.ge.nx .or.
     *          jlow .lt. 1 .or. jlow.ge.ny .or.
     *          klow .lt. 1 .or. klow.ge.nz ) then
               if(ltrace2)write(*,99)myid,'  NO   i/j/k:',ilow,jlow,klow
            else
               !
               if(ltrace2)write(*,99)myid,' YES   i/j/k:',ilow,jlow,klow
               ! Fractional distance point is from lower left point:
               fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
               fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
               fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg
               !
#ifdef GH
               expansion = interp(  q(gr_rad_exp),
     *                     fracx, fracy, fracz,
     *                     ilow,jlow,klow,
     *                     nx,ny,nz)
#endif
               q(horizon_theta(kk)+index) = expansion
!              if (expansion .lt. 0.d0) then
!                 write(*,98) myid,'grid_horizon_comp: Expansion neg',
!    *              x,y,z
!              end if
               !
            end if
            !
         end do
      end do

      if (ltrace) then
         write(*,99) myid, 'grid_horizon_comp: Done.   ',gi,kk,iteration
      end if
      iteration = iteration + 1

 98   format('[',I3,'] ',A,3F13.7)
 99   format('[',I3,'] ',A,3I5)
      return
      end         ! END: grid_horizon_comp

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  compute_exp:                                                              cc
cc                     compute the expansion pointwise                        cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
      subroutine compute_exp( n1,n2,n3,d1n1,d2n2,d3n3,d1n2,d1n3,d2n3,
     *           mask, g00, g01, g02, g03,
     *           g11, g12, g13, g22, g23, g33,
     *           K11, K12, K13, K22, K23, K33,
     *           D100, D101, D102, D103,
     *           D111, D112, D113, D122, D123, D133,
     *           D200, D201, D202, D203,
     *           D211, D212, D213, D222, D223, D233,
     *           D300, D301, D302, D303,
     *           D311, D312, D313, D322, D323, D333, 
     *           rad_speed, rad_exp, khole )
          
      implicit   none
      real(kind=8)   n1, n2, n3, d1n1, d2n2, d3n3, d1n2, d1n3, d2n3, 
     *               g11, g12, g13, g22, g23, g33, g00, g01, g02,
     *               g03, mask, d1g11, d1g12, d1g13,
     *               D100,D101,D102,D103,D111,D112,D113,D122,D123,D133,
     *               D200,D201,D202,D203,D211,D212,D213,D222,D223,D233,
     *               D300,D301,D302,D303,D311,D312,D313,D322,D323,D333,
     *               K11, K12, K13, K22, K23, K33,
     *               rad_speed, rad_exp
      integer        khole
      
      ! Local vars:
      real(kind=8)   x,y,z,r,h11,h12,h13,h22,h23,h33,
     *               huu11,huu12,huu13,huu22,huu23,huu33,
     *               deth, b1, b2, b3, alp, alp2,
     *               eK11,eK12,eK13,eK22,eK23,eK33,etrK,     
     *               nr1, nr2, nr3, nru1, nru2, nru3, nrnorm,
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
     *               d2nr3,d3nr1,d3nr2,d3nr3, besafe, tempval

      logical     ltrace
      parameter ( ltrace = .false.)
       
      rad_speed = 0.0
      rad_exp = 0.0
      besafe = 1.
               
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
         besafe = -1.
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
        besafe = -1.
      end if 
      alp  = sqrt(alp2)
               
      !
      ! compute the normal vector
      !
      nr1    = n1
      nr2    = n2
      nr3    = n3

      nrnorm = sqrt( huu11*nr1*nr1 + huu22*nr2*nr2
     *      + huu33*nr3*nr3 + 2.0d0*(huu12*nr1*nr2 + huu13*nr1*nr3
     *      + huu23*nr2*nr3))
       
      if(nrnorm.lt.1.e-8) then
         nrnorm = 1.e-8
         besafe = -1.
      end if
 
      nr1 = nr1/nrnorm
      nr2 = nr2/nrnorm
      nr3 = nr3/nrnorm      

      nru1 = huu11*nr1 + huu12*nr2 + huu13*nr3
      nru2 = huu12*nr1 + huu22*nr2 + huu23*nr3
      nru3 = huu13*nr1 + huu23*nr2 + huu33*nr3     
      !
      ! its unnormalized derivative
      !
      d1nr1 = d1n1; d2nr2 = d2n2; d3nr3 = d3n3
      d1nr2 = d1n2; d1nr3 = d1n3; d2nr3 = d2n3
      d2nr1 = d1n2; d3nr1 = d1n3; d3nr2 = d2n3
  
      !
      !
      ! compute the expansion
      !
     
      !  the extrinsic curvature
               
      eK11 = -(-alp*K11 - D101 - D101 
     *     + b1*(D111 + D111) + b2*(D112 + D112) 
     *     + b3*(D113 + D113))/(2.0*alp)
      eK12 = -(-alp*K12 - D102 - D201 
     *     + b1*(D112 + D211) + b2*(D122 + D212) 
     *     + b3*(D123 + D213))/(2.0*alp)
      eK13 = -(-alp*K13 - D103 - D301 
     *     + b1*(D113 + D311) + b2*(D123 + D312) 
     *     + b3*(D133 + D313))/(2.0*alp)
      eK22 = -(-alp*K22 - D202 - D202 
     *     + b1*(D212 + D212) + b2*(D222 + D222) 
     *     + b3*(D223 + D223))/(2.0*alp)
      eK23 = -(-alp*K23 - D203 - D302 
     *     + b1*(D213 + D312) + b2*(D223 + D322)
     *     + b3*(D233 + D323))/(2.0*alp)
      eK33 = -(-alp*K33 - D303 - D303 
     *     + b1*(D313 + D313) + b2*(D323 + D323) 
     *     + b3*(D333 + D333))/(2.0*alp)
      
      etrK = huu11*eK11 + huu12*eK12 + huu13*eK13
     *     + huu12*eK12 + huu22*eK22 + huu23*eK23
     *     + huu13*eK13 + huu23*eK23 + huu33*eK33

      ! the Christophel symbols with the first index up      
            
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

      ! contractions with n^i
      
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

      rad_exp = tempval

      return
      end subroutine compute_exp

      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  cartospher:   Convert cartesian coordiantes to spherical.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cartospher(r,theta,phi, x,y,z)
      implicit     none
      real(kind=8) r, theta, phi
      real(kind=8) x, y, z
      include     'largesmall.inc'
      !
      real(kind=8) cpi
      parameter  ( cpi  =  3.14159 26535 89793 23846 26433 83279 d0)

      r = sqrt( x**2 + y**2 + z**2 )

      if (r .gt. SMALLNUMBER) then
      !if (r .gt. 0) then
         theta = acos( z / r )
      else
         theta = 0.d0
      end if

      if (x.gt.SMALLNUMBER) then
      !if (x.gt.0) then
         phi = atan( y / x )
      else if (x.lt.-SMALLNUMBER) then
      !else if (x.lt.0) then
         phi = atan( y / x ) + cPi
      else if (y.lt.-SMALLNUMBER) then
      !else if (y.lt.0) then
         phi = 1.5d0*cPi
      else
         phi = 0.5d0*cPi
      end if

      !
      ! Make sure angle is in interval [0,2Pi]
      !
      !if (phi .lt. 0) phi   = phi + 2.d0*cPi
      if (phi .lt. -SMALLNUMBER) then
         phi = phi + 2.d0*cPi
      else if (phi .lt. SMALLNUMBER) then
         phi = 0.d0
      end if
      if (phi .gt. 2.d0*cPi+SMALLNUMBER) then
         phi = phi - 2.d0*cPi
      else if (phi .gt. 2.d0*cPi-SMALLNUMBER) then
         phi = 2.d0*cPi
      end if

      return
      end subroutine cartospher

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp2dNU:                                                               cc
cc              Bilinearly interpolate in 2d on a nonuniform                  cc
cc              four point stencil.                                           cc
cc              NOTE: this turns out to be a somewhat difficult problem       cc
cc              known as irregular interpolation. Since we don't care too     cc
cc              much about accuracy (since we'll be looking for a new horizon cc
cc              again anyway), just using "max" (since horizon moves in).     cc
cc                Points labeled as:                                          cc
cc                                                                            cc
cc                        f3          f4                                      cc
cc                                                                            cc
cc                        f1          f2                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interp2dNU( x, y,
     *                                      f4, x4, y4,
     *                                      f3, x3, y3,
     *                                      f2, x2, y2,
     *                                      f1, x1, y1     )
      implicit     none
      real(kind=8) x,y, f4, x4, y4, f3,x3,y3, f2,x2,y2, f1,x1,y1

      logical     ltrace
      parameter ( ltrace = .false. )

!     interp2d =
!    *      + ( 1.d0 - (x-x1)/x1) * (1.d0 - (y-y1)/y1) * f1
!    *      + ( 1.d0 - (x-x1)/x1) * (1.d0 - (y-y1)/y1) * f3
!    *      + (1.d0-fx) * (1.d0-fy) * f1
!    *      +       fx  * (1.d0-fy) * f2
!    *      + (1.d0-fx) *       fy  * f3
!    *      +       fx  *       fy  * f4
      interp2dNU = max(f1,f2,f3,f4)

      return
      end    ! END: interp2dNU

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  horizon_recenter:                                                         cc
cc              Determine the new R(theta,phi) when a new center is used.     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine horizon_recenter( newR, oldR, newC, oldC, dtheta, dphi,
     *                     ntheta, nphi)
      implicit     none
      integer      ntheta, nphi
      real(kind=8) newR(ntheta,nphi), oldR(ntheta,nphi)
      real(kind=8) newC(3), oldC(3), dtheta, dphi
      !
      include     'largesmall.inc'
      !
      integer      j,k, myid, num
      integer      nk,nj, ni
      real(kind=8) p,t, p1, t1, r_t1, r_p1, r_p1t1
      real(kind=8) x_p1t1, y_p1t1, z_p1t1, x_p1, y_p1, z_p1,
     *             x_t1, y_t1, z_t1, nx,ny,nz, nx_t1,ny_t1,nz_t1,
     *             nx_p1,ny_p1,nz_p1, nx_p1t1,ny_p1t1,nz_p1t1,
     *             nr,nt,np, nr_p1,nr_t1,nr_p1t1, np1,nt1,
     *             r, x, y, z, nth, nph
      real(kind=8) cpi
      parameter  ( cpi  =  3.14159 26535 89793 23846 26433 83279 d0)
      integer      proc_return_myid
      external     proc_return_myid
      real(kind=8) interp2dNU
      external     interp2dNU
      logical      ltrace
      parameter  ( ltrace  = .false. )
      logical      ltrace2            ! output on all procs
      parameter  ( ltrace2 = .false. )

      myid = proc_return_myid()
      if (ltrace2.or.(myid.eq.0.and.ltrace)) then
         write(*,98)myid,'horizon_recenter: Enter'
         write(*,98)myid,'horizon_recenter:newC',newC(1),newC(2),newC(3)
         write(*,98)myid,'horizon_recenter:oldC',oldC(1),oldC(2),oldC(3)
      end if

      ! So we can tell if any points didn't get re-set:
      if(ltrace)call load_scal1D(newR,LARGENUMBER, ntheta*nphi)

      ! Loop over interpolation square:
      do k = 1, nphi-1
      do j = 1, ntheta-1
         !
         ! For each square, define the four points:
         !   (r_t1,p,t1), (r_p1t1,p1,t1)
         !   (r,   p,t),  (r_p1,  p1,t)
         !
         p     = (k-1)*dphi
         p1    = (k-1+1)*dphi
         !
         t     = (j-1)*dtheta
         t1    = (j-1+1)*dtheta
         !
         r     = oldR(j,  k  )
         r_p1  = oldR(j,  k+1)
         r_t1  = oldR(j+1,k  )
         r_p1t1= oldR(j+1,k+1)
         ! Using conventions of Pretorius: gr-qc/0407110
         x     = oldC(1) + r     *cos(p )*sin(t )
         y     = oldC(2) + r     *sin(p )*sin(t )
         z     = oldC(3) + r     *        cos(t )
         x_p1t1= oldC(1) + r_p1t1*cos(p1)*sin(t1)
         y_p1t1= oldC(2) + r_p1t1*sin(p1)*sin(t1)
         z_p1t1= oldC(3) + r_p1t1*        cos(t1)
         x_p1  = oldC(1) + r_p1  *cos(p1)*sin(t )
         y_p1  = oldC(2) + r_p1  *sin(p1)*sin(t )
         z_p1  = oldC(3) + r_p1  *        cos(t )
         x_t1  = oldC(1) + r_t1  *cos(p )*sin(t1)
         y_t1  = oldC(2) + r_t1  *sin(p )*sin(t1)
         z_t1  = oldC(3) + r_t1  *        cos(t1)
         ! Express point (x,y,z) in terms of new center:
         !    (All such descriptors associated with
         !     the new center have a prefix of "n")
         nx      = x      - newC(1)
         nx_p1   = x_p1   - newC(1)
         nx_t1   = x_t1   - newC(1)
         nx_p1t1 = x_p1t1 - newC(1)
         ny      = y      - newC(2)
         ny_p1   = y_p1   - newC(2)
         ny_t1   = y_t1   - newC(2)
         ny_p1t1 = y_p1t1 - newC(2)
         nz      = z      - newC(3)
         nz_p1   = z_p1   - newC(3)
         nz_t1   = z_t1   - newC(3)
         nz_p1t1 = z_p1t1 - newC(3)
         ! Find new spherical coordinates w/r/t new center:
         call cartospher(nr,     nt, np,  nx,     ny,     nz     )
         call cartospher(nr_p1,  nt, np1, nx_p1,  ny_p1,  nz_p1  )
         call cartospher(nr_t1,  nt1,np,  nx_t1,  ny_t1,  nz_t1  )
         call cartospher(nr_p1t1,nt1,np1, nx_p1t1,ny_p1t1,nz_p1t1)
         if (np1 .lt. np) np1 = np1 + 2.d0*cPi
         !
         if (ltrace2) then
            !write(*,*)'horizon_recenter:x/y/z     ',x,     y,     z
            !write(*,*)'horizon_recenter:x/y/z_p1  ',x_p1,  y_p1,  z_p1
            !write(*,*)'horizon_recenter:x/y/z_t1  ',x_t1,  y_t1,  z_t1
            !write(*,*)'horizon_recenter:x/y/z_p1t1',x_p1t1,y_p1t1,z_p1t1
            !write(*,*) 'horizon_recenter: From square:'
            !write(*,*) 'horizon_recenter:  p-->p1:  ', p,  p1
            !write(*,*) 'horizon_recenter:  t-->t1:  ', t,  t1
            !write(*,*) 'horizon_recenter:  p,  t:  ', p,  t
            !write(*,*) 'horizon_recenter:  p1, t:  ', p1, t
            !write(*,*) 'horizon_recenter:  p,  t1: ', p,  t1
            !write(*,*) 'horizon_recenter:  p1, t1: ', p1, t1
            write(*,*) 'horizon_recenter: Interpolating in square:'
            !write(*,*) 'horizon_recenter: np, nt:  ',np, nt
            !write(*,*) 'horizon_recenter: np1,nt:  ',np1,nt
            !write(*,*) 'horizon_recenter: np, nt1: ',np, nt1
            !write(*,*) 'horizon_recenter: np1,nt1: ',np1,nt1
            !write(*,*) 'horizon_recenter: '
            write(*,*) 'horizon_recenter: nt-->nt1:  ', nt,  nt1
            write(*,*) 'horizon_recenter: np-->np1:  ', np,  np1
            write(*,*) 'horizon_recenter: '
         end if
         !
         ! Interpolate recentered radii from
         ! points within this square
         !
         ! Find indices for first point to be interpolated: (nj,nk)
         !
         nk =  np / dphi   + 1
         nph= (nk-1)*dphi
         if ( nph .lt. np-SMALLNUMBER ) nk = nk + 1
  11     nj =  nt / dtheta + 1
         nth= (nj-1)*dtheta
         if ( nth .lt. nt-SMALLNUMBER ) nj = nj + 1
         !
         nph= (nk-1)*dphi
         ! Check that point's angles (nth,nph) fall within the square:
         if ( (nph.ge.np-SMALLNUMBER) .and.(nph.le.np1+SMALLNUMBER))then
  12        nth= (nj-1)*dtheta
            if((nth.ge.nt-SMALLNUMBER).and.(nth.le.nt1+SMALLNUMBER))then
               if(ltrace2)write(*,96) 'Set:nj/k:',nj,nk,nth,nph
               newR(nj,nk) = interp2dNU(nth, nph,
     *                                  nr_p1t1, nt1, np1,
     *                                  nr_t1,   nt,  np1,
     *                                  nr_p1,   nt1, np,
     *                                  nr,      nt,  np    )
               ! go to next point
               nj = nj + 1
               if (nj.le.ntheta) goto 12
            else if (ltrace2) then
               write(*,96) 'NopeJ:nj/k:',nj,nk,nth,nph
            end if
            nk = nk + 1
            if (nk.le.nphi) goto 11
         else if (ltrace2) then
               write(*,96) 'NopeK:nj/k:',nj,nk,nth,nph
         end if
         !
      end do
      end do

      ! Enforce boundary conditions:
      k = nphi
      do j = 1, ntheta
         ! identify azimuthal sides
         newR(j,nphi) = newR(j,1)
      end do
      ! Extrapolate to poles if not interpolated:
      if (newR(ntheta,1).ge.LARGE_M_SMALL) then
         newR(ntheta,1) = 2.d0*newR(ntheta-1,1)-newR(ntheta-2,1)
      end if
      if (newR(1,1)     .ge.LARGE_M_SMALL) then
         newR(1,1)      = 2.d0*newR(2,1)       -newR(3,1)
      end if
      ! Identify all points along poles as the same point
      do k = 2, nphi
         ! top:
         newR(ntheta,k) = newR(ntheta,1)
         ! bottom:
         newR(1,     k) = newR(1,1)
      end do

      if (ltrace.or. .true.) then
         ! Check if all points got set:
         if(ltrace)write(*,*)'horizon_recenter:ntheta/phi: ',ntheta,nphi
         do k = 1, nphi
         do j = 1, ntheta
            if (newR(j,k).ge.LARGE_M_SMALL) then
               write(*,96) 'Unset j/k',j,k,(j-1)*dtheta,(k-1)*dphi
               ! If not set, try to extrapolate/average
               num = 0
               newR(j,k) = 0.d0
               if (j.lt.ntheta) then
                  if (newR(j+1,k).lt.LARGE_M_SMALL) then
                     newR(j,k) = newR(j,k) + newR(j+1,k)
                     num = num + 1
                  end if
               end if
               if (j.gt.1) then
                  if (newR(j-1,k).lt.LARGE_M_SMALL) then
                     newR(j,k) = newR(j,k) + newR(j-1,k)
                     num = num + 1
                  end if
               end if
               if (k.lt.nphi) then
                  if (newR(j,k+1).lt.LARGE_M_SMALL) then
                     newR(j,k) = newR(j,k) + newR(j,k+1)
                     num = num + 1
                  end if
               end if
               if (k.gt.1) then
                  if (newR(j,k-1).lt.LARGE_M_SMALL) then
                     newR(j,k) = newR(j,k) + newR(j,k-1)
                     num = num + 1
                  end if
               end if
               if (num.gt.0) then
                     newR(j,k) = newR(j,k)/num
               else
                     write(*,*) 'Problem, could not extrap: ',j,k
               end if
            end if
         end do 
         end do 
      end if

      if (ltrace2.or.(myid.eq.0.and.ltrace)) then
         write(*,99) myid, 'horizon_recenter: Done'
      end if

 96   format(A,2I5,3F13.5)
 97   format('[',I3,'] ',A,I5,3F13.5)
 98   format('[',I3,'] ',A,3F13.5)
 99   format('[',I3,'] ',A,2I5)

      return
      end    ! END: horizon_recenter
