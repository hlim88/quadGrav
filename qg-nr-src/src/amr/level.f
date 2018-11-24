cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_refine:                                                             cc
cc                Refine level by creating grids at level+1.                  cc
cc                If level==-1, then we're creating the coarse level.         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_refine( level )
      implicit none
      integer       level
      include      'grid_methods.inc'
      include      'grid.inc'
      include      'glob.inc'
      include      'action.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      include      'output.inc'
      include      'param.inc'
      include      'largesmall.inc'
      integer      nxl,nyl,nzl, gi, gi_tmp
      integer      flagl, flagk
      integer      owner, proc, nextproc, action, rc
      real(kind=8) minx, maxx, miny, maxy, minz, maxz, hl
      real(kind=8) minxc,      minyc,      minzc
      real(kind=8) bbox(6),    time,       tmp_ref_level
      integer      nx,         ny,         nz
      integer      mini,       minj,       mink
      integer      i, numbox, sizemem
      integer      nxc,nyc,nzc
      real(kind=8)  myl2norm3d
      external      myl2norm3d
      ! If user chose ref_level<0, then user wants FMR:
      logical      usingfmr

      !
      ! Bounding boxes of sub grids in refinement process:
      !
      integer  bbox_minx(maxnumchildren), bbox_maxx(maxnumchildren),
     *         bbox_miny(maxnumchildren), bbox_maxy(maxnumchildren),
     *         bbox_minz(maxnumchildren), bbox_maxz(maxnumchildren)
      integer  child_procassignments(maxnumchildren)

      ! To store "recv_action" parameters:
      integer      gparent, gnx, gny, gnz, gax,gay,gaz,
     *             gfrom,   glevel, ggi, gowner

      logical      proc_havemem
      external     proc_havemem
      integer      proc_return_myid,   proc_pick,
     *             mem_alloc
      external     proc_return_myid,   proc_pick,
     *             mem_alloc
      integer      tmp_numprocs
      ! Determine if complete overlap for many levels;
      integer      tmplev
      logical      totalrefining

      !
      ! When one runs out of memory, try and rescue:
      !
      ! Factor by which to increase ethreshold when not enough mem:
      real(kind=8) INCREASEFACTOR
      parameter (  INCREASEFACTOR = 2.d0 )
      ! Number of times to try increasing ethreshold
      integer         numretries    
      integer      MAXNUMRETRIES    
      parameter (  MAXNUMRETRIES  = 5    )
      !If you don't want any retries:
      !parameter (  MAXNUMRETRIES  = -10    )
      ! If you don't want any backouts and just want to quit
      ! if not enough memory
      logical      nokeepalive
      parameter (  nokeepalive  = .false.  )
      logical      nogridsonmaster
      parameter  ( nogridsonmaster  = .false. )

      logical      ltrace
      parameter (  ltrace  = .false.)
      logical      ltrace2
      parameter (  ltrace2 = .false.)
      logical      outputerr
      parameter (  outputerr = .false.)


      numretries = 0

      if (ltrace) then
         write(*,*)'level_refine: Refining level: ',level,'***********'
         write(*,*)'level_refine: Refining level: allowedl=',allowedl
      end if
      if (nogridsonmaster) then
         write(*,*)'level_refine: running with nogridsonmaster set'
      end if

      !
      ! Don't refine if we've hit our level limit:
      !
      !if (level .ge. allowedl) return
      !if (level .eq. allowedl) return
      ! SLL: After checkpointing, allowedl may be decreased by the user,
      !      in which case a finer level may exist.
      !      For example, say levels 0..9 exit, but then allowedl is
      !      changed to 6. In that case, one would not re-refine
      !      on level 9, but one would on other levels.
      !      And so we bypass re-refinement only if we're at
      !      the allowedl limit *and* no finer level already exists
      !if (level .ge. allowedl .and. .not.level_exists(level+1)) return
      if (level .ge. allowedl      .and.
     *    level .eq. level_return_finest(Levelp,maxlev)) return
      if (level .ge. allowedl) then
         write(*,*)'level_refine: Finer levels than allowed'
         write(*,*)'level_refine: Zeroing out TRE'
         ref_level(level) = 1.d-14
      end if

      !
      ! Find physical bounds of level:
      !
      if (level.eq.-1) then
         !
         ! If level==-1, we're creating the coarse grid:
         !
         minx = minx0
         miny = miny0
         minz = minz0
         maxx = maxx0
         maxy = maxy0
         maxz = maxz0
         gi   =  -1
         hl   = h * refine_factor
         time = 0.d0
      else
         if (ltrace2) call level_tree_dump()
         !
         ! Find bounds of level:
         !
         call level_find_bounds(level,minx,maxx,miny,maxy,minz,maxz)
         !
         ! Find beginning of level
         !
         gi = level_return_start(level)
         !
         ! Grid spacing for the level:
         !
         hl    = grid_return_h(gi)
         !
         !  Time on level (just for tracing/output)
         !
         time         = grid_return_time(level_return_start(level))
      end if

      !
      ! Get my processor id:
      !
      myid  = proc_return_myid()

      !
      ! Dimensions of bounding box for entire level:
      !
      nxl   = NINT( (maxx-minx)/hl ) + 1
      nyl   = NINT( (maxy-miny)/hl ) + 1
      nzl   = NINT( (maxz-minz)/hl ) + 1
      if (ltrace) then
         write(*,*) 'level_refine: refine_factor  =',refine_factor
         write(*,*) 'level_refine: ethreshold     =',ethreshold
         write(*,*) 'level_refine: numretries     =',numretries
         write(*,*) 'level_refine: shadow         =',shadow
         write(*,*) 'level_refine: h  =',h
         write(*,*) 'level_refine: hl =',hl
         write(*,*) 'level_refine: max/inx=',maxx,minx
         write(*,*) 'level_refine: max/iny=',maxy,miny
         write(*,*) 'level_refine: max/inz=',maxz,minz
         write(*,*) 'level_refine: nx/y/zl=',nxl,nyl,nzl
      end if

      !
      ! Allocate memory for flag array:
      !    flagl  is used to store the error values from all the grids,
      !           and after that it is used for temporary storage that
      !           some of the routines need
      !    flagk  is used as the flag array for the entire level
      !
      !    NB: these both get deallocated before the return statement
      !
      sizemem = nxl*nyl*nzl
      if (6*maxnumchildren .gt. sizemem) then
         if (ltrace) write(*,*) 'level_refine: Using 6*memchildren',
     *                  sizemem, 6*maxnumchildren
         sizemem = 6*maxnumchildren
      else
         if (ltrace) write(*,*) 'level_refine: Using nxl*nyl*nzl',
     *                  sizemem, 6*maxnumchildren
      end if
      !
      if (.not. proc_havemem(myid, sizemem)) then
         write(*,*) 'level_refine: Not enough memory to'
         write(*,*) '              (re)refine this level: ',level
         write(*,*) 'level_refine: nx/y/zl: ',nxl,nyl,nzl
         write(*,*) 'level_refine: Need to allocate sizemem: ',sizemem
         call mem_stat()
         call proc_dump_wkld()
         write(*,*) 'level_refine: Returning...'
         return
      end if
      if (ltrace)write(*,*)'level_refine: Allocating sizemem: ',sizemem
      flagl   = mem_alloc( sizemem     )
      call proc_addmem(myid, sizemem    )
      !
      if (.not. proc_havemem(myid, nxl*nyl*nzl)) then
         write(*,*) 'level_refine: Not enough memory to'
         write(*,*) '              (re)refine this level: ',level
         write(*,*) 'level_refine: nx/y/zl: ',nxl,nyl,nzl
         write(*,*) 'level_refine: Need to allocate nxl*nyl*nzl: ',
     *                  nxl*nyl*nzl
         call mem_stat()
         call proc_dump_wkld()
         write(*,*) 'level_refine: Returning...'
         return
      end if
      if(ltrace)write(*,*)'level_refine:Allocat nxl*nyl*nzl',nxl*nyl*nzl
      flagk   = mem_alloc( nxl*nyl*nzl )
      call proc_addmem(myid, nxl*nyl*nzl)

      !
      ! If user wants to read in cluster info from a previous run:
      !
      if (clusterreadwrite .eq. 1) then
         if (ltrace) write(*,*) 'level_refine: Reading in cluster info'
         call level_bbox_read(level, bbox_minx,bbox_maxx,
     *                               bbox_miny,bbox_maxy,
     *                               bbox_minz,bbox_maxz,
     *                               minx,miny,minz,hl,
     *                               time, numbox, assume_symmetry,rc)
         if(ltrace)write(*,*)'level_refine: rc from bbox_read:',rc
         ! Check status code:
         !   ...if it found an entry, then proceed normally
         if (rc.gt.0) call level_mkall_dead(level+1)
         !   ...else, clean things and return w/ no change to hierarchy
         goto 20
      end if


      ! If we end up running out of memory, we retry here with
      ! a different error threshold:
 33   continue

      !
      ! Initialize to negative value (since all errors should be
      ! positive): (negative value is arbitrary)
      !
      call load_scal3d( q(flagl), -2.0001d0, nxl, nyl, nzl)

      if (ltrace2) then
         write(*,*) 'level_refine:*** Refining level:',level
         write(*,*) 'level_refine:    minx / maxx = ',minx,maxx
         write(*,*) 'level_refine:    miny / maxy = ',miny,maxy
         write(*,*) 'level_refine:    minz / maxz = ',minz,maxz
         write(*,*) 'level_refine:    hl          = ',hl
         write(*,*) 'level_refine:    nxl,nyl,nzl = ',nxl,nyl,nzl
         write(*,*) 'level_refine:    flagl array = ',flagl
         write(*,*) 'level_refine:    flagk array = ',flagk
      end if

      !
      ! Mark all grids in next level as DEAD
      !     The clustering process will tell us
      !     which grids we should have. We might
      !     bring some of the DEAD ones back to life,
      !     but whichever ones stay dead after this,
      !     will then be freed. Only the local proc
      !     is aware that they are dead. The other procs
      !     only find out about which ones should be freed.
      !
      if (ltrace2) write(*,*) 'level_refine: Killing childgrids locally'
      !
      call level_mkall_dead(level+1)

      !
      ! For creating coarse grid or for FMR, can skip much of this stuff:
      !
      usingfmr = .false.
      if (level.ge.0) then
         if(ltrace2)write(*,*)'level_refine:ref_level=',ref_level(level)
         usingfmr = ref_level(level) .le. 0
         if(ltrace2)write(*,*)'level_refine: usingfmr ',usingfmr
      else
         usingfmr = .false.
         if(ltrace2)write(*,*)'level_refine: Not usingfmr ',usingfmr
      end if
      totalrefining = .true.
      if (level.eq.-1) goto 133
      tmplev        = level
 132  if (totalrefining  .and. ref_level(tmplev).eq.-1.0) then
         tmplev = tmplev - 1
         if (tmplev .ge. 0) goto 132
      else
         if(ltrace2)write(*,*)'level_refine: not totalrefining:',tmplev,
     *                      ref_level(tmplev)
         totalrefining = .false.
      end if
      if(ltrace)write(*,*)'level_refine: totalrefining:',totalrefining
 133  if (level.eq.-1 .or.
     *    totalrefining .or.
     *    (level.eq.0.and.shadow.ne.0.and. .not.usingfmr) ) then
         if (ltrace2) write(*,*) 'level_refine: Refine whole domain'
         numbox = 1
         bbox_minx(1) = 1
         bbox_maxx(1) = nxl
         bbox_miny(1) = 1
         bbox_maxy(1) = nyl
         bbox_minz(1) = 1
         bbox_maxz(1) = nzl
         goto 15
      end if

      !
      ! Loop over all grids in level 
      !    and get its error data:
      !
      gi_tmp = level_return_start(level)
      call send_actionAll(GETERROR,gi_tmp,0,level,0,nx,ny,nz,0,0,0)
 10   if (grid_return_existence(gi_tmp)) then
         !
         ! Find "anchor" points to level field:
         !
         mini = 1 + NINT( (gr_minx(gi_tmp)-minx)/hl )
         minj = 1 + NINT( (gr_miny(gi_tmp)-miny)/hl )
         mink = 1 + NINT( (gr_minz(gi_tmp)-minz)/hl )
         nx   = gr_nx(gi_tmp)
         ny   = gr_ny(gi_tmp)
         nz   = gr_nz(gi_tmp)
         !
         owner  = grid_return_owner(gi_tmp)
         if (owner .eq. myid) then
            if(ltrace)write(*,*)'level_refine:get error locally:',gi_tmp
            call grid_mask_error(gi_tmp)
            call load_pointers(gi_tmp)
            call level_combine(q(flagl), q(gr_error),
     *                         mini,minj,mink,
     *                         nxl,nyl,nzl,
     *                         nx, ny, nz  )
         else
            if(ltrace)write(*,*)'level_refine:get error remotel:',gi_tmp
            call MPI_Recv(q(flagk), nx*ny*nz, MPI_DOUBLE_PRECISION,
     *               owner, TAG_ERROR_DATA, MPI_COMM_WORLD, status,ierr)
            call level_combine(q(flagl), q(flagk),
     *                         mini,minj,mink,
     *                         nxl,nyl,nzl,
     *                         nx, ny, nz  )
         end if
         gi_tmp = grid_return_sibling(gi_tmp)
         goto 10
      end if

      !
      ! Output error array?
      !
      if (ltrace) then
         write(*,*) 'level_refine: Outputting error at time: ',time
         call field_out3d(q(flagl),time,'amrperror',minx,maxx,
     *                    miny,maxy,minz,maxz,nxl,nyl,nzl,0)
      end if

      !
      !if(ltrace)write(*,*) 'level_refine: Smoothing error....'
      ! This smoothing introduces a strange assymmetry and isn't
      ! necessary, so I'm commenting it out (sll 11/3/09)
      !call smooth3D(q(flagl),q(flagk),nxl,nyl,nzl)
     
      if (level.eq.
     *    level_return_finest(Levelp,maxlev)) then
         ! Output l2norm/min/max of error:
         call field_dump_infob(q(flagl),nxl,nyl,nzl,'Error: ')
      end if
      if (myl2norm3d(q(flagl),nxl,nyl,nzl) .gt.10.d0) then
         write(*,*)'level_refine: something wrong with the error'
         write(*,*)'level_refine: Very likely the code has  NaNs'
         write(*,*)'level_refine: if so, then the clusterer is'
         write(*,*)'level_refine: almost certainly going to mess up'
         write(*,*)'level_refine: Fix the evolution first!'
         write(*,*)'level_refine: Below is the norm of the error:'
         call field_dump_infob(q(flagl),nxl,nyl,nzl,'Error: ')
      end if
      if (ltrace.or.(outputerr.and. .not.usingfmr)) then
         write(*,*) 'level_refine: Level/time: ',level,time
         call field_dump_infob(q(flagl),nxl,nyl,nzl,'Error: ')
         call field_out2d(q(flagl+((nzl+1)/2-1)*nyl*nxl),time,
     *                    'amrerror2D',minx,maxx,
     *                    miny,maxy,nxl,nyl,0)
         call field_out3d(q(flagl),time,'amrerror',minx,maxx,
     *                    miny,maxy,minz,maxz,nxl,nyl,nzl,0)
      end if


      !
      ! Make flag array:
      !
 17   if (ltrace2) write(*,*) 'level_refine: Making flag array'
      if (usingfmr) then
         !
         ! Implement simple FMR by flagging central
         ! fraction of area: 1/-ref_level(level)
         !
         if(ltrace)write(*,*)'level_refine:Flagging fixed central regio'
         if(ltrace)write(*,*)'level_refine:ref_level: ',ref_level(level)
         call load_scal3d(q(flagk), FLAG_DISALLOW, nxl,nyl,nzl )
         if (ref_level(level) .ge. -0.5d0) then
            !
            ! This was done for Jan to resolve jet-like features along z:
            write(*,*)'level_refine: specialized FMR ',ref_level(level)
            !call load_scalallz( q(flagk), -ref_level(level),
            call load_scalallz( q(flagk), level,
     *             hl, minx,miny,minz, FLAG_REFINE, nxl, nyl, nzl,
     *             nx0,ny0,nz0)
            !
         else if (ref_level(level) .gt. -1) then
            !
            ! This was done for Eric *not* to resolve jet-like features along z:
            write(*,*)'level_refine: specialized FMR ',ref_level(level)
            call load_scalallz2( q(flagk), -ref_level(level),
     *             hl, minx,miny,minz, FLAG_REFINE, nxl, nyl, nzl)
            !
         else
            !
            tmp_ref_level = -ref_level(level)
            do i=1,numretries
               tmp_ref_level = tmp_ref_level*INCREASEFACTOR
                write(*,*)'level_refine:tmp_ref_level=',tmp_ref_level
            end do
            if(ltrace2)write(*,*)'level_refine:calling load_scalfrac3d'
            ! 2 coarse grids give us an extended boundary of four points
            call load_scalfrac3d( q(flagk), tmp_ref_level,
            !call load_scalfrac3d( q(flagk), -ref_level(level),
     *             hl, minx,miny,minz,
     *             FLAG_REFINE, assume_symmetry, 2, nxl, nyl, nzl)
            if(ltrace2)write(*,*)'level_refine:call load_scal_bound3dW'
            call load_scal_bound3dW( q(flagk), FLAG_DISALLOW,
     *             ghostwidth, assume_symmetry,nxl, nyl, nzl)
            !
         end if
      else
         !
         ! Normal: flag points where error above threshold
         !
         if(ltrace)write(*,*)'level_refine:Rescale by:',ref_level(level)
         call load_scal_mult3d( q(flagl), q(flagl), ref_level(level),
     *                                                    nxl,nyl,nzl)
         if(ltrace)write(*,*)'level_refine:Flagging error exceede',
     *         ethreshold
         call level_makeflag(q(flagk),q(flagl),level, minx,miny,minz,hl,
     *                       nxl,nyl,nzl,time)
      end if


      if (ltrace.or.outputerr) then
         call field_dump_infob(q(flagk),nxl,nyl,nzl,'FlagK: ')
         call field_out3d(q(flagk),time,'amrflag',minx,maxx,
     *                    miny,maxy,minz,maxz,nxl,nyl,nzl,0)
      end if

      !
      ! Using flagl for extra storage for the clusterer:
      !
      if (4*(nxl+nyl+nzl) .gt. nxl*nyl*nzl ) then
         write(*,*) 'level_refine: Memory problem!',nxl,nyl,nzl
         call my_exit('level_refine: memory problem')
      end if
      if(ltrace)write(*,*)'level_refine: Calling level_cluster()'
      call level_cluster(q(flagk),
     *                   q(flagl                  ),
     *                   q(flagl+  nxl            ),
     *                   q(flagl+  nxl+  nyl      ),
                         !
     *                   q(flagl+  nxl+  nyl+  nzl),
     *                   q(flagl+2*nxl+  nyl+  nzl),
     *                   q(flagl+2*nxl+2*nyl+  nzl),
                         !
     *                   q(flagl+2*nxl+2*nyl+2*nzl),
     *                   q(flagl+3*nxl+2*nyl+2*nzl),
     *                   q(flagl+3*nxl+3*nyl+2*nzl),
                         !
     *                   q(flagl+3*nxl+3*nyl+3*nzl),
     *                   q(flagl+4*nxl+3*nyl+3*nzl),
     *                   q(flagl+4*nxl+4*nyl+3*nzl),
     *                   time,
     *                   bbox_minx,bbox_maxx,
     *                   bbox_miny,bbox_maxy,
     *                   bbox_minz,bbox_maxz,
     *                   minx,maxx,
     *                   miny,maxy,
     *                   minz,maxz,
     *                   numbox, nxl,nyl,nzl )


 15   continue
      if (ltrace2) then
         write(*,*) 'level_refine: Done clustering, level: ',level
         write(*,*) 'level_refine: Done clustering, numbox:',numbox
         do i = 1, numbox
            write(*,*) 'level_refine: i:',i
            write(*,*) 'level_refine: x:',bbox_minx(i),bbox_maxx(i)
            write(*,*) 'level_refine: y:',bbox_miny(i),bbox_maxy(i)
            write(*,*) 'level_refine: z:',bbox_minz(i),bbox_maxz(i)
            if (bbox_maxx(i) .lt. bbox_minx(i)) then
               call my_exit('level_refine: Bad box in x')
            end if
            if (bbox_maxy(i) .lt. bbox_miny(i)) then
               call my_exit('level_refine: Bad box in y')
            end if
            if (bbox_maxz(i) .lt. bbox_minz(i)) then
               call my_exit('level_refine: Bad box in z')
            end if
         end do
      end if

      !
      ! If user wants to write out cluster info 
      !    (NB: write out clusters before decomposing them)
      !
      if (clusterreadwrite .eq. 2) then
         if (ltrace) write(*,*) 'level_refine: Writing cluster info'
         call level_bbox_write(level,bbox_minx,bbox_maxx,
     *                               bbox_miny,bbox_maxy,
     *                               bbox_minz,bbox_maxz,
     *                               minx,miny,minz,hl,
     *                               time, numbox)
      end if

      !
      ! The following decomposes the grids called for above.
      ! There are cases where we may not want that, at least
      ! for the coarse grid.
      !    clusterDD == 0   cluster as per usual
      !    clusterDD == 1   do not cluster the coarse grid
      !    clusterDD >  1   do not cluster at any level
      !
 20   if ( clusterDD.le.0 .or. (clusterDD.eq.1.and.level.ne.-1) ) then
!20   if ( clusterDD.eq.0 .or. (clusterDD.eq.1.and.level.ne.-1) .or.
!    *    (clusterDD.eq.2.and.level.eq.-1) ) then
         if (clusterDD .lt. 0) then
            tmp_numprocs = abs(clusterDD)
            if(.true.)write(*,*)'level_refine: Decomposing into ',
     *                       tmp_numprocs
         else
            tmp_numprocs = numprocs
            if (nogridsonmaster.and.numprocs.gt.1)
     *                  tmp_numprocs = numprocs-1
         end if
         if(ltrace)write(*,*)'level_refine: calling level_clusterdd'
         call level_clusterDD(
     *             q(flagl),                  q(flagl+maxnumchildren),
     *             q(flagl+2*maxnumchildren), q(flagl+3*maxnumchildren),
     *             q(flagl+4*maxnumchildren), q(flagl+5*maxnumchildren),
     *             bbox_minx,bbox_maxx,
     *             bbox_miny,bbox_maxy,
     *             bbox_minz,bbox_maxz,
     *             numbox, tmp_numprocs, maxnumchildren )
c    *             numbox, numprocs, maxnumchildren )
c    *             numbox, 4, maxnumchildren )
      if (ltrace2.or.ltrace) then
         write(*,*) 'level_refine: Done clusterDDing, level: ',level
         write(*,*) 'level_refine: Done clusterDDing, numbox:',numbox
         do i = 1, numbox
            write(*,*) 'level_refine: i:',i
            write(*,*) 'level_refine: x:',bbox_minx(i),bbox_maxx(i)
            write(*,*) 'level_refine: y:',bbox_miny(i),bbox_maxy(i)
            write(*,*) 'level_refine: z:',bbox_minz(i),bbox_maxz(i)
            if (bbox_maxx(i) .lt. bbox_minx(i)) then
               call my_exit('level_refine: Bad box in x')
            end if
            if (bbox_maxy(i) .lt. bbox_miny(i)) then
               call my_exit('level_refine: Bad box in y')
            end if
            if (bbox_maxz(i) .lt. bbox_minz(i)) then
               call my_exit('level_refine: Bad box in z')
            end if
         end do
      end if
       end if

      !
      ! Figure out what processors to put grids on:
      !
c     if (numbox .gt. 0) call proc_pick_group(
c    *             bbox_minx,bbox_maxx,
c    *             bbox_miny,bbox_maxy,
c    *             bbox_minz,bbox_maxz,
c    *             child_procassignments,
c    *             minx,miny,minz, hl,
c    *             numbox, level)

      !
      ! Loop and create grids:
      !
      if (ltrace2)write(*,*) 'level_refine: Create # grids: ',numbox
      do i = 1, numbox
         !
         ! Properties of grid to be created:
         !
         nxc    = (bbox_maxx(i)-bbox_minx(i))*refine_factor+1
         nyc    = (bbox_maxy(i)-bbox_miny(i))*refine_factor+1
         nzc    = (bbox_maxz(i)-bbox_minz(i))*refine_factor+1
         minxc  = minx + (bbox_minx(i)-1) * hl
         minyc  = miny + (bbox_miny(i)-1) * hl
         minzc  = minz + (bbox_minz(i)-1) * hl
         if(ltrace)write(*,*)'level_refine: grid nx/y/zc: ',nxc,nyc,nzc
         !
         ! Look for existing, matching grid:
         !
         if (ltrace2)write(*,*)'Lookfor matching grid on level:',level+1
         gi_tmp = level_return_start(level+1)
 30      if ( grid_return_existence(gi_tmp) ) then
         if (ltrace2)write(*,*) ' Does gi_tmp match? ',gi_tmp
            if (       minxc .eq. gr_minx(gi_tmp)
     *           .and. minyc .eq. gr_miny(gi_tmp)
     *           .and. minzc .eq. gr_minz(gi_tmp)
     *           .and. nxc   .eq. gr_nx(gi_tmp)
     *           .and. nyc   .eq. gr_ny(gi_tmp)
     *           .and. nzc   .eq. gr_nz(gi_tmp) ) then
               if (ltrace)write(*,*) ' Matching grid found: ',gi_tmp
               call grid_mk_pending(gi_tmp)
               !call grid_mk_alive(gi_tmp)
               goto 60
            else if (ltrace2) then
               write(*,*) ' This grid is not a match:'
               write(*,*) '       ', minxc,gr_minx(gi_tmp)
               write(*,*) '       ', minyc,gr_miny(gi_tmp)
               write(*,*) '       ', minzc,gr_minz(gi_tmp)
               write(*,*) '       ', nxc,  gr_nx(gi_tmp)
               write(*,*) '       ', nyc,  gr_ny(gi_tmp)
               write(*,*) '       ', nzc,  gr_nz(gi_tmp)
            end if
            gi_tmp = grid_return_sibling(gi_tmp)
            goto 30
         end if
         !
         ! Assign grid number and pick a processor for it:
         !
         gi     = grid_return_number_norec()
c        proc   = child_procassignments(i)
         proc   = proc_pick(gi,nxc*nyc*nzc)
         if (proc .lt. 0) then
            write(*,*) 'level_refine: Not enough memory on level=',level
            write(*,*) 'level_refine: to create grid:            ',gi
            if (level.lt.0.or.clusterreadwrite.eq.1) then
               call my_exit('level_refine:Not enough memory for coarse')
            else if (nokeepalive) then
               call my_exit('level_refine:Not enough mem, nokeepalive')
            else if (numretries .le. MAXNUMRETRIES+1) then
               ! Change error threshold:
               if (numretries .eq. 0) then
                  write(*,*) '......-------------------------..........'
                  write(*,*) '......    WARNING!!!!          ..........'
                  write(*,*) '...... By all rights, this run ..........'
                  write(*,*) '...... should be dead, but I   ..........'
                  write(*,*) '...... will try to revive.     ..........'
                  write(*,*) '......-------------------------..........'
                  ethreshold_orig = ethreshold
               end if
               !call level_tree_dump()
               write(*,*)'level_refine:Increasing ethreshold',ethreshold
               if (numretries .eq.MAXNUMRETRIES+1) then
                  ! Last ditch effort, only points that *have* to be refined:
                  ethreshold = LARGENUMBER
               else
                  ethreshold = INCREASEFACTOR * ethreshold
               end if
               write(*,*)'level_refine:to:                  ',ethreshold
               ! Try again:
               numretries = numretries + 1
               write(*,*)'level_refine: numretries = ',numretries
               !
               ! Revert to pre-existing structure by:
               !    1) Resetting pending grids to dead
               !    2) Freeing all new grids (ones marked alive)
               !
               call level_mkpending_dead(level+1)
               call level_freeall_alive(level+1)
               !
               goto 33
!           else if (.not. level_return_existence(level+2)) then
!              write(*,*)'level_refine: INCREASEFACTOR =',INCREASEFACTOR
!              write(*,*)'level_refine: MAXNUMRETRIES  =',MAXNUMRETRIES 
!              write(*,*)'level_refine: No more retrying.'
!              goto 70
!           else 
!              ! If higher resolution levels exist, then we have to 
!              ! quit if we cannot create this level otherwise
!              ! the code will crash anyway:
!              call my_exit('level_refine:No way around more memory')
            else 
               ! I think we should just be able to back ourselves
               ! out of here retaining whatever resolution we had
               ! before we entered:
               call level_freeall_alive(level+1)
               call level_mkall_alive(level+1)
               write(*,*)'level_refine: No more retry. Try to backout'
               goto 70
            end if
         end if
         !
         if (ltrace) then
            write(*,*) 'level_refine: New grid: ',gi
            write(*,*) 'level_refine: on proc:  ',proc
            write(*,*) 'level_refine: x:',bbox_minx(i),bbox_maxx(i)
            write(*,*) 'level_refine: y:',bbox_miny(i),bbox_maxy(i)
            write(*,*) 'level_refine: z:',bbox_minz(i),bbox_maxz(i)
            write(*,*) 'level_refine: minx/y/zc: ',minxc,minyc,minzc
            write(*,*) 'level_refine: nx/y/zc: ',nxc,nyc,nzc
         end if
         if (ltrace) write(*,*) '   Sending all NEWGRID '
         !write(*,*) '   Sending all NEWGRID '
         call send_actionAll(NEWGRID,gi,proc,level+1,
     *                       gi,nxc,nyc,nzc,
     *                       bbox_minx(i),bbox_miny(i),bbox_minz(i))
         if (proc.ne.myid) then
            !
            ! Add grid to hierarchy:
            !
            if (ltrace) write(*,*) '   Adding grid locally:',gi,proc
            call grid_newa(gi,proc,level+1,
     *                     nxc,nyc,nzc,
     *                     minxc,minyc,minzc)
            ! Hear back from process which created grid:
            !    (because it may have run out of memory)
            if (ltrace) write(*,*) 'level_refine:Waiting for word',proc
            call recv_action(action,ggi,gowner,glevel,gparent,
     *                            gnx,gny,gnz,gax,gay,gaz,gfrom)
            if (ltrace2)write(*,*) 'level_refine:Word from:',gfrom
            if (action .eq. NEWGRID) then
               if (ltrace) write(*,*) '  That is what we expected'
            else if (action .eq. QUIT) then
               if (ltrace) write(*,*) '  Must have run out of memory'
               write(*,*)'level_refine: Proc ran out of memory: ',gfrom
               call my_exit('Processor ran out of memory creating grid')
            else
               write(*,*) 'level_refine: action ne NEWGRID: ',action
               write(*,*) 'level_refine: gfrom: ',gfrom
               call my_exit('level_refine: action ne NEWGRID')
            end if
         else
            !
            ! Create grid locally:
            !
            if (ltrace)write(*,*) ' Creating grid locally: ',gi
            call grid_createA(gi,myid,level+1,
     *                       nxc,nyc,nzc,
     *                       minxc,minyc,minzc     )
         end if
         if (ltrace2) call grid_dump_info(gi)
         !
         ! Init grid either by:
         !            -- if (t=0), 'exact' initial data
         !            -- interpolation from coarse grid parent
         !               & transfer values from an old fine grid

         ! Send message to initialize grid:

         if (ltrace) write(*,*) '   Sending all GRIDINIT ',gi
         !write(*,*) '   Sending all GRIDINIT ',gi
         call send_actionAll(GRIDINIT,gi,proc,level+1,
     *                       gi,nxc,nyc,nzc,
     *                       bbox_minx(i),bbox_miny(i),bbox_minz(i))
         call grid_initialize(gi)
         if(ltrace) write(*,*) 'level_refine: Finished grid init',gi
 60      continue
      end do

      !
      ! Any grids on the next level still marked as DEAD,
      ! should now be freed:
      !
 70   if(ltrace)write(*,*)'level_refine:Freeing dead grids lev=',level+1
      call level_freeall_dead(level+1)
      call level_mkpending_alive(level+1)

      if (ltrace) write(*,*) 'level_refine: Freeing allocated memory'
      call mem_dealloc( flagl, sizemem     )
      call mem_dealloc( flagk, nxl*nyl*nzl )
      call proc_submem( myid,  sizemem     )
      call proc_submem( myid,  nxl*nyl*nzl )

      if (ltrace) call level_tree_dump()
      if (ltrace)write(*,*)'level_refine:***Done refining level ',level

      !
      ! If we changed ethreshold because of memory, let's change it back:
      !
      if (numretries .gt. 0   .and.
     *       level_return_finest(Levelp,maxlev).eq.level
     *    .or. ethreshold .gt. LARGE_M_SMALL) then
         write(*,*)'level_refine: Resetting ethreshold to:',
     *                ethreshold_orig
         ethreshold = ethreshold_orig
         write(*,*)'level_refine: Resetting numretries...'
         numretries = 0
      end if

      return
      end      ! END: level_refine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_combine:                                                            cc
cc                 Take the error field from a given grid and combine         cc
cc                 it with the error field for the entire level in the        cc
cc                 appropriate way.                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_combine(level_error, grid_error,
     *                         mini,minj,mink,
     *                         nxl,nyl,nzl,
     *                         nxg,nyg,nzg)
      implicit none
      integer       nxl,nyl,nzl, nxg,nyg,nzg, mini,minj,mink
      real(kind=8)  level_error(nxl,nyl,nzl), grid_error(nxg,nyg,nzg)
      real(kind=8)  temp
      integer      i, j, k
      integer      il,jl,kl
      real(kind=8)  myl2norm3d, errornorm
      external      myl2norm3d

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_combine: Combining error:'
         write(*,*) 'level_combine: mini/j/k: ', mini,minj,mink
         write(*,*) 'level_combine: nx/y/zl:  ', nxl,nyl,nzl
         write(*,*) 'level_combine: nx/y/zg:  ', nxg,nyg,nzg
         write(*,*) 'level_combine: level_error:'
         call field_dump_infob(level_error,nxl,nyl,nzl,'level_error')
         call field_dump_infoB(grid_error, nxg,nyg,nzg,'griderror')
         call field_out3d(grid_error,0.d0,'errorZ',-1.d0,1.d0,
     *                     -1.d0,1.d0,-1.d0,1.d0,nxg,nyg,nzg,0)
      end if

      errornorm = myl2norm3d(grid_error,nxg,nyg,nzg)
      if (errornorm .gt. 1.0d9) then 
         write(*,*) 'level_combine:  * * *'
         write(*,*) 'level_combine:  Error/TRE is HUGE'
         write(*,*) 'level_combine: Refinement will almost certainly'
         write(*,*) 'level_combine: Fail because of this, and this'
         write(*,*) 'level_combine: would seem to indicate the run'
         write(*,*) 'level_combine: has blown up.'
         write(*,*) 'level_combine: You can look at the error in'
         write(*,*) 'level_combine: file 0errorZ.sdf'
         write(*,*) 'level_combine:  * * *'
      end if

      do k = 1, nzg
         kl = k + mink - 1
         do j = 1, nyg
            jl = j + minj - 1
            do i = 1, nxg
               il = i + mini - 1
               temp = grid_error(i,j,k)
               if (level_error(il,jl,kl) .lt. temp ) then
                  level_error(il,jl,kl)  =   temp
               end if
            end do
         end do
      end do

      if (ltrace) then
         !write(*,*) 'level_combine: level_error:'
         !call field_dump_info(level_error,nxl,nyl,nzl)
         !write(*,*) 'level_combine: grid_error:'
         !call field_dump_info(grid_error, nxg,nyg,nzg)
         call field_dump_infob(level_error,nxl,nyl,nzl,'level_error')
         call field_dump_infoB(grid_error, nxg,nyg,nzg,'griderror')
         write(*,*) 'level_combine: Done.'
      end if

      return
      end        ! END: level_combine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_remove_grid:                                                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_remove_grid(grid)
      implicit none
      integer       grid
      include      'mpif.h'
      include      'mpi_stuff.inc'
      include      'glob.inc'
      include      'grid.inc'
      include      'grid_methods.inc'
      integer       gi, gi_prev, gi_next
      integer       level

      logical      ltrace
      parameter (  ltrace = .false.)
      logical      ltrace2
      parameter (  ltrace2 = .false.)

      if (ltrace) then
         write(*,99) myid,'level_remove_grid: Remove grid: ',grid
c        call level_tree_dump()
      end if

      level    = grid_return_level(grid)

      gi_prev  = -1
      gi       = level_return_start(level)

  15  if (grid_return_existence(gi)) then
         gi_next = grid_return_sibling(gi)
         if (gi .eq. grid) then
            if (ltrace2) then
              write(*,99)myid,'level_remove_grid: Freeing  gi: ',gi
              write(*,99)myid,'level_remove_grid:     gi_prev: ',gi_prev
              write(*,99)myid,'level_remove_grid:     gi_next: ',gi_next
            end if
            if (gi_prev .eq. -1) then
               Levelp(level) = gi_next
            else
               call grid_set_sibling(gi_prev,gi_next)
            end if
            call grid_free(gi)
            call proc_subtractwork(grid_return_owner(gi),
     *                             grid_return_work(gi),
     *                             gr_nx(gi)*gr_ny(gi)*gr_nz(gi) )
            if (ltrace) call level_tree_dump()
            return
         else
            gi_prev = gi
         end if
         gi      = gi_next
         goto 15
      end if

      if (ltrace) then
         call level_tree_dump()
         write(*,99) myid, 'level_remove_grid: Grid not found'
      end if

  99  format('[',I3,'] ',A,3I5)

      return
      end      ! END: level_remove_grid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mkpending_dead:                                                     cc
cc                    All grids marked as PENDING change to DEAD.             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mkpending_dead( level)
      implicit none
      integer       level
      include      'grid_methods.inc'
      integer       gi

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_mkpending_dead: Make dead level: ',level
      end if

      gi       = level_return_start(level)

  15  if (grid_return_existence(gi)) then
         if (ltrace) write(*,*) 'leve_mkpending_dead: gi: ',gi
         if (grid_pending(gi)) call grid_mk_dead(gi)
         gi      = grid_return_sibling(gi)
         goto 15
      end if

      if (ltrace) then
         call level_tree_dump()
         write(*,*) 'level_mkpending_dead: dead level: ',level
      end if

      return
      end      ! END: level_mkpending_dead

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mkpending_alive:                                                    cc
cc                    All grids marked as PENDING change to ALIVE.            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mkpending_alive( level)
      implicit none
      integer       level
      include      'grid_methods.inc'
      integer       gi

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_mkpending_alive: Make alive level: ',level
      end if

      gi       = level_return_start(level)

  15  if (grid_return_existence(gi)) then
         if (ltrace) write(*,*) 'leve_mkpending_alive: gi: ',gi
         if (grid_pending(gi)) call grid_mk_alive(gi)
         gi      = grid_return_sibling(gi)
         goto 15
      end if

      if (ltrace) then
         call level_tree_dump()
         write(*,*) 'level_mkpending_alive: alive level: ',level
      end if

      return
      end      ! END: level_mkpending_alive

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mkall_alive:                                                        cc
cc                    Mark all grids on a given level as alive.               cc
cc                    This routine is to back out from refinement when        cc
cc                    we run out of memory to re-refine that level.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mkall_alive( level)
      implicit none
      integer       level
      !include      'glob.inc'
      include      'grid_methods.inc'
      integer       gi

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_mkall_alive: Make alive level: ',level
      end if

      gi       = level_return_start(level)

  15  if (grid_return_existence(gi)) then
         if (ltrace) write(*,*) 'leve_mkall_alive: Consider gi: ',gi
         call grid_mk_alive(gi)
         gi      = grid_return_sibling(gi)
         goto 15
      end if

      if (ltrace) then
         call level_tree_dump()
         write(*,*) 'level_mkall_alive: alive level: ',level
      end if

      return
      end      ! END: level_mkall_alive


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mkall_dead:                                                         cc
cc                    Mark all grids on a given level as dead.                cc
cc                    A dead grid has not been freeed yet, because we         cc
cc                    might either resue it (and hence make it alive), or     cc
cc                    use its data to initialize a new grid. Any grids        cc
cc                    still marked dead after refinement, will then be freed. cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mkall_dead( level)
      implicit none
      integer       level
      !include      'glob.inc'
      include      'grid_methods.inc'
      integer       gi

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_mkall_dead: Make dead level: ',level
      end if

      gi       = level_return_start(level)

  15  if (grid_return_existence(gi)) then
         if (ltrace) write(*,*) 'leve_mkall_dead: Consider gi: ',gi
         call grid_mk_dead(gi)
         gi      = grid_return_sibling(gi)
         goto 15
      end if

      if (ltrace) then
         call level_tree_dump()
         write(*,*) 'level_mkall_dead: Dead level: ',level
      end if

      return
      end      ! END: level_mkall_dead

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_freeall_dead:                                                       cc
cc                      Any grids in a given level still marked               cc
cc                      DEAD, will be removed from the hierarchy.             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_freeall_dead( level)
      implicit none
      integer       level
      include      'glob.inc'
      include      'grid.inc'
      include      'grid_methods.inc'
      include      'action.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer       gi, gi_prev, gi_next
      integer       action, ggi, gowner, glevel, gparent, gnx,gny,gnz
      integer       gax, gay, gaz, gfrom

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_freeall_dead: Free level: ',level
c        call level_tree_dump()
      end if

      gi_prev  = -1
      gi       = level_return_start(level)

  15  if (grid_return_existence(gi)) then
         if (ltrace) write(*,*) 'level_freeall_dead: Consider gi: ',gi
         gi_next = grid_return_sibling(gi)
         if (grid_dead(gi)) then
         !if (.not.grid_alive(gi)) then
            call send_actionAll(FREEGRID,gi,myid,0,0,0,0,0,0,0,0)
            if (ltrace) then
               write(*,*) 'level_freeall_dead: Freeing  gi: ',gi
               write(*,*) 'level_freeall_dead:     gi_prev: ',gi_prev
               write(*,*) 'level_freeall_dead:     gi_next: ',gi_next
            end if
            if (gi_prev .eq. -1) then
               Levelp(level) = gi_next
            else
               call grid_set_sibling(gi_prev,gi_next)
            end if
            call grid_free(gi)
            call proc_subtractwork(grid_return_owner(gi),
     *                             grid_return_work(gi),
     *                             gr_nx(gi)*gr_ny(gi)*gr_nz(gi) )
         else
            if (ltrace) write(*,*) 'level_freeall_dead: Alive!'
            gi_prev = gi
         end if
         gi      = gi_next
         goto 15
      end if

      if (ltrace) then
         call level_tree_dump()
         write(*,*) 'level_freeall_dead: Freed level: ',level
      end if

      return
      end      ! END: level_freeall_dead

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_freeall_alive:                                                      cc
cc                      Used to revert back to original grid structure.       cc
cc                      Any grids in a given level marked ALIVE               cc
cc                      are removed from teh hierarchy.                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_freeall_alive( level)
      implicit none
      integer       level
      include      'glob.inc'
      include      'grid.inc'
      include      'grid_methods.inc'
      include      'action.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer       gi, gi_prev, gi_next
      integer       action, ggi, gowner, glevel, gparent, gnx,gny,gnz
      integer       gax, gay, gaz, gfrom

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_freeall_alive: Free level: ',level
c        call level_tree_dump()
      end if

      gi_prev  = -1
      gi       = level_return_start(level)

  15  if (grid_return_existence(gi)) then
         if (ltrace) write(*,*) 'level_freeall_alive: Consider gi: ',gi
         gi_next = grid_return_sibling(gi)
         if (grid_alive(gi)) then
         !if (.not.grid_alive(gi)) then
            call send_actionAll(FREEGRID,gi,myid,0,0,0,0,0,0,0,0)
            if (ltrace) then
               write(*,*) 'level_freeall_alive: Freeing  gi: ',gi
               write(*,*) 'level_freeall_alive:     gi_prev: ',gi_prev
               write(*,*) 'level_freeall_alive:     gi_next: ',gi_next
            end if
            if (gi_prev .eq. -1) then
               Levelp(level) = gi_next
            else
               call grid_set_sibling(gi_prev,gi_next)
            end if
            call grid_free(gi)
            call proc_subtractwork(grid_return_owner(gi),
     *                             grid_return_work(gi),
     *                             gr_nx(gi)*gr_ny(gi)*gr_nz(gi) )
         else
            if (ltrace) write(*,*) 'level_freeall_alive: Alive!'
            gi_prev = gi
         end if
         gi      = gi_next
         goto 15
      end if

      if (ltrace) then
         call level_tree_dump()
         write(*,*) 'level_freeall_alive: Freed level: ',level
      end if

      return
      end      ! END: level_freeall_alive

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_makeflag:                                                           cc
cc                 Make flag array for level:                                 cc
cc                     (1) Initialize flag array to DISALLOW                  cc
cc                     (2) If error exceeds threshold, flag                   cc
cc                     (3) Everywhere a refined grid exists, flag             cc
cc                     (4) If any points are flagged, flag where masked       cc
cc                     (5) Add a buffer region of flagged points              cc
cc                     (6) Disallow ghostregion of level                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_makeflag( flag, error, level,
     *                           minx,miny,minz, h,nx,ny,nz,time)
      implicit none
      integer       nx,ny,nz, level
      real(kind=8)  minx,miny,minz, h, time
      real(kind=8)  flag(nx,ny,nz), error(nx,ny,nz)
      include      'glob.inc'
      include      'grid.inc'
      include      'grid_methods.inc'
      include      'mask.inc'
      include      'largesmall.inc'
      real(kind=8)  maxx,maxy,maxz
      integer       gi, numpoints
      integer       i, j, k, l
      integer       mini,maxi, minj,maxj, mink,maxk
      integer       bi,ei, bj,ej, bk,ek
      integer       resolution, izero,jzero, kzero
      logical       inmaskedregion, unresolved, cancellevel
      logical       childlevelexists
      real(kind=8)  r1_2, x,y,z
      logical       double_equal
      external      double_equal
      real(kind=8)  myl2norm3d, errornorm
      external      myl2norm3d

      !
      ! Minimum number of grid points in any direction for a masked region:
      !    (in units of number of grid points)
      !    (set to 0 to turn off this "feature")
      !
      integer       MINRESOLUTION
      parameter (   MINRESOLUTION = 13 )
      !parameter (   MINRESOLUTION = 8 )
      !
      ! Width of buffer region around a masked region to refine:
      !    (in units of number of grid points)
      !
      integer       MASKBUFFER
      parameter (   MASKBUFFER = 5 )

      ! Allow for refined regions to touch the outer boundary:
      !     (enabled for Mariana and her flux tubes)
      !     (now set if using periodicBC, but default is false)
      logical       ALLOWBOUND
      !parameter (   ALLOWBOUND = .false. )

      logical      ltrace
      parameter (  ltrace  = .false.)
      logical      ltrace2
      parameter (  ltrace2 = .false.)


      if (ltrace) write(*,*) 'level_makeflag: Initializing flag: ',level
      if(ltrace2) call field_dump_info(error,nx,ny,nz)
      errornorm = myl2norm3d(error,nx,ny,nz)
      if (errornorm .gt. 1.0d9) then
         write(*,*) 'level_makeflag:  * * *'
         write(*,*) 'level_makeflag: The norm of the TRE is HUGE!'
         write(*,*) 'level_makeflag: Refinement will almost certainly'
         write(*,*) 'level_makeflag: Fail because of this, and this'
         write(*,*) 'level_makeflag: would seem to indicate the run'
         write(*,*) 'level_makeflag: has blown up.'
         write(*,*) 'level_makeflag: You can look at the error in'
         write(*,*) 'level_makeflag: file 0errorB.sdf'
         write(*,*) 'level_makeflag:  * * *'
      end if
      if (ltrace2) then
          call field_out3d(error,0.d0,'errorB',-1.d0,1.d0,
     *                     -1.d0,1.d0,-1.d0,1.d0,nx,ny,nz,0)
      end if

      call load_scal3d(flag, FLAG_DISALLOW, nx,ny,nz )

      !
      ! Keep track of bounding box for FLAG_REFINE points:
      !
      numpoints = 0

      maxx = minx + (nx-1)*h
      maxy = miny + (ny-1)*h
      maxz = minz + (nz-1)*h
      if(ltrace)write(*,*)'level_makeflag:maxx/y/z:',maxx,maxy,maxz

      !
      ! If certain bad conditions, then we will not
      ! create this level
      !
      cancellevel      = .false.
      childlevelexists = .false.

      ! Find bounds of region requiring refinement:
      mini = nx
      minj = ny
      mink = nz
      maxi = 1
      maxj = 1
      maxk = 1

      if (ltrace) write(*,*) 'level_makeflag: Flagging where high error'
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (error(i,j,k) .ge. ethreshold-SMALLNUMBER) then
                  flag(i,j,k) = FLAG_REFINE
                  numpoints   = numpoints + 1
                  if(.false.)write(*,*)'level_makeflag: Flagged: ',
     *                                                i,j,k,error(i,j,k)
                  if (i.lt.mini) mini = i
                  if (j.lt.minj) minj = j
                  if (k.lt.mink) mink = k
                  if (i.gt.maxi) maxi = i
                  if (j.gt.maxj) maxj = j
                  if (k.gt.maxk) maxk = k
               else if (error(i,j,k) .ge. 0-SMALLNUMBER) then
                  flag(i,j,k) = FLAG_NOREFINE
               else if (ltrace2) then
                  if (NINT(error(i,j,k)).ne.-1) then
                     write(*,*) 'level_makeflag: Disallowed pt: ',
     *                                         error(i,j,k),i,j,k
                  end if
               end if
            end do
         end do
      end do
      if (ltrace2) then
         write(*,*) 'level_makeflag:   ethreshold: ',ethreshold
         write(*,*) 'level_makeflag:   numpoints flagged: ',numpoints
         write(*,*) 'level_makeflag: Bounds of region requiring refine:'
         write(*,*) 'level_makeflag:   mini/j/k:',mini,minj,mink
         write(*,*) 'level_makeflag:   maxi/j/k:',maxi,maxj,maxk
         write(*,*) 'level_makeflag:   nx/y/z:  ',nx,ny,nz
         call field_dump_info(error,nx,ny,nz,'error')
         call field_dump_infoB(error,nx,ny,nz,'error')
      end if
      bi = mini
      ei = maxi
      bj = minj
      ej = maxj
      bk = mink
      ek = maxk

      if (ltrace2) then
          call field_out3d(flag,time,'flagB',-1.d0,1.d0,
     *                     -1.d0,1.d0,-1.d0,1.d0,nx,ny,nz,0)
             write(*,*) 'level_makeflag:  lower x:',minx+(mini-1)*h
             write(*,*) 'level_makeflag:  upper x:',minx+(maxi-1)*h
      end if

      if (ltrace) write(*,*) 'level_makeflag: Flagging where refined'
      gi = level_return_start(level+2)
      if (ltrace2)write(*,*) 'level_makeflag: Starting w/ gi: ',gi
      !
 10   if (grid_return_existence(gi)) then
         childlevelexists = .true.
         if (ltrace) write(*,*) 'level_makeflag:   ...grid ',gi
         !
         mini = 1 + NINT( (gr_minx(gi)-minx)/h )
         minj = 1 + NINT( (gr_miny(gi)-miny)/h )
         mink = 1 + NINT( (gr_minz(gi)-minz)/h )
         ! Because resolutions are different, need to make sure
         ! that the full grid is "covered":
         if (minx+(mini-1)*h .gt. gr_minx(gi).and.mini.gt.1) mini=mini-1
         if (miny+(minj-1)*h .gt. gr_miny(gi).and.minj.gt.1) minj=minj-1
         if (minz+(mink-1)*h .gt. gr_minz(gi).and.mink.gt.1) mink=mink-1
         !
         maxi = 1 + NINT( (gr_maxx(gi)-minx)/h )
         maxj = 1 + NINT( (gr_maxy(gi)-miny)/h )
         maxk = 1 + NINT( (gr_maxz(gi)-minz)/h )
         if (minx+(maxi-1)*h .lt. gr_maxx(gi).and.maxi.lt.nx)maxi=maxi+1
         if (miny+(maxj-1)*h .lt. gr_maxy(gi).and.maxj.lt.ny)maxj=maxj+1
         if (minz+(maxk-1)*h .lt. gr_maxz(gi).and.maxk.lt.nz)maxk=maxk+1
         !
         if (ltrace) then
             write(*,*) 'level_makeflag:   ...looping'
             write(*,*) 'level_makeflag:   mini/j/k:',mini,minj,mink
             write(*,*) 'level_makeflag:   maxi/j/k:',maxi,maxj,maxk
             write(*,*) 'level_makeflag: h: ',h
             write(*,*) 'level_makeflag:    minx/y/z:',minx,
     *                  miny,minz
             write(*,*) 'level_makeflag: gr_minx/y/z:',gr_minx(gi),
     *                  gr_miny(gi),gr_minz(gi)
             write(*,*) 'level_makeflag: gr_maxx/y/z:',gr_maxx(gi),
     *                  gr_maxy(gi),gr_maxz(gi)
             write(*,*) 'level_makeflag:  lower x:',minx+(mini-1)*h
             write(*,*) 'level_makeflag:  upper x:',minx+(maxi-1)*h
          end if
         do k = mink, maxk
            do j = minj, maxj
               do i = mini, maxi
                  if( NINT(flag(i,j,k)).ne.FLAG_REFINE)
     *                       numpoints   = numpoints + 1
                  flag(i,j,k) = FLAG_REFINE
                  if (i.lt.bi) bi = i
                  if (i.gt.ei) ei = i
                  if (j.lt.bj) bj = j
                  if (j.gt.ej) ej = j
                  if (k.lt.bk) bk = k
                  if (k.gt.ek) ek = k
               end do
            end do
         end do
         !
         gi = grid_return_sibling(gi)
         goto 10
      end if
      if (ltrace2) then
          call field_out3d(flag,time,'flagC',-1.d0,1.d0,
     *                     -1.d0,1.d0,-1.d0,1.d0,nx,ny,nz,0)
          write(*,*) 'level_makeflag:   bi/j/k: ',bi,bj,bk
          write(*,*) 'level_makeflag:   ei/j/k: ',ei,ej,ek
      end if

      !
      ! Force refinement in masked region(s) if:
      !     ---any points have already been flagged
      !     ---if resolution of masked region too poor
      !
      if (ltrace)write(*,*) 'level_makeflag: Are there masked regions?'
      if (ltrace)write(*,*) 'level_makeflag: num_masks = ',num_masks
      if (num_masks.gt.0) then
         if (ltrace)write(*,*) 'level_makeflag: Masked region(s) exist'
         ! Only flag the masked region if the level has flagged points:
         !
         unresolved = .false.
         do l = 1, max_num_masks
            if (bh_true(l)) then
               resolution = min(
     *                     NINT(mask_coords(2,l)-mask_coords(1,l))/h,
     *                     NINT(mask_coords(4,l)-mask_coords(3,l))/h,
     *                     NINT(mask_coords(6,l)-mask_coords(5,l))/h
     *                          )
               unresolved = unresolved .or. resolution.lt.MINRESOLUTION
               if (numpoints.eq.0.and.unresolved) then
               !if (ltrace.or. (numpoints.eq.0.and.unresolved)) then
               write(*,*) 'level_makeflag: Forcing new level to resolve'
               write(*,*) 'level_makeflag: Hole #:       ',l
               !write(*,*) 'level_makeflag: numpoints:    ',numpoints
               write(*,*) 'level_makeflag: resolution:   ',resolution
               !write(*,*) 'level_makeflag: unresolved:   ',unresolved
               write(*,*) 'level_makeflag: MINRESOLUTION:',MINRESOLUTION
               end if
            end if
         end do
         !
         if ( numpoints .gt. 0 .or. unresolved) then
            if(ltrace)write(*,*)'level_makeFlagged pnts exist',numpoints
            do k = 1, nz
               z = minz + h * (k-1)
            do j = 1, ny
               y = miny + h * (j-1)
            do i = 1, nx
               x = minx + h * (i-1)
               inmaskedregion = .false.
               do l = 1, max_num_masks
                  if (bh_true(l)) then
                     !
                     ! Add buffer region around mask:
                     !
                     inmaskedregion = inmaskedregion .or.
     *          (      ( x .ge. (mask_coords(1,l)-MASKBUFFER*h) ) .and.
     *                 ( x .le. (mask_coords(2,l)+MASKBUFFER*h) ) .and.
     *                 ( y .ge. (mask_coords(3,l)-MASKBUFFER*h) ) .and.
     *                 ( y .le. (mask_coords(4,l)+MASKBUFFER*h) ) .and.
     *                 ( z .ge. (mask_coords(5,l)-MASKBUFFER*h) ) .and.
     *                 ( z .le. (mask_coords(6,l)+MASKBUFFER*h) ) )
                     !
                  end if
               end do
               if (inmaskedregion) then
                  ! This point should       be refined:
                  if( NINT(flag(i,j,k)).ne.FLAG_REFINE)
     *                       numpoints   = numpoints + 1
                  flag(i,j,k) = 1.d0 * FLAG_REFINE
                  !
                  ! If the masked region would occur too close
                  ! to the boundary of a new level then 
                  ! let us not even create the level:
                  !
                  if ( i.le.ghostwidth .or. i.gt.nx-ghostwidth .or.
     *                 j.le.ghostwidth .or. j.gt.ny-ghostwidth .or.
     *                 k.le.ghostwidth .or. k.gt.nz-ghostwidth ) then
                     cancellevel = .true.
                  end if
               end if
            end do
            end do
            end do
            if (ltrace)
     *      write(*,*)'level_makeflag: Total num flagged pts ',numpoints
         end if
      else
         if (ltrace)write(*,*) 'level_makeflag: No masked regions'
      end if

      if (ltrace) write(*,*) 'level_makeflag: Buffering flag array'
      if (ltrace2)write(*,*) 'level_makeflag:     buffer = ',buffer
      !
      ! Use "error" array as temporary storage for this routine:
      !
      if (ltrace2) then
          call field_out3d(flag,time,'flagpreb',-1.d0,1.d0,
     *                     -1.d0,1.d0,-1.d0,1.d0,nx,ny,nz,0)
      end if
      call mat_buffer( flag, error,  buffer, nx,ny,nz)
      if (ltrace2) then
          call field_out3d(flag,time,'flagpostb',-1.d0,1.d0,
     *                     -1.d0,1.d0,-1.d0,1.d0,nx,ny,nz,0)
      end if

      !
      ! Disallow clusters at the boundaries of the level
      !  NB:  this is where the level gets boundary data from its parent
      !  NB2: 
      !
      if (ltrace) write(*,*) 'level_makeflag: NOrefining boundaries'
      ALLOWBOUND = .false.
      if (periodicBC.gt.0) then
        if(ltrace)write(*,*)'level_makeflag: allowing for refinement'
        if(ltrace)write(*,*)'level_makeflag: at coarse grid boundaries'
        if(ltrace)write(*,*)'level_makeflag: because of periodic BCs'
        ALLOWBOUND = .true.
      end if

      ! x-boundaries
!     if(  mini.gt.0  .or.
!    *     (assume_symmetry.ne.1.and.assume_symmetry.ne.6) )then
      if(  (assume_symmetry.ne.1.and.assume_symmetry.ne.6) )then
         if (.not.double_equal(minx,minx0).or. .not.ALLOWBOUND) then
            if(ltrace)write(*,*)'level_makeflag: Disallowing minx bndry'
            do k = 1, nz
               do j = 1, ny
                  do i = 1, ghostwidth
                     flag(i,j,k) = FLAG_DISALLOW
                  end do
               end do
            end do
         else
            if(ltrace)write(*,*)'level_makeflag: Allowing minx bndry'
         end if
      else
         if (ltrace) then
         write(*,*)'level_makeflag: Allowing min x boundary'
         write(*,*)'level_makeflag: assume_symmetry=',assume_symmetry
         end if
         ! Allow flagging only two coarse grid points less than zero
         !        (NB: 4 fine grid points for extended boundary)
         izero   = -NINT(minx/h) + 1
         do k = 1, nz
            do j = 1, ny
               do i = 1,izero-3
                  flag(i,j,k) = FLAG_DISALLOW
               end do
            end do
         end do
      end if
      if (.not.double_equal(maxx,maxx0) .or. .not.ALLOWBOUND) then
         if(ltrace)write(*,*)'level_makeflag: Disallowing maxx bndry'
         do k = 1, nz
            do j = 1, ny
               do i = nx-ghostwidth+1, nx
                  flag(i,j,k) = FLAG_DISALLOW
               end do
            end do
         end do
      else
         if(ltrace)write(*,*)'level_makeflag: Allowing maxx bndry'
      end if

      !
      ! y-boundaries
!     if(  minj.gt.0  .or.
!    *     (assume_symmetry.ne.2.and.assume_symmetry.ne.6) )then
      if(  (assume_symmetry.ne.2.and.assume_symmetry.ne.6) )then
         if (.not.double_equal(miny,miny0).or. .not.ALLOWBOUND) then
            if(ltrace)write(*,*)'level_makeflag: Disallowing miny bndry'
            do k = 1, nz
               do i = 1, nx
                  do j = 1, ghostwidth
                     flag(i,j,k) = FLAG_DISALLOW
                  end do
               end do
            end do
         else
            if(ltrace)write(*,*)'level_makeflag: Allowing min y bndry'
         end if
      else
         if (ltrace) then
         write(*,*)'level_makeflag: Allowing min y boundary'
         write(*,*)'level_makeflag: assume_symmetry=',assume_symmetry
         end if
         ! Allow flagging only two coarse grid points less than zero
         !        (NB: 4 fine grid points for extended boundary)
         jzero   = -NINT(miny/h) + 1
         do k = 1, nz
            do i = 1, nx
               do j = 1,jzero-3
                  flag(i,j,k) = FLAG_DISALLOW
               end do
            end do
         end do
      end if
      !
      if (.not.double_equal(maxy,maxy0).or. .not.ALLOWBOUND) then
         if(ltrace)write(*,*)'level_makeflag: Disallowing outer y bndry'
         do k = 1, nz
            do i = 1, nx
               do j = ny-ghostwidth+1, ny
                  flag(i,j,k) = FLAG_DISALLOW
               end do
            end do
         end do
      else
         if(ltrace)write(*,*)'level_makeflag: Allowing outer y bndry'
      end if

      !
      ! z-boundaries
      !
      !    NB: If using reflection symmetry about z
      !        we do not want to disallow:
      !if (.not.(assume_symmetry.eq.3.and.minz .le. 0.d0) )then
!     if(  mink.gt.0  .or.
!    *     (assume_symmetry.ne.3.and.assume_symmetry.ne.6) )then
      if(  (assume_symmetry.ne.3.and.assume_symmetry.ne.6) )then
         if (.not.double_equal(minz,minz0).or. .not.ALLOWBOUND) then
            if(ltrace)write(*,*)'level_makeflag: Disallowing minz bndry'
            do j = 1, ny
               do i = 1, nx
                  do k = 1, ghostwidth
                     flag(i,j,k) = FLAG_DISALLOW
                  end do
               end do
            end do
         else
            if(ltrace)write(*,*)'level_makeflag: Allowing min.  z bndry'
         end if
      else
         if (ltrace) then
         write(*,*)'level_makeflag: Allowing min z boundary'
         write(*,*)'level_makeflag: assume_symmetry=',assume_symmetry
         write(*,*)'level_makeflag: minz           =',minz
         end if
         ! Allow flagging only two coarse grid points less than zero
         !        (NB: 4 fine grid points for extended boundary)
         kzero   = -NINT(minz/h) + 1
         do j = 1, ny
            do i = 1, nx
               do k = 1,kzero-3
                  flag(i,j,k) = FLAG_DISALLOW
               end do
            end do
         end do
      end if
      !
      if(ltrace2)write(*,*)'level_makeflag: maxz/0:',maxz,maxz0
      if (.not.double_equal(maxz,maxz0).or. .not.ALLOWBOUND) then
         if (ltrace) write(*,*)'level_makeflag: Disallowing max z bndry'
         do j = 1, ny
            do i = 1, nx
               do k = nz-ghostwidth+1, nz
                  flag(i,j,k) = FLAG_DISALLOW
               end do
            end do
         end do
      else
         if (ltrace) write(*,*)'level_makeflag: Allowing outer  z bndry'
      end if


      if (ltrace2) then
!         call field_out3d(flag,0.d0,'flagZ',-1.d0,1.d0,
!    *                     -1.d0,1.d0,-1.d0,1.d0,nx,ny,nz,0)
          call field_out3d(flag,time,'flagZ',minx,miny,minz,
     *                            maxx,maxy,maxz,nx,ny,nz,0)
      end if


      if (cancellevel) then
         if (childlevelexists) then
            write(*,*)'level_makeflag: --------------------------'
            write(*,*)'level_makeflag: level = ',level
            write(*,*)'level_makeflag: Cannot cancel level because'
            write(*,*)'level_makeflag: child level exists, but a masked'
            write(*,*)'level_makeflag: region is too close to boundary'
            write(*,*)'level_makeflag: --------------------------'
         else
            write(*,*) 'level_makeflag: Cancelling, level: ',level
            call load_scal1D(flag,1.d0*FLAG_NOREFINE,nx*ny*nz)
         end if
      end if

      if (ltrace) write(*,*) 'level_makeflag: Done on level ',level
      if (ltrace2) call field_dump_info(flag,nx,ny,nz)

      return
      end        ! END: level_makeflag


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_cluster:                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_cluster(flag,sigi,sigj,sigk, lapi,lapj,lapk,
     *                         asigi,asigj,asigk,alapi,alapj,alapk,
     *                         time, b_minx,b_maxx, b_miny,b_maxy,
     *                         b_minz,b_maxz,
     *                         minx,maxx,
     *                         miny,maxy,
     *                         minz,maxz,
     *                         numbox, nx,ny,nz )
      implicit     none
      integer      numbox, nx,ny,nz
      integer      b_minx(*), b_maxx(*),
     *             b_miny(*), b_maxy(*),
     *             b_minz(*), b_maxz(*)
      real(kind=8) flag(nx,ny,nz), sigi(nx),  sigj(ny),  sigk(nz),
     *                             lapi(nx),  lapj(ny),  lapk(nz),
     *                             asigi(nx), asigj(ny), asigk(nz),
     *                             alapi(nx), alapj(ny), alapk(nz),
     *                             minx, maxx, miny,maxy, minz,maxz,
     *                             time
      include     'glob.inc'
      include     'largesmall.inc'
      !
      !  sigi_min  ---- minimum value of Signature X within the window
      !  sigi_min_i---- index of "i" at which Signature X reaches minimum
      !  sigi_min_s---- span...
      !
      ! Could conserve memory here by using some temporary array for 
      ! the sig's and lap's (these are the so-called Laplacian's from
      ! the B&R article on their clustering).
      !
      real(kind=8)      sigi_min,       sigj_min,       sigk_min,
     *            lapi_max,       lapj_max,       lapk_max,
     *            alapi_max,      alapj_max,      alapk_max,
     *            lapmax,         efficiency,
     *            new1efficiency, new2efficiency, netefficiency
      real(kind=8)      tmpdp1, tmpdp2
      integer     i,              j,              k,   ii, l,
     *            sigi_min_i,     sigj_min_j,     sigk_min_k,
     *            sigi_min_s,     sigj_min_s,     sigk_min_s,
     *            lapi_max_i,     lapj_max_j,     lapk_max_k,
     *            alapi_max_i,    alapj_max_j,    alapk_max_k,
     *            lengthx,        lengthy,        lengthz,
     *            maxlength,      minlength,      numdisallowed,
     *            xcut_quality,   ycut_quality,   zcut_quality,
     *            tmpi
      integer     tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
      logical     backcutout
      logical     double_equal
      external    double_equal
      real(kind=8)getabsmax
      external    getabsmax
      !
      ! To store the original bbox before a cut:
      !
      integer      o_minx, o_maxx,
     *             o_miny, o_maxy,
     *             o_minz, o_maxz
       integer    cuttype,        cutaxis,        cutpos, already

       !
       ! Constants for code to decide on type of cut:
       !
       integer     NOCUT,         LAPCUT,       ZEROCUT,     DISCUT
       parameter ( NOCUT  = 0,    LAPCUT  = 1,  ZEROCUT = 2, DISCUT = 3)
       integer     XAXIS,         YAXIS,        ZAXIS
       parameter ( XAXIS  = 0,    YAXIS   = 1,  ZAXIS   = 2)
       ! Constant to remember which axis we have tried to cut:
       integer     NILL
       parameter ( NILL= -1)

      !
      ! Threshold by which an inflection cut must
      ! increase the net efficiency (between 0 and 1)
      ! or else be backed out
      !
      real(kind=8) ETHRESH
      !parameter  ( ETHRESH = 0.1d0 )
      !parameter  ( ETHRESH = 0.008d0 )
      parameter  ( ETHRESH = 0.02d0 )
      !
      ! Use cuts based on inflection points  as well as zero signatures ala B&R '91:
      !
      logical     inflection
      !
      ! If a box is found to contain any disallowed points (points
      ! not "covered" by its parent level), then we need to cut things up:
      !
      logical     anydisallowed

      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !
      ! Let user restrict clustering if needed:
      !
      if (clusterstyle .eq. 1) then
         if (ltrace) write(*,*) 'level_cluster: Turning off inflection'
         inflection = .false.
      else
         inflection = .true.
      end if

      !
      ! initialize first bounding box:
      !
      numbox    = 1
      i         = 1
      b_minx(i) = 1
      b_miny(i) = 1
      b_minz(i) = 1
      b_maxx(i) = nx
      b_maxy(i) = ny
      b_maxz(i) = nz
      ALREADY   = NILL

      if (ltrace) then
         write(*,*) '********************************'
         write(*,*) 'level_cluster: Clustering at time: ', time
         write(*,*) '                        nx,ny,nz:', nx,ny,nz
         write(*,*) '                      inflection:', inflection
         write(*,*) '                          mindim:', mindim
         write(*,*) '                   minefficiency:', minefficiency
         write(*,*) '                         ETHRESH:', ETHRESH
      end if

      !
      ! compute bounding box, and signatures:
      !
 15   continue
      if (ltrace) write(*,*) '   Shrinkwrapping box i = ',i
      call find_bbox( flag,
     *                sigi,           sigj,           sigk,
     *                b_minx(i),      b_miny(i),      b_minz(i),
     *                b_maxx(i),      b_maxy(i),      b_maxz(i),
     *                nx,             ny,             nz,
     *                efficiency )
      if ( efficiency .eq. 0 ) then
         if (ltrace) write(*,*) 'level_cluster: no flagged pts'
         numbox = 0
         return
      else if(double_equal(minx,minx0).and.double_equal(maxx,maxx0).and.
     *        double_equal(miny,miny0).and.double_equal(maxy,maxy0).and.
     *        double_equal(minz,minz0).and.double_equal(maxz,maxz0).and.
     *        (b_minx(i).eq.1).and. (b_maxx(i).eq.nx)  .and.
     *        (b_miny(i).eq.1).and. (b_maxy(i).eq.ny)  .and.
     *        (b_minz(i).eq.1).and. (b_maxz(i).eq.nz)  .and.
     *        abs(efficiency-1.d0) .lt. SMALLNUMBER              ) then
            if(ltrace)write(*,*)'level_cluster:CompleteCover',efficiency
            return
      else
         if (ltrace) write(*,*) 'level_cluster: efficiency=',efficiency
      end if
      !
      lengthx = b_maxx(i) - b_minx(i) + 1
      lengthy = b_maxy(i) - b_miny(i) + 1
      lengthz = b_maxz(i) - b_minz(i) + 1
      !
      ! compute asig (signature of the disallowed points, if any)
      !
      call compute_disallowed( flag,
     *                asigi,          asigj,          asigk,
     *                alapi,          alapj,          alapk,
     *                b_minx(i),      b_miny(i),      b_minz(i),
     *                b_maxx(i),      b_maxy(i),      b_maxz(i),
     *                nx,             ny,             nz,
     *                anydisallowed )
      if (anydisallowed) then
         if (ltrace) write(*,*) 'level_cluster: Disallowed points ',i
       call find_inflect(alapi(b_minx(i)),alapi_max,alapi_max_i,lengthx)
       call find_inflect(alapj(b_miny(i)),alapj_max,alapj_max_j,lengthy)
       call find_inflect(alapk(b_minz(i)),alapk_max,alapk_max_k,lengthz)
      else
         if (ltrace) write(*,*) 'level_cluster: Good, no disa. pts'
      end if

      !
      ! Avoid trying to cut where it would produce small grids:
      !     (0.1d0 is just an arbitrary, nonzero number)
      !
      call load_scal1d(sigi(b_minx(i)),0.1d0, mindim)
      call load_scal1d(sigj(b_miny(i)),0.1d0, mindim)
      call load_scal1d(sigk(b_minz(i)),0.1d0, mindim)
      call load_scal1d(sigi(b_minx(i)+lengthx-mindim),0.1d0, mindim)
      call load_scal1d(sigj(b_miny(i)+lengthy-mindim),0.1d0, mindim)
      call load_scal1d(sigk(b_minz(i)+lengthz-mindim),0.1d0, mindim)

      !
      ! Compute signatures and find minimums:
      !
      call find_min1d( sigi(b_minx(i)),
     *                 sigi_min, sigi_min_i, sigi_min_s, lengthx)
      call find_min1d( sigj(b_miny(i)),
     *                 sigj_min, sigj_min_j, sigj_min_s, lengthy)
      call find_min1d( sigk(b_minz(i)),
     *                 sigk_min, sigk_min_k, sigk_min_s, lengthz)

      !
      ! Compute laplacians (as per B&R algorithm):
      !
      call  compute_lap(lapi(b_minx(i)), sigi(b_minx(i)), lengthx)
      call  compute_lap(lapj(b_miny(i)), sigj(b_miny(i)), lengthy)
      call  compute_lap(lapk(b_minz(i)), sigk(b_minz(i)), lengthz)

      ! Zero out regions near beginning and ends of this box
      ! since we do not want to split there else we would get
      ! too small clusters:
      !    (needs to go one more point than with the signatures)
      !
      call load_scal1d(lapi(b_minx(i)),0.d0, mindim+1)
      call load_scal1d(lapj(b_miny(i)),0.d0, mindim+1)
      call load_scal1d(lapk(b_minz(i)),0.d0, mindim+1)
      call load_scal1d(lapi(b_minx(i)+lengthx-mindim-1),0.d0, mindim+1)
      call load_scal1d(lapj(b_miny(i)+lengthy-mindim-1),0.d0, mindim+1)
      call load_scal1d(lapk(b_minz(i)+lengthz-mindim-1),0.d0, mindim+1)

      !
      ! Look for maximum inflection points in laplacians:
      !
      call find_inflect(lapi(b_minx(i)), lapi_max, lapi_max_i, lengthx)
      call find_inflect(lapj(b_miny(i)), lapj_max, lapj_max_j, lengthy)
      call find_inflect(lapk(b_minz(i)), lapk_max, lapk_max_k, lengthz)

      if ( ltrace ) then
          write(*,*) '   ************'
          write(*,*) '   Box     i = ', i
          write(*,*) '       b_minx(i), mxx(i) = ',b_minx(i),b_maxx(i)
          write(*,*) '       b_miny(i), mxy(i) = ',b_miny(i),b_maxy(i)
          write(*,*) '       b_minz(i), mxz(i) = ',b_minz(i),b_maxz(i)
          write(*,*) '   mindim: ',mindim
          write(*,*) '  lengths: ',lengthx,lengthy,lengthz
          write(*,96)'  sigi_min, sigj_min, sigk_min: ',
     *                               sigi_min,sigj_min,sigk_min
          write(*,*) '  efficiency: ',efficiency
          write(*,*) 'lapi_max, lapi_max_i = ',lapi_max, lapi_max_i
          write(*,*) 'lapj_max, lapj_max_j = ',lapj_max, lapj_max_j
          write(*,*) 'lapk_max, lapk_max_k = ',lapk_max, lapk_max_k
            if (ltrace2) then
              if (.not.anydisallowed) then
              write(*,*)' * * *   i:   sigi(i):                lapi(i):'
               do j = b_minx(i), b_minx(i)+lengthx-1
                  write(*,*) j, sigi(j),lapi(j)
               end do
              write(*,*)' * * *   j:   sigj(j):                lapj(j):'
               do j = b_miny(i), b_miny(i)+lengthy-1
                  write(*,*) j, sigj(j),lapj(j)
               end do
              write(*,*)' * * *   k:   sigk(k):                lapk(k):'
               do j = b_minz(i), b_minz(i)+lengthz-1
                  write(*,*) j, sigk(j),lapk(j)
               end do
              else
                write(*,*)' Cluster contains disallowed points:',i
                write(*,*)'alapi_max,alapi_max_i=',alapi_max,alapi_max_i
                write(*,*)'alapj_max,alapj_max_j=',alapj_max,alapj_max_j
                write(*,*)'alapk_max,alapk_max_k=',alapk_max,alapk_max_k
              write(*,*)' * * *   j:   asigi(j):              alapi(j):'
                  do j = b_minx(i), b_minx(i)+lengthx-1
                     write(*,*) j, asigi(j),alapi(j)
                  end do
              write(*,*)' * * *   j:   asigj(j):              alapj(j):'
                  do j = b_miny(i), b_miny(i)+lengthy-1
                     write(*,*) j, asigj(j),alapj(j)
                  end do
              write(*,*)' * * *   j:   asigk(j):              alapk(j):'
                  do j = b_minz(i), b_minz(i)+lengthz-1
                     write(*,*) j, asigk(j),alapk(j)
                  end do
              end if
            end if
      end if

      !
      ! Cut bbox into simply connected pieces:
      !
      if (efficiency .lt. minefficiency .or. anydisallowed) then
         if (ltrace)write(*,*)'  Efficiency too low:',efficiency
         !
         ! Decide best place to split (if possible to split):
         !   (1) If can do a split at a zero of sig, do the
         !       "best" zero split first
         !   (2) If no zero splits possible, look for best
         !       possible inflection point split
         !   (3) If no split possible, simply move on to next
         !       box, if there is one
         !
         cuttype = NOCUT
         if (ltrace) write(*,*) '   Looking for zero cuts first'
         !
         ! Find the best possible zero cut:
         !
         xcut_quality = min( sigi_min_i, lengthx - sigi_min_i)
         ycut_quality = min( sigj_min_j, lengthy - sigj_min_j)
         zcut_quality = min( sigk_min_k, lengthz - sigk_min_k)
         if (ltrace) then
            write(*,*) '   Looking for zero cuts:'
            write(*,*) '      xcut_quality = ',xcut_quality
            write(*,*) '      ycut_quality = ',ycut_quality
            write(*,*) '      zcut_quality = ',zcut_quality
         end if
         if ( sigi_min .eq. 0 .and. xcut_quality .ge. mindim) then
            !
            ! Zero split possible only if:
            !                disconnected flagged regions
            !                "left"  region would have minimum dimension
            !                "right" region would have minimum dimension
            !
            if (ltrace) write(*,*) '    X zero cut possible'
            cuttype = ZEROCUT
            cutaxis = XAXIS
            cutpos  = sigi_min_i
         end if
         if ( sigj_min .eq. 0 .and. ycut_quality .ge. mindim) then
            if (ltrace) write(*,*) '    Y zero cut possible'
            !
            ! Pick cut with highest quality
            !
            if (sigi_min.eq.0 .and. xcut_quality.gt.ycut_quality) then
               if (ltrace) write(*,*) '    Go with X cut instead'
            else
               cuttype = ZEROCUT
               cutaxis = YAXIS
               cutpos  = sigj_min_j
            endif
         end if
         if ( sigk_min .eq. 0 .and. zcut_quality .ge. mindim) then
            if (ltrace) write(*,*) '    Z zero cut possible'
            !
            ! Pick cut with highest quality
            !
            if ( (sigi_min.eq.0.and.xcut_quality.gt.zcut_quality) .or.
     *           (sigj_min.eq.0.and.ycut_quality.gt.zcut_quality) ) then
               if (ltrace) write(*,*) '    Go with other cut instead'
            else
               cuttype = ZEROCUT
               cutaxis = ZAXIS
               cutpos  = sigk_min_k
            end if
         end if
         !
         ! Cuts based on disallowed points:
         !
         if ( cuttype .eq. NOCUT .and. anydisallowed ) then
            xcut_quality = min(alapi_max_i, lengthx -alapi_max_i)
            ycut_quality = min(alapj_max_j, lengthy -alapj_max_j)
            zcut_quality = min(alapk_max_k, lengthz -alapk_max_k)
            if (ltrace) then
               write(*,*) '   Looking for DISallowed cuts:'
               write(*,*) '      xcut_quality = ',xcut_quality
               write(*,*) '      ycut_quality = ',ycut_quality
               write(*,*) '      zcut_quality = ',zcut_quality
               write(*,*) '      alapi_max_i  = ',alapi_max_i
               write(*,*) '      alapj_max_j  = ',alapj_max_j
               write(*,*) '      alapk_max_k  = ',alapk_max_k
               write(*,*) '      lengthx      = ',lengthx
               write(*,*) '      lengthy      = ',lengthy
               write(*,*) '      lengthz      = ',lengthz
            end if
            !
            lapmax = 0.d0
            if (alapk_max .gt. 0.d0) then
                if (ltrace) write(*,*) '    zcut possible'
                cuttype = DISCUT
                cutaxis = ZAXIS
                cutpos  = alapk_max_k
                lapmax  = alapk_max
            end if
            if (alapj_max .gt. lapmax) then
                if (ltrace) write(*,*) '    ycut preferable'
                cuttype = DISCUT
                cutaxis = YAXIS
                cutpos  = alapj_max_j
                lapmax  = alapj_max
            end if
            if (alapi_max .gt. lapmax ) then
                if (ltrace) write(*,*) '    xcut most preferred'
                cuttype = DISCUT
                cutaxis = XAXIS
                cutpos  = alapi_max_i
                lapmax  = alapi_max
            end if
         end if
         if ( cuttype .eq. NOCUT .and. inflection ) then
            !
            ! Decide based on:
            ! (1) New boxes must have minimum dimensions
            ! (2) lapX_max must be non zero, else not an inflection point
            ! (3) Pick axis with the maximum value of the inflection point
            !
            !
            xcut_quality = min( lapi_max_i, lengthx - lapi_max_i)
            ycut_quality = min( lapj_max_j, lengthy - lapj_max_j)
            zcut_quality = min( lapk_max_k, lengthz - lapk_max_k)
            if (ltrace) then
               write(*,*) '   Looking for inflection pts:'
               write(*,*) '      xcut_quality = ',xcut_quality
               write(*,*) '      ycut_quality = ',ycut_quality
               write(*,*) '      zcut_quality = ',zcut_quality
            end if
            !
  500       lapmax = 0.d0
            if ( lapk_max .gt. 0.d0
     *          .and. ALREADY .ne. ZAXIS
     *          .and. zcut_quality .gt. mindim       ) then
                if (ltrace) write(*,*) '    zcut possible'
                cuttype = LAPCUT
                cutaxis = ZAXIS
                cutpos  = lapk_max_k
                lapmax  = lapk_max
            end if
            if ( lapj_max .gt. lapmax
     *          .and. ALREADY .ne. YAXIS
     *                    .and. ycut_quality .gt. mindim       ) then
                if (ltrace) write(*,*) '    ycut preferable'
                cuttype = LAPCUT
                cutaxis = YAXIS
                cutpos  = lapj_max_j
                lapmax  = lapj_max
            end if
            if ( lapi_max .gt. lapmax
     *          .and. ALREADY .ne. XAXIS
     *                    .and. xcut_quality .gt. mindim       ) then
                if (ltrace) write(*,*) '    xcut most preferred'
                cuttype = LAPCUT
                cutaxis = XAXIS
                cutpos  = lapi_max_i
                lapmax  = lapi_max
            end if
         end if
         !
         ! Make the cut:
         !
         if (cuttype .ne. NOCUT) then
            if (ltrace) then
               write(*,*) ' A cut has been chosen for box:',i
               write(*,*) '      cuttype = ',cuttype
               write(*,*) '      cutaxis = ',cutaxis
               write(*,*) '      cutpos  = ',cutpos
            end if
            !
            ! Store original box:
            !
            o_minx         = b_minx(i)
            o_miny         = b_miny(i)
            o_minz         = b_minz(i)
            o_maxx         = b_maxx(i)
            o_maxy         = b_maxy(i)
            o_maxz         = b_maxz(i)
            !
            ! Add new box, and copy bounding box:
            !
            numbox         = numbox + 1
            b_minx(numbox) = b_minx(i)
            b_miny(numbox) = b_miny(i)
            b_minz(numbox) = b_minz(i)
            b_maxx(numbox) = b_maxx(i)
            b_maxy(numbox) = b_maxy(i)
            b_maxz(numbox) = b_maxz(i)
            if      (cutaxis .eq. XAXIS) then
               b_maxx(i)      = b_minx(i) + cutpos - 1
               b_minx(numbox) = b_maxx(i) + 1
            else if (cutaxis .eq. YAXIS) then
               b_maxy(i)      = b_miny(i) + cutpos - 1
               b_miny(numbox) = b_maxy(i) + 1
            else if (cutaxis .eq. ZAXIS) then
               b_maxz(i)      = b_minz(i) + cutpos - 1
               b_minz(numbox) = b_maxz(i) + 1
            end if
            call find_bbox( flag,
     *                sigi,           sigj,           sigk,
     *                b_minx(i),      b_miny(i),      b_minz(i),
     *                b_maxx(i),      b_maxy(i),      b_maxz(i),
     *                nx,             ny,             nz,
     *                new1efficiency )
            call find_bbox( flag,
     *                sigi,           sigj,           sigk,
     *                b_minx(numbox), b_miny(numbox), b_minz(numbox),
     *                b_maxx(numbox), b_maxy(numbox), b_maxz(numbox),
     *                nx,             ny,             nz,
     *                new2efficiency )
            netefficiency = (
     *                (new1efficiency*(b_maxx(i)     -b_minx(i)+1)
     *                               *(b_maxy(i)     -b_miny(i)+1)
     *                               *(b_maxz(i)     -b_minz(i)+1))
     *              + (new2efficiency*(b_maxx(numbox)-b_minx(numbox)+1)
     *                               *(b_maxy(numbox)-b_miny(numbox)+1)
     *                               *(b_maxz(numbox)-b_minz(numbox)+1))
     *                       )/(
     *                                (b_maxx(i)     -b_minx(i)+1)
     *                               *(b_maxy(i)     -b_miny(i)+1)
     *                               *(b_maxz(i)     -b_minz(i)+1)
     *              +                 (b_maxx(numbox)-b_minx(numbox)+1)
     *                               *(b_maxy(numbox)-b_miny(numbox)+1)
     *                               *(b_maxz(numbox)-b_minz(numbox)+1))
            if (ltrace) then
              write(*,96)
     *         'old,new1,new2:',efficiency,new1efficiency,new2efficiency
             write(*,*)'netefficiency=',netefficiency
             write(*,*)'After cut: Box     i = ', i
             write(*,*)'  b_min/maxx(i)= ',b_minx(i),b_maxx(i),
     *                                     b_maxx(i)-b_minx(i)+1
             write(*,*)'  b_min/maxy(i)= ',b_miny(i),b_maxy(i),
     *                                     b_maxy(i)-b_miny(i)+1
             write(*,*)'  b_min/maxz(i)= ',b_minz(i),b_maxz(i),
     *                                     b_maxz(i)-b_minz(i)+1
             write(*,*)'After cut: Box     i = ', numbox
            write(*,*)'  b_min/maxx(i)= ',b_minx(numbox),b_maxx(numbox),
     *                                   b_maxx(numbox)-b_minx(numbox)+1
            write(*,*)'  b_min/maxy(i)= ',b_miny(numbox),b_maxy(numbox),
     *                                   b_maxy(numbox)-b_miny(numbox)+1
            write(*,*)'  b_min/maxz(i)= ',b_minz(numbox),b_maxz(numbox),
     *                                   b_maxz(numbox)-b_minz(numbox)+1
            end if
 96         format(A,3F15.5)
            !
            ! Back the cut out?
            !
            if (cuttype .eq. LAPCUT) then
               if (      (b_maxx(i)-b_minx(i) .le.   mindim) .or.
     *                   (b_maxy(i)-b_miny(i) .le.   mindim) .or.
     *                   (b_maxz(i)-b_minz(i) .le.   mindim) ) then
                  !
                  ! If new box is too small....
                  !
                  if (ltrace) write(*,*) ' Too small box i=',i
                  backcutout = .true.
               else if ( 
     *             (b_maxx(numbox)-b_minx(numbox) .le.  mindim).or.
     *             (b_maxy(numbox)-b_miny(numbox) .le.  mindim).or.
     *             (b_maxz(numbox)-b_minz(numbox) .le.  mindim))then
                  !
                  ! If new box is too small....
                  !
                  if (ltrace) write(*,*) ' Too small box i=',numbox
                  backcutout = .true.
               else if ( (netefficiency-efficiency).le.ETHRESH) then
                  !
                  ! If efficiency did not increase above threshold...
                  !
                  if (ltrace) write(*,*) ' No more efficient'
                  backcutout = .true.
               else
                  backcutout = .false.
               end if
            else
               backcutout = .false.
            end if

            if ( backcutout ) then
               !
               ! Restore original box:
               !
               if (ltrace) write(*,*) ' Backing cut out!'
               b_minx(i) = o_minx
               b_miny(i) = o_miny
               b_minz(i) = o_minz
               b_maxx(i) = o_maxx
               b_maxy(i) = o_maxy
               b_maxz(i) = o_maxz
               ! Remove added box:
               numbox    = numbox - 1
               !
               ! Try another cut (if we have not already):
               !
               if (cuttype .eq. LAPCUT .and. ALREADY.eq.NILL) then
                  !
                  ! Remember which axis we tried already:
                  !
                  if (ltrace) write(*,*) ' Try a different axis'
                  ALREADY = cutaxis
                  cuttype = NOCUT
                  goto 500
               else
                  ! Move on to consider the next box:
                  if (ltrace) write(*,*) ' Moving on to next box'
                  i         = i + 1
                  ALREADY   = NILL
               end if
            else
               if (ltrace) write(*,*) ' Cut is good'
               ALREADY = NILL
            end if
         else
            if (ltrace) write(*,*) ' No cut possible'
            i = i + 1
         end if
      else
         !
         ! Box of sufficient efficiency, look at next one.
         !
         if (ltrace)write(*,*)'Box has sufficient efficiency',efficiency
         i = i + 1
      end if
      !
      if (ltrace) write(*,*) 'Moving on: i, numbox = ',i,numbox
      if (numbox .eq. maxnumchildren) then
         write(*,*) ' We have reached maxnumchildren, no more boxes'
      else if (i.le.numbox) then 
         if (ltrace) write(*,*) ' '
         goto 15         ! more boxes to process
      end if

      !
      ! Clustering is now done, just need to do some extra stuff
      !
      if (ltrace .and. .true.) then
      !if (ltrace .and. .false.) then
          !
          ! For debugging purposes we want to output
          ! the flag array such that we can examine
          ! the clusters it produced:
          !    (1) Recompute signatures
          !    (2) For every flagged point in each
          !        box, replace with the grid number
          !    (3) Along the edges, place the signatures
          !
          tmp1 = 1
          tmp2 = 1
          tmp3 = 1
          tmp4 = nx
          tmp5 = ny
          tmp6 = nz
          call find_bbox( flag,
     *                sigi,           sigj,           sigk,
     *                tmp1, tmp2, tmp3, tmp4, tmp5, tmp6,
     *                nx,             ny,             nz,
     *                efficiency )
         call  compute_lap(lapi, sigi, nx)
         call  compute_lap(lapj, sigj, ny)
         call  compute_lap(lapk, sigk, nz)
          write(*,*) 'level_cluster: Outputting flag array'
          !
          ! (2) Show into which cluster flagged points fall
          !
          do l = 1, numbox
             write(*,*) 'level_cluster:   Marking box l=',l
             do k = b_minz(l), b_maxz(l)
             do j = b_miny(l), b_maxy(l)
             do i = b_minx(l), b_maxx(l)
                if (flag(i,j,k) .eq. FLAG_REFINE) then
                   flag(i,j,k) = 1.d0*l
                end if
             end do
             end do
             end do
          end do
          !
          ! (3) Add normalized signature information:
          !      ( This can cause problems when clusters
          !        are checked for negative values. )
          !
          if (.false.) then
          tmpdp1 = getabsmax(sigi,nx)/numbox
          tmpdp2 = getabsmax(lapi,nx)/numbox
          do k = 1, nz
          do i = 1, nx
             flag(i,1, k) = sigi(i)/tmpdp1
             flag(i,ny,k) = lapi(i)/tmpdp2
          end do
          end do
          tmpdp1 = getabsmax(sigj,ny)/numbox
          tmpdp2 = getabsmax(lapj,ny)/numbox
          do k = 1, nz
          do j = 1, ny
             flag(1, j,k) = sigj(j)/tmpdp1
             flag(nx,j,k) = lapj(j)/tmpdp2
          end do
          end do
          end if
          call field_out3d(flag,time,'flags',minx, maxx,
     *                     miny,maxy,minz,maxz,nx,ny,nz,0)
      end if

      !
      ! Add ghostwidth to clusters:
      !    (ghostwidth is defined on the fine grid, 
      !     so need to figure out how many coarse grid points)
      !
      tmpi = NINT( (1.d0*ghostwidth) / refine_factor )
      if (ltrace) write(*,*) '   Adding tmpi points: ',tmpi
      if (.false.) then
         ! This is replaced by routine below, but until
         ! that routine is more fully tested, I'm leaving
         ! this code in case anyone needs to return to it
         do i = 1, numbox
            if (ltrace) write(*,*) '   Growing box i=',i
            call grow_bbox(flag, FLAG_DISALLOW, tmpi,
     *               b_minx(i),b_maxx(i),
     *               b_miny(i),b_maxy(i),
     *               b_minz(i),b_maxz(i),
     *               nx, ny, nz )
         end do
      else
         if (ltrace) write(*,*) '   Growing box i=',i
         call grow_bboxall(flag, FLAG_DISALLOW, tmpi,
     *               b_minx,b_maxx,
     *               b_miny,b_maxy,
     *               b_minz,b_maxz,
     *               nx, ny, nz, numbox )
      end if


      !
      ! Done making boxes, now check them:
      !
      if (ltrace) write(*,*) '   Check bounds of boxes'
      do i = 1, numbox
         !
         if (ltrace) then
            write(*,*) '   Box     i = ', i
            write(*,*) '      b_minx(i), mxx(i) = ',b_minx(i),b_maxx(i),
     *                                             b_maxx(i)-b_minx(i)+1
            write(*,*) '      b_miny(i), mxy(i) = ',b_miny(i),b_maxy(i),
     *                                             b_maxy(i)-b_miny(i)+1
            write(*,*) '      b_minz(i), mxz(i) = ',b_minz(i),b_maxz(i),
     *                                             b_maxz(i)-b_minz(i)+1
         end if


         if ( ( (b_maxx(i)-b_minx(i)+1) .lt. mindim ) .or.
     *        ( (b_maxy(i)-b_miny(i)+1) .lt. mindim ) .or.
     *        ( (b_maxz(i)-b_minz(i)+1) .lt. mindim ) ) then
            write(*,*) 'level_cluster: Small box created'
            write(*,*) 'level_cluster: Presuming cut was correct,'
            write(*,*) 'level_cluster: but box shrank afterward.'
            write(*,*) '   nx,ny,nz  = ',nx,ny,nz
            write(*,*) '   mindim    = ',mindim
            write(*,*) '   numbox    = ',numbox
            write(*,*) '   Box     i = ', i
            write(*,*) '       b_minx(i), mxx(i) = ',b_minx(i),b_maxx(i)
            write(*,*) '       b_miny(i), mxy(i) = ',b_miny(i),b_maxy(i)
            write(*,*) '       b_minz(i), mxz(i) = ',b_minz(i),b_maxz(i)
            if (numbox .eq. 1) then
               write(*,*) 'level_cluster: Simply removing bad box:'
               write(*,*) 'level_cluster: Continuing...'
               numbox = 0
               return
            else
               write(*,*) 'level_cluster: More than one box=',numbox
               write(*,*) 'level_cluster: So cannot remove'
               write(*,*) '**** There very likely is a problem'
               write(*,*) '**** with the clusters produced, but'
               write(*,*) '**** I will not quit so maybe you '
               write(*,*) '**** can see what is wrong'
               write(*,*) ''
               !write(*,*) 'level_cluster: Stopping...'
               !call my_exit('level_cluster: Too small a cluster')
            end if
            write(*,*)'level_cluster: Outputting flag array: flagDEBUG'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: Please include this file',
     *                           ' for debugging'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            call field_out3d(flag,time,'flagDEBUG',minx,maxx,
     *                     miny,maxy,minz,maxz,nx,ny,nz,0)
         end if
         !
         ! Ensure sub grids aren't at interior grid boundaries:
         !
         if (  (b_minx(i).lt. 1) .or.  (b_maxx(i).gt. nx) .or.
     *         (b_miny(i).lt. 1) .or.  (b_maxy(i).gt. ny) .or.
     *         (b_minz(i).lt. 1) .or.  (b_maxz(i).gt. nz) ) then
            write(*,*) 'level_cluster: Problem:'
            write(*,*) '   nx,ny,nz  = ',nx,ny,nz
            write(*,*) '   mindim    = ',mindim
            write(*,*) '   Box     i = ', i
            write(*,*) '       b_minx(i), mxx(i) = ',b_minx(i),b_maxx(i)
            write(*,*) '       b_miny(i), mxy(i) = ',b_miny(i),b_maxy(i)
            write(*,*) '       b_minz(i), mxz(i) = ',b_minz(i),b_maxz(i)
            write(*,*) 'level_cluster: Bbox extends to/past parent bdy'
            write(*,*)'level_cluster: Outputting flag array: flagDEBUG'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: Please include this file',
     *                           ' for debugging'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            call field_out3d(flag,time,'flagDEBUG',minx,maxx,
     *                     miny,maxy,minz,maxz,nx,ny,nz,0)
            call my_exit('Problem grid_cluster Box extends to/past bdy')
         end if
         !
         ! Ensure bboxes are properly formed:
         !
         if ( (b_minx(i) .ge. b_maxx(i) ) .or.
     *        (b_miny(i) .ge. b_maxy(i) ) .or.
     *        (b_minz(i) .ge. b_maxz(i) )      ) then
            if (numbox .eq. 1) then
               numbox = 0
               if (ltrace) write(*,*) 'level_cluster: '
            end if
            write(*,*) 'level_cluster: Box is malformed'
            write(*,*) 'level_cluster: Box i: ', i
            write(*,*) '       b_minx(i), mxx(i) = ',b_minx(i),b_maxx(i)
            write(*,*) '       b_miny(i), mxy(i) = ',b_miny(i),b_maxy(i)
            write(*,*) '       b_minz(i), mxz(i) = ',b_minz(i),b_maxz(i)
            write(*,*) 'level_cluster: Quitting...'
            write(*,*)'level_cluster: Outputting flag array: flagDEBUG'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: Please include this file',
     *                           ' for debugging'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            call field_out3d(flag,time,'flagDEBUG',minx,maxx,
     *                     miny,maxy,minz,maxz,nx,ny,nz,0)
            call my_exit('level_cluster: Malformed box')
         end if
         !
         ! Check that no disallowed points (points that would NOT be
         ! fully nested in the parent level) are included:
         !
         numdisallowed = 0
         do k  = b_minz(i), b_maxz(i)
         do j  = b_miny(i), b_maxy(i)
         do ii = b_minx(i), b_maxx(i)
            if (double_equal(flag(ii,j,k),FLAG_DISALLOW)) then
               numdisallowed = numdisallowed + 1
               write(*,*) 'level_cluster: at point ii,j,k:',ii,j,k
            end if
         end do
         end do
         end do
         if (numdisallowed.gt.0) then
            write(*,*) 'level_cluster: Disallowed point found'
            write(*,*) 'level_cluster: Cluster: ',i
            write(*,*) '       b_minx(i), mxx(i) = ',b_minx(i),b_maxx(i)
            write(*,*) '       b_miny(i), mxy(i) = ',b_miny(i),b_maxy(i)
            write(*,*) '       b_minz(i), mxz(i) = ',b_minz(i),b_maxz(i)
            write(*,*) '   nx,ny,nz  = ',nx,ny,nz
            write(*,*) '   numbox    = ',numbox
            write(*,*) 'level_cluster: numdisallowed: ',numdisallowed
            write(*,*) 'level_cluster: outputting flagdis'
             !do l = numbox, 1, -1
             do l = 1, numbox
             write(*,*) 'level_cluster:   Marking box l=',l
            write(*,*) '       b_minx(i), mxx(i) = ',b_minx(l),b_maxx(l)
            write(*,*) '       b_miny(i), mxy(i) = ',b_miny(l),b_maxy(l)
            write(*,*) '       b_minz(i), mxz(i) = ',b_minz(l),b_maxz(l)
             do k  = b_minz(l), b_maxz(l)
             do j  = b_miny(l), b_maxy(l)
             do ii = b_minx(l), b_maxx(l)
                if (flag(ii,j,k) .eq. FLAG_REFINE) then
                   flag(ii,j,k) = 1.d0*l
                end if
             end do
             end do
             end do
            end do
            write(*,*)'level_cluster: Outputting flag array: flagdis'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: Please include this file',
     *                           ' for debugging'
            write(*,*)'level_cluster: *******'
            write(*,*)'level_cluster: *******'
            call field_out3d(flag,time,'flagdis',minx, maxx,
     *                     miny,maxy,minz,maxz,nx,ny,nz,0)
c        call level_output_flag(   flag, 0.d0,
c    *                    b_minx,b_maxx,
c    *                    b_miny,b_maxy,
c    *                    b_minz,b_maxz,
c    *                    'flag',nx,ny,nz,numbox)
            call my_exit('level_cluster: Disallowed pt clustered')
         end if
      end do

 800  if (ltrace) then
         write(*,*) 'level_cluster: DONE Clustering:'
         do i = 1, numbox
           write(*,*) '   Box     i = ', i
           write(*,*) '       b_minx(i), mxx(i) = ',b_minx(i),b_maxx(i),
     *                                             b_maxx(i)-b_minx(i)+1
           write(*,*) '       b_miny(i), mxy(i) = ',b_miny(i),b_maxy(i),
     *                                             b_maxy(i)-b_miny(i)+1
           write(*,*) '       b_minz(i), mxz(i) = ',b_minz(i),b_maxz(i),
     *                                             b_maxz(i)-b_minz(i)+1
         end do
         write(*,*) '********************************'
      end if

      !write(*,*) 'h'
      return
      end      ! END: level_cluster


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_clusterDD:                                                          cc
cc                  Take as input a the set of subgrids to create             cc
cc                  And output that same set but domain decomposed            cc
cc                  according to the number of processors.                    cc
cc                                                                            cc
cc                  NB: to help ensure that the chopped up grids              cc
cc                      get placed onto all different processors,             cc
cc                      the list should be kept in order instead of           cc
cc                      putting the new blocks at the end of the list.        cc
cc                      Hence the need for temp storage. Reals are used       cc
cc                      because such space is readily available in            cc
cc                      level_refine().                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_clusterDD(tmp_mini, tmp_maxi, tmp_minj, tmp_maxj,
     *                           tmp_mink, tmp_maxk, b_mini,   b_maxi,
     *                           b_minj,   b_maxj,   b_mink,   b_maxk,
     *                           numbox, nump, maxnum)
      implicit     none
      integer      numbox, nump, maxnum
      real(kind=8) tmp_mini(maxnum), tmp_maxi(maxnum),
     *             tmp_minj(maxnum), tmp_maxj(maxnum),
     *             tmp_mink(maxnum), tmp_maxk(maxnum)
      integer      b_mini(maxnum),   b_maxi(maxnum),
     *             b_minj(maxnum),   b_maxj(maxnum),
     *             b_mink(maxnum),   b_maxk(maxnum)
      include      'glob.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer      ii,   i, j, k, dims(3), gwc
      integer      lengthi, lengthj, lengthk
      integer      numx, numy, numz
      ! Offsets to overlap decomposed grids:
      integer      minus, plus
      ! Minimum dimension into which to cut:
      integer      minsize
      integer      boxcount, fixcount
      integer      fixminx, fixminy, fixminz, fixmaxx, fixmaxy, fixmaxz
      logical      fixshow

      logical      ltrace
      parameter (  ltrace  = .false. )
      
      !
      ! Each cluster which we are trying to split among
      ! the procs should have at least ghostwidth+mindim
      ! points in each direction. If we split, then each
      !
      minsize = ghostwidth / refine_factor + mindim

      ! The ghostwidth as it appears on the parent level:
      gwc  = ghostwidth / refine_factor + 1

      if (ltrace) then
         write(*,*) '********************************'
         write(*,*) 'level_clusterDD:  numbox     = ',numbox
         write(*,*) 'level_clusterDD:  mindim     = ',mindim
         write(*,*) 'level_clusterDD:  ghostwidth = ',ghostwidth
         write(*,*) 'level_clusterDD:  gwc        = ',gwc       
         write(*,*) 'level_clusterDD:  minsize    = ',minsize
      end if

      !            
      ! Copy entries to temp storage:
      !            
      do ii = 1, numbox
         tmp_mini(ii) = b_mini(ii)
         tmp_maxi(ii) = b_maxi(ii)
         tmp_minj(ii) = b_minj(ii)
         tmp_maxj(ii) = b_maxj(ii)
         tmp_mink(ii) = b_mink(ii)
         tmp_maxk(ii) = b_maxk(ii)
      end do

      if (nump .gt. 64) then
         dims(1) = NINT(nump**(0.3333))
         dims(3) = nump/dims(1)**2
         dims(2) = nump/dims(1)/dims(3)
      else if (nump .ne. 64) then
         dims(1) = 0      !  x: LEFT   and RIGHT
         dims(2) = 0      !  y: BACK   and FRONT
         dims(3) = 0      !  z: BOTTOM and TOP
         call MPI_DIMS_CREATE( nump, 3, dims, ierr )
      else
         ! For some strange reason above function yields 8x4x2
         ! instead of 4x4x4, so bypass
         dims(1) = 4
         dims(2) = 4
         dims(3) = 4
      end if
      if (ltrace) write(*,*) 'level_clusterDD: dims1/2/3:',
     *                             dims(1),dims(2),dims(3)


      !            
      ! For each box, break up into nump pieces:
      !            
      boxcount = 0
      do ii = 1, numbox
         ! Keep track if clustering needs to be "fixed":
         fixminx = 0
         fixminy = 0
         fixminz = 0
         fixmaxx = 0
         fixmaxy = 0
         fixmaxz = 0
         fixcount= boxcount
         lengthi = NINT( (tmp_maxi(ii) - tmp_mini(ii) ) / dims(1) )
         lengthj = NINT( (tmp_maxj(ii) - tmp_minj(ii) ) / dims(2) )
         lengthk = NINT( (tmp_maxk(ii) - tmp_mink(ii) ) / dims(3) )
         if (ltrace) then
            write(*,*) 'level_clusterDD:  Splitting box: ',ii
            write(*,*) 'level_clusterDD:  i: ',tmp_mini(ii),tmp_maxi(ii)
            write(*,*) 'level_clusterDD:  j: ',tmp_minj(ii),tmp_maxj(ii)
            write(*,*) 'level_clusterDD:  k: ',tmp_mink(ii),tmp_maxk(ii)
            write(*,*) 'level_clusterDD:    lengthi/j/k: ',
     *                  lengthi, lengthj, lengthk
         end if
         !
         ! Decide how we should chop up:
         !    (1) Chop as per MPI_DIMS_Create() says if possible
         !        (ie if that would produce sufficient size boxes)
         !    (2) Else chop into "minsize" pieces
         !
         if (lengthi .lt. minsize .and. dims(1).gt.1) then
         !if (lengthi .lt. mindim .and. dims(1).gt.1) then
            numx = (tmp_maxi(ii)-tmp_mini(ii)) / minsize
            if (numx .le. 0) numx = 1
            lengthi = NINT( (tmp_maxi(ii) - tmp_mini(ii) ) / numx )
         else
            numx = dims(1)
         end if
         if (lengthj .lt. minsize .and. dims(2).gt.1) then
            numy = (tmp_maxj(ii)-tmp_minj(ii)) / minsize
            if (numy .le. 0) numy = 1
            lengthj = NINT( (tmp_maxj(ii) - tmp_minj(ii) ) / numy )
         else
            numy = dims(2)
         end if
         if (lengthk .lt. minsize .and. dims(3).gt.1) then
            numz = (tmp_maxk(ii)-tmp_mink(ii)) / minsize
            if (numz .le. 0) numz = 1
            if (numx*numy*numz.gt.nump) numz = nump / numx / numy
            lengthk = NINT( (tmp_maxk(ii) - tmp_mink(ii) ) / numz )
         else
            numz = dims(3)
         end if
         if (ltrace) write(*,*) '   split numbers:',
     *                       numx,numy,numz,numx*numy*numz,nump
         !
         ! If the boxes need splitting and we have enough
         ! children to accomodate:
         !
         if ( (numx.gt.1.or.numy.gt.1.or.numz.gt.1) 
     *                       .and. 
     *         (maxnum .ge. (boxcount+numx*numy*numz)) ) then
         !if (numx.gt.1.or.numy.gt.1.or.numz.gt.1) then
            !
            ! Split up box "ii":
            !
            do k = 1, numz
            do j = 1, numy
            do i = 1, numx
               boxcount         = boxcount + 1
               b_mini(boxcount) = tmp_mini(ii) + (i-1) * lengthi
               b_maxi(boxcount) = tmp_mini(ii) +   i   * lengthi
               b_minj(boxcount) = tmp_minj(ii) + (j-1) * lengthj
               b_maxj(boxcount) = tmp_minj(ii) +   j   * lengthj
               b_mink(boxcount) = tmp_mink(ii) + (k-1) * lengthk
               b_maxk(boxcount) = tmp_mink(ii) +   k   * lengthk
               !
               ! Force overlap:
               !    NB: The amount of overlap among the grids
               !    depends on the stencils used and how
               !    the parameter bwidth is set in util.f.
               !    Also, keep in mind that we're working
               !    in the coarse grid space in which the
               !    grids are defined but the grids to be
               !    completed will have greater resolution.
               !    Assuming a minimum refinement factor of 2,
               !    then we need to overlap with 4 of these
               !    coarse points so that the grids themselves
               !    share 2*bwidth points. 
               !
               plus  = bound_width / 2
               if (mod(bound_width,2) .eq. 0) then
                  minus = plus
               else
                  minus = plus + 1
               end if
               b_mini(boxcount) = b_mini(boxcount) - minus
               b_maxi(boxcount) = b_maxi(boxcount) + plus
               b_minj(boxcount) = b_minj(boxcount) - minus
               b_maxj(boxcount) = b_maxj(boxcount) + plus
               b_mink(boxcount) = b_mink(boxcount) - minus
               b_maxk(boxcount) = b_maxk(boxcount) + plus
               !
               ! In case, dimensions don't divide evenly
               ! and to overrule the overlap at the bounds,
               ! make sure we cover the whole box:
               !
               if (i.eq. 1    ) b_mini(boxcount) = tmp_mini(ii)
               if (j.eq. 1    ) b_minj(boxcount) = tmp_minj(ii)
               if (k.eq. 1    ) b_mink(boxcount) = tmp_mink(ii)
               if (i.eq. numx ) b_maxi(boxcount) = tmp_maxi(ii)
               if (j.eq. numy ) b_maxj(boxcount) = tmp_maxj(ii)
               if (k.eq. numz ) b_maxk(boxcount) = tmp_maxk(ii)
               !
               ! Must ensure that grids do not terminate within
               ! the ghostregion for the level...that is, each
               ! grid that does not touch the amr boundary on
               ! one of its faces, should not end *close* to
               ! the amr boundary. If one does, one gets points
               ! that are both injected to the parent, and
               ! those same points on the parent are used to
               ! reset the boundary.
               !
               if(ltrace) then
                  write(*,*) '---boxcount= ',boxcount
                  write(*,*) '    i: ',b_mini(boxcount),b_maxi(boxcount)
                  write(*,*) '    j: ',b_minj(boxcount),b_maxj(boxcount)
                  write(*,*) '    k: ',b_mink(boxcount),b_maxk(boxcount)
                  write(*,*) '  . . . '
               end if
               if (b_mini(boxcount).ne.tmp_mini(ii) .and.
     *             b_mini(boxcount).le.tmp_mini(ii)+gwc-1 ) then
                  if(ltrace)write(*,*)'level_clusterDD: fixing small i',
     *                   b_mini(boxcount),tmp_mini(ii)+gwc
                  if ((tmp_mini(ii)+gwc-b_mini(boxcount)) .gt. fixminx)
     *               fixminx = tmp_mini(ii)+gwc-b_mini(boxcount)
                  b_mini(boxcount) = tmp_mini(ii)+gwc
               end if
               if (b_maxi(boxcount).ne.tmp_maxi(ii) .and.
     *             b_maxi(boxcount).ge.tmp_maxi(ii)-gwc+1 ) then
                  if(ltrace)write(*,*)'level_clusterDD: fixing large i',
     *                   b_maxi(boxcount),tmp_maxi(ii)-gwc
                  if ((b_maxi(boxcount)- (tmp_maxi(ii)-gwc).gt.fixmaxx))
     *               fixmaxx = b_maxi(boxcount)- (tmp_maxi(ii)-gwc)
                  b_maxi(boxcount) = tmp_maxi(ii)-gwc
               end if
               if (b_minj(boxcount).ne.tmp_minj(ii) .and.
     *             b_minj(boxcount).le.tmp_minj(ii)+gwc-1 ) then
                  if(ltrace)write(*,*)'level_clusterDD: fixing small j',
     *                   b_minj(boxcount),tmp_minj(ii)+gwc
                  if ( (tmp_minj(ii)+gwc-b_minj(boxcount)) .gt. fixminy)
     *               fixminy = tmp_minj(ii)+gwc-b_minj(boxcount)
                  b_minj(boxcount) = tmp_minj(ii)+gwc
               end if
               if (b_maxj(boxcount).ne.tmp_maxj(ii) .and.
     *             b_maxj(boxcount).ge.tmp_maxj(ii)-gwc+1 ) then
                  if(ltrace)write(*,*)'level_clusterDD: fixing large j',
     *                   b_maxj(boxcount),tmp_maxj(ii)-gwc
                  if ((b_maxj(boxcount)- (tmp_maxj(ii)-gwc).gt.fixmaxy))
     *               fixmaxy = b_maxj(boxcount)- (tmp_maxj(ii)-gwc)
                  b_maxj(boxcount) = tmp_maxj(ii)-gwc
               end if
               if (b_mink(boxcount).ne.tmp_mink(ii) .and.
     *             b_mink(boxcount).le.tmp_mink(ii)+gwc-1 ) then
                  if(ltrace)write(*,*)'level_clusterDD: fixing small k',
     *                   b_mink(boxcount),tmp_mink(ii)+gwc
                  if ( (tmp_mink(ii)+gwc-b_mink(boxcount)) .gt. fixminz)
     *               fixminz = tmp_mink(ii)+gwc-b_mink(boxcount)
                  b_mink(boxcount) = tmp_mink(ii)+gwc
               end if
               if (b_maxk(boxcount).ne.tmp_maxk(ii) .and.
     *             b_maxk(boxcount).ge.tmp_maxk(ii)-gwc+1 ) then
                  if(ltrace)write(*,*)'level_clusterDD: fixing large k',
     *                   b_maxk(boxcount),tmp_maxk(ii)-gwc
                  if ((b_maxk(boxcount)- (tmp_maxk(ii)-gwc).gt.fixmaxz))
     *               fixmaxz = b_maxk(boxcount)- (tmp_maxk(ii)-gwc)
                  b_maxk(boxcount) = tmp_maxk(ii)-gwc
               end if
               !
               if (ltrace) then
                  write(*,*) 'bound_width = ',bound_width
                  write(*,*) 'plus        = ',plus
                  write(*,*) 'minus       = ',minus
                  write(*,*) 'gwc         = ',gwc
                  write(*,*) '   boxcount= ',boxcount
                  write(*,*) '    i: ',b_mini(boxcount),b_maxi(boxcount)
                  write(*,*) '    j: ',b_minj(boxcount),b_maxj(boxcount)
                  write(*,*) '    k: ',b_mink(boxcount),b_maxk(boxcount)
                  write(*,*) '-----------------------------'
               end if
            end do
            end do
            end do
            if (ltrace) then
            !if (ltrace.or. .true.) then
               write(*,*)'level_clusterDD:fixminx/y/z:',fixminx,fixminy,
     *                        fixminz
               write(*,*)'level_clusterDD:fixmaxx/y/z:',fixmaxx,fixmaxy,
     *                        fixmaxz
            end if
            !
            ! Fix the boxes on the border:
            !
            do k = fixcount+1,boxcount
               fixshow = .false.
               if (ltrace) then
                  write(*,*) 'level_clusterDD: Fixing box: ',k
                  write(*,*) '    i: ',b_mini(k),b_maxi(k)
                  write(*,*) '    j: ',b_minj(k),b_maxj(k)
                  write(*,*) '    k: ',b_mink(k),b_maxk(k)
               end if
               if (fixminx.gt.0 .and. b_mini(k).eq.tmp_mini(ii)) then
                  b_maxi(k) = b_maxi(k) + fixminx
                  if (b_maxi(k).gt.tmp_maxi(ii)) b_maxi(k)=tmp_maxi(ii)
                  if(ltrace)write(*,*)'level_clusterDD: fixed minx'
                  fixshow = .true.
               end if
               if (fixmaxx.gt.0 .and. b_maxi(k).eq.tmp_maxi(ii)) then
                  b_mini(k) = b_mini(k) - fixmaxx
                  if (b_mini(k).lt.tmp_mini(ii)) b_mini(k)=tmp_mini(ii)
                  if(ltrace)write(*,*)'level_clusterDD: fixed maxx'
                  fixshow = .true.
               end if
               if (fixminy.gt.0 .and. b_minj(k).eq.tmp_minj(ii)) then
                  b_maxj(k) = b_maxj(k) + fixminy
                  if (b_maxj(k).gt.tmp_maxj(ii)) b_maxj(k)=tmp_maxj(ii)
                  if(ltrace)write(*,*)'level_clusterDD: fixed miny'
                  fixshow = .true.
               end if
               if (fixmaxy.gt.0 .and. b_maxj(k).eq.tmp_maxj(ii)) then
                  b_minj(k) = b_minj(k) - fixmaxy
                  if (b_minj(k).lt.tmp_minj(ii)) b_minj(k)=tmp_minj(ii)
                  if(ltrace)write(*,*)'level_clusterDD: fixed maxy'
                  fixshow = .true.
               end if
               if (fixminz.gt.0 .and. b_mink(k).eq.tmp_mink(ii)) then
                  b_maxk(k) = b_maxk(k) + fixminz
                  if (b_maxk(k).gt.tmp_maxk(ii)) b_maxk(k)=tmp_maxk(ii)
                  if(ltrace)write(*,*)'level_clusterDD: fixed minz'
                  fixshow = .true.
               end if
               if (fixmaxz.gt.0 .and. b_maxk(k).eq.tmp_maxk(ii)) then
                  b_mink(k) = b_mink(k) - fixmaxz
                  if (b_mink(k).lt.tmp_mink(ii)) b_mink(k)=tmp_mink(ii)
                  if(ltrace)write(*,*)'level_clusterDD: fixed maxz'
                  fixshow = .true.
               end if
               if (ltrace .and. fixshow) then
                  write(*,*) 'level_clusterDD: After FIX box: ',k
                  write(*,*) '    i: ',b_mini(k),b_maxi(k)
                  write(*,*) '    j: ',b_minj(k),b_maxj(k)
                  write(*,*) '    k: ',b_mink(k),b_maxk(k)
                  write(*,*) ' ---------------------'
               end if
            end do
         else
            if (ltrace) write(*,*) ' Not cutting box'
            if ( maxnum .lt. (boxcount+numx*numy*numz)) then
               write(*,*) ' Not decomposing anymore grids because of '
               write(*,*) ' hardcoded limit to number of child grids.'
               write(*,*) ' Current limit: maxnum =',maxnum
               write(*,*) ' Consider increasing in had/include/glob.inc'
               write(*,*) '       *  *  *  *  *  *  *'
            end if
            boxcount         = boxcount + 1
            b_mini(boxcount) = tmp_mini(ii)
            b_maxi(boxcount) = tmp_maxi(ii)
            b_minj(boxcount) = tmp_minj(ii)
            b_maxj(boxcount) = tmp_maxj(ii)
            b_mink(boxcount) = tmp_mink(ii)
            b_maxk(boxcount) = tmp_maxk(ii)
         end if
      end do

      numbox = boxcount

      if (ltrace) then
         write(*,*) 'level_clusterDD:  numbox = ',numbox
         write(*,*) 'level_clusterDD: DONE.'
         write(*,*) '********************************'
      end if

      return
      end      ! END: level_clusterDD

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_output_flag:                                                        cc
cc                 Pretty prints flag array for debugging of clusterer(s).    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_output_flag(flag,time,bminx,bmaxx, bminy,bmaxy,
     *                             bminz,bmaxz,name,nx,ny,nz,numbox)
      implicit       none
      integer        nx,ny,nz, numbox
      integer        bminx(numbox), bmaxx(numbox),
     *               bminy(numbox), bmaxy(numbox),
     *               bminz(numbox), bmaxz(numbox)
      real(kind=8)   flag(nx,ny,nz), time
      character*(*)  name
      character(60)  line
      character(2)   m_string
      include       'glob.inc'
      integer        i,j,k,   l, m, num, sigi
      integer        sigj(60)
      integer        mystringlength
      external       mystringlength
      logical        marker

      logical      ltrace
      parameter (  ltrace = .false.)

      if (ltrace) then
         write(*,*) 'level_output_flag: Outputting flag array: ',name
         write(*,*) 'level_output_flag: at time:  ',time
         write(*,*) 'level_output_flag: nx,ny,nz: ',nx,ny,nz
      end if

      !
      ! Initialize signatures in y direction:
      !
      do j = 1, ny
         sigj(j) = 0
      end do

      if (.true.) then
         !
         ! ASCII output:
         !
         !call int2str(ny,line)
         !num = mystringlength(line)
         !
         k = (nz+1) / 2
         do j = 1, ny
            line = ''
            sigi = 0
            do i = 1, nx
               !
               ! Mark BEGINS of clusters:
               !
               marker = .false.
               do m = 1, numbox
                  if ( (bminx(m).eq.i) .and.
     *                 (j.ge.bminy(m)) .and. (j.le.bmaxy(m)) .and.
     *                 (k.ge.bminz(m)) .and. (k.le.bmaxz(m)) ) then
                     call int2str(m, m_string)
                     l          = mystringlength(line)
                     line(l+1:) = m_string
                     marker     = .true.
                  end if
               end do
               !if (.not.marker) then
               !   l          = mystringlength(line)
               !   line(l+1:) = '_'
               !end if
               !
               !
               !
               l = mystringlength(line)
               if (flag(i,j,k) .eq. FLAG_DISALLOW) then
                  line(l+1:) = 'O'
               else if (flag(i,j,k) .eq. FLAG_REFINE) then
                  line(l+1:) = 'X'
                  sigi       = sigi    + 1
                  sigj(i)    = sigj(i) + 1
               else if (flag(i,j,k) .eq. FLAG_NOREFINE) then
                  line(l+1:) = '.'
               else 
                  write(*,*) 'level_output_flag: Unknown flag value'
                  stop
               end if
               !
               ! Mark ENDS of clusters:
               !
               marker = .false.
               do m = 1, numbox
                  if ( (bmaxx(m).eq.i) .and.
     *                 (j.ge.bminy(m)) .and. (j.le.bmaxy(m)) .and.
     *                 (k.ge.bminz(m)) .and. (k.le.bmaxz(m)) ) then
                     call int2str(m, m_string)
                     l          = mystringlength(line)
                     line(l+1:) = m_string
                     marker     = .true.
                  end if
               end do
               !if (.not.marker) then
               !   l          = mystringlength(line)
               !   line(l+1:) = '_'
               !end if
            end do
            write(*,90) line, sigi
 90         format(A,I5)
         end do
      end if

      !
      ! Output y signatures:
      !
c     line = ''
c     do i = 1, nx
c        l = mystringlength(line)
c        call int2str(sigj(i),m_string)
c        line(l+1:) = m_string
c        l = mystringlength(line)
c        line(l+1:) = '_'
c     end do
c     write(*,*) line
c     write(*,*) 

      return
      end      ! END: level_output_flag
