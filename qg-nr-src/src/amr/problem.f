cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  check_for_collapse:                                                       cc
cc                   Problem-specific check for collapse.                     cc
cc                Checks global variables which are computed in               cc
cc                grid_comp_deriv() to see if conditions met for              cc
cc                termination of code. This routine pulls logic that          cc
cc                used to be in grid_comp_deriv() and puts it in its          cc
cc                own routine called from main.f. The reason being            cc
cc                that previously, the code would be terminated before        cc
cc                appropriate output was done. Now, the code outputs          cc
cc                grid functions, then checks for collapse (singularity or BH)cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_for_collapse( )
      implicit    none
      include    'glob.inc'     ! need "potent"
      include    'param.inc'
      integer     fileunit
       include     'mpif.h'
       include     'mpi_stuff.inc'

      logical     ltrace
      parameter ( ltrace = .false. )
      !
      ! Check if user selected threshold has been exceeded:
      ! Exit if collapse detected.
      ! Outputs "TX" (texture) if collapse detected to 
      !   facilitate automatic critical search (ala "mysearch").
      !
c     Don't exit if max occurred in initial data:
c

      if (ltrace) then
         if ( maxchi .gt. maxchi_thresh) write(*,*) ' max > thresh'
         if ( maxchi_thresh .gt. 0 ) write(*,*) ' thres > 0'
         if ( maxchi_t .gt. maxchi_minctime  ) then
            write(*,*) ' t >minctime'
         else
            write(*,*) ' t !>minctime', maxchi_t, maxchi_minctime
         end if
      end if

      if (       (maxchi        .gt. maxchi_thresh   )
     *     .and. (maxchi_thresh .gt. 0               )
     *     .and. (maxchi_t      .gt. maxchi_minctime ) ) then
         write(*,99) myid,'check_for_collapse: maxchi = ',maxchi
      end if

      if (       (maxchi        .gt. maxchi_thresh   )
     *     .and. (maxchi_thresh .gt. 0               )
     *     .and. (maxchi_t      .gt. maxchi_minctime ) ) then
         fileunit = 7
         open(fileunit, file='TX',form='formatted')
         write(fileunit,*) maxchi,maxchi_t,maxchi_l
         close(fileunit)
         write(*,*) 'check_for_collapse: Chi threshold exceeded.'
         write(*,*) 'check_for_collapse: maxchi_thresh = ',maxchi_thresh
         write(*,*) 'check_for_collapse: maxchi        = ',maxchi
         write(*,*) 'check_for_collapse: maxchi_t      = ',maxchi_t
         write(*,*) 'check_for_collapse: maxchi_l      = ',maxchi_l
         write(*,*) 'check_for_collapse: maxchi_i,j,k  = ',maxchi_i,
     *                                        maxchi_j, maxchi_k
         write(*,*) 'check_for_collapse: Quitting...'
         call my_exit('Singularity Formation.  ')
      end if

 99      format('[',I3,'] ',A,3I5)

      return
      end    ! END: check_for_collapse


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_comp_deriv:                                                          cc
cc                NB: this has basically become a wrapper routine to          cc
cc                    interface to project-specific computation.              cc
cc                   Computes derived quantities on a grid.                   cc
cc                These quantities are not needed for the evolution           cc
cc                but are maintained and ouput. May want a more               cc
cc                sophisticated method such that energy, for example,         cc
cc                is computed using fine grid information.                    cc
cc                   Intefaces to RNPL generated subroutine.                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_comp_deriv( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'param.inc'
      include    'glob.inc'
      integer    nx, ny, nz, numfields,
     *           tmp_maxchi_i, tmp_maxchi_j, tmp_maxchi_k
      integer    proc_return_myid, myid
      real*8     hg,  tmp_maxchi, tmp
      logical     ltrace
      parameter ( ltrace = .false. )
      integer    i,j,k,index, l


      myid = proc_return_myid()
      if (ltrace)  then
         write(*,99) myid,'grid_comp_deriv:Compute err',gridnum
         write(*,99) myid,'grid_comp_deriv: Calling grid_testA:',gridnum
         call grid_test(gridnum)
      end if

      call load_pointers(gridnum)

      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      hg   = gr_h(gridnum)

      !
      ! Compute the "amr error" for this grid as per 
      ! the instructions specific to each project:
      !   (NB: note that this file resides in the project
      !        directory and gets copied here to facilitate
      !        easy compilation)
      !
      call comp_amr_error(gridnum)

      if(ltrace)then
         write(*,99) myid,'grid_comp_deriv: Calling grid_test:',gridnum
         call grid_test(gridnum)
         write(*,99) myid,'grid_comp_deriv: done.',gridnum
      end if

 99      format('[',I4,'] ',A,3I5)

      return
      end    ! END: grid_comp_deriv

