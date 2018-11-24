c.......................................................................
c
c
c globs.inc:   Global variables to the distributed AMR code
c
c
c.......................................................................

      !
      ! Maximum number of Levels possible
      !
      integer     maxlev
      parameter ( maxlev         =  25        )

      !
      ! Various other globals:
      !
      ! epsilon for K-O dissipation (e.g. 0.05)
      real(kind=8)  epsdis,
     *        lambda,           ! CFL ratio dt/dx (e.g. 0.3)
     *        minx0, miny0, minz0,! User selected bounds on coarse grid
     *        maxx0, maxy0, maxz0,
     *        ref_level(0:maxlev-1), ! refinement factor for level
     *        ell_epsilon,      ! Threshold for residual in elliptic
     *        refine_deltat,    ! how often to refine in units of time
     *        mask_period_dt,   ! how often to find a horizon after transition
     *        horizon_thresh,   ! threshold on expansion for horizon 
     *        horizon_growth,   ! factor by which initial guess is increased
     *        output_f1_lb1,    ! lower bound point for family 1 output
     *        output_f1_lb2,    ! lower bound point for family 1 output
     *        output_f1_lb3,    ! lower bound point for family 1 output
     *        output_f1_ub1,    ! upper bound point for family 1 output
     *        output_f1_ub2,    ! upper bound point for family 1 output
     *        output_f1_ub3,    ! upper bound point for family 1 output
     *        output_f2_lb1,    ! lower bound point for family 2 output
     *        output_f2_lb2,    ! lower bound point for family 2 output
     *        output_f2_lb3,    ! lower bound point for family 2 output
     *        output_f2_ub1,    ! upper bound point for family 2 output
  	  ! Variable timestep modification by Dominic Marcello 2009-08-04
     *        cfl_lambda,        ! cfl factor for variable timestep
	  ! --------------end modification by Dominic Marcello 2009-08-04
     *        output_f2_ub2,    ! upper bound point for family 2 output
     *        output_f2_ub3,    ! upper bound point for family 2 output
     *        output_f3_lb1,    ! lower bound point for family 3 output
     *        output_f3_lb2,    ! lower bound point for family 3 output
     *        output_f3_lb3,    ! lower bound point for family 3 output
     *        output_f3_ub1,    ! upper bound point for family 3 output
     *        output_f3_ub2,    ! upper bound point for family 3 output
     *        output_f3_ub3,    ! upper bound point for family 3 output
     *        output_f4_lb1,    ! lower bound point for family 4 output
     *        output_f4_lb2,    ! lower bound point for family 4 output
     *        output_f4_lb3,    ! lower bound point for family 4 output
     *        output_f4_ub1,    ! upper bound point for family 4 output
     *        output_f4_ub2,    ! upper bound point for family 4 output
     *        output_f4_ub3     ! upper bound point for family 4 output

      integer trace_level,     ! determines how much output tracing
                               !   0      minimal time step output
                               ! ( 0,10]  output from main() only
     *        steps,           ! loop counter for coarse grid evolution
     *        refine_factor,   ! ratio dx_1/dx_2 of grid resolutions
                               ! must be even
     *        weno_interp,      ! Determine if using weno or not                         
     *        refine_period,   ! number of time steps between refinement
     *        refine_period_ctrl, ! controls whether period is in units of t or steps
     *        simpleFMR,       ! flag to use simple volume based FMR
     *        clusterDD,       ! controls whether domain decomposition 
     *        clusterstyle,    ! restricts clustering
     *        diss_afterinj,   ! add dissipation after injection?(if!=0)
     *        Levelp(0:maxlev-1),! pointer to 'first' grid @ given level
                               ! negative if no grids @ that level
     *        out_period(0:maxlev-1),! output period for each level
     *        lev_count(0:maxlev-1), ! step counter for various levels
     *        output_style,  ! output by: level(1), together(2), both(3)
     *        output_dim,    ! dimension of output
     *        output_f1_type,   ! type  of output for family 1
     *        output_f1_level,  ! level of output for family 1
     *        output_f1_period, ! # of time steps between output for family 1
     *        output_f2_type,   ! type  of output for family 2
     *        output_f2_level,  ! level of output for family 2
     *        output_f2_period, ! # of time steps between output for family 2
     *        output_f3_type,   ! type  of output for family 3
     *        output_f3_level,  ! level of output for family 3
     *        output_f3_period, ! # of time steps between output for family 3
     *        output_f4_type,   ! type  of output for family 4
     *        output_f4_level,  ! level of output for family 4
     *        output_f4_period, ! # of time steps between output for family 4
     *        linearbounds,  ! use linear interp if > 0
     *        shadow,        ! shadow hierarchy
     *        allowedl,      ! maximum # of levels allowed by user
     *        num_evol_iters,! # of iterations for an evol step
     *        bound_width,   ! width of overlap region for DD
     *        num_injfields, ! number of fields which get injected
     *        num_synfields, ! number of fields which get synced
     *        num_anafields, ! number of fields requiring analysis
     *        chkpt_period,  ! how often to dump state for checkpoint
     *        chkpt_readstate,! how often to dump state for checkpoint
     *        chkpt_control, ! whether to dump state in turns or at same time
     *        chkpt_fnum,    ! number/index to chkpt_fname array
     *        findhorizon,   ! how to find a horizon (or use static)
     *        mask_period,   ! how often to find a horizon
     *        mask_usemin,   ! use min of lapse to determine mask guess?
     *        mask_minfield, ! number of field to take minimum of
     *        horizon_ntheta, !number of points on horizon surface
     *        horizon_nphi,   !number of points on horizon surface
     *        horizon_nholes, !number of masked regions to search for
     *        horizon_recenterp, !period for recentering surface (default=0 never)
     *        emulate_proc,   !processor to emulate (-1 to turn off)
     *        assume_symmetry,!option to enforce a symmetry on domain (default=0 off)
  	  ! Variable timestep modification by Dominic Marcello 2009-08-04
     *        variable_timestep, ! whether or not to use a variable timestep
	  ! --------------end modification by Dominic Marcello 2009-08-04
     *        flush_period,   !how often to close SDF files to flush (default=0 off)
     *        periodicBC      !option to enforce periodic boundary condtions (default=0 off)

       character(4) tag      ! prefix for output sdf files
       character(128) chkpt_fname(2) ! name of file for saving state

      !
      ! Clustering parameters:
      !
       real(kind=8)  ethreshold,! error threshold to flag for rfnmnt
     *         ethreshold_orig, ! original value so that it can be reset
     *         minefficiency ! minimum efficiency for new refined grids
       integer buffer,  ! number of coarse grid pts around flagged
     *         mindim,  ! smallest # of coarse grid pts for subgrid
     *         window,  ! # of coarse grid pts from bdy in which
                                     !   clusterer looks to slice box
     *         ghostwidth,        ! Width of region set  by parent
     *         update_scheme,     ! choose RK3 or ICN
     *         amrbound_prepost,  ! enforce subgrd bnds before evolve?
     *         amrbound_timealign,! enforce subgrd only timealigned?
     *         elliptic_solve,    ! ==0 no elliptic solve, >0 then do so
     *         nvcycle,           ! # of vcycles to carry-out
     *         preswp,            ! # of pre-sweeps  of relax. in vcycle
     *         pstswp,            ! # of post-sweeps of relax. in vcycle
     *         maxsweeps,         ! max # of relaxation sweeps for solve
     *         clusterreadwrite   ! Flag:
                                  !   0 --- normal clustering
                                  !   1 --- read  in  hierarchy information
                                  !   2 --- write out hierarchy information

       real(kind=8) total_work! diagnostic measure of integrated "work"
                              !   done relative to coarse grid.

      !
      ! Store maximum value of chi
      !   to detect singularity formation:
      !
      real(kind=8)  maxchi,   ! maximum (absolute) value of chi
     *        maxchi_t        ! time at which max is achieved
      integer maxchi_i, maxchi_j, maxchi_k, ! gridpoint of maximum
     *        maxchi_l,       ! level of maximum
     *        maxchi_reach    ! =0 if not reached collapsed time

      !
      ! Interrupt handling:
      !
      integer interrupt_set,  !  set >0 if interrupt detected
     *        interrupt_rc    !  return code for signal() call

      !
      ! File unit number for output of runtime statistics
      !
      integer fileout

      !
      ! Bounding boxes for refinement: (not in a common block)
      !

      integer     maxnumchildren ! total # of children for any grid
      parameter ( maxnumchildren = 16384       )
      !parameter ( maxnumchildren = 2048       )
      !parameter ( maxnumchildren = 512        )

      !
      ! Flagging constants:
      !
      !     ...points     to be refined
      real(kind=8)   FLAG_REFINE
      parameter    ( FLAG_REFINE   =  1.d0 )
      !     ...points NOT to be refined
      real(kind=8)   FLAG_NOREFINE
      parameter    ( FLAG_NOREFINE =  0.d0 )
      !     ...points which do not exist on the level
      real(kind=8)   FLAG_DISALLOW
      parameter    ( FLAG_DISALLOW = -1.d0 )

      common  / myglobals_doubles /
     *        ref_level,     epsdis,        lambda,
     *        minx0,         miny0,         minz0,
     *        maxx0,         maxy0,         maxz0,
     *        ethreshold,    ethreshold_orig, minefficiency, total_work,
     *        maxchi,        maxchi_t,      ell_epsilon,
     *        refine_deltat, mask_period_dt, horizon_thresh,
     *        horizon_growth,
     *        output_f1_lb1, output_f1_lb2, output_f1_lb3,
     *        output_f1_ub1, output_f1_ub2, output_f1_ub3,
     *        output_f2_lb1, output_f2_lb2, output_f2_lb3,
  	  ! Variable timestep modification by Dominic Marcello 2009-08-04
     *        cfl_lambda, 
	  ! --------------end modification by Dominic Marcello 2009-08-04
     *        output_f2_ub1, output_f2_ub2, output_f2_ub3,
     *        output_f3_lb1, output_f3_lb2, output_f3_lb3,
     *        output_f3_ub1, output_f3_ub2, output_f3_ub3,
     *        output_f4_lb1, output_f4_lb2, output_f4_lb3,
     *        output_f4_ub1, output_f4_ub2, output_f4_ub3


      common  / myglobals_ints /
     *        out_period,    lev_count,      Levelp,
     *        steps,         ghostwidth,     update_scheme,
     *        amrbound_prepost,              amrbound_timealign,
     *        elliptic_solve, clusterreadwrite,
     *        nvcycle, preswp, pstswp, maxsweeps,
     *        window,        buffer,         mindim,
     *        trace_level,   refine_factor,  linearbounds, shadow,
     *        weno_interp,
     *        allowedl,num_evol_iters, bound_width, num_injfields,
     *        num_synfields, num_anafields,
     *        chkpt_period, chkpt_readstate, chkpt_control, chkpt_fnum,
     *        findhorizon,  mask_period, mask_usemin, mask_minfield,
     *        horizon_ntheta, horizon_nphi, horizon_nholes,
     *        horizon_recenterp, emulate_proc, assume_symmetry,

  	  ! Variable timestep modification by Dominic Marcello 2009-08-04
     *        variable_timestep, 
	  ! --------------end modification by Dominic Marcello 2009-08-04

     *        flush_period,  periodicBC,
     *        refine_period, refine_period_ctrl,
     *        clusterDD, clusterstyle, simpleFMR,
     *        diss_afterinj,  output_style, output_dim,
     *        output_f1_type, output_f1_level, output_f1_period,
     *        output_f2_type, output_f2_level, output_f2_period,
     *        output_f3_type, output_f3_level, output_f3_period,
     *        output_f4_type, output_f4_level, output_f4_period,
     *        maxchi_i,      maxchi_j,       maxchi_k,
     *        maxchi_l,      maxchi_reach,   
     *        interrupt_set,  fileout

      common  / myglobl_strings / tag, chkpt_fname
