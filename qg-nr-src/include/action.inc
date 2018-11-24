c----------------------------------------------------------------------
c     Constants characterizing various operations defined on grids
c        or levels:
c
c        Used by level_apply and for inter-processor communication.
c
c----------------------------------------------------------------------
      integer  STEP,

  	  ! Variable timestep modification by Dominic Marcello 2009-08-03
     *         CFL,
	  ! --------------end modification by Dominic Marcello 2009-08-03

     *         INJECT,
     *         BOUNDS,
     *         SWAP_N,
     *         OUTPUT,
     *         REFINE,
     *         COMPUT,
     *         NEWGRID,
     *         DELGRID,
     *         CREATEGRID,
     *         DONE,
     *         NONE,
     *         DUMP,
     *         BOUNDSM,
     *         INJECTM,
     *         LEVINC,
     *         KILLKIDS,
     *         MKALIVE,
     *         FIXGRID,
     *         GRIDINIT,
     *         S_INIT,
     *         QUIT,
     *         RESETMSK,
     *         STEP2,
     *         SYNCBND,
     *         SETMSK,
     *         BCASTPAR,
     *         GETFRALL,
     *         GETERROR,
     *         NEWGRIDA,
     *         CREATEGRIDA,
     *         FREEGRID,
     *         RELAX,
     *         SYNCBNDI,
     *         MG_CGCCOMPUTE,
     *         MG_CGCTRANSFER,
     *         MG_CGCTRANSFERM,
     *         MG_PSTCOMPUTE,
     *         MG_PSTTRANSFER,
     *         MG_PSTTRANSFERM,
     *         MY_BARRIER,
     *         BHMASK,
     *         INIT,
     *         RESETTIME,
     *         SAVESTATE,
     *         READ_STATE,
     *         ELLRHS,
     *         ELLPOST,
     *         SURFCOMP,
     *         ELLSOLVE,
     *         TRACERINIT

      parameter (
     *         STEP          =  0,
  	  ! Variable timestep modification by Dominic Marcello 2009-08-03
     *         CFL           = 50,
	  ! --------------end modification by Dominic Marcello 2009-08-03
     *         INJECT        =  1,
     *         BOUNDS        =  2,
     *         SWAP_N        =  3,
     *         OUTPUT        =  4,
     *         REFINE        =  5,
     *         COMPUT        =  6,
     *         NEWGRID       =  7,
     *         DELGRID       =  8,
     *         CREATEGRID    =  9,
     *         DONE          =  10,
     *         NONE          =  11,
     *         DUMP          =  12,
     *         BOUNDSM       =  13,
     *         INJECTM       =  14,
     *         LEVINC        =  15,
     *         KILLKIDS      =  16,
     *         MKALIVE       =  17,
     *         FIXGRID       =  18,
     *         GRIDINIT      =  19,
     *         S_INIT        =  20,
     *         QUIT          =  21,
     *         RESETMSK      =  22,
     *         STEP2         =  23,
     *         SYNCBND       =  24,
     *         SETMSK        =  25,
     *         BCASTPAR      =  26,
     *         GETFRALL      =  27,
     *         GETERROR      =  28,
     *         NEWGRIDA      =  29,
     *         CREATEGRIDA   =  30,
     *         FREEGRID      =  31,
     *         RELAX         =  32,
     *         SYNCBNDI      =  33,
     *         MG_CGCCOMPUTE =  34,
     *         MG_CGCTRANSFER=  35,
     *         MG_CGCTRANSFERM= 36,
     *         MG_PSTCOMPUTE =  37,
     *         MG_PSTTRANSFER=  38,
     *         MG_PSTTRANSFERM= 39,
     *         MY_BARRIER    =  40,
     *         BHMASK        =  41,
     *         INIT          =  42,
     *         RESETTIME     =  43,
     *         SAVESTATE     =  44,
     *         READ_STATE    =  45,
     *         ELLRHS        =  46,
     *         ELLPOST       =  47,
     *         SURFCOMP      =  48,
     *         ELLSOLVE      =  49,
     *         TRACERINIT    =  50
     *         )


  !.................................................................
  !
  ! INJECT     ---     1-1 message passing from parent grid
  !                                        to   child  grid
  !
  ! CREATEGRID ---     1-1 message passing from parent grid
  !                                        to   child  grid
  !
  ! REFINE     ---     1-1 message passing from master processor
  !                                        to   child  grid
  !
  ! NEWGRID    ---     collective passing. Loops through
  !                                        processors. Starts with
  !                                        owner of grid, looping
  !                                        through all processors,
  !                                        but skips the processor
  !                                        which owns the parent grid
  !                                        because that processor is 
  !                                        already aware of grid.
  !
  ! LEVINC     ---     collective passing. Loops through all procs
  !                                        returning to originator proc.
  !
  ! DONE       ---     collective passing. Loops through all procs
  !                                        returning to originator proc.
  !
  ! SYNCBNDI   ---     same as SYNCBND but has no role in evolution,
  !                                        just the Initial data
  !
  !.................................................................
