
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! These fields and variables are for the tracer particles.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer     tracers_number 
      parameter(  tracers_number = 1000 )
      real(kind=8)tracers_t
      real(kind=8)tracers_x(tracers_number)
      real(kind=8)tracers_y(tracers_number)
      real(kind=8)tracers_z(tracers_number)
      real(kind=8)tracers_xprev(tracers_number)
      real(kind=8)tracers_yprev(tracers_number)
      real(kind=8)tracers_zprev(tracers_number)
      !
      real(kind=8)tracers_Vx(tracers_number)
      real(kind=8)tracers_Vy(tracers_number)
      real(kind=8)tracers_Vz(tracers_number)
      real(kind=8)tracers_Vxprev(tracers_number)
      real(kind=8)tracers_Vyprev(tracers_number)
      real(kind=8)tracers_Vzprev(tracers_number)
      !
      ! Temporary storage for MPI transfers...does not need to be stored
      real(kind=8)tracers_tmp(tracers_number)
      ! will need to include various fluid characteristics that
      ! will be interpolated and stored for each particle, such as rho, Ye, etc


      !
      !
      !            How many coarse level time steps in between updates
      !               of the tracer particles? Set by parameter
      integer      tracers_period
      !
      !
      !            How to initialize the tracer particles?
      !               ==0 [default] do nothing, continue from checkpoint
      !               ==1 set in low density regions
      !               > 1 unusued so far
      integer      tracers_initial
      !
      !
      !            Defines the update scheme for tracer particles
      !                (such as: leap-frog, PIC, null-geodesic, etc) unused so far
      integer      tracers_scheme
      common  / tracers /  tracers_x,
     *                     tracers_y,
     *                     tracers_z,
     *                     tracers_xprev,
     *                     tracers_yprev,
     *                     tracers_zprev,
     *                     tracers_Vx,
     *                     tracers_Vy,
     *                     tracers_Vz,
     *                     tracers_Vxprev,
     *                     tracers_Vyprev,
     *                     tracers_Vzprev,
     *                     tracers_t,
     *                     tracers_period,
     *                     tracers_initial,
     *                     tracers_scheme
