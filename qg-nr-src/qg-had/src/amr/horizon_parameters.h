c-----------------------------------------------------------------------
c	Header file for Apparent Horizon locator.
c	Defines global variables used in the program.
c-----------------------------------------------------------------------
c		Include file for parameters.
c	
c			Maximum sizes for surface mesh arrays
		integer  NSTM,NSPM

c			Maximum sizes for Cartesian grid. (Cartesian grid version)

		integer	NXM,NYM,NZM

		parameter(NSTM=50)
		parameter(NSPM=50)

		parameter(NXM=160)
		parameter(NYM=160)
		parameter(NZM=160)

c			Number of holes. (See schwarzschild.f)
		integer	nholes

      real*8   Mass

c     Value of the outgoing divergence to solve to 
      real*8   kappa0

c			Radius of initial starting surface.
c				(Could be ignored if surface data is read in from file.) 
      real*8   radius

c			Finite Difference molecule spacing.
      real*8   h

c			Spherical mesh spacings
      real*8   dphi
      real*8   dtheta

      real*8   Pi

c			Offsets
      real*8   x0,y0,z0

		real*8	tau

		real*8	center(3)

		real*8	M1,M2
		
c		Jacobian generation. \delta \chi
		real*8	epsilon
	
		real*8	mom(2,3), spin(2,3), cntr(2,3), ahole(2)
	
		logical  VERBOSE

c		Stopping Criterion

		real*8	newtst, ilucgstp

c		Stuff for scott's black box extrap...
		integer soc(4)
        data soc/4,10,20,35/


c  Declare common variables
      common   /Eval/kappa0

      common   /Toggles/VERBOSE

      common   /Parameter1/h,dphi,dtheta,Pi,Mass,radius

      common   /Solver/epsilon

      common   /Constants/x0,y0,z0,tau

		common	/Coordinates/center

		common	/Holes/M1,M2,mom,spin, cntr, ahole, nholes

		common	/Stopping/newtst, ilucgstp


c		common/expol/soc
