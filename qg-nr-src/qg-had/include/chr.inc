      !----------------------------------------------------------------------
      !     Constants characterizing grid points
      !
      !     CHR_interior   --- these points get updated by evolution routines
      !
      !     CHR_amr_bdy    --- these points lie on boundaries of grids which are
      !                        interior to the coarse grid.                        
      !
      !     CHR_Xmin_bdy   --- for X=[x,y,z] are points lying at the *true* outer
      !     CHR_Xmax_bdy       boundary (ie the boundary of the coarse level)
      !
      !     CHR_Xmax_bdy   --- boundary (ie the boundary of the coarse level)
      !
      !     CHR_deco_bdy   --- Similar to CHR_amr_bdy, these points are at (or
      !                        next to) boundary points of one grid that overlap
      !                        with another grid and fall in the interior of the
      !                        overlapping grid.  This type occurs when using     
      !                        domain decomposition in which the boundaries of
      !                        fine grids are either true boundary points or
      !                        overlap into other grids (in which case the boundary
      !                        and one more point next to the boundary point are
      !                        so delineated).
      !
      !----------------------------------------------------------------------

      real(kind=8)         CHR_interior
      real(kind=8)         CHR_amr_bdy
      real(kind=8)         CHR_xmin_bdy
      real(kind=8)         CHR_xmax_bdy
      real(kind=8)         CHR_ymin_bdy
      real(kind=8)         CHR_ymax_bdy
      real(kind=8)         CHR_zmin_bdy
      real(kind=8)         CHR_zmax_bdy
      real(kind=8)         CHR_deco_bdy
      !
      real(kind=8)         CHR_xmin_ymin_bdy
      real(kind=8)         CHR_xmin_ymax_bdy
      real(kind=8)         CHR_xmax_ymin_bdy
      real(kind=8)         CHR_xmax_ymax_bdy
      !
      real(kind=8)         CHR_xmin_zmin_bdy
      real(kind=8)         CHR_xmin_zmax_bdy
      real(kind=8)         CHR_xmax_zmin_bdy
      real(kind=8)         CHR_xmax_zmax_bdy
      !
      real(kind=8)         CHR_ymin_zmin_bdy
      real(kind=8)         CHR_ymin_zmax_bdy
      real(kind=8)         CHR_ymax_zmin_bdy
      real(kind=8)         CHR_ymax_zmax_bdy
      !
      ! Corners:
      !
      real(kind=8)         CHR_xmin_ymin_zmin_bdy
      real(kind=8)         CHR_xmin_ymin_zmax_bdy
      real(kind=8)         CHR_xmin_ymax_zmin_bdy
      real(kind=8)         CHR_xmin_ymax_zmax_bdy
      real(kind=8)         CHR_xmax_ymin_zmin_bdy
      real(kind=8)         CHR_xmax_ymin_zmax_bdy
      real(kind=8)         CHR_xmax_ymax_zmin_bdy
      real(kind=8)         CHR_xmax_ymax_zmax_bdy
      !
      ! Reflection Boundary
      !
      real(kind=8)         CHR_Refl_bdy

      parameter (         CHR_interior          =  0.d0)
      parameter (         CHR_amr_bdy           =  1.d0)
      parameter (         CHR_xmin_bdy          =  2.d0)
      parameter (         CHR_xmax_bdy          =  3.d0)
      parameter (         CHR_ymin_bdy          =  4.d0)
      parameter (         CHR_ymax_bdy          =  5.d0)
      parameter (         CHR_zmin_bdy          =  6.d0)
      parameter (         CHR_zmax_bdy          =  7.d0)
      parameter (         CHR_deco_bdy          =  8.d0)
      !
      ! Edges:
      !
      parameter (         CHR_xmin_ymin_bdy     =  9.d0)
      parameter (         CHR_xmin_ymax_bdy     = 10.d0)
      parameter (         CHR_xmax_ymin_bdy     = 11.d0)
      parameter (         CHR_xmax_ymax_bdy     = 12.d0)
      !
      parameter (         CHR_xmin_zmin_bdy     = 13.d0)
      parameter (         CHR_xmin_zmax_bdy     = 14.d0)
      parameter (         CHR_xmax_zmin_bdy     = 15.d0)
      parameter (         CHR_xmax_zmax_bdy     = 16.d0)
      !
      parameter (         CHR_ymin_zmin_bdy     = 17.d0)
      parameter (         CHR_ymin_zmax_bdy     = 18.d0)
      parameter (         CHR_ymax_zmin_bdy     = 19.d0)
      parameter (         CHR_ymax_zmax_bdy     = 20.d0)
      !
      ! Corners:
      !
      parameter (         CHR_xmin_ymin_zmin_bdy     = 21.d0)
      parameter (         CHR_xmin_ymin_zmax_bdy     = 22.d0)
      parameter (         CHR_xmin_ymax_zmin_bdy     = 23.d0)
      parameter (         CHR_xmin_ymax_zmax_bdy     = 24.d0)
      parameter (         CHR_xmax_ymin_zmin_bdy     = 25.d0)
      parameter (         CHR_xmax_ymin_zmax_bdy     = 26.d0)
      parameter (         CHR_xmax_ymax_zmin_bdy     = 27.d0)
      parameter (         CHR_xmax_ymax_zmax_bdy     = 28.d0)
      !
      ! Reflection Boundary
      !
      parameter (         CHR_Refl_bdy               = 29.d0)

