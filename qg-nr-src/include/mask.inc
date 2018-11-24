c
c     ! Mask variables
c

c     ! Maximum number of masks anticipated:
      integer     max_num_masks
      parameter ( max_num_masks  =  2 )

c     ! Masking constants:
      integer     BHMASK_UNMASKED
      integer     BHMASK_MASKED
      integer     BHMASK_REPOPULATE
      parameter ( BHMASK_UNMASKED      =  +1  )
      parameter ( BHMASK_MASKED        =  -1  )
      parameter ( BHMASK_REPOPULATE    =  +10 )

c     ! Number of masks in use at runtime:
      integer     num_masks

c     ! existence or not of an inflow boundary
      logical bh_true(max_num_masks)

c     ! Store whether two black holes have been merged:
      logical      forcedmerge
c     ! Store time at which a transition to one black 
c     ! hole has occurred (either 2 BHs--> 1 BH or 0-->1BH)
      real(kind=8) transition_time
 
      !
      ! For determining if it is time to 
      ! do the mask:
      !
      integer     mask_level, mask_iteration

c     !
c     ! Parameterization of the mask locations:
c     !  NB: likely to change when we have horizon finder
c     !  NB: the three is for the up-to-three coordinates
c     !    
      real(kind=8) mask_center(3,max_num_masks)
      real(kind=8) mask_radius(3,max_num_masks)
      real(kind=8) mask_left(3,max_num_masks)
      real(kind=8) mask_right(3,max_num_masks)
      real(kind=8) local_left(3,max_num_masks)
      real(kind=8) local_right(3,max_num_masks)

      !
      ! Moving to a rectangular mask stored:
      !    mask_bbox(1)   minx
      !    mask_bbox(2)   maxx
      !    mask_bbox(3)   miny
      !    mask_bbox(4)   maxy
      !    mask_bbox(5)   minz
      !    mask_bbox(6)   maxz
      !
      real(kind=8)     mask_bbox(6,max_num_masks)
      real(kind=8) old_mask_bbox(6,max_num_masks)

c     !
c     ! Store the coordinates of what is actually
c     ! excised...separate from what the finder
c     ! *wants* excised. Store (minx,max,miny,maxy...)
c     !
      real(kind=8) mask_coords(6,max_num_masks)

c     !
c     !  Memory to store old masks
c     !     (no need to store in common block)
c     !
      real(kind=8) old_mask_center(3,max_num_masks)
      real(kind=8) old_mask_radius(1,max_num_masks)

c     ! Fuzz-factor so that different levels
c     ! agree on whether a point is in the mask or not
      real(kind=8)  FUZZ
      parameter   ( FUZZ = 1.d-9 )


      common      / mask_common /
     *              mask_coords,
     *              mask_center,
     *              mask_radius,
     *              mask_bbox,
     *              mask_left,
     *              mask_right,
     *              local_left,
     *              local_right,
     *              transition_time,
     *              num_masks,
     *              mask_level,
     *              mask_iteration,
     *              bh_true,
     *              forcedmerge


      !                _____
      !
      !            Apparent Horizon Stuff
      !                _____
      !
      ! Using angle conventions of Pretorius: http://arxiv.org/abs/gr-qc/0407110
      !
      !  For the flow finder (findhorizon==2), one evolves
      !  the surface in some arbitrary time. One can do
      !  either do an explicit update, or the pseudo-implicit iCN:
      logical        HORIZON_USEICN
      parameter    ( HORIZON_USEICN = .true. )

      ! Pointers to the fields defined on the horizon:
      !     horizon_theta   --- the expansion of null geodesics
      !     horizon_R       --- the location of the horizon
      !
      integer        horizon_theta(max_num_masks),
     *               horizon_R(max_num_masks)
      ! Only allocated if using iCN for flow evolution:
      integer        horizon_theta_np1(max_num_masks),
     *               horizon_R_np1(max_num_masks)
      integer        horizon_shape(2)
      real(kind=8)   horizon_center(3,max_num_masks)

      common  / horizon_ints    / horizon_shape,
     *                            horizon_R,       horizon_theta,
     *                            horizon_R_np1,   horizon_theta_np1
      common  / horizon_doubles / horizon_center


