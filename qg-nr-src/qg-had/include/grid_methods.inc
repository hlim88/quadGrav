c----------------------------------------------------------------------
c
c     Functions/methods defined which act on grids or levels:
c
c----------------------------------------------------------------------

      ! 
      !  INTEGER functions:
      ! 
      integer  grid_return_finest,
     *         grid_return_level,
     *         grid_return_parent,
     *         grid_return_child,
     *         grid_return_grandchild,
     *         grid_return_sibling,
     *         grid_return_number,
     *         grid_return_number_norec,
     *         grid_return_owner,
     *         grid_return_count,
     *         grid_kill_children,
     *         grid_return_numgfuncs,
     *         grid_return_ana,
     *         level_return_start,
     *         level_return_finest_ever,
     *         level_return_finest,
     *         level_find_first_child

      external grid_return_finest,
     *         grid_return_level,
     *         grid_return_parent,
     *         grid_return_child,
     *         grid_return_grandchild,
     *         grid_return_sibling,
     *         grid_return_number,
     *         grid_return_number_norec,
     *         grid_return_owner,
     *         grid_return_count,
     *         grid_kill_children,
     *         grid_return_numgfuncs,
     *         grid_return_ana,
     *         level_return_start,
     *         level_return_finest_ever,
     *         level_return_finest,
     *         level_find_first_child

      ! 
      !  LOGICAL functions:
      ! 
      logical  grid_is_refinedgk,
     *         grid_is_refined,
     *         grid_is_childof,
     *         grid_is_contained_within,
     *         grid_return_existence,
     *         grid_alive,
     *         grid_dead,
     *         grid_pending,
     *         grid_at_t0,
     *         grid_is_local,
     *         grid_intersect,
     *         grid_contains_pt,
     *         level_fix_parents,
     *         level_grid_exists_within,
     *         level_return_existence,
     *         level_time_aligned,
     *         field_is_FV

      external grid_is_refinedgk,
     *         grid_is_refined,
     *         grid_is_childof,
     *         grid_is_contained_within,
     *         grid_return_existence,
     *         grid_alive,
     *         grid_at_t0,
     *         grid_is_local,
     *         grid_intersect,
     *         grid_contains_pt,
     *         level_fix_parents,
     *         level_grid_exists_within,
     *         level_return_existence,
     *         level_time_aligned,
     *         field_is_FV

      ! 
      !  REAL*8 functions:
      ! 
      real(kind=8)   grid_return_time,
     *         grid_return_work,
     *         grid_return_h,   
     *         grid_return_resnrm,   
     *         grid_return_minx,
     *         grid_return_miny,
     *         grid_return_minz,
     *         grid_intfind,
     *         level_return_resnrm,   
     *         level_return_work

      external grid_return_time,
     *         grid_return_work,
     *         grid_return_h,   
     *         grid_return_resnrm,   
     *         grid_return_minx,
     *         grid_return_miny,
     *         grid_return_minz,
     *         grid_intfind,
     *         level_return_resnrm,   
     *         level_return_work


