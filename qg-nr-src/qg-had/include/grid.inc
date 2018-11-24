      !
      ! Grid Numbers:
      !
      ! total # of possible grids
      integer     maxnumgrids
      !parameter ( maxnumgrids    = 1024       )
      parameter ( maxnumgrids    = 32768       )

      !
      ! Grid CONSTANTS:
      !
      integer     ALIVE,     DEAD,     PENDING
      parameter ( ALIVE = 1, DEAD = 0, PENDING = -1 )

      include    'fields.inc'
      include    'mem.inc'

      integer     numgrids,                     ! # grids in use
     *            max_numgrids,                 ! max # grids used
     *            num_freegrids,                ! # of free grid numers
     *            max_num_freegrids,            ! max # free at 1 time
     *            free_grid_numbers(maxnumgrids)! array free grid nmbrs

      !
      ! Quantities and fields that live in a grid:
      !
      real(kind=8)      gr_h(maxnumgrids),
     *            gr_t(maxnumgrids),
     *            gr_work(maxnumgrids),
     *            gr_maxx(maxnumgrids),
     *            gr_maxy(maxnumgrids),
     *            gr_maxz(maxnumgrids),
     *            gr_minx(maxnumgrids),
     *            gr_miny(maxnumgrids),
     *            gr_minz(maxnumgrids)
      integer
     *            gr_own(maxnumgrids),
     *            gr_count(maxnumgrids),
     *            gr_level(maxnumgrids),
     *            gr_parent(maxnumgrids),
     *            gr_child(maxnumgrids),
     *            gr_sibling(maxnumgrids),
     *            gr_size(maxnumgrids),
     *            gr_nx(maxnumgrids),
     *            gr_ny(maxnumgrids),
     *            gr_nz(maxnumgrids),
     *            gr_x(maxnumgrids),
     *            gr_y(maxnumgrids),
     *            gr_z(maxnumgrids),
     *            gr_alive(maxnumgrids),
     *            gr_iter(maxnumgrids)

      common      / grid_structure /
     *            gr_h,       gr_t,       gr_work,
     *            gr_maxx,    gr_maxy,    gr_maxz,
     *            gr_minx,    gr_miny,    gr_minz,
     *            gr_own,     gr_count,   gr_level,
     *            gr_parent,  gr_child,   gr_sibling,
     *            gr_size,
     *            gr_nx,      gr_ny,      gr_nz,
     *            gr_x,       gr_y,       gr_z,
     *            gr_alive,   gr_iter,
     *            free_grid_numbers,
     *            numgrids, max_numgrids,
     *            num_freegrids,
     *            max_num_freegrids
