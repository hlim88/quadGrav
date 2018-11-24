      !
      ! Output functions found in bbhutil:
      !
      integer     gft_out_full,     gft_out_bbox
      external    gft_out_full,     gft_out_bbox
      !
      integer     gft_write_id_int, gft_write_id_float, gft_write_id_str
      external    gft_write_id_int, gft_write_id_float, gft_write_id_str
      integer     mygft_write_id_int, mygft_write_id_float
      external    mygft_write_id_int, mygft_write_id_float
      integer     mygft_read_id_int, mygft_read_id_float
      external    mygft_read_id_int, mygft_read_id_float
      !
      integer     gft_read_id_int,  gft_read_id_float,  gft_read_id_str
      external    gft_read_id_int,  gft_read_id_float,  gft_read_id_str
      !
      integer     gft_read_id_gf,   gft_write_id_gf
      external    gft_read_id_gf,   gft_write_id_gf


      !
      ! Variables to use the above routines:
      !
      integer       gft_rank
      integer       gft_shape(3)
      character(64) gft_cnames
      character(20) gft_name
      integer       gft_rc

      parameter  (  gft_rank   = 3       )
      parameter  (  gft_cnames = 'x|y|z' )

