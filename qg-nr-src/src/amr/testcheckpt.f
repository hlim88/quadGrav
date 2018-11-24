      program test

      implicit none

      integer  i, j
      real(kind=8) x, y
      real(kind=8) xarray(30), yarray(30)

      integer  gft_rc
      integer  mygft_write_id_int
      external mygft_write_id_int
      integer  mygft_read_id_int
      external mygft_read_id_int
      integer  mygft_write_id_float
      external mygft_write_id_float
      integer  mygft_read_id_float
      external mygft_read_id_float
      character(12) fname


      fname  = '.mytestfile'

      i = 1111
      gft_rc = mygft_write_id_int(fname,'i',i,1)
      x = 1111.d0
      gft_rc = mygft_write_id_float(fname,'x',x,1)

      do i = 1, 30
         xarray(i) = 1.d0*i*i
      end do
      gft_rc = mygft_write_id_float(fname,'xarray',xarray,30)
      call gft_close_all()  

      gft_rc = mygft_read_id_int(fname,'i',j,1)
      gft_rc = mygft_read_id_float(fname,'x',y,1)

      write(*,*) ' write i: ',i,' and read into j: ',j
      write(*,*) ' write x: ',x,' and read into y: ',y

      gft_rc = mygft_read_id_float(fname,'xarray',yarray,30)
      do i = 1, 30
          write(*,*) i,' write x(i): ',xarray(i),' and read into y(i):',
     .                yarray(i)
      end do
      end     ! END: masks_return_num
