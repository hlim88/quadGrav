c-----------------------------------------------------------------------
c  2 Jun 2008 
c  This routine gathers all the designated field for the specified level
c  into a new array which is a pointer into the arena.  Used 
c  with Mijan's horizon finder.
c-----------------------------------------------------------------------
        subroutine level_gather_field(level,fieldnum,field,nx,ny,nz)
        implicit none
        include     'grid_methods.inc'
        include     'grid.inc'
        !include     'glob.inc'
        include     'mpif.h'
        include     'mpi_stuff.inc'

c        character(*) fieldname
        integer fieldnum
        integer nx,ny,nz,level
        real(kind=8) field(nx,ny,nz)

c       local variables
        integer gi,owner
        real(kind=8) minx, maxx, miny, maxy, minz, maxz
        integer      nxl,nyl,nzl
        integer      nxg,nyg,nzg, length
        integer      i,j,k,index
        integer      mini,minj,mink
        real*8, pointer, dimension(:,:,:) :: tmp

        logical     ltrace
        parameter ( ltrace = .false. )

 99   format('[',I3,'] ',A,3I5)

      if (ltrace) then
         write(*,99)myid,'level_gather_field: Entering',fieldnum
      end if

      !
      !  loop over all grids on the level
      !
      gi = level_return_start(level)
 10   if (grid_return_existence(gi)) then
          owner  = grid_return_owner(gi)
          nxg    = gr_nx(gi)
          nyg    = gr_ny(gi)
          nzg    = gr_nz(gi)
          length = nxg*nyg*nzg
          if(ltrace)write(*,99)myid,'level_gather_field: gi=',gi,owner
          if(ltrace)write(*,99)myid,'level_gather_field:nxg',nxg,nyg,nzg
          if (myid.eq.master) then
             allocate(tmp(nxg,nyg,nzg))
             if (owner .ne. master) then
                !
                !...MPI recv data into tmp...
                !
                if(ltrace)write(*,99)myid,'level_gather_field: Recving'
                if(ltrace)write(*,99)myid,'length = ',length
                if(ltrace)write(*,99)myid,'owner = ',owner
                if(ltrace)write(*,99)myid,'TAG_GATHER = ',TAG_GATHER
                call MPI_Recv(tmp,length,
     &           MPI_DOUBLE_PRECISION,owner,TAG_GATHER,
     &           MPI_COMM_WORLD,status,ierr)
             else
                !
                !...Grid is local, just copy...
                !
                if(ltrace)write(*,99)myid,'level_gather_field: Copying'
                !call load_pointers(gi)
                do k=1,nzg
                  do j=1,nyg
                    do i=1,nxg
                      index = (k-1)*nxg*nyg + (j-1)*nxg + (i-1)
                      tmp(i,j,k) = q(gfunc_pointer(fieldnum,gi)+index)
                    end do
                  end do
                end do
             end if
             !
             !...Combine tmp data into allocated arrays...
             !
             if(ltrace)write(*,99)myid,'level_gather_field: Combingin'
             !
             if(ltrace)write(*,99)myid,'level_gather_field:Getting bnds'
             call level_find_bounds(level,minx,maxx,miny,maxy,minz,maxz)
             !
             mini = NINT( (gr_minx(gi)-minx)/gr_h(gi) )
             minj = NINT( (gr_miny(gi)-miny)/gr_h(gi) )
             mink = NINT( (gr_minz(gi)-minz)/gr_h(gi) )
             if(ltrace)write(*,99)myid,'level_gather:min',mini,minj,mink
             do k=1,nzg
               do j=1,nyg
                 do i=1,nxg
                   index = (mink+k-1)*nx*ny+(minj+j-1)*nx+(mini+i-1)
                   field(mini+i,minj+j,mink+k) = tmp(i,j,k)
                 end do
               end do
             end do
             deallocate(tmp)
          else if (myid.eq.owner) then
             !
             !...Send data to master...
             !
             if(ltrace)write(*,99)myid,'level_gather_field: Sending'
             if(ltrace)write(*,99)myid,'fieldnum = ',fieldnum
             if(ltrace)write(*,99)myid,'gi = ',gi
             if(ltrace)write(*,99)myid,'length = ',length
             if(ltrace)write(*,99)myid,'master = ',master
             if(ltrace)write(*,99)myid,'TAG_GATHER = ',TAG_GATHER
             if(ltrace)write(*,*)'gfunc_ptr=',gfunc_pointer(fieldnum,gi)
             call MPI_Send(q(gfunc_pointer(fieldnum,gi)), length,
     &         MPI_DOUBLE_PRECISION,
     &         master, TAG_GATHER, MPI_COMM_WORLD, ierr)
          else
             !...do nothing...
             if(ltrace)write(*,99)myid,' doing nothing'
          end if
          gi = grid_return_sibling(gi)
          if(ltrace)write(*,99)myid,'level_gather_field: Next grid',gi
          goto 10
        end if

        end subroutine
