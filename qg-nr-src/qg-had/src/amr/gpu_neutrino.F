c-----------------------------------------------------------------------
c  This routine is called by level_lib to interface with 
c  the gpu neutrino transport code
c-----------------------------------------------------------------------
        subroutine gpu_neutrino(level)
        implicit none
        integer  gi
        include     'grid_methods.inc'
        include     'grid.inc'
        include     'glob.inc'
        include     'mask.inc'
        include     'mpif.h'
        include     'mpi_stuff.inc'
        include     'param.inc'
        integer level,i,j,k,flag
        real*8, pointer, dimension(:,:,:) :: etemp
        real*8, pointer, dimension(:,:,:) :: edens
        real*8, pointer, dimension(:) :: xx
        real*8, pointer, dimension(:) :: yy
        real*8, pointer, dimension(:) :: zz
        real(kind=8) hl
        real(kind=8) minx, maxx, miny, maxy, minz, maxz
        integer      nxl,nyl,nzl
        integer      readflag
        !
        character(5) fieldname
        logical      ltrace
        parameter (  ltrace = .false. )


 99   format('[',I3,'] ',A,3I5)

c find the bounds of the level being used for horizon finding
        call level_find_bounds(level,minx,maxx,miny,maxy,minz,maxz)

        gi = level_return_start(level)
        hl    = grid_return_h(gi)
        nxl   = NINT( (maxx-minx)/hl ) + 1
        nyl   = NINT( (maxy-miny)/hl ) + 1
        nzl   = NINT( (maxz-minz)/hl ) + 1

      if (ltrace) then
         write(*,99) myid, 'gpu_neutrino: Called on level=',level
         write(*,99) myid, 'gpu_neutrino: nx/y/zl:',nxl,nyl,nzl
      end if

c allocate memory big enough to hold what we need
        if (myid.eq.master) then
          if(ltrace)write(*,99)myid,'mijans_horizon: Allocating memory'
          allocate(etemp(nxl,nyl,nzl))
          allocate(edens(nxl,nyl,nzl))
          allocate(xx(nxl))
          allocate(yy(nyl))
          allocate(zz(nzl))
        end if

c       gather the adm variables from the finest mesh onto 
c       the master processor
        if (ltrace) write(*,*) myid, 'gpu_neutrino: Gathering'
c the fieldnumber is found in had/include/fields.inc
        call level_gather_field(level,525,etemp,nxl,nyl,nzl) 
        call level_gather_field(level,505,edens,nxl,nyl,nzl) 
        if (ltrace) write(*,*) myid, 'gpu_neutrino: Done Gathering'

        if (myid.eq.master) then
c  initialize the coordinate arrays
          do i=1,nxl
            xx(i) = minx + hl*(i-1)
          end do
          do j=1,nyl
            yy(j) = miny + hl*(j-1)
          end do
          do k=1,nzl
            zz(k) = minz + hl*(k-1)
          end do

c  call the neutrino transport code
#ifdef CUDA
          call Neutrino(etemp,edens,xx,yy,zz)
#endif

c write it out to an sdf file from the master processor --debug
          if (ltrace) call field_out3d(etemp,0.0,"etemp",
     &             0.0,1.0,0.0,1.0,0.0,1.0,
     &            nxl,nyl,nzl,0)

          deallocate(etemp)
          deallocate(edens)
          deallocate(xx)
          deallocate(yy)
          deallocate(zz)
        end if

        end subroutine
