cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  update1:                                                                  cc
cc             Compute derived quantities                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update1( chi, chit, edens, error, dx, nx,ny,nz)
      implicit none
      integer  nx,ny,nz
      real(kind=8)   chi(nx,ny,nz),        chit(nx,ny,nz),
     *         edens(nx,ny,nz),      error(nx,ny,nz),
     *         dx
      real(kind=8)   chi_x,  chi_y, chi_z
      integer  i,j,k

      logical     ltrace
      parameter ( ltrace = .false. )

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
            if (i.eq.1) then
               chi_x = (  - 3.d0 * chi(i,  j,k) 
     *                    + 4.d0 * chi(i+1,j,k)
     *                    -        chi(i+2,j,k)  ) / (2.d0*dx)
            else if (i.eq.nx) then
               chi_x = (  + 3.d0 * chi(i,  j,k) 
     *                    - 4.d0 * chi(i-1,j,k)
     *                    +        chi(i-2,j,k)  ) / (2.d0*dx)
            else
               chi_x = (           chi(i+1,j,k) 
     *                    -        chi(i-1,j,k)  ) / (2.d0*dx)
            end if
            if (j.eq.1) then
               chi_y = (  - 3.d0 * chi(i,j  ,k) 
     *                    + 4.d0 * chi(i,j+1,k)
     *                    -        chi(i,j+2,k)  ) / (2.d0*dx)
            else if (j.eq.ny) then
               chi_y = (  + 3.d0 * chi(i,j  ,k) 
     *                    - 4.d0 * chi(i,j-1,k)
     *                    +        chi(i,j-2,k)  ) / (2.d0*dx)
            else
               chi_y = (           chi(i,j+1,k) 
     *                    -        chi(i,j-1,k)  ) / (2.d0*dx)
            end if
            if (k.eq.1) then
               chi_z = (  - 3.d0 * chi(i,j,k  ) 
     *                    + 4.d0 * chi(i,j,k+1)
     *                    -        chi(i,j,k+2)  ) / (2.d0*dx)
            else if (k.eq.nz) then
               chi_z = (  + 3.d0 * chi(i,j,k  ) 
     *                    - 4.d0 * chi(i,j,k-1)
     *                    +        chi(i,j,k-2)  ) / (2.d0*dx)
            else
               chi_z = (           chi(i,j,k+1) 
     *                    -        chi(i,j,k-1)  ) / (2.d0*dx)
            end if
         edens(i,j,k) = 0.5d0 * (    chit(i,j,k)**2
     *                            +  chi_x**2 + chi_y**2 + chi_z**2 )
         error(i,j,k) = dx**2 * edens(i,j,k)
      end do
      end do
      end do

      return
      end      ! END: update1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  calc_rhs:                                                                 cc
cc             Calculate RHS's for Klein Gordon equation                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_rhs( chi_rhs, chit_rhs, chi, chit, CHR,
     *                     r, xcoords, ycoords, zcoords, dx,nx,ny,nz)
      implicit none
      integer  nx,ny,nz
      real(kind=8)   chi_rhs(nx,ny,nz),    chit_rhs(nx,ny,nz),
     *         chi(nx,ny,nz),        chit(nx,ny,nz),
     *         CHR(nx,ny,nz),        r(nx,ny,nz),
     *         xcoords(nx),          ycoords(ny),
     *         zcoords(nz),          dx
      include 'chr.inc'
      real(kind=8)   chit_x,  chit_y, chit_z
      integer  i,j,k

      logical     ltrace
      parameter ( ltrace = .false. )

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (CHR(i,j,k) .eq. CHR_amr_bdy .or.
     *       CHR(i,j,k) .eq. CHR_deco_bdy    ) then
            if (ltrace) write(*,*) 'Amr/Deco bdy: ',i,j,k
            !
            ! Artificial/finegrid boundaries:
            !
            chi_rhs(i,j,k)  =  0.d0
            chit_rhs(i,j,k) =  0.d0
         else if (CHR(i,j,k) .eq. CHR_interior) then
c           if (ltrace) write(*,*) 'Interior:     ',i,j,k
            !
            ! Interior points:
            !
            chi_rhs(i,j,k)  =  chit(i,j,k)
            chit_rhs(i,j,k) = 
     *              ( chi(i+1,j,k)-2.d0*chi(i,j,k)+chi(i-1,j,k) )/dx**2
     *            + ( chi(i,j+1,k)-2.d0*chi(i,j,k)+chi(i,j-1,k) )/dx**2
     *            + ( chi(i,j,k+1)-2.d0*chi(i,j,k)+chi(i,j,k-1) )/dx**2
         else
c           if (ltrace) write(*,*) 'Exterior bdy: ',i,j,k
            !
            ! Coarse grid boundary ("physical boundary):
            !
            chi_rhs(i,j,k)  =  chit(i,j,k)
            if (r(i,j,k).eq.0) then
               write(*,*) 'calc_rhs: r cannot be zero on boundary!'
               call my_exit('r value zero on boundary')
            end if
            !
            ! Take care of forwards and backwards differencing at bdy's
            !
            if (i.eq.1) then
               chit_x = ( - 3.d0 * chit(i,  j,k) 
     *                    + 4.d0 * chit(i+1,j,k)
     *                    -        chit(i+2,j,k)  ) / (2.d0*dx)
            else if (i.eq.nx) then
               chit_x = ( + 3.d0 * chit(i,  j,k) 
     *                    - 4.d0 * chit(i-1,j,k)
     *                    +        chit(i-2,j,k)  ) / (2.d0*dx)
            end if
            if (j.eq.1) then
               chit_y = ( - 3.d0 * chit(i,j  ,k) 
     *                    + 4.d0 * chit(i,j+1,k)
     *                    -        chit(i,j+2,k)  ) / (2.d0*dx)
            else if (j.eq.ny) then
               chit_y = ( + 3.d0 * chit(i,j  ,k) 
     *                    - 4.d0 * chit(i,j-1,k)
     *                    +        chit(i,j-2,k)  ) / (2.d0*dx)
            end if
            if (k.eq.1) then
               chit_z = ( - 3.d0 * chit(i,j,k  ) 
     *                    + 4.d0 * chit(i,j,k+1)
     *                    -        chit(i,j,k+2)  ) / (2.d0*dx)
            else if (k.eq.nz) then
               chit_z = ( + 3.d0 * chit(i,j,k  ) 
     *                    - 4.d0 * chit(i,j,k-1)
     *                    +        chit(i,j,k-2)  ) / (2.d0*dx)
            end if
            chit_rhs(i,j,k) =  
     *            - ( xcoords(i) / r(i,j,k) ) * chit_x
     *            - ( ycoords(j) / r(i,j,k) ) * chit_y
     *            - ( zcoords(k) / r(i,j,k) ) * chit_z
     *            -  chit(i,j,k) / r(i,j,k)
            if (ltrace) write(*,*) 'chit_rhs: ',i,j,k,chit_rhs(i,j,k)
            if (ltrace) write(*,*) 'chit_rhs: x,y,z: ',
     *                          xcoords(i),ycoords(j),zcoords(k)
            if (ltrace) write(*,*) 'chit_rhs: derivs:',chit_x,
     *                          chit_y,chit_z
         end if
      end do
      end do
      end do

      return
      end      ! END: calc_rhs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  rk_advance:                                                               cc
cc               Advance a field according to a generalized Runge_kutte method.c
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rk_advance(field_n,field_np1,field_rhs,CHR,
     *                      factor,nx,ny,nz)
      implicit none
      integer  nx,ny,nz
      real(kind=8)   field_n(nx,ny,nz),   field_np1(nx,ny,nz),
     *         field_rhs(nx,ny,nz), CHR(nx,ny,nz),
     *         factor
      include 'chr.inc'
      integer  i,j,k

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if ( CHR(i,j,k) .ne. CHR_amr_bdy .and.
     *        CHR(i,j,k) .ne. CHR_deco_bdy      ) then
            field_np1(i,j,k) = field_n(i,j,k) + factor* field_rhs(i,j,k)
         end if
      end do
      end do
      end do

      return
      end      ! END: rk_advance
