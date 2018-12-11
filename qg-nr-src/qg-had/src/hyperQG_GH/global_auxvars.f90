!-----------------------------------------------------------------
!
!  $Id: global_auxvars.f90,v 1.4 2007-04-03 21:20:53 carlos Exp $
!
!-----------------------------------------------------------------

#include "cctk.h"
	 
MODULE HYPER_AUXVARS
    use params
    use GF
    use auxvars
  contains

!-----------------------------------------------------------------
!
! SUBROUTINE hyperauxvars
!
!-----------------------------------------------------------------
  subroutine hyperauxvars(v, w, u, par)

    implicit none

    type(gridfunction), dimension(NU)   :: u
    type(gridfunction), dimension(NV)   :: v
    type(gridfunction), dimension(NW)   :: w
    CCTK_REAL, dimension(:)             :: par
    
    CCTK_REAL, dimension(NU)            :: u_pt
    CCTK_REAL, dimension(NV)            :: v_pt
    CCTK_REAL, dimension(NW)            :: w_pt
    CCTK_INT                            :: ii,jj,kk, nx, ny, nz, m
    CCTK_REAL                           :: myl2norm3d
    CCTK_INT			        :: use_mask
    CCTK_REAL, pointer, dimension(:,:,:):: mask
    logical, parameter                  :: ltrace2 = .false.

    mask => w(H_MASK)%d
    nx = nint(par(P_NX))
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))
    use_mask = nint(par(P_USE_MASK))
   
!   add a parameter here to compute auxvars or not    
    if (ltrace2) write(0,*)'HYPERAUXVARS: calling define_aux_vars'
             
    do kk = 1, nz
       do jj = 1, ny
          do ii = 1, nx

           if ( (use_mask .EQ. 0) .OR. & 
	       ((use_mask .EQ. 1).AND.(mask(ii,jj,kk) .GT. 0.0)) ) then
             do m = 1, NU
                u_pt(m) = u(m)%d(ii,jj,kk)
             end do
             do m = 1, P_NV_REAL_SIZE
                v_pt(m) = v(m)%d(ii,jj,kk)
             end do
             do m = 1, NW
                w_pt(m) = w(m)%d(ii,jj,kk)
             end do
			      
             call define_aux_vars(v_pt, w_pt, u_pt, par)	     

             do m=1, P_NV_REAL_SIZE
                v(m)%d(ii,jj,kk) = v_pt(m)
             end do
           end if
	   
          end do
       end do
    end do
	      
    if (ltrace2) write(0,*)'HYPERAUXVARS: return'

    return
  end subroutine hyperauxvars
end MODULE HYPER_AUXVARS
