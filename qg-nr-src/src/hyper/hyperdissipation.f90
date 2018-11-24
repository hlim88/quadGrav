
#include "cctk.h"

MODULE HYPER_DISSIPATION

  use params
  use GF
  use HYPER_DIS_21
 
contains

!-----------------------------------------------------------------
!
! SUBROUTINE HYPERDISSIPATION
!
!-----------------------------------------------------------------
  subroutine hyperdissipation(dxu, dyu, dzu, u, w, imask, par)
    implicit none
 
    CCTK_REAL, dimension(:,:,:)         :: dxu, dyu, dzu, u
    type(gridfunction), dimension(NW)   :: w
    CCTK_INT,          dimension(:,:,:) :: imask
    CCTK_REAL, dimension(NPAR)             :: par

    CCTK_INT                            :: direc, use_mask, dissipation
    CCTK_INT                            :: order
    logical, parameter                  :: ltrace = .false.
    integer :: myid, proc_return_myid

    if (ltrace) myid = proc_return_myid()

   if(ltrace) write(*,*) myid,']  entering hyperdissipation'

    use_mask    = nint(par(P_USE_MASK   ))
    order       = nint(par(P_DERIV_ORDER))
    dissipation = nint(par(P_DISSIPATION))


    if (      (order .le.   4)  &  
         .or. (order .eq.  42)  & 
         .or. (order .eq.  44)  & 
         .or. (order .eq. 642)  & 
         .or. (order .eq. 666)     ) then

      if (use_mask .eq. 1) then
        !---------------------------------------
        ! MASKED DERIVATIVES
        !---------------------------------------
        if      (dissipation .eq. 1) then
           direc = 1
           call Dis21(dxu,u,direc,imask,par)
  
           direc = 2
           call Dis21(dyu,u,direc,imask,par)
  
           direc = 3
           call Dis21(dzu,u,direc,imask,par)
  
        else if (dissipation .eq. 2) then
           direc = 1
           call KODis21(dxu,u,direc,imask,par)
  
           direc = 2
           call KODis21(dyu,u,direc,imask,par)
   
           direc = 3
           call KODis21(dzu,u,direc,imask,par)
   
        else if (dissipation .eq. 4) then
           direc = 1
           call KODis42(dxu,u,direc,imask,par)
   
           direc = 2
           call KODis42(dyu,u,direc,imask,par)
   
           direc = 3
           call KODis42(dzu,u,direc,imask,par)
 
        else
           write(0,*)'masked derivatives not implemented for dissipation '&
                     ,dissipation
        end if
  
      else
        !---------------------------------------
        ! UNMASKED DERIVATIVES
        !---------------------------------------
        if (dissipation .eq. 1) then
           direc = 1
           call Dis21NoMask(dxu,u,direc,par)
  
           direc = 2
           call Dis21NoMask(dyu,u,direc,par)
  
           direc = 3
           call Dis21NoMask(dzu,u,direc,par)
  
        else if (dissipation .eq. 2) then
           direc = 1
           call KODis21NoMask(dxu,u,direc,par)
  
           direc = 2
           call KODis21NoMask(dyu,u,direc,par)
  
           direc = 3
           call KODis21NoMask(dzu,u,direc,par)
  
        else if (dissipation .eq. 4) then
           direc = 1
           call KODis42NoMask(dxu,u,direc,par)
  
           direc = 2
           call KODis42NoMask(dyu,u,direc,par)
  
           direc = 3
           call KODis42NoMask(dzu,u,direc,par)
  
        else if (dissipation .eq. 42 .or. dissipation .eq. 44) then
           direc = 1
           call KODis42NoMask(dxu,u,direc,par)
  
           direc = 2
           call KODis42NoMask(dyu,u,direc,par)
  
           direc = 3
           call KODis42NoMask(dzu,u,direc,par)
  
        else if (dissipation .eq. 642 .or. dissipation .eq. 666) then
           direc = 1
           call KODis642NoMask(dxu,u,direc,par)
  
           direc = 2
           call KODis642NoMask(dyu,u,direc,par)
  
           direc = 3
           call KODis642NoMask(dzu,u,direc,par)
  
        else 
           write(0,*)'unmasked derivatives not implemented for dissipation '&
                    ,dissipation
        end if
  
      end if
      else
         write(0,*)'Currently only 2nd order derivatives are available ',order
         write(0,*)'Well... not really.  Someone should really figure '
         write(0,*)'this out.  '
  
    end if
   if(ltrace) write(*,*) myid,']  leaving hyperdissipation'
  
    return
  end subroutine hyperdissipation

end MODULE HYPER_DISSIPATION
