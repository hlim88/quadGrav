#include "cctk.h"

!--------------- R E A D   M E ------------------------------------------
! A more usueful analysis routine
!------------------------------------------------------------------------


  subroutine HyperAnalysis2(u0, v, w, par)
  use params
  use GF
  use UTILEQS
  use MOD_ANALYSIS2

  implicit none
  
  type(gridfunction), dimension(NU)    :: u0
  type(gridfunction), dimension(NV)    :: v
  type(gridfunction), dimension(NW)    :: w
  CCTK_REAL, dimension(NPAR)           :: par
  CCTK_INT, dimension(:,:,:), allocatable :: imask
  logical, parameter                   :: ltrace = .false.
  

  call analysis2(u0, v, w, par)


  if (ltrace) then
     write(0,*)'>>> HYPER_ANALYSIS2: return from analysis'
  end if
  
  return
end subroutine HyperAnalysis2
    
  
  
