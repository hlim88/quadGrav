        MODULE GF
        implicit none
        type gridfunction
        real(kind=8), dimension(:,:,:), pointer :: d
        integer :: take_dx, take_dy, take_dz, stiff
        integer :: dissipation, zsym, ysym, xsym
        real(kind=8) :: excise_val, diss_factor
        end type gridfunction
        end MODULE GF
