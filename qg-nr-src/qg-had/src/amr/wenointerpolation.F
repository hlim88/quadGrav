#include "cctk.h"

c----------------------------------------------------------------------
c
c     function weno_interp_centered
c
c     This routine is called from interp_from_parentW, which is the had
c     routine that sets ups the interpolation stencil for WENO. 
c
c----------------------------------------------------------------------
      CCTK_REAL function weno_interp_centered(u)
        implicit none
        CCTK_REAL, dimension(6) :: u

        ! local vars
        CCTK_REAL               :: a, b, c, d, e
        CCTK_REAL, external     :: interpmiddle

        a = u(1)
        b = u(2)
        c = u(3)
        d = u(4)
        e = u(5)
      
        weno_interp_centered = interpmiddle(a,b,c,d,e)

      end function weno_interp_centered

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      CCTK_REAL function weno_interp_1d(xi,u,x)
        implicit none
        CCTK_REAL               :: xi
        CCTK_REAL, dimension(5) :: u, x

        ! local vars
        CCTK_REAL               :: ans, a, b, c, d, e
        CCTK_REAL, external     :: interpedge
        CCTK_REAL, external     :: interpmiddle

        a = u(1)
        b = u(2)
        c = u(3)
        d = u(4)
        e = u(5)
      
        ans = 0.0d0
        !if      (x(1) .lt. xi .and. xi .lt. x(2)) then
        !  ans = interpedge(e,d,c,b,a)
        !else 
        if (x(2) .lt. xi .and. xi .lt. x(3)) then
          ans = interpmiddle(e,d,c,b,a)
        else if (x(3) .lt. xi .and. xi .lt. x(4)) then
          ans = interpmiddle(a,b,c,d,e)
        !else if (x(4) .lt. xi .and. xi .lt. x(5)) then
        !  ans = interpedge(a,b,c,d,e)
        else
          write(0,*)'weno_interp_1d: invalid point'
          write(0,*)'weno_interp_1d: xi = ',xi
          write(0,*)'weno_interp_1d: x  = ',x
          write(0,*)'weno_interp_1d: u  = ',u
        end if

        weno_interp_1d = ans
      end function weno_interp_1d


c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      CCTK_REAL function weno_interp_1d_edge(xi,u,x)
        implicit none
        CCTK_REAL               :: xi
        CCTK_REAL, dimension(5) :: u, x

        ! local vars
        CCTK_REAL               :: ans, a, b, c, d, e
        CCTK_REAL, external     :: interpedge
        CCTK_REAL, external     :: interpmiddle

        a = u(1)
        b = u(2)
        c = u(3)
        d = u(4)
        e = u(5)

        ans = 0.0d0
        if      (x(1) .lt. xi .and. xi .lt. x(2)) then
          ans = interpedge(e,d,c,b,a)
        else if (x(4) .lt. xi .and. xi .lt. x(5)) then
          ans = interpedge(a,b,c,d,e)
        else
          write(0,*)'weno_interp_1d_edge: invalid point'
          write(0,*)'weno_interp_1d_edge: xi = ',xi
          write(0,*)'weno_interp_1d_edge: x  = ',x
          write(0,*)'weno_interp_1d_edge: u  = ',u
        end if

        weno_interp_1d_edge = ans
      end function weno_interp_1d_edge

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      CCTK_REAL function interpmiddle(a,b,c,d,e)
        implicit none
        CCTK_REAL, intent(IN) :: a,b,c,d,e
      
        CCTK_REAL :: ans, p(3), w(3)
        integer :: i
      
        ans = 0.0d0
      
        call weno_polynomials(a,b,c,d,e,p)
        call weno_weights(a,b,c,d,e,w)
        do i=1,3
          ans = ans + w(i)*p(i)
        end do
        interpmiddle = ans
      end function interpmiddle

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      CCTK_REAL function interpedge(a,b,c,d,e)
        implicit none
        CCTK_REAL, intent(IN) :: a,b,c,d,e
      
        CCTK_REAL  :: ans, p(3), w(3)
        integer :: i

        ans = 0.0d0

        call weno_edge_polynomials(a,b,c,d,e,p)
        call weno_edge_weights(a,b,c,d,e,w)

        if(w(3)>w(2)) then
          ans = p(3)
        else
          ans = p(2)
        end if

        interpedge = ans
      end function interpedge

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      subroutine weno_polynomials(a,b,c,d,e, p)
        implicit none
        CCTK_REAL, intent(IN)  :: a,b,c,d,e
        CCTK_REAL, intent(OUT) :: p(3)
      
        p(1) = (3.0d0*a-10.0d0*b+15.0d0*c)/8.0d0
        p(2) = (-b+6.0d0*c+3.0d0*d)/8.0d0
        p(3) = (3.0d0*c+6.0d0*d-e)/8.0d0
      end subroutine weno_polynomials
      
c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      subroutine weno_edge_polynomials(a,b,c,d,e, p)
        implicit none
        CCTK_REAL, intent(IN)  :: a,b,c,d,e
        CCTK_REAL, intent(OUT) :: p(3)
      
        p(1) = (15.0d0*a-42.0d0*b+35.0d0*c)/8.0d0
        p(2) = (3.0d0*b-10.0d0*c+15.0d0*d)/8.0d0
        p(3) = (-c+6.0d0*d+3.0d0*e)/8.0d0
      end subroutine weno_edge_polynomials

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      subroutine weno_smoothness(a,b,c,d,e, s)
        implicit none
        CCTK_REAL, intent(IN)  :: a,b,c,d,e
        CCTK_REAL, intent(OUT) :: s(3)

        s(3)=(10.0d0*c*c-31.0d0*c*d+25.0d0*d*d+11.0d0*c*e
     &                   -19.0d0*d*e+4.0d0*e*e)/3.0d0
        s(2)=(4.0d0*b*b-13.0d0*c*b+13.0d0*c*c+5.0d0*b*d
     &                   -13.0d0*c*d+4.0d0*d*d)/3.0d0
        s(1)=(4.0d0*a*a-19.0d0*a*b+25.0d0*b*b+11.0d0*c*a
     &                   -31.0d0*b*c+10.0d0*c*c)/3.0d0

      end subroutine weno_smoothness
      
c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      subroutine weno_mjb_smoothness(a,b,c,d,e,s)
        implicit none
        CCTK_REAL, intent(IN) :: a,b,c,d,e
        CCTK_REAL, intent(OUT):: s(3)
        CCTK_REAL, external   :: mjb_smoothness
      
        s(1)=mjb_smoothness(a,b,c)
        s(2)=mjb_smoothness(b,c,d)
        s(3)=mjb_smoothness(c,d,e)

      end subroutine weno_mjb_smoothness

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      CCTK_REAL function mjb_smoothness(a,b,c)
        implicit none
        CCTK_REAL, intent(IN) :: a,b,c
      
        mjb_smoothness = (c-2.0d0*b+a)**2+(a*a+2.0d0*b*b
     &                 + c*c-2.0d0*b*(a+c))/2.0d0
      end function mjb_smoothness

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      subroutine weno_weights(a,b,c,d,e,w)
        implicit none
        CCTK_REAL, intent(IN)  :: a,b,c,d,e
        CCTK_REAL, intent(OUT) :: w(3)
      
        CCTK_REAL              :: s(3), total
        CCTK_REAL, parameter   :: epsilon = 1.0d-5
        integer                :: i


        call weno_smoothness(a,b,c,d,e,s)
      
        ! I set the d values such that if all the smoothnesses are the same
        ! the rusult will be a interpolated value based on all 5 points
        w(1)=1.0d0/(16.0d0*((epsilon+s(1))**3))
        w(2)=10.0d0/(16.0d0*((epsilon+s(2))**3))
        w(3)=5.0d0/(16.0d0*((epsilon+s(3))**3))

        total=w(1)+w(2)+w(3)

        do i=1,3
          w(i)=w(i)/total
        end do
      end subroutine weno_weights

c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
      subroutine weno_edge_weights(a,b,c,d,e,w)
        implicit none
        CCTK_REAL, intent(IN)  :: a,b,c,d,e
        CCTK_REAL, intent(OUT) :: w(3)

        CCTK_REAL              :: s(3), total
        CCTK_REAL, parameter   :: epsilon = 1.0d-5
        integer                :: i

        call weno_smoothness(a,b,c,d,e,s)

        w(1)=1.0d0/((epsilon+s(1))**2)
        w(2)=5.0d0/((epsilon+s(2))**2)
        w(3)=10.0d0/((epsilon+s(3))**2)
      
        total=w(1)+w(2)+w(3)

        do i=1,3
          w(i)=w(i)/total
        end do
      end subroutine weno_edge_weights

