	MODULE MPSI4

	CONTAINS

	subroutine frame(l, n, m_R, m_I, x, y, z, alpha, beta, g, g_uu, detg)
	implicit none
	
	real*8 :: x, y, z, g(4,4), alpha, beta(3)
	real*8 :: l(4), n(4), m_R(4), m_I(4), g_uu(4,4)
	real*8 :: u(4), rr(4), t(4), f(4), sq2
	
!local quantities for intermediate steps
	real*8 :: v(3,3),tw1,tw2,tw3, vn(3,3), fud
	real*8 :: eps(3,3,3),r, temp, detg,w(3,3)
	real*8 :: t2,t4,t7,t9,t12,t15,t20,t24,t30
	integer :: li,lj,lk,ld
	logical, parameter :: ltrace =.false.
	
!define the total antysmmetric tensor
       	eps = 0.0d0
	eps(1,2,3) = 1.0d0; eps(2,3,1) = 1.0d0; eps(3,1,2) = 1.0d0
	eps(1,3,2) = -1.0d0; eps(2,1,3) = -1.0d0; eps(3,2,1) = -1.0d0
	
	r = sqrt(x**2+y**2+z**2)
	fud = 0.0
        if( (abs(y).lt.1.e-5) .and. (abs(x).lt.1.e-5) ) then
	  if (ltrace) print*, "we are on the z axis; tetrad set to zero"
	  l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
	else

	v(1,1) = -y; v(1,2) = x; v(1,3) = 0.0d0
	v(2,1) = x; v(2,2) = y; v(2,3) = z
	
	do li=1,3	
	   temp=0.0
	   do lj=1,3
	     do lk=1,3
	       do ld=1,3
	         temp =temp + sqrt(detg)*g_uu(li,lj)*eps(lj,lk,ld)*v(1,lk)*v(2,ld)
	       end do
	     end do
	   end do	
	   v(3,li)=temp
	end do
		  
	  
!now redefine things
! vector 1

!now get the inner product to orthonormalize wrt the metric

        temp = 0.0d0
	do lk=1,3
	  do ld=1,3
	     temp = temp + g(lk,ld)*v(1,lk)*v(1,ld)	
	  end do
	end do
	w(1,1) = temp 
	      
	if(w(1,1).lt.1.e-10) then
          l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
	  if (ltrace) print*, 'trouble with w11 in frame',tw1,v(1,1),v(1,2),x,y,z
	  GOTO 100	  
!	  STOP
	  tw1 = 1.0
	else
          tw1 = sqrt(w(1,1))
	  tw1 = 1.0d0/tw1
	end if
	  
	v(1,:) = v(1,:) * tw1
	
! vector 2
	
!now get the inner product to orthonormalize wrt the metric
	do li=1,2
	  temp = 0.0d0
	  do lk=1,3
	     do ld=1,3
	       temp = temp + g(lk,ld)*v(li,lk)*v(2,ld)	
	     end do
	  end do
	  w(li,2) = temp
	end do
	
	if(w(2,2).lt.1.e-10) then
          l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
	  if (ltrace) print*, 'trouble with w22 in frame',tw2,v(2,1),v(2,2)
	  GOTO 100
	  tw2 = 1.0
	else
          tw2 = sqrt(w(2,2))
	  tw2 = 1.0d0/tw2		  
	end if
	
	v(2,:) = (v(2,:) - v(1,:)*w(1,2)) * tw2
	
! vector 3	
!now get the inner product to orthonormalize wrt the metric
	do li=1,3
	  temp = 0.0d0
	  do lk=1,3
	    do ld=1,3
	      temp = temp + g(lk,ld)*v(li,lk)*v(3,ld)	
	    end do
	  end do
	  w(li,3) = temp	
	end do

	if(w(3,3).lt.1.e-10) then
          l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
	  if (ltrace) print*, 'trouble with w33 in frame', x,y,z,w(3,3)
	  GOTO 100	  
	  tw3 = 1.0
	else
          tw3 = sqrt(w(3,3))
	  tw3 = 1.0d0/tw3		  
	end if	
	
	v(3,:) = (v(3,:) - v(1,:)*w(1,3) - v(2,:)*w(2,3)) * tw3
	
	
!wit these constrcut the tetrad	
	rr(1) = v(2,1); rr(2) = v(2,2); rr(3) = v(2,3); rr(4) = 0.0
	t(1) = v(3,1); t(2) = v(3,2); t(3) = v(3,3); t(4) = 0.0
	f(1) = v(1,1); f(2) = v(1,2); f(3) = v(1,3); f(4) = 0.0
	u(1) = -1./alpha*beta(1);u(2)=-1./alpha*beta(2)
	u(3) = -1./alpha*beta(3); u(4) = 1./alpha

	sq2 = 1./sqrt(2.0d0)
	
	do li=1,4
	  l(li) = sq2 * (u(li) + rr(li))
	  n(li) = sq2 * (u(li) - rr(li))
	  m_r(li) = sq2 * t(li)
	  m_i(li) = sq2 * f(li) 
	end do

100 	CONTINUE

	end if
	
	end subroutine frame	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine rieman(R, Chd, guu, g, dg, d_dg, dime)
	implicit none
	
	real*8, dimension(:,:,:,:) :: R
	real*8, dimension(:,:) :: g, guu
	real*8, dimension(:,:,:) :: dg, Chd
	real*8, dimension(:,:,:,:) :: d_dg
	integer :: dime
	
        !local vars
	real*8, dimension(:,:,:),allocatable::Chu

#       include "temp_varsrieman.inc"
        
        allocate(Chu(4,4,4))

        ! compute the Christophel with the index down Chd and up Chp
#        include "Ch1.inc"

#       include "R4down1.inc"

        deallocate(Chu)

        end subroutine rieman


	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FROM 3D to 4D Rieman
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine rieman3to4(R4, R3, Ch, g_uu, g, K, D_K)
	implicit none
	
	real*8, dimension(:,:,:,:) :: R4, R3
	real*8, dimension(:,:,:) :: Ch, D_K
	real*8, dimension(:,:) :: g, K, g_uu
	
	integer :: lj,lk,ll,li
	real*8,dimension(3,3) :: Ric3
	real*8:: temp, trK
	
!get Ricci
	do li=1,3
	do lj=1,3
	
	 temp = 0.0d0
	  do lk=1,3
	  do ll=1,3	
	  temp = temp + R3(li,lk,lj,ll)*g_uu(lk,ll)
	  end do 
	  end do
	
	 Ric3(li,lj) = temp
	  
	end do
	end do
	
	
!get traceK
	trK = 0.0d0
	
	do li=1,3
	 do lj=1,3	
	  trK = trK + K(li,lj)*g_uu(li,lj)
	 end do
	end do
	
	
!fill spatial comps for 4 Rieman
	do li=1,3
	do lj=1,3
	do lk=1,3
	do ll=1,3
	R4(li,lj,lk,ll) = R3(li,lj,lk,ll) + &
	&                 K(li,lk)*K(ll,lj)-K(li,ll)*K(lk,lj)
	
	end do
	end do
	end do
	end do

!fill time-spa-time-spa for 4 Rieman	
	
	do lj=1,3
	do ll=1,3
	  temp = 0.0d0
	  do lk = 1,3
	  do li = 1,3
	    temp = temp + K(lj,li)*K(ll,lk)*g_uu(li,lk)
	  end do
	  end do	
	
	R4(4,lj,4,ll) = Ric3(lj,ll) - temp + trK * K(lj,ll)
	
	end do
	end do

!fill time-spa-spa-spa 4Riem
	do lj=1,3
	do lk=1,3
	do ll=1,3
	  temp = 0.0d0
	   do li=1,3
	     temp = temp + Ch(li,lj,lk)*K(ll,li)-Ch(li,lj,ll)*K(lk,li)
	   end do	
	   
	 R4(4,lj,lk,ll) = - ( D_K(ll,lj,lk)-D_K(lk,lj,ll) + temp )  
	end do
	end do
	end do
	
!apply symmetry	
        do lj=1,3
        do lk=1,3
        do ll=1,3
           R4(ll,4,lj,ll) = -R4(4,ll,lj,ll)
           R4(ll,lj,4,lk) =  R4(4,lk,ll,lj)
           R4(ll,lj,lk,4) = -R4(4,lk,ll,lj)
        end do
           R4(lj,4,lk,4) =  R4(4,lj,4,lk)
	   R4(lj,4,4,lk) = -R4(4,lj,4,lk)
	   R4(4,lj,lk,4) = -R4(4,lj,4,lk)
        end do
        end do
	
	return
	end subroutine rieman3to4
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	subroutine calc_psi4(PSI4_R, PSI4_I, R4, n, m_R, m_I)
	implicit none
	real*8, dimension(:,:,:,:) :: R4
	real*8, dimension(:) :: n, m_R, m_I
	real*8 :: PSI4_R, PSI4_I	
	
        !local vars
	integer:: la,lb,lc,ld
	real*8 :: PP(2)

#       include "temp_psi4.inc"
	
	PSI4_R = 0.0d0
	PSI4_I = 0.0d0

#       include "psi4.inc"

        PSI4_R = PP(1)
        PSI4_I = PP(2)

!	do la=1,4
!	do lb=1,4
!	do lc=1,4
!	do ld=1,4

!	 PSI4_R = PSI4_R + R4(la,lb,lc,ld)* &
!	&     ( n(la)*m_R(lb)*n(lc)*m_R(ld) - n(la)*m_I(lb)*n(lc)*m_I(ld) )
	
!	 PSI4_I = PSI4_I - R4(la,lb,lc,ld)* &
!	&     ( n(la)*m_I(lb)*n(lc)*m_R(ld) + n(la)*m_R(lb)*n(lc)*m_I(ld) )
		
!	end do
!	end do
!	end do
!	end do

!	print*, "psi4R"
!        print*, PSI4_R
!        print*, PSI4_I	
!	stop

	
	return
	
	end subroutine calc_psi4

	

      !--------------------------------------------------------------------
      !
      !   calculate the 4-curvature
      !
      !--------------------------------------------------------------------
      subroutine cal_curvature(curvature, R4, g_uu)
      implicit   none

      real*8                     :: curvature
      real*8, dimension(4,4,4,4) ::  R4
      real*8, dimension(4,4)     :: g_uu

      ! local vars
      integer :: la, lb, lc, ld, le, lf, lg, lh
      real*8 CC(1)

#     include "temp_curv1.inc"   

      curvature = 0.0d0

!      do la=1,4; do lb=1,4; do lc=1,4; do ld=1,4
!         do le=1,4; do lf=1,4; do lg=1,4; do lh=1,4
!             curvature = curvature + &
!     &                  (R(la,lb,lc,ld)*R(le,lf,lg,lh))*g_uu(la,le) &
!     &                     *g_uu(lb,lf)*g_uu(lc,lg)*g_uu(ld,lh)
!           end do; end do; end do; end do
!         end do; end do; end do; end do

#       include "curv1.inc"
        curvature = CC(1)


      return
      end subroutine cal_curvature

	END MODULE MPSI4
