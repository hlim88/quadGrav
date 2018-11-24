c----------------------------------------------------------------------
c    sphere.f
c    Contains routines to construct an initial surface 
c    Modified 5-21-93    MFH
c        - Switched to Spherical coordinates
c        - S(theta,phi,3) -> 1 -> radial, 2 -> theta, 3-> phi
c        - Indices for theta-phi mesh on surface
c    Modified 6-30-93    FH
c        - Switched to non-staggered grid
c    Modified Sat Sep 18 18:14:27 CDT 1993 MFH
c        - Center of coordinate system center(3)
c        - Offset of surface from origin x0,y0,z0
c----------------------------------------------------------------------

      subroutine GenSur(S,Sc,Sf,NDIM,NST,NSP)
        implicit none
        include 'horizon_parameters.h'

        integer NDIM
        integer NST
        integer NSP
        integer np,nt
        
        real*8 S(NST,NSP,NDIM)
        real*8 Sc(NST,NSP,NDIM)
        real*8 Sf(NST,NSP)



c    Local variables

        integer i,j,k
        
        real*8 x,y,z
        real*8 t,p
        real*8 temp

c        write(6,*)'Radius ',radius
c        write(6,*)'Number ',NSS*(NSS+1)
c        write(0,*)x0,y0,z0,tau
c        write(0,*)'CENTER',center(1),center(2),center(3)
c        write(2,*)'CENTER',center(1),center(2),center(3)

c    Define PI
        Pi = 4.0e0*atan(1.0e0)

c    Do loop to calculate points on a sphere

c    Store values as r,theta,phi
        do i=1,NST
            t = float(i-1)*dtheta
            do j=1,NSP
                 p = float(j-1)*dphi

                 x = radius*cos(p)*sin(t) + x0 + center(1)
                 y = radius*sin(p)*sin(t) + y0 + center(2)
                 z = radius*cos(t)        + z0 + center(3)

c                 Sf(i,j)  = 0.0e0
c                 Sf(i,j)  = sqrt(x**2 + y**2 + z**2)
c                 Sf(i,j)  = radius
c                 S(i,j,1) = radius

                 Sf(i,j)  = sqrt((x-center(1))**2 + (y-center(2))**2 +
     #                                 (z-center(3))**2)
c            Spherical coordinates
c                 S(i,j,1)  = sqrt((x-center(1))**2 + (y-center(2))**2 +
c     #                                 (z-center(3))**2)
                 S(i,j,1) = Sf(i,j)
                 S(i,j,2) = t
                 S(i,j,3) = p

c            Cartesian Coordinates
                 Sc(i,j,1) = x
                 Sc(i,j,2) = y
                 Sc(i,j,3) = z

                
c                 write(66,*)S(i,j,1),S(i,j,2),S(i,j,3)
c                 write(67,*)Sc(i,j,1),Sc(i,j,2),Sc(i,j,3)
c                 write(68,*)x,y,z
c                 write(3,*)i,j,Sf(i,j)
            enddo
c            write(3,*)
c            write(68,*)
        enddo
        return
        end
c----------------------------------------------------------------------
        

        
        
    

