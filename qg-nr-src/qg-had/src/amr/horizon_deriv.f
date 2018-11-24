C-----------------------------------------------------------------------
c    deriv.f
c
c    created 4-14-93    Mijan Huq
c
c-----------------------------------------------------------------------
        subroutine Get1DDeriv(GrF,h,NG,Result)
            implicit none

            integer    NG
            
            real*8        GrF(NG)
            real*8        h
            real*8        Result(NG)

c        Local Variables

            integer     i
            
            real*8        diff

c        Use a forward differencing scheme until we reach the end

            do i=1,NG
                
                if(i .lt. NG .and. i .gt. 1) then
                    diff = 0.5e0*(GrF(i+1) - GrF(i-1))
                    Result(i) = diff/h
                endif
                if(i .eq. NG) then
                    diff = 2.0e0*GrF(i-1) - 3.0e0*GrF(i-2) + GrF(i)
                    Result(i) = 0.25e0*diff/h
                endif
                if(i .eq. 1) then
                    diff = 2.0e0*GrF(i+1) - 3.0e0*GrF(i) + GrF(i+2)
                    Result(i) = 0.25e0*diff/h
                endif
            
            enddo

            end
c-----------------------------------------------------------------------    

        subroutine ddx(f,h,df,NP)
            implicit none

            integer    NP

            real*8    f(NP)
            real*8    h
            real*8    df


            df = 0.5e0*(f(3) - f(1))/h

        return
        end
c-----------------------------------------------------------------------
        subroutine Deriv1D2nd(f,d2fdx2,h,NDIM)
            implicit none

            integer    NDIM

            real*8    f(NDIM)
            real*8    d2fdx2
            real*8    h

            d2fdx2=(f(3) - 2*f(2) + f(1))/h**2    

            return

            end
c-----------------------------------------------------------------------
c        Fri Sep 10 18:38:59 CDT 1993
c        The following routine takes a function defined on a 3D cartesian
c        grid and does finite differencing to second order.
c        Note that this routine returns zeros for the boundary points
c        so as to avoid all the boundary headaches with endpoint
c        differencing.
c

        subroutine GrSc3D2D(GrF,GrDF,hx,hy,hz,NX,NY,NZ)
            implicit none

            integer    NX,NY,NZ

            real*8    hx,hy,hz
            real*8    GrF(NX,NY,NZ)
            real*8    GrDF(NX,NY,NZ,3)

c        Local Variables

            integer    i,j,k

            real*8    divhx    
            real*8    divhy    
            real*8    divhz    

c        Calculate half of inverse h[x-z] ahead of time.

            divhx = 0.5e0/hx
            divhy = 0.5e0/hy
            divhz = 0.5e0/hz

c        Zero out the derivative array

            do i=1,NX
                do j=1,NY
                    do k=1,NZ

                        GrDF(i,j,k,3) = 0.0e0
                    
                    enddo
                enddo
            enddo

c        Loop through the X direction skipping the boundaries
        
            do i=2,NX-1

c            Loop through the Y direction skipping the boundaries
        
                do j=2,NY-1

c                Loop through the Z direction skipping the boundaries
        
                    do k=2,NZ-1

c                    Centered Differencing

                     GrDF(i,j,k,1) = divhx*(GrF(i+1,j,k) - GrF(i-1,j,k))
                     GrDF(i,j,k,2) = divhy*(GrF(i,j+1,k) - GrF(i,j-1,k))
                     GrDF(i,j,k,3) = divhz*(GrF(i,j,k+1) - GrF(i,j,k-1))

c                    Done

                    enddo
                enddo
            enddo

            return
        end
c-----------------------------------------------------------------------
