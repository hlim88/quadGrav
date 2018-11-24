c-----------------------------------------------------------------------
c    Interpolation Routines
c

c=======================================================================
      real*8 function L4irr(ft,x,tx)
         implicit none

         real*8   ft(4)
         real*8   x(4)
         real*8   tx


         L4irr = (tx - x(2))*(tx - x(3))*(tx - x(4))*ft(1)/
     .            ((x(1) - x(2))*(x(1) - x(3))*(x(1) - x(4))) +

     .            (tx - x(1))*(tx - x(3))*(tx - x(4))*ft(2)/
     .             ((x(2) - x(1))*(x(2) - x(3))*(x(2) - x(4)))+

     .            (tx - x(1))*(tx - x(2))*(tx - x(4))*ft(3)/
     .            ((x(3) - x(1))*(x(3) - x(2))*(x(3) - x(4)))+

     .            (tx - x(1))*(tx - x(2))*(tx - x(3))*ft(4)/
     .            ((x(4) - x(1))*(x(4) - x(2))*(x(4) - x(3)))
         return
      end


c=======================================================================
        real*8 function L3intp(f,x,xi)

            implicit none
        
            real*8    f(3), x(3), xi, fi

c            Carry out lagrange interpolation

            fi = (xi - x(2))*(xi - x(3))*f(1)/((x(1)-x(2))*(x(1)-x(3)))
     *     +  (xi - x(1))*(xi - x(3))*f(2)/((x(2)-x(1))*(x(2)-x(3)))
     *     +  (xi - x(1))*(xi - x(2))*f(3)/((x(3)-x(1))*(x(3)-x(2)))

            L3intp = fi

            return
        end

c=======================================================================
        real*8 function L4intp(f,x,xi)

            implicit none
        
            real*8    f(4), x(4), xi, fi

c            Carry out lagrange interpolation

            fi=  (xi-x(2))*(xi-x(3))*(xi-x(4))*f(1)
     *        /((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4)))
     *     +  (xi - x(1))*(xi - x(3))*(xi-x(4))*f(2)
     *        /((x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4)))
     *     +  (xi - x(1))*(xi - x(2))*(xi-x(4))*f(3)
     *        /((x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4)))
     *     +  (xi - x(1))*(xi - x(2))*(xi-x(3))*f(4)
     *        /((x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3)))

            L4intp = fi

            return
        end
c=======================================================================
c        4th order interpolator
      real*8 function L4interp(ft,x,h,tx)
         implicit none

         real*8   ft(4)
         real*8   x(4)
         real*8   h
         real*8   tx

c        write(20,*)'P ',tx
c        write(20,*)'X ',x(1),x(2),x(3),x(4)

         L4interp = -(tx-x(2))*(tx - x(3))*(tx -x(4))*
     #      ft(1)/6.0d0 +
     #   (tx-x(1))*(tx - x(3))*(tx -x(4))*ft(2)/2.0d0-
     #   (tx-x(1))*(tx - x(2))*(tx -x(4))*ft(3)/2.0d0+
     #   (tx-x(1))*(tx - x(2))*(tx -x(3))*ft(4)/6.0d0
         L4interp = L4interp/(h**3)

         return
      end
c=======================================================================
    
c    3D Interpolation routine.  Created Thu Mar 17 15:39:55 CST 1994 MFH
c        3rd Order


        subroutine DINT3D3(ftmpl,cx,cy,cz,xi,fintp,hx,hy,hz)

            implicit none

c            Function value cube
            real*8    ftmpl(3,3,3)

c            Coordinate value cube
            real*8    cx(3), cy(3), cz(3)

c            Point to interpolate to...
            real*8    xi(3)

c            Return value for interpolated value for function
            real*8    fintp

c            Spacings...
            real*8    hx, hy, hz

c        Local Variables

            integer    i, j, k,
     *            ix, iy, iz

            real*8   fx(3), fy(3), fz(3)

c        Function prototypes
            
            real*8    L3intp


c            Loop over the z-direction
            do k=1,3

c                Loop over the y-direction
                do j=1,3

c                    Load in the function values in the x-direction
                    do i=1,3
                        fx(i) = ftmpl(i,j,k)
                    enddo

c                    Interpolate in the x-direction
                    fy(j) = L3intp(fx,cx,xi(1))
                
                end do

c                Now interpolate in the y-direction
                fz(k) = L3intp(fy,cy,xi(2))    

            end do

c            Now interpolate in the z-direction
            fintp = L3intp(fz,cz,xi(3)) 
                    
            return
        end
c=======================================================================
c    Interpolate Scalar function Thu Mar 17 17:03:50 CST 1994 MFH
c    3rd Order...
c
c    Routine returns the value at point ix
c
c    GrF       := Grid function
c    Gr[X-Z] := Coordinates on grid
c    ix          := Point at which interpolation is to be done
c    h[x-z]  := grid spacings
c    N[X-Z]  := Grid sizes
c

        real*8 function Intp3Sc3(GrF,GrX,GrY,GrZ,ix,hx,hy,hz,NX,NY,NZ)
    
            implicit none

            integer    NX, NY, NZ

            real*8    GrF( NX, NY, NZ),
     *               GrX(NX),  GrY(NY),  GrZ(NZ),
     *            ix(3),
     *                hx, hy, hz
        

c        Local Variables

            integer    i, j, k, xint, yint, zint

            real*8    ftmpl(3,3,3),
     *                cx(3), cy(3), cz(3)

c        Functions


c        Search for neighboring points for the interpolation.
c        Note that interpolation at the edges is not carried out here.

c            Calculate integer coordinates for interpolation point

            xint = int((ix(1)  - GrX(1))/hx) + 1
            yint = int((ix(2)  - GrY(1))/hy) + 1
            zint = int((ix(3)  - GrZ(1))/hz) + 1

c            Create the interpolation template

            do i=-1,1
                do j=-1,1
                    do k=-1,1

                        ftmpl(i+2,j+2,k+2) = GrF(xint+i, yint+j, zint+k)

                    end do
                end do
            end do
            
c            Coordinates...

            do i=-1,1

                cx(i+2) = GrX(xint+i)    
                cy(i+2) = GrY(yint+i)    
                cz(i+2) = GrZ(zint+i)    

            end do

c            Interpolate
            call DINT3D3(ftmpl,cx,cy,cz,ix,Intp3Sc3,hx,hy,hz)

            return
        end
c=======================================================================
c
c    Tensor Interpolation routine.
c    Interpolates rank 2 tensors to 3rd order
c    Arguments
c    GrT         := Grid Tensor
c    Gr[X-Z]    := Coordinates on grid
c    ix            := Interpolation point
c    it            := Return value for tensor
c    N[C,R]    := Number of components
c    N[X-Z]    := Number of grid points
c
        subroutine Int2Tsr3(GrT,GrX,GrY,GrZ,ix,it,hx,hy,hz,
     #                    NR,NC,NX,NY,NZ)
            implicit none

            integer    NR, NC,
     *                NX, NY, NZ

            real*8    GrT(NX, NY, NZ, NR, NC),
     *                GrX(NX), GrY(NY), GrZ(NZ),
     *                ix(3),    it(NR,NC),
     *                hx, hy, hz

c        Local Variables

            integer i, j, k, xint, yint, zint,
     *                p, q

            real*8   ftmpl(3,3,3),
     *            cx(3), cy(3), cz(3)

c     Functions


c     Search for neighboring points for the interpolation.
c     Note that interpolation at the edges is not carried out here.

c        Calculate integer coordinates for interpolation point

         xint = int((ix(1)  - GrX(1))/hx) + 1
         yint = int((ix(2)  - GrY(1))/hy) + 1
         zint = int((ix(3)  - GrZ(1))/hz) + 1

c            Coordinates...

            do i=-1,1

                cx(i+2) = GrX(xint+i)    
                cy(i+2) = GrY(yint+i)    
                cz(i+2) = GrZ(zint+i)    

            end do

c        Create the interpolation template for each of the components
c            of the tensor

            do p=1,NR
                do q=1,NC
                    do i=-1,1
                        do j=-1,1
                            do k=-1,1

                  ftmpl(i+2,j+2,k+2) = GrT(xint+i, yint+j, zint+k,
     *                                                p, q)

                            end do
                        end do
                    end do
                    call DINT3D3(ftmpl,cx,cy,cz,ix,it(p,q),hx,hy,hz)
                end do
            end do

            return
        end
c=======================================================================
c
c    Tensor Interpolation routine.
c    Interpolates rank 3 tensors to 3rd order
c    Arguments
c    GrT         := Grid Tensor
c    Gr[X-Z]    := Coordinates on grid
c    ix            := Interpolation point
c    it            := Return value for tensor
c    N[C,R]    := Number of components
c    N[X-Z]    := Number of grid points
c
        subroutine Int3Tsr3(GrT,GrX,GrY,GrZ,ix,it,hx,hy,hz,
     #                    N1,N2,N3,NX,NY,NZ)
            implicit none

            integer    N1, N2, N3,
     *                NX, NY, NZ

            real*8    GrT(NX, NY, NZ, N1, N2, N3),
     *                GrX(NX), GrY(NY), GrZ(NZ),
     *                ix(3),    it(N1,N2,N3),
     *                hx, hy, hz

c        Local Variables

            integer i, j, k, xint, yint, zint,
     *                p, q, r

            real*8   ftmpl(3,3,3),
     *            cx(3), cy(3), cz(3)

c     Functions


c     Search for neighboring points for the interpolation.
c     Note that interpolation at the edges is not carried out here.

c        Calculate integer coordinates for interpolation point

         xint = int((ix(1)  - GrX(1))/hx) + 1
         yint = int((ix(2)  - GrY(1))/hy) + 1
         zint = int((ix(3)  - GrZ(1))/hz) + 1

c            Coordinates...

            do i=-1,1

                cx(i+2) = GrX(xint+i)    
                cy(i+2) = GrY(yint+i)    
                cz(i+2) = GrZ(zint+i)    

            end do

c        Create the interpolation template for each of the components
c            of the tensor

            do p=1,N1
                do q=1,N2
                    do r=1,N3
                        do i=-1,1
                            do j=-1,1
                                do k=-1,1

                        ftmpl(i+2,j+2,k+2) =GrT(xint+i, yint+j, zint+k,
     *                                            p, q, r)

                                end do
                            end do
                        end do
                      call DINT3D3(ftmpl,cx,cy,cz,ix,it(p,q,r),hx,hy,hz)
                    end do
                end do
            end do

            return
        end
c=======================================================================
    
c    3D Interpolation routine.  Created Tue Apr  5 14:00:39 CDT 1994 MFH
c        4th Order


        subroutine DINT4D3(ftmpl,cx,cy,cz,xi,fintp,hx,hy,hz)

            implicit none

c            Function value cube
            real*8    ftmpl(4,4,4)

c            Coordinate value cube
            real*8    cx(4), cy(4), cz(4)

c            Point to interpolate to...
            real*8    xi(3)

c            Return value for interpolated value for function
            real*8    fintp

c            Spacings...
            real*8    hx, hy, hz

c        Local Variables

            integer    i, j, k,
     *            ix, iy, iz

            real*8   fx(4), fy(4), fz(4)

c        Function prototypes
            
            real*8    L4intp


c            Loop over the z-direction
            do k=1,4

c                Loop over the y-direction
                do j=1,4

c                    Load in the function values in the x-direction
                    do i=1,4
                        fx(i) = ftmpl(i,j,k)
                    enddo

c                    Interpolate in the x-direction
                    fy(j) = L4intp(fx,cx,xi(1))
                
                end do

c                Now interpolate in the y-direction
                fz(k) = L4intp(fy,cy,xi(2))    

            end do

c            Now interpolate in the z-direction
            fintp = L4intp(fz,cz,xi(3)) 
                    
            return
        end
c=======================================================================
c    Interpolate Scalar function Tue Apr  5 14:17:04 CDT 1994 MFH
c    4th Order...
c
c    Routine returns the value at point ix
c
c    GrF       := Grid function
c    Gr[X-Z] := Coordinates on grid
c    ix          := Point at which interpolation is to be done
c    h[x-z]  := grid spacings
c    N[X-Z]  := Grid sizes
c

        real*8 function Intp4Sc3(GrF,GrX,GrY,GrZ,ix,hx,hy,hz,NX,NY,NZ)
    
            implicit none

            integer    NX, NY, NZ

            real*8    GrF( NX, NY, NZ),
     *               GrX(NX),  GrY(NY),  GrZ(NZ),
     *            ix(3),
     *                hx, hy, hz
        

c        Local Variables

            integer    i, j, k, xint, yint, zint

            real*8    ftmpl(4,4,4),
     *                cx(4), cy(4), cz(4)

c        Functions


c        Search for neighboring points for the interpolation.
c        Note that interpolation at the edges is not carried out here.

c            Calculate integer coordinates for interpolation point

            xint = int((ix(1)  - GrX(1))/hx) + 1
            yint = int((ix(2)  - GrY(1))/hy) + 1
            zint = int((ix(3)  - GrZ(1))/hz) + 1

c            Create the interpolation template

            do i=-1,2
                do j=-1,2
                    do k=-1,2

                        ftmpl(i+2,j+2,k+2) = GrF(xint+i, yint+j, zint+k)

                    end do
                end do
            end do
            
c            Coordinates...

            do i=-1,2

                cx(i+2) = GrX(xint+i)    
                cy(i+2) = GrY(yint+i)    
                cz(i+2) = GrZ(zint+i)    

            end do

c            Interpolate
            call DINT4D3(ftmpl,cx,cy,cz,ix,Intp4Sc3,hx,hy,hz)
c            call DINT4ADPT(ftmpl,cx,cy,cz,ix,Intp4Sc3,hx,hy,hz)

            return
        end
c=======================================================================
c
c    Tensor Interpolation routine.
c    Interpolates rank 2 tensors to 4th order
c    Arguments
c    GrT         := Grid Tensor
c    Gr[X-Z]    := Coordinates on grid
c    ix            := Interpolation point
c    it            := Return value for tensor
c    N[C,R]    := Number of components
c    N[X-Z]    := Number of grid points
c
        subroutine Int2Tsr4(GrT,GrX,GrY,GrZ,ix,it,hx,hy,hz,
     #                    NR,NC,NX,NY,NZ)
            implicit none

            integer    NR, NC,
     *                NX, NY, NZ

            real*8    GrT(NX, NY, NZ, NR, NC),
     *                GrX(NX), GrY(NY), GrZ(NZ),
     *                ix(3),    it(NR,NC),
     *                hx, hy, hz

c        Local Variables

            integer i, j, k, xint, yint, zint,
     *                p, q

            real*8   ftmpl(4,4,4),
     *            cx(4), cy(4), cz(4)

c     Functions


c     Search for neighboring points for the interpolation.
c     Note that interpolation at the edges is not carried out here.

c        Calculate integer coordinates for interpolation point

         xint = int((ix(1)  - GrX(1))/hx) + 1
         yint = int((ix(2)  - GrY(1))/hy) + 1
         zint = int((ix(3)  - GrZ(1))/hz) + 1

c            Coordinates...

            do i=-1,2

                cx(i+2) = GrX(xint+i)    
                cy(i+2) = GrY(yint+i)    
                cz(i+2) = GrZ(zint+i)    

            end do

c        Create the interpolation template for each of the components
c            of the tensor

            do p=1,NR
                do q=1,NC
                    do i=-1,2
                        do j=-1,2
                            do k=-1,2

                  ftmpl(i+2,j+2,k+2) = GrT(xint+i, yint+j, zint+k,
     *                                                p, q)

                            end do
                        end do
                    end do
                    call DINT4D3(ftmpl,cx,cy,cz,ix,it(p,q),hx,hy,hz)
c                    call DINT4ADPT(ftmpl,cx,cy,cz,ix,it(p,q),hx,hy,hz)
                end do
            end do

            return
        end
c=======================================================================
c
c    Tensor Interpolation routine.
c    Interpolates rank 3 tensors to 4th order
c    Arguments
c    GrT         := Grid Tensor
c    Gr[X-Z]    := Coordinates on grid
c    ix            := Interpolation point
c    it            := Return value for tensor
c    N[C,R]    := Number of components
c    N[X-Z]    := Number of grid points
c
        subroutine Int3Tsr4(GrT,GrX,GrY,GrZ,ix,it,hx,hy,hz,
     #                    N1,N2,N3,NX,NY,NZ)
            implicit none

            integer    N1, N2, N3,
     *                NX, NY, NZ

            real*8    GrT(NX, NY, NZ, N1, N2, N3),
     *                GrX(NX), GrY(NY), GrZ(NZ),
     *                ix(3),    it(N1,N2,N3),
     *                hx, hy, hz

c        Local Variables

            integer i, j, k, xint, yint, zint,
     *                p, q, r

            real*8   ftmpl(4,4,4),
     *            cx(4), cy(4), cz(4)

c     Functions


c     Search for neighboring points for the interpolation.
c     Note that interpolation at the edges is not carried out here.

c        Calculate integer coordinates for interpolation point

         xint = int((ix(1)  - GrX(1))/hx) + 1
         yint = int((ix(2)  - GrY(1))/hy) + 1
         zint = int((ix(3)  - GrZ(1))/hz) + 1

c            Coordinates...

            do i=-1,2

                cx(i+2) = GrX(xint+i)    
                cy(i+2) = GrY(yint+i)    
                cz(i+2) = GrZ(zint+i)    

            end do

c        Create the interpolation template for each of the components
c            of the tensor

            do p=1,N1
                do q=1,N2
                    do r=1,N3
                        do i=-1,2
                            do j=-1,2
                                do k=-1,2

                         ftmpl(i+2,j+2,k+2) =GrT(xint+i, yint+j, zint+k,
     *                                            p, q, r)

                                end do
                            end do
                        end do
                      call DINT4D3(ftmpl,cx,cy,cz,ix,it(p,q,r),hx,hy,hz)
c                        call DINT4ADPT(ftmpl,cx,cy,cz,ix,it(p,q,r),hx,hy,hz)
                    end do
                end do
            end do

            return
        end
c=======================================================================
c    The following routine carries out cubic interpolation with half points
        subroutine DINT3ADPT(ftmpl,cx,cy,cz,xi,fintp,hx,hy,hz)

            implicit none

c            Function value cube
            real*8    ftmpl(3,3,3)

c            Coordinate value cube
            real*8    cx(3), cy(3), cz(3)

c            Point to interpolate to...
            real*8    xi(3)

c            Return value for interpolated value for function
            real*8    fintp

c            Spacings...
            real*8    hx, hy, hz

c        Local Variables

            integer    i, j, k,
     *            ix, iy, iz

            real*8   fx(3), fy(3), fz(3),
     *                fhalf(3,3,3), chx(3),chy(3),chz(3),
     *                tx(3)

c        Function prototypes
            
            real*8    L3intp

c            Set up center point
            chx(2) = cx(2)
            chy(2) = cy(2)
            chz(2) = cz(2)

c            Set up rest of the half points
            chx(1) = chx(2) - hx*0.5d0
            chy(1) = chy(2) - hy*0.5d0
            chz(1) = chz(2) - hz*0.5d0

            chx(3) = chx(2) + hx*0.5d0
            chy(3) = chy(2) + hy*0.5d0
            chz(3) = chz(2) + hz*0.5d0

            do i=1,3
                do j=1,3
                    do k=1,3

                        tx(1) = chx(i)
                        tx(2) = chy(j)
                        tx(3) = chz(k)

                   call DINT3D3(ftmpl,cx,cy,cz,tx,fhalf(i,j,k),hx,hy,hz)

                    end do
                end do
            end do
            

                                    
c            Loop over the z-direction
            do k=1,3

c                Loop over the y-direction
                do j=1,3

c                    Load in the function values in the x-direction
                    do i=1,3
                        fx(i) = fhalf(i,j,k)
                    enddo

c                    Interpolate in the x-direction
                    fy(j) = L3intp(fx,chx,xi(1))
                
                end do

c                Now interpolate in the y-direction
                fz(k) = L3intp(fy,chy,xi(2))    

            end do

c            Now interpolate in the z-direction
            fintp = L3intp(fz,chz,xi(3)) 
                    
            return
        end
c=======================================================================
c    The following routine carries out fourth order interpolation 
        subroutine DINT4ADPT(ftmpl,cx,cy,cz,xi,fintp,hx,hy,hz)

            implicit none

c            Function value cube
            real*8    ftmpl(4,4,4)

c            Coordinate value cube
            real*8    cx(4), cy(4), cz(4)

c            Point to interpolate to...
            real*8    xi(3)

c            Return value for interpolated value for function
            real*8    fintp

c            Spacings...
            real*8    hx, hy, hz

c        Local Variables

            integer    i, j, k,
     *            ix, iy, iz

            real*8   fx(4), fy(4), fz(4),
     *                fhalf(4,4,4), chx(4),chy(4),chz(4),
     *                tx(3)

c        Function prototypes
            
            real*8    L4intp

c            Set up center point
            chx(2) = cx(2)
            chy(2) = cy(2)
            chz(2) = cz(2)

c            Set up rest of the half points
            chx(1) = chx(2) - hx*0.5d0
            chy(1) = chy(2) - hy*0.5d0
            chz(1) = chz(2) - hz*0.5d0

            chx(3) = chx(2) + hx*0.5d0
            chy(3) = chy(2) + hy*0.5d0
            chz(3) = chz(2) + hz*0.5d0

            chx(4) = chx(2) + hx
            chy(4) = chy(2) + hy
            chz(4) = chz(2) + hz

            do i=1,4
                do j=1,4
                    do k=1,4

                        tx(1) = chx(i)
                        tx(2) = chy(j)
                        tx(3) = chz(k)

                   call DINT4D3(ftmpl,cx,cy,cz,tx,fhalf(i,j,k),hx,hy,hz)

                    end do
                end do
            end do
            

                                    
c            Loop over the z-direction
            do k=1,4

c                Loop over the y-direction
                do j=1,4

c                    Load in the function values in the x-direction
                    do i=1,4
                        fx(i) = fhalf(i,j,k)
                    enddo

c                    Interpolate in the x-direction
                    fy(j) = L4intp(fx,chx,xi(1))
                
                end do

c                Now interpolate in the y-direction
                fz(k) = L4intp(fy,chy,xi(2))    

            end do

c            Now interpolate in the z-direction
            fintp = L4intp(fz,chz,xi(3)) 
                    
            return
        end
c=======================================================================
c        Grid to Grid interpolation 
c    
c    Given a spherical grid of given size interpolate to a new spherical
c    Grid of either smaller or larger size.
c
c    Input:
c        Sc1,S1,Sf1 (Cartesian,Spherical coordinates and function)1
c        Sizes in NST1 and NSP1 of Grid 1
c        Sizes in NST2 and NSP2 of Grid 2
c
c    Output:
c        Sc2,S2,Sf2 (Cartesian,Spherical coordinates and function)2
c
        subroutine S2Grid2GridIntp(S1,Sf1,S2,Sf2,tmpl,
     *                                    NST1,NSP1,NST2,NSP2)
            implicit none
            
            integer    NST1, NSP1, NST2, NSP2,
     *                tmpl(NST1,NSP1,4,4,2)

            real*8    S1(NST1,NSP1,3), Sf1(NST1,NSP1),
     *                S2(NST2,NSP2,3), Sf2(NST2,NSP2)

c        Local Variables

            integer    i,j

            real*8    dtheta1,dphi1, dtheta2, dphi2, Pi,
     *                r1,t,p, temp, ff

            Pi = 4.0d0*atan(1.0d0)

            dtheta1 = Pi / float(NST1 - 1)
            dphi1 = 2.0d0*Pi / float(NSP1)
            dtheta2 = Pi / float(NST2 - 1)
            dphi2 = 2.0d0*Pi / float(NSP2)

        
c        Loop over the new grid points

c        First interpolate r
            do i=1,NST2
                t = dtheta2*float(i-1)
                do j=1,NSP2
                    p = dphi2*float(j-1)
    
                    S2(i,j,2) = t
                    S2(i,j,3) = p

                        
            call Intp2Dtmpl(S1,Sf1,p,t,ff,temp,tmpl,4,
     *         dtheta1,dphi1,NST1,NSP1,3)

                    Sf2(i,j) = ff


                end do
            end do

            return
        end





c=======================================================================


        subroutine Intp2Dtmpl(S,Sf,phi,theta,iff,df,tmpl,order,
     *            dtheta,dphi,NST,NSP,NDIM)
            implicit none

            integer NST
            integer NSP
            integer NDIM
            integer order

            real*8    S(NST,NSP,NDIM)
            real*8    Sf(NST,NSP)
            real*8    phi
            real*8    theta
            real*8    iff
            real*8    df
            integer    tmpl(NST,NSP,4,4,2)
            real*8    dphi
            real*8    dtheta

    
c    Local Variables

            integer    i,j,k
            integer     nt,np

            integer     across

            integer    flag
            
            real*8    rtt(4),rpp(4)
            real*8    rff(4,4)
            real*8    tff(4)
            real*8    ff(4)

            real*8    prem,trem
            real*8    newp,newt


            logical    PRINT_TEMPLATE
    
c    Function Calls

            real*8    L4Interp

c    Start by finding nearest integer coordinates of phi and theta

            np = int(phi/dphi) + 1
            nt = int(theta/dtheta) + 1




c    np = 1     -> phi = 0      ;   nt=1   -> theta = 0
c    np = NSP -> phi = 2Pi ;   nt=NST -> theta = Pi

c    Find the remainders.
        
            prem = phi - S(nt,np,3)
            trem = theta - S(nt,np,2)

c            if(nt .lt. 0 .or. nt .gt. NST)then
c                write(6,*)'Error in MyIntp: Improper value of nt'
c                write(6,*)'Value = ',nt
c                stop
c            endif


c  Set up the function template and phi and theta

         do i=1,4
            rpp(i) = float(i)
            rtt(i) = float(i)
            do j=1,4
               rff(i,j) = Sf(tmpl(nt,np,i,j,1),tmpl(nt,np,i,j,2))
            enddo
         enddo

            
            newp = 2.0d0 + prem/dphi  
            newt = 2.0d0 + trem/dtheta

       do i=1,4
                do j=1,4
                    ff(j) = rff(j,i)
                end do
                tff(i) = L4Interp(ff,rtt,1.0d0,newt)
       enddo
            iff = L4Interp(tff,rpp,1.0d0,newp)



            


            
            return

            end
c=======================================================================
c     Creates the interpolation templates
c
      subroutine CreateTemplate(tmpl,NST,NSP)
         implicit none

         integer  NST,NSP
         integer   tmpl(NST,NSP,4,4,2)

c     Local Variables

         integer  nt,np,i,j

         do nt=1,NST
            do np=1,NSP
               do i=1,4
                  do j=1,4

                     tmpl(nt,np,i,j,1)=nt + (i-2)
                     tmpl(nt,np,i,j,2)=np + (j-2)


c  Theta's
                     if(tmpl(nt,np,i,j,1) .le. 0)then
                        tmpl(nt,np,i,j,1) = abs(tmpl(nt,np,i,j,1))+2
                    tmpl(nt,np,i,j,2) = tmpl(nt,np,i,j,2)+int(NSP/2)
                        if(tmpl(nt,np,i,j,2) .gt. NSP)then
                    tmpl(nt,np,i,j,2) = tmpl(nt,np,i,j,2) - NSP
                        endif
                     endif
                     if(tmpl(nt,np,i,j,1) .gt. NST)then
                        tmpl(nt,np,i,j,1) = 2*NST-tmpl(nt,np,i,j,1)
                  tmpl(nt,np,i,j,2) = tmpl(nt,np,i,j,2)+int(NSP/2)
                        if(tmpl(nt,np,i,j,2) .gt. NSP)then
                     tmpl(nt,np,i,j,2) = tmpl(nt,np,i,j,2) - NSP
                       endif
                     endif

c  Phi's
                     if(tmpl(nt,np,i,j,2) .gt. NSP)then
                        tmpl(nt,np,i,j,2)=tmpl(nt,np,i,j,2)-NSP
                     endif
                     if(tmpl(nt,np,i,j,2) .le. 0)then
                        tmpl(nt,np,i,j,2)=tmpl(nt,np,i,j,2)+NSP
                     endif

                  enddo
               enddo
            enddo
         enddo

         return
      end
c=======================================================================
c    Routine to calculate interpolation coeffients for 4th order Lagrange 
c    interpolation. 
c    Input:
c        pt:= point at which coefficient is desired.
c        list:= List of points that define the interpolation neighborhood
c                in 1D.
c        Ensure that the point pt is the one closest to the 2nd point
c        in the list. Preferably
c        x  x o x  x
c    
c    Convergence tested for 1D interpolation problem.
c    See test622a.f in tests/ directory.
c    Wed Jun 22 21:14:19 CDT 1994 MFH
c
        subroutine lgrncf4(pt,list,cfs)
            implicit none

            real*8    pt, list(4), cfs(4)

c        Local Variables
            
            integer    i

            real*8    dd

c    Evaluate the coefficients

            dd = (list(1) - list(2))*(list(1) - list(3))*
     *          (list(1) - list(4))

            cfs(1) = ((pt - list(2))*(pt - list(3))*(pt - list(4)))/dd

            dd = (list(2) - list(1))*(list(2) - list(3))*
     *          (list(2) - list(4))

            cfs(2) = ((pt - list(1))*(pt - list(3))*(pt - list(4)))/dd

            dd = (list(3) - list(1))*(list(3) - list(2))*
     *          (list(3) - list(4))

            cfs(3) = ((pt - list(1))*(pt - list(2))*(pt - list(4)))/dd

            dd = (list(4) - list(1))*(list(4) - list(2))*
     *          (list(4) - list(3))

            cfs(4) = ((pt - list(1))*(pt - list(2))*(pt - list(3)))/dd

            return
        end
c-----------------------------------------------------------------------
c    Routine to calculate the derivative of the interpolation
c    coefficients.
c    Arguments:
c        pt:= Point at which it is to be evaluated.
c        list:= List of points embedding the interpolation point.
c        dcfs:= Returns an array of derivatives of the interpolation coeffs.
c
c        Ensure that the point pt is the one closest to the 2nd point
c        in the list. Preferably
c        x  x o x  x
c
        subroutine lgrndcf4(pt,list,dcfs)
            implicit none
    
            real*8    pt, list(4), dcfs(4)

c        Local Variables
    
            real*8 dd,sum

            dd = (list(1) - list(2))*(list(1) - list(3))*
     *          (list(1) - list(4))

            sum = (pt - list(3))*(pt - list(4)) +
     *            (pt - list(2))*(pt - list(4)) +
     *            (pt - list(2))*(pt - list(3))

            dcfs(1) = sum / dd
        
            dd = (list(2) - list(1))*(list(2) - list(3))*
     *          (list(2) - list(4))

            sum = (pt - list(3))*(pt - list(4)) +
     *            (pt - list(1))*(pt - list(4)) +
     *            (pt - list(1))*(pt - list(3))

            dcfs(2) = sum / dd
            
            dd = (list(3) - list(1))*(list(3) - list(2))*
     *          (list(3) - list(4))

            sum = (pt - list(2))*(pt - list(4)) +
     *            (pt - list(1))*(pt - list(4)) +
     *            (pt - list(1))*(pt - list(2))

            dcfs(3) = sum / dd

            dd = (list(4) - list(1))*(list(4) - list(2))*
     *          (list(4) - list(3))

            sum = (pt - list(2))*(pt - list(3)) +
     *            (pt - list(1))*(pt - list(3)) +
     *            (pt - list(1))*(pt - list(2))

            dcfs(4) = sum / dd


            return

        end
c-----------------------------------------------------------------------
c  Loads in interpolation template from template.dat
c
      subroutine LoadTemplate(tmpl,NST,NSP)
         implicit none

         integer  NST,NSP
         integer  nt,np,i,j
         integer  counter

         integer  tmpl(NST,NSP,4,4,2)


         open(unit=97,file='template.dat')

         counter = 0
100      continue
         read(97,1212,end=200)nt,np,i,j,tmpl(nt,np,i,j,1),
     #               tmpl(nt,np,i,j,2)

         counter = counter + 1
         goto 100

200      continue

         write(0,*)'Read',counter,'lines'

1212     format(6I4)
         return
      end


c-----------------------------------------------------------------------
