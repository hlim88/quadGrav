c=======================================================================
c    Compute the derivatives of the 3-metric
c  Tested for diagonal terms 2-21-94
c    
        subroutine GrGTDeriv(GrG,GrDG,hx,hy,hz,temp,temp2,NX,NY,NZ)
            implicit none

            integer    NX,NY,NZ

            real*8    GrG(NX,NY,NZ,3,3)
            real*8    GrDG(NX,NY,NZ,3,3,3)
            real*8    temp(NX),temp2(NX)

            real*8    hx,hy,hz

c    Local Variables

            integer    i,j,k,m,p,q

            real*8    divhx,divhy,divhz

c    Declare temporary arrays



c    Zero out derivative array
            do i=1,NX
                do j=1,NY
                    do k=1,NZ

                        do m=1,3
                            do p=1,3
                                do q=1,3
                                    GrDG(i,j,k,m,p,q) = 0.0e0
                                enddo
                            enddo
                        enddo

                    enddo
                enddo
            enddo

c    Compute the X derivatives
c    1,1 - 3 components
            do m=1,3
                do i=1,NY
                    do j=1,NZ
                        do k=1,NX

                            temp(k)=GrG(k,i,j,1,m)

                        enddo

                        call DVD2Q1(temp,temp2,hx,NX)

                        do k=1,NX

                            GrDG(k,i,j,1,m,1) = temp2(k)

                        enddo
                    enddo
                enddo
            enddo

c    2,2-3 components
            do m=2,3
                do i=1,NY
                    do j=1,NZ
                        do k=1,NX

                            temp(k)=GrG(k,i,j,2,m)

                        enddo

                        call DVD2Q1(temp,temp2,hx,NX)

                        do k=1,NX

                            GrDG(k,i,j,2,m,1) = temp2(k)

                        enddo
                    enddo
                enddo
            enddo
                            
c    3,3 components
            do i=1,NY
                do j=1,NZ
                    do k=1,NX
            
                        temp(k) = GrG(k,i,j,3,3)
        
                    enddo

                    call DVD2Q1(temp,temp2,hx,NX)

                    do k=1,NX

                        GrDG(k,i,j,3,3,1) = temp2(k)

                    enddo
                enddo
            enddo
        
c    Y derivatives                                        
c    1,1 - 3 components
            do m=1,3
                do i=1,NX
                    do j=1,NZ
                        do k=1,NY

                            temp(k)=GrG(i,k,j,1,m)

                        enddo

                        call DVD2Q1(temp,temp2,hy,NY)

                        do k=1,NY

                            GrDG(i,k,j,1,m,2) = temp2(k)

                        enddo
                    enddo
                enddo
            enddo

c    2,2-3 components
            do m=2,3
                do i=1,NX
                    do j=1,NZ
                        do k=1,NY

                            temp(k)=GrG(i,k,j,2,m)

                        enddo

                        call DVD2Q1(temp,temp2,hy,NY)

                        do k=1,NY

                            GrDG(i,k,j,2,m,2) = temp2(k)

                        enddo
                    enddo
                enddo
            enddo
                            
c    3,3 components
            do i=1,NX
                do j=1,NZ
                    do k=1,NY
            
                        temp(k) = GrG(i,k,j,3,3)
        
                    enddo

                    call DVD2Q1(temp,temp2,hy,NY)

                    do k=1,NY

                        GrDG(i,k,j,3,3,2) = temp2(k)

                    enddo
                enddo
            enddo
        
c    Z Derivatives
c    1,1 - 3 components
            do m=1,3
                do i=1,NX
                    do j=1,NY
                        do k=1,NZ

                            temp(k)=GrG(i,j,k,1,m)

                        enddo

                        call DVD2Q1(temp,temp2,hz,NZ)

                        do k=1,NZ

                            GrDG(i,j,k,1,m,3) = temp2(k)

                        enddo
                    enddo
                enddo
            enddo

c    2,2-3 components
            do m=2,3
                do i=1,NX
                    do j=1,NY
                        do k=1,NZ

                            temp(k)=GrG(i,j,k,2,m)

                        enddo

                        call DVD2Q1(temp,temp2,hz,NZ)

                        do k=1,NZ

                            GrDG(i,j,k,2,m,3) = temp2(k)

                        enddo
                    enddo
                enddo
            enddo
                            
c    3,3 components
            do i=1,NX
                do j=1,NY
                    do k=1,NZ
            
                        temp(k) = GrG(i,j,k,3,3)
        
                    enddo

                    call DVD2Q1(temp,temp2,hz,NZ)

                    do k=1,NZ

                        GrDG(i,j,k,3,3,3) = temp2(k)

                    enddo
                enddo
            enddo
        
c    Copy  to the lower diagonals

            do i=1,NX
                do j=1,NY
                    do k=1,NZ

                        GrDG(i,j,k,2,1,1) = GrDG(i,j,k,1,2,1)
                        GrDG(i,j,k,3,1,1) = GrDG(i,j,k,1,3,1)
                        GrDG(i,j,k,3,2,1) = GrDG(i,j,k,2,3,1)
                        GrDG(i,j,k,2,1,2) = GrDG(i,j,k,1,2,2)
                        GrDG(i,j,k,3,1,2) = GrDG(i,j,k,1,3,2)
                        GrDG(i,j,k,3,2,2) = GrDG(i,j,k,2,3,2)
                        GrDG(i,j,k,2,1,3) = GrDG(i,j,k,1,2,3)
                        GrDG(i,j,k,3,1,3) = GrDG(i,j,k,1,3,3)
                        GrDG(i,j,k,3,2,3) = GrDG(i,j,k,2,3,3)
                    enddo
                enddo
            enddo

c    Done

            return

        end
c-----------------------------------------------------------------------

c    Tue Aug 17 18:27:42 CDT 1993 MFH
c
        subroutine GrGDeriv(GrG,DGrGD,hx,hy,hz,NX,NY,NZ)
            implicit none

            integer    NX,NY,NZ

            real*8    GrG(NX,NY,NZ,3,3)
            real*8    DGrGD(NX,NY,NZ,3,3,3)

            real*8    hx,hy,hz
    
c    Local Variables

            integer    i,j,k,p,q

            real*8    divhx
            real*8    divhy
            real*8    divhz


            divhx = 0.5e0/hx
            divhy = 0.5e0/hy
            divhz = 0.5e0/hz
c


c    Zero out the Derivatives Array
            do i=1,NX
                do j=1,NY
                    do k=1,NZ
                        do p=1,3
                        do q=1,3
                        DGrGD(i,j,k,p,q,1) = 0.0e0
                        DGrGD(i,j,k,p,q,2) = 0.0e0
                        DGrGD(i,j,k,p,q,3) = 0.0e0
                        end do
                        end do
                    enddo
                enddo
            enddo

c    Do the interior
            do i=2,NX-1
                do j=2,NY-1
                    do k=2,NZ-1
                        do p=1,3
                        do q=1,3

                        DGrGD(i,j,k,p,q,1) =divhx*(GrG(i+1,j,k,p,q) 
     *                                        - GrG(i-1,j,k,p,q))
                        DGrGD(i,j,k,p,q,2) =divhy*(GrG(i,j+1,k,p,q) 
     *                                        - GrG(i,j-1,k,p,q))
                        DGrGD(i,j,k,p,q,3) =divhz*(GrG(i,j,k+1,p,q) 
     *                                        - GrG(i,j,k-1,p,q))
                        
                        end do
                        end do
                    enddo
                enddo
            enddo



            return
        end
c=======================================================================
c
c    Routine to generate metric for tests
c
        subroutine GenMet(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
     $      k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,NX,NY,NZ)
            implicit none

            integer    NX, NY, NZ

            real*8    hx, hy, hz,
     $    g11(NX,NY,NZ),g12(NX,NY,NZ),g13(NX,NY,NZ),
     $    g22(NX,NY,NZ),g23(NX,NY,NZ),g33(NX,NY,NZ),
     $    k11(NX,NY,NZ),k12(NX,NY,NZ),k13(NX,NY,NZ),
     $    k22(NX,NY,NZ),k23(NX,NY,NZ),k33(NX,NY,NZ),
     $    msk(NX,NY,NZ),GrX(NX), GrY(NY), GrZ(NZ)

c        Local Variables

            integer    i, j, k
            real*8 psi, rr, max, min

            max = 2.0e0
            min = -2.0e0

            hx = (max - min) / float(NX-1) 
            hy = hx
            hz = hx

c        set up coordinates
    
            do i = 1, NX
                GrX(i) = hx * float(i-1) + min
                GrY(i) = hx * float(i-1) + min
                GrZ(i) = hx * float(i-1) + min
c                write(6,*)GrX(i),GrX(i), GrY(i)
            end do

c        Now the metric

            do k = 1, NZ
                do j = 1, NY
                    do i = 1,NX
                        rr = sqrt(grx(i)**2 + gry(j)**2 + grz(k)**2)
                        if(rr .ge. 0.50e0)then
                            msk(i,j,k) = 1.0
                            
                            g11(i,j,k) = (1.0e0 + 0.5e0/rr)**4
                            g12(i,j,k) = 0.0e0
                            g13(i,j,k) = 0.0e0
                            g22(i,j,k) = g11(i,j,k)
                            g23(i,j,k) = 0.0e0
                            g33(i,j,k) = g11(i,j,k)

                            k11(i,j,k) = 0.0e0
                            k12(i,j,k) = 0.0e0
                            k13(i,j,k) = 0.0e0
                            k22(i,j,k) = 0.0e0
                            k23(i,j,k) = 0.0e0
                            k33(i,j,k) = 0.0e0

                        else
                            msk(i,j,k) = 0.0e0
                        end if
                    end do
                end do
            end do

            return
        end
c-----------------------------------------------------------------------
c    The following routine will generate the 3-metric and extrinsic 
c    curvature data along with coordinates and characteristic function.
c
        subroutine GenEF(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
     $      k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,NX,NY,NZ)
            implicit none

            integer    NX, NY, NZ

            real*8    hx, hy, hz,
     $    g11(NX,NY,NZ),g12(NX,NY,NZ),g13(NX,NY,NZ),
     $    g22(NX,NY,NZ),g23(NX,NY,NZ),g33(NX,NY,NZ),
     $    k11(NX,NY,NZ),k12(NX,NY,NZ),k13(NX,NY,NZ),
     $    k22(NX,NY,NZ),k23(NX,NY,NZ),k33(NX,NY,NZ),
     $    msk(NX,NY,NZ),GrX(NX), GrY(NY), GrZ(NZ)
c    Local Variables

            integer    i,j,k,m,p
            real*8 rr, max, min, ccx,ccy,ccz,mx,my,mz,mskrad, x,y,z
            
            max = 2.0e0
            min = -2.0e0

c    Read in grid properties

            write(0,*)'Reading in eddfink.dat'
            open(unit=1,file='eddfink.dat')
            read(1,*)min,max
            read(1,*)ccx,ccy,ccz

c            mask origin and radius
            read(1,*)mskrad
            read(1,*)mx,my,mz
            close(unit=1)
            
            write(0,*)'Read in min, max ',min,max
            write(0,*)'Read in coordinates of hole',ccx,ccy,ccz
            write(2,*)'Read in min, max ',min,max
            write(2,*)'Read in coordinates of hole',ccx,ccy,ccz
            write(2,*)'Mask radius ',mskrad
            write(2,*)'Mask location ',mx,my,mz

            hx = (max - min) / float(NX-1) 
            hy = hx
            hz = hx

            write(2,*)'Cartesian grid size hx=hy=hz=',hx

c        set up coordinates
    
            do i = 1, NX
                GrX(i) = hx * float(i-1) + min
                GrY(i) = hx * float(i-1) + min
                GrZ(i) = hx * float(i-1) + min
            end do

c        set up the mask
            do i = 1, NX
                do j = 1, NY
                    do k = 1, NZ
                        rr = sqrt((GrX(i)-mx)**2 + (GrY(j)-my)**2 +
     .                     (GrZ(k)-mz)**2)
                        if(rr .gt. mskrad)then
                            msk(i,j,k) = 1.0
                        else
                            msk(i,j,k) = 0.0
                        end if
                    end do
                end do
            end do

            do i = 1, NX
                do j = 1, NY
                    do k = 1, NZ
                    
                    if(msk(i,j,k) .eq. 0)then
                        g11(i,j,k) = 0.0
                        g12(i,j,k) = 0.0
                        g13(i,j,k) = 0.0
                        g22(i,j,k) = 0.0
                        g23(i,j,k) = 0.0
                        g33(i,j,k) = 0.0
                        k11(i,j,k) = 0.0
                        k12(i,j,k) = 0.0
                        k13(i,j,k) = 0.0
                        k22(i,j,k) = 0.0
                        k23(i,j,k) = 0.0
                        k33(i,j,k) = 0.0
                    else

                        x = grx(i) - ccx
                        y = gry(j) - ccy
                        z = grz(k) - ccz

c                metric functions first

                    g11(i,j,k) = (y**2*sqrt(x**2+y**2+z**2)+2*x**2+
     #            sqrt(x**2+y**2+z**2)*z**2+sqrt(x**2+y**2+z**2)
     #            *x**2)/sqrt(x**2+y**2+z**2)**3

                    g12(i,j,k) = 2*x*y/sqrt(x**2+y**2+z**2)**3

                    g13(i,j,k) = 2*z*x/sqrt(x**2+y**2+z**2)**3

                    g22(i,j,k) = (2*y**2+y**2*sqrt(x**2+y**2+z**2)+
     #        sqrt(x**2+y**2+z**2)*z**2+sqrt(x**2+y**2+z**2)*x**2)
     #        /sqrt(x**2+y**2+z**2)**3

                    g23(i,j,k) = 2*z*y/sqrt(x**2+y**2+z**2)**3

                    g33(i,j,k) = (y**2*sqrt(x**2+y**2+z**2)+2*z**2+
     #        sqrt(x**2+y**2+z**2)*z**2+sqrt(x**2+y**2+z**2)*x**2)/
     #        sqrt(x**2+y**2+z**2)**3



                    k11(i,j,k) = 2*(y**2*sqrt(x**2+y**2+z**2)-x**2+
     #        sqrt(x**2+y**2+z**2)*z**2-sqrt(x**2+y**2+z**2)*x**2)/
     #        sqrt((sqrt(x**2+y**2+z**2)+2)/sqrt(x**2+y**2+z**2))/
     #        sqrt(x**2+y**2+z**2)**5
        
                    k12(i,j,k) = -2*y*x*(2*sqrt(x**2+y**2+z**2)+1)/
     #        sqrt((sqrt(x**2+y**2+z**2)+2)/sqrt(x**2+y**2+z**2))/
     #        sqrt(x**2+y**2+z**2)**5

                    k13(i,j,k) = -2*(2*sqrt(x**2+y**2+z**2)+1)*z*x/
     #        sqrt((sqrt(x**2+y**2+z**2)+2)/sqrt(x**2+y**2+z**2))/
     #        sqrt(x**2+y**2+z**2)**5

                    k22(i,j,k) = -2*(y**2*sqrt(x**2+y**2+z**2)+y**2-
     #        sqrt(x**2+y**2+z**2)*x**2-sqrt(x**2+y**2+z**2)*z**2)/
     #        sqrt((sqrt(x**2+y**2+z**2)+2)/sqrt(x**2+y**2+z**2))/
     #        sqrt(x**2+y**2+z**2)**5

                    k23(i,j,k) = -2*(2*sqrt(x**2+y**2+z**2)+1)*y*z/
     #        sqrt((sqrt(x**2+y**2+z**2)+2)/sqrt(x**2+y**2+z**2))/
     #        sqrt(x**2+y**2+z**2)**5

                    k33(i,j,k) = 2*(y**2*sqrt(x**2+y**2+z**2)-z**2+
     #        sqrt(x**2+y**2+z**2)*x**2-sqrt(x**2+y**2+z**2)*z**2)/
     #        sqrt((sqrt(x**2+y**2+z**2)+2)/sqrt(x**2+y**2+z**2))/
     #        sqrt(x**2+y**2+z**2)**5

                    end if

                    end do
                end do
            end do

            return
        end
c-----------------------------------------------------------------------
