c-----------------------------------------------------------------------
c    metricdata.f
c    Created Mon Jun 13 17:49:36 CDT 1994 MFH
c
c    This file contains routines that will interpolate or set metric and
c    extrinsic curvature data.
c-----------------------------------------------------------------------
        subroutine MetIntp(S,Sc,Sg,Sdg,Sk,Skt,GrX,GrY,GrZ,GrG,GrDG,
     *                GrKt,hx,hy,hz,th0,thF,ph0,phF,
     *                NX,NY,NZ,NST,NSP)
            implicit none

            include 'horizon_parameters.h' 

            integer    NX,NY,NZ,NST,NSP, th0,thF,ph0,phF

c    Surface Grid Functions

            real*8    S(NST,NSP,3),Sc(NST,NSP,3),
     *                Sg(NST,NSP,3,3), Sdg(NST,NSP,3,3,3),
     *                Sk(NST,NSP), Skt(NST,NSP,3,3)

c    Cartesian Grid functions

            real*8    GrX(NX), GrY(NY), GrZ(NZ), 
     *                GrG(NX,NY,NZ,3,3), GrDG(NX,NY,NZ,3,3,3),
     *                GrKt(NX,NY,NZ,3,3)

            real*8    hx,hy,hz

c        Local Variables

            integer    i,j,k,m,p

            real*8    ix(3),t1(3,3),t2(3,3,3),t3
            real*8    invg(NST,NSP,3,3),detg(NST,NSP)
            real*8    Intp4Sc3

c        Loop over points that are to be interpolated over
            do i=th0,thF
                do j=ph0,phF
                    
c            Load in the coordinates of the point
                    ix(1) = Sc(i,j,1)    
                    ix(2) = Sc(i,j,2)    
                    ix(3) = Sc(i,j,3)    

c            Interpolate the metric
                    call Int2Tsr3(GrG,GrX,GrY,GrZ,ix,t1,hx,hy,hz,
     #               3,3,NX,NY,NZ)

c            Copy over the metric 

                    do k=1,3
                        do m=1,3
                            Sg(i,j,k,m) = t1(k,m)
                        end do
                    end do
c                    write(71,*)ix(1),ix(2),ix(3),Sg(i,j,1,1)
c                    write(72,*)ix(1),ix(2),ix(3),Sg(i,j,2,2)
c                    write(73,*)ix(1),ix(2),ix(3),Sg(i,j,3,3)
                end do
            end do
c            write(0,*)'Got here'

            do i=th0,thF
                do j=ph0,phF
                    
c            Load in the coordinates of the point
                    ix(1) = Sc(i,j,1)    
                    ix(2) = Sc(i,j,2)    
                    ix(3) = Sc(i,j,3)    

c            Interpolate the derivatives of the metric

                    call Int3Tsr3(GrDG,GrX,GrY,GrZ,ix,t2,hx,hy,hz,
     #               3,3,3,NX,NY,NZ)
        
c            Copy over the derivatives of the metric 

                    do k=1,3
                        do m=1,3
                            Sdg(i,j,k,m,1) = t2(k,m,1)
                            Sdg(i,j,k,m,2) = t2(k,m,2)
                            Sdg(i,j,k,m,3) = t2(k,m,3)
                        end do
                    end do
                
c            Interpolate the Extrinsic curvature

                    call Int2Tsr3(GrKt,GrX,GrY,GrZ,ix,t1,hx,hy,hz,
     #               3,3,NX,NY,NZ)                

c            Copy over the components of the extrinsic curvature tensor

                    do k=1,3
                        do m=1,3

                            Skt(i,j,k,m) = t1(k,m)
        
                        end do
                    end do

c        Calculate the inverse of the 3-metric first.

c        Determinant of the metric


             detg(i,j) = Sg(i,j,1,1)*Sg(i,j,2,2)*Sg(i,j,3,3) -
     #                Sg(i,j,1,1)*Sg(i,j,2,3)*Sg(i,j,3,2)  -
     #                Sg(i,j,2,1)*Sg(i,j,1,2)*Sg(i,j,3,3)  +
     #                Sg(i,j,2,1)*Sg(i,j,1,3)*Sg(i,j,3,2)  +
     #                Sg(i,j,3,1)*Sg(i,j,1,2)*Sg(i,j,2,3)  -
     #                Sg(i,j,3,1)*Sg(i,j,1,3)*Sg(i,j,2,2)

c        Calculate the components

               invg(i,j,1,1) = (Sg(i,j,2,2)*Sg(i,j,3,3) -
     #                  Sg(i,j,2,3)*Sg(i,j,3,2))/detg(i,j)

               invg(i,j,1,2) = (Sg(i,j,1,3)*Sg(i,j,3,2) -
     #                  Sg(i,j,1,2)*Sg(i,j,3,3))/detg(i,j)

               invg(i,j,1,3) = (Sg(i,j,1,2)*Sg(i,j,2,3) -
     #                  Sg(i,j,1,3)*Sg(i,j,2,2))/detg(i,j)

               invg(i,j,2,1) = (Sg(i,j,2,3)*Sg(i,j,3,1) -
     #                  Sg(i,j,2,1)*Sg(i,j,3,3))/detg(i,j)

               invg(i,j,2,2) = (Sg(i,j,1,1)*Sg(i,j,3,3) -
     #               Sg(i,j,1,3)*Sg(i,j,3,1))/detg(i,j)

               invg(i,j,2,3) = (Sg(i,j,1,3)*Sg(i,j,2,1) -
     #               Sg(i,j,1,1)*Sg(i,j,2,3))/detg(i,j)

               invg(i,j,3,1) = (Sg(i,j,2,1)*Sg(i,j,3,2) -
     #               Sg(i,j,2,2)*Sg(i,j,3,1))/detg(i,j)

               invg(i,j,3,2) = (Sg(i,j,1,2)*Sg(i,j,3,1) -
     #               Sg(i,j,1,1)*Sg(i,j,3,2))/detg(i,j)

               invg(i,j,3,3) = (Sg(i,j,1,1)*Sg(i,j,2,2) -
     #               Sg(i,j,1,2)*Sg(i,j,2,1))/detg(i,j)

                    Sk(i,j)= 0.0e0
                    do k=1,3
                        do m=1,3

                          Sk(i,j) = Sk(i,j) + invg(i,j,k,m)*Skt(i,j,k,m)

                        end do
                    end do

                end do
            end do


            return
        end
c-----------------------------------------------------------------------

            

