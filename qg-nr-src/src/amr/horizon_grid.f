c-----------------------------------------------------------------------
c grid.f 
c Written by Mijan Huq 
c Univ. of Texas at Austin
c
c Routines to carry out extrapolation/interpolation of 3D grid functions
c for the purposes of an apparent horizon tracker.
c
c Preliminary use of Scott Klasky's black box interpolation routines.
c
c Modification Tue Aug  8 10:06:31 CDT 1995 MFH
c        Put in new black box interpolation routines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c    Routine to interface to apparent horizon locator. 
c
c    Inputs:
c        g11,g12,g13,g22,g23,g33 of size (NX,NY,NZ)
c        k11,k12,k13,k22,k23,k33 of size (NX,NY,NZ)
c
c    Outputs: 
c        Sg(NT,NP,3,3), Sdg(NT,NP,3), Skt(NT,NP,3,3)
c
c    Work Arrays necessary:
c        wkspc(soc(order)*(soc(order) + 4)*NT*NP
c        func(NT,NP), tfunc(NT*NP)

        subroutine metextp3(grx,gry,grz,g11,g12,g13,g22,g23,g33,
     .                k11,k12,k13,k22,k23,k33, Sg, Sdg, Skt,
     .                Sc, lx, ly, lz, msk,
     .                func, tfunc,order, NT,NP,NX,NY,NZ,success)

            implicit none

            include 'horizon_terms.inc'

            integer    NX, NY, NZ, NT, NP, order

            real*8    g11(NX,NY,NZ), g12(NX,NY,NZ), g13(NX,NY,NZ),
     .            g22(NX,NY,NZ), g23(NX,NY,NZ), g33(NX,NY,NZ),
     .            k11(NX,NY,NZ), k12(NX,NY,NZ), k13(NX,NY,NZ),
     .            k22(NX,NY,NZ), k23(NX,NY,NZ), k33(NX,NY,NZ),
     .            grx(NX), gry(NY),grz(NZ)

            real*8    Sg(NT,NP,3,3), Sdg(NT,NP,3,3,3), Skt(NT,NP,3,3),
     .                    Sc(NT,NP,3)

            real*8    func(NT,NP), tfunc(NT*NP)    

            real*8 lx(NT*NP), ly(NT*NP), lz(NT*NP), msk(NX,NY,NZ)
            logical success

c        Local variables

            integer    i, j, k, ierr, d4lintxyzcd, d4lintxyzc,
     .        d2lintxyzcd, d2lintxyzc

            integer iwk(20*NT*NP),swk(5*NT*NP)
            integer iwkx(20*NT*NP),swkx(5*NT*NP)
            integer iwky(20*NT*NP),swky(5*NT*NP)
            integer iwkz(20*NT*NP),swkz(5*NT*NP)

            integer diditx,didity,diditz,diditf

            diditf = 0
            diditx = 0
            didity = 0
            diditz = 0

c        First begin by setting up the coordinate list to extrapolate to

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    lx(k) = Sc(i,j,1)
                    ly(k) = Sc(i,j,2)
                    lz(k) = Sc(i,j,3)
c                    write(6,*)lx(k), ly(k), lz(k)
                end do
            end do

c                write(6,*)'****metextp begin****'
c            do i = 1, NX
c                write(6,*)grx(i), gry(i), grz(i)
c            end do
c                write(6,*)'****metextp end****'



c        Evaluate g11

            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g11,tfunc,msk,iwk,swk,diditf,' ')
c         write(6,*)'ERROR CODE g11',ierr
            if (ierr .gt. 0) then
c             the horizon finder has failed -- exit out cleanly
              success = .false.
              return
            end if
            
            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sg(i,j,1,1) = tfunc(k)
                end do
            end do

            
c        Evaluate g12
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g12,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sg(i,j,1,2) = tfunc(k)
                    Sg(i,j,2,1) = tfunc(k)
                end do
            end do

c        Evaluate g13
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g13,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sg(i,j,1,3) = tfunc(k)
                    Sg(i,j,3,1) = tfunc(k)
                end do
            end do


            
c        Evaluate g22
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g22,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sg(i,j,2,2) = tfunc(k)
                end do
            end do


c        Evaluate g23
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g23,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sg(i,j,2,3) = tfunc(k)
                    Sg(i,j,3,2) = tfunc(k)
                end do
            end do


            
c        Evaluate g33
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g33,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sg(i,j,3,3) = tfunc(k)
                end do
            end do


c        Derivatives of the metric 
c
c    Start with g11_i

c        Evaluate 11i
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g11,tfunc,msk,iwkx,swkx,diditx,'x' )

            

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,1,1) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g11,tfunc,msk,iwky,swky,didity,'y' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,1,2) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g11,tfunc,msk,iwkz,swkz,diditz,'z' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,1,3) = tfunc(k)
                end do
            end do
                
c        Evaluate 12i
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g12,tfunc,msk,iwkx,swkx,diditx,'x' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,2,1) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g12,tfunc,msk,iwky,swky,didity,'y' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,2,2) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g12,tfunc,msk,iwkz,swkz,diditz,'z' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,2,3) = tfunc(k)
                end do
            end do

c        Evaluate 13i
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g13,tfunc,msk,iwkx,swkx,diditx,'x' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,3,1) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g13,tfunc,msk,iwky,swky,didity,'y' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,3,2) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g13,tfunc,msk,iwkz,swkz,diditz,'z' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,1,3,3) = tfunc(k)
                end do
            end do

c        Evaluate 22i
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g22,tfunc,msk,iwkx,swkx,diditx,'x' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,2,2,1) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g22,tfunc,msk,iwky,swky,didity,'y' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,2,2,2) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g22,tfunc,msk,iwkz,swkz,diditz,'z' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,2,2,3) = tfunc(k)
                end do
            end do
c        Evaluate 23i
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g23,tfunc,msk,iwkx,swkx,diditx,'x' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,2,3,1) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g23,tfunc,msk,iwky,swky,didity,'y' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,2,3,2) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g23,tfunc,msk,iwkz,swkz,diditz,'z' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,2,3,3) = tfunc(k)
                end do
            end do
c        Evaluate 33i
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g33,tfunc,msk,iwkx,swkx,diditx,'x' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,3,3,1) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g33,tfunc,msk,iwky,swky,didity,'y' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,3,3,2) = tfunc(k)
                end do
            end do
                
            ierr = d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,g33,tfunc,msk,iwkz,swkz,diditz,'z' )

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Sdg(i,j,3,3,3) = tfunc(k)
                end do
            end do

            do k = 1, 3
                do i = 1, NT
                    do j = 1, NP
                        Sdg(i,j,2,1,k) = Sdg(i,j,1,2,k)
                        Sdg(i,j,3,1,k) = Sdg(i,j,1,3,k)
                        Sdg(i,j,3,2,k) = Sdg(i,j,2,3,k)
                    end do
                end do
            end do

c        Next the Extrinsic curvature components

c        Evaluate k11

            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,k11,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Skt(i,j,1,1) = tfunc(k)
                end do
            end do

            
c        Evaluate k12
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,k12,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Skt(i,j,1,2) = tfunc(k)
                    Skt(i,j,2,1) = tfunc(k)
                end do
            end do

c        Evaluate k13
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,k13,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Skt(i,j,1,3) = tfunc(k)
                    Skt(i,j,3,1) = tfunc(k)
                end do
            end do


            
c        Evaluate k22
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,k22,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Skt(i,j,2,2) = tfunc(k)
                end do
            end do


c        Evaluate k23
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,k23,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Skt(i,j,2,3) = tfunc(k)
                    Skt(i,j,3,2) = tfunc(k)
                end do
            end do


            
c        Evaluate k33
            ierr= d4lintxyzcd( NX,NY,NZ,grx,gry,grz,NT*NP,
     .                    lx,ly,lz,k33,tfunc,msk,iwk,swk,diditf,' ')

            k = 0 
            do i = 1, NT
                do j = 1, NP
                    k = k + 1
                    Skt(i,j,3,3) = tfunc(k)
                end do
            end do


            
            success = .true.
            return
            end
c-----------------------------------------------------------------------
