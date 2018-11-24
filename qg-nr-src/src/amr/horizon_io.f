c-----------------------------------------------------------------------
c io.f
c
c Creation Date: Fri Mar 25 16:47:48 CST 1994 
c
        subroutine Dump2Arr(Arr,N1,N2,Unum)
            implicit none

            integer    N1,N2,Unum

            real*8    Arr(N1,N2)

c            Local Variables
        
            integer    i,j,k

            do i=1,N1
                do j=1,N2

                        write(Unum,*)i,j,Arr(i,j)
    
                enddo
                write(Unum,*)
            enddo

            return
        end
c-----------------------------------------------------------------------
        subroutine Dump3Arr(Arr,N1,N2,N3,Unum)
            implicit none

            integer    N1,N2,N3,Unum

            real*8    Arr(N1,N2,N3)

c            Local Variables
        
            integer    i,j,k

            do i=1,N1
                do j=1,N2
                    do k=1,N3
                        write(Unum,101)i,j,k,Arr(i,j,k)
                    enddo    
                enddo
                write(Unum,*)
            enddo
101        format(3I4,F10.6)
            return
        end
c-----------------------------------------------------------------------
        subroutine Dump3CoordCart(Arr,N1,N2,Unum)
            implicit none

            integer    N1,N2,Unum

            real*8    Arr(N1,N2,3)

c        Local Variables

            integer    i,j,k

            do i=1,N1
                do j=1,N2

                    write(Unum,*)Arr(i,j,1),Arr(i,j,2),Arr(i,j,3)

                end do
                write(Unum,*)
            end do

            return
        end
c-----------------------------------------------------------------------
        subroutine Dump3Arr3Col(Arr,N1,N2,N3,Unum)
            implicit none

            integer    N1,N2,N3,Unum

            real*8    Arr(N1,N2,N3)

c        Local variables

            integer    i,j,k

            if(N3 .ne. 3)then
                write(0,*)'**FILE NOT WRITTEN**'
             write(0,*)'Cannot use Dump3Arr3Col with N3 not equal to 3.'
                write(0,*)'    Use Dump3Arr instead.'
            else
                do i=1,N1
                    do j=1,N1

                     write(Unum,201)i,j,Arr(i,j,1),Arr(i,j,2),Arr(i,j,3)
            
                    end do 
                    end do
            end if

201        format(2I4,3F10.6)
            return
        end
c-----------------------------------------------------------------------
        subroutine Dump4Arr(Arr,N1,N2,N3,N4,Unum)
            implicit none

            integer    N1,N2,N3,N4,Unum

            real*8    Arr(N1,N2,N3,N4)

c            Local Variables
        
            integer    i,j,k,m

            do i=1,N1
                do j=1,N2
                    do k=1,N3
                        do m=1,N4
                            write(Unum,101)i,j,k,m,Arr(i,j,k,m)
                        enddo    
                    enddo    
                enddo
                write(Unum,*)
            enddo
101        format(4I4,F10.6)
            return
        end
c-----------------------------------------------------------------------
        subroutine ReadS2Cart(fname,Array,NST,NSP)
            implicit none

            integer    NST,NSP

            real*8    Array(NST,NSP,3)
    
            character*64    fname

c        Local Variables

            integer    i,j,k

            open(unit=1,file=fname)
        
            read(1,*)NST,NSP

            k=1
            do i=1,NST
                do j=1,NSP

                    read(1,*)Array(i,j,1),Array(i,j,2),Array(i,j,3)

        
                end do
            end do

            close(unit=1)
            return
        end
c-----------------------------------------------------------------------
        subroutine    DumpS2CartFunc(fname,Coord,Array,NST,NSP)
            implicit none
            
            integer    NST,NSP

            real*8    Coord(NST,NSP,3)
            real*8    Array(NST,NSP)

            character*64    fname

c        Local Variables

            integer    i,j,k

            open(unit=11,file=fname)
    
            write(11,*)NST,NSP
            do i=1,NST
                do j=1,NSP

                    write(11,301)Coord(i,j,1),Coord(i,j,2),Coord(i,j,3),
     *            Array(i,j)

                end do
            end do

301    format(4F15.10)

            return
        end
                    

c-----------------------------------------------------------------------
c    This routine will copy over information from an one dimensional 
c    array to a 3D array
      subroutine Tr1D3D(oned,threed,N1,N2,N3)
         implicit none

         integer  N1,N2,N3

         real*8   oned(1),threed(N1,N2,N3)

c     Local Variables

         integer  i,j,k,m

         k=1
c            write(0,*)N1,N2,N3
         do i=1,N1
            do j=1,N2
                    do m=0,N3-1

                        threed(i,j,m+1) = oned(k+m)
c                        threed(i,j,1) = oned(k)
c                        threed(i,j,2) = oned(k+1)
c                        threed(i,j,3) = oned(k+2)



                    end do
                    k = k + N3
            end do
         end do

         return
      end
c-----------------------------------------------------------------------
c    The following routine reads in spherical surface function data.
c    S2FC stands for S2 function in cartesian coordinates.
c    There will be a S2FS for spherical coordinates.
        subroutine ReadS2FC(fname,Coord,Func,NST,NSP)
            implicit none

            integer    NST,NSP

            real*8    Coord(1), Func(1)

            character*64    fname

c        Local Variables

            integer    i,j,ck,fk

            open(unit=99,file=fname)

            read(99,*)NST,NSP

            if(NST .eq. 0 .or. NSP .eq. 0)then
                write(0,*)'Error: ReadS2FC'
                write(0,*)'Zero size in file',fname
                stop
            end if

c    ck:= coordinate array counter. Updated by 3 each time
c    fk:= function array counter. Updated by 1 each time
            ck = 1
            fk = 1
            do i=1,NST
                do j=1,NSP

                    read(99,*)Coord(ck),Coord(ck+1),Coord(ck+2),Func(fk)

                    ck = ck + 3
                    fk = fk + 1

        
                end do
            end do


            close(unit = 99)
            return
        end
c-----------------------------------------------------------------------
c    The following routine reads in spherical surface function data.
        subroutine ReadS2F(fname,Coord,Func,NST,NSP)
            implicit none

            integer    NST,NSP

            real*8    Coord(NST,NSP,3), Func(NST,NSP)

            character*64    fname

c        Local Variables

            integer    i,j,ck,fk

            open(unit=99,file=fname)

c            write(6,*)NST,NSP
c            read(99,*)NST,NSP

            if(NST .eq. 0 .or. NSP .eq. 0)then
                write(0,*)'Error: ReadS2FC'
                write(0,*)'Zero size in file',fname
                stop
            end if

c    ck:= coordinate array counter. Updated by 3 each time
c    fk:= function array counter. Updated by 1 each time
            ck = 1
            fk = 1
            do i=1,NST
                do j=1,NSP
        
            read(99,*)Coord(i,j,1),Coord(i,j,2),
     *                    Coord(i,j,3),Func(i,j)
c                    ck = ck + 3
c                    fk = fk + 1
                end do
            end do
c1002        format(4F15.10)

        


            close(unit = 99)
            return
        end
c-----------------------------------------------------------------------
c    The following routine reads in the cartesian coordinates of 
c    a S2 surface.
c    S2C stands for S2 coordinates in cartesian coordinates.
c    There will be a S2S for spherical coordinates.
        subroutine ReadS2C(fname,Coord1,Coord2,Func,center,NST,NSP)
            implicit none

            integer    NST,NSP

            real*8    Coord1(NST,NSP,3),Coord2(NST,NSP,3),
     *                Func(NST,NSP), center(3)

            character*64    fname

c        Local Variables

            integer    i,j,ck


            open(unit=99,file=fname)

            read(99,*)NST,NSP
            read(99,*)center(1),center(2),center(3)

            if(NST .eq. 0 .or. NSP .eq. 0)then
                write(0,*)'Error: ReadS2C'
                write(0,*)'Zero size in file',fname
                stop
            end if


c        First write the cartesian coordinates...

            do i=1,NST
                do j=1,NSP

                    read(99,*)Coord1(i,j,1),Coord1(i,j,2),
     *                        Coord1(i,j,3)

c                    ck = ck + 3

        
                end do
            end do


c        now write the spherical coordinates

            do i=1,NST
                do j=1,NSP

                    read(99,*)Coord2(i,j,1),Coord2(i,j,2),
     *                        Coord2(i,j,3)


        
                end do
            end do

c        Write the function

            do i=1,NST
                do j=1,NSP

                    read(99,*)Func(i,j)

c                    ck = ck + 3

        
                end do
            end do




            close(unit = 99)
            return
        end
c-----------------------------------------------------------------------
c    The following routine reads in the cartesian coordinates of 
c    a S2 surface.
c    S2C stands for S2 coordinates in cartesian coordinates.
c    There will be a S2S for spherical coordinates.
        subroutine WriteS2(fname,Coord1,Coord2,Func,center,NST,NSP)
            implicit none

            integer    NST,NSP

            real*8    Coord1(NST,NSP,3), Coord2(NST,NSP,3), 
     .                Func(NST,NSP), center(3)

            character*64    fname

c        Local Variables

            integer    i,j,ck


            if(NST .eq. 0 .or. NSP .eq. 0)then
                write(0,*)'Error: WriteS2'
                write(0,*)'Zero size in file',fname
                stop
            end if

            open(unit=99,file=fname)

            write(99,*)NST,NSP
            write(99,*)center(1),center(2),center(3)


c        First write the cartesian coordinates...

            do i=1,NST
                do j=1,NSP

                    write(99,*)Coord1(i,j,1),Coord1(i,j,2),
     *                        Coord1(i,j,3)

c                    ck = ck + 3

        
                end do
            end do


c        now write the spherical coordinates

            do i=1,NST
                do j=1,NSP

                    write(99,*)Coord2(i,j,1),Coord2(i,j,2),
     *                        Coord2(i,j,3)


        
                end do
            end do

c        Write the function

            do i=1,NST
                do j=1,NSP

                    write(99,*)Func(i,j)

c                    ck = ck + 3

        
                end do
            end do


            close(unit = 99)
            return
        end
c
c-----------------------------------------------------------------------
        subroutine WriteS2F(fname,Coord,Func,NST,NSP)
            implicit none

            integer    NST,NSP

            real*8    Coord(NST,NSP,3), Func(NST,NSP)

            character*64    fname

c        Local Variables

            integer    i,j,ck,fk

            open(unit=99,file=fname)


            if(NST .eq. 0 .or. NSP .eq. 0)then
                write(0,*)'Error: ReadS2FC'
                write(0,*)'Zero size in file',fname
                stop
            end if

            write(99,*)NST,NSP
c    ck:= coordinate array counter. Updated by 3 each time
c    fk:= function array counter. Updated by 1 each time
            ck = 1
            fk = 1
            do i=1,NST
                do j=1,NSP

            write(99,100)Coord(i,j,1),Coord(i,j,2),Coord(i,j,3),
     *            Func(i,j)

c                    ck = ck + 3
c                    fk = fk + 1

        
                end do
            end do


            close(unit = 99)
100        format(4F15.10)
            return
        end
c-----------------------------------------------------------------------
c    ReadSSph
c    
c    Reads in a a file containing Spherical coordinate surface information.
c    Specifically for the apparent horizon locator
        subroutine ReadSSph(fname,S,Sf,Sc,NST,NSP)
            implicit none

            integer    NST, NSP

            real*8    S(NST,NSP,3), Sf(NST,NSP),Sc(NST,NSP,3)

            character*64    fname

c        Local Variables

            integer    i,j, it,ip

            open(unit = 98, file = fname)

            read(98,*)it,ip

            if(it .ne. NST .or. ip .ne. NSP)then
                write(0,*)'Error in ReadSSph: Sizes not compatible'
                stop
            end if

            do i=1,NST
                do j=1,NSP

                    read(98,*)S(i,j,1),S(i,j,2),S(i,j,3)
                    Sf(i,j) = S(i,j,1)

                    Sc(i,j,1) = S(i,j,1)*sin(S(i,j,2))*cos(S(i,j,3))
                    Sc(i,j,2) = S(i,j,1)*sin(S(i,j,2))*sin(S(i,j,3))
                    Sc(i,j,3) = S(i,j,1)*cos(S(i,j,2))

                end do
            end do
            close(unit = 98)


            return
        end
c-----------------------------------------------------------------------
                    
