c-----------------------------------------------------------------------
c    main.f
c    Apparent horizon locator via the numerical perturbation method.
c
c    Modification Date: Tue Jul 26 14:01:21 CDT 1994 MFH
c        . Uses Richard Matzner's extrinsic curvature code to generate 
c          extrinsic curvature given momenta and spins for the black holes.
c          The routine is supposed to move the holes, compute the initial
c          data and track the apparent horizon(s).
c
c    Modification Wed Nov  2 16:44:48 CST 1994
c        . Added in a area calculation routine that works with the final
c          converged answer. 
c
c-----------------------------------------------------------------------
        subroutine mijan_horizon_interface(
     &              g11,g12,g13,
     &              g22,g23,g33,
     &              k11,k12,k13,
     &              k22,k23,k33,
     &              grx,gry,grz,
     &              msk,nx,ny,nz,nst,nsp,hflag,bhcen,
     &              final_area,final_ct,final_cp,
     &              bh_mass,bh_radius,flag)

        implicit none

        include 'horizon_parameters.h'

        integer    NST, NSP, NEQN, NRED

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        integer  itmpl(NSTM,NSPM,4,4,2)

        real*8  bh_mass,bh_radius
        real*8    S(NSTM,NSPM,3), Sc(NSTM,NSPM,3),Sf(NSTM,NSPM),
     *            Sn(NSTM,NSPM,3), Sb(NSTM,NSPM,3,3),
     *            Sfmol(NSTM,NSPM,3,3,3),
     *            Sg(NSTM,NSPM,3,3),Sdg(NSTM,NSPM,3,3,3),
     *            Skt(NSTM,NSPM,3,3), Sk(NSTM,NSPM),
     *            AH(NSTM,NSPM), AHP(NSTM,NSPM),
     *            SP(NSTM,NSPM,3),ScP(NSTM,NSPM,3), SfP(NSTM,NSPM),
     *            SnP(NSTM,NSPM,3),SbP(NSTM,NSPM,3,3),
     *         SgP(NSTM,NSPM,3,3),SdgP(NSTM,NSPM,3,3,3),
     *         SktP(NSTM,NSPM,3,3), SkP(NSTM,NSPM),
     *         SgU(NSTM,NSPM,3,3),SdgU(NSTM,NSPM,3,3,3),
     *         SktU(NSTM,NSPM,3,3), SkU(NSTM,NSPM)

        
        real*8    divh,divh2, xc, yc, zc,final_area

c        Cartesian Metric data

                integer level
        integer    NX, NY, NZ, order, flag
        integer hflag

        real*8 bhcen(3)
        real*8    GrX(NX),GrY(NY),GrZ(NZ),msk(NX,NY,NZ),
     .        g11(NX,NY,NZ), g12(NX,NY,NZ),g13(NX,NY,NZ),
     .        g22(NX,NY,NZ), g23(NX,NY,NZ),g33(NX,NY,NZ),
     .        k11(NX,NY,NZ), k12(NX,NY,NZ),k13(NX,NY,NZ),
     .        k22(NX,NY,NZ), k23(NX,NY,NZ),k33(NX,NY,NZ)

c        Temporaray arrays

        real*8    Intg(NSPM, NSTM),
     .        tfunc(NSTM*NSPM),func(NSTM,NSPM), 
     .    lx(NSTM*NSPM), ly(NSTM*NSPM),lz(NSTM*NSPM)

        real*8    area, cp, ct, a,v, aharea, aharea2, PropArea1
        real*8 final_ct,final_cp

c        Other Cartesian grid parameters

        real*8    hx,hy,hz, fx, fy, fz
      real*8  rmax(3)
      real*8  rmin(3)
      real*8  rpos1(3)
      real*8  rpos2(3)
      real*8  mom1(3)
      real*8  mom2(3)
      real*8  spin1(3)
      real*8  spin2(3)
      real*8  mass1,mass2
      logical success



c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Mass = bh_mass
      radius = bh_radius
      if ( Mass .lt. 0.05 ) then
! searching for a black hole of zero mass
! won't work
        hflag = 0
        return
      end if
!      Mass = 1.0
!      radius = 2.0
      epsilon = 1.d-7
      newtst = 1.d-7
      ilucgstp = -256.0d0
      kappa0 = 0.0
      level = 1
      center(1) = bhcen(1)
      center(2) = bhcen(2)
      center(3) = bhcen(3)
      hx = grx(2) - grx(1)
      hy = gry(2) - gry(1)
      hz = grz(2) - grz(1)

!      write(0,*)NHOLES,NX,NST,NSP,Mass,radius,epsilon,newtst,
!     .       ilucgstp,flag,kappa0,level,center


101    format(I4)
102    format(F15.10)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Pi = 4.0d0 * atan(1.0e0)

        dtheta = Pi / float(NST - 1)
        dphi = 2.0d0 * Pi / float(NSP)

        h = 1.0d0/float(NST)
        divh = 0.5d0 / h
        divh2= 1.0d0 / h**2

        cntr(1,1) = 0.0d0
        cntr(1,2) = 0.0d0
        cntr(1,3) = 0.0d0

        ahole(1) = 0.5d0
            

        x0 = 0.0d0
        y0 = 0.0d0
        z0 = 0.0d0

      order = 4

c        epsilon = 1.0e-6

c                tinitial = time()

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c    Compute interpolation stencils on 2d mesh

        call CreateTemplate(itmpl,NST,NSP)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c    Either generates a two-sphere or reads 2-surface from file.
        if(flag .eq. 0)then
            call GenSur(S,Sc,SF,3,NST,NSP)
        else
         write(*,*) 'Using the previous horizon to initialize finder'
c        call ReadS2C('initial.dat',Sc,S,Sf,center,NST,NSP)
          call ReadS2C('solution.dat',Sc,S,Sf,center,NST,NSP)
        end if

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c        Load in the numerical metric data...

c            First load in the control file and figure out boundaries and such

c         call LoadBoundary(rmax,rmin,rpos1,rpos2,mom1,mom2,
c     $        spin1,spin2,mass1,mass2,hx,hy,hz,NX,NY,NZ)

c            Calculate what the coordinates based on the boundaries...

c         call CalcCoord(GrX,GrY,GrZ,hx,hy,hz,rmax,rmin,NX,NY,NZ)


c            Load in the metric and extrinsic curvature data from hdf files.

c         call HDFLoadMetric(GrG,GrKt,GrX,GrY,GrZ,temp,NX,NY,NZ)
c          call  HDFMetricData(nx, ny, nz, level,grx, gry, grz,
c     .           g11, g12, g13, g22, g23, g33,
c     .           k11, k12, k13, k22, k23, k33,
c     .           msk,hx,hy,hz)

c        Set up metric data locally

c            call GenMet(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
c     $        k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,NX,NY,NZ)
c            call GenKerr(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
c     $      k11,k12,k13,k22,k23,k33, msk,a,hx,hy,hz,Mass,NX,NY,NZ)
c            call GenKerrBst(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
c     $    k11,k12,k13,k22,k23,k33, msk,a,v,0.0,hx,hy,hz,NX,NY,NZ)
c            call GenEF2(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
c     $      k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,NX,NY,NZ)
c            call GenBstIso(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
c     $      k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,0.0,NX,NY,NZ)


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c    Evaluator Mode

        goto 110

        call GenSur(S,Sc,SF,3,NST,NSP)

c        Interpolate metric data to surface

        call metextp3(grx,gry,grz,g11,g12,g13,g22,g23,g33,
     .            k11,k12,k13,k22,k23,k33, Sg, Sdg, Skt,
     .            Sc, lx, ly, lz, msk,
     .            func, tfunc,order, NST,NSP,NX,NY,NZ,success)


c        call MetData(S,Sc,Sg,Sdg,Skt,Sk,NST,NSP)


        call DumpSg(Sg,Sdg,Sk,Skt,NST,NSP)

c        Create finite difference molecules

      call IntpMol(S,Sc,Sf,Sfmol,itmpl,
     *      1,1,NST,NSP,NST,NSP)

c        Finite difference

      call Derivatives(S,Sfmol,Sn,Sb,divh,divh2,1,1,
     $               NST,NSP,NST,NSP)

c        Evaluate apparent horizon equation

      call Evaluator(S,Sg,Sdg,Skt,Sk,AH,Sn,Sb,Sfmol,
     #                  1,1,NST,NSP,NST,NSP)


c     call WriteS2F('eval.dat',Sc,AH,NST,NSP)
      open(unit=1,file='eval1.dat')
      call Dp(Sc,AH,NST,NSP)
      close(unit=1)
c      open(unit=1,file='eval2.dat')
c      call Dp(Sc,AH,NST,NSP)
c      close(unit=1)
c      call Diff(Ah,AHP,NST,NSP)

110    continue

            

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c        goto 111

c        write(6,*)'Calling solver'
c                tinitial = time()
      call  SlvAH(S,Sc,Sf,Sfmol,Sn,Sb,AH,AHP,
     *      Sg,Sdg,Skt,Sk,itmpl,
     .      SP,ScP,SfP,SnP,SbP,SgP,SdgP,SktP,SkP,SgU,SdgU,SktU,SkU,
     .      g11, g12, g13, g22, g23, g33, k11, k12, k13, k22, k23, k33,
     .      GrX,GrY,GrZ,msk,func,tfunc,lx,ly,lz,
     .      order,hx,hy,hz,NX,NY,NZ,
     *      NST,NSP,success)
c                tfinal = time()

c                write(6,*)'time taken for find = ',tfinal-tinitial
c                write(2,*)'time taken for find = ',tfinal-tinitial


      if ( success ) then
        call WriteS2('solution.dat',Sc,S,Sf,center,NST,NSP)
        open(unit=1,file='solnAH.dat')
        call Dp(Sc,AH,NST,NSP)
        close(unit=1)

c        call kerrdump(S,Sc,a,center,NST,NSP,19)
c            open(unit=59,file='solncrt.dat')
c        call Dump3CoordCart(Sc,NST,NSP,59)
c        call Dp3(Sc,NST,NSP,59)
c        close(unit=59)




c        Projection approach
         final_area = 0.0d0
          area = AHArea(nst, nsp, dtheta, dphi,
     .                          S, Sc, Sn, Sg, Intg)

        final_area = final_area + area
        write(0,*)'Final Area of Apparent horizon (measure 1) = ',area
c        write(2,*)'Final Area of Apparent horizon (measure 1)= ',area
c        write(6,*)'Final Area of Apparent horizon (measure 1)= ',area 
c        Approx approach by considering tangents
        area =  PropArea1(nst, nsp, Sc, Sg)
        final_area = final_area + area
        write(0,*)'Final Area of Apparent horizon (measure 2) = ',area
c        write(2,*)'Final Area of Apparent horizon (measure 2)= ',area

c        Approach with rewriting line elements
        area = AHArea2(nst, nsp, dtheta, dphi,
     .                          S, Sc, Sn, Sg, Intg)

        final_area = final_area + area
        write(0,*)'Final Area of Apparent horizon (measure 3) = ',area
c        write(2,*)'Final Area of Apparent horizon (measure 3)= ',area
c        write(6,*)'Final Area of Apparent horizon (measure 3)= ',area   
c        Ordinary coordinate transformation
c  Set up the integrand

        call Horizon_Integrand(S,Sc,Sg,Intg,NST,NSP)

c  Integrate using the trapezoidal rule
        call SurfArea(Intg,area,dtheta,dphi,NST,NSP)
        final_area = final_area + area
        write(0,*)'Final Area of Apparent horizon (measure 4) = ',area
c        write(2,*)'Final Area of Apparent horizon (measure 4)= ',area

        call CalcCirc(S,Sc,Sg,cp,ct,NST,NSP)

        write(0,*)'Polar (theta) Circumference',ct*2.0e0
        write(0,*)'Equatorial (phi) Circumference',cp

c        write(2,*)'Polar (theta) Circumference',ct*2.0e0
c        write(2,*)'Equatorial (phi) Circumference',cp

c       call compbst(0.0,Sc,ScP,S,SP,Sf,SfP,NST,NSP)
        hflag = 1
        final_area = 0.25d0*final_area
        final_ct = ct*2.0d0
        final_cp = cp
      else
        hflag = 0
        final_area = 0.0d0
        final_ct = 0.0d0
        final_cp = 0.0d0
      end if





c        center(1) = 0.0e0
c        center(2) = 0.0e0
c        center(3) = 0.0e0

c        x0 = cntr(1,1)
c        y0 = cntr(1,2)
c        z0 = cntr(1,3)

c        radius = ahole(1)
c        call GenSur(S,Sc,SF,3,NST,NSP)

c        open(unit=59,file='offsph.dat')
c        call Dump3CoordCart(S,NST,NSP,59)
c        close(unit=59)
c        open(unit=59,file='offcrt.dat')
c        call Dump3CoordCart(Sc,NST,NSP,59)
c        close(unit=59)
111    continue

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


c            stop
        end subroutine
c=======================================================================
        subroutine Dp(Sc,AH,NST,NSP)

            implicit none

            integer    NST, NSP

            real*8    Sc(NST,NSP,3), AH(NST,NSP)

            integer i,j

            do i=1,NST
                do j=1,NSP

                    write(1,1122)i,j,AH(i,j)

                end do
                write(1,*)
            end do

1122    format(2I4,F20.15)

            return
        end
c-----------------------------------------------------------------------
      subroutine Diff(Ah,AHP,NST,NSP)
         implicit none

         integer  NST,NSP

         real*8  AH(NST,NSP), AHP(NST,NSP)

         integer  i,j

         do i=1,NST
            do j=1,NSP

               write(1,10101)i,j,AHP(i,j)-AH(i,j)

            enddo
            write(1,*)
         enddo

10101    format(2I5,F20.15)

         return
      end
c-----------------------------------------------------------------------
         subroutine Perturb(S,Sc,Sf,it,ip,NST,NSP)
         implicit none

         include 'horizon_parameters.h'

         integer  NSP, NST, it,ip

         real*8  S(NST,NSP,3)
         real*8  Sc(NST,NSP,3)
         real*8  Sf(NST,NSP)

         integer  i,j

         i=it
         j=ip
            if(i .eq. 1 .or. i .eq. NST)then
                do j=1,NSP
                Sf(i,j) = Sf(i,j) + epsilon

                S(i,j,1) = Sf(i,j)

                Sc(i,j,1)= S(i,j,1)*cos(S(i,j,3))*sin(S(i,j,2))
                Sc(i,j,2)= S(i,j,1)*sin(S(i,j,3))*sin(S(i,j,2))
                Sc(i,j,3)= S(i,j,1)*cos(S(i,j,2))
                end do
            else

                Sf(i,j) = Sf(i,j) + epsilon

                S(i,j,1) = Sf(i,j)

                Sc(i,j,1)= S(i,j,1)*cos(S(i,j,3))*sin(S(i,j,2))
                Sc(i,j,2)= S(i,j,1)*sin(S(i,j,3))*sin(S(i,j,2))
                Sc(i,j,3)= S(i,j,1)*cos(S(i,j,2))
            end if
         return
      end

c-----------------------------------------------------------------------
        subroutine DumpSg(Sg,Sdg,Sk,Skt,NST,NSP)
            implicit none

            integer    NST, NSP

            real*8    Sg(NST,NSP,3,3)
            real*8    Sdg(NST,NSP,3,3,3)
            real*8    Sk(NST,NSP)
            real*8    Skt(NST,NSP,3,3)

c        Local
    
        integer    i,j

        do i=1,NST
            do j=1,NSP
    
                write(92,*)Sg(i,j,1,1),Sg(i,j,1,2),Sg(i,j,1,3)
                write(92,*)Sg(i,j,2,2),Sg(i,j,2,3)
                write(92,*)Sg(i,j,3,3)

                write(93,*)Sdg(i,j,1,1,1),Sdg(i,j,2,2,1),Sdg(i,j,3,3,1)
                write(93,*)Sdg(i,j,1,1,2),Sdg(i,j,2,2,2),Sdg(i,j,3,3,2)
                write(93,*)Sdg(i,j,1,1,3),Sdg(i,j,2,2,3),Sdg(i,j,3,3,3)
                
                write(94,*)Sk(i,j),Sk(i,j),Sk(i,j)
                write(94,*)Sk(i,j),Sk(i,j)
                write(94,*)Sk(i,j)

                write(95,*)Skt(i,j,1,1),Skt(i,j,1,2),Skt(i,j,1,3)
                write(95,*)Skt(i,j,2,2),Skt(i,j,2,3)
                write(95,*)Skt(i,j,3,3)
            end do
        end do

        return
        end

c-----------------------------------------------------------------------

        subroutine Horizon_Integrand(S,Sc,Sg,Intg,NST,NSP)
            implicit none

            integer    NST, NSP

            real*8 S(NST,NSP,3), Sc(NST,NSP,3), 
     *            Sg(NST,NSP,3,3),Intg(NST,NSP)

c  Local Variables

            integer    i, j, k,m, p, q

            real*8        detg(NST,NSP)
            real*8        lambda(3,3), gg(3,3)

            do i = 1, NST
                do j = 1,NSP

c        X to r, theta, phi
                     lambda(1,1) = sin(S(i,j,2))*cos(S(i,j,3))
                     lambda(1,2) = S(i,j,1)*cos(S(i,j,2))*cos(S(i,j,3))
                     lambda(1,3) = -S(i,j,1)*sin(S(i,j,2))*sin(S(i,j,3))

c        Y to r, theta, phi
                     lambda(2,1) = sin(S(i,j,2))*sin(S(i,j,3))
                     lambda(2,2) = S(i,j,1)*cos(S(i,j,2))*sin(S(i,j,3))
                     lambda(2,3) = S(i,j,1)*sin(S(i,j,2))*cos(S(i,j,3))

c        Z to r, theta, phi
                     lambda(3,1) = cos(S(i,j,2))
                     lambda(3,2) = -S(i,j,1)*sin(S(i,j,2))
                     lambda(3,3) = 0.0e0

                     do k = 1, 3
                        do m = 1,3

                            gg(k,m) = 0.0e0

                            do p =1, 3
                                do q = 1,3
                                    gg(k,m) = gg(k,m) + 
     *            lambda(p,k)*lambda(q,m)*Sg(i,j,p,q)
                                end do
                            end do

                        end do
                     end do

                    

c                 detg(i,j) = Sg(i,j,1,1)*Sg(i,j,2,2)*Sg(i,j,3,3) -
c     #            Sg(i,j,1,1)*Sg(i,j,2,3)*Sg(i,j,3,2)  -
c     #            Sg(i,j,2,1)*Sg(i,j,1,2)*Sg(i,j,3,3)  +
c     #            Sg(i,j,2,1)*Sg(i,j,1,3)*Sg(i,j,3,2)  +
c     #            Sg(i,j,3,1)*Sg(i,j,1,2)*Sg(i,j,2,3)  -
c     #            Sg(i,j,3,1)*Sg(i,j,1,3)*Sg(i,j,2,2)

                    detg(i,j) = gg(2,2)*gg(3,3) - gg(2,3)*gg(3,2)


c                    Intg(i,j) = S(i,j,1)*S(i,j,1) * sin(S(i,j,2)) *
c     #                            sqrt(detg(i,j))
                    Intg(i,j) = sqrt(detg(i,j))

c                    write(98,*)Intg(i,j) - S(i,j,1)*S(i,j,1) * 
c     #                sin(S(i,j,2)) * Sg(i,j,1,1)

c                    write(99,*)Sg(i,j,1,1),Sg(i,j,2,2),Sg(i,j,3,3)

    
                end do
            end do

            return
        end

c-----------------------------------------------------------------------
      subroutine SurfArea(Itemp,II,dx,dy,NX,NY)
         implicit none

         integer  NX,NY

         real*8   Itemp(NX,NY), II,
     *            dx, dy

c     Local Variables

         integer  i,j

         real*8   dii

c     Initialize the total integral to be zero
         II = 0.0d0

c     Loop over patches on the surface and carry out the trapezoidal
c     method

         do i=1,NX-1
            do j=1,NY-1

               dii = dx*dy*(Itemp(i,j)+Itemp(i+1,j)+Itemp(i,j+1) +
     *                     Itemp(i+1,j+1))*0.25d0

               II = II + dii
            end do
c        Fix the last segment
            dii = dx*dy*(Itemp(i,NY) + Itemp(i+1,NY)+Itemp(i,1) +
     *            Itemp(i+1,1))*0.25d0
            II = II + dii
         end do




         return
      end
c-----------------------------------------------------------------------

        subroutine CalcCirc(S,Sc,Sg,Ce,Cp,NST,NSP)
            implicit none

            integer    NST, NSP
            
          real*8    S(NST,NSP,3), Sc(NST,NSP,3), Sg(NST,NSP,3,3), Ce, Cp

c        Local Variables

            integer    i, j, k, m

            real*8 dr, dx(3)

c        First calculate Ce (Equatorial Circumference. theta=pi/2)
        
            Ce = 0.0
            i = (NST/2)
            do j = 1, NSP-1

                dr = 0.0
                do k=1,3
                do m=1,3
                    dr = dr + Sg(i,j,k,m)*(Sc(i,j+1,k) - Sc(i,j,k))*
     .            (Sc(i,j+1,m) - Sc(i,j,m))

                end do
                enddo
                dr = sqrt(dr)
                Ce = Ce + dr
            end do

c        Patch last point b/w 2pi-dphi to 2pi

            dr = 0.0
            do k=1,3
                do m=1,3
                    dr = dr + Sg(i,NSP,k,m)*(Sc(i,1,k) - Sc(i,NSP,k))*
     .            (Sc(i,1,m) - Sc(i,NSP,m))

                end do
            enddo
            Ce = Ce + dr
                
c        Calculate Cp (Polar) phi = 0
            
            Cp = 0.0
            j = 1
            do i = 1, NST-1
                dr = 0.0
                do k=1,3
                do m=1,3
                    dr = dr + Sg(i,j,k,m)*(Sc(i+1,j,k) - Sc(i,j,k))*
     .            (Sc(i+1,j,m) - Sc(i,j,m))

                end do
                enddo
                dr = sqrt(dr)
                Cp = Cp + dr
            end do
        
            return
        end
c-----------------------------------------------------------------------    

c    Subroutine to locate a new center for the calculation.
c    The necessity for the is routine comes out of having problems with 
c    offset black holes that have different origins from coordinate
c    system used. Locates an average center point.

        subroutine LocateCenter2(S,Sc,newcntr,NST,NSP)
            implicit none

            include 'horizon_parameters.h'

            integer    NST, NSP

            real*8    S(NST,NSP,3), Sc(NST,NSP,3), 
     *                newcntr(3)

c    Local Variables

            integer    i, j, k

            real*8    ncntr(3), min(3),max(3)

            ncntr(1) = 0.0e0
            ncntr(2) = 0.0e0
            ncntr(3) = 0.0e0



            min(1) = 99999.0
            min(2) = 99999.0
            min(3) = 99999.0
            max(1) = -99999.0
            max(2) = -99999.0
            max(3) = -99999.0

            do i = 1, NST
                do j = 1, NSP
                    if(Sc(i,j,1) .lt. min(1))min(1) = Sc(i,j,1)
                    if(Sc(i,j,1) .gt. max(1))max(1) = Sc(i,j,1)

                    if(Sc(i,j,2) .lt. min(2))min(2) = Sc(i,j,2)
                    if(Sc(i,j,2) .gt. max(2))max(2) = Sc(i,j,2)

                    if(Sc(i,j,3) .lt. min(3))min(3) = Sc(i,j,3)
                    if(Sc(i,j,3) .gt. max(3))max(3) = Sc(i,j,3)
                end do
            end do


c                newcntr(1) = newcntr(1) + ncntr(1)
c                newcntr(2) = newcntr(2) + ncntr(2)
c                newcntr(3) = newcntr(3) + ncntr(3)

c                write(0,*)min(1),max(1)
c                write(0,*)min(2),max(2)
c                write(0,*)min(3),max(3)

                newcntr(1) =( max(1) + min(1) )/2.0
                newcntr(2) = (max(2) + min(2) )/2.0
                newcntr(3) = (max(3) + min(3) )/2.0
            return
        end
c=======================================================================
c=======================================================================
c        
        subroutine findrad(radius,S,NST,NSP)
            implicit none
    
            integer    NST, NSP

            real*8    radius, S(NST,NSP,3)

c        Local variables

            integer    i, j, k

            real*8    radmax,radmin

            radmin = 99999.9
            do i =1 ,NST
                do j = 1, NSP
                    if(S(i,j,1) .lt. radmin)then
                        radmin = S(i,j,1)
                    end if
                    if(S(i,j,1) .gt. radmax)then
                        radmax = S(i,j,1)
                    end if
                end do
            end do

            radius = 0.5*(radmin + radmax)
            return
        end
c-----------------------------------------------------------------------
        subroutine compbst(t,Sc,ScP,S,SP,Sf,SfP,NST,NSP)
            implicit none

            include 'horizon_parameters.h'

            integer    NST, NSP

            real*8        Sc(NST,NSP,3), S(NST,NSP,3), Sf(NST,NSP),
     .            ScP(NST,NSP,3), SP(NST,NSP,3), SfP(NST,NSP), t

c        Local Variables

            integer    i, j, k 

            real*8    gamma, dx,dy,dz, drad

            real*8    min,max,ccx,ccy,ccz, mskrad,mx,my,mz,v

         open(unit=1,file='bstiso.dat')
         read(1,*)min,max
         read(1,*)ccx,ccy,ccz

c        mask origin and radius
         read(1,*)mskrad
         read(1,*)mx,my,mz
         read(1,*)v
         close(unit=1)


            open(unit=56,file='cmprad.dat')
            open(unit=55,file='cmpcrt.dat')

c        First generate a sphere of radius 1.0 centered at center(3)

            radius = 1.0
            
            call GenSur(SP,ScP,SFP,3,NST,NSP)

            gamma  = 1.0 / sqrt(1.0 - v**2)
            do i = 1, NST
                do j = 1, NSP
                    dx =  (ScP(i,j,1) - v*t)*gamma
                    dy =  ScP(i,j,2)
                    dz =  ScP(i,j,3)
            
                    drad = 1.0 - sqrt(dx**2 + dy**2 + dz**2)
                    
                    write(55,231)dx,dy,dz
                    write(56,232)drad
                end do
                write(55,*)
            end do

            close(unit=55)
            close(unit=56)

            
231    format(3F20.15)
232    format(F20.15)
            return
        end subroutine
