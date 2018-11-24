c=======================================================================
c jacobian.f
c Created Mon May 23 13:40:31 CDT 1994 MFH
c
c    Mijan Huq
c    Center for Relativity,
c    University of Texas at Austin.
c
c Version 2
c    Modified June 30th, 1994 MFH
c     . Put in EvalTerms from version 1 since that was conv. tested.
c   . Making use of precomputed interpolation coefficients.
c   . Precompute the interpolation neighborhoods.
c   . Poles still tricky. Since poles supply the "boundaries" in our
c   . problem there would be a singular jacobian from that if not correct.
c    Modified July 11th, 1994 MFH
c     . Added in r dependence.
c    Modified Tue Jul 26 14:08:25 CDT 1994 MFH
c     .    Changed routines to numerical perturbation routines. Analytic
c     . jacobian techniques still under development under ~/Projects/Aljac
c    Modified Tue Aug  8 13:40:02 CDT 1995 MFH
c     . Changed metexp to use new interpolation routines
c
c=======================================================================
c    Routine: PAnlJac
c    Created Mon Jul 18 18:24:57 CDT 1994 Mijan Huq
c    
c    Pseudo-Analytic Jacobian generation. Generated the functional 
c    derivatives of the first and second derivatives via a numerical
c    perturbation approach like the completely numerical approach done
c    previously for the apparent horizon equation.
c
c    Input:
c   S(*)     := Spherical coordinates for surface wrt to center
c   SP(*)    := Perturbed Spherical coordinates for surface wrt to center
c   Sc(*)    := Cartesian coordinates for surface
c   ScP(*)   := Perturbed Cartesian coordinates for surface
c   Sf(*)    := Surface function values on mesh
c   SfP(*)   := Perturbed Surface function values on mesh
c   Sfmol(*) := Finite difference molecules at each point.
c   Sn(*)    := First derivatives of Sf at each grid point.
c   SnP(*)   := Perturbed First derivatives of Sf at each grid point.
c   Sb(*)    := Second and mixed derivatives at each grid point.
c   SbP(*)   := Perturbed Second and mixed derivatives at each grid point.
c   Sg(*)    := Metric at surface grid points.
c   Sdg(*)   := Derivatives of metric at surface grid points.
c   itmpl(*) := Template for interpolation in surface.
c   NST/NSP  := Number of grid points in theta-phi mesh.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        subroutine PAnlJac(jac,S,Sc,Sf,Sfmol,Sn,Sb,Sg,Sdg,Skt,Sk,
     .        SP,ScP,SfP,SnP,SbP, AH,AHP,SgP,SdgP,SktP,SkP,
     .        SgU,SdgU,SktU,SkU,itmpl,
     .        theta_indx,phi_indx,theta_red,phi_red, cindx,
     .       g11, g12, g13, g22, g23, g33, k11, k12, k13, k22, k23, k33,
     .        GrX,GrY,GrZ,msk,func,tfunc,lx,ly,lz,order,
     .        hx,hy,hz,NX,NY,NZ,
     .        NEQN,NRED,NST,NSP,success)
            implicit none
            
            include 'horizon_parameters.h'

            integer    NST, NSP, NX, NY, NZ, order

            integer    itmpl(NST,NSP,4,4,2),
     .            NRED, NEQN,
     .            theta_indx(NSTM*NSPM),phi_indx(NSTM*NSPM),
     .            theta_red(2*NSPM), phi_red(2*NSPM),
     .            cindx(NSTM,NSPM)
            
            real*8    jac(NST*NSP,NST*NSP)

        real*8  S(NST,NSP,3), Sc(NST,NSP,3), Sf(NST,NSP),
     .          SP(NST,NSP,3), ScP(NST,NSP,3), SfP(NST,NSP),
     .          Sfmol(NST,NSP,3,3,3), 
     .          Sn(NST,NSP,3), Sb(NST,NSP,3,3),
     .          SnP(NST,NSP,3), SbP(NST,NSP,3,3),
     .          Sg(NST,NSP,3,3), Sdg(NST,NSP,3,3,3),
     .          Skt(NST,NSP,3,3), Sk(NST,NSP),
     .          SgP(NST,NSP,3,3), SdgP(NST,NSP,3,3,3),
     .          SktP(NST,NSP,3,3), SkP(NST,NSP),
     .          SgU(NST,NSP,3,3), SdgU(NST,NSP,3,3,3),
     .          SktU(NST,NSP,3,3), SkU(NST,NSP),
     .          AH(NST,NSP),AHP(NST,NSP)
        logical success

c     Cartesian Grid Functions

       real*8  GrX(NX),GrY(NY),GrZ(NZ),msk(NX,NY,NZ),
     .      g11(NX,NY,NZ), g12(NX,NY,NZ),g13(NX,NY,NZ),
     .      g22(NX,NY,NZ), g23(NX,NY,NZ),g33(NX,NY,NZ),
     .      k11(NX,NY,NZ), k12(NX,NY,NZ),k13(NX,NY,NZ),
     .      k22(NX,NY,NZ), k23(NX,NY,NZ),k33(NX,NY,NZ)

c     Temporaray arrays

       real*8   tfunc(NST*NSP),func(NST,NSP),
     .   lx(NST*NSP), ly(NST*NSP),lz(NST*NSP)

       real*8  hx, hy, hz



c        Local Variables

       integer  i,j,k,m,p,q,r,t,mu,nu,nt,np,
     *          th0,ph0,thF,phF

       real*8   sgn, divh, divh2, divh4, diveps, tmp


        divh = 0.5e0 / h

        divh2 = 1.0e0 / h**2

c        write(6,*)'EPSILON = ',epsilon
        diveps = 1.0e0 / epsilon

c        Zero out the jacobian matrix to start with
        do nu=1,NEQN
          do mu=1,NEQN
            jac(mu,nu) = 0.0e0
          enddo
        enddo

c        Perturb everyone

        do i = 1, NST
          do j = 1, NSP

            SfP(i,j) = Sf(i,j) + epsilon
            SP(i,j,1) = S(i,j,1) + epsilon
            SP(i,j,2) = S(i,j,2)
            SP(i,j,3) = S(i,j,3)
            ScP(i,j,1) = SfP(i,j)*cos(S(i,j,3))*sin(S(i,j,2))
            ScP(i,j,2) = SfP(i,j)*sin(S(i,j,3))*sin(S(i,j,2))
            ScP(i,j,3) = SfP(i,j)*cos(S(i,j,2))
            ScP(i,j,1) = ScP(i,j,1) + center(1)
            ScP(i,j,2) = ScP(i,j,2) + center(2)
            ScP(i,j,3) = ScP(i,j,3) + center(3)

          end do
        end do

c        Calculate perturbed and unperturbed metric


c      call MetData(SP,ScP,SgP,SdgP,SktP,SkP,1,1,NST,NSP,NST,NSP)
c         call MetIntp(SP,ScP,SgP,SdgP,SkP,SktP,GrX,GrY,GrZ,GrG,GrDG,
c     *            GrKt,hx,hy,hz,1,NST,1,NSP,
c     *            NX,NY,NZ,NST,NSP)
c      call metextp(grx,gry,grz,g11,g12,g13,g22,g23,g33,
c     .            k11,k12,k13,k22,k23,k33, SgU, SdgU, SktU,
c     .            Sc, lx, ly, lz, msk,
c     .            wkspc, func, tfunc,order, NST,NSP,NX,NY,NZ)

c      call metextp(grx,gry,grz,g11,g12,g13,g22,g23,g33,
c     .            k11,k12,k13,k22,k23,k33, SgP, SdgP, SktP,
c     .            ScP, lx, ly, lz, msk,
c     .            wkspc, func, tfunc,order, NST,NSP,NX,NY,NZ)
        call metextp3(grx,gry,grz,g11,g12,g13,g22,g23,g33,
     .            k11,k12,k13,k22,k23,k33, SgU, SdgU, SktU,
     .            Sc, lx, ly, lz, msk,
     .            func, tfunc,order, NST,NSP,NX,NY,NZ,success)
        if (.not. success) then
          return
        end if
        call metextp3(grx,gry,grz,g11,g12,g13,g22,g23,g33,
     .            k11,k12,k13,k22,k23,k33, SgP, SdgP, SktP,
     .            ScP, lx, ly, lz, msk,
     .            func, tfunc,order, NST,NSP,NX,NY,NZ,success)
        if (.not. success) then
          return
        end if

c      call MetData(S,Sc,SgU,SdgU,SktU,SkU,1,1,NST,NSP,NST,NSP)
c         call MetIntp(S,Sc,SgU,SdgU,SkU,SktU,GrX,GrY,GrZ,GrG,GrDG,
c     *            GrKt,hx,hy,hz,1,NST,1,NSP,
c     *            NX,NY,NZ,NST,NSP)

c        Copy over unperturbed functions

        do i = 1, NST
          do j = 1, NSP
            SfP(i,j) = Sf(i,j)
            SP(i,j,1) = S(i,j,1)
            SP(i,j,2) = S(i,j,2)
            SP(i,j,3) = S(i,j,3)

            ScP(i,j,1) = Sc(i,j,1)
            ScP(i,j,2) = Sc(i,j,2)
            ScP(i,j,3) = Sc(i,j,3)

            Sk(i,j) = SkU(i,j)

          end do
        end do


        do m=1,3
          do k=1,3
            do j=1,NSP
              do i = 1,NST
                Sg(i,j,k,m) = SgU(i,j,k,m)
                Skt(i,j,k,m) = SktU(i,j,k,m)
              end do
            end do
          end do
        end do

        do p=1,3
          do m=1,3
            do k=1,3
              do j=1,NSP
                do i=1,NST
                  Sdg(i,j,k,m,p) = SdgU(i,j,k,m,p)
                end do
              end do
            end do
          end do
        end do

c    Interpolate and create finite difference molecule.
        call IntpMol(S,Sc,Sf,Sfmol,itmpl,
     *      1,1,NST,NSP,NST,NSP)

c    Calculate derivatives
        call Derivatives(S,Sfmol,Sn,Sb,divh,divh2,1,1,
     $               NST,NSP,NST,NSP)

c    Evaluate the apparent horizon equation.
        call Evaluator(S,Sg,Sdg,Skt,Sk,AH,Sn,Sb,Sfmol,
     #                  1,1,NST,NSP,NST,NSP)

c        write(0,*)'Looping'
c    Loop over all the points

        do mu = 1, NEQN
          i = theta_indx(mu)
          j = phi_indx(mu)
c Perturb point(s)
c Ensure that the pole points are set up appropriately
          if(i .eq. 1 .or. i .eq. NST)then
c        If at the pole perturb along phi
            do k=1,NSP
              SfP(i,k) = SfP(i,k) + epsilon
              SP(i,k,1) = SfP(i,k)

              ScP(i,k,1)= SP(i,k,1)*cos(S(i,k,3))*sin(S(i,k,2))
              ScP(i,k,2)= SP(i,k,1)*sin(S(i,k,3))*sin(S(i,k,2))
              ScP(i,k,3)= SP(i,k,1)*cos(S(i,k,2))
              ScP(i,k,1)= ScP(i,k,1) + center(1)
              ScP(i,k,2)= ScP(i,k,2) + center(2)
              ScP(i,k,3)= ScP(i,k,3) + center(3)

            enddo
          else
            SfP(i,j) =  SfP(i,j) + epsilon
            SP(i,j,1) = SfP(i,j)
            ScP(i,j,1)= SP(i,j,1)*cos(S(i,j,3))*sin(S(i,j,2))
            ScP(i,j,2)= SP(i,j,1)*sin(S(i,j,3))*sin(S(i,j,2))
            ScP(i,j,3)= SP(i,j,1)*cos(S(i,j,2))
            ScP(i,j,1)= ScP(i,j,1) + center(1)
            ScP(i,j,2)= ScP(i,j,2) + center(2)
            ScP(i,j,3)= ScP(i,j,3) + center(3)
          endif
c     Set up evaluation points (ie: Set up neighborhood over which to
c                                            evaluate.)
          if(i .le. 5 .or. i .gt. NST-5)then
            if(i .le. 5)then
              th0 = 1
              thF = 7
            else
              th0 = NST - 7
              thF = NST
            end if
c               th0 = 1
c               thF = NST
            ph0 = 1
            phF = NSP
          else
            th0 = i - 4
            thF = i + 4
            ph0 = j - 3
            phF = j + 3
          end if

c    write(6,*)'=>',i,j
c    write(6,*)'==>',th0,ph0,thF,phF
c        Metric functions at perturbed point

          if( i .ne. 1 .and. i .ne. NST)then
            do m=1,3
              do k=1,3
                Sg(i,j,k,m) = SgP(i,j,k,m)
                Skt(i,j,k,m) = SktP(i,j,k,m)
                Sdg(i,j,k,m,1) = SdgP(i,j,k,m,1)
                Sdg(i,j,k,m,2) = SdgP(i,j,k,m,2)
                Sdg(i,j,k,m,3) = SdgP(i,j,k,m,3)
              end do
            end do
            Sk(i,j) = SkP(i,j)
c              call MetData(SP,ScP,Sg,Sdg,Skt,Sk,i,j,i,j,NST,NSP)
          else
            do p=1,NSP
              do m=1,3
                do k=1,3
                  Sg(i,p,k,m) = SgP(i,p,k,m)
                  Skt(i,p,k,m) = SktP(i,p,k,m)
                  Sdg(i,p,k,m,1) = SdgP(i,p,k,m,1)
                  Sdg(i,p,k,m,2) = SdgP(i,p,k,m,2)
                  Sdg(i,p,k,m,3) = SdgP(i,p,k,m,3)
                end do
              end do
              Sk(i,p) = SkP(i,p)
            end do
c           call MetData(SP,ScP,Sg,Sdg,Skt,Sk,i,1,i,NSP,NST,NSP)
          endif

          if(ph0 .lt. 1 .or. phF .gt. NST)then
            ph0 = 1
            phF = 5

            call IntpMol(SP,ScP,SfP,Sfmol,itmpl,
     *                th0,ph0,thF,phF,NST,NSP)

            call Derivatives(SP,Sfmol,SnP,SbP,divh,divh2,
     $               th0,ph0,thF,phF,NST,NSP)

            call Evaluator(SP,Sg,Sdg,Skt,Sk,AHP,SnP,SbP,Sfmol,
     #                          th0,ph0,thF,phF,NST,NSP)

c        Calculate jacobian matrix component

            do k =th0,thF
              do m=ph0,phF
                nu = cindx(k,m)
                jac(nu,mu) = diveps*(AHP(k,m) - AH(k,m))
              end do
            end do

            do k = th0, thF
              do m = ph0, phF
                AHP(k,m) = AH(k,m)
                SnP(k,m,1) = Sn(k,m,1)
                SnP(k,m,2) = Sn(k,m,2)
                SnP(k,m,3) = Sn(k,m,3)

                SbP(k,m,1,1) = Sb(k,m,1,1)
                SbP(k,m,2,2) = Sb(k,m,2,2)
                SbP(k,m,3,3) = Sb(k,m,3,3)
                SbP(k,m,1,2) = Sb(k,m,1,2)
                SbP(k,m,1,3) = Sb(k,m,1,3)
                SbP(k,m,2,3) = Sb(k,m,2,3)
                SbP(k,m,2,1) = Sb(k,m,2,1)
                SbP(k,m,3,1) = Sb(k,m,3,1)
                SbP(k,m,3,2) = Sb(k,m,3,2)
              end do
            end do

            ph0 = NSP - 5
            phF = NSP

            call IntpMol(SP,ScP,SfP,Sfmol,itmpl,
     *          th0,ph0,thF,phF,NST,NSP)

            call Derivatives(SP,Sfmol,SnP,SbP,divh,divh2,
     $               th0,ph0,thF,phF,NST,NSP)


            call Evaluator(SP,Sg,Sdg,Skt,Sk,AHP,SnP,SbP,Sfmol,
     #               th0,ph0,thF,phF,NST,NSP)

            do k =th0,thF
              do m=ph0,phF
                nu = cindx(k,m)
                jac(nu,mu) = diveps*(AHP(k,m) - AH(k,m))
              end do
            end do

            do k = th0, thF
              do m = ph0, phF
                AHP(k,m) = AH(k,m)
                SnP(k,m,1) = Sn(k,m,1)
                SnP(k,m,2) = Sn(k,m,2)
                SnP(k,m,3) = Sn(k,m,3)

                SbP(k,m,1,1) = Sb(k,m,1,1)
                SbP(k,m,2,2) = Sb(k,m,2,2)
                SbP(k,m,3,3) = Sb(k,m,3,3)
                SbP(k,m,1,2) = Sb(k,m,1,2)
                SbP(k,m,1,3) = Sb(k,m,1,3)
                SbP(k,m,2,3) = Sb(k,m,2,3)
                SbP(k,m,2,1) = Sb(k,m,2,1)
                SbP(k,m,3,1) = Sb(k,m,3,1)
                SbP(k,m,3,2) = Sb(k,m,3,2)
              end do
            end do
          else
            call IntpMol(SP,ScP,SfP,Sfmol,itmpl,
     *                  th0,ph0,thF,phF,NST,NSP)

            call Derivatives(SP,Sfmol,SnP,SbP,divh,divh2,
     $               th0,ph0,thF,phF,NST,NSP)

            call Evaluator(SP,Sg,Sdg,Skt,Sk,AHP,SnP,SbP,Sfmol,
     #                  th0,ph0,thF,phF,NST,NSP)

c        Calculate jacobian matrix component

            do k =th0,thF
              do m=ph0,phF
                nu = cindx(k,m)
                jac(nu,mu) = diveps*(AHP(k,m) - AH(k,m))
              end do
            end do

            do k = th0, thF
              do m = ph0, phF
                AHP(k,m) = AH(k,m)
                SnP(k,m,1) = Sn(k,m,1)
                SnP(k,m,2) = Sn(k,m,2)
                SnP(k,m,3) = Sn(k,m,3)

                SbP(k,m,1,1) = Sb(k,m,1,1)
                SbP(k,m,2,2) = Sb(k,m,2,2)
                SbP(k,m,3,3) = Sb(k,m,3,3)
                SbP(k,m,1,2) = Sb(k,m,1,2)
                SbP(k,m,1,3) = Sb(k,m,1,3)
                SbP(k,m,2,3) = Sb(k,m,2,3)
                SbP(k,m,2,1) = Sb(k,m,2,1)
                SbP(k,m,3,1) = Sb(k,m,3,1)
                SbP(k,m,3,2) = Sb(k,m,3,2)
              end do
            end do
          end if

          if(i .eq. 1 .or. i .eq. NST)then
            do k=1,NSP
              SfP(i,k) = Sf(i,k) 
              SP(i,k,1) = Sf(i,k)

              ScP(i,k,1)= Sc(i,k,1)
              ScP(i,k,2)= Sc(i,k,2)
              ScP(i,k,3)= Sc(i,k,3)
            enddo
          else
            SfP(i,j) =  Sf(i,j)
            SP(i,j,1) = Sf(i,j)
            ScP(i,j,1)= Sc(i,j,1)
            ScP(i,j,2)= Sc(i,j,2)
            ScP(i,j,3)= Sc(i,j,3)
          endif

          if( i .ne. 1 .and. i .ne. NST)then
            do m=1,3
              do k=1,3
                Sg(i,j,k,m) = SgU(i,j,k,m)
                Skt(i,j,k,m) = SktU(i,j,k,m)
                Sdg(i,j,k,m,1) = SdgU(i,j,k,m,1)
                Sdg(i,j,k,m,2) = SdgU(i,j,k,m,2)
                Sdg(i,j,k,m,3) = SdgU(i,j,k,m,3)
              end do
            end do
            Sk(i,j) = SkU(i,j)
          else
            do p=1,NSP
              do m=1,3
                do k=1,3
                  Sg(i,p,k,m) = SgU(i,p,k,m)
                  Skt(i,p,k,m) = SktU(i,p,k,m)
                  Sdg(i,p,k,m,1) = SdgU(i,p,k,m,1)
                  Sdg(i,p,k,m,2) = SdgU(i,p,k,m,2)
                  Sdg(i,p,k,m,3) = SdgU(i,p,k,m,3)
                end do
              end do
              Sk(i,p) = SkU(i,p)
            end do
          endif

        end do

        success = .true.
        return
        end
