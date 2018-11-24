c=======================================================================
c  solver.f
c
c    Created Fri Jul  1 13:39:59 CDT 1994 Mijan Huq
c
c  This file contains routines for calculating the jacobian and
c  using it to solve the apparent horizon equation.
c
c    Modified Tue Jul 26 14:05:42 CDT 1994 MFH
c        .    Removed the analytic jacobian calculation routines and arrays.
c            Converted to numerical perturbation approach.
c=======================================================================
      subroutine SlvAH(S,Sc,Sf,Sfmol,Sn,Sb,AH,AHP,
     *      Sg,Sdg,Skt,Sk,itmpl,
     .        SP,ScP,SfP,SnP,SbP,SgP,SdgP,SktP,SkP,SgU,SdgU,SktU,SkU,
     .      g11, g12, g13, g22, g23, g33, k11, k12, k13, k22, k23, k33,
     .      GrX,GrY,GrZ,msk,func,tfunc,lx,ly,lz,
     .      order,hx,hy,hz,NX,NY,NZ,
     *      NST,NSP,success)
         implicit none

         include 'horizon_parameters.h'

         integer  NST, NSP, NX, NY, NZ,order,
     *            itmpl(NST,NSP,4,4,2)

c  Unperturbed Surface Grid Functions
         real*8   S(NST,NSP,3),  Sc(NST,NSP,3), Sf(NST,NSP),
     *         Sfmol(NST,NSP,3,3,3),
     *         Sn(NST,NSP,3), Sb(NST,NSP,3,3),AH(NST,NSP)

c    Perturbed Surface Grid Functions
            real*8    SP(NST,NSP,3),ScP(NST,NSP,3), SfP(NST,NSP),
     *            SnP(NST,NSP,3),SbP(NST,NSP,3,3)

         logical success

c  Surface Metric Grid Functions
      real*8   Sg(NST,NSP,3,3), Sdg(NST,NSP,3,3,3),
     *         Skt(NST,NSP,3,3), Sk(NST,NSP), AHP(NST,NSP)
      real*8   SgP(NST,NSP,3,3), SdgP(NST,NSP,3,3,3),
     *         SktP(NST,NSP,3,3), SkP(NST,NSP)
      real*8   SgU(NST,NSP,3,3), SdgU(NST,NSP,3,3,3),
     *         SktU(NST,NSP,3,3), SkU(NST,NSP)

c     Cartesian Grid Functions

      real*8  GrX(NX),GrY(NY),GrZ(NZ),msk(NX,NY,NZ),
     .      g11(NX,NY,NZ), g12(NX,NY,NZ),g13(NX,NY,NZ),
     .      g22(NX,NY,NZ), g23(NX,NY,NZ),g33(NX,NY,NZ),
     .      k11(NXM,NYM,NZM), k12(NXM,NYM,NZM),k13(NXM,NYM,NZM),
     .      k22(NXM,NYM,NZM), k23(NXM,NYM,NZM),k33(NXM,NYM,NZM)

c     Temporaray arrays

      real*8  tfunc(NST*NSP),func(NST,NSP),
     .   lx(NST*NSP), ly(NST*NSP),lz(NST*NSP)

      real*8  hx, hy, hz, sss
      integer  ti2,tf2



c  Local Variables- - - - - - - - - - - - - - - - - - - -

      integer  i,j,k, ii,jj, iter,NEQN, NRED, NNZ

c  Indexing for conversion between two type of indexing used here.

      integer theta_indx(NSTM*NSPM),phi_indx(NSTM*NSPM),
     *    cindx(NSTM,NSPM),theta_red(2*NSTM), phi_red(2*NSTM)

c  Jacobian and related

      integer  ia(NSTM*NSPM),ja(NSTM*NSPM*NSTM),
     *            itemp(NSTM*NSPM*NSTM*NSPM)

c  Modification -MFH
c     removed jac array and used rtemp as temp array
      real*8     bb(NSTM*NSPM),
     *         hfunc(NSTM*NSPM), jac2(NSTM*NSPM*NSTM*NSPM),
     *         L2AH,L2H,rcond,
     *            divh, divh2, dcntr(3),
     *            rtemp(NSTM*NSPM*NSTM*NSPM)

        real*8    Intg(NSTM,NSPM), area, pl2h, pl2ah

c  Functions Used

      real*8   DNRM2
      integer  time
      external time

        character*64 stem

        logical stat, silucg

        integer    ie, it, nn, ifail
        real*8        eps,jnorm, ti, tf

        divh = 0.5e0 / h
        divh2= 1.0e0 /h**2

c    Set up the index arrays

      call SetUpIndex(theta_indx, phi_indx,theta_red,phi_red,
     #          cindx,NEQN,NRED,NST,NSP)

c        open(unit=99,file='band.dat')

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c    Start carrying out Newton Iterations.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        nn = NEQN
        ifail = 0

c        call second(sss)
c                sss =time()
c    write(0,*)'time at start of newton iterations ',sss
        do iter = 1, 20

c            write(stem,1043)iter
c            call WriteS2F(stem,S,AH,NST,NSP)


            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c        Call the jacobian generation routine. (See jacobian.f)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c            call second(ti)
c             ti = time()
          call PAnlJac(rtemp,S,Sc,Sf,Sfmol,Sn,Sb,Sg,Sdg,Skt,Sk,
     .      SP,ScP,SfP,SnP,SbP, AH,AHP,SgP,SdgP,SktP,SkP,
     .      SgU,SdgU,SktU,SkU,itmpl,
     .      theta_indx,phi_indx,theta_red,phi_red, cindx,
     .      g11, g12, g13, g22, g23, g33, k11, k12, k13, k22, k23, k33,
     .      GrX,GrY,GrZ,msk,func,tfunc,lx,ly,lz,order,
     .      hx,hy,hz,NX,NY,NZ,
     .      NEQN,NRED,NST,NSP,success)
          if (.not. success) then
            return
          end if
c            call second(tf)
c             tf = time()
c            write(0,*)'Time taken to generate jacobian ',tf-ti

c            stop



c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c        Copy over the jacobian from jac to jac2 using the format needed
c        by ilucg solver. This part of the code needs proper modification
c        once a final solver is chosen. Need to optimize for memory here.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

          call tosparse(rtemp,jac2,ia,ja,jnorm,NEQN,NST,NSP)


c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c        Copy over the RHS from AH to bb which is a NEQN sized array
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        do ii=1,NEQN
c               i = theta_indx(ii)
c               j = phi_indx(ii)
        bb(ii) = -AH(theta_indx(ii),phi_indx(ii))
        end do

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c        If this is the first iteration put in a value for the starting
c        hfunc that the ILUCG solver can use.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        if(iter .eq. 1)then
          pl2h = 99999.0
      pl2ah= 999999.0
      do ii=1,NEQN
c       i = theta_indx(ii)
c       j = phi_indx(ii)
        hfunc(ii) = 0.001e0
        bb(ii) = -AH(theta_indx(ii),phi_indx(ii))
          end do
      end if

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c           Calculate the L2 norm for the RHS
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            L2AH = DNRM2(NEQN,bb,1)/sqrt(float(NEQN))
            write(6,*)'L_2 norm of AH before ILUCG call',L2AH

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c                Call the ILUCG solver.
c                See ilucg/ilucg.f for comments on what the values of 
c                eps and it mean.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                eps = ilucgstp
                it = 0
c                call second(sss)
c                                ti2 = time()
        stat = silucg(NEQN,ia,ja,jac2,bb,hfunc,itemp,rtemp,eps,it,ie)  
c                call second(sss)
c                                tf2 = time()
c                                write(6,*)'time taken for silucg',tf2-ti2
c                                write(2,*)'time taken for silucg',tf2-ti2
c                write(0,*)'time after ',sss

c                write (6,*) 'return status = ',stat, '          ie =',ie

c            write(6,*)'L_2 norm of AH after ILUCG call',L2AH

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c           Calculate the L2-norm for the solution
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            L2H = DNRM2(NEQN,hfunc,1)/sqrt(float(NEQN))

                write(0,*)'Iteration Number ',iter
            write(0,*)'Norm of RHS before solve ',L2AH
            write(0,*)'Norm of solution ',L2H

c            write(2,*)'rcond =',rcond
c                write(2,*)'Iteration Number ',iter
c                write (2,*) 'return status = ',stat, '          ie =',ie
c            write(2,*)'Norm of RHS before solve ',L2AH
c            write(2,*)'Norm of solution ',L2H
c

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c        If the L2-norm of the apparent horizon equation is a factor of 10
c        greater than the last iteration then stop
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
c            if(L2H .gt. 5.0 * pl2h .or. L2AH .gt. 2.0)then
c                ifail = ifail + 1
c                write(0,*)'Changing epsilon by factor of ten'
c                write(2,*)'Changing epsilon by factor of ten'
c                epsilon = epsilon * 0.1
c                write(0,*)'New value of epsilon = ',epsilon
c                if(ifail .gt. 4)then
c                    write(0,*)'Exceeded max number of failures'
c                    write(0,*)'current value of ifail = ',ifail
c                    write(0,*)'Exiting: NO CONVERGENCE'
c                    write(2,*)'Exceeded max number of failures'
c                    write(2,*)'current value of ifail = ',ifail
c                    write(2,*)'Exiting: NO CONVERGENCE'
c                    goto 1040
c                end if
c            else
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                Update Sf by hfunc
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c                write(stem,1042)iter
c                open(unit=76,file=stem)

c                do i=1,NEQN
c                    write(76,*)i,hfunc(i)
c                end do
c                close(unit = 76)


                ii=1
                        do j=1,NSP
                 Sf(theta_indx(ii),j) = Sf(theta_indx(ii),j) + hfunc(ii)

                      S(theta_indx(ii),j,1) = Sf(theta_indx(ii),j)

                      Sc(theta_indx(ii),j,1) = S(theta_indx(ii),j,1)*
     #             cos(S(theta_indx(ii),j,3))*sin(S(theta_indx(ii),j,2))
     #                  + center(1)
                      Sc(theta_indx(ii),j,2)= S(theta_indx(ii),j,1)*
     #             sin(S(theta_indx(ii),j,3))*sin(S(theta_indx(ii),j,2))
     #                  + center(2)
                      Sc(theta_indx(ii),j,3)= S(theta_indx(ii),j,1)*
     #                    cos(S(theta_indx(ii),j,2))
     #                        + center(3)

                        end do
                ii=NEQN
                        do j=1,NSP
                 Sf(theta_indx(ii),j) = Sf(theta_indx(ii),j) + hfunc(ii)

                     S(theta_indx(ii),j,1) = Sf(theta_indx(ii),j)

                      Sc(theta_indx(ii),j,1) = S(theta_indx(ii),j,1)*
     #             cos(S(theta_indx(ii),j,3))*sin(S(theta_indx(ii),j,2))
     #                  + center(1)
                      Sc(theta_indx(ii),j,2)= S(theta_indx(ii),j,1)*
     #             sin(S(theta_indx(ii),j,3))*sin(S(theta_indx(ii),j,2))
     #                  + center(2)
                      Sc(theta_indx(ii),j,3)= S(theta_indx(ii),j,1)*
     #                    cos(S(theta_indx(ii),j,2))
     #                        + center(3)

                        end do

                do ii=2,NEQN-1
                   Sf(theta_indx(ii),phi_indx(ii)) = Sf(theta_indx(ii),
     #              phi_indx(ii)) + hfunc(ii)

                  S(theta_indx(ii),phi_indx(ii),1) = Sf(theta_indx(ii),
     #                phi_indx(ii))

                  Sc(theta_indx(ii),phi_indx(ii),1) = 
     #            S(theta_indx(ii),phi_indx(ii),1)*cos(S(theta_indx(ii),
     #            phi_indx(ii),3))*sin(S(theta_indx(ii),phi_indx(ii),2))
     #                       + center(1)

                  Sc(theta_indx(ii),phi_indx(ii),2)= 
     #            S(theta_indx(ii),phi_indx(ii),1)*
     #            sin(S(theta_indx(ii),phi_indx(ii),3))*
     #            sin(S(theta_indx(ii),phi_indx(ii),2))
     #                       + center(2)

                  Sc(theta_indx(ii),phi_indx(ii),3)= 
     #            S(theta_indx(ii),phi_indx(ii),1)*
     #            cos(S(theta_indx(ii),phi_indx(ii),2))
     #                    +    center(3)

            end do


c                write(stem,1041)iter
c
c                call WriteS2F(stem,S,Sf,NST,NSP)
c            open(unit=59,file=stem)
c            call Dp3(Sc,NST,NSP,59)
c            close(unit=59)
c            write(stem,1044)iter
c            open(unit=59,file=stem)
c            call Dp3(S,NST,NSP,59)
c            close(unit=59)
c            write(stem,1045)iter
c            open(unit=59,file=stem)
c            call Dp2(AH,NST,NSP,59)
c            close(unit=59)

c  Set up the integrand

      call Horizon_Integrand(S,Sc,Sg,Intg,NST,NSP)

c  Integrate using the trapezoidal rule

      call SurfArea(Intg,area,dtheta,dphi,NST,NSP)

      write(0,*)'Area of Apparent horizon = ',area
c      write(2,*)'Area of Apparent horizon = ',area
c            Store current as old
            pl2ah = l2ah
            pl2h = l2h
c            end if

c        call findcntr(Sc,S,center,NST,NSP)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c            If L2-norm of apparent horizon equation is less than the
c            tolerable amount then stop newton iteration.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c              if(iter .le. 2 .and. pl2h .ge. 1.0e-2)then
c                 write(0,*)'Repositioning '
c                 call LocateCenter2(S,Sc,center,NST,NSP)    
c                 call  findrad(radius,S,NST,NSP)
c             call GenSur(S,Sc,SF,3,NST,NSP)
c              end if
           if(L2aH .lt. newtst)then
               goto 1040
            endif

         end do
1040     continue
c        call second(sss)
c                sss = time()
c        write(0,*)'time at end of newton iterations ',sss

                    write(0,*)'l2-norm of jacobian =',jnorm
c                    write(2,*)'l2-norm of jacobian =',jnorm
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c            Dump solution to files.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            write(0,*)'Ending value of epsilon = ',epsilon
c            write(2,*)'Ending value of epsilon = ',epsilon
            open(unit=59,file='solncrt.dat')
c            call Dump3CoordCart(Sc,NST,NSP,59)
            call Dp3(Sc,NST,NSP,59)
            close(unit=59)
            open(unit=59,file='solnsph.dat')
c            call Dump3CoordCart(S,NST,NSP,59)
            call Dp3(S,NST,NSP,59)
            close(unit=59)
c            close(unit=99)

1041        format('Step.',I2)
1042        format('delh.',I1)
1043        format('eval.',I1)
1044        format('Stepsph.',I2)
1045        format('Residual.',I2)

         success = .true.
         return
      end
c=======================================================================
c        Routines to dump data to files.
c            Dumps rank 3 arrays
c
                subroutine Dp3(Sc,NST,NSP,uu)

                        implicit none

                                include 'horizon_parameters.h'
                        integer NST, NSP, uu

                        real*8    Sc(NST,NSP,3)

                        integer i,j

                        do i=1,NST
                                do j=1,NSP

                          write(uu,110)Sc(i,j,1),Sc(i,j,2),Sc(i,j,3)

                                end do
                                write(uu,*)
                        end do

110            format(3F20.15)

                        return
                end

c        Dumps rank 2 arrays            

                subroutine Dp2(Sf,NST,NSP,uu)

                        implicit none
                                include 'horizon_parameters.h'

                        integer NST, NSP, uu

                        real*8    Sf(NST,NSP)

                        integer i,j

                        do i=1,NST
                                do j=1,NSP

                          write(uu,111)i,j,Sf(i,j)

                                end do
                                write(uu,*)
                        end do
                close(unit = 99)

111            format(2I4,F20.15)
                        return
                end
c-----------------------------------------------------------------------
c    Routine to
c        Copy over the jacobian from jac to jac2 using the format needed
c        by ilucg solver. This part of the code needs proper modification
c        once a final solver is chosen. Need to optimize for memory here.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        subroutine tosparse(jac,jac2,ia,ja,jacnorm,NEQN,NST,NSP)
            implicit none
            include 'horizon_parameters.h'
    
            integer    NST, NSP, NEQN

            integer  ia(NSTM*NSPM),ja(NSTM*NSPM*NSTM*NSPM)

            real*8 jac(NST*NSP,NST*NSP), jac2(NSTM*NSPM*NSTM*NSPM),
     .        jacnorm

c        Local variables 

            integer ii,jj,i,j, count

            jacnorm = 0.0
            do j =1, NEQN
                do i =1, NEQN
                    jacnorm = jacnorm + jac(i,j)**2
                end do
            end do
            jacnorm = sqrt(jacnorm)/float(NEQN)

                open(unit=87,file='count.dat')
            count =0
            ii = 1
            jj = 1
            do i=1,NEQN
                ia(jj) = ii
                jj = jj + 1
                do j=1,NEQN
c                    Check if the term is non-zero
                    if(jac(i,j) .ne. 0.0e0)then
                        count = count + 1
                        ja(ii) = j
                        jac2(ii) = jac(i,j)
                        ii = ii + 1
                    end if
                end do
                    write(87,'(2I6)')i, count
                    count = 0
            end do
            ia(jj) = ii
            close(unit=87)

            return
        end
c-----------------------------------------------------------------------
c
        subroutine findcntr(Sc,S,center,NST,NSP)
            implicit none

            integer    NST, NSP

            real*8    Sc(NST,NSP,3), S(NST,NSP,3),center(3)

c        Local variables

            integer    i, j, k

            real*8    sumx,sumy,sumz

            sumx = 0.0
            sumy = 0.0
            sumz = 0.0

            do i = 1, NST
                do j = 1, NSP
                    sumx = sumx + Sc(i,j,1)
                    sumy = sumy + Sc(i,j,2)
                    sumz = sumz + Sc(i,j,3)
                end do
            end do

            center(1) = sumx/(NST*NSP)
            center(2) = sumy/(NST*NSP)
            center(3) = sumz/(NST*NSP)

            write(0,*)'New center',center(1),center(2),center(3)

            return
        end
