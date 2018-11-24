c-----------------------------------------------------------------------
c   apparent.f
c   
c   Written by 
c      Mijan Huq 
c      Center for Relativity,
c      University of Texas at Austin.
c
c   This file contains routines that will evaluate the apparent horizon
c   equation given a surface in S, metric data and extrinsic curvature
c   data.
c
c   S(NST,NSP,3)        -   Surface coordintes (Spherical)
c   Sc(NST,NSP,3)       -   Surface coordinates (Cartesian)
c   Sf(NST,NSP)           -    Surface function. Same as S(,,1)
c   Sg(NST,NSP,3,3)     -    Metric tensor
c   Sdg(NST,NSP,3,3,3)  -   Derivatives of the above
c  Skt(NST,NSP,3,3)    -   Extrinsic Curvature tensor
c   Sk(NST,NSP)         -   Extrinsic Curvature scalar
c   Sn(NST,NSP,3)       -   Normal vectors
c   AH(NST,NSP)         -   Apparent horizon equation evaluation
c
c   Need to pass a work array
c   Sfmol(NST,NSP,3,3,3)-   Surface Molecule function values.
c   Sb(NST,NSP,3,3)     -   Second and cross derivatives
c
c
c   Evaluator routine 
c   Additions to that are:
c
c   BNT   :=   Begining theta
c   BNP   := Begining phi
c   NSST  := End theta
c   NSSP  := End phi
c

      subroutine EvalAH(S,Sc,Sf,Sg,Sdg,Skt,Sk,AH,Sn,Sb,Sfmol,Simol,
     #   itmpl,BNT,BNP,NSST,NSSP,NST,NSP)
         implicit none
         include 'horizon_parameters.h'

         integer   NST,NSP
         integer   BNT,BNP
         integer   NSST,NSSP

               
         real*8   S(NST,NSP,3)
         real*8   Sc(NST,NSP,3)
         real*8   Sf(NST,NSP)
         real*8   Sg(NST,NSP,3,3)
         real*8   Sdg(NST,NSP,3,3,3)
         real*8   Skt(NST,NSP,3,3)
         real*8   Sk(NST,NSP)

         real*8   AH(NST,NSP)
         real*8   Sn(NST,NSP,3)
         real*8   Sb(NST,NSP,3,3)
         real*8   Sfmol(NST,NSP,3,3,3)
         real*8   Simol(NST,NSP,3,3,3,3)
         integer   itmpl(NST,NSP,4,4,2)


c   Local Variables

         integer   i,j,k,m,p,q,r,
     *            iref(4,4,2)

         real*8   xx(3),tx(3),ax(3)

         real*8   tmpl(3,3,3,3)

         real*8   divh, divh2
         real*8   temp

         real*8   at(4),ap(4)





c   Calculate half inverse h and inverse h square

         divh = 0.5e0/h
         divh2= 1.0e0/h**2

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Calculate derivatives of Sf and store in Sn and Sb
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         call Derivatives(S,Sfmol,Sn,Sb,divh,divh2,BNT,BNP,
     $               NSST,NSSP,NST,NSP)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Evaluate the apparent horizon and store in AH
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


         call Evaluator(S,Sg,Sdg,Skt,Sk,AH,Sn,Sb,Sfmol,
     #                  BNT,BNP,NSST,NSSP,NST,NSP)


         return
      end
c

c=======================================================================
c   Routine to calculate the derivatives via centered differencing.
c
c   Input:
c      S       := Spherical coordinate locations of surface mesh points.
c      Sfmol   := Finite difference molecule function values.
c      divh  := 1/(2h) where h is the finite difference molecule size
c      divh2 := 1/h^2  where h is the same as above...
c      BNT   := Begining theta mesh coordinate.
c      BNP   := Begining phi mesh coordinate.
c      NSST   := Ending theta mesh coordinate.
c      NSSP   := Ending phi mesh coordinate.
c
c   Output:
c      Sn      := First derivatives.
c      Sb      := Second derivatives.
c
      subroutine QDerivatives(S,Sfmol,Sn,Sb,divh,divh2,BNT,BNP,
     $                  NSST,NSSP,NST,NSP)

         implicit none
         include 'horizon_parameters.h'
         
         integer    NST,NSP
         integer    BNT,BNP
         integer    NSST,NSSP


         real*8   S(NST,NSP,3)
         real*8   Sfmol(NST,NSP,3,3,3)
         real*8   Sn(NST,NSP,3)
         real*8   Sb(NST,NSP,3,3)

         real*8   divh
         real*8   divh2
c      Local variables
         
         integer   i,j,k,TOT,imap(nstm*nspm),jmap(nstm*nspm)

         real*8   invr

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Begin the loop over mesh points
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         TOT = 0
         do i=BNT,NSST
            do j=BNP,NSSP
               TOT = TOT + 1
               imap(TOT) = i   
               jmap(TOT) = j   
            end do
         end do
                  

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Calculate the first derivatives
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c   X derivative

CDIR$ IVDEP
         do i=1,TOT
            Sn(imap(i),jmap(i),1) = Sfmol(imap(i),jmap(i),3,2,2) - 
     $         Sfmol(imap(i),jmap(i),1,2,2)
         end do
CDIR$ IVDEP
         do i=1,TOT
            Sn(imap(i),jmap(i),1) = Sn(imap(i),jmap(i),1)*divh
         end do
c   Y derivative

CDIR$ IVDEP
         do i=1,TOT
            Sn(imap(i),jmap(i),2) = Sfmol(imap(i),jmap(i),2,3,2) - 
     $               Sfmol(imap(i),jmap(i),2,1,2)
         end do
CDIR$ IVDEP
         do i=1,TOT
            Sn(imap(i),jmap(i),2) = Sn(imap(i),jmap(i),2)*divh
         end do

c   Z derivative

CDIR$ IVDEP
         do i=1,TOT
            Sn(imap(i),jmap(i),3) = Sfmol(imap(i),jmap(i),2,2,3)
     $                - Sfmol(imap(i),jmap(i),2,2,1)
         end do

CDIR$ IVDEP
         do i=1,TOT
            Sn(imap(i),jmap(i),3) = Sn(imap(i),jmap(i),3) *divh
         end do
c               write(81,*)i,j,Sn(i,j,1)
c               write(82,*)i,j,Sn(i,j,2)
c               write(83,*)i,j,Sn(i,j,3)


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Calculate the second derivatives
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c   XX derivative
   
CDIR$ IVDEP
         do i=1,TOT
            Sb(imap(i),jmap(i),1,1)=Sfmol(imap(i),jmap(i),1,2,2)+
     $               Sfmol(imap(i),jmap(i),3,2,2)
         end do

CDIR$ IVDEP
         do i=1,TOT
              Sb(imap(i),jmap(i),1,1) = Sb(imap(i),jmap(i),1,1) 
     $               - Sfmol(imap(i),jmap(i),2,2,2)*2.0e0
         end do

CDIR$ IVDEP
         do i=1,TOT
              Sb(imap(i),jmap(i),1,1) = Sb(imap(i),jmap(i),1,1)*divh2
         end do

c   YY derivative
CDIR$ IVDEP
         do i = 1, TOT
            Sb(imap(i),jmap(i),2,2)=Sfmol(imap(i),jmap(i),2,1,2)+
     $               Sfmol(imap(i),jmap(i),2,3,2)
         end do

CDIR$ IVDEP
         do i = 1, TOT
              Sb(imap(i),jmap(i),2,2) = Sb(imap(i),jmap(i),2,2) - 
     $               Sfmol(imap(i),jmap(i),2,2,2)*2.0e0
         end do

CDIR$ IVDEP
         do i = 1, TOT
            Sb(imap(i),jmap(i),2,2) =  Sb(imap(i),jmap(i),2,2) * divh2
         end do

c   ZZ derivative
CDIR$ IVDEP
         do i = 1, TOT
            Sb(imap(i),jmap(i),3,3)=Sfmol(imap(i),jmap(i),2,2,1)+
     $               Sfmol(imap(i),jmap(i),2,2,3)
         end do
CDIR$ IVDEP
         do i=1,TOT
              Sb(imap(i),jmap(i),3,3)= Sb(imap(i),jmap(i),3,3)-
     $                Sfmol(imap(i),jmap(i),2,2,2)*2.0e0
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Calculate cross derivatives
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c   XY derivative
CDIR$ IVDEP
         do i=1,TOT
               Sb(imap(i),jmap(i),1,2)= Sfmol(imap(i),jmap(i),3,3,2) -
     $   Sfmol(imap(i),jmap(i),3,1,2) - Sfmol(imap(i),jmap(i),1,3,2) +
     $                Sfmol(imap(i),jmap(i),1,1,2)
         end do
CDIR$ IVDEP
         do i=1,TOT
            Sb(imap(i),jmap(i),1,2) = Sb(imap(i),jmap(i),1,2)
     $              *0.25e0*divh2
         enddo

c   XY derivative
CDIR$ IVDEP
         do i = 1, TOT
               Sb(imap(i),jmap(i),2,1) = Sb(imap(i),jmap(i),1,2)
         end do

c   XZ derivative
CDIR$ IVDEP
         do i = 1, TOT
            Sb(imap(i),jmap(i),1,3)= Sfmol(imap(i),jmap(i),3,2,3) -
     $   Sfmol(imap(i),jmap(i),3,2,1) - Sfmol(imap(i),jmap(i),1,2,3) +
     $                Sfmol(imap(i),jmap(i),1,2,1)
         end do

CDIR$ IVDEP
         do i = 1, TOT
            Sb(imap(i),jmap(i),1,3) = Sb(imap(i),jmap(i),1,3)
     $       * 0.25e0*divh2
         end do

c   ZX derivative
         do i=1,TOT
               Sb(imap(i),jmap(i),3,1) = Sb(imap(i),jmap(i),1,3)
         end do

c   YZ derivative
CDIR$ IVDEP
         do i = 1, TOT
               Sb(imap(i),jmap(i),2,3)= Sfmol(imap(i),jmap(i),2,3,3) -
     $   Sfmol(imap(i),jmap(i),2,3,1) - Sfmol(imap(i),jmap(i),2,1,3) + 
     $               Sfmol(imap(i),jmap(i),2,1,1)
         end do

CDIR$ IVDEP
         do i=1,TOT
          Sb(imap(i),jmap(i),2,3)=  Sb(imap(i),jmap(i),2,3)*0.25e0*divh2
         end do

c   ZY derivative
         do i=1,TOT
               Sb(imap(i),jmap(i),3,2) = Sb(imap(i),jmap(i),2,3)
         end do


         return
      end
c=======================================================================
c   This routine evaluates the apparent horizon equation.
c
c  Input:
c     S    := Spherical coorinates of mesh points.
c     Sg   := 3-metric tensor.
c     Sdg  := Spacial Derivatives of the 3-metric.
c     Skt  := Extrinsic curvature tensor.
c     Sk   := Trace of the extrinsic curvature.
c     Sn   := First derivative of the surface function.
c     Sb   := Second and Cross derivatives of the surface function.
c      BNT   := Begining theta mesh coordinate.
c      BNP   := Begining phi mesh coordinate.
c      NSST   := Ending theta mesh coordinate.
c      NSSP   := Ending phi mesh coordinate.
c     
c  Output:
c     AH   :=  LHS of the apparent horizon equation. F[\varphi]
c
      subroutine Evaluator(S,Sg,Sdg,Skt,Sk,AH,Sn,Sb,Sfmol,
     #                  BNT,BNP,NSST,NSSP,NST,NSP)

         implicit none

         include 'horizon_parameters.h'
   
         integer   BNT,BNP
         integer   NSST,NSSP
         integer    NST,NSP

         real*8   S(NST,NSP,3)
         real*8   Sg(NST,NSP,3,3)
         real*8   Sdg(NST,NSP,3,3,3)
         real*8   Skt(NST,NSP,3,3)
         real*8   Sk(NST,NSP)

         real*8   AH(NST,NSP)
         real*8   Sn(NST,NSP,3)
         real*8   Sb(NST,NSP,3,3)
         real*8   Sfmol(NST,NSP,3,3,3)

c   Local Variables
   
         integer   i,j,k,m,p,q,t,TOT
         integer   imap(nstm*nspm),jmap(nspm*nstm)

         real*8   detg(nstm*nspm)
         real*8   invg(nstm*nspm,3,3)
         real*8   dinvg(nstm*nspm,3,3,3)
         real*8   mag2(nstm*nspm)
         real*8   term1(nstm*nspm)
         real*8   term2(nstm*nspm)
         real*8   term3(nstm*nspm)
         real*8   term4(nstm*nspm)
         real*8   term5(nstm*nspm)
         real*8   term6(nstm*nspm)

c         write(0,*)'NST = ',NST
c         write(0,*)'NSP = ',NSP
c         write(0,*)'BNT = ',BNT
c         write(0,*)'BNP = ',BNP
c         write(0,*)'NSST = ',NSST
c         write(0,*)'NSSP = ',NSSP

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c   Set up the integer map
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         TOT = 0
         do i=BNT,NSST
            do j=BNP,NSSP
               TOT = TOT + 1
               imap(TOT) = i
               jmap(TOT) = j
            end do
         end do
      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Begin the loop over mesh points
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c         do i=BNT,NSST
c            do j=BNP,NSSP

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Calculate the inverse 3-metric 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c            First calculate the determinant
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CDIR$ IVDEP
         do k=1,TOT
             detg(k) = Sg(imap(k),jmap(k),1,1)*
     #               Sg(imap(k),jmap(k),2,2)*Sg(imap(k),jmap(k),3,3) 
         end do
CDIR$ IVDEP
         do k=1,TOT
            detg(k) = detg(k) -
     #                 Sg(imap(k),jmap(k),1,1)*
     #               Sg(imap(k),jmap(k),2,3)*Sg(imap(k),jmap(k),3,2) 
         end do
CDIR$ IVDEP
         do k=1,TOT
            detg(k) = detg(k) -
     #                 Sg(imap(k),jmap(k),2,1)*
     #               Sg(imap(k),jmap(k),1,2)*Sg(imap(k),jmap(k),3,3)  
         end do
CDIR$ IVDEP
         do k=1,TOT
            detg(k) = detg(k) +
     #                 Sg(imap(k),jmap(k),2,1)*
     #               Sg(imap(k),jmap(k),1,3)*Sg(imap(k),jmap(k),3,2)  
         end do
CDIR$ IVDEP
         do k=1,TOT
            detg(k) = detg(k) +
     #                 Sg(imap(k),jmap(k),3,1)*
     #               Sg(imap(k),jmap(k),1,2)*Sg(imap(k),jmap(k),2,3)  
         end do
CDIR$ IVDEP
         do k=1,TOT
            detg(k) = detg(k) -
     #               Sg(imap(k),jmap(k),3,1)*
     #               Sg(imap(k),jmap(k),1,3)*Sg(imap(k),jmap(k),2,2)
         end do

CDIR$ IVDEP
         do k =1,TOT
             detg(k) = 1.0e0 / detg(k)
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c            Calculate the components of the inverse 3-metric
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
         do k=1,TOT

            invg(k,1,1) = (Sg(imap(k),jmap(k),2,2)*
     &           Sg(imap(k),jmap(k),3,3) - 
     #      Sg(imap(k),jmap(k),2,3)*Sg(imap(k),jmap(k),3,2))*detg(k)

         end do
CDIR$ IVDEP
         do k=1,TOT

            invg(k,1,2) = (Sg(imap(k),jmap(k),1,3)*
     &           Sg(imap(k),jmap(k),3,2) - 
     #      Sg(imap(k),jmap(k),1,2)*Sg(imap(k),jmap(k),3,3))*detg(k)

         end do
CDIR$ IVDEP
         do k=1,TOT
            invg(k,1,3) = (Sg(imap(k),jmap(k),1,2)*
     &           Sg(imap(k),jmap(k),2,3) - 
     #      Sg(imap(k),jmap(k),1,3)*Sg(imap(k),jmap(k),2,2))*detg(k)

         end do
CDIR$ IVDEP
         do k=1,TOT

            invg(k,2,1) = invg(k,1,2)

         end do
CDIR$ IVDEP
         do k=1,TOT

            invg(k,2,2) = (Sg(imap(k),jmap(k),1,1)*
     &           Sg(imap(k),jmap(k),3,3) - 
     #      Sg(imap(k),jmap(k),1,3)*Sg(imap(k),jmap(k),3,1))*detg(k)

         end do
CDIR$ IVDEP
         do k=1,TOT

            invg(k,2,3) = (Sg(imap(k),jmap(k),1,3)*
     &            Sg(imap(k),jmap(k),2,1) - 
     #      Sg(imap(k),jmap(k),1,1)*Sg(imap(k),jmap(k),2,3))*detg(k)

         end do
CDIR$ IVDEP
         do k=1,TOT

            invg(k,3,1) = invg(k,1,3)

         end do
CDIR$ IVDEP
         do k=1,TOT

            invg(k,3,2) = invg(k,2,3)

         end do
CDIR$ IVDEP
         do k=1,TOT

            invg(k,3,3) = (Sg(imap(k),jmap(k),1,1)*
     &           Sg(imap(k),jmap(k),2,2) - 
     #      Sg(imap(k),jmap(k),1,2)*Sg(imap(k),jmap(k),2,1))*detg(k)

         end do
c         do k=1,TOT
c            write(83,*)invg(k,1,1),invg(k,1,2),invg(k,1,3)
c            write(83,*)invg(k,2,2),invg(k,2,3),invg(k,3,3)
c            write(83,*)'***********************'
c             write(83,*)detg(k)

c         end do


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Raise the indices of the derivative of the metric
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

               do k=1,3
                  do m=1,3
                     do p=1,3

                     do i=1,TOT
                        dinvg(i,k,m,p) = 0.0e0
                     end do

                        do q=1,3
                           do t=1,3
                        do i=1,TOT
                        dinvg(i,k,m,p) = dinvg(i,k,m,p) - 
     #   invg(i,k,q)*invg(i,m,t)*Sdg(imap(i),jmap(i),q,t,p)
                        enddo
                        enddo
                        enddo
                     enddo
                  enddo
               enddo

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Evaluate the apparent horizon equation term by term.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Calculate the magnitude squared of the first derviatives.
c      ie: Calculate the normalization
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
            do i=1,TOT
               mag2(i) = 0.0e0
            end do

            do k=1,3
               do m=1,3
                  do i=1,TOT
               mag2(i) = mag2(i) + invg(i,k,m)*
     #               Sn(imap(i),jmap(i),k)*Sn(imap(i),jmap(i),m)
                  enddo
               enddo
c               write(0,*)imap(i),jmap(i),mag2(i)
            end do

CDIR$ IVDEP
            do i=1,TOT
c               write(81,*)mag2(i)
               mag2(i) = 1.0e0/mag2(i)
            end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      First term. 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c            write(0,*)'First term ',i,j
CDIR$ IVDEP
            do i=1,TOT
               term1(i) = 0.0e0
            end do


            do k=1,3
               do m=1,3
                  do i=1,TOT

            term1(i) = term1(i) + invg(i,k,m)*Sb(imap(i),jmap(i),k,m)

                  enddo
               enddo
            end do

c         Divide term 1 as in equation: Remember mag2 has been inverted.
CDIR$ IVDEP
         do i=1,TOT
            term1(i) = term1(i) * sqrt(mag2(i))
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Second term. 
c         Derivative of metric term from cov deriv.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c            write(0,*)'Second term ',i,j
CDIR$ IVDEP
         do i=1,TOT
            term2(i)=0.0e0
         end do
            
         do k=1,3
            do m=1,3
               do i=1,TOT
         term2(i) = term2(i) + dinvg(i,k,m,k)*Sn(imap(i),jmap(i),m)
               enddo
            enddo
         end do

c         Divide as in equation: Remember mag2 has been inverted.

CDIR$ IVDEP
         do i=1,TOT
            term2(i) = term2(i) * sqrt(mag2(i))
         enddo


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Third term. 
c         Derivative of the normalizing factor: first one
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c            write(0,*)'Third term ',i,j
CDIR$ IVDEP
         do i=1,TOT
               term3(i) = 0.0e0
         enddo

         do k=1,3
            do m=1,3
               do p=1,3
                  do q=1,3
                     do i=1,TOT
                        
               term3(i) = term3(i) + invg(i,k,m)*dinvg(i,p,q,k)*
     #   Sn(imap(i),jmap(i),m)*Sn(imap(i),jmap(i),p)*
     &   Sn(imap(i),jmap(i),q)

                        enddo
                     enddo
                  enddo
               enddo   
            end do
            
c      Multiply and divide as in equation: Remember mag2 has been inverted.

         do i=1,TOT
            term3(i) = -0.5e0 * term3(i) * (mag2(i)**(1.5))
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Fourth term. 
c         Derivative of the normalizing factor: second one
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c            write(0,*)'Fourth term ',i,j
CDIR$ IVDEP
         do i=1,TOT
            term4(i) = 0.0e0
         end do
         do k=1,3
            do m=1,3
               do p=1,3
                  do q=1,3

                  do i=1,TOT

               term4(i) = term4(i) + invg(i,k,m)*invg(i,p,q)*
     #   Sn(imap(i),jmap(i),m)*Sb(imap(i),jmap(i),k,p)*
     &   Sn(imap(i),jmap(i),q)

                  enddo
                  enddo
               enddo
            enddo
         enddo
      
c      Divide by the appropriate factor again: Remember mag2 has been inverted.

CDIR$ IVDEP
         do i=1,TOT
            term4(i) = -term4(i) * (mag2(i)**(1.5))
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Fifth term. 
c         Christoffel symbol term
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
         do i=1,TOT
            term5(i) = 0.0e0
         end do
         do k=1,3
            do m=1,3
               do p=1,3
                  do q=1,3
   
                  do i=1,TOT
            term5(i) = term5(i) + invg(i,k,m)*invg(i,p,q)*
     #      Sdg(imap(i),jmap(i),k,m,q)*Sn(imap(i),jmap(i),p)

                  enddo
                  enddo
               enddo
            enddo
         end do

c         Divide the appropriate terms : Remember mag2 has been inverted.

CDIR$ IVDEP
         do i=1,TOT
            term5(i) = term5(i) * sqrt(mag2(i))
         end do
CDIR$ IVDEP
         do i=1,TOT
            term5(i) = term5(i)*0.5e0
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Sixth term. 
c         Extrinsic curvature tensor term 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
         do i=1,TOT
               term6(i) = 0.0e0
         end do

         do k=1,3
            do m=1,3
               do p=1,3
                  do q=1,3
      

                  do i=1,TOT
                  term6(i) = term6(i) + Skt(imap(i),jmap(i),k,p)*
     #   invg(i,k,m)*invg(i,p,q)*Sn(imap(i),jmap(i),m)*
     &   Sn(imap(i),jmap(i),q)

                  end do
                  enddo
               enddo
            enddo
         enddo

c   Divide by the magnitude squared for this term: Remember mag2 has been inverted.
CDIR$ IVDEP
         do i=1,TOT                     
            term6(i) = term6(i) * mag2(i)
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Since the last term is just the trace of the extrinsic  curvature 
c  we put it in as is.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
         do i = 1, TOT
            Sk(imap(i),jmap(i)) = 0.0e0
         end do
         do k = 1,3
            do m = 1, 3
               do i = 1, TOT
                  Sk(imap(i),jmap(i)) = Sk(imap(i),jmap(i)) +
     .       invg(i,k,m)*Skt(imap(i),jmap(i),k,m)
               end do
            end do
         end do

CDIR$ IVDEP
         do i=1,TOT
            AH(imap(i),jmap(i)) = term1(i) + term2(i) + term3(i) + 
     #      term4(i) +term5(i) + term6(i) - Sk(imap(i),jmap(i))  - 
     #      kappa0
         end do

         return
      end
c=======================================================================


      subroutine Derivatives(S,Sfmol,Sn,Sb,divh,divh2,BNT,BNP,
     $                  NSST,NSSP,NST,NSP)

         implicit none
      
         include 'horizon_parameters.h'
         
         integer    NST,NSP
         integer    BNT,BNP
         integer    NSST,NSSP


         real*8   S(NST,NSP,3)
         real*8   Sfmol(NST,NSP,3,3,3)
         real*8   Sn(NST,NSP,3)
         real*8   Sb(NST,NSP,3,3)

         real*8   divh
         real*8   divh2
c      Local variables
         
         integer   i,j,k,TOT,imap(nstm*nspm),jmap(nstm*nspm)

         real*8   invr

         TOT = 0
         do i=BNT,NSST
            do j=BNP,NSSP
               TOT = TOT + 1
               imap(TOT) = i   
               jmap(TOT) = j   
            end do
         end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c      Calculate the first derivatives
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CDIR$ IVDEP
c         do i=BNT,NSST
c            do j=BNP,NSSP
         do i = 1, TOT
            
c   X derivative

               Sn(imap(i),jmap(i),1) = divh*(
     $   Sfmol(imap(i),jmap(i),3,2,2) - Sfmol(imap(i),jmap(i),1,2,2))
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   Y derivative

               Sn(imap(i),jmap(i),2) = divh*(
     $   Sfmol(imap(i),jmap(i),2,3,2) - Sfmol(imap(i),jmap(i),2,1,2))
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   Z derivative

               Sn(imap(i),jmap(i),3) = divh*
     $   (Sfmol(imap(i),jmap(i),2,2,3) - Sfmol(imap(i),jmap(i),2,2,1))
c            enddo
         enddo

c               write(81,*)imap(i),jmap(i),Sn(imap(i),jmap(i),1)
c               write(82,*)imap(i),jmap(i),Sn(imap(i),jmap(i),2)
c               write(83,*)imap(i),jmap(i),Sn(imap(i),jmap(i),3)


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Calculate the second derivatives
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   XX derivative
   
               Sb(imap(i),jmap(i),1,1)=divh2*(
     $   Sfmol(imap(i),jmap(i),1,2,2)+Sfmol(imap(i),jmap(i),3,2,2)
     $            - 2.0e0*Sfmol(imap(i),jmap(i),2,2,2))
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   YY derivative
               Sb(imap(i),jmap(i),2,2)=divh2*(
     $    Sfmol(imap(i),jmap(i),2,1,2)+Sfmol(imap(i),jmap(i),2,3,2)
     $            - 2.0e0*Sfmol(imap(i),jmap(i),2,2,2))
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   ZZ derivative
               Sb(imap(i),jmap(i),3,3)=divh2*
     $   (Sfmol(imap(i),jmap(i),2,2,1)+Sfmol(imap(i),jmap(i),2,2,3)
     $            - 2.0e0*Sfmol(imap(i),jmap(i),2,2,2))
c            enddo
         enddo

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Calculate cross derivatives
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   XY derivative
               Sb(imap(i),jmap(i),1,2)= 0.25e0*divh2*
     $   (Sfmol(imap(i),jmap(i),3,3,2) -
     $   Sfmol(imap(i),jmap(i),3,1,2) - Sfmol(imap(i),jmap(i),1,3,2) 
     $   + Sfmol(imap(i),jmap(i),1,1,2))
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   XY derivative
               Sb(imap(i),jmap(i),2,1) = Sb(imap(i),jmap(i),1,2)
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   XZ derivative
               Sb(imap(i),jmap(i),1,3)= 0.25e0*divh2*
     $   (Sfmol(imap(i),jmap(i),3,2,3) -
     $   Sfmol(imap(i),jmap(i),3,2,1) - Sfmol(imap(i),jmap(i),1,2,3) 
     $   + Sfmol(imap(i),jmap(i),1,2,1))
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   ZX derivative
               Sb(imap(i),jmap(i),3,1) = Sb(imap(i),jmap(i),1,3)
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   YZ derivative
               Sb(imap(i),jmap(i),2,3)= 0.25e0*divh2*
     $   (Sfmol(imap(i),jmap(i),2,3,3) -
     $   Sfmol(imap(i),jmap(i),2,3,1) - Sfmol(imap(i),jmap(i),2,1,3) 
     $      + Sfmol(imap(i),jmap(i),2,1,1))
c            enddo
         enddo

CDIR$ IVDEP
         do i = 1, TOT
c         do i=BNT,NSST
c            do j=BNP,NSSP
c   ZY derivative
               Sb(imap(i),jmap(i),3,2) = Sb(imap(i),jmap(i),2,3)

c            enddo
         enddo

         return
      end
