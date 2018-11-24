c=======================================================================
c molecule.f
c
c   Written by
c   Mijan Huq
c   Center for Relativity,
c   University of Texas at Austin.
c
c Created Wed Jun  8 14:33:14 CDT 1994 MFH
c   from jacobian.f. 
c Modified Wed Jun 22 22:51:54 CDT 1994 MFH
c   Changing finite difference molecule generation routine. Now called
c   GenMolCf.
c
c   Added in routine to calculate interpolation coefficients. GenMolCf
c   
c
c   Routines to setup a finite difference molecule at each point on a
c   given 2-surface and to calculate interpolation coefficient given
c  the point coordinates on the molecule.
c   
c=======================================================================
c   The following routine generates finite difference molecules at points
c   on a mesh in spherical coordinates.
c   Fri Apr 22 13:43:56 CDT 1994 MFH
c   Modified Thu Jun 23 09:51:12 CDT 1994 MFH
c      Changed the original routine Gen2DMol to only generate the 
c      theta-phi coordinates of the finite difference points.
c      See GenMolCf(...) for generation of interpolation coefficients.
c-----------------------------------------------------------------------   
c   S(NST,NSP,3) := Coordinates of points in spherical coordinates.
c   Sc(NST,NSP,3):= Coordinates of points in cartesian coordinates.
c   Simol(*) := theta-phi coordinates of Finite difference molecule.
c   center(3)    := Coordinates of spherical coordinate center
c   dtheta       := Grid spacing in theta direction.
c   dphi           := Grid spacing in phi direction.
c   h             := Grid spacing for cartesian finite difference molcule
c   NST & NSP    := Grid size for theta-phi mesh.

      subroutine GenMol(S,Sc,Simol,center,dtheta,dphi,h,
     *      BNT,BNP,NSST,NSSP,NST,NSP)
         implicit none

         integer   BNP,BNT,NSST,NSSP,NST, NSP
      
         real*8   S(NST,NSP,3), Sc(NST,NSP,3), 
     *            Simol(NST,NSP,3,3,3,3),
     *            center(3),dtheta, dphi, h

c      Local Variables

         integer   i,j,k,m,p, cin

         real*8   xx(3),tx(3),ax(3),
     *            tmpl(3,3,3,3)

         real*8  Pi

c         write(0,*)NST,NSP
c      Define the value for Pi
         Pi = 4.0e0 * atan(1.0e0)

c      Set up the finite difference molecule

c         call DLOADVEC(tmpl,0.0e0,81)

         do k=1,3
            do m=1,3

c      Neighboring x coordinates
               tmpl(1,1,k,m) =  -1.0e0*h
               tmpl(1,2,k,m) =  0.0e0
               tmpl(1,3,k,m) =   h

c     Neighboring y coordinates
               tmpl(2,k,1,m) =  -1.0e0*h
               tmpl(2,k,2,m) =  0.0e0
               tmpl(2,k,3,m) =   h

c     Neighboring z coordinates
               tmpl(3,k,m,1) =  -1.0e0*h
               tmpl(3,k,m,2) =  0.0e0
               tmpl(3,k,m,3) =   h
            
            enddo
         enddo

c      Loop over all grid points

         do i=BNT,NSST
            do j=BNP,NSSP
c         Copy the values of cartesian coordinates to a local variable
c         Subtract off the center of the coordinate system

               xx(1) = Sc(i,j,1) - center(1)
               xx(2) = Sc(i,j,2) - center(2)
               xx(3) = Sc(i,j,3) - center(3)

c               xx(1)=S(i,j,1)*cos(S(i,j,3))*sin(S(i,j,2))
c               xx(2)=S(i,j,1)*sin(S(i,j,3))*sin(S(i,j,2))
c               xx(3)=S(i,j,1)*cos(S(i,j,2))

c         Set up the finite difference molecule

               do k=1,3
                  do m=1,3
                     do p=1,3



c            Find each of the coordinates of the finite difference points.
                        tx(1) = xx(1) + tmpl(1,k,m,p)
                        tx(2) = xx(2) + tmpl(2,k,m,p)
                        tx(3) = xx(3) + tmpl(3,k,m,p)

c            Convert these coordinates to polar spherical coordinates so
c            as to be able to interpolate in the 2-surface.

                        ax(2) =atan2(sqrt(tx(1)**2 + tx(2)**2),tx(3))

c                        if(tx(1) .eq. 0.0e0 .and. tx(2) .eq. 0.0e0)then
c                           ax(3) = 0.0e0
c                        else
                          ax(3) =2.0e0*Pi + atan2(tx(2),tx(1))
c                        endif


c                        if( ax(3) .ge. 2.0e0*Pi)then
c                           ax(3) = ax(3) - 2.0e0*Pi
c                        endif
                         ax(3) = mod(ax(3),2.0e0*Pi)


                        ax(1) = sqrt(tx(1)**2 + tx(2)**2  + tx(3)**2)

                        


                        Simol(i,j,k,m,p,1) = ax(1)
                        Simol(i,j,k,m,p,2) = ax(2)
                        Simol(i,j,k,m,p,3) = ax(3)

c                        if( k*m*p .eq. 8)then
c                           write(84,*)'C',ax(2),ax(3)
c                           write(84,*)'A',S(i,j,2),S(i,j,3)
                           
c            Added in this set of commands because at theta =0,pi get
c            0 for calculated value of phi vs true value causing problems
c            in the jacobian. Get off-diagonals due to this. 7-17-94.
c
c                           Simol(i,j,k,m,p,1) = S(i,j,1)
c                           Simol(i,j,k,m,p,2) = S(i,j,2)
c                           Simol(i,j,k,m,p,3) = S(i,j,3)
c                        end if




                     enddo
                  enddo
               enddo
            enddo
         enddo

         return
      end         
c=======================================================================
c   Routine GenMolCf
c   Created Thu Jun 23 09:55:49 CDT 1994 MFH
c   This routine will take in finite difference molecule spherical 
c   coordinates and calculate interpolation coefficients for each 
c   of the points.
c-----------------------------------------------------------------------
c   Arguments
c   S(NST,NSP,3) := Coordinates of points in spherical coordinates.
c   Simol(*) := theta-phi coordinates of Finite difference molecule.
c   ipcf(*) := Interpolation coefficients to be returned by subroutine.
c   center(3)    := Coordinates of spherical coordinate center
c   dtheta       := Grid spacing in theta direction.
c   dphi           := Grid spacing in phi direction.
c   h             := Grid spacing for cartesian finite difference molcule
c   NST & NSP    := Grid size for theta-phi mesh.
c      Coefficients checked. 1) Convergence to 4th order. 2) In place diff.
c

      subroutine GenMolCf(S,Simol,ipcf,center,dtheta,dphi,h,
     *      BNT,BNP,NSST,NSSP,NST,NSP)
         implicit none

         integer   NST, NSP,BNT,BNP,NSST,NSSP
      
         real*8   S(NST,NSP,3),
     *            Simol(NST,NSP,3,3,3,3),
     *            ipcf(NST,NSP,27,4,2),
     *            center(3),dtheta, dphi, h

c      Local Variables

         integer   i,j,k,m,p, mi,mm,
     *            nt, np

         real*8   at(4), ap(4), list(4), trem, prem, pt

c      Initialize the interpolation points list
      list(1) = 1.0e0
      list(2) = 2.0e0
      list(3) = 3.0e0
      list(4) = 4.0e0

c      Loop over surface mesh points

         do i=BNT,NSST
          do j=BNP,NSSP

c         Loop through finite difference molecule points.
           mi = 0
           do k=1,3
            do m=1,3
             do p=1,3

              mi = mi + 1   

c            Calculate integer coordinates for point in question.
              nt = int(Simol(i,j,k,m,p,2)/dtheta) + 1
              np = int(Simol(i,j,k,m,p,3)/dphi) + 1

c            Calculate the remainder 
              trem = Simol(i,j,k,m,p,2) - S(nt,np,2)
              prem = Simol(i,j,k,m,p,3) - S(nt,np,3)

c               write(6,*)Simol(i,j,k,m,p,2),S(nt,np,2)

c            Start with theta direction.
              pt = 2.0e0 + trem/dtheta
              call lgrncf4(pt,list,at)
c            if(k*m*p .eq. 8)then
c               ap(1) = 0.0e0
c               ap(2) = 1.0e0
c               ap(3) = 0.0e0
c               ap(4) = 0.0e0
c               at(1) = 0.0e0
c               at(2) = 1.0e0
c               at(3) = 0.0e0
c               at(4) = 0.0e0
c            end if

c            Copy over to coefficient array
               ipcf(i,j,mi,1,1) = at(1)
               ipcf(i,j,mi,2,1) = at(2)
               ipcf(i,j,mi,3,1) = at(3)
               ipcf(i,j,mi,4,1) = at(4)


c            Phi direction.
              pt = 2.0e0 + prem/dphi
              call lgrncf4(pt,list,ap)

c            Copy over to coefficient array
               ipcf(i,j,mi,1,2) = ap(1)
               ipcf(i,j,mi,2,2) = ap(2)
               ipcf(i,j,mi,3,2) = ap(3)
               ipcf(i,j,mi,4,2) = ap(4)

c            if( (i .eq. 2 .and. j .eq. 7) .and. k*m*p .eq. 8)then
c                  write(80,*)1,at(1), ap(1)
c                  write(80,*)2,at(2), ap(2)
c                  write(80,*)3,at(3), ap(3)
c                  write(80,*)4,at(4), ap(4)
c            end if


              end do
             end do
            end do
c         End loop over mesh points
             end do
         end do

         return
      end
               
               

            
c=======================================================================
c
      subroutine SetUpIndex(theta_indx, phi_indx,theta_red,phi_red,
     #          cindx,NEQN,NRED,NST,NSP)
         implicit none

         include 'horizon_parameters.h'

         integer theta_indx(NSTM*NSPM)
         integer phi_indx(NSTM*NSPM)
         integer theta_red(2*NSPM)
         integer phi_red(2*NSPM)

         integer cindx(NSTM,NSPM)

         integer  NEQN
         integer  NRED
         integer  NST,NSP

c  Local Variables

         integer  i,j,k
         integer  kk
         integer  rr


         kk = 1
         rr = 1
         NRED = 0
         do i=1,NST
            do j=1,NSP
               if((i .eq. 1 .or. i .eq. NST)
     #                  .and. (j .ne. int(NST/2)))then
c     #               )then
                  theta_red(rr) = i
                  phi_red(rr) = j

c                 write(83,*)i,j,rr

                  NRED = rr
                  rr = rr + 1

               else

                  theta_indx(kk) = i
                  phi_indx(kk) = j
                  NEQN = kk
c                 write(82,*)i,j,kk

                  kk = kk + 1
               endif
            enddo
         enddo

c      Set up the inverse mapping
         do kk=1,NEQN
            i = theta_indx(kk)
            j = phi_indx(kk)
            if(i .eq. 1 .or. i .eq. NST)then
               do k=1,NSP
                  cindx(i,k) = kk
               end do
            else
               cindx(i,j) = kk
            end if
         end do

c         do i=1,NST
c            do j=1,NSP
c               write(0,*)i,j,cindx(i,j)
c            end do
c         end do

         return
      end
c=======================================================================
c   Routine: IntpCf2D4
c
c   Evaluates the function values at the finite difference molecule
c   points and stores them in Sfmol.
c
c   Input:
c      Sf(*) := Function values defined on mesh.
c      Simol(*) := r,theta,phi coordinates of the FDM points.
c      ipcf(*)  := Interpolation coeficients.
c      NST/NSP  := Mesh size.
c      BNT/NSST := Begin/End integer coordinates of mesh. (Theta)
c      BNP/NSSP := Begin/End integer coordinates of mesh. (Phi)
c   
c   Output:
c      Sfmol(*) := Function values at FDM points.
c
c   Routine tested at level of evaluation of the apparent horizon
c   equation. Got 2nd order convergence. See notes 6-23-94.
c


      subroutine IntpCf2D4(Sf,Sfmol,Simol,ipcf, NiMol,
     *   BNT,BNP,NSST,NSSP,NST, NSP)
         implicit none

         integer   NST, NSP, BNT, BNP, NSST, NSSP,
     *            NiMol(NST,NSP,27,4,4,2)

         real*8   Sf(NST,NSP), ipcf(NST,NSP,27,4,2),
     *            Sfmol(NST,NSP,3,3,3), Simol(NST,NSP,3,3,3,3)

c      Local Variables

         integer   i,j,k,m,p,q,r,mi,nt,np

         real*8   sum,at(4),ap(4)

         do i=BNT,NSST
            do j=BNP,NSSP

               mi = 0   

               do k=1,3
                do m=1,3
                 do p=1,3

                  mi = mi + 1

                  Sfmol(i,j,k,m,p) = 0.0e0

                  do q=1,4
                     at(q) = ipcf(i,j,mi,q,1)
                     ap(q) = ipcf(i,j,mi,q,2)
                  end do
         
                  do q=1,4
                     do r=1,4
         
c                        nt = itmpl(i,j,q,r,1)
c                        np = itmpl(i,j,q,r,2)

                        nt = NiMol(i,j,mi,q,r,1)
                        np = NiMol(i,j,mi,q,r,2)

                        Sfmol(i,j,k,m,p) = Sfmol(i,j,k,m,p) +
     *               at(q)*ap(r)*Sf(nt,np)
                     end do
                  end do

c                     if(k*m*p .eq. 8)then
c                        write(0,*)i,j,Sfmol(i,j,k,m,p),Sf(i,j)
c                        Sfmol(i,j,k,m,p)=Sf(i,j)
c                     endif

c         Subtract F = r - h(t,p)
                  Sfmol(i,j,k,m,p) = Simol(i,j,k,m,p,1) - 
     *               Sfmol(i,j,k,m,p)

                 end do                                  
                end do                                  
               end do                                  

            end do
         end do

         return
      end      
c=======================================================================
c   Routine: GenNbhdFDM
c
c   Copies the intepolation neighborhood from itmpl to an array for 
c   each point on the finite difference molecule.

      subroutine GenNbhdFDM(S,Simol,Nimol,itmpl,dtheta,dphi,
     .                  BNT,BNP,NSST,NSSP,NST,NSP)
         implicit none

         integer   NST,NSP, BNT,BNP,NSST,NSSP,
     .      itmpl(NST,NSP,4,4,2), NiMol(NST,NSP,27,4,4,2)

         real*8   S(NST,NSP,3), Simol(NST,NSP,3,3,3,3),
     .      dtheta, dphi

c      Local Variables

         integer   i,j,k,m,p,q,nt,np,it,ip

         do i = BNT, NSST
          do j = BNP, NSSP
           q = 1
           do k=1,3
            do m=1,3
             do p=1,3
           
            nt = int(Simol(i,j,k,m,p,2)/dtheta) + 1
            np = int(Simol(i,j,k,m,p,3)/dphi) + 1
         
              do it=1,4
               do ip=1,4

            NiMol(i,j,q,it,ip,1) = itmpl(nt,np,it,ip,1)
            NiMol(i,j,q,it,ip,2) = itmpl(nt,np,it,ip,2)

c               if ( i .eq. 5 .and. j .eq. 17 .and.
c     *         k*m*p .eq. 8)then
c                  write(81,*)NiMol(i,j,q,it,ip,1),
c     *            NiMol(i,j,q,it,ip,2)
c                  write(81,*)nt,np,it,ip,q
c               end if


               end do
              end do

              q = q + 1


             end do
            end do
           end do
         
          end do
         end do
         return
      end
c=======================================================================
c   The following routine generates finite difference molecules at points
c   on a mesh in spherical coordinates.
c   Fri Apr 22 13:43:56 CDT 1994 MFH
c   Modified Thu Jun 23 09:51:12 CDT 1994 MFH
c      Changed the original routine Gen2DMol to only generate the 
c      theta-phi coordinates of the finite difference points.
c      See GenMolCf(...) for generation of interpolation coefficients.
c-----------------------------------------------------------------------   
c   S(NST,NSP,3) := Coordinates of points in spherical coordinates.
c   Sc(NST,NSP,3):= Coordinates of points in cartesian coordinates.
c   Simol(*) := theta-phi coordinates of Finite difference molecule.
c   center(3)    := Coordinates of spherical coordinate center
c   dtheta       := Grid spacing in theta direction.
c   dphi           := Grid spacing in phi direction.
c   h             := Grid spacing for cartesian finite difference molcule
c   NST & NSP    := Grid size for theta-phi mesh.

      subroutine IntpMolOLD(S,Sc,Sf,Sfmol,itmpl,
     *      BNT,BNP,NSST,NSSP,NST,NSP)
         implicit none

         include 'horizon_parameters.h'

         integer   BNP,BNT,NSST,NSSP,NST, NSP,
     *            itmpl(NST,NSP,4,4,2)
      
         real*8   S(NST,NSP,3), Sc(NST,NSP,3), Sf(NST,NSP),
     *            Sfmol(NST,NSP,3,3,3)

c      Local Variables

         integer   i,j,k,m,p, cin

         real*8   xx(3),tx(3),ax(3),df,
     *            tmpl(3,3,3,3)


c         write(0,*)NST,NSP
c      Define the value for Pi
         Pi = 4.0e0 * atan(1.0e0)

c      Set up the finite difference molecule

c         call DLOADVEC(tmpl,0.0e0,81)

         do k=1,3
            do m=1,3

c      Neighboring x coordinates
               tmpl(1,1,k,m) =  -1.0e0*h
               tmpl(1,2,k,m) =  0.0e0
               tmpl(1,3,k,m) =   h

c     Neighboring y coordinates
               tmpl(2,k,1,m) =  -1.0e0*h
               tmpl(2,k,2,m) =  0.0e0
               tmpl(2,k,3,m) =   h

c     Neighboring z coordinates
               tmpl(3,k,m,1) =  -1.0e0*h
               tmpl(3,k,m,2) =  0.0e0
               tmpl(3,k,m,3) =   h
            
            enddo
         enddo

c      Loop over all grid points
         do i=BNT,NSST
            do j=BNP,NSSP
               do k=1,3
                  do m=1,3
                     do p=1,3

c         Copy the values of cartesian coordinates to a local variable
c         Subtract off the center of the coordinate system

               xx(1) = Sc(i,j,1) - center(1)
               xx(2) = Sc(i,j,2) - center(2)
               xx(3) = Sc(i,j,3) - center(3)

c               xx(1)=S(i,j,1)*cos(S(i,j,3))*sin(S(i,j,2))
c               xx(2)=S(i,j,1)*sin(S(i,j,3))*sin(S(i,j,2))
c               xx(3)=S(i,j,1)*cos(S(i,j,2))

c         Set up the finite difference molecule
c
c               Sfmol(i,j,2,2,2) = S(i,j,1) - Sf(i,j)




c            Find each of the coordinates of the finite difference points.
                        tx(1) = xx(1) + tmpl(1,k,m,p)
                        tx(2) = xx(2) + tmpl(2,k,m,p)
                        tx(3) = xx(3) + tmpl(3,k,m,p)

c            Convert these coordinates to polar spherical coordinates so
c            as to be able to interpolate in the 2-surface.

                        ax(2) =atan2(sqrt(tx(1)**2 + tx(2)**2),tx(3))

                        if(tx(1) .eq. 0.0e0 .and. tx(2) .eq. 0.0e0)then
                           ax(3) = 0.0e0
                        else
                          ax(3) =2.0e0*Pi + atan2(tx(2),tx(1))
                        endif


                        if( ax(3) .ge. 2.0e0*Pi)then
                           ax(3) = ax(3) - 2.0e0*Pi
                        endif
c                         ax(3) = mod(ax(3),2.0e0*Pi)


                        ax(1) = sqrt(tx(1)**2 + tx(2)**2  + tx(3)**2)

                        
                        call Intp2Dtmpl(S,Sf,ax(3),ax(2),
     *         Sfmol(i,j,k,m,p),df,itmpl,4,dtheta,dphi,NST,NSP,3)


                     
c           Don't forget that Sf is h(theta,phi). So construct
c           F = r - h(theta,phi)
c
                        Sfmol(i,j,k,m,p) = ax(1) - Sfmol(i,j,k,m,p)



                     enddo
                  enddo
               enddo
            enddo
         enddo

         return
      end         
c=======================================================================

        subroutine IntpMol(S,Sc,Sf,Sfmol,itmpl,
     *        BNT,BNP,NSST,NSSP,NST,NSP)
        implicit none

         include 'horizon_parameters.h'

        integer        BNP,BNT,NSST,NSSP,NST, NSP,
     *                 itmpl(NST,NSP,4,4,2)
                
        real*8        S(NST,NSP,3), Sc(NST,NSP,3), Sf(NST,NSP),
     *              Sfmol(NST,NSP,3,3,3)

        integer        i,j,k,m,p, cin, TOT

        integer        ii,jj,np(nstm*nspm),nt(nstm*nspm),
     *                  imap(nstm*nspm),jmap(nstm*nspm)

        real*8        xx(nstm*nspm,3),tx(nstm*nspm,3),ax(nstm*nspm,3),
     *              df,temp,
     *              tmpl(3,3,3,3)
        real*8       rpp(nstm*nspm,4),rtt(nstm*nspm,4),rff(nstm*nspm,4),
     *                  tff(nstm*nspm,4)
        real*8        newp(nstm*nspm),newt(nstm*nspm),prem(nstm*nspm),
     *               trem(nstm*nspm)

        real*8  reci


        reci = 1.0e0 / 6.0e0
         k=1
         do i=BNT,NSST
         do j=BNP,NSSP
            imap(k) = i
            jmap(k) = j
            k = k + 1
         enddo 
         enddo
         TOT =  k - 1

        do i=1,4
         do j=1,TOT
          rpp(j,i) = float(i)
          rtt(j,i) = float(i)
         end do
        end do

        do k=1,3
          do m=1,3

            tmpl(1,k,m,1) =  -1.0e0*h
            tmpl(2,k,m,1) =  0.0e0
            tmpl(3,k,m,1) =   h

            tmpl(k,1,m,2) =  -1.0e0*h
            tmpl(k,2,m,2) =  0.0e0
            tmpl(k,3,m,2) =   h

            tmpl(k,m,1,3) =  -1.0e0*h
            tmpl(k,m,2,3) =  0.0e0
            tmpl(k,m,3,3) =   h
                                
          enddo
        enddo

CDIR$ IVDEP
         do k=1,TOT

            xx(k,1) = Sc(imap(k),jmap(k),1) - center(1)
        end do
CDIR$ IVDEP
         do k=1,TOT
            xx(k,2) = Sc(imap(k),jmap(k),2) - center(2)
        end do
CDIR$ IVDEP
         do k=1,TOT
            xx(k,3) = Sc(imap(k),jmap(k),3) - center(3)
        end do

            do k=1,3
              do m=1,3
                do p=1,3

                 
CDIR$ IVDEP
               do i=1,TOT

                tx(i,1) = xx(i,1) + tmpl(k,m,p,1)
                tx(i,2) = xx(i,2) + tmpl(k,m,p,2)
                tx(i,3) = xx(i,3) + tmpl(k,m,p,3)
               
               end do

c      The following set of lines causes problems with offsets in 
c      the x-direction. No problems in the y or z directions though!
c      Thu Oct 20 00:01:36 CDT 1994 MFH
c      
CDIR$ IVDEP
               do i=1,TOT
                if(tx(i,3) .eq. 0.0e0)then
                  ax(i,2) = Pi*0.5e0
                else

                  ax(i,2) =atan2(sqrt(tx(i,1)**2 + tx(i,2)**2),
     &                             tx(i,3))

                end if
               end do

CDIR$ IVDEP
               do i=1,TOT
                if(tx(i,1) .eq. 0.0e0.and.tx(i,2) .eq. 0.0e0)then
                  ax(i,3) = Pi
                else
                  ax(i,3) =2.0e0*Pi + atan2(tx(i,2),tx(i,1))
                endif
               end do

c
c   Added the -dphi into avoid the case where ax(i,j,3) is 
c   very close to 2Pi. There the x-derivatives have problems
c   since 322 and so forth  are not interpolated correctly.
c   Wed Sep  7 00:50:38 EDT 1994 MFH
cc
c      Changed back got negative ax(i,3)
CDIR$ IVDEP
               do i=1,TOT
                if( ax(i,3) .gt. 2.0e0*Pi - 1.0e-10)then
                  ax(i,3) = ax(i,3) - 2.0e0*Pi
                endif
                if( ax(i,3) .ge. 2.0e0*Pi - 1.0e-10)then
                  ax(i,3) = ax(i,3) - 2.0e0*Pi
                endif
                end do

CDIR$ IVDEP
               do i=1,TOT
                ax(i,1) = sqrt(tx(i,1)**2 
     &                         + tx(i,2)**2  
     &                         + tx(i,3)**2)
                end do

CDIR$ IVDEP
               do i=1,TOT

                np(i)= int(ax(i,3)/dphi) +1
                nt(i)= int(ax(i,2)/dtheta) +1

               end do
   
CDIR$ IVDEP
               do i=1,TOT
                   prem(i) = ax(i,3) - S(nt(i),np(i),3)
                   trem(i) = ax(i,2) - S(nt(i),np(i),2)
               end do

CDIR$ IVDEP
               do i=1,TOT
                newp(i) = prem(i)/dphi + 2.0e0
                newt(i) = trem(i)/dtheta + 2.0e0

               end do

                jj=1
                do ii=1,4
CDIR$ IVDEP
               do i=1,TOT
                rff(i,ii) = Sf(itmpl(nt(i),np(i),ii,jj,1),
     &                     itmpl(nt(i),np(i),ii,jj,2))
                end do
               end do

CDIR$ IVDEP
               do i=1,TOT

         tff(i,jj) = 
     #  -(newt(i)-rtt(i,2))*(newt(i) - rtt(i,3))
     #   *(newt(i) -rtt(i,4))*rff(i,1)*reci
               end do
CDIR$ IVDEP
               do i=1,TOT
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,3))*
     #   (newt(i) -rtt(i,4))*rff(i,2)*0.5e0
               end do
CDIR$ IVDEP
               do i=1,TOT
         tff(i,jj) = tff(i,jj) -
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,4))*rff(i,3)*0.5e0
               end do
CDIR$ IVDEP
               do i=1,TOT
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,3))*rff(i,4)*reci

               end do

                jj=2
                do ii=1,4
CDIR$ IVDEP
               do i=1,TOT
                rff(i,ii) = Sf(itmpl(nt(i),np(i),ii,jj,1),
     &                  itmpl(nt(i),np(i),ii,jj,2))
                end do
               end do
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = 
     #  -(newt(i)-rtt(i,2))*(newt(i) - rtt(i,3))*
     #   (newt(i) -rtt(i,4))*rff(i,1)*reci
         enddo
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,3))*
     #   (newt(i) -rtt(i,4))*rff(i,2)*0.5e0
         enddo
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj) -
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,4))*rff(i,3)*0.5e0
         enddo
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,3))*rff(i,4)*reci
         end do

                jj=3
                do ii=1,4
CDIR$ IVDEP
         do i=1,TOT
                rff(i,ii) = Sf(itmpl(nt(i),np(i),ii,jj,1),
     &                  itmpl(nt(i),np(i),ii,jj,2))
                end do
          end do

CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = 
     #  -(newt(i)-rtt(i,2))*(newt(i) - rtt(i,3))*
     #   (newt(i) -rtt(i,4))*rff(i,1)*reci
         end do
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,3))*
     #   (newt(i) -rtt(i,4))*rff(i,2)*0.5e0
         end do
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj) -
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,4))*rff(i,3)*0.5e0
         end do
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,3))*rff(i,4)*reci
         end do
         

                jj=4
                do ii=1,4
CDIR$ IVDEP
         do i=1,TOT
      
                rff(i,ii) = Sf(itmpl(nt(i),np(i),ii,jj,1),
     &               itmpl(nt(i),np(i),ii,jj,2))
         end do   
                end do

CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = 
     #  -(newt(i)-rtt(i,2))*(newt(i) - rtt(i,3))*
     #   (newt(i) -rtt(i,4))*rff(i,1)*reci
         end do
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,3))*
     #   (newt(i) -rtt(i,4))*rff(i,2)*0.5e0
         end do
CDIR$ IVDEP
         do i=1,TOT
         tff(i,jj) = tff(i,jj)-
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,4))*rff(i,3)*0.5e0
         end do
CDIR$ IVDEP
         do i=1,TOT 
         tff(i,jj) = tff(i,jj) +
     #   (newt(i)-rtt(i,1))*(newt(i) - rtt(i,2))*
     #   (newt(i) -rtt(i,3))*rff(i,4)*reci
         end do

CDIR$ IVDEP
         do i=1,TOT
         Sfmol(imap(i),jmap(i),k,m,p) = 
     #  -(newp(i)-rpp(i,2))*(newp(i) - rpp(i,3))*
     #   (newp(i) -rpp(i,4))*tff(i,1)*reci
         end do
CDIR$ IVDEP
         do i=1,TOT
         Sfmol(imap(i),jmap(i),k,m,p) = Sfmol(imap(i),jmap(i),k,m,p)
     #   +(newp(i)-rpp(i,1))*(newp(i) - rpp(i,3))*
     #   (newp(i) -rpp(i,4))*tff(i,2)*0.5e0
         end do
CDIR$ IVDEP
         do i=1,TOT
         Sfmol(imap(i),jmap(i),k,m,p) = Sfmol(imap(i),jmap(i),k,m,p)
     #   -(newp(i)-rpp(i,1))*(newp(i) - rpp(i,2))*
     #   (newp(i) -rpp(i,4))*tff(i,3)*0.5e0
         end do
CDIR$ IVDEP
         do i=1,TOT
         Sfmol(imap(i),jmap(i),k,m,p) = Sfmol(imap(i),jmap(i),k,m,p)
     #   +(newp(i)-rpp(i,1))*(newp(i) - rpp(i,2))*
     #   (newp(i) -rpp(i,3))*tff(i,4)*reci
         end do

CDIR$ IVDEP
               do i=1,TOT

                Sfmol(imap(i),jmap(i),k,m,p) = ax(i,1) -
     *                      Sfmol(imap(i),jmap(i),k,m,p)


              enddo
            enddo
          enddo
        enddo

        return
        end  

