c-----------------------------------------------------------------------
c area.f
c-----------------------------------------------------------------------
c  Routines to compute the area of a two surface given the 3-metric on
c  the surface and the normals to the surface.
c-----------------------------------------------------------------------
c  Uses approach due to Matzner. 
c    Triangulate surface using existing mesh structure.
c    Construct tangents to surface using three metric.
c    Compute cross-products of tangents to get area elements
c    Sum up area elements
c
      real*8 function PropArea1(nst, nsp, Sc, Sg)
         implicit none

         integer   nst, nsp

         real*8   Sc(nst, nsp, 3), Sg(nst, nsp, 3, 3)

c      Local Variables
   
         integer i, j, k
         
         real*8   u(3), v(3), totarea, area1, area2, 
     .            lu, lv, luv, theta

         totarea = 0.0d0

         do i = 1, nst - 1
           do j = 1, nsp - 1
c            Construct upper tangents
             u(1) = Sc(i+1,j,1) - Sc(i,j,1)
             u(2) = Sc(i+1,j,2) - Sc(i,j,2)
             u(3) = Sc(i+1,j,3) - Sc(i,j,3)

             v(1) = Sc(i,j+1,1) - Sc(i,j,1)
             v(2) = Sc(i,j+1,2) - Sc(i,j,2)
             v(3) = Sc(i,j+1,3) - Sc(i,j,3)

c             write(2,*)'u',u(1), u(2), u(3)
c             write(2,*)'v',v(1), v(2), v(3)
c             write(14,*)'Sc',Sc(i,j,1), Sc(i,j,2),Sc(i,j,3)

c            Compute lengths
             lu = sqrt( Sg(i,j,1,1)*u(1)**2 + 
     .                  Sg(i,j,2,2)*u(2)**2 +
     .                  Sg(i,j,3,3)*u(3)**2 +
     .                  Sg(i,j,1,2)*u(1)*u(2)*2.0d0 +
     .                  Sg(i,j,1,3)*u(1)*u(3)*2.0d0 +
     .                  Sg(i,j,2,3)*u(2)*u(3)*2.0d0 )

             lv = sqrt( Sg(i,j,1,1)*v(1)**2 + 
     .                  Sg(i,j,2,2)*v(2)**2 +
     .                  Sg(i,j,3,3)*v(3)**2 +
     .                  Sg(i,j,1,2)*v(1)*v(2)*2.0d0 +
     .                  Sg(i,j,1,3)*v(1)*v(3)*2.0d0 +
     .                  Sg(i,j,2,3)*v(2)*v(3)*2.0d0 )

c            Compute dot-product
             luv = Sg(i,j,1,1) * u(1) * v(1) +
     .             Sg(i,j,1,2) * u(1) * v(2) + 
     .             Sg(i,j,1,3) * u(1) * v(3) + 
     .             Sg(i,j,2,1) * u(2) * v(1) +
     .             Sg(i,j,2,2) * u(2) * v(2) + 
     .             Sg(i,j,2,3) * u(2) * v(3) +
     .             Sg(i,j,3,1) * u(3) * v(1) +
     .             Sg(i,j,3,2) * u(3) * v(2) + 
     .             Sg(i,j,3,3) * u(3) * v(3) 

c            Compute angle
      
c               Zero out pole point contributions since v=0 there...
            if(lu .eq. 0.0d0 .or. lv .eq. 0.0d0)then
               theta = 0.0d0
            else
              theta = acos(luv / (lu*lv))
            end if

c            Compute area element

            area1 = 0.5d0 * lu * lv * sin(theta)

c            Construct lower tangents
             u(1) = Sc(i+1,j+1,1) - Sc(i,j+1,1)
             u(2) = Sc(i+1,j+1,2) - Sc(i,j+1,2)
             u(3) = Sc(i+1,j+1,3) - Sc(i,j+1,3)

             v(1) = Sc(i+1,j,1) - Sc(i+1,j+1,1)
             v(2) = Sc(i+1,j,2) - Sc(i+1,j+1,2)
             v(3) = Sc(i+1,j,3) - Sc(i+1,j+1,3)

c            Compute lengths
             lu = sqrt( Sg(i+1,j+1,1,1)*u(1)**2 + 
     .                  Sg(i+1,j+1,2,2)*u(2)**2 +
     .                  Sg(i+1,j+1,3,3)*u(3)**2 +
     .                  Sg(i+1,j+1,1,2)*u(1)*u(2)*2.0d0 +
     .                  Sg(i+1,j+1,1,3)*u(1)*u(3)*2.0d0 +
     .                  Sg(i+1,j+1,2,3)*u(2)*u(3)*2.0d0 )

             lv = sqrt( Sg(i+1,j+1,1,1)*v(1)**2 + 
     .                  Sg(i+1,j+1,2,2)*v(2)**2 +
     .                  Sg(i+1,j+1,3,3)*v(3)**2 +
     .                  Sg(i+1,j+1,1,2)*v(1)*v(2)*2.0d0 +
     .                  Sg(i+1,j+1,1,3)*v(1)*v(3)*2.0d0 +
     .                  Sg(i+1,j+1,2,3)*v(2)*v(3)*2.0d0 )

c            Compute dot-product
             luv = Sg(i+1,j+1,1,1) * u(1) * v(1) +
     .             Sg(i+1,j+1,1,2) * u(1) * v(2) + 
     .             Sg(i+1,j+1,1,3) * u(1) * v(3) + 
     .             Sg(i+1,j+1,2,1) * u(2) * v(1) +
     .             Sg(i+1,j+1,2,2) * u(2) * v(2) + 
     .             Sg(i+1,j+1,2,3) * u(2) * v(3) +
     .             Sg(i+1,j+1,3,1) * u(3) * v(1) +
     .             Sg(i+1,j+1,3,2) * u(3) * v(2) + 
     .             Sg(i+1,j+1,3,3) * u(3) * v(3) 

c            Compute angle
      
c               Zero out pole point contributions since v=0 there...
            if(lu .eq. 0.0d0 .or. lv .eq. 0.0d0)then
               theta = 0.0d0
            else
              theta = acos(luv / (lu*lv))
            end if

c            Compute area element

            area2 = 0.5d0 * lu * lv * sin(theta)
c            write(1,*)'areas',area1, area2
c            write(1,*)'lu,lv,theta',lu, lv, theta

c            Add to total area

            totarea = totarea + area1 + area2   
           end do
         end do

         PropArea1 = totarea


         return
      end
c-----------------------------------------------------------------------
c The following approach to computing the area uses the normals to the
c surface to project the spherical coordinate 3-metric onto the 2-surface
c and then extracts the measure on the surface. Having done that we integrate
c using trapezoidal integration.
c
      real*8 function AHArea(nst, nsp, dtheta, dphi,
     .                          S, Sc, Sn, Sg, Intg)
         implicit none

         integer   nst, nsp

         real*8   S(nst, nsp, 3), Sc(nst, nsp, 3), 
     .            Sn(nst, nsp, 3), Sg(nst, nsp, 3, 3),
     .            Intg(nst, nsp), dtheta, dphi
 

c      Local Variables
   
         integer i, j, k, m, p, q

         real*8   lambda(3,3), gg(3,3), detg, mag2, invg(3,3),
     .            pp(3,3), sup(3), sdown(3), gp(3,3), trapint,
     .            ssup(3), ssdown(3)

c         write(0,*)'NX/NY', nst, nsp
c         write(0,*)'dtheta/dphi', dtheta, dphi

c      Normalize the derivatives of varphi
         do i = 1, nst
           do j = 1, nsp
c            Determinant of the 3-metric
             detg =  Sg(i,j,1,1)*Sg(i,j,2,2)*Sg(i,j,3,3) -
     #               Sg(i,j,1,1)*Sg(i,j,2,3)*Sg(i,j,3,2)  -
     #               Sg(i,j,2,1)*Sg(i,j,1,2)*Sg(i,j,3,3)  +
     #               Sg(i,j,2,1)*Sg(i,j,1,3)*Sg(i,j,3,2)  +
     #               Sg(i,j,3,1)*Sg(i,j,1,2)*Sg(i,j,2,3)  -
     #               Sg(i,j,3,1)*Sg(i,j,1,3)*Sg(i,j,2,2)

c            write(0,*)'chkpt 1: Determinant detg=',detg
            detg = 1.0d0 / detg

c            Inverse metric components
            invg(1,1) = (Sg(i,j,2,2)*Sg(i,j,3,3) -
     #                   Sg(i,j,2,3)*Sg(i,j,3,2))*detg

            invg(1,2) = (Sg(i,j,1,3)*Sg(i,j,3,2) -
     #                   Sg(i,j,1,2)*Sg(i,j,3,3))*detg

            invg(1,3) = (Sg(i,j,1,2)*Sg(i,j,2,3) -
     #                   Sg(i,j,1,3)*Sg(i,j,2,2))*detg

            invg(2,1) = invg(1,2)

            invg(2,2) = (Sg(i,j,1,1)*Sg(i,j,3,3) -
     #                   Sg(i,j,1,3)*Sg(i,j,3,1))*detg

            invg(2,3) = (Sg(i,j,1,3)*Sg(i,j,2,1) -
     #                   Sg(i,j,1,1)*Sg(i,j,2,3))*detg

            invg(3,1) = invg(1,3)

            invg(3,2) = invg(2,3)

            invg(3,3) = (Sg(i,j,1,1)*Sg(i,j,2,2) -
     #                   Sg(i,j,1,2)*Sg(i,j,2,1))*detg

c            Compute the magnitude

            mag2 = (invg(1,1)*Sn(i,j,1)**2 +
     .             invg(2,2)*Sn(i,j,2)**2 +
     .             invg(3,3)*Sn(i,j,3)**2 +
     .             invg(1,2)*Sn(i,j,1)*Sn(i,j,2) * 2.0d0+
     .             invg(1,3)*Sn(i,j,1)*Sn(i,j,3) * 2.0d0+
     .             invg(2,3)*Sn(i,j,2)*Sn(i,j,3) * 2.0d0)

c            Normalize and store

            sdown(1) = Sn(i,j,1) / sqrt(mag2)
            sdown(2) = Sn(i,j,2) / sqrt(mag2)
            sdown(3) = Sn(i,j,3) / sqrt(mag2)

c            raise an index
   
            sup(1) =  invg(1,1)*sdown(1) + 
     .                invg(1,2)*sdown(2) +
     .                invg(1,3)*sdown(3) 

            sup(2) =  invg(2,1)*sdown(1) + 
     .                invg(2,2)*sdown(2) +
     .                invg(2,3)*sdown(3) 

            sup(3) =  invg(3,1)*sdown(1) + 
     .                invg(3,2)*sdown(2) +
     .                invg(3,3)*sdown(3) 
            

c             construct the transformation matrix
c     X to r, theta, phi
                lambda(1,1) = sin(S(i,j,2))*cos(S(i,j,3))
                lambda(1,2) = S(i,j,1)*cos(S(i,j,2))*cos(S(i,j,3))
                lambda(1,3) = -S(i,j,1)*sin(S(i,j,2))*sin(S(i,j,3))

c     Y to r, theta, phi
                lambda(2,1) = sin(S(i,j,2))*sin(S(i,j,3))
                lambda(2,2) = S(i,j,1)*cos(S(i,j,2))*sin(S(i,j,3))
                lambda(2,3) = S(i,j,1)*sin(S(i,j,2))*cos(S(i,j,3))

c     Z to r, theta, phi
                lambda(3,1) = cos(S(i,j,2))
                lambda(3,2) = -S(i,j,1)*sin(S(i,j,2))
                lambda(3,3) = 0.0d0

c            Now construct the projection tensor

               pp(1,1) = 1.0d0 - sup(1)*sdown(1)
               pp(1,2) =       - sup(1)*sdown(2)
               pp(1,3) =       - sup(1)*sdown(3)

               pp(2,1) =       - sup(2)*sdown(1)
               pp(2,2) = 1.0d0 - sup(2)*sdown(2)
               pp(2,3) =       - sup(2)*sdown(3)

               pp(3,1) =       - sup(3)*sdown(1)
               pp(3,2) =       - sup(3)*sdown(2)
               pp(3,3) = 1.0d0 - sup(3)*sdown(3)

c            Project the cartesian 3-metric

               do k = 1, 3
                 do m = 1, 3
                   gp(k,m) = 0.0d0

                   do p = 1, 3
                     do q = 1,3
                        gp(k,m) = gp(k,m) + 
     .                           pp(p,k)*pp(q,m)*Sg(i,j,p,q)
                     end do
                   end do
                 end do   
               end do

c               write(55,*)'************gp************',i,j
c               write(55,*)gp(1,1), gp(1,2), gp(1,3)
c               write(55,*)gp(2,1), gp(2,2), gp(2,3)
c               write(55,*)gp(3,1), gp(3,2), gp(3,3)

c            Next transform the resulting 
c           Cartesian 3-metric to spherical coordinates
                do k = 1, 3
                  do m = 1,3

                     gg(k,m) = 0.0d0

                     do p =1, 3
                        do q = 1,3
                           gg(k,m) = gg(k,m) +
     *         lambda(p,k)*lambda(q,m)*gp(p,q)
                        end do
                     end do

                  end do
                end do

c               write(56,*)'************gg************',i,j
c               write(56,*)gg(1,1), gg(1,2), gg(1,3)
c               write(56,*)gg(2,1), gg(2,2), gg(2,3)
c               write(56,*)gg(3,1), gg(3,2), gg(3,3)

c            Compute the determinant of the 2-metric 
               detg = gg(2,2)*gg(3,3) - gg(2,3)*gg(3,2)
c            write(0,*)'chkpt 2: Determinant detg=',detg

c               write(99,*)sqrt(detg) - sin(S(i,j,2))

               Intg(i,j) = sqrt(detg)
             
           end do
         end do

c         Now integrate
         AHArea = trapInt(Intg,dtheta,dphi,nst,nsp)

         return
      end
c-----------------------------------------------------------------------
c  The following method is equivalent to the above projection approach.
c  It involves expanding the 3-line element with dr = drho(theta,phi)
      real*8 function AHArea2(nst, nsp, dtheta, dphi,
     .                          S, Sc, Sn, Sg, Intg)
         implicit none

         integer   nst, nsp

         real*8   S(nst, nsp, 3), Sc(nst, nsp, 3), 
     .            Sn(nst, nsp, 3), Sg(nst, nsp, 3, 3),
     .            Intg(nst, nsp), dtheta, dphi
 

c      Local Variables
   
         integer i, j, k, m, p, q,bb,aa

         real*8   lambda(3,3), gg(2,2), detg, omega(3,2), 
     .            dthrho, dphrho, trapInt, gt(3,3)

         do i = 1, nst
            do j = 1, nsp
c      First define the transformation matrix

c            construct the transformation matrix
c     X to r, theta, phi
                lambda(1,1) = sin(S(i,j,2))*cos(S(i,j,3))
                lambda(1,2) = S(i,j,1)*cos(S(i,j,2))*cos(S(i,j,3))
                lambda(1,3) = -S(i,j,1)*sin(S(i,j,2))*sin(S(i,j,3))

c     Y to r, theta, phi
                lambda(2,1) = sin(S(i,j,2))*sin(S(i,j,3))
                lambda(2,2) = S(i,j,1)*cos(S(i,j,2))*sin(S(i,j,3))
                lambda(2,3) = S(i,j,1)*sin(S(i,j,2))*cos(S(i,j,3))

c     Z to r, theta, phi
                lambda(3,1) = cos(S(i,j,2))
                lambda(3,2) = -S(i,j,1)*sin(S(i,j,2))
                lambda(3,3) = 0.0d0

c               write(0,*)'drhodr =',lambda(1,1)*Sn(i,j,1) +
c     .                  lambda(2,1)*Sn(i,j,2) + lambda(3,1)*Sn(i,j,3)

               dthrho = lambda(1,2)*Sn(i,j,1) + 
     .                  lambda(2,2)*Sn(i,j,2) + lambda(3,2)*Sn(i,j,3) 


c               write(0,*)'dthrho =',dthrho

               dphrho = lambda(1,3)*Sn(i,j,1) + 
     .                  lambda(2,3)*Sn(i,j,2) + lambda(3,3)*Sn(i,j,3) 

 
c               write(0,*)'dphrho =',dphrho


c               write(1,*)dthrho, dphrho
c               write(2,*)lambda(1,2),lambda(2,2),lambda(3,2)
c               write(2,*)lambda(1,3),lambda(2,3),lambda(3,3)
c               write(2,*)'************'

c       Compute omega
      
c               omega(1,1) = Sc(i,j,1) / S(i,j,1) * dthrho + lambda(1,2)
c               omega(2,1) = Sc(i,j,2) / S(i,j,1) * dthrho + lambda(2,2)
c               omega(3,1) = Sc(i,j,3) / S(i,j,1) * dthrho + lambda(3,2)

c               omega(1,2) = Sc(i,j,1) / S(i,j,1) * dphrho + lambda(1,3)
c               omega(2,2) = Sc(i,j,2) / S(i,j,1) * dphrho + lambda(2,3)
c               omega(3,2) = Sc(i,j,3) / S(i,j,1) * dphrho + lambda(3,3)

               omega(1,1) = (lambda(1,1) * dthrho + lambda(1,2))
               omega(2,1) = (lambda(2,1) * dthrho + lambda(2,2))
               omega(3,1) = (lambda(3,1) * dthrho + lambda(3,2))

               omega(1,2) = (lambda(1,1) * dphrho + lambda(1,3))
               omega(2,2) = (lambda(2,1) * dphrho + lambda(2,3))
               omega(3,2) = (lambda(3,1) * dphrho + lambda(3,3))

c               omega(1,1) =  dthrho  / S(i,j,1)
c               omega(2,1) = 1.0
c               omega(3,1) = 0.0

c               omega(1,2) = dphrho  / S(i,j,1) / sin(S(i,j,2) + 0.000000001)
c               omega(2,2) = 0.0
c               omega(3,2) = 1.0

c         transform the 3-metric into spherical coordinates
c              do k = 1, 3
c                do m = 1, 3
c                  gt(k,m) = 0.0d0

c                  do p = 1, 3
c                    do q = 1,3
c                       gt(k,m) = gt(k,m) + 
c    .                           lambda(p,k)*lambda(q,m)*Sg(i,j,p,q)
c                    end do
c                  end do
c                end do   
c              end do

c          now project via the basis one-forms
c              do p = 1, 2
c                do q = 1, 2
c                  gg(p,q) = 0.0
c                  do k = 1, 3
c                    do m = 1, 3
c                      gg(p,q) = gg(p,q) + omega(k,p)*omega(m,q)*gt(k,m)
c                    end do
c                  end do
c                end do
c              end do

c      Compute 2-metric

c              gg(1,1) = Sg(i,j,1,1)*omega(1,1) * omega(1,1)
c    .                 + Sg(i,j,2,2)*omega(2,1) * omega(2,1)
c    .                 + Sg(i,j,3,3)*omega(3,1) * omega(3,1)
c    .                 + Sg(i,j,1,2)*omega(1,1) * omega(2,1)*2.0d0
c    .                 + Sg(i,j,1,3)*omega(1,1) * omega(3,1)*2.0d0
c    .                 + Sg(i,j,2,3)*omega(2,1) * omega(3,1)*2.0d0

c              gg(2,2) = Sg(i,j,1,1)*omega(1,2) * omega(1,2)
c    .                 + Sg(i,j,2,2)*omega(2,2) * omega(2,2)
c    .                 + Sg(i,j,3,3)*omega(3,2) * omega(3,2)
c    .                 + Sg(i,j,1,2)*omega(1,2) * omega(2,2)*2.0d0
c    .                 + Sg(i,j,1,3)*omega(1,2) * omega(3,2)*2.0d0
c    .                 + Sg(i,j,2,3)*omega(2,2) * omega(3,2)*2.0d0

c              gg(1,2) = Sg(i,j,1,1)*omega(1,1) * omega(1,2)
c    .                 + Sg(i,j,1,2)*omega(1,1) * omega(2,2)
c    .                 + Sg(i,j,1,3)*omega(1,1) * omega(3,2)
c    .                 + Sg(i,j,2,1)*omega(2,1) * omega(1,2)
c    .                 + Sg(i,j,2,2)*omega(2,1) * omega(2,2)
c    .                 + Sg(i,j,2,3)*omega(2,1) * omega(3,2)
c    .                 + Sg(i,j,3,1)*omega(3,1) * omega(1,2)
c    .                 + Sg(i,j,3,2)*omega(3,1) * omega(2,2)
c    .                 + Sg(i,j,3,3)*omega(3,1) * omega(3,2)

c              gg(2,1) = gg(1,2)

              
                do p = 1, 2
                  do q = 1, 2
                  gg(p,q) = 0.0
                    do aa=1,3
                      do bb=1,3
                        gg(p,q) = gg(p,q) + omega(aa,p)*omega(bb,q) *
     .                                     Sg(i,j,aa,bb)
                      end do
                    end do
                  end do
                end do


c      Compute the determinant

               detg = gg(1,1)*gg(2,2) - gg(1,2)**2

               Intg(i,j) = sqrt(detg)
               
            end do
         end do

c         Now integrate
         AHArea2 = trapInt(Intg,dtheta,dphi,nst,nsp)

         return
      end


c-----------------------------------------------------------------------
      real*8 function trapInt(Itemp,dx,dy,NX,NY)
         implicit none

         integer  NX,NY

         real*8   Itemp(NX,NY), 
     *            dx, dy

c     Local Variables

         integer  i,j

         real*8   dii, II

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

         trapInt = II

         return
      end
c-----------------------------------------------------------------------
