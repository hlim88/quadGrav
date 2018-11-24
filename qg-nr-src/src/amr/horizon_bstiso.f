c-----------------------------------------------------------------------
c    Routine to generate metric and extrinsic curvature data for a boosted
c    black hole. Coordinates used are isotropic coordinates with the boost
c    carried out in the z-direction.
c
      subroutine GenBstIso(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
     $      k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,tau,NX,NY,NZ)
         implicit none

         integer  NX, NY, NZ

         real*8  hx, hy, hz, tau,
     $   g11(NX,NY,NZ),g12(NX,NY,NZ),g13(NX,NY,NZ),
     $   g22(NX,NY,NZ),g23(NX,NY,NZ),g33(NX,NY,NZ),
     $   k11(NX,NY,NZ),k12(NX,NY,NZ),k13(NX,NY,NZ),
     $   k22(NX,NY,NZ),k23(NX,NY,NZ),k33(NX,NY,NZ),
     $   msk(NX,NY,NZ),GrX(NX), GrY(NY), GrZ(NZ)

c  Local Variables

        integer    i, j, k

      real*8 invr,r, max, min, ccx,ccy,ccz,mx,my,mz,mskrad, x,y,z,
     $        v, t, gamma, chi, psi, dpsi(3), aa, ccm

         t = tau


c  Read in grid properties

         write(0,*)'Reading in bstiso.dat'
         open(unit=1,file='bstiso.dat')
         read(1,*)min,max
         read(1,*)ccx,ccy,ccz

c        mask origin and radius
         read(1,*)mskrad
         read(1,*)mx,my,mz
            read(1,*)v
         close(unit=1)

         write(0,*)'Read in min, max ',min,max
         write(0,*)'Read in coordinates of hole',ccx,ccy,ccz
         write(2,*)'Read in min, max ',min,max
         write(2,*)'Read in coordinates of hole',ccx,ccy,ccz
         write(2,*)'Mask radius ',mskrad
         write(2,*)'Mask location ',mx,my,mz

            write(0,*)'Boost velocity = ',v
            write(2,*)'Boost velocity = ',v

            write(0,*)'Time slice ', t
            write(2,*)'Time slice ', t

         hx = (max - min) / float(NX-1)
         hy = hx
         hz = hx

         aa =1.0
         write(2,*)'Cartesian grid size hx=hy=hz=',hx


c     set up coordinates

         do i = 1, NX
            GrX(i) = hx * float(i-1) + min
            GrY(i) = hx * float(i-1) + min
            GrZ(i) = hx * float(i-1) + min
         end do

c     set up the mask
         do i = 1, NX
            do j = 1, NY
               do k = 1, NZ
                  r = sqrt((GrX(i)-mx)**2 + (GrY(j)-my)**2 +
     .                (GrZ(k)-mz)**2)
                  if(r .gt. mskrad)then
                     msk(i,j,k) = 1.0
                  else
                     msk(i,j,k) = 0.0
                  end if
               end do
            end do
         end do
          call BstIso(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
     $   k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,t,v,ccx,ccy,ccz,aa,
     $    NX,NY,NZ)
            return
        end
c
c-----------------------------------------------------------------------
c    Routine to generate metric and extrinsic curvature data for a boosted
c    black hole. Coordinates used are isotropic coordinates with the boost
c    carried out in the z-direction.
c
      subroutine BstIso(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
     $      k11,k12,k13,k22,k23,k33, msk,hx,hy,hz,tau,v,ccx,ccy,ccz,aa,
     $        NX,NY,NZ)
         implicit none

         integer  NX, NY, NZ

         real*8  hx, hy, hz, tau,
     $   g11(NX,NY,NZ),g12(NX,NY,NZ),g13(NX,NY,NZ),
     $   g22(NX,NY,NZ),g23(NX,NY,NZ),g33(NX,NY,NZ),
     $   k11(NX,NY,NZ),k12(NX,NY,NZ),k13(NX,NY,NZ),
     $   k22(NX,NY,NZ),k23(NX,NY,NZ),k33(NX,NY,NZ),
     $   msk(NX,NY,NZ),GrX(NX), GrY(NY), GrZ(NZ)

c  Local Variables

        integer    i, j, k

      real*8 invr,r, max, min, ccx,ccy,ccz,mx,my,mz,mskrad, x,y,z,
     $        v, t, gamma, chi, dpsi(3), aa, ccm

        real*8    psi

         t = tau

         gamma = 1.0e0 / sqrt(1.0e0 - v*v)

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
                  z = (grz(k) - ccz - v*t)*gamma

                   r = sqrt(x**2 + y**2 + z**2)

                        if(r .eq. 0.0)write(0,*)'zero detected'



                   invr = 1.0e0 / r

c     Define phi and psi ..

                   psi = 1.0e0 + aa * invr
                   chi = 1.0e0 - aa * invr

c     Define derivatives of phi and psi ..

                   dpsi(1) = -aa * x * invr**3
                   dpsi(2) = -aa * y * invr**3
                   dpsi(3) = -aa * z * invr**3


c     Define the 3-metric ..

               g11(i,j,k) = psi**4
                    g12(i,j,k)  =0.0
                    g13(i,j,k)  =0.0
               g22(i,j,k) = psi**4
                    g23(i,j,k)  =0.0
               g33(i,j,k) = (gamma**2)*(psi**4 - v*v*(chi/psi)**2)


c     Define the extrinsic curvature tensor


               ccm = 1.0e0/sqrt((1.0e0 + v*v)*(chi**2 *v*v + psi**6))

                        K11(i,j,k)= 2.0e0*ccm*v*psi*chi*dpsi(3)
                        K12(i,j,k) = 0.0e0
                        K13(i,j,k)= -gamma**2 *ccm*v*psi*dpsi(1)*
     .                        (3.0e0*chi + psi)

                        K22(i,j,k) = 2.0e0*ccm*v*psi*chi*dpsi(3)
                        K23(i,j,k)= -gamma**2 *ccm*v*psi*dpsi(2)*
     .                        (3.0e0*chi + psi)

                   K33(i,j,k) = gamma**2 *dpsi(3)*v*ccm*
     .               (chi**3 *v*v + chi**2 *v*v *psi - 2.0e0*psi**7
     .               -4.0e0*psi**6 *chi)/psi**5

                        end if

                        end do 
                    end do
                end do


         return
      end

