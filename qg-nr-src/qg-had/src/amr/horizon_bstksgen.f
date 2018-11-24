c-----------------------------------------------------------------------
c Routine to construct 3-metric and extrinsic curvature on a given
c time slice for a single boosted Kerr black hole in the Kerr-Schild
c slicing. The boost is in an arbitrary direction.
c-----------------------------------------------------------------------
      subroutine GenKerrBst(GrX,GrY,GrZ,g11,g12,g13,g22,g23,g33,
     $   k11,k12,k13,k22,k23,k33, msk,a,v,t,hx,hy,hz,NX,NY,NZ)
         implicit none

         integer  NX, NY, NZ

         real*8  hx, hy, hz, a, v,t,
     $   g11(NX,NY,NZ),g12(NX,NY,NZ),g13(NX,NY,NZ),
     $   g22(NX,NY,NZ),g23(NX,NY,NZ),g33(NX,NY,NZ),
     $   k11(NX,NY,NZ),k12(NX,NY,NZ),k13(NX,NY,NZ),
     $   k22(NX,NY,NZ),k23(NX,NY,NZ),k33(NX,NY,NZ),
     $   msk(NX,NY,NZ),GrX(NX), GrY(NY), GrZ(NZ)

c     Local variables

         integer  i,j,k,p

         real*8 invr,rr, max, min, ccx,ccy,ccz,mx,my,mz,mskrad, x,y,z,
     .      LTB, LXB, LYB, LZB,LT, LX,LY,LZ, HH, NN, DXH, DYH, DZH,
     .      vx, vy, vz,
     .      dtR, dxR, dyR, dzR,
     .      DXLX, DYLX, DZLX, 
     .      DXLY, DYLY, DZLY, 
     .      DXLZ, DYLZ, DZLZ,
     .      DTLX, DTLY, DTLZ, DTH, DTLT,
     .      DXLT, DYLT, DZLT,
     .      DXLXB, DYLXB, DZLXB, DXLTB,
     .      DXLYB, DYLYB, DZLYB, DYLTB,
     .      DXLZB, DYLZB, DZLZB ,DZLTB,
     .      DTLXB, DTLYB, DTLZB, DTLTB,
     .      drx,dry,drz, M

         real*8 t1 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,
     .          t2 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t29 ,
     .          t3 ,t30 ,t31 ,t32 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t39 ,
     .          t4 ,t40 ,t41 ,t42 ,t43 ,t44 ,t45 ,t46 ,t47 ,t48 ,t49 ,
     .          t5 ,t50 ,t52 ,t55 ,t57 ,t59 ,t6 ,t60 ,t61 ,t63 ,t64 ,
     .          t66 ,t67 ,t7 ,t72 ,t77 ,t8 ,t9 



c     Functions
         real*8 ksg11,ksg12, ksg13, ksg22, ksg23, ksg33
         real*8 ksk11,ksk12, ksk13, ksk22, ksk23, ksk33

c  Read in grid properties

         write(0,*)'Reading in ksbst.dat'
         open(unit=1,file='ksbst.dat')
         read(1,*)min,max
         read(1,*)ccx,ccy,ccz

c        mask origin and radius
         read(1,*)mskrad
         read(1,*)mx,my,mz

c        angular momentum and boost velocity
         read(1,*)a      ,v

c                       Direction
         read(1,*)vx,vy,vz
         close(unit=1)

         write(0,*)'Read in min, max ',min,max
         write(0,*)'Read in coordinates of hole',ccx,ccy,ccz
         write(2,*)'Read in min, max ',min,max
         write(2,*)'Read in coordinates of hole',ccx,ccy,ccz
         write(2,*)'Mask radius ',mskrad
         write(2,*)'Mask location ',mx,my,mz
         write(0,*)'Angular momentum ',a
         write(2,*)'Angular momentum ',a
         write(0,*)'Boost Velocity ',v
         write(2,*)'Boost Velocity ',v
         write(0,*)'Direction:',vx,vy,vz
         write(2,*)'Direction:',vx,vy,vz

         hx = (max - min) / float(NX-1)
         hy = hx
         hz = hx

c   Set the mass
         M=1.0

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
                  rr = sqrt((GrX(i)-mx)**2 + (GrY(j)-my)**2 +
     .                (GrZ(k)-mz)**2)
                  if(rr .gt. mskrad)then
                     msk(i,j,k) = 1.0
                  else
                     msk(i,j,k) = 0.0
                  end if
               end do
            end do
         end do
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
                  z = grz(k) - ccz


c        Define the spherioidal radial function RR
      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t6 = 1.D0-t5
      t7 = sqrt(t6)
      t9 = x*t7
      t10 = vx*vy
      t12 = y*t7
      t14 = vx*vz
      t16 = z*t7
      t19 = (v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16)**2.D0
      t20 = 1.D0/t6
      t21 = t19*t20
      t26 = vy**2
      t27 = y*t26
      t29 = vz*vy
      t33 = (vy*v*t-t10*x+t10*t9-t27+t27*t7-t12-t29*z+t29*t16)**2.D0
      t34 = t33*t20
      t41 = vz**2
      t42 = z*t41
      t45 = (vz*v*t-t14*x+t14*t9-t29*y+t29*t12-t42+t42*t7-t16)**2.D0
      t46 = t45*t20
      t47 = a**2
      t49 = (t21+t34+t46-t47)**2.D0
      RR = sqrt(2.D0*t21+2.D0*t34+2.D0*t46-2.D0*t47+2.D0*sqrt(t49+4.D0*
     #t46*t47))/2.D0

c        t,x,y,z,Derivatives of spheriodal RR
      t1 = v*vx
      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t6 = 1.D0-t5
      t7 = sqrt(t6)
      t9 = x*t7
      t10 = vx*vy
      t12 = y*t7
      t14 = vx*vz
      t16 = z*t7
      t18 = t1*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t19 = t18**2
      t20 = 1.D0/t6
      t21 = t19*t20
      t22 = vy*v
      t26 = vy**2
      t27 = y*t26
      t29 = vz*vy
      t32 = t22*t-t10*x+t10*t9-t27+t27*t7-t12-t29*z+t29*t16
      t33 = t32**2
      t34 = t33*t20
      t35 = vz*v
      t41 = vz**2
      t42 = z*t41
      t44 = t35*t-t14*x+t14*t9-t29*y+t29*t12-t42+t42*t7-t16
      t45 = t44**2
      t46 = t45*t20
      t47 = a**2
      t48 = t21+t34+t46-t47
      t49 = t48**2
      t52 = sqrt(t49+4.D0*t46*t47)
      t57 = t18*t20*t1
      t59 = t32*t20*t22
      t60 = t44*t20
      t61 = t60*t35
      dtR = 1.D0/sqrt(2.D0*t21+2.D0*t34+2.D0*t46-2.D0*t47+2.D0*t52)*(4.D
     #0*t57+4.D0*t59+4.D0*t61+1.D0/t52*(2.D0*t48*(2.D0*t57+2.D0*t59+2.D0
     #*t61)+8.D0*t60*t47*vz*v))/4.D0

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t6 = 1.D0-t5
      t7 = sqrt(t6)
      t9 = x*t7
      t10 = vy*vx
      t12 = y*t7
      t14 = vz*vx
      t16 = z*t7
      t18 = v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t19 = t18**2
      t20 = 1.D0/t6
      t21 = t19*t20
      t26 = vy**2
      t27 = y*t26
      t29 = vz*vy
      t32 = vy*v*t-t10*x+t10*t9-t27+t27*t7-t12-t29*z+t29*t16
      t33 = t32**2
      t34 = t33*t20
      t41 = vz**2
      t42 = z*t41
      t44 = vz*v*t-t14*x+t14*t9-t29*y+t29*t12-t42+t42*t7-t16
      t45 = t44**2
      t46 = t45*t20
      t47 = a**2
      t48 = t21+t34+t46-t47
      t49 = t48**2
      t52 = sqrt(t49+4.D0*t46*t47)
      t59 = t18*t20*(-t3+t3*t7-t7)
      t63 = t32*t20*(-t10+t10*t7)
      t64 = t44*t20
      t66 = -t14+t14*t7
      t67 = t64*t66
      dxR = 1.D0/sqrt(2.D0*t21+2.D0*t34+2.D0*t46-2.D0*t47+2.D0*t52)*(4.D
     #0*t59+4.D0*t63+4.D0*t67+1.D0/t52*(2.D0*t48*(2.D0*t59+2.D0*t63+2.D0
     #*t67)+8.D0*t64*t47*t66))/4.D0

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t6 = 1.D0-t5
      t7 = sqrt(t6)
      t9 = x*t7
      t10 = vy*vx
      t12 = y*t7
      t14 = vx*vz
      t16 = z*t7
      t18 = v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t19 = t18**2
      t20 = 1.D0/t6
      t21 = t19*t20
      t26 = vy**2
      t27 = y*t26
      t29 = vz*vy
      t32 = vy*v*t-t10*x+t10*t9-t27+t27*t7-t12-t29*z+t29*t16
      t33 = t32**2
      t34 = t33*t20
      t41 = vz**2
      t42 = z*t41
      t44 = vz*v*t-t14*x+t14*t9-t29*y+t29*t12-t42+t42*t7-t16
      t45 = t44**2
      t46 = t45*t20
      t47 = a**2
      t48 = t21+t34+t46-t47
      t49 = t48**2
      t52 = sqrt(t49+4.D0*t46*t47)
      t59 = t18*t20*(-t10+t10*t7)
      t63 = t32*t20*(-t26+t26*t7-t7)
      t64 = t44*t20
      t66 = -t29+t29*t7
      t67 = t64*t66
      dyR = 1.D0/sqrt(2.D0*t21+2.D0*t34+2.D0*t46-2.D0*t47+2.D0*t52)*(4.D
     #0*t59+4.D0*t63+4.D0*t67+1.D0/t52*(2.D0*t48*(2.D0*t59+2.D0*t63+2.D0
     #*t67)+8.D0*t64*t47*t66))/4.D0

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t6 = 1.D0-t5
      t7 = sqrt(t6)
      t9 = x*t7
      t10 = vx*vy
      t12 = y*t7
      t14 = vz*vx
      t16 = z*t7
      t18 = v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t19 = t18**2
      t20 = 1.D0/t6
      t21 = t19*t20
      t26 = vy**2
      t27 = y*t26
      t29 = vz*vy
      t32 = vy*v*t-t10*x+t10*t9-t27+t27*t7-t12-t29*z+t29*t16
      t33 = t32**2
      t34 = t33*t20
      t41 = vz**2
      t42 = z*t41
      t44 = vz*v*t-t14*x+t14*t9-t29*y+t29*t12-t42+t42*t7-t16
      t45 = t44**2
      t46 = t45*t20
      t47 = a**2
      t48 = t21+t34+t46-t47
      t49 = t48**2
      t52 = sqrt(t49+4.D0*t46*t47)
      t59 = t18*t20*(-t14+t14*t7)
      t63 = t32*t20*(-t29+t29*t7)
      t64 = t44*t20
      t66 = -t41+t41*t7-t7
      t67 = t64*t66
      dzR = 1.D0/sqrt(2.D0*t21+2.D0*t34+2.D0*t46-2.D0*t47+2.D0*t52)*(4.D
     #0*t59+4.D0*t63+4.D0*t67+1.D0/t52*(2.D0*t48*(2.D0*t59+2.D0*t63+2.D0
     #*t67)+8.D0*t64*t47*t66))/4.D0

c     HH as a function of boosted coordaintes
      t1 = RR**2
      t4 = t1**2
      t7 = vz*vx
      t9 = v**2
      t10 = 1.D0-t9
      t11 = sqrt(t10)
      t14 = vz*vy
      t18 = vz**2
      t19 = z*t18
      t23 = (vz*v*t-t7*x+t7*x*t11-t14*y+t14*y*t11-t19+t19*t11-z*t11)**2.
     #D0
      t26 = a**2
      HH = M*t1*RR/(t4+t23/t10*t26)

c    Unboosted LT,LX,LY,LZ as a function of boosted coordinates
      LT = 1.0

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t9 = x*t7
      t10 = vx*vy
      t12 = y*t7
      t14 = vx*vz
      t16 = z*t7
      t20 = 1.D0/t7
      t26 = vy**2
      t27 = y*t26
      t29 = vz*vy
      t36 = RR**2
      t37 = a**2
      LX = (-RR*(v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16)*t20-a*
     #(vy*v*t-t10*x+t10*t9-t27+t27*t7-t12-t29*z+t29*t16)*t20)/(t36+t37)

      t3 = vy*vx
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t8 = x*t7
      t10 = vy**2
      t11 = y*t10
      t13 = y*t7
      t14 = vz*vy
      t16 = z*t7
      t20 = 1.D0/t7
      t24 = vx**2
      t25 = x*t24
      t29 = vx*vz
      t36 = RR**2
      t37 = a**2
      LY = (-RR*(vy*v*t-t3*x+t3*t8-t11+t11*t7-t13-t14*z+t14*t16)*t20+a*
     #(v*vx*t-t25+t25*t7-t8-t3*y+t3*t13-t29*z+t29*t16)*t20)/(t36+t37)

      t3 = vz*vx
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t10 = vz*vy
      t14 = vz**2
      t15 = z*t14
      LZ = -(vz*v*t-t3*x+t3*x*t7-t10*y+t10*y*t7-t15+t15*t7-z*t7)/t7/RR

c      Derivatives of unboosted quantities
      t1 = RR**2
      t3 = t1**2
      t6 = vz*vx
      t8 = v**2
      t9 = 1.D0-t8
      t10 = sqrt(t9)
      t13 = vz*vy
      t17 = vz**2
      t18 = z*t17
      t21 = vz*v*t-t6*x+t6*x*t10-t13*y+t13*y*t10-t18+t18*t10-z*t10
      t22 = t21**2
      t23 = 1.D0/t9
      t25 = a**2
      t27 = t3+t22*t23*t25
      t31 = t1*RR
      t33 = t27**2
      dtH = 3.D0*M*t1/t27*dtR-M*t31/t33*(4.D0*t31*dtR+2.D0*t21*t23*t25*v
     #z*v)

      t1 = RR**2
      t3 = t1**2
      t6 = vz*vx
      t8 = v**2
      t9 = 1.D0-t8
      t10 = sqrt(t9)
      t13 = vz*vy
      t17 = vz**2
      t18 = z*t17
      t21 = vz*v*t-t6*x+t6*x*t10-t13*y+t13*y*t10-t18+t18*t10-z*t10
      t22 = t21**2
      t23 = 1.D0/t9
      t25 = a**2
      t27 = t3+t22*t23*t25
      t31 = t1*RR
      t33 = t27**2
      dxH = 3.D0*M*t1/t27*dxR-M*t31/t33*(4.D0*t31*dxR+2.D0*t21*t23*t25*(
     #-t6+t6*t10))

      t1 = RR**2
      t3 = t1**2
      t6 = vz*vx
      t8 = v**2
      t9 = 1.D0-t8
      t10 = sqrt(t9)
      t13 = vz*vy
      t17 = vz**2
      t18 = z*t17
      t21 = vz*v*t-t6*x+t6*x*t10-t13*y+t13*y*t10-t18+t18*t10-z*t10
      t22 = t21**2
      t23 = 1.D0/t9
      t25 = a**2
      t27 = t3+t22*t23*t25
      t31 = t1*RR
      t33 = t27**2
      dyH = 3.D0*M*t1/t27*dyR-M*t31/t33*(4.D0*t31*dyR+2.D0*t21*t23*t25*(
     #-t13+t13*t10))

      t1 = RR**2
      t3 = t1**2
      t6 = vz*vx
      t8 = v**2
      t9 = 1.D0-t8
      t10 = sqrt(t9)
      t13 = vz*vy
      t17 = vz**2
      t18 = z*t17
      t21 = vz*v*t-t6*x+t6*x*t10-t13*y+t13*y*t10-t18+t18*t10-z*t10
      t22 = t21**2
      t23 = 1.D0/t9
      t25 = a**2
      t27 = t3+t22*t23*t25
      t31 = t1*RR
      t33 = t27**2
      dzH = 3.D0*M*t1/t27*dzR-M*t31/t33*(4.D0*t31*dzR+2.D0*t21*t23*t25*(
     #-t17+t17*t10-t10))

      dtLT = 0.0
      dxLT = 0.0
      dyLT = 0.0
      dzLT = 0.0

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t9 = x*t7
      t10 = vx*vy
      t12 = y*t7
      t14 = vx*vz
      t16 = z*t7
      t18 = v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t20 = 1.D0/t7
      t29 = RR**2
      t30 = a**2
      t31 = t29+t30
      t40 = vy**2
      t41 = y*t40
      t43 = vz*vy
      t50 = t31**2
      dtLX =(-dtR*t18*t20-RR*v*vx*t20-a*vy*v*t20)/t31-2.D0*(-RR*t18*t20-
     #a*(vy*v*t-t10*x+t10*t9-t41+t41*t7-t12-t43*z+t43*t16)*t20)/t50*RR*d
     #tR

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t9 = x*t7
      t10 = vy*vx
      t12 = y*t7
      t14 = vx*vz
      t16 = z*t7
      t18 = v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t20 = 1.D0/t7
      t31 = RR**2
      t32 = a**2
      t33 = t31+t32
      t42 = vy**2
      t43 = y*t42
      t45 = vz*vy
      t52 = t33**2
      dxLX= (-dxR*t18*t20-RR*(-t3+t3*t7-t7)*t20-a*(-t10+t10*t7)*t20)/t33
     #-2.D0*(-RR*t18*t20-a*(vy*v*t-t10*x+t10*t9-t43+t43*t7-t12-t45*z+t45
     #*t16)*t20)/t52*RR*dxR

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t9 = x*t7
      t10 = vy*vx
      t12 = y*t7
      t14 = vx*vz
      t16 = z*t7
      t18 = v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t20 = 1.D0/t7
      t26 = vy**2
      t32 = RR**2
      t33 = a**2
      t34 = t32+t33
      t43 = y*t26
      t45 = vz*vy
      t52 = t34**2
      dyLX= (-dyR*t18*t20-RR*(-t10+t10*t7)*t20-a*(-t26+t26*t7-t7)*t20)/t
     #34-2.D0*(-RR*t18*t20-a*(vy*v*t-t10*x+t10*t9-t43+t43*t7-t12-t45*z+t
     #45*t16)*t20)/t52*RR*dyR

      t3 = vx**2
      t4 = x*t3
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t9 = x*t7
      t10 = vx*vy
      t12 = y*t7
      t14 = vz*vx
      t16 = z*t7
      t18 = v*vx*t-t4+t4*t7-t9-t10*y+t10*t12-t14*z+t14*t16
      t20 = 1.D0/t7
      t26 = vz*vy
      t32 = RR**2
      t33 = a**2
      t34 = t32+t33
      t43 = vy**2
      t44 = y*t43
      t52 = t34**2
      dzLX= (-dzR*t18*t20-RR*(-t14+t14*t7)*t20-a*(-t26+t26*t7)*t20)/t34-
     #2.D0*(-RR*t18*t20-a*(vy*v*t-t10*x+t10*t9-t44+t44*t7-t12-t26*z+t26*
     #t16)*t20)/t52*RR*dzR

      t3 = vy*vx
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t8 = x*t7
      t10 = vy**2
      t11 = y*t10
      t13 = y*t7
      t14 = vz*vy
      t16 = z*t7
      t18 = vy*v*t-t3*x+t3*t8-t11+t11*t7-t13-t14*z+t14*t16
      t20 = 1.D0/t7
      t29 = RR**2
      t30 = a**2
      t31 = t29+t30
      t38 = vx**2
      t39 = x*t38
      t43 = vx*vz
      t50 = t31**2
      dtLY= (-dtR*t18*t20-RR*vy*v*t20+a*v*vx*t20)/t31-2.D0*(-RR*t18*t20+
     #a*(v*vx*t-t39+t39*t7-t8-t3*y+t3*t13-t43*z+t43*t16)*t20)/t50*RR*dtR

      t3 = vy*vx
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t8 = x*t7
      t10 = vy**2
      t11 = y*t10
      t13 = y*t7
      t14 = vz*vy
      t16 = z*t7
      t18 = vy*v*t-t3*x+t3*t8-t11+t11*t7-t13-t14*z+t14*t16
      t20 = 1.D0/t7
      t26 = vx**2
      t32 = RR**2
      t33 = a**2
      t34 = t32+t33
      t41 = x*t26
      t45 = vx*vz
      t52 = t34**2
      dxLY= (-dxR*t18*t20-RR*(-t3+t3*t7)*t20+a*(-t26+t26*t7-t7)*t20)/t34
     #-2.D0*(-RR*t18*t20+a*(v*vx*t-t41+t41*t7-t8-t3*y+t3*t13-t45*z+t45*t
     #16)*t20)/t52*RR*dxR

      t3 = vy*vx
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t8 = x*t7
      t10 = vy**2
      t11 = y*t10
      t13 = y*t7
      t14 = vz*vy
      t16 = z*t7
      t18 = vy*v*t-t3*x+t3*t8-t11+t11*t7-t13-t14*z+t14*t16
      t20 = 1.D0/t7
      t31 = RR**2
      t32 = a**2
      t33 = t31+t32
      t40 = vx**2
      t41 = x*t40
      t45 = vx*vz
      t52 = t33**2
      dyLY= (-dyR*t18*t20-RR*(-t10+t10*t7-t7)*t20+a*(-t3+t3*t7)*t20)/t33
     #-2.D0*(-RR*t18*t20+a*(v*vx*t-t41+t41*t7-t8-t3*y+t3*t13-t45*z+t45*t
     #16)*t20)/t52*RR*dyR

      t3 = vy*vx
      t5 = v**2
      t7 = sqrt(1.D0-t5)
      t8 = x*t7
      t10 = vy**2
      t11 = y*t10
      t13 = y*t7
      t14 = vz*vy
      t16 = z*t7
      t18 = vy*v*t-t3*x+t3*t8-t11+t11*t7-t13-t14*z+t14*t16
      t20 = 1.D0/t7
      t26 = vz*vx
      t32 = RR**2
      t33 = a**2
      t34 = t32+t33
      t41 = vx**2
      t42 = x*t41
      t52 = t34**2
      dzLY= (-dzR*t18*t20-RR*(-t14+t14*t7)*t20+a*(-t26+t26*t7)*t20)/t34-
     #2.D0*(-RR*t18*t20+a*(v*vx*t-t42+t42*t7-t8-t3*y+t3*t13-t26*z+t26*t1
     #6)*t20)/t52*RR*dzR

      t1 = vz*v
      t2 = v**2
      t4 = sqrt(1.D0-t2)
      t5 = 1.D0/t4
      t10 = vz*vx
      t14 = vz*vy
      t18 = vz**2
      t19 = z*t18
      t24 = RR**2
      dtLZ= -t1*t5/RR+(t1*t-t10*x+t10*x*t4-t14*y+t14*y*t4-t19+t19*t4-z*t
     #4)*t5/t24*dtR

      t1 = vz*vx
      t2 = v**2
      t4 = sqrt(1.D0-t2)
      t7 = 1.D0/t4
      t16 = vz*vy
      t20 = vz**2
      t21 = z*t20
      t26 = RR**2
      dxLZ= -(-t1+t1*t4)*t7/RR+(vz*v*t-t1*x+t1*x*t4-t16*y+t16*y*t4-t21+t
     #21*t4-z*t4)*t7/t26*dxR

      t1 = vz*vy
      t2 = v**2
      t4 = sqrt(1.D0-t2)
      t7 = 1.D0/t4
      t13 = vz*vx
      t20 = vz**2
      t21 = z*t20
      t26 = RR**2
      dyLZ= -(-t1+t1*t4)*t7/RR+(vz*v*t-t13*x+t13*x*t4-t1*y+t1*y*t4-t21+t
     #21*t4-z*t4)*t7/t26*dyR

      t1 = vz**2
      t2 = v**2
      t4 = sqrt(1.D0-t2)
      t7 = 1.D0/t4
      t13 = vz*vx
      t17 = vz*vy
      t21 = z*t1
      t26 = RR**2
      dzLZ= -(-t1+t1*t4-t4)*t7/RR+(vz*v*t-t13*x+t13*x*t4-t17*y+t17*y*t4-
     #t21+t21*t4-z*t4)*t7/t26*dzR

c     Boosted LT,LX, LY,LZ components

      LTB =-(-LT+v*vx*LX+v*vy*LY+v*vz*LZ)/sqrt(1-v**2)

      LXB= -(v*vx*LT-LX*vx**2+LX*vx**2*sqrt(1
     #-v**2)-LX*sqrt(1-v**2)-vx*vy*LY+vx*vy*LY
     #*sqrt(1-v**2)-vx*vz*LZ+vx*vz*LZ*sqrt(1-v**2))/sqrt(1-v**2)

      LYB= -(v*vy*LT-vy*vx*LX+vy*vx*LX*sqrt(1-v**2)-LY*vy**2+
     #LY*vy**2*sqrt(1-v**2)-LY*sqrt(1-v**2)-vz*vy*LZ+vz*vy*LZ
     #*sqrt(1-v**2))/sqrt(1-v**2)

      LZB= -(vz*v*LT-vz*vx*LX+vz*vx*LX*sqrt(1-v**2)-vz*vy*LY+vz*vy
     #*LY*sqrt(1-v**2)-LZ*vz**2+LZ*vz**2*sqrt(1-v**2)-LZ*
     #sqrt(1-v**2))/sqrt(1-v**2)

c      Derivatives of boosted LT, LX, LY, LZ

      t14 = v**2
      dtLTB = -(-dtLT+v*vx*dtLX+v*vy*dtLY+v*
     #vz*dtLZ)/sqrt(1.D0-t14)

      t14 = v**2
      dxLTB = -(-dxLT+v*vx*dxLX+v*vy*dxLY+v*
     #vz*dxLZ)/sqrt(1.D0-t14)

      t14 = v**2
      dyLTB = -(-dyLT+v*vx*dyLX+v*vy*dyLY+v*
     #vz*dyLZ)/sqrt(1.D0-t14)

      t14 = v**2
      dzLTB= -(-dzLT+v*vx*dzLX+v*vy*dzLY+v*
     #vz*dzLZ)/sqrt(1.D0-t14)

      t5 = vx**2
      t6 = dtLX*t5
      t7 = v**2
      t9 = sqrt(1.D0-t7)
      t12 = vx*vy
      t14 = dtLY
      t18 = vx*vz
      t20 = dtLZ
      dtLXB = -(v*vx*dtLT-t6+t6*t9-dtLX*t9-t12*t14+t12*t14*
     #t9-t18*t20+t18*t20*t9)/t9

      t5 = vx**2
      t6 = dxLX*t5
      t7 = v**2
      t9 = sqrt(1.D0-t7)
      t12 = vx*vy
      t14 = dxLY
      t18 = vx*vz
      t20 = dxLZ
      dxLXB = -(v*vx*dxLT-t6+t6*t9-dxLX*t9-t12*t14+t12*t14*
     #t9-t18*t20+t18*t20*t9)/t9

      t5 = vx**2
      t6 = dyLX*t5
      t7 = v**2
      t9 = sqrt(1.D0-t7)
      t12 = vx*vy
      t14 = dyLY
      t18 = vx*vz
      t20 = dyLZ
      dyLXB = -(v*vx*dyLT-t6+t6*t9-dyLX*t9-t12*t14+t12*t14*
     #t9-t18*t20+t18*t20*t9)/t9

      t5 = vx**2
      t6 = dzLX*t5
      t7 = v**2
      t9 = sqrt(1.D0-t7)
      t12 = vx*vy
      t14 = dzLY
      t18 = vx*vz
      t20 = dzLZ
      dzLXB = -(v*vx*dzLT-t6+t6*t9-dzLX*t9-t12*t14+t12*t14*
     #t9-t18*t20+t18*t20*t9)/t9

      t5 = vy*vx
      t7 = dtLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vy**2
      t15 = dtLY*t14
      t18 = vz*vy
      t20 = dtLZ
      dtLYB = -(v*vy*dtLT-t5*t7+t5*t11*t7-t15+t15*t11-dtLY*
     #t11-t18*t20+t18*t20*t11)/t11

      t5 = vy*vx
      t7 = dxLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vy**2
      t15 = dxLY*t14
      t18 = vz*vy
      t20 = dxLZ
      dxLYB = -(v*vy*dxLT-t5*t7+t5*t11*t7-t15+t15*t11-dxLY*
     #t11-t18*t20+t18*t20*t11)/t11

      t5 = vy*vx
      t7 = dyLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vy**2
      t15 = dyLY*t14
      t18 = vz*vy
      t20 = dyLZ
      dyLYB = -(v*vy*dyLT-t5*t7+t5*t11*t7-t15+t15*t11-dyLY*
     #t11-t18*t20+t18*t20*t11)/t11

      t5 = vy*vx
      t7 = dzLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vy**2
      t15 = dzLY*t14
      t18 = vz*vy
      t20 = dzLZ
      dzLYB = -(v*vy*dzLT-t5*t7+t5*t11*t7-t15+t15*t11-dzLY*
     #t11-t18*t20+t18*t20*t11)/t11

      t5 = vz*vx
      t7 = dtLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vz*vy
      t16 = dtLY
      t20 = vz**2
      t21 = dtLZ*t20
      dtLZB = -(vz*v*dtLT-t5*t7+t5*t11*t7-t14*t16+t14*t16*t
     #11-t21+t21*t11-dtLZ*t11)/t11

      t5 = vz*vx
      t7 = dxLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vz*vy
      t16 = dxLY
      t20 = vz**2
      t21 = dxLZ*t20
      dxLZB = -(vz*v*dxLT-t5*t7+t5*t11*t7-t14*t16+t14*t16*t
     #11-t21+t21*t11-dxLZ*t11)/t11

      t5 = vz*vx
      t7 = dyLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vz*vy
      t16 = dyLY
      t20 = vz**2
      t21 = dyLZ*t20
      dyLZB = -(vz*v*dyLT-t5*t7+t5*t11*t7-t14*t16+t14*t16*t
     #11-t21+t21*t11-dyLZ*t11)/t11

      t5 = vz*vx
      t7 = dzLX
      t9 = v**2
      t11 = sqrt(1.D0-t9)
      t14 = vz*vy
      t16 = dzLY
      t20 = vz**2
      t21 = dzLZ*t20
      dzLZB = -(vz*v*dzLT-t5*t7+t5*t11*t7-t14*t16+t14*t16*t
     #11-t21+t21*t11-dzLZ*t11)/t11


c      Get the metric components

                 g11(i,j,k) = ksg11(HH, lxb, lyb, lzb)
                 g12(i,j,k) = ksg12(HH, lxb, lyb, lzb)
                 g13(i,j,k) = ksg13(HH, lxb, lyb, lzb)
                 g22(i,j,k) = ksg22(HH, lxb, lyb, lzb)
                 g23(i,j,k) = ksg23(HH, lxb, lyb, lzb)
                 g33(i,j,k) = ksg33(HH, lxb, lyb, lzb)

                 k11(i,j,k) = ksk11(HH, ltB, lxB, lyB, lzB,
     &                        dtH, dxH, dyH, dzH,
     &                        dtltB, dxltB, dyltB, dzltB,
     &                        dtlxB, dxlxB, dylxB, dzlxB,
     &                        dtlyB, dxlyB, dylyB, dzlyB,
     &                        dtlzB, dxlzB, dylzB, dzlzB)

                 k12(i,j,k) = ksk12(HH, ltB, lxB, lyB, lzB,
     &                        dtH, dxH, dyH, dzH,
     &                        dtltB, dxltB, dyltB, dzltB,
     &                        dtlxB, dxlxB, dylxB, dzlxB,
     &                        dtlyB, dxlyB, dylyB, dzlyB,
     &                        dtlzB, dxlzB, dylzB, dzlzB)

                 k13(i,j,k) = ksk13(HH, ltB, lxB, lyB, lzB,
     &                        dtH, dxH, dyH, dzH,
     &                        dtltB, dxltB, dyltB, dzltB,
     &                        dtlxB, dxlxB, dylxB, dzlxB,
     &                        dtlyB, dxlyB, dylyB, dzlyB,
     &                        dtlzB, dxlzB, dylzB, dzlzB)

                 k22(i,j,k) = ksk22(HH, ltB, lxB, lyB, lzB,
     &                        dtH, dxH, dyH, dzH,
     &                        dtltB, dxltB, dyltB, dzltB,
     &                        dtlxB, dxlxB, dylxB, dzlxB,
     &                        dtlyB, dxlyB, dylyB, dzlyB,
     &                        dtlzB, dxlzB, dylzB, dzlzB)

                 k23(i,j,k) = ksk23(HH, ltB, lxB, lyB, lzB,
     &                        dtH, dxH, dyH, dzH,
     &                        dtltB, dxltB, dyltB, dzltB,
     &                        dtlxB, dxlxB, dylxB, dzlxB,
     &                        dtlyB, dxlyB, dylyB, dzlyB,
     &                        dtlzB, dxlzB, dylzB, dzlzB)

                 k33(i,j,k) = ksk33(HH, ltB, lxB, lyB, lzB,
     &                        dtH, dxH, dyH, dzH,
     &                        dtltB, dxltB, dyltB, dzltB,
     &                        dtlxB, dxlxB, dylxB, dzlxB,
     &                        dtlyB, dxlyB, dylyB, dzlyB,
     &                        dtlzB, dxlzB, dylzB, dzlzB)


               end if

               end do
            end do
         end do



         return
      end 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Routines to return 3-metric and extrinsic curvature given the components
c of the null vector in the Kerr-Schild metric: returns a scalar
c
c 3-metric components
c
      function ksg11(H, lx, ly, lz)
        implicit none
        real*8  ksg11
        real*8  H, lx, ly, lz
        
        ksg11 = 1.0 + 2.0*H*lx*lx
      end 

      function ksg12(H, lx, ly, lz)
        implicit none
        real*8  ksg12
        real*8  H, lx, ly, lz
        
        ksg12 = 2.0*H*lx*ly
      end 

      function ksg13(H, lx, ly, lz)
        implicit none
        real*8  ksg13
        real*8  H, lx, ly, lz

        ksg13 = 2.0*H*lx*lz
      end 

      function ksg22(H, lx, ly, lz)
        implicit none
        real*8  ksg22
        real*8  H, lx, ly, lz
        
        ksg22 = 1.0 + 2.0*H*ly*ly
      end 

      function ksg23(H, lx, ly, lz)
        implicit none
        real*8  ksg23
        real*8  H, lx, ly, lz

        ksg23 = 2.0*H*ly*lz
      end 

      function ksg33(H, lx, ly, lz)
        implicit none
        real*8  ksg33
        real*8  H, lx, ly, lz

        ksg33 = 1.0 + 2.0*H*lz*lz
      end 
c
c extrinsic curvature components
c
      function ksk11(HH, lt, lx, ly, lz,          
     &     dtH, dxH, dyH, dzH,       
     &     dtlt, dxlt, dylt, dzlt,   
     &     dtlx, dxlx, dylx, dzlx,  
     &     dtly, dxly, dyly, dzly,   
     &     dtlz, dxlz, dylz, dzlz)
        implicit none
        real*8    ksk11
        real*8    HH, lt, lx, ly, lz,            
     &       dtH, dxH, dyH, dzH,     
     &       dtlt, dxlt, dylt, dzlt,
     &       dtlx, dxlx, dylx, dzlx,
     &       dtly, dxly, dyly, dzly,
     &       dtlz, dxlz, dylz, dzlz

        ksk11 = 
     # (-4*HH**2*LX**3*DTLX-2*HH*LY**2*DTH*LX**2-4*HH**2*LY**2*LX*DT
     #LX+2*HH*LX**3*DXH*LT+4*HH**2*LX**3*DXLT-4*HH**2*LZ*LT*LX*DXLZ+2*HH
     #*LZ*LT*DZH*LX**2+4*HH**2*LZ*LT*LX*DZLX+2*HH*LY*LT*DYH*LX**2+4*HH**
     #2*LY*LT*LX*DYLX-4*HH**2*LY*LT*LX*DXLY-2*HH*LX*DTLX+2*HH*LX*DXLT+2*
     #HH*DXLX*LT+4*HH**2*LZ**2*LX*DXLT-2*HH*LZ**2*DTH*LX**2-4*HH**2*LZ**
     #2*LX*DTLX-2*HH*LX**4*DTH+4*HH**2*LY**2*LX*DXLT+2*DXH*LX*LT-DTH*LX*
     #*2)/(2*HH*LZ**2+1+2*HH*LY**2+2*HH*LX**2-2*HH*LT**2)/sqrt(1+2*HH*LT
     #**2)
      end 

      function ksk12(HH, lt, lx, ly, lz, 
     &     dtH, dxH, dyH, dzH,    
     &     dtlt, dxlt, dylt, dzlt,
     &     dtlx, dxlx, dylx, dzlx,
     &     dtly, dxly, dyly, dzly,
     &     dtlz, dxlz, dylz, dzlz)
        implicit none
        real*8    ksk12
        real*8    HH, lt, lx, ly, lz,            
     &       dtH, dxH, dyH, dzH,     
     &       dtlt, dxlt, dylt, dzlt,
     &       dtlx, dxlx, dylx, dzlx,
     &       dtly, dxly, dyly, dzly,
     &       dtlz, dxlz, dylz, dzlz


        ksk12 = 
     # (HH*LX*DYLT+2*HH**2*LY**3*DXLT-DTH*LX*LY-2*HH**2*LX**3*DTLY+2
     #*HH**2*LX**3*DYLT+DXH*LY*LT+HH*LY*DXLT-2*HH**2*LY**3*DTLX-HH*DTLX*
     #LY+HH*DXLY*LT+HH*DYLX*LT+DYH*LX*LT-HH*LX*DTLY+2*HH**2*LY**2*LX*DYL
     #T-2*HH**2*LY**2*LX*DTLY-2*HH**2*LY**2*DXLY*LT+2*HH**2*LY**2*DYLX*L
     #T+2*HH**2*LX**2*DXLY*LT+2*HH*LY**2*DYH*LX*LT-2*HH**2*LX**2*DTLX*LY
     #+2*HH**2*LX**2*LY*DXLT+2*HH**2*LZ*LT*DZLX*LY+2*HH*LZ*LT*DZH*LX*LY+
     #2*HH**2*LZ*LT*LX*DZLY-2*HH**2*LX**2*DYLX*LT-2*HH**2*LZ*LT*LX*DYLZ+
     #2*HH*LX**2*DXH*LY*LT-2*HH**2*LZ*LT*LY*DXLZ-2*HH**2*LZ**2*LX*DTLY+2
     #*HH**2*LZ**2*LX*DYLT+2*HH**2*LZ**2*LY*DXLT-2*HH**2*LZ**2*DTLX*LY-2
     #*HH*LY**3*DTH*LX-2*HH*LX**3*DTH*LY-2*HH*LZ**2*DTH*LX*LY)/(2*HH*LZ*
     #*2+1+2*HH*LY**2+2*HH*LX**2-2*HH*LT**2)/sqrt(1+2*HH*LT**2)
      end 

      function ksk13(HH, lt, lx, ly, lz,        
     &     dtH, dxH, dyH, dzH,      
     &     dtlt, dxlt, dylt, dzlt, 
     &     dtlx, dxlx, dylx, dzlx,
     &     dtly, dxly, dyly, dzly,
     &     dtlz, dxlz, dylz, dzlz)
        implicit none
        real*8    ksk13
        real*8    HH, lt, lx, ly, lz,
     &       dtH, dxH, dyH, dzH,   
     &       dtlt, dxlt, dylt, dzlt,
     &       dtlx, dxlx, dylx, dzlx,
     &       dtly, dxly, dyly, dzly, 
     &       dtlz, dxlz, dylz, dzlz

        ksk13 = 
     #(-2*HH**2*LY*LT*DXLY*LZ+2*HH**2*LX**2*DXLZ*LT-2*HH**2*LX**2*D
     #ZLX*LT-2*HH**2*LZ**2*DXLZ*LT+2*HH**2*LZ**2*DZLX*LT-2*HH**2*LY*LT*L
     #X*DZLY+2*HH**2*LY*LT*LX*DYLZ+2*HH*LZ**2*DZH*LX*LT+2*HH*LY*LT*DYH*L
     #X*LZ+2*HH*LX**2*DXH*LZ*LT+2*HH**2*LZ**3*DXLT-DTH*LX*LZ+HH*DZLX*LT+
     #HH*DXLZ*LT+DZH*LX*LT+DXH*LZ*LT-2*HH**2*LY**2*LX*DTLZ+HH*LX*DZLT-HH
     #*LX*DTLZ-2*HH**2*LX**3*DTLZ-2*HH**2*LZ**3*DTLX-HH*DTLX*LZ+HH*LZ*DX
     #LT+2*HH**2*LY**2*LX*DZLT+2*HH**2*LY*LT*DYLX*LZ-2*HH**2*LZ**2*LX*DT
     #LZ+2*HH**2*LY**2*LZ*DXLT-2*HH**2*LY**2*DTLX*LZ+2*HH**2*LX**2*LZ*DX
     #LT-2*HH*LY**2*DTH*LX*LZ+2*HH**2*LZ**2*LX*DZLT+2*HH**2*LX**3*DZLT-2
     #*HH**2*LX**2*DTLX*LZ-2*HH*LX**3*DTH*LZ-2*HH*LZ**3*DTH*LX)/(2*HH*LZ
     #**2+1+2*HH*LY**2+2*HH*LX**2-2*HH*LT**2)/sqrt(1+2*HH*LT**2)
      end 

      function ksk22(HH, lt, lx, ly, lz,            
     &     dtH, dxH, dyH, dzH,       
     &     dtlt, dxlt, dylt, dzlt,   
     &     dtlx, dxlx, dylx, dzlx,   
     &     dtly, dxly, dyly, dzly,   
     &     dtlz, dxlz, dylz, dzlz)
        implicit none
        real*8    ksk22
        real*8    HH, lt, lx, ly, lz,            
     &       dtH, dxH, dyH, dzH,       
     &       dtlt, dxlt, dylt, dzlt,   
     &       dtlx, dxlx, dylx, dzlx,   
     &       dtly, dxly, dyly, dzly,  
     &       dtlz, dxlz, dylz, dzlz

        ksk22 = 
     #(4*HH**2*LY**3*DYLT-2*HH*LY**2*DTH*LX**2-4*HH**2*LY*LT*LX*DYL
     #X+4*HH**2*LY*LT*LX*DXLY+2*HH*LY**2*DXH*LX*LT-DTH*LY**2-4*HH**2*LZ*
     #LT*LY*DYLZ+2*HH*LZ*LT*DZH*LY**2+4*HH**2*LZ*LT*LY*DZLY+4*HH**2*LZ**
     #2*LY*DYLT-2*HH*LZ**2*DTH*LY**2-4*HH**2*LZ**2*LY*DTLY+2*HH*LY**3*DY
     #H*LT+4*HH**2*LX**2*LY*DYLT-4*HH**2*LX**2*LY*DTLY-2*HH*LY**4*DTH-4*
     #HH**2*LY**3*DTLY+2*DYH*LY*LT+2*HH*DYLY*LT+2*HH*LY*DYLT-2*HH*LY*DTL
     #Y)/(2*HH*LZ**2+1+2*HH*LY**2+2*HH*LX**2-2*HH*LT**2)/sqrt(1+2*HH*LT*
     #*2)
      end 

      function ksk23(HH, lt, lx, ly, lz,        
     &     dtH, dxH, dyH, dzH,      
     &     dtlt, dxlt, dylt, dzlt,   
     &     dtlx, dxlx, dylx, dzlx,   
     &     dtly, dxly, dyly, dzly,   
     &     dtlz, dxlz, dylz, dzlz)
        implicit none
        real*8    ksk23
        real*8    HH, lt, lx, ly, lz,            
     &       dtH, dxH, dyH, dzH,       
     &       dtlt, dxlt, dylt, dzlt,   
     &       dtlx, dxlx, dylx, dzlx,   
     &       dtly, dxly, dyly, dzly,   
     &       dtlz, dxlz, dylz, dzlz

        ksk23 = 
     #(-2*HH**2*LY**3*DTLZ+2*HH**2*LX*LT*DXLY*LZ-DTH*LY*LZ-2*HH*LZ*
     #*3*DTH*LY-2*HH**2*LZ**2*LY*DTLZ-2*HH**2*LX*LT*DYLX*LZ+2*HH**2*LY**
     #3*DZLT+2*HH*LZ**2*DZH*LY*LT+2*HH**2*LZ**2*DZLY*LT+2*HH*LX*LT*DXH*L
     #Y*LZ+2*HH**2*LY**2*DYLZ*LT+2*HH*LY**2*DYH*LZ*LT+HH*LZ*DYLT+HH*DYLZ
     #*LT+HH*DZLY*LT+DZH*LY*LT+2*HH**2*LZ**3*DYLT-2*HH**2*LX**2*LY*DTLZ-
     #2*HH**2*LY**2*DTLY*LZ+2*HH**2*LX**2*LZ*DYLT-2*HH*LY**3*DTH*LZ+2*HH
     #**2*LZ**2*LY*DZLT+2*HH**2*LY**2*LZ*DYLT-2*HH*LX**2*DTH*LY*LZ-2*HH*
     #*2*LX**2*DTLY*LZ-2*HH**2*LX*LT*DZLX*LY+2*HH**2*LX**2*LY*DZLT-2*HH*
     #*2*LZ**2*DYLZ*LT+DYH*LZ*LT-HH*LY*DTLZ+HH*LY*DZLT-HH*DTLY*LZ-2*HH**
     #2*LZ**3*DTLY-2*HH**2*LY**2*DZLY*LT+2*HH**2*LX*LT*LY*DXLZ)/(2*HH*LZ
     #**2+1+2*HH*LY**2+2*HH*LX**2-2*HH*LT**2)/sqrt(1+2*HH*LT**2)

      end 

      function ksk33(HH, lt, lx, ly, lz, 
     &     dtH, dxH, dyH, dzH,     
     &     dtlt, dxlt, dylt, dzlt,
     &     dtlx, dxlx, dylx, dzlx,
     &     dtly, dxly, dyly, dzly,
     &     dtlz, dxlz, dylz, dzlz)
        implicit none
        real*8    ksk33
        real*8    HH, lt, lx, ly, lz,            
     &       dtH, dxH, dyH, dzH,   
     &       dtlt, dxlt, dylt, dzlt, 
     &       dtlx, dxlx, dylx, dzlx,  
     &       dtly, dxly, dyly, dzly,   
     &       dtlz, dxlz, dylz, dzlz

        ksk33 = 
     #(-DTH*LZ**2-2*HH*LZ*DTLZ+4*HH**2*LZ*LT*LY*DYLZ-4*HH**2*LZ*LT*
     #LY*DZLY+2*HH*LZ**2*DYH*LY*LT-2*HH*LZ**2*DTH*LY**2+2*DZH*LZ*LT+2*HH
     #*DZLZ*LT+2*HH*LZ*DZLT+4*HH**2*LZ**3*DZLT-2*HH*LZ**4*DTH-4*HH**2*LZ
     #**3*DTLZ+2*HH*LZ**3*DZH*LT+4*HH**2*LY**2*LZ*DZLT-4*HH**2*LY**2*LZ*
     #DTLZ+4*HH**2*LX**2*LZ*DZLT-2*HH*LZ**2*DTH*LX**2-4*HH**2*LX**2*LZ*D
     #TLZ-4*HH**2*LZ*LT*LX*DZLX+2*HH*LX*LT*DXH*LZ**2+4*HH**2*LZ*LT*LX*DX
     #LZ)/(2*HH*LZ**2+1+2*HH*LY**2+2*HH*LX**2-2*HH*LT**2)/sqrt(1+2*HH*LT
     #**2)
      end 
