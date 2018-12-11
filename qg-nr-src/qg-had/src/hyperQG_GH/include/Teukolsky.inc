
!     the funcions frr,.. as function of theta and phi
      if (mult .EQ. 0) then
        frr  = 2.0d0 - 3.0d0*(dsin(theta))**2
	frt  = - 3.0d0*dsin(theta)*dcos(theta)
	frp  = 0.0d0
	ftt1 = 3.0d0*(dsin(theta))**2
	ftt2 = -1.0d0
	ftp  = 0.0d0
	fpp1 = - ftt1
	fpp2 = 3.0d0*(dsin(theta))**2 - 1.0d0
      elseif (mult .EQ. 2) then
        write(*,*) "this multipole is not implemented"
	stop      
      else
        write(*,*) "this multipole is not implemented"
	stop
      end if	 

      
!     the functions A,B,C at t=0 as function of eps,r,t0,lambda
      A = 3*amp*exp(-(t0+eps*r)**2/lambda**2)*(-6*r**2*t0*lambda**2+4*r&
     &**2*t0**3+12*r**3*t0**2*eps+12*r**4*t0*eps**2-6*r**3*eps*lambda**2&
     &+4*r**5*eps**3+6*eps*r*lambda**2*t0**2+12*eps**2*r**2*lambda**2*t0&
     &+6*eps**3*r**3*lambda**2+3*lambda**4*t0)/lambda**6/r**5+0.d0
      B = -2*amp*exp(-(t0+eps*r)**2/lambda**2)*(-6*eps*r**3*lambda**4+6&
     &*eps*r**3*lambda**2*t0**2-6*eps**2*r**4*lambda**2*t0-6*eps**3*r**5&
     &*lambda**2+4*eps*r**3*t0**4+16*eps**2*r**4*t0**3+24*eps**3*r**5*t0&
     &**2+16*eps**4*r**6*t0+4*eps**5*r**7-9*lambda**4*r**2*t0+6*lambda**&
     &2*r**2*t0**3+6*eps*lambda**4*r*t0**2+12*eps**2*lambda**4*r**2*t0+6&
     &*eps**3*lambda**4*r**3+3*lambda**6*t0)/lambda**8/r**5+0.d0
      C = amp*exp(-(t0+eps*r)**2/lambda**2)*(-80*r**4*lambda**2*t0**3+1&
     &60*r**7*t0**2*eps**3+16*r**4*t0**5+16*r**9*eps**5+21*lambda**8*t0-&
     &240*r**5*eps*lambda**2*t0**2-240*r**6*eps**2*lambda**2*t0+60*r**4*&
     &lambda**4*t0-12*r**5*eps**3*lambda**4+16*r**7*eps**5*lambda**2-80*&
     &r**7*eps**3*lambda**2+60*r**5*eps*lambda**4+80*r**8*t0*eps**4+160*&
     &r**6*t0**3*eps**2+80*r**5*t0**4*eps+42*eps**3*lambda**6*r**3+36*la&
     &mbda**4*r**2*t0**3-54*lambda**6*r**2*t0-42*r**3*eps*lambda**6+60*r&
     &**3*eps*lambda**4*t0**2+16*r**3*eps*lambda**2*t0**4+64*r**4*eps**2&
     &*lambda**2*t0**3+12*r**4*eps**2*lambda**4*t0+96*r**5*eps**3*lambda&
     &**2*t0**2+64*r**6*eps**4*lambda**2*t0+42*eps*lambda**6*r*t0**2+84*&
     &eps**2*lambda**6*r**2*t0)/lambda**10/r**5/4+0.d0
      dAdt = -3*amp*exp(-(t0+eps*r)**2/lambda**2)*(6*lambda**4*r**2-24*&
     &r**2*lambda**2*t0**2-48*r**3*lambda**2*t0*eps-24*r**4*lambda**2*ep&
     &s**2+8*r**2*t0**4+32*r**3*t0**3*eps+48*r**4*t0**2*eps**2+32*r**5*t&
     &0*eps**3+8*r**6*eps**4-6*eps*r*lambda**4*t0+12*eps*r*lambda**2*t0*&
     &*3+36*eps**2*r**2*lambda**2*t0**2+36*eps**3*r**3*lambda**2*t0-12*e&
     &ps**2*r**2*lambda**4+12*eps**4*r**4*lambda**2-3*lambda**6+6*lambda&
     &**4*t0**2)/lambda**8/r**5+0.d0
      dBdt = 2*amp*exp(-(t0+eps*r)**2/lambda**2)*(-3*lambda**8-42*eps*r&
     &**3*lambda**4*t0+8*eps*r**3*lambda**2*t0**3-48*eps**2*r**4*lambda*&
     &*2*t0**2-72*eps**3*r**5*lambda**2*t0-6*eps*lambda**6*r*t0+12*eps*l&
     &ambda**4*r*t0**3+36*eps**2*lambda**4*r**2*t0**2+36*eps**3*lambda**&
     &4*r**3*t0+8*eps**6*r**8+9*lambda**6*r**2+6*lambda**6*t0**2+8*eps*r&
     &**3*t0**5+40*eps**2*r**4*t0**4+80*eps**3*r**5*t0**3+80*eps**4*r**6&
     &*t0**2+40*eps**5*r**7*t0-6*eps**2*r**4*lambda**4-28*eps**4*r**6*la&
     &mbda**2-36*lambda**4*r**2*t0**2+12*lambda**2*r**2*t0**4-12*eps**2*&
     &lambda**6*r**2+12*eps**4*lambda**4*r**4)/lambda**10/r**5+0.d0
      s1 = -amp*exp(-(t0+eps*r)**2/lambda**2)/4+0.d0
      s3 = 480*r**8*t0**2*eps**4-96*r**4*eps**2*lambda**6+192*r**5*t0**5&
     &*eps+640*r**7*t0**3*eps**3+360*r**6*eps**2*lambda**4+32*r**4*t0**6&
     &+54*lambda**8*r**2-60*lambda**6*r**4+42*lambda**8*t0**2+32*r**10*e&
     &ps**6+192*r**9*t0*eps**5+480*r**6*t0**4*eps**2-88*r**6*eps**4*lamb&
     &da**4-240*r**4*lambda**2*t0**4+360*r**4*lambda**4*t0**2-216*lambda&
     &**6*r**2*t0**2+84*eps**4*lambda**6*r**4-240*r**8*eps**4*lambda**2+&
     &32*r**8*eps**6*lambda**2+72*lambda**4*r**2*t0**4-84*eps**2*lambda*&
     &*8*r**2+720*r**5*eps*lambda**4*t0-960*r**5*eps*lambda**2*t0**3-144&
     &0*r**6*eps**2*lambda**2*t0**2-960*r**7*eps**3*lambda**2*t0-312*r**&
     &3*eps*lambda**6*t0+128*r**3*eps*lambda**4*t0**3-48*r**4*eps**2*lam&
     &bda**4*t0**2-192*r**5*eps**3*lambda**4*t0+32*r**3*eps*lambda**2*t0&
     &**5+160*r**4*eps**2*lambda**2*t0**4+320*r**5*eps**3*lambda**2*t0**&
     &3+320*r**6*eps**4*lambda**2*t0**2+160*r**7*eps**5*lambda**2*t0-42*&
     &eps*lambda**8*r*t0+84*eps*lambda**6*r*t0**3+252*eps**2*lambda**6*r&
     &**2*t0**2+252*eps**3*lambda**6*r**3*t0-21*lambda**10+0.d0
      s4 = 1/lambda**12/r**5+0.d0
      s2 = s3*s4+0.d0
      dCdt = s1*s2




!     the metric and curvature components as function of A,B,C and frr,.. 
      g11 = 1+A*frr
      g12 = r*B*frt
      g13 = r*sin(theta)*B*frp
      g22 = r**2*(1+C*ftt1+A*ftt2)
      g23 = r**2*sin(theta)*ftp*(A-2*C)
      g33 = r**2*sin(theta)**2*(1+C*fpp1+A*fpp2)
      K11 = -0.5D0*dAdt*frr
      K12 = -0.5D0*r*dBdt*frt
      K13 = -0.5D0*r*sin(theta)*dBdt*frp
      K22 = -0.5D0*r**2*(dCdt*ftt1+dAdt*ftt2)
      K23 = -0.5D0*r**2*sin(theta)*ftp*(dAdt-2.D0*dCdt)
      K33 = -0.5D0*r**2*sin(theta)**2*(dCdt*fpp1+dAdt*fpp2)
            
      