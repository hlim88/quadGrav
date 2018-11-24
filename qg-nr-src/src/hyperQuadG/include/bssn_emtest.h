      t1 = abs(chi)
      inv_chi = 0.1D1/t1
      t1 = gtd(1,1)
      t2 = gtd(2,2)
      t4 = gtd(3,3)
      t6 = gtd(2,3)
      t7 = t6**2
      t9 = gtd(1,2)
      t10 = t9**2
      t12 = gtd(1,3)
      t16 = t12**2
      detgtd = t1*t2*t4+2*t9*t12*t6-t1*t7-t10*t4-t16*t2
      idetgtd = 0.1D1/detgtd
      t2 = gtd(2,3)**2
      t14 = gtd(1,3)**2
      t22 = gtd(1,2)**2
      gtu(1,1) = idetgtd*(gtd(2,2)*gtd(3,3)-t2)
      gtu(1,2) = idetgtd*(-gtd(1,2)*gtd(3,3)+gtd(1,3)*gtd(2,3))
      gtu(1,3) = idetgtd*(gtd(1,2)*gtd(2,3)-gtd(1,3)*gtd(2,2))
      gtu(2,1) = gtu(1,2)
      gtu(2,2) = idetgtd*(gtd(1,1)*gtd(3,3)-t14)
      gtu(2,3) = idetgtd*(-gtd(1,1)*gtd(2,3)+gtd(1,2)*gtd(1,3))
      gtu(3,1) = gtu(1,3)
      gtu(3,2) = gtu(2,3)
      gtu(3,3) = idetgtd*(gtd(1,1)*gtd(2,2)-t22)
      t1 = Alpha**2
      t2 = abs(chi)
      t3 = 1/t2
      t4 = Betau(1)**2
      t6 = gtd(1,2)*Betau(2)
      t9 = gtd(1,3)*Betau(3)
      t12 = Betau(2)**2
      t14 = gtd(2,3)*Betau(3)
      t17 = Betau(3)**2
      t26 = 0.1D1*t3*(gtd(1,1)*Betau(1)+t6+t9)
      t31 = 0.1D1*t3*(gtd(1,2)*Betau(1)+gtd(2,2)*Betau(2)+t14)
      t37 = 0.1D1*t3*(gtd(1,3)*Betau(1)+gtd(2,3)*Betau(2)+gtd(3,3)*Betau
     #(3))
      t41 = 0.1D1*t3*gtd(1,2)
      t43 = 0.1D1*t3*gtd(1,3)
      t47 = 0.1D1*t3*gtd(2,3)
      g4d(0,0) = -t1+0.1D1*t3*(gtd(2,2)*t12+2*t14*Betau(2)+gtd(3,3)*t17+
     #gtd(1,1)*t4+2*t6*Betau(1)+2*t9*Betau(1))
      g4d(0,1) = t26
      g4d(0,2) = t31
      g4d(0,3) = t37
      g4d(1,0) = t26
      g4d(1,1) = 0.1D1*t3*gtd(1,1)
      g4d(1,2) = t41
      g4d(1,3) = t43
      g4d(2,0) = t31
      g4d(2,1) = t41
      g4d(2,2) = 0.1D1*t3*gtd(2,2)
      g4d(2,3) = t47
      g4d(3,0) = t37
      g4d(3,1) = t43
      g4d(3,2) = t47
      g4d(3,3) = 0.1D1*t3*gtd(3,3)
      t1 = Alpha**2
      t2 = 1/t1
      t3 = t2*Betau(1)
      t4 = t2*Betau(2)
      t5 = t2*Betau(3)
      t7 = Betau(1)**2
      t12 = chi*gtu(1,2)-t3*Betau(2)
      t15 = chi*gtu(1,3)-t3*Betau(3)
      t17 = Betau(2)**2
      t22 = chi*gtu(2,3)-t4*Betau(3)
      t24 = Betau(3)**2
      g4u(0,0) = -t2
      g4u(0,1) = t3
      g4u(0,2) = t4
      g4u(0,3) = t5
      g4u(1,0) = t3
      g4u(1,1) = chi*gtu(1,1)-t2*t7
      g4u(1,2) = t12
      g4u(1,3) = t15
      g4u(2,0) = t4
      g4u(2,1) = t12
      g4u(2,2) = chi*gtu(2,2)-t2*t17
      g4u(2,3) = t22
      g4u(3,0) = t5
      g4u(3,1) = t15
      g4u(3,2) = t22
      g4u(3,3) = chi*gtu(3,3)-t2*t24
      d_phiR4d(0) = -Alpha*piR+Betau(1)*d_phiR(1)+Betau(2)*d_phiR(2)+Bet
     #au(3)*d_phiR(3)
      d_phiR4d(1) = d_phiR(1)
      d_phiR4d(2) = d_phiR(2)
      d_phiR4d(3) = d_phiR(3)
      t5 = -Alpha*piR+Betau(1)*d_phiR(1)+Betau(2)*d_phiR(2)+Betau(3)*d_p
     #hiR(3)
      d_phiR4u(0) = g4u(0,0)*t5+g4u(0,1)*d_phiR(1)+g4u(0,2)*d_phiR(2)+g4
     #u(0,3)*d_phiR(3)
      d_phiR4u(1) = g4u(0,1)*t5+g4u(1,1)*d_phiR(1)+g4u(1,2)*d_phiR(2)+g4
     #u(1,3)*d_phiR(3)
      d_phiR4u(2) = g4u(0,2)*t5+g4u(1,2)*d_phiR(1)+g4u(2,2)*d_phiR(2)+g4
     #u(2,3)*d_phiR(3)
      d_phiR4u(3) = g4u(0,3)*t5+g4u(1,3)*d_phiR(1)+g4u(2,3)*d_phiR(2)+g4
     #u(3,3)*d_phiR(3)
      d_phiI4d(0) = -Alpha*piI+Betau(1)*d_phiI(1)+Betau(2)*d_phiI(2)+Bet
     #au(3)*d_phiI(3)
      d_phiI4d(1) = d_phiI(1)
      d_phiI4d(2) = d_phiI(2)
      d_phiI4d(3) = d_phiI(3)
      t5 = -Alpha*piI+Betau(1)*d_phiI(1)+Betau(2)*d_phiI(2)+Betau(3)*d_p
     #hiI(3)
      d_phiI4u(0) = g4u(0,0)*t5+g4u(0,1)*d_phiI(1)+g4u(0,2)*d_phiI(2)+g4
     #u(0,3)*d_phiI(3)
      d_phiI4u(1) = g4u(0,1)*t5+g4u(1,1)*d_phiI(1)+g4u(1,2)*d_phiI(2)+g4
     #u(1,3)*d_phiI(3)
      d_phiI4u(2) = g4u(0,2)*t5+g4u(1,2)*d_phiI(1)+g4u(2,2)*d_phiI(2)+g4
     #u(2,3)*d_phiI(3)
      d_phiI4u(3) = g4u(0,3)*t5+g4u(1,3)*d_phiI(1)+g4u(2,3)*d_phiI(2)+g4
     #u(3,3)*d_phiI(3)
      dphiR4sq = d_phiR4u(0)*d_phiR4d(0)+d_phiR4u(1)*d_phiR4d(1)+d_phiR4
     #u(2)*d_phiR4d(2)+d_phiR4u(3)*d_phiR4d(3)
      dphiI4sq = d_phiI4u(0)*d_phiI4d(0)+d_phiI4u(1)*d_phiI4d(1)+d_phiI4
     #u(2)*d_phiI4d(2)+d_phiI4u(3)*d_phiI4d(3)
      dphiIR4sq = d_phiI4u(0)*d_phiR4d(0)+d_phiI4u(1)*d_phiR4d(1)+d_phiI
     #4u(2)*d_phiR4d(2)+d_phiI4u(3)*d_phiR4d(3)
      dphiRI4sq = d_phiR4u(0)*d_phiI4d(0)+d_phiR4u(1)*d_phiI4d(1)+d_phiR
     #4u(2)*d_phiI4d(2)+d_phiR4u(3)*d_phiI4d(3)
      t1 = dil_mass**2
      t2 = phiR**2
      sfVR = t1*t2
      t3 = axn_mass**2
      t4 = phiI**2
      sfVI = t3*t4
      dVdphiR2 = t1
      dVdphiI2 = t3
      t1 = d_phiR4u(0)**2
      t2 = dphiR4sq+sfVR
      t7 = exp(4*axn_alpha*phiR)
      t8 = d_phiI4u(0)**2
      t9 = dphiI4sq+sfVI
      t25 = d_phiR4u(0)*d_phiR4u(1)-0.5D0*g4u(0,1)*t2+0.25D0*t7*(d_phiI4
     #u(0)*d_phiI4u(1)-0.5D0*g4u(0,1)*t9)
      t35 = d_phiR4u(0)*d_phiR4u(2)-0.5D0*g4u(0,2)*t2+0.25D0*t7*(d_phiI4
     #u(0)*d_phiI4u(2)-0.5D0*g4u(0,2)*t9)
      t45 = d_phiR4u(0)*d_phiR4u(3)-0.5D0*g4u(0,3)*t2+0.25D0*t7*(d_phiI4
     #u(0)*d_phiI4u(3)-0.5D0*g4u(0,3)*t9)
      t46 = d_phiR4u(1)**2
      t49 = d_phiI4u(1)**2
      t65 = d_phiR4u(1)*d_phiR4u(2)-0.5D0*g4u(1,2)*t2+0.25D0*t7*(d_phiI4
     #u(1)*d_phiI4u(2)-0.5D0*g4u(1,2)*t9)
      t75 = d_phiR4u(1)*d_phiR4u(3)-0.5D0*g4u(1,3)*t2+0.25D0*t7*(d_phiI4
     #u(1)*d_phiI4u(3)-0.5D0*g4u(1,3)*t9)
      t76 = d_phiR4u(2)**2
      t79 = d_phiI4u(2)**2
      t95 = d_phiR4u(2)*d_phiR4u(3)-0.5D0*g4u(2,3)*t2+0.25D0*t7*(d_phiI4
     #u(2)*d_phiI4u(3)-0.5D0*g4u(2,3)*t9)
      t96 = d_phiR4u(3)**2
      t99 = d_phiI4u(3)**2
      Tsfu(0,0) = t1-0.5D0*g4u(0,0)*t2+0.25D0*t7*(t8-0.5D0*g4u(0,0)*t9)
      Tsfu(0,1) = t25
      Tsfu(0,2) = t35
      Tsfu(0,3) = t45
      Tsfu(1,0) = t25
      Tsfu(1,1) = t46-0.5D0*g4u(1,1)*t2+0.25D0*t7*(t49-0.5D0*g4u(1,1)*t9
     #)
      Tsfu(1,2) = t65
      Tsfu(1,3) = t75
      Tsfu(2,0) = t35
      Tsfu(2,1) = t65
      Tsfu(2,2) = t76-0.5D0*g4u(2,2)*t2+0.25D0*t7*(t79-0.5D0*g4u(2,2)*t9
     #)
      Tsfu(2,3) = t95
      Tsfu(3,0) = t45
      Tsfu(3,1) = t75
      Tsfu(3,2) = t95
      Tsfu(3,3) = t96-0.5D0*g4u(3,3)*t2+0.25D0*t7*(t99-0.5D0*g4u(3,3)*t9
     #)
      t5 = Alpha*(gtd(1,1)*emEu(1)+gtd(1,2)*emEu(2)+gtd(1,3)*emEu(3))
      t10 = Alpha*(gtd(1,2)*emEu(1)+gtd(2,2)*emEu(2)+gtd(2,3)*emEu(3))
      t15 = Alpha*(gtd(1,3)*emEu(1)+gtd(2,3)*emEu(2)+gtd(3,3)*emEu(3))
      t16 = sqrt(inv_chi)
      t17 = t16**2
      t18 = t17*t16
      t19 = t18*emBu(3)
      t20 = t18*emBu(2)
      t21 = t18*emBu(1)
      Femd(0,0) = 0
      Femd(0,1) = -t5
      Femd(0,2) = -t10
      Femd(0,3) = -t15
      Femd(1,0) = t5
      Femd(1,1) = 0
      Femd(1,2) = t19
      Femd(1,3) = -t20
      Femd(2,0) = t10
      Femd(2,1) = -t19
      Femd(2,2) = 0
      Femd(2,3) = t21
      Femd(3,0) = t15
      Femd(3,1) = t20
      Femd(3,2) = -t21
      Femd(3,3) = 0
      t2 = g4u(0,1)*Femd(0,1)
      t3 = g4u(0,2)*Femd(0,2)
      t4 = g4u(0,3)*Femd(0,3)
      t27 = g4u(1,2)*Femd(1,2)
      t28 = g4u(1,3)*Femd(1,3)
      t51 = g4u(2,3)*Femd(2,3)
      Femud(0,0) = g4u(0,0)*Femd(0,0)-t2-t3-t4
      Femud(0,1) = g4u(0,0)*Femd(0,1)+g4u(0,1)*Femd(1,1)-g4u(0,2)*Femd(1
     #,2)-g4u(0,3)*Femd(1,3)
      Femud(0,2) = g4u(0,0)*Femd(0,2)+g4u(0,1)*Femd(1,2)+g4u(0,2)*Femd(2
     #,2)-g4u(0,3)*Femd(2,3)
      Femud(0,3) = g4u(0,0)*Femd(0,3)+g4u(0,1)*Femd(1,3)+g4u(0,2)*Femd(2
     #,3)+g4u(0,3)*Femd(3,3)
      Femud(1,0) = g4u(0,1)*Femd(0,0)-g4u(1,1)*Femd(0,1)-g4u(1,2)*Femd(0
     #,2)-g4u(1,3)*Femd(0,3)
      Femud(1,1) = g4u(1,1)*Femd(1,1)+t2-t27-t28
      Femud(1,2) = g4u(0,1)*Femd(0,2)+g4u(1,1)*Femd(1,2)+g4u(1,2)*Femd(2
     #,2)-g4u(1,3)*Femd(2,3)
      Femud(1,3) = g4u(0,1)*Femd(0,3)+g4u(1,1)*Femd(1,3)+g4u(1,2)*Femd(2
     #,3)+g4u(1,3)*Femd(3,3)
      Femud(2,0) = g4u(0,2)*Femd(0,0)-g4u(1,2)*Femd(0,1)-g4u(2,2)*Femd(0
     #,2)-g4u(2,3)*Femd(0,3)
      Femud(2,1) = g4u(0,2)*Femd(0,1)+g4u(1,2)*Femd(1,1)-g4u(2,2)*Femd(1
     #,2)-g4u(2,3)*Femd(1,3)
      Femud(2,2) = g4u(2,2)*Femd(2,2)+t27+t3-t51
      Femud(2,3) = g4u(0,2)*Femd(0,3)+g4u(1,2)*Femd(1,3)+g4u(2,2)*Femd(2
     #,3)+g4u(2,3)*Femd(3,3)
      Femud(3,0) = g4u(0,3)*Femd(0,0)-g4u(1,3)*Femd(0,1)-g4u(2,3)*Femd(0
     #,2)-g4u(3,3)*Femd(0,3)
      Femud(3,1) = g4u(0,3)*Femd(0,1)+g4u(1,3)*Femd(1,1)-g4u(2,3)*Femd(1
     #,2)-g4u(3,3)*Femd(1,3)
      Femud(3,2) = g4u(0,3)*Femd(0,2)+g4u(1,3)*Femd(1,2)+g4u(2,3)*Femd(2
     #,2)-g4u(3,3)*Femd(2,3)
      Femud(3,3) = g4u(3,3)*Femd(3,3)+t28+t4+t51
      t10 = g4u(0,1)*Femud(0,0)+g4u(1,1)*Femud(0,1)+g4u(1,2)*Femud(0,2)+
     #g4u(1,3)*Femud(0,3)
      t15 = g4u(0,2)*Femud(0,0)+g4u(1,2)*Femud(0,1)+g4u(2,2)*Femud(0,2)+
     #g4u(2,3)*Femud(0,3)
      t20 = g4u(0,3)*Femud(0,0)+g4u(1,3)*Femud(0,1)+g4u(2,3)*Femud(0,2)+
     #g4u(3,3)*Femud(0,3)
      t30 = g4u(0,2)*Femud(1,0)+g4u(1,2)*Femud(1,1)+g4u(2,2)*Femud(1,2)+
     #g4u(2,3)*Femud(1,3)
      t35 = g4u(0,3)*Femud(1,0)+g4u(1,3)*Femud(1,1)+g4u(2,3)*Femud(1,2)+
     #g4u(3,3)*Femud(1,3)
      t45 = g4u(0,3)*Femud(2,0)+g4u(1,3)*Femud(2,1)+g4u(2,3)*Femud(2,2)+
     #g4u(3,3)*Femud(2,3)
      Femu(0,0) = g4u(0,0)*Femud(0,0)+g4u(0,1)*Femud(0,1)+g4u(0,2)*Femud
     #(0,2)+g4u(0,3)*Femud(0,3)
      Femu(0,1) = t10
      Femu(0,2) = t15
      Femu(0,3) = t20
      Femu(1,0) = -t10
      Femu(1,1) = g4u(0,1)*Femud(1,0)+g4u(1,1)*Femud(1,1)+g4u(1,2)*Femud
     #(1,2)+g4u(1,3)*Femud(1,3)
      Femu(1,2) = t30
      Femu(1,3) = t35
      Femu(2,0) = -t15
      Femu(2,1) = -t30
      Femu(2,2) = g4u(0,2)*Femud(2,0)+g4u(1,2)*Femud(2,1)+g4u(2,2)*Femud
     #(2,2)+g4u(2,3)*Femud(2,3)
      Femu(2,3) = t45
      Femu(3,0) = -t20
      Femu(3,1) = -t35
      Femu(3,2) = -t45
      Femu(3,3) = g4u(0,3)*Femud(3,0)+g4u(1,3)*Femud(3,1)+g4u(2,3)*Femud
     #(3,2)+g4u(3,3)*Femud(3,3)
      t2 = chi**(-0.3D1)
      t3 = sqrt(t2)
      t5 = 1/Alpha/t3
      t7 = 0.1D1*t5*Femd(2,3)
      t9 = 0.1D1*t5*Femd(1,3)
      t11 = 0.1D1*t5*Femd(1,2)
      t13 = 0.1D1*t5*Femd(0,3)
      t15 = 0.1D1*t5*Femd(0,2)
      t17 = 0.1D1*t5*Femd(0,1)
      dualFemu(0,0) = 0.D0
      dualFemu(0,1) = t7
      dualFemu(0,2) = -t9
      dualFemu(0,3) = t11
      dualFemu(1,0) = -t7
      dualFemu(1,1) = 0.D0
      dualFemu(1,2) = t13
      dualFemu(1,3) = -t15
      dualFemu(2,0) = t9
      dualFemu(2,1) = -t13
      dualFemu(2,2) = 0.D0
      dualFemu(2,3) = t17
      dualFemu(3,0) = -t11
      dualFemu(3,1) = t15
      dualFemu(3,2) = -t17
      dualFemu(3,3) = 0.D0
      Fsq = Femu(0,0)*Femd(0,0)+2*Femu(0,1)*Femd(0,1)+2*Femu(0,2)*Femd(0
     #,2)+2*Femu(0,3)*Femd(0,3)+Femu(1,1)*Femd(1,1)+2*Femu(1,2)*Femd(1,2
     #)+2*Femu(1,3)*Femd(1,3)+Femu(2,2)*Femd(2,2)+2*Femu(2,3)*Femd(2,3)+
     #Femu(3,3)*Femd(3,3)
      t14 = -Femud(0,0)*Femu(0,1)+Femud(0,1)*Femu(1,1)+Femud(0,2)*Femu(1
     #,2)+Femud(0,3)*Femu(1,3)-0.25D0*g4u(0,1)*Fsq
      t21 = -Femud(0,0)*Femu(0,2)-Femud(0,1)*Femu(1,2)+Femud(0,2)*Femu(2
     #,2)+Femud(0,3)*Femu(2,3)-0.25D0*g4u(0,2)*Fsq
      t28 = -Femud(0,0)*Femu(0,3)-Femud(0,1)*Femu(1,3)-Femud(0,2)*Femu(2
     #,3)+Femud(0,3)*Femu(3,3)-0.25D0*g4u(0,3)*Fsq
      t42 = -Femud(1,0)*Femu(0,2)-Femud(1,1)*Femu(1,2)+Femud(1,2)*Femu(2
     #,2)+Femud(1,3)*Femu(2,3)-0.25D0*g4u(1,2)*Fsq
      t49 = -Femud(1,0)*Femu(0,3)-Femud(1,1)*Femu(1,3)-Femud(1,2)*Femu(2
     #,3)+Femud(1,3)*Femu(3,3)-0.25D0*g4u(1,3)*Fsq
      t63 = -Femud(2,0)*Femu(0,3)-Femud(2,1)*Femu(1,3)-Femud(2,2)*Femu(2
     #,3)+Femud(2,3)*Femu(3,3)-0.25D0*g4u(2,3)*Fsq
      Temu(0,0) = Femud(0,0)*Femu(0,0)+Femud(0,1)*Femu(0,1)+Femud(0,2)*F
     #emu(0,2)+Femud(0,3)*Femu(0,3)-0.25D0*g4u(0,0)*Fsq
      Temu(0,1) = t14
      Temu(0,2) = t21
      Temu(0,3) = t28
      Temu(1,0) = t14
      Temu(1,1) = -Femud(1,0)*Femu(0,1)+Femud(1,1)*Femu(1,1)+Femud(1,2)*
     #Femu(1,2)+Femud(1,3)*Femu(1,3)-0.25D0*g4u(1,1)*Fsq
      Temu(1,2) = t42
      Temu(1,3) = t49
      Temu(2,0) = t21
      Temu(2,1) = t42
      Temu(2,2) = -Femud(2,0)*Femu(0,2)-Femud(2,1)*Femu(1,2)+Femud(2,2)*
     #Femu(2,2)+Femud(2,3)*Femu(2,3)-0.25D0*g4u(2,2)*Fsq
      Temu(2,3) = t63
      Temu(3,0) = t28
      Temu(3,1) = t49
      Temu(3,2) = t63
      Temu(3,3) = -Femud(3,0)*Femu(0,3)-Femud(3,1)*Femu(1,3)-Femud(3,2)*
     #Femu(2,3)+Femud(3,3)*Femu(3,3)-0.25D0*g4u(3,3)*Fsq
      t3 = exp(-0.2D1*dil_alpha*phiR)
      t7 = 1/four_pi_G
      t12 = (Tsfu(0,1)+0.2D1*t3*Temu(0,1))*t7
      t16 = (Tsfu(0,2)+0.2D1*t3*Temu(0,2))*t7
      t20 = (Tsfu(0,3)+0.2D1*t3*Temu(0,3))*t7
      t28 = (Tsfu(1,2)+0.2D1*t3*Temu(1,2))*t7
      t32 = (Tsfu(1,3)+0.2D1*t3*Temu(1,3))*t7
      t40 = (Tsfu(2,3)+0.2D1*t3*Temu(2,3))*t7
      Tu(0,0) = (Tsfu(0,0)+0.2D1*t3*Temu(0,0))*t7
      Tu(0,1) = t12
      Tu(0,2) = t16
      Tu(0,3) = t20
      Tu(1,0) = t12
      Tu(1,1) = (Tsfu(1,1)+0.2D1*t3*Temu(1,1))*t7
      Tu(1,2) = t28
      Tu(1,3) = t32
      Tu(2,0) = t16
      Tu(2,1) = t28
      Tu(2,2) = (Tsfu(2,2)+0.2D1*t3*Temu(2,2))*t7
      Tu(2,3) = t40
      Tu(3,0) = t20
      Tu(3,1) = t32
      Tu(3,2) = t40
      Tu(3,3) = (Tsfu(3,3)+0.2D1*t3*Temu(3,3))*t7
      t1 = Alpha**2
      rho_ADM = t1*Tu(0,0)
      t2 = Betau(1)*Tu(0,0)+Tu(0,1)
      t5 = Betau(2)*Tu(0,0)+Tu(0,2)
      t8 = Betau(3)*Tu(0,0)+Tu(0,3)
      Jtd_ADM(1) = Alpha*(t2*gtd(1,1)+t5*gtd(1,2)+t8*gtd(1,3))
      Jtd_ADM(2) = Alpha*(t2*gtd(1,2)+t5*gtd(2,2)+t8*gtd(2,3))
      Jtd_ADM(3) = Alpha*(t2*gtd(1,3)+t5*gtd(2,3)+t8*gtd(3,3))
      t4 = gtd(1,1)*Betau(1)+gtd(1,2)*Betau(2)+gtd(1,3)*Betau(3)
      t5 = t4**2
      t10 = gtd(1,1)*Tu(0,1)+gtd(1,2)*Tu(0,2)+gtd(1,3)*Tu(0,3)
      t13 = gtd(1,1)**2
      t15 = gtd(1,2)**2
      t17 = gtd(1,3)**2
      t19 = gtd(1,1)*gtd(1,2)
      t22 = gtd(1,1)*gtd(1,3)
      t25 = gtd(1,2)*gtd(1,3)
      t32 = gtd(1,2)*Betau(1)+gtd(2,2)*Betau(2)+gtd(2,3)*Betau(3)
      t38 = gtd(1,2)*Tu(0,1)+gtd(2,2)*Tu(0,2)+gtd(2,3)*Tu(0,3)
      t45 = gtd(1,1)*gtd(2,3)
      t48 = gtd(1,2)*gtd(2,2)
      t50 = gtd(1,2)*gtd(2,3)
      t52 = gtd(1,3)*gtd(2,2)
      t54 = gtd(1,3)*gtd(2,3)
      t56 = t4*t32*Tu(0,0)+gtd(1,1)*gtd(2,2)*Tu(1,2)+t32*t10+t15*Tu(1,2)
     #+t19*Tu(1,1)+t25*Tu(1,3)+t4*t38+t45*Tu(1,3)+t48*Tu(2,2)+t50*Tu(2,3
     #)+t52*Tu(2,3)+t54*Tu(3,3)
      t60 = gtd(1,3)*Betau(1)+gtd(2,3)*Betau(2)+gtd(3,3)*Betau(3)
      t66 = gtd(1,3)*Tu(0,1)+gtd(2,3)*Tu(0,2)+gtd(3,3)*Tu(0,3)
      t76 = gtd(1,2)*gtd(3,3)
      t79 = gtd(1,3)*gtd(3,3)
      t81 = t4*t60*Tu(0,0)+gtd(1,1)*gtd(3,3)*Tu(1,3)+t60*t10+t17*Tu(1,3)
     #+t22*Tu(1,1)+t25*Tu(1,2)+t4*t66+t45*Tu(1,2)+t50*Tu(2,2)+t54*Tu(2,3
     #)+t76*Tu(2,3)+t79*Tu(3,3)
      t82 = t32**2
      t87 = gtd(2,2)**2
      t89 = gtd(2,3)**2
      t95 = gtd(2,2)*gtd(2,3)
      t112 = gtd(2,3)*gtd(3,3)
      t114 = t32*t60*Tu(0,0)+gtd(2,2)*gtd(3,3)*Tu(2,3)+t112*Tu(3,3)+t25*
     #Tu(1,1)+t32*t66+t60*t38+t50*Tu(1,2)+t52*Tu(1,2)+t54*Tu(1,3)+t76*Tu
     #(1,3)+t89*Tu(2,3)+t95*Tu(2,2)
      t115 = t60**2
      t121 = gtd(3,3)**2
      pTtd_ADM(1,1) = t5*Tu(0,0)+0.2D1*t4*t10+t13*Tu(1,1)+t15*Tu(2,2)+t1
     #7*Tu(3,3)+0.2D1*t19*Tu(1,2)+0.2D1*t22*Tu(1,3)+0.2D1*t25*Tu(2,3)
      pTtd_ADM(1,2) = t56
      pTtd_ADM(1,3) = t81
      pTtd_ADM(2,1) = pTtd_ADM(1,2)
      pTtd_ADM(2,2) = t82*Tu(0,0)+0.2D1*t32*t38+t15*Tu(1,1)+t87*Tu(2,2)+
     #t89*Tu(3,3)+0.2D1*t48*Tu(1,2)+0.2D1*t50*Tu(1,3)+0.2D1*t95*Tu(2,3)
      pTtd_ADM(2,3) = t114
      pTtd_ADM(3,1) = pTtd_ADM(1,3)
      pTtd_ADM(3,2) = pTtd_ADM(2,3)
      pTtd_ADM(3,3) = t115*Tu(0,0)+0.2D1*t60*t66+t17*Tu(1,1)+t89*Tu(2,2)
     #+t121*Tu(3,3)+0.2D1*t54*Tu(1,2)+0.2D1*t79*Tu(1,3)+0.2D1*t112*Tu(2,
     #3)
      t2 = gtu(1,2)*Atd(1,2)
      t3 = gtu(1,3)*Atd(1,3)
      t18 = gtu(2,3)*Atd(2,3)
      Atud(1,1) = gtu(1,1)*Atd(1,1)+t2+t3
      Atud(1,2) = gtu(1,1)*Atd(1,2)+gtu(1,2)*Atd(2,2)+gtu(1,3)*Atd(2,3)
      Atud(1,3) = gtu(1,1)*Atd(1,3)+gtu(1,2)*Atd(2,3)+gtu(1,3)*Atd(3,3)
      Atud(2,1) = gtu(1,2)*Atd(1,1)+gtu(2,2)*Atd(1,2)+gtu(2,3)*Atd(1,3)
      Atud(2,2) = gtu(2,2)*Atd(2,2)+t18+t2
      Atud(2,3) = gtu(1,2)*Atd(1,3)+gtu(2,2)*Atd(2,3)+gtu(2,3)*Atd(3,3)
      Atud(3,1) = gtu(1,3)*Atd(1,1)+gtu(2,3)*Atd(1,2)+gtu(3,3)*Atd(1,3)
      Atud(3,2) = gtu(1,3)*Atd(1,2)+gtu(2,3)*Atd(2,2)+gtu(3,3)*Atd(2,3)
      Atud(3,3) = gtu(3,3)*Atd(3,3)+t18+t3
      Atu(1,1) = Atud(1,1)*gtu(1,1)+Atud(1,2)*gtu(1,2)+Atud(1,3)*gtu(1,3
     #)
      Atu(1,2) = Atud(1,1)*gtu(1,2)+Atud(1,2)*gtu(2,2)+Atud(1,3)*gtu(2,3
     #)
      Atu(1,3) = Atud(1,1)*gtu(1,3)+Atud(1,2)*gtu(2,3)+Atud(1,3)*gtu(3,3
     #)
      Atu(2,1) = Atu(1,2)
      Atu(2,2) = Atud(2,1)*gtu(1,2)+Atud(2,2)*gtu(2,2)+Atud(2,3)*gtu(2,3
     #)
      Atu(2,3) = Atud(2,1)*gtu(1,3)+Atud(2,2)*gtu(2,3)+Atud(2,3)*gtu(3,3
     #)
      Atu(3,1) = Atu(1,3)
      Atu(3,2) = Atu(2,3)
      Atu(3,3) = Atud(3,1)*gtu(1,3)+Atud(3,2)*gtu(2,3)+Atud(3,3)*gtu(3,3
     #)
      Ctd(1,1,1) = d_gtd(1,1,1)
      Ctd(1,1,2) = d_gtd(2,1,1)
      Ctd(1,1,3) = d_gtd(3,1,1)
      Ctd(1,2,1) = Ctd(1,1,2)
      Ctd(1,2,2) = 2*d_gtd(2,1,2)-d_gtd(1,2,2)
      Ctd(1,2,3) = d_gtd(2,1,3)+d_gtd(3,1,2)-d_gtd(1,2,3)
      Ctd(1,3,1) = Ctd(1,1,3)
      Ctd(1,3,2) = Ctd(1,2,3)
      Ctd(1,3,3) = 2*d_gtd(3,1,3)-d_gtd(1,3,3)
      Ctd(2,1,1) = 2*d_gtd(1,1,2)-d_gtd(2,1,1)
      Ctd(2,1,2) = d_gtd(1,2,2)
      Ctd(2,1,3) = d_gtd(1,2,3)+d_gtd(3,1,2)-d_gtd(2,1,3)
      Ctd(2,2,1) = Ctd(2,1,2)
      Ctd(2,2,2) = d_gtd(2,2,2)
      Ctd(2,2,3) = d_gtd(3,2,2)
      Ctd(2,3,1) = Ctd(2,1,3)
      Ctd(2,3,2) = Ctd(2,2,3)
      Ctd(2,3,3) = 2*d_gtd(3,2,3)-d_gtd(2,3,3)
      Ctd(3,1,1) = 2*d_gtd(1,1,3)-d_gtd(3,1,1)
      Ctd(3,1,2) = d_gtd(1,2,3)+d_gtd(2,1,3)-d_gtd(3,1,2)
      Ctd(3,1,3) = d_gtd(1,3,3)
      Ctd(3,2,1) = Ctd(3,1,2)
      Ctd(3,2,2) = 2*d_gtd(2,2,3)-d_gtd(3,2,2)
      Ctd(3,2,3) = d_gtd(2,3,3)
      Ctd(3,3,1) = Ctd(3,1,3)
      Ctd(3,3,2) = Ctd(3,2,3)
      Ctd(3,3,3) = d_gtd(3,3,3)
      Ct(1,1,1) = gtu(1,1)*Ctd(1,1,1)+gtu(1,2)*Ctd(2,1,1)+gtu(1,3)*Ctd(3
     #,1,1)
      Ct(1,1,2) = gtu(1,1)*Ctd(1,1,2)+gtu(1,2)*Ctd(2,1,2)+gtu(1,3)*Ctd(3
     #,1,2)
      Ct(1,1,3) = gtu(1,1)*Ctd(1,1,3)+gtu(1,2)*Ctd(2,1,3)+gtu(1,3)*Ctd(3
     #,1,3)
      Ct(1,2,1) = Ct(1,1,2)
      Ct(1,2,2) = gtu(1,1)*Ctd(1,2,2)+gtu(1,2)*Ctd(2,2,2)+gtu(1,3)*Ctd(3
     #,2,2)
      Ct(1,2,3) = gtu(1,1)*Ctd(1,2,3)+gtu(1,2)*Ctd(2,2,3)+gtu(1,3)*Ctd(3
     #,2,3)
      Ct(1,3,1) = Ct(1,1,3)
      Ct(1,3,2) = Ct(1,2,3)
      Ct(1,3,3) = gtu(1,1)*Ctd(1,3,3)+gtu(1,2)*Ctd(2,3,3)+gtu(1,3)*Ctd(3
     #,3,3)
      Ct(2,1,1) = gtu(1,2)*Ctd(1,1,1)+gtu(2,2)*Ctd(2,1,1)+gtu(2,3)*Ctd(3
     #,1,1)
      Ct(2,1,2) = gtu(1,2)*Ctd(1,1,2)+gtu(2,2)*Ctd(2,1,2)+gtu(2,3)*Ctd(3
     #,1,2)
      Ct(2,1,3) = gtu(1,2)*Ctd(1,1,3)+gtu(2,2)*Ctd(2,1,3)+gtu(2,3)*Ctd(3
     #,1,3)
      Ct(2,2,1) = Ct(2,1,2)
      Ct(2,2,2) = gtu(1,2)*Ctd(1,2,2)+gtu(2,2)*Ctd(2,2,2)+gtu(2,3)*Ctd(3
     #,2,2)
      Ct(2,2,3) = gtu(1,2)*Ctd(1,2,3)+gtu(2,2)*Ctd(2,2,3)+gtu(2,3)*Ctd(3
     #,2,3)
      Ct(2,3,1) = Ct(2,1,3)
      Ct(2,3,2) = Ct(2,2,3)
      Ct(2,3,3) = gtu(1,2)*Ctd(1,3,3)+gtu(2,2)*Ctd(2,3,3)+gtu(2,3)*Ctd(3
     #,3,3)
      Ct(3,1,1) = gtu(1,3)*Ctd(1,1,1)+gtu(2,3)*Ctd(2,1,1)+gtu(3,3)*Ctd(3
     #,1,1)
      Ct(3,1,2) = gtu(1,3)*Ctd(1,1,2)+gtu(2,3)*Ctd(2,1,2)+gtu(3,3)*Ctd(3
     #,1,2)
      Ct(3,1,3) = gtu(1,3)*Ctd(1,1,3)+gtu(2,3)*Ctd(2,1,3)+gtu(3,3)*Ctd(3
     #,1,3)
      Ct(3,2,1) = Ct(3,1,2)
      Ct(3,2,2) = gtu(1,3)*Ctd(1,2,2)+gtu(2,3)*Ctd(2,2,2)+gtu(3,3)*Ctd(3
     #,2,2)
      Ct(3,2,3) = gtu(1,3)*Ctd(1,2,3)+gtu(2,3)*Ctd(2,2,3)+gtu(3,3)*Ctd(3
     #,2,3)
      Ct(3,3,1) = Ct(3,1,3)
      Ct(3,3,2) = Ct(3,2,3)
      Ct(3,3,3) = gtu(1,3)*Ctd(1,3,3)+gtu(2,3)*Ctd(2,3,3)+gtu(3,3)*Ctd(3
     #,3,3)
      div_Beta = d_Betau(1,1)+d_Betau(2,2)+d_Betau(3,3)
      d_div_Beta(1) = dd_Betau(1,1,1)+dd_Betau(1,2,2)+dd_Betau(1,3,3)
      d_div_Beta(2) = dd_Betau(1,2,1)+dd_Betau(2,2,2)+dd_Betau(2,3,3)
      d_div_Beta(3) = dd_Betau(1,3,1)+dd_Betau(2,3,2)+dd_Betau(3,3,3)
      CalGamt(1) = 0.5D0*gtu(1,1)*Ct(1,1,1)+0.5D0*gtu(2,2)*Ct(1,2,2)+0.5
     #D0*gtu(3,3)*Ct(1,3,3)+gtu(1,2)*Ct(1,1,2)+gtu(1,3)*Ct(1,1,3)+gtu(2,
     #3)*Ct(1,2,3)
      CalGamt(2) = 0.5D0*gtu(1,1)*Ct(2,1,1)+0.5D0*gtu(2,2)*Ct(2,2,2)+0.5
     #D0*gtu(3,3)*Ct(2,3,3)+gtu(1,2)*Ct(2,1,2)+gtu(1,3)*Ct(2,1,3)+gtu(2,
     #3)*Ct(2,2,3)
      CalGamt(3) = 0.5D0*gtu(1,1)*Ct(3,1,1)+0.5D0*gtu(2,2)*Ct(3,2,2)+0.5
     #D0*gtu(3,3)*Ct(3,3,3)+gtu(1,2)*Ct(3,1,2)+gtu(1,3)*Ct(3,1,3)+gtu(2,
     #3)*Ct(3,2,3)
      t2 = d_chi(1)**2
      t24 = inv_chi*d_chi(3)
      t35 = d_chi(2)**2
      t56 = d_chi(3)**2
      Rpd_1(1,1) = 0.5D0*dd_chi(1,1)-0.25D0*t2*inv_chi-0.25D0*Ct(1,1,1)*
     #d_chi(1)-0.25D0*Ct(2,1,1)*d_chi(2)-0.25D0*Ct(3,1,1)*d_chi(3)
      Rpd_1(1,2) = 0.5D0*dd_chi(1,2)-0.25D0*inv_chi*d_chi(2)*d_chi(1)-0.
     #25D0*Ct(1,1,2)*d_chi(1)-0.25D0*Ct(2,1,2)*d_chi(2)-0.25D0*Ct(3,1,2)
     #*d_chi(3)
      Rpd_1(1,3) = 0.5D0*dd_chi(1,3)-0.25D0*t24*d_chi(1)-0.25D0*Ct(1,1,3
     #)*d_chi(1)-0.25D0*Ct(2,1,3)*d_chi(2)-0.25D0*Ct(3,1,3)*d_chi(3)
      Rpd_1(2,1) = Rpd_1(1,2)
      Rpd_1(2,2) = 0.5D0*dd_chi(2,2)-0.25D0*inv_chi*t35-0.25D0*Ct(1,2,2)
     #*d_chi(1)-0.25D0*Ct(2,2,2)*d_chi(2)-0.25D0*Ct(3,2,2)*d_chi(3)
      Rpd_1(2,3) = 0.5D0*dd_chi(2,3)-0.25D0*t24*d_chi(2)-0.25D0*Ct(1,2,3
     #)*d_chi(1)-0.25D0*Ct(2,2,3)*d_chi(2)-0.25D0*Ct(3,2,3)*d_chi(3)
      Rpd_1(3,1) = Rpd_1(1,3)
      Rpd_1(3,2) = Rpd_1(2,3)
      Rpd_1(3,3) = 0.5D0*dd_chi(3,3)-0.25D0*inv_chi*t56-0.25D0*Ct(1,3,3)
     #*d_chi(1)-0.25D0*Ct(2,3,3)*d_chi(2)-0.25D0*Ct(3,3,3)*d_chi(3)
      t18 = Ct(1,1,1)*Ctd(1,1,2)
      t19 = Ct(1,1,2)*Ctd(1,1,1)
      t21 = Ct(2,1,1)*Ctd(2,1,2)
      t22 = Ct(2,1,2)*Ctd(1,1,2)
      t24 = Ct(3,1,1)*Ctd(3,1,2)
      t25 = Ct(3,1,2)*Ctd(1,1,3)
      t30 = Ct(1,1,1)*Ctd(1,1,3)
      t31 = Ct(1,1,3)*Ctd(1,1,1)
      t33 = Ct(2,1,1)*Ctd(2,1,3)
      t34 = Ct(2,1,3)*Ctd(1,1,2)
      t36 = Ct(3,1,1)*Ctd(3,1,3)
      t37 = Ct(3,1,3)*Ctd(1,1,3)
      t52 = Ct(1,1,2)*Ctd(1,1,2)
      t54 = Ct(2,1,2)*Ctd(1,2,2)
      t56 = Ct(2,1,2)*Ctd(2,1,2)
      t57 = Ct(3,1,2)*Ctd(1,2,3)
      t59 = Ct(3,1,2)*Ctd(3,1,2)
      t63 = Ct(1,1,2)*Ctd(1,1,3)
      t64 = Ct(1,1,3)*Ctd(1,1,2)
      t66 = Ct(2,1,2)*Ctd(2,1,3)
      t67 = Ct(2,1,3)*Ctd(1,2,2)
      t69 = Ct(3,1,2)*Ctd(3,1,3)
      t70 = Ct(3,1,3)*Ctd(1,2,3)
      t87 = Ct(2,1,2)*Ctd(1,2,3)
      t89 = Ct(2,1,3)*Ctd(2,1,2)
      t90 = Ct(3,1,2)*Ctd(1,3,3)
      t92 = Ct(3,1,3)*Ctd(3,1,2)
      t96 = Ct(1,1,3)*Ctd(1,1,3)
      t98 = Ct(2,1,3)*Ctd(1,2,3)
      t100 = Ct(2,1,3)*Ctd(2,1,3)
      t101 = Ct(3,1,3)*Ctd(1,3,3)
      t103 = Ct(3,1,3)*Ctd(3,1,3)
      t125 = 0.25D0*gtu(2,3)*(2*t63+t64+2*t87+t89+2*t90+t92)+0.25D0*gtu(
     #3,3)*(3*t96+2*t98+t100+2*t101+t103)-0.5D0*gtu(1,1)*dd_gtd(1,1,1,1)
     #-0.1D1*gtu(1,2)*dd_gtd(1,2,1,1)-0.1D1*gtu(1,3)*dd_gtd(1,3,1,1)-0.5
     #D0*gtu(2,2)*dd_gtd(2,2,1,1)-0.1D1*gtu(2,3)*dd_gtd(2,3,1,1)-0.5D0*g
     #tu(3,3)*dd_gtd(3,3,1,1)+0.1D1*d_Gamt(1,1)*gtd(1,1)+0.1D1*d_Gamt(1,
     #2)*gtd(1,2)+0.1D1*d_Gamt(1,3)*gtd(1,3)
      t143 = Ct(1,1,2)*Ctd(2,1,1)
      t145 = Ct(2,1,1)*Ctd(2,2,2)
      t148 = Ct(3,1,2)*Ctd(2,1,3)
      t153 = Ct(1,1,1)*Ctd(1,2,3)
      t154 = Ct(1,1,3)*Ctd(2,1,1)
      t155 = Ct(1,2,3)*Ctd(1,1,1)
      t156 = Ct(2,1,1)*Ctd(2,2,3)
      t157 = Ct(2,2,3)*Ctd(1,1,2)
      t158 = Ct(3,1,1)*Ctd(3,2,3)
      t159 = Ct(3,1,3)*Ctd(2,1,3)
      t160 = Ct(3,2,3)*Ctd(1,1,3)
      t170 = Ct(1,1,2)*Ctd(1,2,2)
      t171 = Ct(1,1,2)*Ctd(2,1,2)
      t172 = Ct(1,2,2)*Ctd(1,1,2)
      t173 = Ct(2,1,2)*Ctd(2,2,2)
      t174 = 2*t173
      t176 = Ct(3,1,2)*Ctd(2,2,3)
      t177 = Ct(3,1,2)*Ctd(3,2,2)
      t182 = Ct(1,1,2)*Ctd(1,2,3)
      t183 = Ct(1,1,3)*Ctd(2,1,2)
      t184 = Ct(1,2,3)*Ctd(1,1,2)
      t185 = Ct(2,1,2)*Ctd(2,2,3)
      t186 = Ct(2,1,3)*Ctd(2,2,2)
      t187 = Ct(2,2,3)*Ctd(1,2,2)
      t188 = Ct(3,1,2)*Ctd(3,2,3)
      t189 = Ct(3,1,3)*Ctd(2,2,3)
      t190 = Ct(3,2,3)*Ctd(1,2,3)
      t199 = Ct(1,1,2)*Ctd(2,1,3)
      t201 = Ct(1,2,2)*Ctd(1,1,3)
      t203 = Ct(3,1,2)*Ctd(2,3,3)
      t209 = Ct(1,1,3)*Ctd(1,2,3)
      t210 = Ct(1,1,3)*Ctd(2,1,3)
      t211 = Ct(1,2,3)*Ctd(1,1,3)
      t212 = Ct(2,1,3)*Ctd(2,2,3)
      t214 = Ct(2,2,3)*Ctd(1,2,3)
      t215 = Ct(3,1,3)*Ctd(2,3,3)
      t216 = Ct(3,1,3)*Ctd(3,2,3)
      t217 = Ct(3,2,3)*Ctd(1,3,3)
      s1 = 0.25D0*CalGamt(1)*(Ctd(1,1,2)+Ctd(2,1,1))+0.25D0*CalGamt(2)*(
     #Ctd(1,2,2)+Ctd(2,1,2))+0.25D0*CalGamt(3)*(Ctd(1,2,3)+Ctd(2,1,3))+0
     #.25D0*gtu(1,1)*(Ct(1,1,1)*Ctd(2,1,1)+Ct(3,1,1)*Ctd(2,1,3)+t18+t19+
     #2*t21+t22+t24+t25)+0.25D0*gtu(1,2)*(Ct(1,1,1)*Ctd(1,2,2)+Ct(1,2,2)
     #*Ctd(1,1,1)+Ct(2,2,2)*Ctd(1,1,2)+Ct(3,1,1)*Ctd(3,2,2)+Ct(3,2,2)*Ct
     #d(1,1,3)+t143+t145+t148+t56)+0.25D0*gtu(1,3)*(t153+t154+t155+t156+
     #t89+t157+t158+t159+t160)
      t221 = s1+0.25D0*gtu(1,2)*(Ct(1,1,1)*Ctd(2,1,2)+Ct(3,1,1)*Ctd(2,2,
     #3)+t145+2*t52+t54+t56+t57+t59)+0.25D0*gtu(2,2)*(Ct(2,2,2)*Ctd(1,2,
     #2)+Ct(3,2,2)*Ctd(1,2,3)+t170+t171+t172+t174+t176+t177)+0.25D0*gtu(
     #2,3)*(t182+t183+t184+t185+t186+t187+t188+t189+t190)+0.25D0*gtu(1,3
     #)*(Ct(1,1,1)*Ctd(2,1,3)+Ct(3,1,1)*Ctd(2,3,3)+t156+t63+t64+t87+t89+
     #t90+t92)+0.25D0*gtu(2,3)*(Ct(1,1,3)*Ctd(1,2,2)+Ct(2,2,2)*Ctd(1,2,3
     #)+Ct(3,1,3)*Ctd(3,2,2)+Ct(3,2,2)*Ctd(1,3,3)+t185+t186+t199+t201+t2
     #03)+0.25D0*gtu(3,3)*(t209+t210+t211+2*t212+t214+t215+t216+t217)
      t246 = -0.5D0*gtu(1,1)*dd_gtd(1,1,1,2)-0.1D1*gtu(1,2)*dd_gtd(1,2,1
     #,2)-0.1D1*gtu(1,3)*dd_gtd(1,3,1,2)-0.5D0*gtu(2,2)*dd_gtd(2,2,1,2)-
     #0.1D1*gtu(2,3)*dd_gtd(2,3,1,2)-0.5D0*gtu(3,3)*dd_gtd(3,3,1,2)+0.5D
     #0*d_Gamt(1,1)*gtd(1,2)+0.5D0*d_Gamt(1,2)*gtd(2,2)+0.5D0*d_Gamt(1,3
     #)*gtd(2,3)+0.5D0*d_Gamt(2,1)*gtd(1,1)+0.5D0*d_Gamt(2,2)*gtd(1,2)+0
     #.5D0*d_Gamt(2,3)*gtd(1,3)
      t263 = Ct(1,1,2)*Ctd(3,1,1)
      t264 = Ct(2,1,2)*Ctd(3,1,2)
      t269 = Ct(1,1,3)*Ctd(3,1,1)
      t272 = Ct(2,1,3)*Ctd(3,1,2)
      t274 = Ct(3,1,1)*Ctd(3,3,3)
      t284 = Ct(1,1,2)*Ctd(3,1,2)
      t285 = Ct(2,1,2)*Ctd(3,2,2)
      t290 = Ct(1,1,2)*Ctd(1,3,3)
      t291 = Ct(1,1,3)*Ctd(3,1,2)
      t293 = Ct(2,1,2)*Ctd(2,3,3)
      t294 = Ct(2,1,3)*Ctd(3,2,2)
      t296 = Ct(3,1,2)*Ctd(3,3,3)
      t307 = Ct(1,1,2)*Ctd(3,1,3)
      t308 = Ct(2,1,2)*Ctd(3,2,3)
      t312 = Ct(1,1,3)*Ctd(1,3,3)
      t313 = Ct(1,1,3)*Ctd(3,1,3)
      t314 = Ct(1,3,3)*Ctd(1,1,3)
      t315 = Ct(2,1,3)*Ctd(2,3,3)
      t316 = Ct(2,1,3)*Ctd(3,2,3)
      t318 = Ct(3,1,3)*Ctd(3,3,3)
      t319 = 2*t318
      s1 = 0.25D0*CalGamt(1)*(Ctd(1,1,3)+Ctd(3,1,1))+0.25D0*CalGamt(2)*(
     #Ctd(1,2,3)+Ctd(3,1,2))+0.25D0*CalGamt(3)*(Ctd(1,3,3)+Ctd(3,1,3))+0
     #.25D0*gtu(1,1)*(Ct(1,1,1)*Ctd(3,1,1)+Ct(2,1,1)*Ctd(3,1,2)+t30+t31+
     #t33+t34+2*t36+t37)+0.25D0*gtu(1,2)*(t153+t263+t155+t156+t264+t157+
     #t158+t69+t160)+0.25D0*gtu(1,3)*(Ct(1,1,1)*Ctd(1,3,3)+Ct(1,3,3)*Ctd
     #(1,1,1)+Ct(2,1,1)*Ctd(2,3,3)+Ct(2,3,3)*Ctd(1,1,2)+Ct(3,3,3)*Ctd(1,
     #1,3)+t103+t269+t272+t274)
      t324 = s1+0.25D0*gtu(1,2)*(Ct(1,1,1)*Ctd(3,1,2)+Ct(2,1,1)*Ctd(3,2,
     #2)+t158+t63+t64+t66+t67+t69+t70)+0.25D0*gtu(2,2)*(t182+t284+t184+t
     #185+t285+t187+2*t188+t190)+0.25D0*gtu(2,3)*(Ct(1,3,3)*Ctd(1,1,2)+C
     #t(2,3,3)*Ctd(1,2,2)+Ct(3,3,3)*Ctd(1,2,3)+t216+t290+t291+t293+t294+
     #t296)+0.25D0*gtu(1,3)*(Ct(1,1,1)*Ctd(3,1,3)+Ct(2,1,1)*Ctd(3,2,3)+t
     #100+t101+t103+t274+2*t96+t98)+0.25D0*gtu(2,3)*(t307+t209+t211+t308
     #+t212+t214+t296+t216+t217)+0.25D0*gtu(3,3)*(Ct(2,3,3)*Ctd(1,2,3)+C
     #t(3,3,3)*Ctd(1,3,3)+t312+t313+t314+t315+t316+t319)
      t349 = -0.5D0*gtu(1,1)*dd_gtd(1,1,1,3)-0.1D1*gtu(1,2)*dd_gtd(1,2,1
     #,3)-0.1D1*gtu(1,3)*dd_gtd(1,3,1,3)-0.5D0*gtu(2,2)*dd_gtd(2,2,1,3)-
     #0.1D1*gtu(2,3)*dd_gtd(2,3,1,3)-0.5D0*gtu(3,3)*dd_gtd(3,3,1,3)+0.5D
     #0*d_Gamt(1,1)*gtd(1,3)+0.5D0*d_Gamt(1,2)*gtd(2,3)+0.5D0*d_Gamt(1,3
     #)*gtd(3,3)+0.5D0*d_Gamt(3,1)*gtd(1,1)+0.5D0*d_Gamt(3,2)*gtd(1,2)+0
     #.5D0*d_Gamt(3,3)*gtd(1,3)
      t365 = Ct(2,2,2)*Ctd(2,1,2)
      t372 = Ct(1,2,3)*Ctd(2,1,1)
      t374 = Ct(2,2,3)*Ctd(2,1,2)
      t376 = Ct(3,2,3)*Ctd(2,1,3)
      t398 = Ct(1,2,2)*Ctd(1,2,3)
      t399 = Ct(1,2,3)*Ctd(2,1,2)
      t401 = Ct(2,2,2)*Ctd(2,2,3)
      t402 = Ct(2,2,3)*Ctd(2,2,2)
      t404 = Ct(3,2,2)*Ctd(3,2,3)
      t405 = Ct(3,2,3)*Ctd(2,2,3)
      t428 = Ct(1,2,3)*Ctd(1,2,3)
      t429 = Ct(1,2,3)*Ctd(2,1,3)
      t431 = Ct(2,2,3)*Ctd(2,2,3)
      t433 = Ct(3,2,3)*Ctd(2,3,3)
      t435 = Ct(3,2,3)*Ctd(3,2,3)
      t457 = 0.25D0*gtu(2,3)*(2*Ct(1,2,2)*Ctd(2,1,3)+Ct(1,2,3)*Ctd(1,2,2
     #)+2*Ct(3,2,2)*Ctd(2,3,3)+Ct(3,2,3)*Ctd(3,2,2)+2*t401+t402)+0.25D0*
     #gtu(3,3)*(t428+2*t429+3*t431+2*t433+t435)-0.5D0*gtu(1,1)*dd_gtd(1,
     #1,2,2)-0.1D1*gtu(1,2)*dd_gtd(1,2,2,2)-0.1D1*gtu(1,3)*dd_gtd(1,3,2,
     #2)-0.5D0*gtu(2,2)*dd_gtd(2,2,2,2)-0.1D1*gtu(2,3)*dd_gtd(2,3,2,2)-0
     #.5D0*gtu(3,3)*dd_gtd(3,3,2,2)+0.1D1*d_Gamt(2,1)*gtd(1,2)+0.1D1*d_G
     #amt(2,2)*gtd(2,2)+0.1D1*d_Gamt(2,3)*gtd(2,3)
      t474 = Ct(3,2,2)*Ctd(3,1,3)
      t478 = Ct(1,2,3)*Ctd(3,1,1)
      t480 = Ct(2,2,3)*Ctd(3,1,2)
      t482 = Ct(3,2,3)*Ctd(3,1,3)
      t498 = Ct(1,2,3)*Ctd(3,1,2)
      t501 = Ct(2,2,3)*Ctd(3,2,2)
      t503 = Ct(3,2,2)*Ctd(3,3,3)
      t508 = Ct(2,2,3)*Ctd(2,1,3)
      t518 = Ct(1,2,3)*Ctd(1,3,3)
      t519 = Ct(1,2,3)*Ctd(3,1,3)
      t521 = Ct(2,2,3)*Ctd(2,3,3)
      t522 = Ct(2,2,3)*Ctd(3,2,3)
      t523 = Ct(2,3,3)*Ctd(2,2,3)
      t524 = Ct(3,2,3)*Ctd(3,3,3)
      t525 = 2*t524
      s1 = 0.25D0*CalGamt(1)*(Ctd(2,1,3)+Ctd(3,1,2))+0.25D0*CalGamt(2)*(
     #Ctd(2,2,3)+Ctd(3,2,2))+0.25D0*CalGamt(3)*(Ctd(2,3,3)+Ctd(3,2,3))+0
     #.25D0*gtu(1,1)*(t63+t263+t154+t66+t264+t89+2*t69+t159)+0.25D0*gtu(
     #1,2)*(Ct(1,2,2)*Ctd(3,1,1)+Ct(2,2,2)*Ctd(3,1,2)+t182+t185+t188+t37
     #2+t374+t376+t474)+0.25D0*gtu(1,3)*(Ct(1,3,3)*Ctd(2,1,1)+Ct(2,3,3)*
     #Ctd(2,1,2)+Ct(3,3,3)*Ctd(2,1,3)+t290+t293+t296+t478+t480+t482)
      t530 = s1+0.25D0*gtu(1,2)*(Ct(2,2,2)*Ctd(2,1,3)+t183+t186+t188+t18
     #9+t201+t284+t285+t474)+0.25D0*gtu(2,2)*(Ct(1,2,2)*Ctd(3,1,2)+Ct(2,
     #2,2)*Ctd(3,2,2)+t398+t399+t401+t402+2*t404+t405)+0.25D0*gtu(2,3)*(
     #Ct(1,2,2)*Ctd(1,3,3)+Ct(1,3,3)*Ctd(2,1,2)+Ct(2,2,2)*Ctd(2,3,3)+Ct(
     #2,3,3)*Ctd(2,2,2)+Ct(3,3,3)*Ctd(2,2,3)+t435+t498+t501+t503)+0.25D0
     #*gtu(1,3)*(t307+t210+t211+t308+t212+t508+t296+t215+t482)+0.25D0*gt
     #u(2,3)*(Ct(1,2,2)*Ctd(3,1,3)+Ct(2,2,2)*Ctd(3,2,3)+t428+t429+2*t431
     #+t433+t435+t503)+0.25D0*gtu(3,3)*(Ct(1,3,3)*Ctd(2,1,3)+Ct(3,3,3)*C
     #td(2,3,3)+t518+t519+t521+t522+t523+t525)
      t555 = -0.5D0*gtu(1,1)*dd_gtd(1,1,2,3)-0.1D1*gtu(1,2)*dd_gtd(1,2,2
     #,3)-0.1D1*gtu(1,3)*dd_gtd(1,3,2,3)-0.5D0*gtu(2,2)*dd_gtd(2,2,2,3)-
     #0.1D1*gtu(2,3)*dd_gtd(2,3,2,3)-0.5D0*gtu(3,3)*dd_gtd(3,3,2,3)+0.5D
     #0*d_Gamt(2,1)*gtd(1,3)+0.5D0*d_Gamt(2,2)*gtd(2,3)+0.5D0*d_Gamt(2,3
     #)*gtd(3,3)+0.5D0*d_Gamt(3,1)*gtd(1,2)+0.5D0*d_Gamt(3,2)*gtd(2,2)+0
     #.5D0*d_Gamt(3,3)*gtd(2,3)
      t579 = Ct(3,3,3)*Ctd(3,1,3)
      t600 = Ct(3,3,3)*Ctd(3,2,3)
      t647 = 0.25D0*gtu(2,3)*(Ct(1,3,3)*Ctd(1,2,3)+2*t519+2*t522+t523+t5
     #25+t600)+0.25D0*gtu(3,3)*(Ct(1,3,3)*Ctd(1,3,3)+2*Ct(1,3,3)*Ctd(3,1
     #,3)+Ct(2,3,3)*Ctd(2,3,3)+2*Ct(2,3,3)*Ctd(3,2,3)+3*Ct(3,3,3)*Ctd(3,
     #3,3))-0.5D0*gtu(1,1)*dd_gtd(1,1,3,3)-0.1D1*gtu(1,2)*dd_gtd(1,2,3,3
     #)-0.1D1*gtu(1,3)*dd_gtd(1,3,3,3)-0.5D0*gtu(2,2)*dd_gtd(2,2,3,3)-0.
     #1D1*gtu(2,3)*dd_gtd(2,3,3,3)-0.5D0*gtu(3,3)*dd_gtd(3,3,3,3)+0.1D1*
     #d_Gamt(3,1)*gtd(1,3)+0.1D1*d_Gamt(3,2)*gtd(2,3)+0.1D1*d_Gamt(3,3)*
     #gtd(3,3)
      s1 = 0.5D0*CalGamt(1)*Ctd(1,1,1)+0.5D0*CalGamt(2)*Ctd(1,1,2)+0.5D0
     #*CalGamt(3)*Ctd(1,1,3)+0.25D0*gtu(1,1)*(3*Ct(1,1,1)*Ctd(1,1,1)+2*C
     #t(2,1,1)*Ctd(1,1,2)+Ct(2,1,1)*Ctd(2,1,1)+2*Ct(3,1,1)*Ctd(1,1,3)+Ct
     #(3,1,1)*Ctd(3,1,1))+0.25D0*gtu(1,2)*(t18+2*t19+t21+2*t22+t24+2*t25
     #)
      Rtd(1,1) = s1+0.25D0*gtu(1,3)*(t30+2*t31+t33+2*t34+t36+2*t37)+0.25
     #D0*gtu(1,2)*(2*Ct(2,1,1)*Ctd(1,2,2)+Ct(2,1,2)*Ctd(2,1,1)+2*Ct(3,1,
     #1)*Ctd(1,2,3)+Ct(3,1,2)*Ctd(3,1,1)+2*t18+t19)+0.25D0*gtu(2,2)*(3*t
     #52+2*t54+t56+2*t57+t59)+0.25D0*gtu(2,3)*(t63+2*t64+t66+2*t67+t69+2
     #*t70)+0.25D0*gtu(1,3)*(2*Ct(2,1,1)*Ctd(1,2,3)+Ct(2,1,3)*Ctd(2,1,1)
     #+2*Ct(3,1,1)*Ctd(1,3,3)+Ct(3,1,3)*Ctd(3,1,1)+2*t30+t31)+t125
      Rtd(1,2) = t221+t246
      Rtd(1,3) = t324+t349
      Rtd(2,1) = Rtd(1,2)
      s1 = 0.5D0*CalGamt(1)*Ctd(2,1,2)+0.5D0*CalGamt(2)*Ctd(2,2,2)+0.5D0
     #*CalGamt(3)*Ctd(2,2,3)+0.25D0*gtu(1,1)*(t52+2*t143+3*t56+2*t148+t5
     #9)+0.25D0*gtu(1,2)*(2*Ct(1,2,2)*Ctd(2,1,1)+2*Ct(3,2,2)*Ctd(2,1,3)+
     #t170+t173+t177+2*t365)
      Rtd(2,2) = s1+0.25D0*gtu(1,3)*(t182+2*t372+t185+2*t374+t188+2*t376
     #)+0.25D0*gtu(1,2)*(Ct(3,2,2)*Ctd(3,1,2)+2*t171+t172+t174+2*t176+t3
     #65)+0.25D0*gtu(2,2)*(Ct(1,2,2)*Ctd(1,2,2)+2*Ct(1,2,2)*Ctd(2,1,2)+3
     #*Ct(2,2,2)*Ctd(2,2,2)+2*Ct(3,2,2)*Ctd(2,2,3)+Ct(3,2,2)*Ctd(3,2,2))
     #+0.25D0*gtu(2,3)*(t398+2*t399+t401+2*t402+t404+2*t405)+0.25D0*gtu(
     #1,3)*(Ct(3,2,3)*Ctd(3,1,2)+t184+2*t185+2*t199+2*t203+t374)+t457
      Rtd(2,3) = t530+t555
      Rtd(3,1) = Rtd(1,3)
      Rtd(3,2) = Rtd(2,3)
      Rtd(3,3) = 0.5D0*CalGamt(1)*Ctd(3,1,3)+0.5D0*CalGamt(2)*Ctd(3,2,3)
     #+0.5D0*CalGamt(3)*Ctd(3,3,3)+0.25D0*gtu(1,1)*(t96+2*t269+t100+2*t2
     #72+3*t103)+0.25D0*gtu(1,2)*(t209+2*t478+t212+2*t480+t216+2*t482)+0
     #.25D0*gtu(1,3)*(2*Ct(1,3,3)*Ctd(3,1,1)+2*Ct(2,3,3)*Ctd(3,1,2)+t312
     #+t315+t318+2*t579)+0.25D0*gtu(1,2)*(2*t291+t211+2*t294+t508+2*t216
     #+t482)+0.25D0*gtu(2,2)*(t428+2*t498+t431+2*t501+3*t435)+0.25D0*gtu
     #(2,3)*(2*Ct(1,3,3)*Ctd(3,1,2)+2*Ct(2,3,3)*Ctd(3,2,2)+t518+t521+t52
     #4+2*t600)+0.25D0*gtu(1,3)*(Ct(2,3,3)*Ctd(2,1,3)+2*t313+t314+2*t316
     #+t319+t579)+t647
      t10 = inv_chi*eight_pi_G
      Psi1(1,1) = chi*(Alpha*Rtd(1,1)-dd_Alpha(1,1)+0.5D0*Ct(1,1,1)*d_Al
     #pha(1)+0.5D0*Ct(2,1,1)*d_Alpha(2)+0.5D0*Ct(3,1,1)*d_Alpha(3))-t10*
     #Alpha*pTtd_ADM(1,1)+Alpha*Rpd_1(1,1)-d_Alpha(1)*d_chi(1)
      Psi1(1,2) = chi*(Alpha*Rtd(1,2)-dd_Alpha(1,2)+0.5D0*Ct(1,1,2)*d_Al
     #pha(1)+0.5D0*Ct(2,1,2)*d_Alpha(2)+0.5D0*Ct(3,1,2)*d_Alpha(3))-t10*
     #Alpha*pTtd_ADM(1,2)+Alpha*Rpd_1(1,2)-0.5D0*d_Alpha(1)*d_chi(2)-0.5
     #D0*d_Alpha(2)*d_chi(1)
      Psi1(1,3) = chi*(Alpha*Rtd(1,3)-dd_Alpha(1,3)+0.5D0*Ct(1,1,3)*d_Al
     #pha(1)+0.5D0*Ct(2,1,3)*d_Alpha(2)+0.5D0*Ct(3,1,3)*d_Alpha(3))-t10*
     #Alpha*pTtd_ADM(1,3)+Alpha*Rpd_1(1,3)-0.5D0*d_Alpha(1)*d_chi(3)-0.5
     #D0*d_Alpha(3)*d_chi(1)
      Psi1(2,1) = Psi1(1,2)
      Psi1(2,2) = chi*(Alpha*Rtd(2,2)-dd_Alpha(2,2)+0.5D0*Ct(1,2,2)*d_Al
     #pha(1)+0.5D0*Ct(2,2,2)*d_Alpha(2)+0.5D0*Ct(3,2,2)*d_Alpha(3))-t10*
     #Alpha*pTtd_ADM(2,2)+Alpha*Rpd_1(2,2)-d_Alpha(2)*d_chi(2)
      Psi1(2,3) = chi*(Alpha*Rtd(2,3)-dd_Alpha(2,3)+0.5D0*Ct(1,2,3)*d_Al
     #pha(1)+0.5D0*Ct(2,2,3)*d_Alpha(2)+0.5D0*Ct(3,2,3)*d_Alpha(3))-t10*
     #Alpha*pTtd_ADM(2,3)+Alpha*Rpd_1(2,3)-0.5D0*d_Alpha(2)*d_chi(3)-0.5
     #D0*d_Alpha(3)*d_chi(2)
      Psi1(3,1) = Psi1(1,3)
      Psi1(3,2) = Psi1(2,3)
      Psi1(3,3) = chi*(Alpha*Rtd(3,3)-dd_Alpha(3,3)+0.5D0*Ct(1,3,3)*d_Al
     #pha(1)+0.5D0*Ct(2,3,3)*d_Alpha(2)+0.5D0*Ct(3,3,3)*d_Alpha(3))-t10*
     #Alpha*pTtd_ADM(3,3)+Alpha*Rpd_1(3,3)-d_Alpha(3)*d_chi(3)
      third_trPsi1 = 0.3333333333333333D0*Psi1(1,1)*gtu(1,1)+0.333333333
     #3333333D0*Psi1(2,2)*gtu(2,2)+0.3333333333333333D0*Psi1(3,3)*gtu(3,
     #3)+0.6666666666666666D0*Psi1(1,2)*gtu(1,2)+0.6666666666666666D0*Ps
     #i1(1,3)*gtu(1,3)+0.6666666666666666D0*Psi1(2,3)*gtu(2,3)
      Psi1TF(1,1) = -third_trPsi1*gtd(1,1)+Psi1(1,1)
      Psi1TF(1,2) = -third_trPsi1*gtd(1,2)+Psi1(1,2)
      Psi1TF(1,3) = -third_trPsi1*gtd(1,3)+Psi1(1,3)
      Psi1TF(2,1) = Psi1TF(1,2)
      Psi1TF(2,2) = -third_trPsi1*gtd(2,2)+Psi1(2,2)
      Psi1TF(2,3) = -third_trPsi1*gtd(2,3)+Psi1(2,3)
      Psi1TF(3,1) = Psi1TF(1,3)
      Psi1TF(3,2) = Psi1TF(2,3)
      Psi1TF(3,3) = -third_trPsi1*gtd(3,3)+Psi1(3,3)
      gtd_rhs(1,1) = 0.2D1*gtd(1,1)*d_Betau(1,1)+0.2D1*gtd(1,2)*d_Betau(
     #1,2)+0.2D1*gtd(1,3)*d_Betau(1,3)-0.6666666666666666D0*gtd(1,1)*div
     #_Beta-0.2D1*Alpha*Atd(1,1)
      gtd_rhs(1,2) = d_Betau(1,1)*gtd(1,2)+d_Betau(1,2)*gtd(2,2)+d_Betau
     #(1,3)*gtd(2,3)+d_Betau(2,1)*gtd(1,1)+d_Betau(2,2)*gtd(1,2)+d_Betau
     #(2,3)*gtd(1,3)-0.2D1*Alpha*Atd(1,2)-0.6666666666666666D0*div_Beta*
     #gtd(1,2)
      gtd_rhs(1,3) = d_Betau(1,1)*gtd(1,3)+d_Betau(1,2)*gtd(2,3)+d_Betau
     #(1,3)*gtd(3,3)+d_Betau(3,1)*gtd(1,1)+d_Betau(3,2)*gtd(1,2)+d_Betau
     #(3,3)*gtd(1,3)-0.2D1*Alpha*Atd(1,3)-0.6666666666666666D0*div_Beta*
     #gtd(1,3)
      gtd_rhs(2,1) = gtd_rhs(1,2)
      gtd_rhs(2,2) = 0.2D1*gtd(1,2)*d_Betau(2,1)+0.2D1*gtd(2,2)*d_Betau(
     #2,2)+0.2D1*gtd(2,3)*d_Betau(2,3)-0.6666666666666666D0*gtd(2,2)*div
     #_Beta-0.2D1*Alpha*Atd(2,2)
      gtd_rhs(2,3) = d_Betau(2,1)*gtd(1,3)+d_Betau(2,2)*gtd(2,3)+d_Betau
     #(2,3)*gtd(3,3)+d_Betau(3,1)*gtd(1,2)+d_Betau(3,2)*gtd(2,2)+d_Betau
     #(3,3)*gtd(2,3)-0.2D1*Alpha*Atd(2,3)-0.6666666666666666D0*div_Beta*
     #gtd(2,3)
      gtd_rhs(3,1) = gtd_rhs(1,3)
      gtd_rhs(3,2) = gtd_rhs(2,3)
      gtd_rhs(3,3) = 0.2D1*gtd(1,3)*d_Betau(3,1)+0.2D1*gtd(2,3)*d_Betau(
     #3,2)+0.2D1*gtd(3,3)*d_Betau(3,3)-0.6666666666666666D0*gtd(3,3)*div
     #_Beta-0.2D1*Alpha*Atd(3,3)
      t15 = Alpha*trK-0.6666666666666667D0*div_Beta
      Atd_rhs(1,1) = 0.2D1*Atd(1,1)*(-Alpha*Atud(1,1)+d_Betau(1,1))+0.2D
     #1*Atd(1,2)*(-Alpha*Atud(2,1)+d_Betau(1,2))+0.2D1*Atd(1,3)*(-Alpha*
     #Atud(3,1)+d_Betau(1,3))+Atd(1,1)*t15+Psi1TF(1,1)
      Atd_rhs(1,2) = Atd(1,1)*d_Betau(2,1)+Atd(1,2)*d_Betau(1,1)+Atd(1,2
     #)*d_Betau(2,2)+Atd(1,3)*d_Betau(2,3)+Atd(2,2)*d_Betau(1,2)+Atd(2,3
     #)*d_Betau(1,3)+Atd(1,2)*t15+Psi1TF(1,2)-0.2D1*Alpha*(Atd(1,1)*Atud
     #(1,2)+Atd(1,2)*Atud(2,2)+Atd(1,3)*Atud(3,2))
      Atd_rhs(1,3) = Atd(1,1)*d_Betau(3,1)+Atd(1,2)*d_Betau(3,2)+Atd(1,3
     #)*d_Betau(1,1)+Atd(1,3)*d_Betau(3,3)+Atd(2,3)*d_Betau(1,2)+Atd(3,3
     #)*d_Betau(1,3)+Atd(1,3)*t15+Psi1TF(1,3)-0.2D1*Alpha*(Atd(1,1)*Atud
     #(1,3)+Atd(1,2)*Atud(2,3)+Atd(1,3)*Atud(3,3))
      Atd_rhs(2,1) = Atd_rhs(1,2)
      Atd_rhs(2,2) = 0.2D1*Atd(1,2)*(-Alpha*Atud(1,2)+d_Betau(2,1))+0.2D
     #1*Atd(2,2)*(-Alpha*Atud(2,2)+d_Betau(2,2))+0.2D1*Atd(2,3)*(-Alpha*
     #Atud(3,2)+d_Betau(2,3))+Atd(2,2)*t15+Psi1TF(2,2)
      Atd_rhs(2,3) = Atd(1,2)*d_Betau(3,1)+Atd(1,3)*d_Betau(2,1)+Atd(2,2
     #)*d_Betau(3,2)+Atd(2,3)*d_Betau(2,2)+Atd(2,3)*d_Betau(3,3)+Atd(3,3
     #)*d_Betau(2,3)+Atd(2,3)*t15+Psi1TF(2,3)-0.2D1*Alpha*(Atd(1,2)*Atud
     #(1,3)+Atd(2,2)*Atud(2,3)+Atd(2,3)*Atud(3,3))
      Atd_rhs(3,1) = Atd_rhs(1,3)
      Atd_rhs(3,2) = Atd_rhs(2,3)
      Atd_rhs(3,3) = 0.2D1*Atd(1,3)*(-Alpha*Atud(1,3)+d_Betau(3,1))+0.2D
     #1*Atd(2,3)*(-Alpha*Atud(2,3)+d_Betau(3,2))+0.2D1*Atd(3,3)*(-Alpha*
     #Atud(3,3)+d_Betau(3,3))+Atd(3,3)*t15+Psi1TF(3,3)
      t30 = inv_chi*eight_pi_G
      t33 = t30*Jtd_ADM(1)+0.6666666666666667D0*d_trK(1)
      t37 = t30*Jtd_ADM(2)+0.6666666666666667D0*d_trK(2)
      t41 = t30*Jtd_ADM(3)+0.6666666666666667D0*d_trK(3)
      t46 = Alpha*inv_chi
      t49 = -0.15D1*t46*d_chi(1)-d_Alpha(1)
      t54 = -0.15D1*t46*d_chi(2)-d_Alpha(2)
      t59 = -0.15D1*t46*d_chi(3)-d_Alpha(3)
      t62 = 0.6666666666666667D0*CalGamt(1)*div_Beta-CalGamt(1)*d_Betau(
     #1,1)-CalGamt(2)*d_Betau(2,1)-CalGamt(3)*d_Betau(3,1)+gtu(1,1)*dd_B
     #etau(1,1,1)+gtu(2,2)*dd_Betau(2,2,1)+gtu(3,3)*dd_Betau(3,3,1)+0.2D
     #1*gtu(1,2)*dd_Betau(1,2,1)+0.2D1*gtu(1,3)*dd_Betau(1,3,1)+0.2D1*gt
     #u(2,3)*dd_Betau(2,3,1)+0.3333333333333333D0*gtu(1,1)*d_div_Beta(1)
     #+0.3333333333333333D0*gtu(1,2)*d_div_Beta(2)+0.3333333333333333D0*
     #gtu(1,3)*d_div_Beta(3)+0.2D1*Alpha*(0.5D0*Ct(1,1,1)*Atu(1,1)+0.5D0
     #*Ct(1,2,2)*Atu(2,2)+0.5D0*Ct(1,3,3)*Atu(3,3)+Ct(1,1,2)*Atu(1,2)+Ct
     #(1,1,3)*Atu(1,3)+Ct(1,2,3)*Atu(2,3)-gtu(1,1)*t33-gtu(1,2)*t37-gtu(
     #1,3)*t41)+0.2D1*Atu(1,1)*t49+0.2D1*Atu(1,2)*t54+0.2D1*Atu(1,3)*t59
      t104 = 0.6666666666666667D0*CalGamt(2)*div_Beta-CalGamt(1)*d_Betau
     #(1,2)-CalGamt(2)*d_Betau(2,2)-CalGamt(3)*d_Betau(3,2)+gtu(1,1)*dd_
     #Betau(1,1,2)+gtu(2,2)*dd_Betau(2,2,2)+gtu(3,3)*dd_Betau(3,3,2)+0.2
     #D1*gtu(1,2)*dd_Betau(1,2,2)+0.2D1*gtu(1,3)*dd_Betau(1,3,2)+0.2D1*g
     #tu(2,3)*dd_Betau(2,3,2)+0.3333333333333333D0*gtu(1,2)*d_div_Beta(1
     #)+0.3333333333333333D0*gtu(2,2)*d_div_Beta(2)+0.3333333333333333D0
     #*gtu(2,3)*d_div_Beta(3)+0.2D1*Alpha*(0.5D0*Ct(2,1,1)*Atu(1,1)+0.5D
     #0*Ct(2,2,2)*Atu(2,2)+0.5D0*Ct(2,3,3)*Atu(3,3)+Ct(2,1,2)*Atu(1,2)+C
     #t(2,1,3)*Atu(1,3)+Ct(2,2,3)*Atu(2,3)-gtu(1,2)*t33-gtu(2,2)*t37-gtu
     #(2,3)*t41)+0.2D1*Atu(1,2)*t49+0.2D1*Atu(2,2)*t54+0.2D1*Atu(2,3)*t5
     #9
      t146 = 0.6666666666666667D0*CalGamt(3)*div_Beta-CalGamt(1)*d_Betau
     #(1,3)-CalGamt(2)*d_Betau(2,3)-CalGamt(3)*d_Betau(3,3)+gtu(1,1)*dd_
     #Betau(1,1,3)+gtu(2,2)*dd_Betau(2,2,3)+gtu(3,3)*dd_Betau(3,3,3)+0.2
     #D1*gtu(1,2)*dd_Betau(1,2,3)+0.2D1*gtu(1,3)*dd_Betau(1,3,3)+0.2D1*g
     #tu(2,3)*dd_Betau(2,3,3)+0.3333333333333333D0*gtu(1,3)*d_div_Beta(1
     #)+0.3333333333333333D0*gtu(2,3)*d_div_Beta(2)+0.3333333333333333D0
     #*gtu(3,3)*d_div_Beta(3)+0.2D1*Alpha*(0.5D0*Ct(3,1,1)*Atu(1,1)+0.5D
     #0*Ct(3,2,2)*Atu(2,2)+0.5D0*Ct(3,3,3)*Atu(3,3)+Ct(3,1,2)*Atu(1,2)+C
     #t(3,1,3)*Atu(1,3)+Ct(3,2,3)*Atu(2,3)-gtu(1,3)*t33-gtu(2,3)*t37-gtu
     #(3,3)*t41)+0.2D1*Atu(1,3)*t49+0.2D1*Atu(2,3)*t54+0.2D1*Atu(3,3)*t5
     #9
      Gamt_rhs(1) = t62
      Gamt_rhs(2) = t104
      Gamt_rhs(3) = t146
      t2 = lambda_f1*Alpha+lambda_f0
      Betau_rhs(1) = 0.75D0*t2*Bu(1)
      Betau_rhs(2) = 0.75D0*t2*Bu(2)
      Betau_rhs(3) = 0.75D0*t2*Bu(3)
      Bu_rhs(1) = -feta*Bu(1)+Gamt_rhs(1)
      Bu_rhs(2) = -feta*Bu(2)+Gamt_rhs(2)
      Bu_rhs(3) = -feta*Bu(3)+Gamt_rhs(3)
      Alpha_rhs = -0.2D1*(lambda_f3*Alpha+lambda_f2)*Alpha*(trK-trK0)
      chi_rhs = -0.6666666666666667D0*chi*(-Alpha*trK+div_Beta)
      t11 = gtu(1,1)
      t14 = gtu(2,2)
      t17 = gtu(3,3)
      t20 = gtu(1,2)
      t24 = gtu(1,3)
      t28 = gtu(2,3)
      tr_pT = inv_chi*(t11*pTtd_ADM(1,1)+t14*pTtd_ADM(2,2)+t17*pTtd_ADM(
     #3,3)+0.2D1*t20*pTtd_ADM(1,2)+0.2D1*t24*pTtd_ADM(1,3)+0.2D1*t28*pTt
     #d_ADM(2,3))
      t54 = trK**2
      t61 = d_Alpha(1)
      t64 = d_Alpha(2)
      t67 = d_Alpha(3)
      t87 = d_chi(1)
      t91 = d_chi(2)
      t95 = d_chi(3)
      trK_rhs = Alpha*(Atd(1,1)*Atu(1,1)+Atd(2,2)*Atu(2,2)+Atd(3,3)*Atu(
     #3,3)+0.2D1*Atd(1,2)*Atu(1,2)+0.2D1*Atd(1,3)*Atu(1,3)+0.2D1*Atd(2,3
     #)*Atu(2,3)+0.3333333333333333D0*t54+four_pi_G*(rho_ADM+tr_pT))-chi
     #*(t11*dd_Alpha(1,1)+t14*dd_Alpha(2,2)+t17*dd_Alpha(3,3)+2*t20*dd_A
     #lpha(1,2)+2*t24*dd_Alpha(1,3)+2*t28*dd_Alpha(2,3)-CalGamt(1)*t61-C
     #alGamt(2)*t64-CalGamt(3)*t67)+0.5D0*t11*t61*t87+0.5D0*t20*t61*t91+
     #0.5D0*t24*t61*t95+0.5D0*t20*t64*t87+0.5D0*t14*t64*t91+0.5D0*t28*t6
     #4*t95+0.5D0*t24*t67*t87+0.5D0*t28*t67*t91+0.5D0*t17*t67*t95
      t5 = d_chi(1)
      t8 = d_chi(2)
      t11 = d_chi(3)
      divE = d_emEu(1,1)+d_emEu(2,2)+d_emEu(3,3)-0.15D1*inv_chi*(emEu(3)
     #*t11+emEu(1)*t5+emEu(2)*t8)
      t19 = emBu(1)
      t21 = emBu(2)
      t23 = emBu(3)
      divB = d_emBu(1,1)+d_emBu(2,2)+d_emBu(3,3)-0.15D1*inv_chi*(t23*t11
     #+t19*t5+t21*t8)
      t28 = t19**2
      t31 = t21**2
      t34 = t23**2
      Bsq = t28*gtd(1,1)+t31*gtd(2,2)+t34*gtd(3,3)+0.2D1*t19*t21*gtd(1,2
     #)+0.2D1*t19*t23*gtd(1,3)+0.2D1*t21*t23*gtd(2,3)
      inv_Bsq = 0.1D1/Bsq
      t50 = abs(chi)
      sqrt_chi = sqrt(t50)
      t10 = exp(0.2D1*dil_alpha*phiR)
      t12 = chi**(-0.3D1)
      t13 = sqrt(t12)
      t15 = 1/Alpha/t13
      Ju(1) = 0.2D1*dil_alpha*(-Femu(0,1)*d_phiR4d(0)+Femu(1,1)*d_phiR4d
     #(1)+Femu(1,2)*d_phiR4d(2)+Femu(1,3)*d_phiR4d(3))-t10*(-0.1D1*t15*F
     #emd(2,3)*d_phiI4d(0)+0.1D1*t15*Femd(0,3)*d_phiI4d(2)-0.1D1*t15*Fem
     #d(0,2)*d_phiI4d(3))
      Ju(2) = 0.2D1*dil_alpha*(-Femu(0,2)*d_phiR4d(0)-Femu(1,2)*d_phiR4d
     #(1)+Femu(2,2)*d_phiR4d(2)+Femu(2,3)*d_phiR4d(3))-t10*(0.1D1*t15*Fe
     #md(1,3)*d_phiI4d(0)-0.1D1*t15*Femd(0,3)*d_phiI4d(1)+0.1D1*t15*Femd
     #(0,1)*d_phiI4d(3))
      Ju(3) = 0.2D1*dil_alpha*(-Femu(0,3)*d_phiR4d(0)-Femu(1,3)*d_phiR4d
     #(1)-Femu(2,3)*d_phiR4d(2)+Femu(3,3)*d_phiR4d(3))-t10*(-0.1D1*t15*F
     #emd(1,2)*d_phiI4d(0)+0.1D1*t15*Femd(0,2)*d_phiI4d(1)-0.1D1*t15*Fem
     #d(0,1)*d_phiI4d(2))
      t19 = exp(0.2D1*axn_alpha*phiR)
      t22 = chi**(-0.3D1)
      t23 = sqrt(t22)
      t25 = 1/Alpha/t23
      qE = -0.2D1*dil_alpha*Alpha*(Femu(0,0)*d_phiR4d(0)-Femu(0,1)*d_phi
     #R4d(1)-Femu(0,2)*d_phiR4d(2)-Femu(0,3)*d_phiR4d(3))+0.2D1*t19*Alph
     #a*(-0.1D1*t25*Femd(2,3)*d_phiI4d(1)+0.1D1*t25*Femd(1,3)*d_phiI4d(2
     #)-0.1D1*t25*Femd(1,2)*d_phiI4d(3))
      t7 = emBu(1)*gtd(1,2)+emBu(2)*gtd(2,2)+emBu(3)*gtd(2,3)
      t8 = Alpha*inv_chi
      t10 = -t8*d_chi(3)+d_Alpha(3)
      t23 = emBu(1)*gtd(1,3)+emBu(2)*gtd(2,3)+emBu(3)*gtd(3,3)
      t25 = -t8*d_chi(2)+d_Alpha(2)
      t52 = emBu(1)*gtd(1,1)+emBu(2)*gtd(1,2)+emBu(3)*gtd(1,3)
      t63 = -t8*d_chi(1)+d_Alpha(1)
      emEu_rhs(1) = -emEu(1)*d_Betau(1,1)-emEu(2)*d_Betau(2,1)-emEu(3)*d
     #_Betau(3,1)+sqrt_chi*(-t7*t10-Alpha*(d_emBu(3,1)*gtd(1,2)+d_emBu(3
     #,2)*gtd(2,2)+d_emBu(3,3)*gtd(2,3)+d_gtd(3,1,2)*emBu(1)+d_gtd(3,2,2
     #)*emBu(2)+d_gtd(3,2,3)*emBu(3))+t23*t25+Alpha*(d_emBu(2,1)*gtd(1,3
     #)+d_emBu(2,2)*gtd(2,3)+d_emBu(2,3)*gtd(3,3)+d_gtd(2,1,3)*emBu(1)+d
     #_gtd(2,2,3)*emBu(2)+d_gtd(2,3,3)*emBu(3)))+Alpha*(trK*emEu(1)-chi*
     #(gtu(1,1)*d_emPsi(1)+gtu(1,2)*d_emPsi(2)+gtu(1,3)*d_emPsi(3))-Ju(1
     #))
      emEu_rhs(2) = -emEu(1)*d_Betau(1,2)-emEu(2)*d_Betau(2,2)-emEu(3)*d
     #_Betau(3,2)+sqrt_chi*(t52*t10+Alpha*(d_emBu(3,1)*gtd(1,1)+d_emBu(3
     #,2)*gtd(1,2)+d_emBu(3,3)*gtd(1,3)+d_gtd(3,1,1)*emBu(1)+d_gtd(3,1,2
     #)*emBu(2)+d_gtd(3,1,3)*emBu(3))-t23*t63-Alpha*(d_emBu(1,1)*gtd(1,3
     #)+d_emBu(1,2)*gtd(2,3)+d_emBu(1,3)*gtd(3,3)+d_gtd(1,1,3)*emBu(1)+d
     #_gtd(1,2,3)*emBu(2)+d_gtd(1,3,3)*emBu(3)))+Alpha*(trK*emEu(2)-chi*
     #(gtu(1,2)*d_emPsi(1)+gtu(2,2)*d_emPsi(2)+gtu(2,3)*d_emPsi(3))-Ju(2
     #))
      emEu_rhs(3) = -emEu(1)*d_Betau(1,3)-emEu(2)*d_Betau(2,3)-emEu(3)*d
     #_Betau(3,3)+sqrt_chi*(-t52*t25-Alpha*(d_emBu(2,1)*gtd(1,1)+d_emBu(
     #2,2)*gtd(1,2)+d_emBu(2,3)*gtd(1,3)+d_gtd(2,1,1)*emBu(1)+d_gtd(2,1,
     #2)*emBu(2)+d_gtd(2,1,3)*emBu(3))+t7*t63+Alpha*(d_emBu(1,1)*gtd(1,2
     #)+d_emBu(1,2)*gtd(2,2)+d_emBu(1,3)*gtd(2,3)+d_gtd(1,1,2)*emBu(1)+d
     #_gtd(1,2,2)*emBu(2)+d_gtd(1,2,3)*emBu(3)))+Alpha*(trK*emEu(3)-chi*
     #(gtu(1,3)*d_emPsi(1)+gtu(2,3)*d_emPsi(2)+gtu(3,3)*d_emPsi(3))-Ju(3
     #))
      t7 = emEu(1)*gtd(1,2)+emEu(2)*gtd(2,2)+emEu(3)*gtd(2,3)
      t8 = Alpha*inv_chi
      t10 = -t8*d_chi(3)+d_Alpha(3)
      t23 = emEu(1)*gtd(1,3)+emEu(2)*gtd(2,3)+emEu(3)*gtd(3,3)
      t25 = -t8*d_chi(2)+d_Alpha(2)
      t52 = emEu(1)*gtd(1,1)+emEu(2)*gtd(1,2)+emEu(3)*gtd(1,3)
      t63 = -t8*d_chi(1)+d_Alpha(1)
      emBu_rhs(1) = -emBu(1)*d_Betau(1,1)-emBu(2)*d_Betau(2,1)-emBu(3)*d
     #_Betau(3,1)-sqrt_chi*(-t7*t10-Alpha*(d_emEu(3,1)*gtd(1,2)+d_emEu(3
     #,2)*gtd(2,2)+d_emEu(3,3)*gtd(2,3)+d_gtd(3,1,2)*emEu(1)+d_gtd(3,2,2
     #)*emEu(2)+d_gtd(3,2,3)*emEu(3))+t23*t25+Alpha*(d_emEu(2,1)*gtd(1,3
     #)+d_emEu(2,2)*gtd(2,3)+d_emEu(2,3)*gtd(3,3)+d_gtd(2,1,3)*emEu(1)+d
     #_gtd(2,2,3)*emEu(2)+d_gtd(2,3,3)*emEu(3)))+Alpha*(trK*emBu(1)+chi*
     #(gtu(1,1)*d_emPhi(1)+gtu(1,2)*d_emPhi(2)+gtu(1,3)*d_emPhi(3)))
      emBu_rhs(2) = -emBu(1)*d_Betau(1,2)-emBu(2)*d_Betau(2,2)-emBu(3)*d
     #_Betau(3,2)-sqrt_chi*(t52*t10+Alpha*(d_emEu(3,1)*gtd(1,1)+d_emEu(3
     #,2)*gtd(1,2)+d_emEu(3,3)*gtd(1,3)+d_gtd(3,1,1)*emEu(1)+d_gtd(3,1,2
     #)*emEu(2)+d_gtd(3,1,3)*emEu(3))-t23*t63-Alpha*(d_emEu(1,1)*gtd(1,3
     #)+d_emEu(1,2)*gtd(2,3)+d_emEu(1,3)*gtd(3,3)+d_gtd(1,1,3)*emEu(1)+d
     #_gtd(1,2,3)*emEu(2)+d_gtd(1,3,3)*emEu(3)))+Alpha*(trK*emBu(2)+chi*
     #(gtu(1,2)*d_emPhi(1)+gtu(2,2)*d_emPhi(2)+gtu(2,3)*d_emPhi(3)))
      emBu_rhs(3) = -emBu(1)*d_Betau(1,3)-emBu(2)*d_Betau(2,3)-emBu(3)*d
     #_Betau(3,3)-sqrt_chi*(-t52*t25-Alpha*(d_emEu(2,1)*gtd(1,1)+d_emEu(
     #2,2)*gtd(1,2)+d_emEu(2,3)*gtd(1,3)+d_gtd(2,1,1)*emEu(1)+d_gtd(2,1,
     #2)*emEu(2)+d_gtd(2,1,3)*emEu(3))+t7*t63+Alpha*(d_emEu(1,1)*gtd(1,2
     #)+d_emEu(1,2)*gtd(2,2)+d_emEu(1,3)*gtd(2,3)+d_gtd(1,1,2)*emEu(1)+d
     #_gtd(1,2,2)*emEu(2)+d_gtd(1,2,3)*emEu(3)))+Alpha*(trK*emBu(3)+chi*
     #(gtu(1,3)*d_emPhi(1)+gtu(2,3)*d_emPhi(2)+gtu(3,3)*d_emPhi(3)))
      emPsi_rhs = Alpha*(-kappa_1*emPsi-divE+qE)
      emPhi_rhs = Alpha*(-kappa_2*emPhi+divB)
      phiR_rhs = -Alpha*piR
      phiI_rhs = -Alpha*piI
      t1 = gtu(1,1)
      t4 = gtu(1,2)
      t8 = gtu(1,3)
      t12 = gtu(2,2)
      t15 = gtu(2,3)
      t19 = gtu(3,3)
      t24 = CalGamt(1)
      t25 = d_phiR(1)
      t27 = CalGamt(2)
      t28 = d_phiR(2)
      t30 = CalGamt(3)
      t31 = d_phiR(3)
      t36 = d_chi(1)
      t46 = d_chi(2)
      t56 = d_chi(3)
      t66 = dVdphiR2*phiR
      t69 = exp(-0.2D1*dil_alpha*phiR)
      t73 = 0.4D1*axn_alpha*phiR
      t74 = exp(t73)
      t76 = piI**2
      t80 = -chi*(t1*dd_phiR(1,1)+t12*dd_phiR(2,2)+2*t15*dd_phiR(2,3)+t1
     #9*dd_phiR(3,3)+2*t4*dd_phiR(1,2)+2*t8*dd_phiR(1,3))+chi*(t24*t25+t
     #27*t28+t30*t31)+0.5D0*t1*t25*t36+0.5D0*t4*t28*t36+0.5D0*t8*t31*t36
     #+0.5D0*t4*t25*t46+0.5D0*t12*t28*t46+0.5D0*t15*t31*t46+0.5D0*t8*t25
     #*t56+0.5D0*t15*t28*t56+0.5D0*t19*t31*t56+trK*piR+t66-dil_alpha*t69
     #*Fsq+0.5D0*dil_alpha*t74*(-t76+dphiI4sq)
      t82 = d_Alpha(1)
      t83 = t1*t82
      t85 = t4*t82
      t87 = t8*t82
      t89 = d_Alpha(2)
      t90 = t4*t89
      t92 = t12*t89
      t94 = t15*t89
      t96 = d_Alpha(3)
      t97 = t8*t96
      t99 = t15*t96
      t101 = t19*t96
      piR_rhs = Alpha*t80-chi*(t101*t31+t83*t25+t90*t25+t97*t25+t85*t28+
     #t92*t28+t99*t28+t87*t31+t94*t31)
      t122 = d_phiI(1)
      t124 = d_phiI(2)
      t126 = d_phiI(3)
      t163 = exp(-t73)
      t164 = emEu(1)
      t165 = emBu(1)
      t169 = emEu(2)
      t171 = gtd(1,2)
      t173 = emEu(3)
      t175 = gtd(1,3)
      t177 = emBu(2)
      t184 = gtd(2,3)
      t186 = emBu(3)
      t198 = -chi*(t1*dd_phiI(1,1)+t12*dd_phiI(2,2)+2*t15*dd_phiI(2,3)+t
     #19*dd_phiI(3,3)+2*t4*dd_phiI(1,2)+2*t8*dd_phiI(1,3))+chi*(t24*t122
     #+t27*t124+t30*t126)+0.5D0*t1*t122*t36+0.5D0*t4*t124*t36+0.5D0*t8*t
     #126*t36+0.5D0*t4*t122*t46+0.5D0*t12*t124*t46+0.5D0*t15*t126*t46+0.
     #5D0*t8*t122*t56+0.5D0*t15*t124*t56+0.5D0*t19*t126*t56+trK*piI+t66+
     #0.4D1*axn_alpha*dphiIR4sq+0.4D1*axn_alpha*piR*piI+4*t163*(t164*t16
     #5*gtd(1,1)+t164*t177*t171+t164*t186*t175+t169*t165*t171+t173*t165*
     #t175+t169*t177*gtd(2,2)+t169*t186*t184+t173*t177*t184+t173*t186*gt
     #d(3,3))*inv_chi
      piI_rhs = Alpha*t198-chi*(t101*t126+t83*t122+t90*t122+t97*t122+t85
     #*t124+t92*t124+t99*t124+t87*t126+t94*t126)
