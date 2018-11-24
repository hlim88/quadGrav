      detgtd = gtd(1,1)*gtd(2,2)*gtd(3,3)-gtd(1,1)*gtd(3,2)**2-gtd(2,1)*
     &*2*gtd(3,3)+2*gtd(2,1)*gtd(3,1)*gtd(3,2)-gtd(3,1)**2*gtd(2,2)
      idetgtd = 0.1D1/detgtd
      gtu(1,1) = idetgtd*(gtd(2,2)*gtd(3,3)-gtd(3,2)**2)
      gtu(1,2) = idetgtd*(-gtd(2,1)*gtd(3,3)+gtd(3,1)*gtd(3,2))
      gtu(1,3) = idetgtd*(gtd(2,1)*gtd(3,2)-gtd(3,1)*gtd(2,2))
      gtu(2,1) = gtu(1,2)
      gtu(2,2) = idetgtd*(gtd(1,1)*gtd(3,3)-gtd(3,1)**2)
      gtu(2,3) = idetgtd*(-gtd(1,1)*gtd(3,2)+gtd(2,1)*gtd(3,1))
      gtu(3,1) = gtu(1,3)
      gtu(3,2) = gtu(2,3)
      gtu(3,3) = idetgtd*(gtd(1,1)*gtd(2,2)-gtd(2,1)**2)
      t3 = Dgtd(2,1,3)+Dgtd(3,1,2)-Dgtd(1,2,3)
      t8 = Dgtd(1,2,3)+Dgtd(3,1,2)-Dgtd(2,1,3)
      t13 = Dgtd(1,2,3)+Dgtd(2,1,3)-Dgtd(3,1,2)
      Ctd(1,1,1) = Dgtd(1,1,1)
      Ctd(1,1,2) = Dgtd(2,1,1)
      Ctd(1,1,3) = Dgtd(3,1,1)
      Ctd(1,2,1) = Dgtd(2,1,1)
      Ctd(1,2,2) = 2*Dgtd(2,1,2)-Dgtd(1,2,2)
      Ctd(1,2,3) = t3
      Ctd(1,3,1) = Dgtd(3,1,1)
      Ctd(1,3,2) = t3
      Ctd(1,3,3) = 2*Dgtd(3,1,3)-Dgtd(1,3,3)
      Ctd(2,1,1) = 2*Dgtd(1,1,2)-Dgtd(2,1,1)
      Ctd(2,1,2) = Dgtd(1,2,2)
      Ctd(2,1,3) = t8
      Ctd(2,2,1) = Dgtd(1,2,2)
      Ctd(2,2,2) = Dgtd(2,2,2)
      Ctd(2,2,3) = Dgtd(3,2,2)
      Ctd(2,3,1) = t8
      Ctd(2,3,2) = Dgtd(3,2,2)
      Ctd(2,3,3) = 2*Dgtd(3,2,3)-Dgtd(2,3,3)
      Ctd(3,1,1) = 2*Dgtd(1,1,3)-Dgtd(3,1,1)
      Ctd(3,1,2) = t13
      Ctd(3,1,3) = Dgtd(1,3,3)
      Ctd(3,2,1) = t13
      Ctd(3,2,2) = 2*Dgtd(2,2,3)-Dgtd(3,2,2)
      Ctd(3,2,3) = Dgtd(2,3,3)
      Ctd(3,3,1) = Dgtd(1,3,3)
      Ctd(3,3,2) = Dgtd(2,3,3)
      Ctd(3,3,3) = Dgtd(3,3,3)
      t8 = gtu(1,1)*Ctd(1,1,2)+gtu(1,2)*Ctd(2,1,2)+gtu(1,3)*Ctd(3,1,2)
      t12 = gtu(1,1)*Ctd(1,1,3)+gtu(1,2)*Ctd(2,1,3)+gtu(1,3)*Ctd(3,1,3)
      t20 = gtu(1,1)*Ctd(1,2,3)+gtu(1,2)*Ctd(2,2,3)+gtu(1,3)*Ctd(3,2,3)
      t32 = gtu(1,2)*Ctd(1,1,2)+gtu(2,2)*Ctd(2,1,2)+gtu(2,3)*Ctd(3,1,2)
      t36 = gtu(1,2)*Ctd(1,1,3)+gtu(2,2)*Ctd(2,1,3)+gtu(2,3)*Ctd(3,1,3)
      t44 = gtu(1,2)*Ctd(1,2,3)+gtu(2,2)*Ctd(2,2,3)+gtu(2,3)*Ctd(3,2,3)
      t56 = gtu(1,3)*Ctd(1,1,2)+gtu(2,3)*Ctd(2,1,2)+gtu(3,3)*Ctd(3,1,2)
      t60 = gtu(1,3)*Ctd(1,1,3)+gtu(2,3)*Ctd(2,1,3)+gtu(3,3)*Ctd(3,1,3)
      t68 = gtu(1,3)*Ctd(1,2,3)+gtu(2,3)*Ctd(2,2,3)+gtu(3,3)*Ctd(3,2,3)
      Ct(1,1,1) = gtu(1,1)*Ctd(1,1,1)+gtu(1,2)*Ctd(2,1,1)+gtu(1,3)*Ctd(3
     &,1,1)
      Ct(1,1,2) = t8
      Ct(1,1,3) = t12
      Ct(1,2,1) = t8
      Ct(1,2,2) = gtu(1,1)*Ctd(1,2,2)+gtu(1,2)*Ctd(2,2,2)+gtu(1,3)*Ctd(3
     &,2,2)
      Ct(1,2,3) = t20
      Ct(1,3,1) = t12
      Ct(1,3,2) = t20
      Ct(1,3,3) = gtu(1,1)*Ctd(1,3,3)+gtu(1,2)*Ctd(2,3,3)+gtu(1,3)*Ctd(3
     &,3,3)
      Ct(2,1,1) = gtu(1,2)*Ctd(1,1,1)+gtu(2,2)*Ctd(2,1,1)+gtu(2,3)*Ctd(3
     &,1,1)
      Ct(2,1,2) = t32
      Ct(2,1,3) = t36
      Ct(2,2,1) = t32
      Ct(2,2,2) = gtu(1,2)*Ctd(1,2,2)+gtu(2,2)*Ctd(2,2,2)+gtu(2,3)*Ctd(3
     &,2,2)
      Ct(2,2,3) = t44
      Ct(2,3,1) = t36
      Ct(2,3,2) = t44
      Ct(2,3,3) = gtu(1,2)*Ctd(1,3,3)+gtu(2,2)*Ctd(2,3,3)+gtu(2,3)*Ctd(3
     &,3,3)
      Ct(3,1,1) = gtu(1,3)*Ctd(1,1,1)+gtu(2,3)*Ctd(2,1,1)+gtu(3,3)*Ctd(3
     &,1,1)
      Ct(3,1,2) = t56
      Ct(3,1,3) = t60
      Ct(3,2,1) = t56
      Ct(3,2,2) = gtu(1,3)*Ctd(1,2,2)+gtu(2,3)*Ctd(2,2,2)+gtu(3,3)*Ctd(3
     &,2,2)
      Ct(3,2,3) = t68
      Ct(3,3,1) = t60
      Ct(3,3,2) = t68
      Ct(3,3,3) = gtu(1,3)*Ctd(1,3,3)+gtu(2,3)*Ctd(2,3,3)+gtu(3,3)*Ctd(3
     &,3,3)
      Gamt(1) = half*(gtu(1,1)*Ct(1,1,1)+gtu(2,2)*Ct(1,2,2)+gtu(3,3)*Ct(
     &1,3,3))+gtu(1,2)*Ct(1,1,2)+gtu(1,3)*Ct(1,1,3)+gtu(2,3)*Ct(1,2,3)
      Gamt(2) = half*(gtu(1,1)*Ct(2,1,1)+gtu(2,2)*Ct(2,2,2)+gtu(3,3)*Ct(
     &2,3,3))+gtu(1,2)*Ct(2,1,2)+gtu(1,3)*Ct(2,1,3)+gtu(2,3)*Ct(2,2,3)
      Gamt(3) = half*(gtu(1,1)*Ct(3,1,1)+gtu(2,2)*Ct(3,2,2)+gtu(3,3)*Ct(
     &3,3,3))+gtu(1,2)*Ct(3,1,2)+gtu(1,3)*Ct(3,1,3)+gtu(2,3)*Ct(3,2,3)
