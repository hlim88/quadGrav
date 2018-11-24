      t1 = gd(1,1)
      t2 = gd(2,2)
      t4 = gd(3,3)
      t6 = gd(2,3)
      t7 = t6**2
      t9 = gd(1,2)
      t10 = t9**2
      t12 = gd(1,3)
      t16 = t12**2
      detgd = t1*t2*t4-t1*t7-t10*t4+2*t9*t12*t6-t16*t2
      idetgd = 0.1D1/detgd
      gu(1,1) = idetgd*(gd(2,2)*gd(3,3)-gd(2,3)**2)
      gu(1,2) = idetgd*(-gd(1,2)*gd(3,3)+gd(1,3)*gd(2,3))
      gu(1,3) = idetgd*(gd(1,2)*gd(2,3)-gd(1,3)*gd(2,2))
      gu(2,1) = idetgd*(-gd(1,2)*gd(3,3)+gd(1,3)*gd(2,3))
      gu(2,2) = idetgd*(gd(1,1)*gd(3,3)-gd(1,3)**2)
      gu(2,3) = idetgd*(-gd(1,1)*gd(2,3)+gd(1,2)*gd(1,3))
      gu(3,1) = idetgd*(gd(1,2)*gd(2,3)-gd(1,3)*gd(2,2))
      gu(3,2) = idetgd*(-gd(1,1)*gd(2,3)+gd(1,2)*gd(1,3))
      gu(3,3) = idetgd*(gd(1,1)*gd(2,2)-gd(1,2)**2)
      chi = idetgd**third
      gtd(1,1) = chi*gd(1,1)
      gtd(1,2) = chi*gd(1,2)
      gtd(1,3) = chi*gd(1,3)
      gtd(2,1) = gtd(1,2)
      gtd(2,2) = chi*gd(2,2)
      gtd(2,3) = chi*gd(2,3)
      gtd(3,1) = gtd(1,3)
      gtd(3,2) = gtd(2,3)
      gtd(3,3) = chi*gd(3,3)
      detgtd = gtd(1,1)*gtd(2,2)*gtd(3,3)-gtd(1,1)*gtd(3,2)**2-gtd(2,1)*
     &*2*gtd(3,3)+2*gtd(2,1)*gtd(3,1)*gtd(3,2)-gtd(3,1)**2*gtd(2,2)
      idetgtd = 0.1D1/detgtd
      trK = gu(1,1)*Kd(1,1)+gu(2,2)*Kd(2,2)+gu(3,3)*Kd(3,3)+two*(gu(1,2)
     &*Kd(1,2)+gu(1,3)*Kd(1,3)+gu(2,3)*Kd(2,3))
      Atd(1,1) = chi*(Kd(1,1)-third*gd(1,1)*TrK)
      Atd(1,2) = chi*(Kd(1,2)-third*gd(1,2)*TrK)
      Atd(1,3) = chi*(Kd(1,3)-third*gd(1,3)*TrK)
      Atd(2,1) = Atd(1,2)
      Atd(2,2) = chi*(Kd(2,2)-third*gd(2,2)*TrK)
      Atd(2,3) = chi*(Kd(2,3)-third*gd(2,3)*TrK)
      Atd(3,1) = Atd(1,3)
      Atd(3,2) = Atd(2,3)
      Atd(3,3) = chi*(Kd(3,3)-third*gd(3,3)*TrK)
