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
