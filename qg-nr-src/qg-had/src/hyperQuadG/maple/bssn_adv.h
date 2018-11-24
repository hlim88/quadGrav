      Gamt_rhs_xxx(1) = Betau(1)*adv_d_Gamt(1,1)+Betau(2)*adv_d_Gamt(2,1
     #)+Betau(3)*adv_d_Gamt(3,1)+Gamt_rhs(1)
      Gamt_rhs_xxx(2) = Betau(1)*adv_d_Gamt(1,2)+Betau(2)*adv_d_Gamt(2,2
     #)+Betau(3)*adv_d_Gamt(3,2)+Gamt_rhs(2)
      Gamt_rhs_xxx(3) = Betau(1)*adv_d_Gamt(1,3)+Betau(2)*adv_d_Gamt(2,3
     #)+Betau(3)*adv_d_Gamt(3,3)+Gamt_rhs(3)
      Atd_rhs_xxx(1,1) = Betau(1)*adv_d_Atd(1,1,1)+Betau(2)*adv_d_Atd(2,
     #1,1)+Betau(3)*adv_d_Atd(3,1,1)+Atd_rhs(1,1)
      Atd_rhs_xxx(1,2) = Betau(1)*adv_d_Atd(1,1,2)+Betau(2)*adv_d_Atd(2,
     #1,2)+Betau(3)*adv_d_Atd(3,1,2)+Atd_rhs(1,2)
      Atd_rhs_xxx(1,3) = Betau(1)*adv_d_Atd(1,1,3)+Betau(2)*adv_d_Atd(2,
     #1,3)+Betau(3)*adv_d_Atd(3,1,3)+Atd_rhs(1,3)
      Atd_rhs_xxx(2,1) = Atd_rhs_xxx(1,2)
      Atd_rhs_xxx(2,2) = Betau(1)*adv_d_Atd(1,2,2)+Betau(2)*adv_d_Atd(2,
     #2,2)+Betau(3)*adv_d_Atd(3,2,2)+Atd_rhs(2,2)
      Atd_rhs_xxx(2,3) = Betau(1)*adv_d_Atd(1,2,3)+Betau(2)*adv_d_Atd(2,
     #2,3)+Betau(3)*adv_d_Atd(3,2,3)+Atd_rhs(2,3)
      Atd_rhs_xxx(3,1) = Atd_rhs_xxx(1,3)
      Atd_rhs_xxx(3,2) = Atd_rhs_xxx(2,3)
      Atd_rhs_xxx(3,3) = Betau(1)*adv_d_Atd(1,3,3)+Betau(2)*adv_d_Atd(2,
     #3,3)+Betau(3)*adv_d_Atd(3,3,3)+Atd_rhs(3,3)
      gtd_rhs_xxx(1,1) = Betau(1)*adv_d_gtd(1,1,1)+Betau(2)*adv_d_gtd(2,
     #1,1)+Betau(3)*adv_d_gtd(3,1,1)+gtd_rhs(1,1)
      gtd_rhs_xxx(1,2) = Betau(1)*adv_d_gtd(1,1,2)+Betau(2)*adv_d_gtd(2,
     #1,2)+Betau(3)*adv_d_gtd(3,1,2)+gtd_rhs(1,2)
      gtd_rhs_xxx(1,3) = Betau(1)*adv_d_gtd(1,1,3)+Betau(2)*adv_d_gtd(2,
     #1,3)+Betau(3)*adv_d_gtd(3,1,3)+gtd_rhs(1,3)
      gtd_rhs_xxx(2,1) = gtd_rhs_xxx(1,2)
      gtd_rhs_xxx(2,2) = Betau(1)*adv_d_gtd(1,2,2)+Betau(2)*adv_d_gtd(2,
     #2,2)+Betau(3)*adv_d_gtd(3,2,2)+gtd_rhs(2,2)
      gtd_rhs_xxx(2,3) = Betau(1)*adv_d_gtd(1,2,3)+Betau(2)*adv_d_gtd(2,
     #2,3)+Betau(3)*adv_d_gtd(3,2,3)+gtd_rhs(2,3)
      gtd_rhs_xxx(3,1) = gtd_rhs_xxx(1,3)
      gtd_rhs_xxx(3,2) = gtd_rhs_xxx(2,3)
      gtd_rhs_xxx(3,3) = Betau(1)*adv_d_gtd(1,3,3)+Betau(2)*adv_d_gtd(2,
     #3,3)+Betau(3)*adv_d_gtd(3,3,3)+gtd_rhs(3,3)
      Betau_rhs_xxx(1) = Betau_rhs(1)+lambda_2*(Betau(1)*adv_d_Betau(1,1
     #)+Betau(2)*adv_d_Betau(2,1)+Betau(3)*adv_d_Betau(3,1))
      Betau_rhs_xxx(2) = Betau_rhs(2)+lambda_2*(Betau(1)*adv_d_Betau(1,2
     #)+Betau(2)*adv_d_Betau(2,2)+Betau(3)*adv_d_Betau(3,2))
      Betau_rhs_xxx(3) = Betau_rhs(3)+lambda_2*(Betau(1)*adv_d_Betau(1,3
     #)+Betau(2)*adv_d_Betau(2,3)+Betau(3)*adv_d_Betau(3,3))
      t2 = 0.1D1-lambda_4
      Bu_rhs_xxx(1) = Bu_rhs(1)+Betau(1)*(t2*adv_d_Gamt(1,1)+lambda_3*ad
     #v_d_Bu(1,1))+Betau(2)*(t2*adv_d_Gamt(2,1)+lambda_3*adv_d_Bu(2,1))+
     #Betau(3)*(t2*adv_d_Gamt(3,1)+lambda_3*adv_d_Bu(3,1))
      Bu_rhs_xxx(2) = Bu_rhs(2)+Betau(1)*(t2*adv_d_Gamt(1,2)+lambda_3*ad
     #v_d_Bu(1,2))+Betau(2)*(t2*adv_d_Gamt(2,2)+lambda_3*adv_d_Bu(2,2))+
     #Betau(3)*(t2*adv_d_Gamt(3,2)+lambda_3*adv_d_Bu(3,2))
      Bu_rhs_xxx(3) = Bu_rhs(3)+Betau(1)*(t2*adv_d_Gamt(1,3)+lambda_3*ad
     #v_d_Bu(1,3))+Betau(2)*(t2*adv_d_Gamt(2,3)+lambda_3*adv_d_Bu(2,3))+
     #Betau(3)*(t2*adv_d_Gamt(3,3)+lambda_3*adv_d_Bu(3,3))
      t1 = Betau(1)
      t4 = Betau(2)
      t7 = Betau(3)
      chi_rhs_xxx = t1*adv_d_chi(1)+t4*adv_d_chi(2)+t7*adv_d_chi(3)+chi_
     #rhs
      trK_rhs_xxx = t1*adv_d_trK(1)+t4*adv_d_trK(2)+t7*adv_d_trK(3)+trK_
     #rhs
      Alpha_rhs_xxx = Alpha_rhs+lambda_1*(t1*adv_d_Alpha(1)+t4*adv_d_Alp
     #ha(2)+t7*adv_d_Alpha(3))
