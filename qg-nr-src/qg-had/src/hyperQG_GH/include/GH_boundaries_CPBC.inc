      
     

 
      !--------------------------------------------------------------
      !---------------------THE EASY PART----------------------------
      !--------------------------------------------------------------

      !------maximally dissipative for the L----------------
      if (CPBC_type .EQ. 0) then
      
        sDt100 = 0.0d0
        sDt101 = 0.0d0
        sDt102 = 0.0d0
        sDt103 = 0.0d0
        sDt111 = 0.0d0
        sDt112 = 0.0d0
        sDt113 = 0.0d0
        sDt122 = 0.0d0
        sDt123 = 0.0d0
        sDt133 = 0.0d0

        sDt200 = 0.0d0
        sDt201 = 0.0d0
        sDt202 = 0.0d0
        sDt203 = 0.0d0
        sDt211 = 0.0d0
        sDt212 = 0.0d0
        sDt213 = 0.0d0
        sDt222 = 0.0d0
        sDt223 = 0.0d0
        sDt233 = 0.0d0

        sDt300 = 0.0d0
        sDt301 = 0.0d0
        sDt302 = 0.0d0
        sDt303 = 0.0d0
        sDt311 = 0.0d0
        sDt312 = 0.0d0
        sDt313 = 0.0d0
        sDt322 = 0.0d0
        sDt323 = 0.0d0
        sDt333 = 0.0d0
                        
        sdtphim1 = 0.0d0
        sdtphim2 = 0.0d0
        sdtphim3 = 0.0d0
        sdtphin1 = 0.0d0
        sdtphin2 = 0.0d0
        sdtphin3 = 0.0d0
        sdtphir1 = 0.0d0
        sdtphir2 = 0.0d0
        sdtphir3 = 0.0d0
        sdtphic1 = 0.0d0
        sdtphic2 = 0.0d0
        sdtphic3 = 0.0d0
      
        sdt1H0 = 0.0d0
        sdt2H0 = 0.0d0
        sdt3H0 = 0.0d0
      
      else      
      
      !- do it before changing the definition of the n --------------
      !--to compute the outgoing/ingoing null vectors ---------------
      
      !--- the first order constraints with 3 indexes--------------
      !CD_iab = d_i g_ab - D_iab
 
      ! it is commented because it is not being used right now
           
!      CD100 = dx_g00 - D100
!      CD101 = dx_g01 - D101
!      CD102 = dx_g02 - D102
!      CD103 = dx_g03 - D103
!      CD111 = dx_g11 - D111
!      CD112 = dx_g12 - D112
!      CD113 = dx_g13 - D113
!      CD122 = dx_g22 - D122
!      CD123 = dx_g23 - D123                                                
!      CD133 = dx_g33 - D133 

!      CD200 = dy_g00 - D200
!      CD201 = dy_g01 - D201
!      CD202 = dy_g02 - D202
!      CD203 = dy_g03 - D203
!      CD211 = dy_g11 - D211
!      CD212 = dy_g12 - D212
!      CD213 = dy_g13 - D213
!      CD222 = dy_g22 - D222
!      CD223 = dy_g23 - D223                                                
!      CD233 = dy_g33 - D233 

!      CD300 = dz_g00 - D300
!      CD301 = dz_g01 - D301
!      CD302 = dz_g02 - D302
!      CD303 = dz_g03 - D303
!      CD311 = dz_g11 - D311
!      CD312 = dz_g12 - D312
!      CD313 = dz_g13 - D313
!      CD322 = dz_g22 - D322
!      CD323 = dz_g23 - D323                                                
!      CD333 = dz_g33 - D333 

                  
      !--- the first order constraints with 4 indexes--------------                 
      !CDD_ijab = di_Djab - dj_Diab
      
      CDD1200 = dx_D200 - dy_D100
      CDD1201 = dx_D201 - dy_D101
      CDD1202 = dx_D202 - dy_D102
      CDD1203 = dx_D203 - dy_D103
      CDD1211 = dx_D211 - dy_D111
      CDD1212 = dx_D212 - dy_D112
      CDD1213 = dx_D213 - dy_D113
      CDD1222 = dx_D222 - dy_D122
      CDD1223 = dx_D223 - dy_D123
      CDD1233 = dx_D233 - dy_D133

      CDD1300 = dx_D300 - dz_D100
      CDD1301 = dx_D301 - dz_D101
      CDD1302 = dx_D302 - dz_D102
      CDD1303 = dx_D303 - dz_D103
      CDD1311 = dx_D311 - dz_D111
      CDD1312 = dx_D312 - dz_D112
      CDD1313 = dx_D313 - dz_D113
      CDD1322 = dx_D322 - dz_D122
      CDD1323 = dx_D323 - dz_D123
      CDD1333 = dx_D333 - dz_D133

      CDD2300 = dy_D300 - dz_D200
      CDD2301 = dy_D301 - dz_D201
      CDD2302 = dy_D302 - dz_D202
      CDD2303 = dy_D303 - dz_D203
      CDD2311 = dy_D311 - dz_D211
      CDD2312 = dy_D312 - dz_D212
      CDD2313 = dy_D313 - dz_D213
      CDD2322 = dy_D322 - dz_D222
      CDD2323 = dy_D323 - dz_D223
      CDD2333 = dy_D333 - dz_D233

      CDD1100 = 0.0; CDD1101 = 0.0; CDD1102 = 0.0; CDD1103 = 0.0; CDD1111 = 0.0
      CDD1112 = 0.0; CDD1113 = 0.0; CDD1122 = 0.0; CDD1123 = 0.0; CDD1133 = 0.0
      CDD2200 = 0.0; CDD2201 = 0.0; CDD2202 = 0.0; CDD2203 = 0.0; CDD2211 = 0.0
      CDD2212 = 0.0; CDD2213 = 0.0; CDD2222 = 0.0; CDD2223 = 0.0; CDD2233 = 0.0
      CDD3300 = 0.0; CDD3301 = 0.0; CDD3302 = 0.0; CDD3303 = 0.0; CDD3311 = 0.0
      CDD3312 = 0.0; CDD3313 = 0.0; CDD3322 = 0.0; CDD3323 = 0.0; CDD3333 = 0.0
                                  
      CDD2100 = - CDD1200; CDD2101 = - CDD1201; CDD2102 = - CDD1202
      CDD2103 = - CDD1203; CDD2111 = - CDD1211; CDD2112 = - CDD1212
      CDD2113 = - CDD1213; CDD2122 = - CDD1222; CDD2123 = - CDD1223
      CDD2133 = - CDD1233                                                      

      CDD3100 = - CDD1300; CDD3101 = - CDD1301; CDD3102 = - CDD1302
      CDD3103 = - CDD1303; CDD3111 = - CDD1311; CDD3112 = - CDD1312
      CDD3113 = - CDD1313; CDD3122 = - CDD1322; CDD3123 = - CDD1323
      CDD3133 = - CDD1333                                                      
      
      CDD3200 = - CDD2300; CDD3201 = - CDD2301; CDD3202 = - CDD2302
      CDD3203 = - CDD2303; CDD3211 = - CDD2311; CDD3212 = - CDD2312
      CDD3213 = - CDD2313; CDD3222 = - CDD2322; CDD3223 = - CDD2323
      CDD3233 = - CDD2333                                                      

!     the n projection -----------------------------
      CDDn100 = nu1*CDD1100 + nu2*CDD2100 + nu3*CDD3100       
      CDDn101 = nu1*CDD1101 + nu2*CDD2101 + nu3*CDD3101       
      CDDn102 = nu1*CDD1102 + nu2*CDD2102 + nu3*CDD3102
      CDDn103 = nu1*CDD1103 + nu2*CDD2103 + nu3*CDD3103
      CDDn111 = nu1*CDD1111 + nu2*CDD2111 + nu3*CDD3111
      CDDn112 = nu1*CDD1112 + nu2*CDD2112 + nu3*CDD3112
      CDDn113 = nu1*CDD1113 + nu2*CDD2113 + nu3*CDD3113
      CDDn122 = nu1*CDD1122 + nu2*CDD2122 + nu3*CDD3122
      CDDn123 = nu1*CDD1123 + nu2*CDD2123 + nu3*CDD3123
      CDDn133 = nu1*CDD1133 + nu2*CDD2133 + nu3*CDD3133

      CDDn200 = nu1*CDD1200 + nu2*CDD2200 + nu3*CDD3200       
      CDDn201 = nu1*CDD1201 + nu2*CDD2201 + nu3*CDD3201       
      CDDn202 = nu1*CDD1202 + nu2*CDD2202 + nu3*CDD3202
      CDDn203 = nu1*CDD1203 + nu2*CDD2203 + nu3*CDD3203
      CDDn211 = nu1*CDD1211 + nu2*CDD2211 + nu3*CDD3211
      CDDn212 = nu1*CDD1212 + nu2*CDD2212 + nu3*CDD3212
      CDDn213 = nu1*CDD1213 + nu2*CDD2213 + nu3*CDD3213
      CDDn222 = nu1*CDD1222 + nu2*CDD2222 + nu3*CDD3222
      CDDn223 = nu1*CDD1223 + nu2*CDD2223 + nu3*CDD3223
      CDDn233 = nu1*CDD1233 + nu2*CDD2233 + nu3*CDD3233

      CDDn300 = nu1*CDD1300 + nu2*CDD2300 + nu3*CDD3300       
      CDDn301 = nu1*CDD1301 + nu2*CDD2301 + nu3*CDD3301       
      CDDn302 = nu1*CDD1302 + nu2*CDD2302 + nu3*CDD3302
      CDDn303 = nu1*CDD1303 + nu2*CDD2303 + nu3*CDD3303
      CDDn311 = nu1*CDD1311 + nu2*CDD2311 + nu3*CDD3311
      CDDn312 = nu1*CDD1312 + nu2*CDD2312 + nu3*CDD3312
      CDDn313 = nu1*CDD1313 + nu2*CDD2313 + nu3*CDD3313
      CDDn322 = nu1*CDD1322 + nu2*CDD2322 + nu3*CDD3322
      CDDn323 = nu1*CDD1323 + nu2*CDD2323 + nu3*CDD3323
      CDDn333 = nu1*CDD1333 + nu2*CDD2333 + nu3*CDD3333
                   

      !CPPij= di_djphi - dj_diphi-----------------------
      CPPm12 = dx_d2phim - dy_d1phim 
      CPPm13 = dx_d3phim - dz_d1phim
      CPPm23 = dy_d3phim - dz_d2phim

      CPPn12 = dx_d2phin - dy_d1phin 
      CPPn13 = dx_d3phin - dz_d1phin
      CPPn23 = dy_d3phin - dz_d2phin
      
      CPPr12 = dx_d2phir - dy_d1phir 
      CPPr13 = dx_d3phir - dz_d1phir
      CPPr23 = dy_d3phir - dz_d2phir

      CPPc12 = dx_d2phic - dy_d1phic 
      CPPc13 = dx_d3phic - dz_d1phic
      CPPc23 = dy_d3phic - dz_d2phic
            
      CPPm11 = 0.0; CPPm22 = 0.0; CPPm33 = 0.0
      CPPn11 = 0.0; CPPn22 = 0.0; CPPn33 = 0.0
      CPPr11 = 0.0; CPPr22 = 0.0; CPPr33 = 0.0
      CPPc11 = 0.0; CPPc22 = 0.0; CPPc33 = 0.0
                  
      CPPm21 = - CPPm12; CPPm31 = - CPPm13; CPPm32 = - CPPm23
      CPPn21 = - CPPn12; CPPn31 = - CPPn13; CPPn32 = - CPPn23
      CPPr21 = - CPPr12; CPPr31 = - CPPr13; CPPr32 = - CPPr23
      CPPc21 = - CPPc12; CPPc31 = - CPPc13; CPPc32 = - CPPc23
      
      !the n projection-------------------------
      CPPmn1 = nu1*CPPm11 + nu2*CPPm21 + nu3*CPPm31
      CPPmn2 = nu1*CPPm12 + nu2*CPPm22 + nu3*CPPm32
      CPPmn3 = nu1*CPPm13 + nu2*CPPm23 + nu3*CPPm33

      CPPnn1 = nu1*CPPn11 + nu2*CPPn21 + nu3*CPPn31
      CPPnn2 = nu1*CPPn12 + nu2*CPPn22 + nu3*CPPn32
      CPPnn3 = nu1*CPPn13 + nu2*CPPn23 + nu3*CPPn33

      CPPrn1 = nu1*CPPr11 + nu2*CPPr21 + nu3*CPPr31
      CPPrn2 = nu1*CPPr12 + nu2*CPPr22 + nu3*CPPr32
      CPPrn3 = nu1*CPPr13 + nu2*CPPr23 + nu3*CPPr33

      CPPcn1 = nu1*CPPc11 + nu2*CPPc21 + nu3*CPPc31
      CPPcn2 = nu1*CPPc12 + nu2*CPPc22 + nu3*CPPc32
      CPPcn3 = nu1*CPPc13 + nu2*CPPc23 + nu3*CPPc33

      !CGGij= di_djHo - dj_diHo-----------------------
      CGGm12 = dx_d2H0 - dy_d1H0 
      CGGm13 = dx_d3H0 - dz_d1H0
      CGGm23 = dy_d3H0 - dz_d2H0
                  
      CGGm11 = 0.0; CGGm22 = 0.0; CGGm33 = 0.0                  
      CGGm21 = - CGGm12; CGGm31 = - CGGm13; CGGm32 = - CGGm23
      
      !the n projection-------------------------
      CGGmn1 = nu1*CGGm11 + nu2*CGGm21 + nu3*CGGm31
      CGGmn2 = nu1*CGGm12 + nu2*CGGm22 + nu3*CGGm32
      CGGmn3 = nu1*CGGm13 + nu2*CGGm23 + nu3*CGGm33
      
            
      !-------------------------------------------
      !--- the 3-dimensional projector------------
      !-------------------------------------------
      
      Pud11 = 1.0d0 - nu1*nd1
      Pud12 = - nu1*nd2
      Pud13 = - nu1*nd3
      Pud21 = - nu2*nd1
      Pud22 = 1.0d0 - nu2*nd2
      Pud23 = - nu2*nd3
      Pud31 = - nu3*nd1
      Pud32 = - nu3*nd2
      Pud33 = 1.0d0 - nu3*nd3
 
      sDt100 = Dt100 - bn*(Pud11*CDDn100 + Pud21*CDDn200 + Pud31*CDDn300)
      sDt101 = Dt101 - bn*(Pud11*CDDn101 + Pud21*CDDn201 + Pud31*CDDn301)
      sDt102 = Dt102 - bn*(Pud11*CDDn102 + Pud21*CDDn202 + Pud31*CDDn302)
      sDt103 = Dt103 - bn*(Pud11*CDDn103 + Pud21*CDDn203 + Pud31*CDDn303)
      sDt111 = Dt111 - bn*(Pud11*CDDn111 + Pud21*CDDn211 + Pud31*CDDn311)
      sDt112 = Dt112 - bn*(Pud11*CDDn112 + Pud21*CDDn212 + Pud31*CDDn312)
      sDt113 = Dt113 - bn*(Pud11*CDDn113 + Pud21*CDDn213 + Pud31*CDDn313)
      sDt122 = Dt122 - bn*(Pud11*CDDn122 + Pud21*CDDn222 + Pud31*CDDn322)
      sDt123 = Dt123 - bn*(Pud11*CDDn123 + Pud21*CDDn223 + Pud31*CDDn323)
      sDt133 = Dt133 - bn*(Pud11*CDDn133 + Pud21*CDDn233 + Pud31*CDDn333)

      sDt200 = Dt200 - bn*(Pud12*CDDn100 + Pud22*CDDn200 + Pud32*CDDn300)
      sDt201 = Dt201 - bn*(Pud12*CDDn101 + Pud22*CDDn201 + Pud32*CDDn301)
      sDt202 = Dt202 - bn*(Pud12*CDDn102 + Pud22*CDDn202 + Pud32*CDDn302)
      sDt203 = Dt203 - bn*(Pud12*CDDn103 + Pud22*CDDn203 + Pud32*CDDn303)
      sDt211 = Dt211 - bn*(Pud12*CDDn111 + Pud22*CDDn211 + Pud32*CDDn311)
      sDt212 = Dt212 - bn*(Pud12*CDDn112 + Pud22*CDDn212 + Pud32*CDDn312)
      sDt213 = Dt213 - bn*(Pud12*CDDn113 + Pud22*CDDn213 + Pud32*CDDn313)
      sDt222 = Dt222 - bn*(Pud12*CDDn122 + Pud22*CDDn222 + Pud32*CDDn322)
      sDt223 = Dt223 - bn*(Pud12*CDDn123 + Pud22*CDDn223 + Pud32*CDDn323)
      sDt233 = Dt233 - bn*(Pud12*CDDn133 + Pud22*CDDn233 + Pud32*CDDn333)

      sDt300 = Dt300 - bn*(Pud13*CDDn100 + Pud23*CDDn200 + Pud33*CDDn300)
      sDt301 = Dt301 - bn*(Pud13*CDDn101 + Pud23*CDDn201 + Pud33*CDDn301)
      sDt302 = Dt302 - bn*(Pud13*CDDn102 + Pud23*CDDn202 + Pud33*CDDn302)
      sDt303 = Dt303 - bn*(Pud13*CDDn103 + Pud23*CDDn203 + Pud33*CDDn303)
      sDt311 = Dt311 - bn*(Pud13*CDDn111 + Pud23*CDDn211 + Pud33*CDDn311)
      sDt312 = Dt312 - bn*(Pud13*CDDn112 + Pud23*CDDn212 + Pud33*CDDn312)
      sDt313 = Dt313 - bn*(Pud13*CDDn113 + Pud23*CDDn213 + Pud33*CDDn313)
      sDt322 = Dt322 - bn*(Pud13*CDDn122 + Pud23*CDDn222 + Pud33*CDDn322)    
      sDt323 = Dt323 - bn*(Pud13*CDDn123 + Pud23*CDDn223 + Pud33*CDDn323)        
      sDt333 = Dt333 - bn*(Pud13*CDDn133 + Pud23*CDDn233 + Pud33*CDDn333)
                        
      sdtphim1 = dtphim1 - bn*(Pud11*CPPmn1 + Pud21*CPPmn2 + Pud31*CPPmn3)
      sdtphim2 = dtphim2 - bn*(Pud12*CPPmn1 + Pud22*CPPmn2 + Pud32*CPPmn3)
      sdtphim3 = dtphim3 - bn*(Pud13*CPPmn1 + Pud23*CPPmn2 + Pud33*CPPmn3)
      sdtphin1 = dtphin1 - bn*(Pud11*CPPnn1 + Pud21*CPPnn2 + Pud31*CPPnn3)
      sdtphin2 = dtphin2 - bn*(Pud12*CPPnn1 + Pud22*CPPnn2 + Pud32*CPPnn3)
      sdtphin3 = dtphin3 - bn*(Pud13*CPPnn1 + Pud23*CPPnn2 + Pud33*CPPnn3)
      sdtphir1 = dtphir1 - bn*(Pud11*CPPrn1 + Pud21*CPPrn2 + Pud31*CPPrn3)
      sdtphir2 = dtphir2 - bn*(Pud12*CPPrn1 + Pud22*CPPrn2 + Pud32*CPPrn3)
      sdtphir3 = dtphir3 - bn*(Pud13*CPPrn1 + Pud23*CPPrn2 + Pud33*CPPrn3)
      sdtphic1 = dtphic1 - bn*(Pud11*CPPcn1 + Pud21*CPPcn2 + Pud31*CPPcn3)
      sdtphic2 = dtphic2 - bn*(Pud12*CPPcn1 + Pud22*CPPcn2 + Pud32*CPPcn3)
      sdtphic3 = dtphic3 - bn*(Pud13*CPPcn1 + Pud23*CPPcn2 + Pud33*CPPcn3)
      
      sdt1H0 = dt1H0 - bn*(Pud11*CGGmn1 + Pud21*CGGmn2 + Pud31*CGGmn3)
      sdt2H0 = dt2H0 - bn*(Pud12*CGGmn1 + Pud22*CGGmn2 + Pud32*CGGmn3)
      sdt3H0 = dt3H0 - bn*(Pud13*CGGmn1 + Pud23*CGGmn2 + Pud33*CGGmn3)
            
      end if      !----- CPBC_type 0
      
                       
      !--------------------------------------------------------------
      !---------------------THE DIFFICOULT PART----------------------
      !--------------------------------------------------------------                  

      !------maximally dissipative for the L----------------
      if (CPBC_type .LE. 1) then

        sL00 = 0.0d0
        sL01 = 0.0d0
        sL02 = 0.0d0
        sL03 = 0.0d0
        sL11 = 0.0d0
        sL12 = 0.0d0
        sL13 = 0.0d0
        sL22 = 0.0d0
        sL23 = 0.0d0
        sL33 = 0.0d0
      
        sSm = 0.0d0
        sSr = 0.0d0
        sSc = 0.0d0
        sG  = 0.0d0
      
      
      !------Sommerfeld boundary conditions for L-----------            
      !------ but maximally dissipative for the gauge ------
      elseif (CPBC_type .EQ. 2) then      

        sL00 = (1.0d0/r - sigma2)*g00
        sL01 = (1.0d0/r - sigma2)*g01
        sL02 = (1.0d0/r - sigma2)*g02
        sL03 = (1.0d0/r - sigma2)*g03
        sL11 = (1.0d0/r - sigma2)*g11
        sL12 = (1.0d0/r - sigma2)*g12
        sL13 = (1.0d0/r - sigma2)*g13
        sL22 = (1.0d0/r - sigma2)*g22
        sL23 = (1.0d0/r - sigma2)*g23
        sL33 = (1.0d0/r - sigma2)*g33
      
        sSm = (1.0d0/r - sigma2)*phim
        sSn = (1.0d0/r - sigma2)*phin
        sSr = (1.0d0/r - sigma2)*phir
        sSc = (1.0d0/r - sigma2)*phic
        sG  = (1.0d0/r - sigma2)*H0


      !------Caltech CPBC for L-----------------------------
      else
      
        call my_exit('boundary problem : this CPBC_type does not exist')
                       
     
      end if 
                  
      
