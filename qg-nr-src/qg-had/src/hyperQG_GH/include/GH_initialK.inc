      
      !---compute the extrinsic curvature from the time derivative of the metric
      !---saved in K and the other quantities-----------------------------------
      
      !---invert the metric -----------------------
      h11 = g_e(1,1); h12 = g_e(1,2); h13 = g_e(1,3)
      h22 = g_e(2,2); h23 = g_e(2,3); h33 = g_e(3,3)
                                          
      huu11 = -h23**2 + h22*h33
      huu12 = h13*h23 - h12*h33
      huu13 = -(h13*h22) + h12*h23
      huu22 = -h13**2 + h11*h33
      huu23 = h12*h13 - h11*h23
      huu33 = -h12**2 + h11*h22
      deth = h11*huu11 + h12*huu12 + h13*huu13
      if(deth.lt.1e-1) deth = 1.e-1
      huu11 = huu11/deth
      huu12 = huu12/deth
      huu13 = huu13/deth
      huu22 = huu22/deth
      huu23 = huu23/deth
      huu33 = huu33/deth
                     
      b1 = g_e(0,1)*huu11 + g_e(0,2)*huu12 + g_e(0,3)*huu13
      b2 = g_e(0,1)*huu12 + g_e(0,2)*huu22 + g_e(0,3)*huu23
      b3 = g_e(0,1)*huu13 + g_e(0,2)*huu23 + g_e(0,3)*huu33
      alp2 = -g_e(0,0) + b1**2*h11 + b2**2*h22 + 2*(b1*(b2*h12 + b3*h13) &
     &     + b2*b3*h23) + b3**2*h33
      if(alp2.lt.0.01) alp2=0.01
      alp = sqrt(alp2)

      guu00 = -(1.0d0/alp2)
      guu01 = b1/alp2
      guu02 = b2/alp2
      guu03 = b3/alp2
      guu11 = -(b1**2/alp2) + huu11
      guu12 = -((b1*b2)/alp2) + huu12
      guu13 = -((b1*b3)/alp2) + huu13
      guu22 = -(b2**2/alp2) + huu22
      guu23 = -((b2*b3)/alp2) + huu23
      guu33 = -(b3**2/alp2) + huu33

      !--- read the derivatives of the metric--------
      D001 = K_e(0,1); D002 = K_e(0,2); D003 = K_e(0,3); D000 = K_e(0,0)      
      D011 = K_e(1,1); D012 = K_e(1,2); D013 = K_e(1,3)
      D022 = K_e(2,2); D023 = K_e(2,3); D033 = K_e(3,3)
      
      D100 = dg_e(1,0,0); D101 = dg_e(1,0,1); D102 = dg_e(1,0,2)
      D103 = dg_e(1,0,3); D111 = dg_e(1,1,1); D112 = dg_e(1,1,2)
      D113 = dg_e(1,1,3); D122 = dg_e(1,2,2); D123 = dg_e(1,2,3)
      D133 = dg_e(1,3,3)
       
      D200 = dg_e(2,0,0); D201 = dg_e(2,0,1); D202 = dg_e(2,0,2)
      D203 = dg_e(2,0,3); D211 = dg_e(2,1,1); D212 = dg_e(2,1,2)
      D213 = dg_e(2,1,3); D222 = dg_e(2,2,2); D223 = dg_e(2,2,3)
      D233 = dg_e(2,3,3)

      D300 = dg_e(3,0,0); D301 = dg_e(3,0,1); D302 = dg_e(3,0,2)
      D303 = dg_e(3,0,3); D311 = dg_e(3,1,1); D312 = dg_e(3,1,2)
      D313 = dg_e(3,1,3); D322 = dg_e(3,2,2); D323 = dg_e(3,2,3)
      D333 = dg_e(3,3,3)
           
      !---compute the D00a ------------------                   
      D000 = (guu11*(D011 - 2.*D101) + guu22*(D022 -&
        & 2.*D202) + 2.*(guu12*D012 + guu13*D013 + guu23*D023 -&
        & guu01*D100 - guu02*D200 - guu12*(D102 + D201) -&
        & guu03*D300 - guu13*(D103 + D301) - guu23*(D203 +&
        & D302)) + guu33*(D033 - 2.*D303))/guu00

      D001 = D100/2. + ((-(guu11*D111) + guu22*D122 + guu33*D133)/2.&
        & - guu12*D211 - guu22*D212 - guu13*D311 + guu23*(D123 - D213 - D312)&
        & - guu33*D313 - guu01*D011 - guu03*D013 + guu02*D102 -&
        & guu02*(D012 + D201) + guu03*(D103 - D301))/guu00

      D002 = D200/2. + (-(guu11*D112) - guu12*D122 - guu13*D123 -&
        & (guu22*D222)/2. + (guu11*D211 + guu33*D233)/2. + guu13*(D213 - D312)&
        & - guu23*D322 - guu33*D323 - guu02*D022 - guu01*(D012 +&
        & D102) + guu01*D201 + guu03*D203 - guu03*(D023 + D302))/guu00

      D003 = D300/2. + (-(guu11*D113) - guu13*D133 - guu12*(D123 +&
        & D213) - guu22*D223 - guu23*D233 + guu12*D312 + (guu11*D311 +&
        & guu22*D322)/2. - (guu33*D333)/2. - guu03*D033 - guu01*(D013 +&
        & D103) - guu02*(D023 + D203) + guu01*D301 + guu02*D302)/guu00

	
!-----we have Hs and D0ab, let us compute the evolved variable K-----------------------

      ! the difference between the evolved Ks and the extrinsic curvature--------------
      X_e(0,0) = D000 + D000 - ( b1*(D001 + D001) + b2*(D002 + D002) + b3*(D003 + D003) )
      X_e(0,1) = D001 + D100 - ( b1*(D011 + D101) + b2*(D012 + D102) + b3*(D013 + D103) )
      X_e(0,2) = D002 + D200 - ( b1*(D012 + D201) + b2*(D022 + D202) + b3*(D023 + D203) )
      X_e(0,3) = D003 + D300 - ( b1*(D013 + D301) + b2*(D023 + D302) + b3*(D033 + D303) )
      X_e(1,1) = D101 + D101 - ( b1*(D111 + D111) + b2*(D112 + D112) + b3*(D113 + D113) )
      X_e(1,2) = D102 + D201 - ( b1*(D112 + D211) + b2*(D122 + D212) + b3*(D123 + D213) )
      X_e(1,3) = D103 + D301 - ( b1*(D113 + D311) + b2*(D123 + D312) + b3*(D133 + D313) )
      X_e(2,2) = D202 + D202 - ( b1*(D212 + D212) + b2*(D222 + D222) + b3*(D223 + D223) )
      X_e(2,3) = D203 + D302 - ( b1*(D213 + D312) + b2*(D223 + D322) + b3*(D233 + D323) )
      X_e(3,3) = D303 + D303 - ( b1*(D313 + D313) + b2*(D323 + D323) + b3*(D333 + D333) )      
      ! the evolved Ks------------------------------------------------------------------                        
      K_e(0,0) = -( D000 - (b1*D100 + b2*D200 + b3*D300) )/alp
      K_e(0,1) = -( D001 - (b1*D101 + b2*D201 + b3*D301) )/alp
      K_e(0,2) = -( D002 - (b1*D102 + b2*D202 + b3*D302) )/alp
      K_e(0,3) = -( D003 - (b1*D103 + b2*D203 + b3*D303) )/alp                  
      K_e(1,1) = -( D011 - (b1*D111 + b2*D211 + b3*D311) )/alp
      K_e(1,2) = -( D012 - (b1*D112 + b2*D212 + b3*D312) )/alp
      K_e(1,3) = -( D013 - (b1*D113 + b2*D213 + b3*D313) )/alp
      K_e(2,2) = -( D022 - (b1*D122 + b2*D222 + b3*D322) )/alp
      K_e(2,3) = -( D023 - (b1*D123 + b2*D223 + b3*D323) )/alp
      K_e(3,3) = -( D033 - (b1*D133 + b2*D233 + b3*D333) )/alp
      ! the extrinsic curvature----------------------------------------------------------
      K_e(0,0) = 0.5d0*( K_e(0,0) + X_e(0,0)/alp )
      K_e(0,1) = 0.5d0*( K_e(0,1) + X_e(0,1)/alp )
      K_e(0,2) = 0.5d0*( K_e(0,2) + X_e(0,2)/alp )
      K_e(0,3) = 0.5d0*( K_e(0,3) + X_e(0,3)/alp )
      K_e(1,1) = 0.5d0*( K_e(1,1) + X_e(1,1)/alp )
      K_e(1,2) = 0.5d0*( K_e(1,2) + X_e(1,2)/alp )
      K_e(1,3) = 0.5d0*( K_e(1,3) + X_e(1,3)/alp )
      K_e(2,2) = 0.5d0*( K_e(2,2) + X_e(2,2)/alp )
      K_e(2,3) = 0.5d0*( K_e(2,3) + X_e(2,3)/alp )
      K_e(3,3) = 0.5d0*( K_e(3,3) + X_e(3,3)/alp )
