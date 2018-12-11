
!-----the first order constraints------------------------------------------

    order1_const_pt = &
          + (D100 - dx_g00)**2 + (D200 - dy_g00)**2 + (D300 - dz_g00)**2 & 
	  + (D101 - dx_g01)**2 + (D201 - dy_g01)**2 + (D301 - dz_g01)**2 & 
	  + (D102 - dx_g02)**2 + (D202 - dy_g02)**2 + (D302 - dz_g02)**2 & 
	  + (D103 - dx_g03)**2 + (D203 - dy_g03)**2 + (D303 - dz_g03)**2 & 
	  + (D111 - dx_g11)**2 + (D211 - dy_g11)**2 + (D311 - dz_g11)**2 & 
	  + (D112 - dx_g12)**2 + (D212 - dy_g12)**2 + (D312 - dz_g12)**2 & 
	  + (D113 - dx_g13)**2 + (D213 - dy_g13)**2 + (D313 - dz_g13)**2 & 
	  + (D122 - dx_g22)**2 + (D222 - dy_g22)**2 + (D322 - dz_g22)**2 & 
	  + (D123 - dx_g23)**2 + (D223 - dy_g23)**2 + (D323 - dz_g23)**2 & 
	  + (D133 - dx_g33)**2 + (D233 - dy_g33)**2 + (D333 - dz_g33)**2     
    order1_const_pt = dsqrt(order1_const_pt/30.0)
	  
    order2_const_pt = &
          + (dx_D200 - dy_D100)**2 + (dx_D300 - dz_D100)**2 &
          + (dy_D300 - dz_D200)**2 &
          + (dx_D201 - dy_D101)**2 + (dx_D301 - dz_D101)**2 &
          + (dy_D301 - dz_D201)**2 &
	  + (dx_D202 - dy_D102)**2 + (dx_D302 - dz_D102)**2 &
          + (dy_D302 - dz_D202)**2 &
	  + (dx_D203 - dy_D103)**2 + (dx_D303 - dz_D103)**2 &
          + (dy_D303 - dz_D203)**2 &
          + (dx_D211 - dy_D111)**2 + (dx_D311 - dz_D111)**2 &
          + (dy_D311 - dz_D211)**2 &
          + (dx_D212 - dy_D112)**2 + (dx_D312 - dz_D112)**2 &
          + (dy_D312 - dz_D212)**2 &
          + (dx_D213 - dy_D113)**2 + (dx_D313 - dz_D113)**2 &
          + (dy_D313 - dz_D213)**2 &
          + (dx_D222 - dy_D122)**2 + (dx_D322 - dz_D122)**2 &
          + (dy_D322 - dz_D222)**2 &
          + (dx_D223 - dy_D123)**2 + (dx_D323 - dz_D123)**2 &
          + (dy_D323 - dz_D223)**2 &
          + (dx_D233 - dy_D133)**2 + (dx_D333 - dz_D133)**2 &
          + (dy_D333 - dz_D233)**2 
    order2_const_pt = dsqrt(order2_const_pt/30.0)


!-----the 3+1 quantities----------------------------------------------------
       
       h11 = g11
       h12 = g12
       h13 = g13
       h22 = g22
       h23 = g23
       h33 = g33
                                                        
       huu11 = -h23**2 + h22*h33
       huu12 = h13*h23 - h12*h33
       huu13 = -(h13*h22) + h12*h23
       huu22 = -h13**2 + h11*h33
       huu23 = h12*h13 - h11*h23
       huu33 = -h12**2 + h11*h22
       deth = h11*huu11 + h12*huu12 + h13*huu13
       if(deth.lt.0.1) deth = 0.1
       huu11 = huu11/deth
       huu12 = huu12/deth
       huu13 = huu13/deth
       huu22 = huu22/deth
       huu23 = huu23/deth
       huu33 = huu33/deth
                     
       b1 = g01*huu11 + g02*huu12 + g03*huu13
       b2 = g01*huu12 + g02*huu22 + g03*huu23
       b3 = g01*huu13 + g02*huu23 + g03*huu33
       alp2 = -g00 + b1**2*h11 + b2**2*h22 + 2*(b1*(b2*h12 + b3*h13) &
     &      + b2*b3*h23) + b3**2*h33
       if(alp2.lt.1e-4) alp2=1e-4
       alp = dsqrt(alp2)

       if (alp2.le.0) then
          !
!          write(0,*) myid,']  analysis: alp2 <=0, Old alp2: ',alp2
!          write(0,*) myid,']  analysis: g00      = ',g00
!          write(0,*) myid,']  analysis: dx       = ',par(p_dx)
          !
          alp2 = 0.001d0 * par(p_dx)        
          alp  = dsqrt(alp2)                
          !
!          write(0,*) myid,']  analysis: New alp2 = ',alp2
!          write(0,*) myid,']  analysis: New alp  = ',alp
          !
       end if
      
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
       
       eK11 = -(-alp*K11 - D101 - D101 &
      &    +b1*(D111 + D111) + b2*(D112 + D112) + b3*(D113 + D113))/(2.0*alp)
       eK12 = -(-alp*K12 - D102 - D201 &
      &    +b1*(D112 + D211) + b2*(D122 + D212) + b3*(D123 + D213))/(2.0*alp)
       eK13 = -(-alp*K13 - D103 - D301 &
      &    +b1*(D113 + D311) + b2*(D123 + D312) + b3*(D133 + D313))/(2.0*alp)
       eK22 = -(-alp*K22 - D202 - D202 &
      &    +b1*(D212 + D212) + b2*(D222 + D222) + b3*(D223 + D223))/(2.0*alp)
       eK23 = -(-alp*K23 - D203 - D302 &
      &    +b1*(D213 + D312) + b2*(D223 + D322) + b3*(D233 + D323))/(2.0*alp)
       eK33 = -(-alp*K33 - D303 - D303 &
      &    +b1*(D313 + D313) + b2*(D323 + D323) + b3*(D333 + D333))/(2.0*alp)


!-----a lot of stuff for the energy constraint---------------------------
      
       Dddu111 = huu11*D111 + huu12*D112 + huu13*D113
       Dddu112 = huu12*D111 + huu22*D112 + huu23*D113
       Dddu113 = huu13*D111 + huu23*D112 + huu33*D113
       Dddu121 = huu11*D112 + huu12*D122 + huu13*D123
       Dddu122 = huu12*D112 + huu22*D122 + huu23*D123
       Dddu123 = huu13*D112 + huu23*D122 + huu33*D123
       Dddu131 = huu11*D113 + huu12*D123 + huu13*D133
       Dddu132 = huu12*D113 + huu22*D123 + huu23*D133
       Dddu133 = huu13*D113 + huu23*D123 + huu33*D133
       Dddu211 = huu11*D211 + huu12*D212 + huu13*D213
       Dddu212 = huu12*D211 + huu22*D212 + huu23*D213
       Dddu213 = huu13*D211 + huu23*D212 + huu33*D213
       Dddu221 = huu11*D212 + huu12*D222 + huu13*D223
       Dddu222 = huu12*D212 + huu22*D222 + huu23*D223
       Dddu223 = huu13*D212 + huu23*D222 + huu33*D223
       Dddu231 = huu11*D213 + huu12*D223 + huu13*D233
       Dddu232 = huu12*D213 + huu22*D223 + huu23*D233
       Dddu233 = huu13*D213 + huu23*D223 + huu33*D233
       Dddu311 = huu11*D311 + huu12*D312 + huu13*D313
       Dddu312 = huu12*D311 + huu22*D312 + huu23*D313
       Dddu313 = huu13*D311 + huu23*D312 + huu33*D313
       Dddu321 = huu11*D312 + huu12*D322 + huu13*D323
       Dddu322 = huu12*D312 + huu22*D322 + huu23*D323
       Dddu323 = huu13*D312 + huu23*D322 + huu33*D323
       Dddu331 = huu11*D313 + huu12*D323 + huu13*D333
       Dddu332 = huu12*D313 + huu22*D323 + huu23*D333
       Dddu333 = huu13*D313 + huu23*D323 + huu33*D333
       
       Dduu111 = huu11**2*D111 + huu12**2*D122 + 2*(huu11*(huu12*D112 +&
     & huu13*D113) + huu12*huu13*D123) + huu13**2*D133
       Dduu112 = huu12**2*D112 + huu11* (huu12*D111 + huu22*D112 +&
     & huu23*D113) + huu12*(huu13*D113 + huu22*D122 + huu23*D123) +&
     & huu13* (huu22*D123 + huu23*D133)
       Dduu113 = huu13**2*D113 + huu11* (huu23*D112 + huu33*D113) +&
     & huu12*(huu23*D122 + huu33*D123) + huu13*(huu11*D111 + huu12*D112&
     & + huu23*D123 + huu33*D133)
       Dduu122 = huu12**2*D111 + huu22**2*D122 + 2*(huu12*(huu22*D112 +&
     & huu23*D113) + huu22*huu23*D123) + huu23**2*D133
       Dduu123 = huu13*(huu12*D111 + huu22*D112) + huu23**2*D123 +&
     & huu33* (huu12*D113 + huu22*D123) + huu23*(huu12*D112 + huu13*D113&
     & + huu22*D122 + huu33*D133)
       Dduu133 = huu13**2*D111 + huu23**2*D122 + 2*(huu13*(huu23*D112 +&
     & huu33*D113) + huu23*huu33*D123) + huu33**2*D133
       Dduu211 = huu11**2*D211 + huu12**2*D222 + 2*(huu11*(huu12*D212 +&
     & huu13*D213) + huu12*huu13*D223) + huu13**2*D233
       Dduu212 = huu12**2*D212 + huu11* (huu12*D211 + huu22*D212 +&
     & huu23*D213) + huu12*(huu13*D213 + huu22*D222 + huu23*D223) +&
     & huu13* (huu22*D223 + huu23*D233)
       Dduu213 = huu13**2*D213 + huu11* (huu23*D212 + huu33*D213) +&
     & huu12*(huu23*D222 + huu33*D223) + huu13*(huu11*D211 + huu12*D212&
     & + huu23*D223 + huu33*D233)
       Dduu222 = huu12**2*D211 + huu22**2*D222 + 2*(huu12*(huu22*D212 +&
     & huu23*D213) + huu22*huu23*D223) + huu23**2*D233
       Dduu223 = huu13*(huu12*D211 + huu22*D212) + huu23**2*D223 +&
     & huu33* (huu12*D213 + huu22*D223) + huu23*(huu12*D212 + huu13*D213&
     & + huu22*D222 + huu33*D233)
       Dduu233 = huu13**2*D211 + huu23**2*D222 + 2*(huu13*(huu23*D212 +&
     & huu33*D213) + huu23*huu33*D223) + huu33**2*D233
       Dduu311 = huu11**2*D311 + huu12**2*D322 + 2*(huu11*(huu12*D312 +&
     & huu13*D313) + huu12*huu13*D323) + huu13**2*D333
       Dduu312 = huu12**2*D312 + huu11* (huu12*D311 + huu22*D312 +&
     & huu23*D313) + huu12*(huu13*D313 + huu22*D322 + huu23*D323) +&
     & huu13* (huu22*D323 + huu23*D333)
       Dduu313 = huu13**2*D313 + huu11* (huu23*D312 + huu33*D313) +&
     & huu12*(huu23*D322 + huu33*D323) + huu13*(huu11*D311 + huu12*D312&
     & + huu23*D323 + huu33*D333)
       Dduu322 = huu12**2*D311 + huu22**2*D322 + 2*(huu12*(huu22*D312 +&
     & huu23*D313) + huu22*huu23*D323) + huu23**2*D333
       Dduu323 = huu13*(huu12*D311 + huu22*D312) + huu23**2*D323 +&
     & huu33* (huu12*D313 + huu22*D323) + huu23*(huu12*D312 + huu13*D313&
     & + huu22*D322 + huu33*D333)
       Dduu333 = huu13**2*D311 + huu23**2*D322 + 2*(huu13*(huu23*D312 +&
     & huu33*D313) + huu23*huu33*D323) + huu33**2*D333
       Dudd111 = huu11*D111 + huu12*D211 + huu13*D311
       Dudd112 = huu11*D112 + huu12*D212 + huu13*D312
       Dudd113 = huu11*D113 + huu12*D213 + huu13*D313
       Dudd122 = huu11*D122 + huu12*D222 + huu13*D322
       Dudd123 = huu11*D123 + huu12*D223 + huu13*D323
       Dudd133 = huu11*D133 + huu12*D233 + huu13*D333
       Dudd211 = huu12*D111 + huu22*D211 + huu23*D311
       Dudd212 = huu12*D112 + huu22*D212 + huu23*D312
       Dudd213 = huu12*D113 + huu22*D213 + huu23*D313
       Dudd222 = huu12*D122 + huu22*D222 + huu23*D322
       Dudd223 = huu12*D123 + huu22*D223 + huu23*D323
       Dudd233 = huu12*D133 + huu22*D233 + huu23*D333
       Dudd311 = huu13*D111 + huu23*D211 + huu33*D311
       Dudd312 = huu13*D112 + huu23*D212 + huu33*D312
       Dudd313 = huu13*D113 + huu23*D213 + huu33*D313
       Dudd322 = huu13*D122 + huu23*D222 + huu33*D322
       Dudd323 = huu13*D123 + huu23*D223 + huu33*D323
       Dudd333 = huu13*D133 + huu23*D233 + huu33*D333
       trD1 = Dddu111 + Dddu122 + Dddu133
       trD2 = Dddu211 + Dddu222 + Dddu233
       trD3 = Dddu311 + Dddu322 + Dddu333
       trE1 = Dudd111 + Dudd212 + Dudd313
       trE2 = Dudd112 + Dudd222 + Dudd323
       trE3 = Dudd113 + Dudd223 + Dudd333
       trDu1 = huu11*trD1 + huu12*trD2 + huu13*trD3
       trDu2 = huu12*trD1 + huu22*trD2 + huu23*trD3
       trDu3 = huu13*trD1 + huu23*trD2 + huu33*trD3
       trEu1 = huu11*trE1 + huu12*trE2 + huu13*trE3
       trEu2 = huu12*trE1 + huu22*trE2 + huu23*trE3
       trEu3 = huu13*trE1 + huu23*trE2 + huu33*trE3

       !with the first index up                                             
       G111 = Dddu111 - Dudd111/2.
       G112 = (Dddu121 + Dddu211 - Dudd112)/2.
       G113 = (Dddu131 + Dddu311 - Dudd113)/2.
       G122 = Dddu221 - Dudd122/2.
       G123 = (Dddu231 + Dddu321 - Dudd123)/2.
       G133 = Dddu331 - Dudd133/2.
       G211 = Dddu112 - Dudd211/2.
       G212 = (Dddu122 + Dddu212 - Dudd212)/2.
       G213 = (Dddu132 + Dddu312 - Dudd213)/2.
       G222 = Dddu222 - Dudd222/2.
       G223 = (Dddu232 + Dddu322 - Dudd223)/2.
       G233 = Dddu332 - Dudd233/2.
       G311 = Dddu113 - Dudd311/2.
       G312 = (Dddu123 + Dddu213 - Dudd312)/2.
       G313 = (Dddu133 + Dddu313 - Dudd313)/2.
       G322 = Dddu223 - Dudd322/2.
       G323 = (Dddu233 + Dddu323 - Dudd323)/2.
       G333 = Dddu333 - Dudd333/2.
                    
       eKud11 = huu11*eK11 + huu12*eK12 + huu13*eK13
       eKud12 = huu11*eK12 + huu12*eK22 + huu13*eK23
       eKud13 = huu11*eK13 + huu12*eK23 + huu13*eK33
       eKud21 = huu12*eK11 + huu22*eK12 + huu23*eK13
       eKud22 = huu12*eK12 + huu22*eK22 + huu23*eK23
       eKud23 = huu12*eK13 + huu22*eK23 + huu23*eK33
       eKud31 = huu13*eK11 + huu23*eK12 + huu33*eK13
       eKud32 = huu13*eK12 + huu23*eK22 + huu33*eK23
       eKud33 = huu13*eK13 + huu23*eK23 + huu33*eK33
       etrK = eKud11 + eKud22 + eKud33       
       
       energy_const = (dx_D122 - dx_D212 - dy_D112 + dy_D211)*huu12**2 -&
     & dz_D113*huu13**2 + (dx_D133 - dx_D313 + dz_D311)*huu13**2 +&
     & ((-2*dx_D123 + dz_D112)*huu11 - dz_D211*huu11 + (-dx_D223 +&
     & 2*dy_D213 + dz_D122)*huu12 - dz_D212*huu12 + (-dx_D323 + dy_D133&
     & + 2*dz_D312)*huu13)*huu23 - dz_D223*huu23**2 + (dy_D233 -&
     & dy_D323 + dz_D322)*huu23**2 + huu13*(-(dx_D312*huu12) -&
     & 2*dy_D213*huu22 + (dx_D223 - dx_D322 + dy_D123 + dy_D312 +&
     & dz_D212)*huu22 + (dx_D233 - dy_D313)*huu23 - dz_D123*huu23) -&
     & huu13*(dz_D112*huu12 + dz_D122*huu22 + dz_D213*huu23) -&
     & dz_D322*huu22*huu33 + ((-dx_D133 + dz_D113)*huu11 + (dy_D313 -&
     & 2*dz_D312)*huu12 + (-dy_D233 + dy_D323 + dz_D223)*huu22)*huu33 +&
     & huu11*((-dx_D122 + dx_D212 + dy_D112)*huu22 + (dx_D213 + dx_D312&
     & + dy_D113)*huu23 + dx_D313*huu33) + huu12*(-(dy_D113*huu13) +&
     & (2*dx_D123 - dx_D213 + dy_D311 + dz_D211)* huu13 + (dx_D322 -&
     & dy_D123)*huu23 - dy_D312*huu23 - dy_D133*huu33 + (-dx_D233 +&
     & dx_D323 + dz_D123 + dz_D213)*huu33) - huu11*(dy_D211*huu22 +&
     & dy_D311*huu23 + dz_D311*huu33) - eKud11**2 - eKud22**2 -&
     & 2*(eKud12*eKud21 + eKud13*eKud31 + eKud23*eKud32) - eKud33**2 -&
     & trE1*trEu1 - trE2*trEu2 - trE3*trEu3 + etrK**2 + trEu1*trD1 -&
     & (trD1*trDu1)/4. + trD2*(trEu2 - trDu2/4.) + trD3*(trEu3 -&
     & trDu3/4.) - Dddu111*Dduu111 - Dddu121*Dduu112 - Dddu131*Dduu113 -&
     & Dddu221*Dduu122 - Dddu231*Dduu123 - Dddu331*Dduu133 -&
     & Dddu112*Dduu211 - Dddu122*Dduu212 - Dddu132*Dduu213 -&
     & Dddu222*Dduu222 - Dddu232*Dduu223 - Dddu332*Dduu233 -&
     & Dddu113*Dduu311 - Dddu123*Dduu312 - Dddu133*Dduu313 -&
     & Dddu223*Dduu322 - Dddu233*Dduu323 - Dddu333*Dduu333 +&
     & (G111*Dduu111 + G122*Dduu122 + G133*Dduu133 + G211*Dduu211 +&
     & G222*Dduu222 + G233*Dduu233 + G311*Dduu311 + G322*Dduu322 +&
     & G333*Dduu333)/ 2. + Dduu111*Dudd111 + Dduu112*(G112 - Dddu211 +&
     & 2*Dudd112) + Dduu113*(G113 - Dddu311 + 2*Dudd113) +&
     & Dduu122*Dudd122 + Dduu123* (G123 - Dddu321 + 2*Dudd123) +&
     & Dduu133*Dudd133 + Dduu211*Dudd211 + Dduu212*(G212 - Dddu212 +&
     & 2*Dudd212) + Dduu213*(G213 - Dddu312 + 2*Dudd213) +&
     & Dduu222*Dudd222 + Dduu223* (G223 - Dddu322 + 2*Dudd223) +&
     & Dduu233*Dudd233 + Dduu311*Dudd311 + Dduu312*(G312 - Dddu213 +&
     & 2*Dudd312) + Dduu313*(G313 - Dddu313 + 2*Dudd313) +&
     & Dduu322*Dudd322 + Dduu323* (G323 - Dddu323 + 2*Dudd323) +&
     & Dduu333*Dudd333
     
!-----the energy density with the sqrt(h) from the volume factor----------     
     
       tu0 = 1.0/alp; tu1 = -b1/alp; tu2 = -b2/alp; tu3 = -b3/alp 
       tau = tu0*tu0*T00 + tu1*tu1*T11 + tu2*tu2*T22 + tu3*tu3*T33 & 
     &      + 2.d0*(tu0*tu1*T01 + tu0*tu2*T02 + tu0*tu3*T03 &
     &            + tu1*tu2*T12 + tu1*tu3*T13 + tu2*tu3*T23)
       		          
       energy_const = energy_const - 16.0*pi*tau
       energy_density = tau*sqrt(deth)
       
!-----the Noether density with the sqrt(h) from the volume factor----------
                   
!      noether_density = 2.0d0*deth*(phir*pic - phic*pir)
      noether_density = sqrt(deth)*(phir*pic - phic*pir)
       
            
!-----the Zs quantities----------------------------------------------------
                  
       D000 = -alp*K00 + b1*D100 + b2*D200 + b3*D300 
       D001 = -alp*K01 + b1*D101 + b2*D201 + b3*D301
       D002 = -alp*K02 + b1*D102 + b2*D202 + b3*D302 
       D003 = -alp*K03 + b1*D103 + b2*D203 + b3*D303 
       D011 = -alp*K11 + b1*D111 + b2*D211 + b3*D311 
       D012 = -alp*K12 + b1*D112 + b2*D212 + b3*D312 
       D013 = -alp*K13 + b1*D113 + b2*D213 + b3*D313 
       D022 = -alp*K22 + b1*D122 + b2*D222 + b3*D322 
       D023 = -alp*K23 + b1*D123 + b2*D223 + b3*D323
       D033 = -alp*K33 + b1*D133 + b2*D233 + b3*D333
       
       Gd000 = D000/2.
       Gd001 = D100/2.
       Gd002 = D200/2.
       Gd003 = D300/2.
       Gd011 = -D011/2. + D101
       Gd012 = (-D012 + D102 + D201)/2.
       Gd013 = (-D013 + D103 + D301)/2.
       Gd022 = -D022/2. + D202
       Gd023 = (-D023 + D203 + D302)/2.
       Gd033 = -D033/2. + D303
       Gd100 = D001 - D100/2.
       Gd101 = D011/2.
       Gd102 = (D012 - D102 + D201)/2.
       Gd103 = (D013 - D103 + D301)/2.
       Gd111 = D111/2.
       Gd112 = D211/2.
       Gd113 = D311/2.
       Gd122 = -D122/2. + D212
       Gd123 = (-D123 + D213 + D312)/2.
       Gd133 = -D133/2. + D313
       Gd200 = D002 - D200/2.
       Gd201 = (D012 + D102 - D201)/2.
       Gd202 = D022/2.
       Gd203 = (D023 - D203 + D302)/2.
       Gd211 = D112 - D211/2.
       Gd212 = D122/2.
       Gd213 = (D123 - D213 + D312)/2.
       Gd222 = D222/2.
       Gd223 = D322/2.
       Gd233 = -D233/2. + D323
       Gd300 = D003 - D300/2.
       Gd301 = (D013 + D103 - D301)/2.
       Gd302 = (D023 + D203 - D302)/2.
       Gd303 = D033/2.
       Gd311 = D113 - D311/2.
       Gd312 = (D123 + D213 - D312)/2.
       Gd313 = D133/2.
       Gd322 = D223 - D322/2.
       Gd323 = D233/2.
       Gd333 = D333/2.

       trG0 = guu00*Gd000 + guu11*Gd011 + guu22*Gd022 + guu33*Gd033 &
      &   +2.0*(guu01*Gd001 + guu02*Gd002 + guu03*Gd003 +&
      &         guu12*Gd012 + guu13*Gd013 + guu23*Gd023) 
      
       trG1 = guu00*Gd100 + guu11*Gd111 + guu22*Gd122 + guu33*Gd133 &
      &   +2.0*(guu01*Gd101 + guu02*Gd102 + guu03*Gd103 +&
      &         guu12*Gd112 + guu13*Gd113 + guu23*Gd123) 

       trG2 = guu00*Gd200 + guu11*Gd211 + guu22*Gd222 + guu33*Gd233 &
      &   +2.0*(guu01*Gd201 + guu02*Gd202 + guu03*Gd203 +&
      &         guu12*Gd212 + guu13*Gd213 + guu23*Gd223) 
 
       trG3 = guu00*Gd300 + guu11*Gd311 + guu22*Gd322 + guu33*Gd333 &
      &   +2.0*(guu01*Gd301 + guu02*Gd302 + guu03*Gd303 +&
      &         guu12*Gd312 + guu13*Gd313 + guu23*Gd323)       
                        
       Z0 = (-H0e - trG0)/2.
       Z1 = (-H1e - trG1)/2.
       Z2 = (-H2e - trG2)/2.
       Z3 = (-H3e - trG3)/2.
             
       physZ_const_pt = Z0**2 + Z1**2 + Z2**2 + Z3**2
       physZ_const_pt = dsqrt(physZ_const_pt/4.0) 
       allconst_pt = order1_const_pt + order2_const_pt + physZ_const_pt   
       
!-----computing the radial outgoing speed -----------------------------
!       nr1 = x; nr2 = y; nr3 = z        
!       nrnorm = dsqrt(huu11*nr1*nr1 + huu22*nr2*nr2 + huu33*nr3*nr3 +& 
!                2.0d0*(huu12*nr1*nr2 + huu13*nr1*nr3 + huu23*nr2*nr3))
!       if (nrnorm .LT. 1E-10) nrnorm = 1E-10
!       nr1 = nr1/nrnorm
!       nr2 = nr2/nrnorm
!       nr3 = nr3/nrnorm      
              
!       rad_speed_pt = - (nr1*b1 + nr2*b2 + nr3*b3) + alp
       
!-----computing the radial expansion -----------------------------       
!       nru1 = huu11*nr1 + huu12*nr2 + huu13*nr3
!       nru2 = huu12*nr1 + huu22*nr2 + huu23*nr3
!       nru3 = huu13*nr1 + huu23*nr2 + huu33*nr3
                            
!       Gndd11 = nr1*G111 + nr2*G211 + nr3*G311
!       Gndd12 = nr1*G112 + nr2*G212 + nr3*G312
!       Gndd13 = nr1*G113 + nr2*G213 + nr3*G313
!       Gndd22 = nr1*G122 + nr2*G222 + nr3*G322
!       Gndd23 = nr1*G123 + nr2*G223 + nr3*G323
!       Gndd33 = nr1*G133 + nr2*G233 + nr3*G333
              
!       rad_exp_pt = ( (huu11 - nru1*nru1) + (huu22 - nru2*nru2)&
!      &      +        (huu33 - nru3*nru3) )/nrnorm &
!      &      - ( huu11*Gndd11 + huu22*Gndd22 + huu33*Gndd33 &
!      &      + 2.0*(huu12*Gndd12 + huu13*Gndd13 + huu23*Gndd23) ) &
!      &      + nru1*nru1*(eK11  + 0.5*(nru1*D111 + nru2*D211 + nru3*D311)) &
!      &      + nru2*nru2*(eK22  + 0.5*(nru1*D122 + nru2*D222 + nru3*D322)) & 
!      &      + nru3*nru3*(eK33  + 0.5*(nru1*D133 + nru2*D233 + nru3*D333)) &       
!      &  + 2.0*nru1*nru2*(eK12  + 0.5*(nru1*D112 + nru2*D212 + nru3*D312)) & 
!      &  + 2.0*nru1*nru3*(eK13  + 0.5*(nru1*D113 + nru2*D213 + nru3*D313)) & 
!      &  + 2.0*nru2*nru3*(eK23  + 0.5*(nru1*D123 + nru2*D223 + nru3*D323)) &
!      &      - etrK
      
      !-- compute the angular momentum from a killing----------
	
       xiu0 = 0.0d0
       xiu1 = -y	
       xiu2 = x 
       xiu3 = 0.0d0
	 
       xinorm = g00*xiu0*xiu0 + g11*xiu1*xiu1 + g22*xiu2*xiu2 + g33*xiu3*xiu3 &
     &      + 2.0d0*( g01*xiu0*xiu1 + g02*xiu0*xiu2 + g03*xiu0*xiu3 &
     &              + g12*xiu1*xiu2 + g13*xiu1*xiu3 + g23*xiu2*xiu3)

       if (dabs(xinorm) .LT. 1E-8) xinorm=1E-8
          
       xiu0 = xiu0/xinorm; xiu1 = xiu1/xinorm; xiu2 = xiu2/xinorm; xiu3 = xiu3/xinorm
      
       kang_momz = - sqrt(deth)*&
     &	   (xiu0*tu0*T00 + xiu0*tu1*T01 + xiu0*tu2*T02 + xiu0*tu3*T03 + &
     &      xiu1*tu0*T01 + xiu1*tu1*T11 + xiu1*tu2*T12 + xiu1*tu3*T13 + &
     &      xiu2*tu0*T02 + xiu2*tu1*T12 + xiu2*tu2*T22 + xiu2*tu3*T23 + &
     &      xiu3*tu0*T03 + xiu3*tu1*T13 + xiu3*tu2*T23 + xiu3*tu3*T33 )
      
!--------------------------------------------------------------------------------
!-----computing the frame, the Riemann and the Psi4 -----------------------------       
!--------------------------------------------------------------------------------
       if ((psi4_analysis .EQ. 1) .OR. (psi4_analysis .EQ. 3)) then

         !---------the 4-metric----------------------------------
  	 g(1,1) = g11; g(1,2) = g12; g(1,3) = g13; g(1,4) = g01
         g(2,1) = g12; g(2,2) = g22; g(2,3) = g23; g(2,4) = g02
         g(3,1) = g13; g(3,2) = g23; g(3,3) = g33; g(3,4) = g03		
         g(4,1) = g01; g(4,2) = g02; g(4,3) = g03; g(4,4) = g00
	
	 !---------the inverse of the 4-metric--------------------
         g_uu(1,1) = guu11; g_uu(1,2) = guu12; g_uu(1,3) = guu13; g_uu(1,4) = guu01
         g_uu(2,1) = guu12; g_uu(2,2) = guu22; g_uu(2,3) = guu23; g_uu(2,4) = guu02
         g_uu(3,1) = guu13; g_uu(3,2) = guu23; g_uu(3,3) = guu33; g_uu(3,4) = guu03
         g_uu(4,1) = guu01; g_uu(4,2) = guu02; g_uu(4,3) = guu03; g_uu(4,4) = guu00

	 !---------the space derivatives of the 4-metric---------
	 dg(1,1,1) = u_pt(H_D1G11); dg(1,1,2) = u_pt(H_D1G12); dg(1,1,3) = u_pt(H_D1G13); dg(1,1,4) = u_pt(H_D1G01)
         dg(1,2,1) = u_pt(H_D1G12); dg(1,2,2) = u_pt(H_D1G22); dg(1,2,3) = u_pt(H_D1G23); dg(1,2,4) = u_pt(H_D1G02)
         dg(1,3,1) = u_pt(H_D1G13); dg(1,3,2) = u_pt(H_D1G23); dg(1,3,3) = u_pt(H_D1G33); dg(1,3,4) = u_pt(H_D1G03)
         dg(1,4,1) = u_pt(H_D1G01); dg(1,4,2) = u_pt(H_D1G02); dg(1,4,3) = u_pt(H_D1G03); dg(1,4,4) = u_pt(H_D1G00)
	
	 dg(2,1,1) = u_pt(H_D2G11); dg(2,1,2) = u_pt(H_D2G12); dg(2,1,3) = u_pt(H_D2G13); dg(2,1,4) = u_pt(H_D2G01)
         dg(2,2,1) = u_pt(H_D2G12); dg(2,2,2) = u_pt(H_D2G22); dg(2,2,3) = u_pt(H_D2G23); dg(2,2,4) = u_pt(H_D2G02)
         dg(2,3,1) = u_pt(H_D2G13); dg(2,3,2) = u_pt(H_D2G23); dg(2,3,3) = u_pt(H_D2G33); dg(2,3,4) = u_pt(H_D2G03)
         dg(2,4,1) = u_pt(H_D2G01); dg(2,4,2) = u_pt(H_D2G02); dg(2,4,3) = u_pt(H_D2G03); dg(2,4,4) = u_pt(H_D2G00)
        
	 dg(3,1,1) = u_pt(H_D3G11); dg(3,1,2) = u_pt(H_D3G12); dg(3,1,3) = u_pt(H_D3G13); dg(3,1,4) = u_pt(H_D3G01)
         dg(3,2,1) = u_pt(H_D3G12); dg(3,2,2) = u_pt(H_D3G22); dg(3,2,3) = u_pt(H_D3G23); dg(3,2,4) = u_pt(H_D3G02)
         dg(3,3,1) = u_pt(H_D3G13); dg(3,3,2) = u_pt(H_D3G23); dg(3,3,3) = u_pt(H_D3G33); dg(3,3,4) = u_pt(H_D3G03)
         dg(3,4,1) = u_pt(H_D3G01); dg(3,4,2) = u_pt(H_D3G02); dg(3,4,3) = u_pt(H_D3G03); dg(3,4,4) = u_pt(H_D3G00)

	 !---------the "time" derivatives of the 4-metric---------
	 dg(4,1,1) = D011; dg(4,1,2) = D012; dg(4,1,3) = D013; dg(4,1,4) = D001
	 dg(4,2,1) = D012; dg(4,2,2) = D022; dg(4,2,3) = D023; dg(4,2,4) = D002
	 dg(4,3,1) = D013; dg(4,3,2) = D023; dg(4,3,3) = D033; dg(4,3,4) = D003
	 dg(4,4,1) = D001; dg(4,4,2) = D002; dg(4,4,3) = D003; dg(4,4,4) = D000

	 !---------the space derivates of the 4-metric------------
	 D_dg = 0.0d0
         !---x direction -------------	
	 D_dg(1,1,1,1) = dxu_pt(H_D1G11); D_dg(1,1,1,2) = dxu_pt(H_D1G12)
	 D_dg(1,1,1,3) = dxu_pt(H_D1G13); D_dg(1,1,1,4) = dxu_pt(H_D1G01)
         D_dg(1,1,2,1) = dxu_pt(H_D1G12); D_dg(1,1,2,2) = dxu_pt(H_D1G22)
	 D_dg(1,1,2,3) = dxu_pt(H_D1G23); D_dg(1,1,2,4) = dxu_pt(H_D1G02)
         D_dg(1,1,3,1) = dxu_pt(H_D1G13); D_dg(1,1,3,2) = dxu_pt(H_D1G23)
	 D_dg(1,1,3,3) = dxu_pt(H_D1G33); D_dg(1,1,3,4) = dxu_pt(H_D1G03)		 
         D_dg(1,1,4,1) = dxu_pt(H_D1G01); D_dg(1,1,4,2) = dxu_pt(H_D1G02)
	 D_dg(1,1,4,3) = dxu_pt(H_D1G03); D_dg(1,1,4,4) = dxu_pt(H_D1G00)

         D_dg(1,2,1,1) = dxu_pt(H_D2G11); D_dg(1,2,1,2) = dxu_pt(H_D2G12)
	 D_dg(1,2,1,3) = dxu_pt(H_D2G13); D_dg(1,2,1,4) = dxu_pt(H_D2G01)
         D_dg(1,2,2,1) = dxu_pt(H_D2G12); D_dg(1,2,2,2) = dxu_pt(H_D2G22)
	 D_dg(1,2,2,3) = dxu_pt(H_D2G23); D_dg(1,2,2,4) = dxu_pt(H_D2G02)
         D_dg(1,2,3,1) = dxu_pt(H_D2G13); D_dg(1,2,3,2) = dxu_pt(H_D2G23)
	 D_dg(1,2,3,3) = dxu_pt(H_D2G33); D_dg(1,2,3,4) = dxu_pt(H_D2G03)		 
         D_dg(1,2,4,1) = dxu_pt(H_D2G01); D_dg(1,2,4,2) = dxu_pt(H_D2G02)
	 D_dg(1,2,4,3) = dxu_pt(H_D2G03); D_dg(1,2,4,4) = dxu_pt(H_D2G00)

         D_dg(1,3,1,1) = dxu_pt(H_D3G11); D_dg(1,3,1,2) = dxu_pt(H_D3G12)
  	 D_dg(1,3,1,3) = dxu_pt(H_D3G13); D_dg(1,3,1,4) = dxu_pt(H_D3G01)
         D_dg(1,3,2,1) = dxu_pt(H_D3G12); D_dg(1,3,2,2) = dxu_pt(H_D3G22)
	 D_dg(1,3,2,3) = dxu_pt(H_D3G23); D_dg(1,3,2,4) = dxu_pt(H_D3G02)
         D_dg(1,3,3,1) = dxu_pt(H_D3G13); D_dg(1,3,3,2) = dxu_pt(H_D3G23)
	 D_dg(1,3,3,3) = dxu_pt(H_D3G33); D_dg(1,3,3,4) = dxu_pt(H_D3G03)		 
         D_dg(1,3,4,1) = dxu_pt(H_D3G01); D_dg(1,3,4,2) = dxu_pt(H_D3G02)
	 D_dg(1,3,4,3) = dxu_pt(H_D3G03); D_dg(1,3,4,4) = dxu_pt(H_D3G00)
	
         !---y direction -------------	
	 D_dg(2,1,1,1) = dyu_pt(H_D1G11); D_dg(2,1,1,2) = dyu_pt(H_D1G12)
	 D_dg(2,1,1,3) = dyu_pt(H_D1G13); D_dg(2,1,1,4) = dyu_pt(H_D1G01)
         D_dg(2,1,2,1) = dyu_pt(H_D1G12); D_dg(2,1,2,2) = dyu_pt(H_D1G22)
	 D_dg(2,1,2,3) = dyu_pt(H_D1G23); D_dg(2,1,2,4) = dyu_pt(H_D1G02)
         D_dg(2,1,3,1) = dyu_pt(H_D1G13); D_dg(2,1,3,2) = dyu_pt(H_D1G23)
	 D_dg(2,1,3,3) = dyu_pt(H_D1G33); D_dg(2,1,3,4) = dyu_pt(H_D1G03)		 
         D_dg(2,1,4,1) = dyu_pt(H_D1G01); D_dg(2,1,4,2) = dyu_pt(H_D1G02)
	 D_dg(2,1,4,3) = dyu_pt(H_D1G03); D_dg(2,1,4,4) = dyu_pt(H_D1G00)

         D_dg(2,2,1,1) = dyu_pt(H_D2G11); D_dg(2,2,1,2) = dyu_pt(H_D2G12)
	 D_dg(2,2,1,3) = dyu_pt(H_D2G13); D_dg(2,2,1,4) = dyu_pt(H_D2G01)
         D_dg(2,2,2,1) = dyu_pt(H_D2G12); D_dg(2,2,2,2) = dyu_pt(H_D2G22)
	 D_dg(2,2,2,3) = dyu_pt(H_D2G23); D_dg(2,2,2,4) = dyu_pt(H_D2G02)
         D_dg(2,2,3,1) = dyu_pt(H_D2G13); D_dg(2,2,3,2) = dyu_pt(H_D2G23)
	 D_dg(2,2,3,3) = dyu_pt(H_D2G33); D_dg(2,2,3,4) = dyu_pt(H_D2G03)		 
         D_dg(2,2,4,1) = dyu_pt(H_D2G01); D_dg(2,2,4,2) = dyu_pt(H_D2G02)
	 D_dg(2,2,4,3) = dyu_pt(H_D2G03); D_dg(2,2,4,4) = dyu_pt(H_D2G00)

         D_dg(2,3,1,1) = dyu_pt(H_D3G11); D_dg(2,3,1,2) = dyu_pt(H_D3G12)
	 D_dg(2,3,1,3) = dyu_pt(H_D3G13); D_dg(2,3,1,4) = dyu_pt(H_D3G01)
         D_dg(2,3,2,1) = dyu_pt(H_D3G12); D_dg(2,3,2,2) = dyu_pt(H_D3G22)
	 D_dg(2,3,2,3) = dyu_pt(H_D3G23); D_dg(2,3,2,4) = dyu_pt(H_D3G02)
         D_dg(2,3,3,1) = dyu_pt(H_D3G13); D_dg(2,3,3,2) = dyu_pt(H_D3G23)
	 D_dg(2,3,3,3) = dyu_pt(H_D3G33); D_dg(2,3,3,4) = dyu_pt(H_D3G03)		 
         D_dg(2,3,4,1) = dyu_pt(H_D3G01); D_dg(2,3,4,2) = dyu_pt(H_D3G02)
	 D_dg(2,3,4,3) = dyu_pt(H_D3G03); D_dg(2,3,4,4) = dyu_pt(H_D3G00)
	
         !---z direction -------------	
	 D_dg(3,1,1,1) = dzu_pt(H_D1G11); D_dg(3,1,1,2) = dzu_pt(H_D1G12)
	 D_dg(3,1,1,3) = dzu_pt(H_D1G13); D_dg(3,1,1,4) = dzu_pt(H_D1G01)
         D_dg(3,1,2,1) = dzu_pt(H_D1G12); D_dg(3,1,2,2) = dzu_pt(H_D1G22)
	 D_dg(3,1,2,3) = dzu_pt(H_D1G23); D_dg(3,1,2,4) = dzu_pt(H_D1G02)
         D_dg(3,1,3,1) = dzu_pt(H_D1G13); D_dg(3,1,3,2) = dzu_pt(H_D1G23)
	 D_dg(3,1,3,3) = dzu_pt(H_D1G33); D_dg(3,1,3,4) = dzu_pt(H_D1G03)		 
         D_dg(3,1,4,1) = dzu_pt(H_D1G01); D_dg(3,1,4,2) = dzu_pt(H_D1G02)
	 D_dg(3,1,4,3) = dzu_pt(H_D1G03); D_dg(3,1,4,4) = dzu_pt(H_D1G00)

         D_dg(3,2,1,1) = dzu_pt(H_D2G11); D_dg(3,2,1,2) = dzu_pt(H_D2G12)
	 D_dg(3,2,1,3) = dzu_pt(H_D2G13); D_dg(3,2,1,4) = dzu_pt(H_D2G01)
         D_dg(3,2,2,1) = dzu_pt(H_D2G12); D_dg(3,2,2,2) = dzu_pt(H_D2G22)
	 D_dg(3,2,2,3) = dzu_pt(H_D2G23); D_dg(3,2,2,4) = dzu_pt(H_D2G02)
         D_dg(3,2,3,1) = dzu_pt(H_D2G13); D_dg(3,2,3,2) = dzu_pt(H_D2G23)
	 D_dg(3,2,3,3) = dzu_pt(H_D2G33); D_dg(3,2,3,4) = dzu_pt(H_D2G03)		 
         D_dg(3,2,4,1) = dzu_pt(H_D2G01); D_dg(3,2,4,2) = dzu_pt(H_D2G02)
	 D_dg(3,2,4,3) = dzu_pt(H_D2G03); D_dg(3,2,4,4) = dzu_pt(H_D2G00)

         D_dg(3,3,1,1) = dzu_pt(H_D3G11); D_dg(3,3,1,2) = dzu_pt(H_D3G12)
	 D_dg(3,3,1,3) = dzu_pt(H_D3G13); D_dg(3,3,1,4) = dzu_pt(H_D3G01)
         D_dg(3,3,2,1) = dzu_pt(H_D3G12); D_dg(3,3,2,2) = dzu_pt(H_D3G22)
	 D_dg(3,3,2,3) = dzu_pt(H_D3G23); D_dg(3,3,2,4) = dzu_pt(H_D3G02)
         D_dg(3,3,3,1) = dzu_pt(H_D3G13); D_dg(3,3,3,2) = dzu_pt(H_D3G23)
	 D_dg(3,3,3,3) = dzu_pt(H_D3G33); D_dg(3,3,3,4) = dzu_pt(H_D3G03)		 
         D_dg(3,3,4,1) = dzu_pt(H_D3G01); D_dg(3,3,4,2) = dzu_pt(H_D3G02)
	 D_dg(3,3,4,3) = dzu_pt(H_D3G03); D_dg(3,3,4,4) = dzu_pt(H_D3G00)	 
	 					
         !---t direction -------------		
	 call compute_rhs_analysis(u_pt, dxu_pt, dyu_pt, dzu_pt, &
     &                             v_pt, dxv_pt, dyv_pt, dzv_pt, par, D_dg)
              
         !--compute the frame--------- 
         call frame(Tl, Tn, Tm_R, Tm_I, x, y, z, alp, beta, g, g_uu, deth)

         !--compute the 4D Rieman-----
         call rieman(R4, Ch, g_uu, g, dg, D_dg, 4)
	 
	 !-- compute the psi 4--------
         call calc_psi4(Psi4R, Psi4I, R4, Tn, Tm_R, Tm_I)
    
	 Psi4R = Psi4R * r
	 Psi4I = Psi4I * r

         !-- compute the 4 curvature --------------------------------
         call cal_curvature(curvature, R4, g_uu)
 

         !---the Bondi gur -----------------------------------
         !----------------------------------------------------
         !--------- r_a = (l_a - n_a)/sqrt(2)-----------------
	 
         gur = 0.0d0
         do la=1,4
           do lb=1,4
	     gur = gur + g_uu(la,lb)*Tn(la)*((Tl(lb)-Tn(lb))/dsqrt(2.0d0))	 
	  end do
	end do

         !-- compute the ADM mass ----------------------------------
	 ADMmass = 0.0         
	 do li=1,3
           do lm=1,3
             ADMmass = ADMmass + &
     &          (dg(li,li,lm)-dg(lm,li,li))*(Tl(lm)-Tn(lm))*r**2
	   end do
	 end do
         ADMmass = ADMmass/(16*pi*dsqrt(2.0d0)) 
	 
         !-- compute the linear and angular momentum-----------------
	 !--the M quantities------------
	 KtrK(1,1) = eK11 - g11*etrK
	 KtrK(1,2) = eK12 - g12*etrK
	 KtrK(1,3) = eK13 - g13*etrK
	 KtrK(2,2) = eK22 - g22*etrK
	 KtrK(2,3) = eK23 - g23*etrK
	 KtrK(3,3) = eK33 - g33*etrK	 

	 KtrK(2,1) = KtrK(1,2)	
	 KtrK(3,1) = KtrK(1,3)
	 KtrK(3,2) = KtrK(2,3)	
	
	 M = 0.0d0
         do li=1,3; do lm=1,3
	   M(li) = M(li) + KtrK(li,lm)*(Tl(lm)-Tn(lm))*r*r
	 end do; end do  	 
	 		         	 	 
	 lin_momx = M(1)/(8*pi*dsqrt(2.0d0)) 
	 lin_momy = M(2)/(8*pi*dsqrt(2.0d0)) 
	 lin_momz = M(3)/(8*pi*dsqrt(2.0d0))  	
	
	 ang_momx = (y*lin_momz - z*lin_momy)
	 ang_momy = (z*lin_momx - x*lin_momz)
	 ang_momz = (x*lin_momy - y*lin_momx)
	 
       end if
       
!--------------------------------------------------------------------------       
!----- compute the 2-curvature on the spheres normal to n------------------
!--------------------------------------------------------------------------       

       if ((psi4_analysis .EQ. 2) .OR. (psi4_analysis .EQ. 3)) then         
	                  	 
         !-- compute the 3 curvature -------------------------------  
         K(1,1) = eK11; K(1,2) = eK12; K(1,3) = eK13 
         K(2,1) = eK12; K(2,2) = eK22; K(2,3) = eK23 
         K(3,1) = eK13; K(3,2) = eK23; K(3,3) = eK33 

         h_uu(1,1) = huu11; h_uu(1,2) = huu12; h_uu(1,3) = huu13
         h_uu(2,1) = huu12; h_uu(2,2) = huu22; h_uu(2,3) = huu23
         h_uu(3,1) = huu13; h_uu(3,2) = huu23; h_uu(3,3) = huu33
  
         R3 = 0.0d0	 
         do li=1,3; do lj=1,3; do lk=1,3; do lm=1,3 
             R3(li,lj,lk,lm) = R4(li,lj,lk,lm) - K(li,lk)*K(lj,lm) + K(li,lm)*K(lk,lj)
         end do; end do; end do; end do
  
         Ricci3 = 0.0d0
         do li=1,3; do lj=1,3; do lk=1,3; do lm=1,3 
            Ricci3(li,lj) = Ricci3(li,lj) + R3(lk,li,lm,lj)*h_uu(lk,lm)
         end do; end do; end do; end do
	          
	 R3scalar = 0.0d0
         do li=1,3; do lj=1,3
            R3scalar = R3scalar + Ricci3(li,lj)*h_uu(li,lj)
         end do; end do
 
        !--other auxiliary variables -------------------	 
         nr1 = x; nr2 = y; nr3 = z             
         nrnorm = dsqrt(huu11*nr1*nr1 + huu22*nr2*nr2 + huu33*nr3*nr3 +& 
     &            2.0d0*(huu12*nr1*nr2 + huu13*nr1*nr3 + huu23*nr2*nr3))
         if (nrnorm .LT. 1E-8) nrnorm=1E-8 
         nr1 = nr1/nrnorm
         nr2 = nr2/nrnorm
         nr3 = nr3/nrnorm      
                                                 
         nru1 = huu11*nr1 + huu12*nr2 + huu13*nr3
         nru2 = huu12*nr1 + huu22*nr2 + huu23*nr3
         nru3 = huu13*nr1 + huu23*nr2 + huu33*nr3
                            
         Rnn = Ricci3(1,1)*nru1*nru1 + Ricci3(2,2)*nru2*nru2 + Ricci3(3,3)*nru3*nru3 &
     &    + 2.0d0*(Ricci3(1,2)*nru1*nru2 + Ricci3(1,3)*nru1*nru3 + Ricci3(2,3)*nru2*nru3)
 
         Gndd11 = nr1*G111 + nr2*G211 + nr3*G311
         Gndd12 = nr1*G112 + nr2*G212 + nr3*G312
         Gndd13 = nr1*G113 + nr2*G213 + nr3*G313
         Gndd22 = nr1*G122 + nr2*G222 + nr3*G322
         Gndd23 = nr1*G123 + nr2*G223 + nr3*G323
         Gndd33 = nr1*G133 + nr2*G233 + nr3*G333
              
         Dnn1 = D111*nru1*nru1 + D122*nru2*nru2 + D133*nru3*nru3 + &
     &    2.0d0*(D112*nru1*nru2 + D113*nru1*nru3 + D123*nru2*nru3) 
         Dnn2 = D211*nru1*nru1 + D222*nru2*nru2 + D233*nru3*nru3 + &
     &    2.0d0*(D212*nru1*nru2 + D213*nru1*nru3 + D223*nru2*nru3) 
         Dnn3 = D311*nru1*nru1 + D322*nru2*nru2 + D333*nru3*nru3 + &
     &    2.0d0*(D312*nru1*nru2 + D313*nru1*nru3 + D323*nru2*nru3) 
       
         D1n1 = (1.0d0 - nr1*nru1)/nrnorm + 0.50d0*Dnn1*nr1 - Gndd11
         D2n1 = (0.0d0 - nr1*nru2)/nrnorm + 0.50d0*Dnn2*nr1 - Gndd12
         D3n1 = (0.0d0 - nr1*nru3)/nrnorm + 0.50d0*Dnn3*nr1 - Gndd13                     
         D1n2 = (0.0d0 - nr2*nru1)/nrnorm + 0.50d0*Dnn1*nr2 - Gndd12
         D2n2 = (1.0d0 - nr2*nru2)/nrnorm + 0.50d0*Dnn2*nr2 - Gndd22
         D3n2 = (0.0d0 - nr2*nru3)/nrnorm + 0.50d0*Dnn3*nr2 - Gndd23              
         D1n3 = (0.0d0 - nr3*nru1)/nrnorm + 0.50d0*Dnn1*nr3 - Gndd13
         D2n3 = (0.0d0 - nr3*nru2)/nrnorm + 0.50d0*Dnn2*nr3 - Gndd23
         D3n3 = (1.0d0 - nr3*nru3)/nrnorm + 0.50d0*Dnn3*nr3 - Gndd33                     

         !-define the 2 projector over the hs----
         Sud11 = 1.0d0 - nru1*nr1
         Sud12 = - nru1*nr2
         Sud13 = - nru1*nr3	
         Sud21 = - nru2*nr1
         Sud22 = 1.0d0 - nru2*nr2
         Sud23 = - nru2*nr3        
         Sud31 = - nru3*nr1
         Sud32 = - nru3*nr2
         Sud33 = 1.0d0 - nru3*nr3

         Suu11 = huu11 - nru1*nru1
         Suu12 = huu12 - nru1*nru2
         Suu13 = huu13 - nru1*nru3	
         Suu22 = huu22 - nru2*nru2
         Suu23 = huu23 - nru2*nru3        
         Suu33 = huu33 - nru3*nru3
         
	 sK11 = D1n1*Sud11**2 + D2n2*Sud21**2 + (D3n1*Sud11 + (D2n3 +&
     &   D3n2)*Sud21)*Sud31 + D3n3*Sud31**2 + Sud11* ((D1n2 + D2n1)*Sud21&
     &   + D1n3*Sud31)
         sK12 = Sud12*(D1n1*Sud11 + D1n2*Sud21 + D1n3*Sud31) +&
     &   Sud22*(D2n1*Sud11 + D2n2*Sud21 + D2n3*Sud31) + (D3n1*Sud11 +&
     &   D3n2*Sud21 + D3n3*Sud31)*Sud32
         sK13 = Sud13*(D1n1*Sud11 + D1n2*Sud21 + D1n3*Sud31) +&
     &   Sud23*(D2n1*Sud11 + D2n2*Sud21 + D2n3*Sud31) + (D3n1*Sud11 +&
     &   D3n2*Sud21 + D3n3*Sud31)*Sud33
         sK22 = D1n1*Sud12**2 + D2n2*Sud22**2 + (D3n1*Sud12 + (D2n3 +&
     &   D3n2)*Sud22)*Sud32 + D3n3*Sud32**2 + Sud12* ((D1n2 + D2n1)*Sud22&
     &   + D1n3*Sud32)
         sK23 = Sud13*(D1n1*Sud12 + D1n2*Sud22 + D1n3*Sud32) +&
     &   Sud23*(D2n1*Sud12 + D2n2*Sud22 + D2n3*Sud32) + (D3n1*Sud12 +&
     &   D3n2*Sud22 + D3n3*Sud32)*Sud33
         sK33 = D1n1*Sud13**2 + D2n2*Sud23**2 + (D3n1*Sud13 + (D2n3 +&
     &   D3n2)*Sud23)*Sud33 + D3n3*Sud33**2 + Sud13* ((D1n2 + D2n1)*Sud23&
     &   + D1n3*Sud33)
           
         sKud11 = Suu11*sK11 + Suu12*sK12 + Suu13*sK13
         sKud12 = Suu11*sK12 + Suu12*sK22 + Suu13*sK23
         sKud13 = Suu11*sK13 + Suu12*sK23 + Suu13*sK33
         sKud21 = Suu12*sK11 + Suu22*sK12 + Suu23*sK13
         sKud22 = Suu12*sK12 + Suu22*sK22 + Suu23*sK23
         sKud23 = Suu12*sK13 + Suu22*sK23 + Suu23*sK33
         sKud31 = Suu13*sK11 + Suu23*sK12 + Suu33*sK13
         sKud32 = Suu13*sK12 + Suu23*sK22 + Suu33*sK23
         sKud33 = Suu13*sK13 + Suu23*sK23 + Suu33*sK33
         strK = sKud11 + sKud22 + sKud33       

         sKK = sKud11*sKud11 + sKud22*sKud22 + sKud33*sKud33 &
     &     + 2.0d0*(sKud12*sKud21 + sKud13*sKud31 + sKud23*sKud32)	 
       
         !---finally, the curvature in the 2-surface ----------
         R2scalar = R3scalar - 2*Rnn + strK*strK - sKK
	 R2scalar = R2scalar*r*r
	  
	 !---the intrinsic metric induced on the 2-sphere ----
	 !--- in spherical coordinates -----------------------
	          
         Sdd11 = (h11 - nr1*nr1)
         Sdd12 = (h12 - nr1*nr2)
         Sdd13 = (h13 - nr1*nr3)
         Sdd22 = (h22 - nr2*nr2)
         Sdd23 = (h23 - nr2*nr3)        
         Sdd33 = (h33 - nr3*nr3)

         rho = sqrt(x*x + y*y)	 
         if ( rho .LT. 1E-6 ) rho = 1E-6  
         costheta = z/r
	 sintheta = rho/r
	 cosphi   = x/rho
	 sinphi   = y/rho
           
         JA(1,1) = cosphi*sintheta
         JA(2,1) = sinphi*sintheta
         JA(3,1) = cosphi

         JA(1,2) = r*cosphi*costheta
         JA(2,2) = r*sinphi*costheta
         JA(3,2) = -r*sintheta
  
         JA(1,3) = -r*sinphi*sintheta
         JA(2,3) = r*cosphi*sintheta
         JA(3,3) = 0.0
     
!         Sdds11 = JA(1,1)*JA(1,1)*Sdd11 + JA(2,1)*JA(2,1)*Sdd22 + JA(3,1)*JA(3,1)*Sdd33 &
!     &     + 2.0*(JA(1,1)*JA(2,1)*Sdd12 + JA(1,1)*JA(3,1)*Sdd13 + JA(2,1)*JA(3,1)*Sdd23) 
         Sdds22 = JA(1,2)*JA(1,2)*Sdd11 + JA(2,2)*JA(2,2)*Sdd22 + JA(3,2)*JA(3,2)*Sdd33 &
     &     + 2.0*(JA(1,2)*JA(2,2)*Sdd12 + JA(1,2)*JA(3,2)*Sdd13 + JA(2,2)*JA(3,2)*Sdd23) 
         Sdds33 = JA(1,3)*JA(1,3)*Sdd11 + JA(2,3)*JA(2,3)*Sdd22 + JA(3,3)*JA(3,3)*Sdd33 &
     &     + 2.0*(JA(1,3)*JA(2,3)*Sdd12 + JA(1,3)*JA(3,3)*Sdd13 + JA(2,3)*JA(3,3)*Sdd23) 
         Sdds23 = JA(1,2)*JA(1,3)*Sdd11 + JA(2,2)*JA(2,3)*Sdd22 + JA(3,2)*JA(3,3)*Sdd33 &
     &     + 2.0*(JA(1,2)*JA(2,3)*Sdd12 + JA(1,2)*JA(3,3)*Sdd13 + JA(2,2)*JA(3,3)*Sdd23) 

         Sdds22 = Sdds22/(r*r)
         Sdds23 = Sdds23/(r*r)               
         Sdds33 = Sdds33/(r*r)
	 	
	
!	 if ((abs(z).LT. 0.001) .AND. (abs(x-2) .LT. 0.001) .AND. (abs(y-2) .LT. 0.001) ) then
!           print*, "======================================="
!	   print*, "x,y,z" , x, y, z
!           
!           print*, "hdd11=", h11
!           print*, "hdd12=", h12
!           print*, "hdd13=", h13
!           print*, "hdd22=", h22
!           print*, "hdd23=", h23
!           print*, "hdd33=", h33
!
!           print*, "n1=", nr1
!           print*, "n2=", nr2
!           print*, "n3=", nr3
!	   
!	   print*, "======================"	   
!	   
!           print*, "Sdd11=", Sdd11
!           print*, "Sdd12=", Sdd12
!           print*, "Sdd13=", Sdd13
!           print*, "Sdd22=", Sdd22
!           print*, "Sdd23=", Sdd23
!           print*, "Sdd33=", Sdd33
!	   
!	   print*, "======================"	   	   
!	   
!	   print*, "Sdds11=", Sdds11
!           print*, "Sdds22=", Sdds22
!           print*, "Sdds23=", Sdds23
!           print*, "Sdds33=", Sdds33
!	   
!           print*, "R3scalar=", R3scalar
!           print*, "Rnn     =", Rnn
!           print*, "strK^2  =", strK*strK
!           print*, "sKK     =", sKK
!           print*, "R2scalar=", R2scalar
!	   	   
!	   print*, "-------RIEMANN 4D-------------------"	   
!           do la=1,4; do lb=1,4; do lc=1,4; do ld=1,4
!              print*, la,lb,lc,ld,R4(la,lb,lc,ld)
!           end do; end do; end do; end do
!	   
!	   print*, "-------RICCI 4D-------------------"	   
!           do la=1,4; do lb=1,4
!	      temp = 0.0d0
!	      do lc=1,4; do ld=1,4
!                 temp = temp + R4(la,lc,lb,ld)*g_uu(lc,ld)
!              end do; end do
!              print*, la,lb,temp	   
!	   end do; end do  
!	   
!	   print*, "-------RIEMANN 3D-------------------"	   
!           do la=1,3; do lb=1,3; do lc=1,3; do ld=1,3
!              print*, la,lb,lc,ld,R3(la,lb,lc,ld)
!           end do; end do; end do; end do
!
!	   print*, "-------3D K-------------------"	   
!           do li=1,3; do lj=1,3
!              print*, li,lj,K(li,lj)
!          end do; end do
!
!	   print*, "-------3D huu-------------------"	   
!           do li=1,3; do lj=1,3
!              print*, li,lj,h_uu(li,lj)
!           end do; end do
!	   
!	   	   	   
!	   print*, "-------RICCI 3D-------------------"	   
!           do li=1,3; do lj=1,3
!              print*, li,lj,Ricci3(li,lj)
!           end do; end do
!	   	   
!	   stop
!         end if
	  
	  
       end if

!     !--FOR TESTING REASONS -----------------------------------------	        
!      !-------l=2,m=2 mode ------------------------         
!         costheta = z/r
!	  sintheta = dsqrt(1.0d0 - costheta**2)	 
!         cosphi   = x/(r*sintheta)
!         sinphi   = y/(r*sintheta)	 	          
!	 Y2p2 = sqrt(5.D0)*(costheta+1)**2/sqrt(0.3141592653589793D1)/8       
!         if (dabs(Y2p2) .LT. 1E5) then
!           Psi4R = Y2p2*(cosphi**2 - sinphi**2)
!           Psi4I = Y2p2*2.0d0*cosphi*sinphi
!	 else
!	   Psi4R = 0.0d0
!	   Psi4I = 0.0d0
!	 end if    
!        !--------l=2,m=0 mode ------------------------
!         costheta = z/r
!	 sintheta = dsqrt(1.0d0 - (z/r)**2)	 	          
!         Y20 = (sintheta**2)*dsqrt(15.0d0/(32.0d0*3.141592653589793))
!         Psi4R = Y20
!         Psi4I = 0.0d0
	 
