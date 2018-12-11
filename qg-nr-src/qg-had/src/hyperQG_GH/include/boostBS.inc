
      Jud(0,0) = 1.0d0
      Jud(0,1) = 0.0d0
      Jud(0,2) = 0.0d0
      Jud(0,3) = 0.0d0      
      Jud(1,0) = 0.0d0
      Jud(1,1) = 1.0d0
      Jud(1,2) = 0.0d0
      Jud(1,3) = 0.0d0
      Jud(2,0) = 0.0d0 
      Jud(2,1) = 0.0d0
      Jud(2,2) = 1.0d0
      Jud(2,3) = 0.0d0
      Jud(3,0) = 0.0d0
      Jud(3,1) = 0.0d0
      Jud(3,2) = 0.0d0
      Jud(3,3) = 1.0d0

      if (v2 .GT. 1E-2) then
       
        Jud(0,0) = cgamma      
        Jud(0,1) =-cgamma*vx
        Jud(0,2) =-cgamma*vy
        Jud(0,3) = 0.0d0      
        Jud(1,0) =-cgamma*vx      
        Jud(1,1) = 1.0d0 + (cgamma -1.0d0)*vx*vx/v2
        Jud(1,2) = (cgamma -1.0d0)*vx*vy/v2
        Jud(1,3) = 0.0d0
        Jud(2,0) =-cgamma*vy 
        Jud(2,1) = (cgamma -1.0d0)*vx*vy/v2
        Jud(2,2) = 1.0d0 + (cgamma -1.0d0)*vy*vy/v2
        Jud(2,3) = 0.0d0
        Jud(3,0) = 0.0d0
        Jud(3,1) = 0.0d0
        Jud(3,2) = 0.0d0
        Jud(3,3) = 1.0d0

        DJud = 0.d0
	
      end if

!-----reading the data ----------------------------

      !--the metric components-------------      
      gold(0,0) = g_e(0,0)      
      gold(0,1) = g_e(0,1)
      gold(0,2) = g_e(0,2)
      gold(0,3) = g_e(0,3)          
      gold(1,0) = g_e(0,1)      
      gold(1,1) = g_e(1,1)
      gold(1,2) = g_e(1,2)
      gold(1,3) = g_e(1,3)
      gold(2,0) = g_e(0,2)      
      gold(2,1) = g_e(1,2)
      gold(2,2) = g_e(2,2)
      gold(2,3) = g_e(2,3)
      gold(3,0) = g_e(0,3)      
      gold(3,1) = g_e(1,3)
      gold(3,2) = g_e(2,3)
      gold(3,3) = g_e(3,3)      
      
      !--the time derivatives------------------
      Dgold(0,0,0) = K_e(0,0)      
      Dgold(0,0,1) = K_e(0,1)
      Dgold(0,0,2) = K_e(0,2)
      Dgold(0,0,3) = K_e(0,3)
      Dgold(0,1,0) = K_e(0,1)      
      Dgold(0,1,1) = K_e(1,1)
      Dgold(0,1,2) = K_e(1,2)
      Dgold(0,1,3) = K_e(1,3)
      Dgold(0,2,0) = K_e(0,2)      
      Dgold(0,2,1) = K_e(1,2)
      Dgold(0,2,2) = K_e(2,2)
      Dgold(0,2,3) = K_e(2,3)
      Dgold(0,3,0) = K_e(0,3)      
      Dgold(0,3,1) = K_e(1,3)      
      Dgold(0,3,2) = K_e(2,3)
      Dgold(0,3,3) = K_e(3,3)      
            
      !--the spatial derivatives-------------      
      Dgold(1:3,0,0) = dg_e(1:3,0,0)      
      Dgold(1:3,0,1) = dg_e(1:3,0,1)
      Dgold(1:3,0,2) = dg_e(1:3,0,2)
      Dgold(1:3,0,3) = dg_e(1:3,0,3)      
      Dgold(1:3,1,0) = dg_e(1:3,0,1)      
      Dgold(1:3,1,1) = dg_e(1:3,1,1)
      Dgold(1:3,1,2) = dg_e(1:3,1,2)
      Dgold(1:3,1,3) = dg_e(1:3,1,3)
      Dgold(1:3,2,0) = dg_e(1:3,0,2)      
      Dgold(1:3,2,1) = dg_e(1:3,1,2)
      Dgold(1:3,2,2) = dg_e(1:3,2,2)
      Dgold(1:3,2,3) = dg_e(1:3,2,3)
      Dgold(1:3,3,0) = dg_e(1:3,0,3)            
      Dgold(1:3,3,1) = dg_e(1:3,1,3)
      Dgold(1:3,3,2) = dg_e(1:3,2,3)
      Dgold(1:3,3,3) = dg_e(1:3,3,3)


      !--the derivative of the SF--------------- 
      Dphirold(0:3) = dphir_e(0:3)            
      Dphicold(0:3) = dphic_e(0:3)            
      
!------------------------------------------------------------------      
!----ow get the stuff in the fixed system of reference-------------
!------------------------------------------------------------------

      !--setting the things to zero----------      
      gnew     = 0.0d0
      Dgnew    = 0.0d0     
      Dphirnew = 0.0d0
      Dphicnew = 0.0d0

      !--converting the vectors --------------------------------      
      do li=0,3; do lk=0,3
         Dphirnew(li) = Dphirnew(li) + Dphirold(lk) * Jud(lk,li)
         Dphicnew(li) = Dphicnew(li) + Dphicold(lk) * Jud(lk,li)	 
      end do; end do
      
            
      !--converting the metric --------------------------------      
      do li=0,3; do lj=0,3; 
	 do lk=0,3; do lm=0,3      
            gnew(li,lj) = gnew(li,lj) + gold(lk,lm) * Jud(lk,li) * Jud(lm,lj)	 
         end do; end do; 
      end do; end do

      !--converting the derivatives of the metric --------------      
      do li=0,3; do lj=0,3; do ll=0,3
	 
	 Temp1 = 0.0
	 Temp2 = 0.0  
	       
         do lk=0,3; do lm=0,3; do lp=0,3
	    Temp1 =  Temp1   &
	     &   + Dgold(lp,lk,lm) * Jud(lp,ll) * Jud(lk,li) * Jud(lm,lj) 
	 end do; end do ; end do
	   
	 do lk=0,3; do lm=0,3
	    Temp2 =  Temp2  &
	     &    + gold(lk,lm) * DJud(lk,li,ll) * Jud(lm,lj) &
	     &    + gold(lk,lm) * Jud(lk,li) * DJud(lm,lj,ll)    
	 end do; end do
	   
	 Dgnew(ll,li,lj) = Temp1 + Temp2	 
	 
	     
       end do; end do; end do
             
                    
!------------------------------------------------------------------      
!------recovering the values --------------------------------------
!------------------------------------------------------------------      
      
      g_e(0,0) = gnew(0,0)
      g_e(0,1) = gnew(0,1) 
      g_e(0,2) = gnew(0,2)
      g_e(0,3) = gnew(0,3)                                         
      g_e(1,1) = gnew(1,1)
      g_e(1,2) = gnew(1,2)
      g_e(1,3) = gnew(1,3)
      g_e(2,2) = gnew(2,2)
      g_e(2,3) = gnew(2,3)
      g_e(3,3) = gnew(3,3)
             
      dg_e(1:3,0,0) = Dgnew(1:3,0,0)
      dg_e(1:3,0,1) = Dgnew(1:3,0,1) 
      dg_e(1:3,0,2) = Dgnew(1:3,0,2)
      dg_e(1:3,0,3) = Dgnew(1:3,0,3)                                         
      dg_e(1:3,1,1) = Dgnew(1:3,1,1)
      dg_e(1:3,1,2) = Dgnew(1:3,1,2)
      dg_e(1:3,1,3) = Dgnew(1:3,1,3)
      dg_e(1:3,2,2) = Dgnew(1:3,2,2)
      dg_e(1:3,2,3) = Dgnew(1:3,2,3)
      dg_e(1:3,3,3) = Dgnew(1:3,3,3)

      K_e(0,0) = Dgnew(0,0,0)
      K_e(0,1) = Dgnew(0,0,1) 
      K_e(0,2) = Dgnew(0,0,2)
      K_e(0,3) = Dgnew(0,0,3)                                         
      K_e(1,1) = Dgnew(0,1,1)
      K_e(1,2) = Dgnew(0,1,2)
      K_e(1,3) = Dgnew(0,1,3)
      K_e(2,2) = Dgnew(0,2,2)
      K_e(2,3) = Dgnew(0,2,3)
      K_e(3,3) = Dgnew(0,3,3)

      dphir_e(0:3)   =  Dphirnew(0:3)
      dphic_e(0:3)   =  Dphicnew(0:3)      
       