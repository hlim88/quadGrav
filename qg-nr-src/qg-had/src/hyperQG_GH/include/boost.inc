
      Jud = 0.0d0 

      if (boost_type .EQ. 0) then 
        Jud(0,0) = 1.0d0
        Jud(1,0) =-vx
        Jud(2,0) =-vy
        Jud(3,0) =-vz      
        Jud(1,1) = 1.0d0      
        Jud(2,2) = 1.0d0
        Jud(3,3) = 1.0d0            
      else
        cgamma = 1.0d0/dsqrt(1.0d0 - vx*vx - vy*vy)
        Jud(0,0) = cgamma
        Jud(0,1) = -vx*cgamma
        Jud(0,2) = -vy*cgamma
        Jud(1,0) = -vx*cgamma
        Jud(2,0) = -vy*cgamma
        Jud(1,1) = cgamma
        Jud(2,2) = cgamma
        Jud(3,3) = 1.0d0      
      end if
                      
      DJud = 0.0d0 

!-----reading the data ----------------------------

      !--setting the things to zero----------      
      gold     = 0.0d0
      Dgold    = 0.0d0     
      Dphirold = 0.0d0
      Dphicold = 0.0d0
      
      !--the metric components-------------      
      gold(0,0) = u(H_G00)%d(i,j,k)
      gold(0,1) = u(H_G01)%d(i,j,k)
      gold(0,2) = u(H_G02)%d(i,j,k)
      gold(0,3) = u(H_G03)%d(i,j,k) 
      gold(1,0) = u(H_G01)%d(i,j,k)
      gold(1,1) = u(H_G11)%d(i,j,k)
      gold(1,2) = u(H_G12)%d(i,j,k)
      gold(1,3) = u(H_G13)%d(i,j,k)
      gold(2,0) = u(H_G02)%d(i,j,k)
      gold(2,1) = u(H_G12)%d(i,j,k)
      gold(2,2) = u(H_G22)%d(i,j,k)
      gold(2,3) = u(H_G23)%d(i,j,k)
      gold(3,0) = u(H_G03)%d(i,j,k)
      gold(3,1) = u(H_G13)%d(i,j,k)
      gold(3,2) = u(H_G23)%d(i,j,k)
      gold(3,3) = u(H_G33)%d(i,j,k)
      
      !--the time derivatives------------------
      Dgold(0,:,:) = 0.0d0
            
      !--the spatial derivatives-------------      
      Dgold(1,0,0) = u(H_D1G00)%d(i,j,k)
      Dgold(1,0,1) = u(H_D1G01)%d(i,j,k)
      Dgold(1,0,2) = u(H_D1G02)%d(i,j,k)
      Dgold(1,0,3) = u(H_D1G03)%d(i,j,k) 
      Dgold(1,1,0) = u(H_D1G01)%d(i,j,k)
      Dgold(1,1,1) = u(H_D1G11)%d(i,j,k)
      Dgold(1,1,2) = u(H_D1G12)%d(i,j,k)
      Dgold(1,1,3) = u(H_D1G13)%d(i,j,k)
      Dgold(1,2,0) = u(H_D1G02)%d(i,j,k)
      Dgold(1,2,1) = u(H_D1G12)%d(i,j,k)
      Dgold(1,2,2) = u(H_D1G22)%d(i,j,k)
      Dgold(1,2,3) = u(H_D1G23)%d(i,j,k)
      Dgold(1,3,0) = u(H_D1G03)%d(i,j,k)
      Dgold(1,3,1) = u(H_D1G13)%d(i,j,k)
      Dgold(1,3,2) = u(H_D1G23)%d(i,j,k)
      Dgold(1,3,3) = u(H_D1G33)%d(i,j,k)

      Dgold(2,0,0) = u(H_D2G00)%d(i,j,k)
      Dgold(2,0,1) = u(H_D2G01)%d(i,j,k)
      Dgold(2,0,2) = u(H_D2G02)%d(i,j,k)
      Dgold(2,0,3) = u(H_D2G03)%d(i,j,k) 
      Dgold(2,1,0) = u(H_D2G01)%d(i,j,k)
      Dgold(2,1,1) = u(H_D2G11)%d(i,j,k)
      Dgold(2,1,2) = u(H_D2G12)%d(i,j,k)
      Dgold(2,1,3) = u(H_D2G13)%d(i,j,k)
      Dgold(2,2,0) = u(H_D2G02)%d(i,j,k)
      Dgold(2,2,1) = u(H_D2G12)%d(i,j,k)
      Dgold(2,2,2) = u(H_D2G22)%d(i,j,k)
      Dgold(2,2,3) = u(H_D2G23)%d(i,j,k)
      Dgold(2,3,0) = u(H_D2G03)%d(i,j,k)
      Dgold(2,3,1) = u(H_D2G13)%d(i,j,k)
      Dgold(2,3,2) = u(H_D2G23)%d(i,j,k)
      Dgold(2,3,3) = u(H_D2G33)%d(i,j,k)

      Dgold(3,0,0) = u(H_D3G00)%d(i,j,k)
      Dgold(3,0,1) = u(H_D3G01)%d(i,j,k)
      Dgold(3,0,2) = u(H_D3G02)%d(i,j,k)
      Dgold(3,0,3) = u(H_D3G03)%d(i,j,k) 
      Dgold(3,1,0) = u(H_D3G01)%d(i,j,k)
      Dgold(3,1,1) = u(H_D3G11)%d(i,j,k)
      Dgold(3,1,2) = u(H_D3G12)%d(i,j,k)
      Dgold(3,1,3) = u(H_D3G13)%d(i,j,k)
      Dgold(3,2,0) = u(H_D3G02)%d(i,j,k)
      Dgold(3,2,1) = u(H_D3G12)%d(i,j,k)
      Dgold(3,2,2) = u(H_D3G22)%d(i,j,k)
      Dgold(3,2,3) = u(H_D3G23)%d(i,j,k)
      Dgold(3,3,0) = u(H_D3G03)%d(i,j,k)
      Dgold(3,3,1) = u(H_D3G13)%d(i,j,k)
      Dgold(3,3,2) = u(H_D3G23)%d(i,j,k)
      Dgold(3,3,3) = u(H_D3G33)%d(i,j,k)
      
      !--the derivative of the SF--------------- 
      !-- in dphim is the fluid 4-velocity -----

      Dphirold(0) = u(H_PIR)%d(i,j,k)                  
      Dphicold(0) = u(H_PIC)%d(i,j,k)                  
      Dphimold(0) = u(H_PIM)%d(i,j,k)                  
            
      Dphirold(1) = u(H_D1PHIR)%d(i,j,k)
      Dphirold(2) = u(H_D2PHIR)%d(i,j,k)
      Dphirold(3) = u(H_D3PHIR)%d(i,j,k)      

      Dphicold(1) = u(H_D1PHIC)%d(i,j,k)
      Dphicold(2) = u(H_D2PHIC)%d(i,j,k)
      Dphicold(3) = u(H_D3PHIC)%d(i,j,k)      

      Dphimold(1) = u(H_D1PHIM)%d(i,j,k)
      Dphimold(2) = u(H_D2PHIM)%d(i,j,k)
      Dphimold(3) = u(H_D3PHIM)%d(i,j,k)      
            
!------------------------------------------------------------------      
!----ow get the stuff in the fixed system of reference-------------
!------------------------------------------------------------------

      !--setting the things to zero----------      
      gnew     = 0.0d0
      Dgnew    = 0.0d0     
      Dphirnew = 0.0d0
      Dphicnew = 0.0d0
      Dphimnew = 0.0d0

      !--converting the vectors --------------------------------      
      call boost_vec(Dphirnew, Dphirold, Jud)
      call boost_vec(Dphicnew, Dphicold, Jud)
      call boost_vec(Dphimnew, Dphimold, Jud)
                  
      !--converting the metric --------------------------------      
      call boost_metric(gnew, gold, Jud)

      !--converting the derivatives of the metric --------------      
      call boost_metric_derivs(Dgnew, Dgold, gold, Jud, DJud)       
                    
!------------------------------------------------------------------      
!------recovering the values --------------------------------------
!------------------------------------------------------------------      
      
      u(H_G00)%d(i,j,k) = gnew(0,0)
      u(H_G01)%d(i,j,k) = gnew(0,1) 
      u(H_G02)%d(i,j,k) = gnew(0,2)
      u(H_G03)%d(i,j,k) = gnew(0,3)                                         
      u(H_G11)%d(i,j,k) = gnew(1,1)
      u(H_G12)%d(i,j,k) = gnew(1,2)
      u(H_G13)%d(i,j,k) = gnew(1,3)
      u(H_G22)%d(i,j,k) = gnew(2,2)
      u(H_G23)%d(i,j,k) = gnew(2,3)
      u(H_G33)%d(i,j,k) = gnew(3,3)
             
      u(H_K00)%d(i,j,k) = Dgnew(0,0,0)
      u(H_K01)%d(i,j,k) = Dgnew(0,0,1) 
      u(H_K02)%d(i,j,k) = Dgnew(0,0,2)
      u(H_K03)%d(i,j,k) = Dgnew(0,0,3)                                         
      u(H_K11)%d(i,j,k) = Dgnew(0,1,1)
      u(H_K12)%d(i,j,k) = Dgnew(0,1,2)
      u(H_K13)%d(i,j,k) = Dgnew(0,1,3)
      u(H_K22)%d(i,j,k) = Dgnew(0,2,2)
      u(H_K23)%d(i,j,k) = Dgnew(0,2,3)
      u(H_K33)%d(i,j,k) = Dgnew(0,3,3)

      u(H_D1G00)%d(i,j,k) = Dgnew(1,0,0)
      u(H_D1G01)%d(i,j,k) = Dgnew(1,0,1) 
      u(H_D1G02)%d(i,j,k) = Dgnew(1,0,2)
      u(H_D1G03)%d(i,j,k) = Dgnew(1,0,3)                                         
      u(H_D1G11)%d(i,j,k) = Dgnew(1,1,1)
      u(H_D1G12)%d(i,j,k) = Dgnew(1,1,2)
      u(H_D1G13)%d(i,j,k) = Dgnew(1,1,3)
      u(H_D1G22)%d(i,j,k) = Dgnew(1,2,2)
      u(H_D1G23)%d(i,j,k) = Dgnew(1,2,3)
      u(H_D1G33)%d(i,j,k) = Dgnew(1,3,3)

      u(H_D2G00)%d(i,j,k) = Dgnew(2,0,0)
      u(H_D2G01)%d(i,j,k) = Dgnew(2,0,1) 
      u(H_D2G02)%d(i,j,k) = Dgnew(2,0,2)
      u(H_D2G03)%d(i,j,k) = Dgnew(2,0,3)                                         
      u(H_D2G11)%d(i,j,k) = Dgnew(2,1,1)
      u(H_D2G12)%d(i,j,k) = Dgnew(2,1,2)
      u(H_D2G13)%d(i,j,k) = Dgnew(2,1,3)
      u(H_D2G22)%d(i,j,k) = Dgnew(2,2,2)
      u(H_D2G23)%d(i,j,k) = Dgnew(2,2,3)
      u(H_D2G33)%d(i,j,k) = Dgnew(2,3,3)

      u(H_D3G00)%d(i,j,k) = Dgnew(3,0,0)
      u(H_D3G01)%d(i,j,k) = Dgnew(3,0,1) 
      u(H_D3G02)%d(i,j,k) = Dgnew(3,0,2)
      u(H_D3G03)%d(i,j,k) = Dgnew(3,0,3)                                         
      u(H_D3G11)%d(i,j,k) = Dgnew(3,1,1)
      u(H_D3G12)%d(i,j,k) = Dgnew(3,1,2)
      u(H_D3G13)%d(i,j,k) = Dgnew(3,1,3)
      u(H_D3G22)%d(i,j,k) = Dgnew(3,2,2)
      u(H_D3G23)%d(i,j,k) = Dgnew(3,2,3)
      u(H_D3G33)%d(i,j,k) = Dgnew(3,3,3)

                  
      u(H_PIR)%d(i,j,k)    = Dphirnew(0)
      u(H_D1PHIR)%d(i,j,k) = Dphirnew(1)      
      u(H_D2PHIR)%d(i,j,k) = Dphirnew(2)      
      u(H_D3PHIR)%d(i,j,k) = Dphirnew(3)            
      
      u(H_PIC)%d(i,j,k)    = Dphicnew(0)
      u(H_D1PHIC)%d(i,j,k) = Dphicnew(1)      
      u(H_D2PHIC)%d(i,j,k) = Dphicnew(2)      
      u(H_D3PHIC)%d(i,j,k) = Dphicnew(3)      

      u(H_PIM)%d(i,j,k)    = Dphimnew(0)
      u(H_D1PHIM)%d(i,j,k) = Dphimnew(1)      
      u(H_D2PHIM)%d(i,j,k) = Dphimnew(2)      
      u(H_D3PHIM)%d(i,j,k) = Dphimnew(3)      
