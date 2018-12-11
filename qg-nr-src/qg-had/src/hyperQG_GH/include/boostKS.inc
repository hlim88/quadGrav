      cgamma = 1/sqrt(1-vx**2-vy**2)
      r2 = cgamma**2*(vx*(x-vx*t)+vy*(y-vy*t))**2+(x-vx*t)**2+(y-vy*t)**2+z**2
      r = sqrt(r2)
      dtr = (2*cgamma**2*(vx*(x-vx*t)+vy*(y-vy*t))*(-vx**2-vy**2)-2*vx*&
            (x-vx*t)-2*vy*(y-vy*t))/r/2
      dxr = (2*cgamma**2*(vx*(x-vx*t)+vy*(y-vy*t))*vx+2*x-2*vx*t)/r/2
      dyr = (2*cgamma**2*(vx*(x-vx*t)+vy*(y-vy*t))*vy+2*y-2*vy*t)/r/2
      dzr = z/r
      
      H = M/r
      nult = cgamma*(1-cgamma*(vx*(x-vx*t)+vy*(y-vy*t))/r)
      nulx = (x-vx*t)/r-vx*cgamma*(1-cgamma*(vx*(x-vx*t)+vy*(y-vy*t))/r)
      nuly = (y-vy*t)/r-vy*cgamma*(1-cgamma*(vx*(x-vx*t)+vy*(y-vy*t))/r)
      nulz = z/r
      dtH = -M/r**2*dtr
      dtnult = cgamma*(-cgamma*(-vx**2-vy**2)/r+cgamma*(vx*(x-vx*t)+vy*&
               (y-vy*t))/r**2*dtr)
      dtnulx = -vx/r-(x-vx*t)/r**2*dtr-vx*cgamma*(-cgamma*(-vx**2-vy**2&
                )/r+cgamma*(vx*(x-vx*t)+vy*(y-vy*t))/r**2*dtr)
      dtnuly = -vy/r-(y-vy*t)/r**2*dtr-vy*cgamma*(-cgamma*(-vx**2-vy**2&
                )/r+cgamma*(vx*(x-vx*t)+vy*(y-vy*t))/r**2*dtr)
      dtnulz = -z/r**2*dtr
      dxH = -M/r**2*dxr
      dxnult = cgamma*(-cgamma*vx/r+cgamma*(vx*(x-vx*t)+vy*(y-vy*t))/r**2*dxr)
      dxnulx = 1/r-(x-vx*t)/r**2*dxr-vx*cgamma*(-cgamma*vx/r+cgamma*(vx&
                *(x-vx*t)+vy*(y-vy*t))/r**2*dxr)
      dxnuly = -(y-vy*t)/r**2*dxr-vy*cgamma*(-cgamma*vx/r+cgamma*(vx*(x&
                -vx*t)+vy*(y-vy*t))/r**2*dxr)
      dxnulz = -z/r**2*dxr
      dyH = -M/r**2*dyr
      dynult = cgamma*(-cgamma*vy/r+cgamma*(vx*(x-vx*t)+vy*(y-vy*t))/r**2*dyr)
      dynulx = -(x-vx*t)/r**2*dyr-vx*cgamma*(-cgamma*vy/r+cgamma*(vx*(x&
                -vx*t)+vy*(y-vy*t))/r**2*dyr)
      dynuly = 1/r-(y-vy*t)/r**2*dyr-vy*cgamma*(-cgamma*vy/r+cgamma*(vx&
                *(x-vx*t)+vy*(y-vy*t))/r**2*dyr)
      dynulz = -z/r**2*dyr
      dzH = -M/r**2*dzr
      dznult = cgamma**2*(vx*(x-vx*t)+vy*(y-vy*t))/r**2*dzr
      dznulx = -(x-vx*t)/r**2*dzr-vx*cgamma**2*(vx*(x-vx*t)+vy*(y-vy*t))/r**2*dzr
      dznuly = -(y-vy*t)/r**2*dzr-vy*cgamma**2*(vx*(x-vx*t)+vy*(y-vy*t))/r**2*dzr
      dznulz = 1/r-z/r**2*dzr

      