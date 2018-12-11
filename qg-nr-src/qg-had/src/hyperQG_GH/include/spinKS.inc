      r2 = x**2/2+y**2/2+z**2/2-a**2/2+sqrt(x**4+2*x**2*y**2+2*x**2*z**&
     &2-2*x**2*a**2+y**4+2*y**2*z**2-2*y**2*a**2+z**4+2*a**2*z**2+a**4)/2
      r = sqrt(r2)
      dtr = 0
      dxr = (x+1/sqrt(x**4+2*x**2*y**2+2*x**2*z**2-2*x**2*a**2+y**4+2*y&
     &**2*z**2-2*y**2*a**2+z**4+2*a**2*z**2+a**4)*(4*x**3+4*x*y**2+4*x*z&
     &**2-4*x*a**2)/4)/r/2
      dyr = (y+1/sqrt(x**4+2*x**2*y**2+2*x**2*z**2-2*x**2*a**2+y**4+2*y&
     &**2*z**2-2*y**2*a**2+z**4+2*a**2*z**2+a**4)*(4*x**2*y+4*y**3+4*y*z&
     &**2-4*y*a**2)/4)/r/2
      dzr = (z+1/sqrt(x**4+2*x**2*y**2+2*x**2*z**2-2*x**2*a**2+y**4+2*y&
     &**2*z**2-2*y**2*a**2+z**4+2*a**2*z**2+a**4)*(4*x**2*z+4*y**2*z+4*z&
     &**3+4*a**2*z)/4)/r/2

      H = M*r**3/(r**4+a**2*z**2)
      nult = 1
      nulx = (r*x+a*y)/(r**2+a**2)
      nuly = (r*y-a*x)/(r**2+a**2)
      nulz = z/r
      dtH = 3*M*r**2/(r**4+a**2*z**2)*dtr-4*M*r**6/(r**4+a**2*z**2)**2*dtr
      dtnult = 0
      dtnulx = dtr*x/(r**2+a**2)-2*(r*x+a*y)/(r**2+a**2)**2*r*dtr
      dtnuly = dtr*y/(r**2+a**2)-2*(r*y-a*x)/(r**2+a**2)**2*r*dtr
      dtnulz = -z/r**2*dtr
      dxH = 3*M*r**2/(r**4+a**2*z**2)*dxr-4*M*r**6/(r**4+a**2*z**2)**2*dxr
      dxnult = 0
      dxnulx = (dxr*x+r)/(r**2+a**2)-2*(r*x+a*y)/(r**2+a**2)**2*r*dxr
      dxnuly = (dxr*y-a)/(r**2+a**2)-2*(r*y-a*x)/(r**2+a**2)**2*r*dxr
      dxnulz = -z/r**2*dxr
      dyH = 3*M*r**2/(r**4+a**2*z**2)*dyr-4*M*r**6/(r**4+a**2*z**2)**2*dyr
      dynult = 0
      dynulx = (dyr*x+a)/(r**2+a**2)-2*(r*x+a*y)/(r**2+a**2)**2*r*dyr
      dynuly = (dyr*y+r)/(r**2+a**2)-2*(r*y-a*x)/(r**2+a**2)**2*r*dyr
      dynulz = -z/r**2*dyr
      dzH = 3*M*r**2/(r**4+a**2*z**2)*dzr-M*r**3/(r**4+a**2*z**2)**2*(4&
     &*r**3*dzr+2*a**2*z)
      dznult = 0
      dznulx = dzr*x/(r**2+a**2)-2*(r*x+a*y)/(r**2+a**2)**2*r*dzr
      dznuly = dzr*y/(r**2+a**2)-2*(r*y-a*x)/(r**2+a**2)**2*r*dzr
      dznulz = 1/r-z/r**2*dzr
       
