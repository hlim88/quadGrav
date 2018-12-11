
      s1 = 1+2*a*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))-2*y*a*sin(c*x*y*z&
     &)*c*x*z*exp(-b*(x**2+y**2+z**2))-4*a*cos(c*x*y*z)*b*x**2*exp(-b*(x&
     &**2+y**2+z**2))+a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2-2&
     &*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)*c*x*&
     &z-4*a**2*cos(c*x*y*z)**2*b*x**2*exp(-b*(x**2+y**2+z**2))**2+y**2*a&
     &**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2+0.d0
      g11 = s1+4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*y*&
     &sin(c*x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b**2*x**4*exp(-b*(x**2+y**&
     &2+z**2))**2+y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*&
     &c**2*z**2+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**2*s&
     &in(c*x*y*z)*c*x*z+4*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x**2+y*&
     &*2+z**2))**2*y**2+z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*&
     &z)**2*c**2*y**2+4*y*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)*b+4*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x&
     &**2+y**2+z**2))**2*z**2+0.d0
      s1 = -a*sin(c*x*y*z)*c*x**2*z*exp(-b*(x**2+y**2+z**2))-4*x*a*cos(c&
     &*x*y*z)*b*y*exp(-b*(x**2+y**2+z**2))-a**2*sin(c*x*y*z)*c*x**2*z*ex&
     &p(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)-4*x*a**2*cos(c*x*y*z)**2*b*&
     &y*exp(-b*(x**2+y**2+z**2))**2+a**2*sin(c*x*y*z)**2*c**2*x**3*z**2*&
     &exp(-b*(x**2+y**2+z**2))**2*y+4*a**2*cos(c*x*y*z)*b*y**2*exp(-b*(x&
     &**2+y**2+z**2))**2*sin(c*x*y*z)*c*x**2*z+2*x**4*a**2*exp(-b*(x**2+&
     &y**2+z**2))**2*cos(c*x*y*z)*b*sin(c*x*y*z)*c*z+4*x**3*a**2*exp(-b*&
     &(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b**2*y+0.d0
      s2 = s1-y**2*a*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*z-y**2*a**2&
     &*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z*cos(c*x*y*z)+y**3*a*&
     &*2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*z**2*x+2*y**4*&
     &a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z*cos(c*x*y*z)*b+0.d0
      g12 = s2+4*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b&
     &**2*x+z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y&
     &*x+2*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*y**2*cos&
     &(c*x*y*z)*b+2*a**2*sin(c*x*y*z)*c*x**2*z**3*exp(-b*(x**2+y**2+z**2&
     &))**2*cos(c*x*y*z)*b+4*x*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+&
     &y**2+z**2))**2*z**2+0.d0
      s1 = -a*sin(c*x*y*z)*c*x**2*y*exp(-b*(x**2+y**2+z**2))-4*x*a*cos(c&
     &*x*y*z)*b*z*exp(-b*(x**2+y**2+z**2))-a**2*sin(c*x*y*z)*c*x**2*y*ex&
     &p(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)-4*x*a**2*cos(c*x*y*z)**2*b*&
     &z*exp(-b*(x**2+y**2+z**2))**2+x**3*y**2*a**2*sin(c*x*y*z)**2*c**2*&
     &exp(-b*(x**2+y**2+z**2))**2*z+4*x**2*y*a**2*cos(c*x*y*z)*b*z**2*ex&
     &p(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c+2*x**4*a**2*exp(-b*(x**2+&
     &y**2+z**2))**2*cos(c*x*y*z)*b*sin(c*x*y*z)*c*y+4*x**3*a**2*exp(-b*&
     &(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b**2*z+0.d0
      s2 = s1+y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2&
     &*z*x+2*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z**2*c&
     &os(c*x*y*z)*b+2*x**2*y**3*a**2*sin(c*x*y*z)*c*exp(-b*(x**2+y**2+z*&
     &*2))**2*cos(c*x*y*z)*b+4*x*y**2*a**2*cos(c*x*y*z)**2*b**2*z*exp(-b&
     &*(x**2+y**2+z**2))**2+0.d0
      g13 = s2-z**2*a*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*y-z**2*a**&
     &2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*y*cos(c*x*y*z)+z**3*a&
     &**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y**2*x+2*z**4&
     &*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*y*cos(c*x*y*z)*b+&
     &4*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b**2*x+0.d0
      s1 = 1+x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*&
     &z**2+4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*y*sin(&
     &c*x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*y**2+2*a*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))-2*y*a*sin(c&
     &*x*y*z)*c*x*z*exp(-b*(x**2+y**2+z**2))-4*a*cos(c*x*y*z)*b*y**2*exp&
     &(-b*(x**2+y**2+z**2))+a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2)&
     &)**2+0.d0
      g22 = s1-2*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x&
     &*y*z)*c*x*z-4*a**2*cos(c*x*y*z)**2*b*y**2*exp(-b*(x**2+y**2+z**2))&
     &**2+y**2*a**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**&
     &2))**2+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**2*sin(&
     &c*x*y*z)*c*x*z+4*a**2*cos(c*x*y*z)**2*b**2*y**4*exp(-b*(x**2+y**2+&
     &z**2))**2+z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*x**2+4*y*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z**2))**2&
     &*cos(c*x*y*z)*b+4*a**2*cos(c*x*y*z)**2*b**2*y**2*exp(-b*(x**2+y**2&
     &+z**2))**2*z**2+0.d0
      s1 = x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*z*&
     &y+2*x**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z**2*cos(&
     &c*x*y*z)*b+2*x**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b*&
     &y**2*sin(c*x*y*z)*c+4*x**2*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*&
     &x*y*z)**2*b**2*y*z-y**2*a*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*&
     &x-4*y*a*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)*b*z-y**2*a**2*exp(-b&
     &*(x**2+y**2+z**2))**2*cos(c*x*y*z)*sin(c*x*y*z)*c*x-4*y*a**2*exp(-&
     &b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b*z+0.d0
      s2 = s1+y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2&
     &*x**2*z+4*y**2*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*z&
     &**2*cos(c*x*y*z)*b+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x&
     &*y*z)*b*sin(c*x*y*z)*c*x+4*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*c&
     &os(c*x*y*z)**2*b**2*z+0.d0
      g23 = s2-z**2*a*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*x-z**2*a**&
     &2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*sin(c*x*y*z)*c*x+z**3*a&
     &**2*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)**2*c**2*x**2+2*z**4&
     &*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*cos(c*x*y*z)*b+&
     &4*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b**2*y+0.d0
      s1 = 1+x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*&
     &y**2+4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*y*sin(&
     &c*x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*z**2+y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2&
     &*c**2*x**2+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**2*&
     &sin(c*x*y*z)*c*x*z+4*a**2*cos(c*x*y*z)**2*b**2*y**2*exp(-b*(x**2+y&
     &**2+z**2))**2*z**2+2*a*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))+0.d0
      g33 = s1-2*y*a*sin(c*x*y*z)*c*x*z*exp(-b*(x**2+y**2+z**2))-4*a*cos&
     &(c*x*y*z)*b*z**2*exp(-b*(x**2+y**2+z**2))+a**2*cos(c*x*y*z)**2*exp&
     &(-b*(x**2+y**2+z**2))**2-2*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**&
     &2))**2*y*sin(c*x*y*z)*c*x*z-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y*&
     &*2+z**2))**2*b*z**2+y**2*a**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-&
     &b*(x**2+y**2+z**2))**2+4*y*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2&
     &+y**2+z**2))**2*cos(c*x*y*z)*b+4*a**2*cos(c*x*y*z)**2*b**2*z**4*ex&
     &p(-b*(x**2+y**2+z**2))**2+0.d0
      
      s2 = 4*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*c**&
     &2*z**2*x+4*a**2*cos(c*x*y*z)**2*b*x**3*exp(-b*(x**2+y**2+z**2))**2&
     &*y**2*c**2*z**2-12*a*cos(c*x*y*z)*b*x*exp(-b*(x**2+y**2+z**2))-2*y&
     &**2*a*cos(c*x*y*z)*c**2*z**2*x*exp(-b*(x**2+y**2+z**2))+8*y*a*sin(&
     &c*x*y*z)*c*x**2*z*b*exp(-b*(x**2+y**2+z**2))+8*a*cos(c*x*y*z)*b**2&
     &*x**3*exp(-b*(x**2+y**2+z**2))-8*y**2*a**2*sin(c*x*y*z)**2*c**2*x*&
     &*3*z**2*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      s1 = s2-4*a*sin(c*x*y*z)*c*y*z*exp(-b*(x**2+y**2+z**2))+4*y*a**2*s&
     &in(c*x*y*z)*c*z**3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b+4*y*&
     &*2*a**2*cos(c*x*y*z)**2*c**2*z**4*x*exp(-b*(x**2+y**2+z**2))**2*b-&
     &24*y*a**2*sin(c*x*y*z)*c*x**2*z**3*exp(-b*(x**2+y**2+z**2))**2*cos&
     &(c*x*y*z)*b**2+8*a**2*cos(c*x*y*z)**2*b**2*x*exp(-b*(x**2+y**2+z**&
     &2))**2*z**2-16*a**2*cos(c*x*y*z)**2*b**3*x**3*exp(-b*(x**2+y**2+z*&
     &*2))**2*z**2+8*a**2*cos(c*x*y*z)**2*b**2*x*exp(-b*(x**2+y**2+z**2)&
     &)**2*y**2-16*a**2*cos(c*x*y*z)**2*b**3*x**3*exp(-b*(x**2+y**2+z**2&
     &))**2*y**2+0.d0
      s2 = s1+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*&
     &y**3*cos(c*x*y*z)-8*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*&
     &y*z)**2*c**2*z**2*b*x+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(&
     &c*x*y*z)*c**3*z**3*cos(c*x*y*z)+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*&
     &(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z-24*a**2*cos(c*x*y*z)*b**2*y*&
     &*3*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x**2*z-16*a**2*cos(c&
     &*x*y*z)**2*b**3*x**5*exp(-b*(x**2+y**2+z**2))**2+32*a**2*cos(c*x*y&
     &*z)**2*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2+0.d0
      D111 = s2+2*y**3*a**2*sin(c*x*y*z)*c**3*x**2*z**3*exp(-b*(x**2+y*&
     &*2+z**2))**2*cos(c*x*y*z)+28*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z&
     &**2))**2*y*sin(c*x*y*z)*c*x**2*z*b-2*a**2*cos(c*x*y*z)**2*exp(-b*(&
     &x**2+y**2+z**2))**2*y**2*c**2*z**2*x-24*a**2*cos(c*x*y*z)*b**2*x**&
     &4*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)*c*z-4*a**2*cos(c*x*y*&
     &z)*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*y*z-12*a**2*cos(c*x*&
     &y*z)**2*exp(-b*(x**2+y**2+z**2))**2*b*x+4*a**2*sin(c*x*y*z)**2*c**&
     &2*y**2*z**2*exp(-b*(x**2+y**2+z**2))**2*x-8*z**4*a**2*exp(-b*(x**2&
     &+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y**2*b*x+0.d0
      s3 = 12*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))&
     &**2*c**2*x**2*z**2+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x&
     &*y*z)*c**3*z**3*cos(c*x*y*z)*x+28*a**2*cos(c*x*y*z)**2*b**2*x**2*e&
     &xp(-b*(x**2+y**2+z**2))**2*y+2*y**2*a**2*sin(c*x*y*z)*c**3*x**3*z*&
     &*3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)+0.d0
      s2 = s3-8*y**3*a**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y*&
     &*2+z**2))**2*b-a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*y*&
     &c**2*x**2*z**2-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*&
     &b*y+2*a*sin(c*x*y*z)*c*x**3*z*b*exp(-b*(x**2+y**2+z**2))+8*a*cos(c&
     &*x*y*z)*b**2*x**2*y*exp(-b*(x**2+y**2+z**2))+0.d0
      s3 = s2-2*a*sin(c*x*y*z)*c*x*z*exp(-b*(x**2+y**2+z**2))+6*y**2*a*s&
     &in(c*x*y*z)*c*x*z*b*exp(-b*(x**2+y**2+z**2))-16*a**2*cos(c*x*y*z)*&
     &*2*b**3*x**4*exp(-b*(x**2+y**2+z**2))**2*y-4*a*cos(c*x*y*z)*b*y*ex&
     &p(-b*(x**2+y**2+z**2))+0.d0
      s1 = s3-y*a*cos(c*x*y*z)*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))-2&
     &*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*z+&
     &4*a**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2*&
     &y+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y**2*sin(c*x*y*&
     &z)*c*x*z*b-6*a**2*sin(c*x*y*z)**2*c**2*x**4*z**2*b*exp(-b*(x**2+y*&
     &*2+z**2))**2*y+2*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z&
     &)**2*c**2*z**2+0.d0
      s3 = s1-2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*z**2*b+4*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*&
     &b**2+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c**2*&
     &y**3*b+4*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+y**2+z**2))**2*z&
     &**2+0.d0
      s2 = s3-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2+y**2+z**2))&
     &**2*z**2*y-y**3*a*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)*c**2*z**2-&
     &y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c**2*z**2+2*&
     &y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c**2*z**2*b+&
     &4*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y&
     &*z)*b+z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y+0.d0
      s3 = s2-2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*y**3*b-16*a**2*cos(c*x*y*z)*b**2*y**4*exp(-b*(x**2+y**2+z**2))*&
     &*2*sin(c*x*y*z)*c*x*z-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x*&
     &*2+y**2+z**2))**2*y**3+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c**3*y**2*cos(c*x*y*z)*x+2*y*a**2*cos(c*x*y*z)**2*c**2*x&
     &**2*z**4*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      D112 = s3-16*y**2*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z*&
     &*2))**2*cos(c*x*y*z)*b**2-6*y*a**2*sin(c*x*y*z)**2*c**2*x**2*z**4*&
     &exp(-b*(x**2+y**2+z**2))**2*b-8*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b&
     &*(x**2+y**2+z**2))**2*z**3*sin(c*x*y*z)*c+2*a**2*cos(c*x*y*z)**2*b&
     &*x**4*exp(-b*(x**2+y**2+z**2))**2*y*c**2*z**2-8*a**2*cos(c*x*y*z)*&
     &b**2*x**5*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z-24*a**2*cos&
     &(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2*y**2*sin(c*x*y*z)*&
     &c*z+0.d0
      s3 = -16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2+y**2+z**2))**&
     &2*z**3+y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*&
     &z-8*a**2*cos(c*x*y*z)*b**2*x**5*exp(-b*(x**2+y**2+z**2))**2*sin(c*&
     &x*y*z)*c*y+2*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))*&
     &*2*c**2*x**2*z-8*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**&
     &2))**2*y**3*sin(c*x*y*z)*c+0.d0
      s2 = s3+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c*x-16*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**2+z**&
     &2))**2*sin(c*x*y*z)*c*x*z**2-2*y**4*a**2*exp(-b*(x**2+y**2+z**2))*&
     &*2*sin(c*x*y*z)**2*c**2*z**3*b-24*a**2*cos(c*x*y*z)*b**2*x**3*exp(&
     &-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)*c*z**2+2*a**2*cos(c*x*y*z)*&
     &*2*b*x**4*exp(-b*(x**2+y**2+z**2))**2*y**2*c**2*z+0.d0
      s3 = s2-8*y**2*a**2*sin(c*x*y*z)**2*c**2*x**2*z**3*exp(-b*(x**2+y*&
     &*2+z**2))**2*b-6*a**2*sin(c*x*y*z)**2*c**2*x**4*y**2*b*exp(-b*(x**&
     &2+y**2+z**2))**2*z+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x&
     &*y*z)*c**3*z**2*cos(c*x*y*z)*x+4*y**2*a**2*cos(c*x*y*z)**2*c**2*x*&
     &*2*z**3*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      s1 = s3+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c*&
     &*2*z**3*b+4*y**2*a**2*cos(c*x*y*z)**2*b**2*z*exp(-b*(x**2+y**2+z**&
     &2))**2-z**3*a*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)*c**2*y**2-z**3&
     &*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c**2*y**2+8*a*co&
     &s(c*x*y*z)*b**2*x**2*z*exp(-b*(x**2+y**2+z**2))-2*a**2*cos(c*x*y*z&
     &)*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*y+0.d0
      s3 = s1-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*b*z+4*a&
     &**2*sin(c*x*y*z)**2*c**2*x**2*y**2*exp(-b*(x**2+y**2+z**2))**2*z-y&
     &**2*a*cos(c*x*y*z)*c**2*x**2*z*exp(-b*(x**2+y**2+z**2))+2*a*sin(c*&
     &x*y*z)*c*x**3*y*b*exp(-b*(x**2+y**2+z**2))+0.d0
      s2 = s3+6*y*a*sin(c*x*y*z)*c*x*z**2*b*exp(-b*(x**2+y**2+z**2))-16*&
     &a**2*cos(c*x*y*z)**2*b**3*x**4*exp(-b*(x**2+y**2+z**2))**2*z-2*a*s&
     &in(c*x*y*z)*c*x*y*exp(-b*(x**2+y**2+z**2))-6*a**2*sin(c*x*y*z)**2*&
     &c**2*x**2*y**4*b*exp(-b*(x**2+y**2+z**2))**2*z-4*a*cos(c*x*y*z)*b*&
     &z*exp(-b*(x**2+y**2+z**2))+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+&
     &z**2))**2*y*sin(c*x*y*z)*c*x*z**2*b+0.d0
      s3 = s2-a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*y**2*c**2&
     &*x**2*z+12*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*si&
     &n(c*x*y*z)*c*y+28*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x**2+y**2&
     &+z**2))**2*z+2*y**3*a**2*sin(c*x*y*z)*c**3*x**3*z**2*exp(-b*(x**2+&
     &y**2+z**2))**2*cos(c*x*y*z)+4*z**3*a**2*exp(-b*(x**2+y**2+z**2))**&
     &2*cos(c*x*y*z)**2*b**2+0.d0
      D113 = s3+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2&
     &*c**2*y**2*b-16*y*a**2*sin(c*x*y*z)*c*x*z**4*exp(-b*(x**2+y**2+z**&
     &2))**2*cos(c*x*y*z)*b**2+2*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*s&
     &in(c*x*y*z)**2*c**2*y**2-2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*s&
     &in(c*x*y*z)**2*c**2*y**2*b+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c**3*y**3*cos(c*x*y*z)*x-16*a**2*cos(c*x*y*z)**2*b**&
     &3*x**2*exp(-b*(x**2+y**2+z**2))**2*y**2*z+0.d0
      s2 = 4*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*c**&
     &2*z**2*x+4*a**2*cos(c*x*y*z)**2*b*x**3*exp(-b*(x**2+y**2+z**2))**2&
     &*y**2*c**2*z**2-4*a*cos(c*x*y*z)*b*x*exp(-b*(x**2+y**2+z**2))-2*y*&
     &*2*a*cos(c*x*y*z)*c**2*z**2*x*exp(-b*(x**2+y**2+z**2))+4*y*a*sin(c&
     &*x*y*z)*c*x**2*z*b*exp(-b*(x**2+y**2+z**2))-8*y**2*a**2*sin(c*x*y*&
     &z)**2*c**2*x**3*z**2*exp(-b*(x**2+y**2+z**2))**2*b-4*a*sin(c*x*y*z&
     &)*c*y*z*exp(-b*(x**2+y**2+z**2))+4*y*a**2*sin(c*x*y*z)*c*z**3*exp(&
     &-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b+0.d0
      s3 = s2+4*y**2*a**2*cos(c*x*y*z)**2*c**2*z**4*x*exp(-b*(x**2+y**2+&
     &z**2))**2*b-16*y*a**2*sin(c*x*y*z)*c*x**2*z**3*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)*b**2+24*a**2*cos(c*x*y*z)**2*b**2*x*exp(-b*(&
     &x**2+y**2+z**2))**2*y**2-16*a**2*cos(c*x*y*z)**2*b**3*x**3*exp(-b*&
     &(x**2+y**2+z**2))**2*y**2+0.d0
      s1 = s3-4*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*z**2*b*x+12*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**&
     &2*sin(c*x*y*z)*c*z-24*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**&
     &2+z**2))**2*sin(c*x*y*z)*c*x**2*z+2*y**3*a**2*sin(c*x*y*z)*c**3*x*&
     &*2*z**3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)+20*a**2*cos(c*x*y&
     &*z)*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)*c*x**2*z*b+0.d0
      s2 = s1-2*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*y**2*c*&
     &*2*z**2*x-16*a**2*cos(c*x*y*z)*b**2*x**4*exp(-b*(x**2+y**2+z**2))*&
     &*2*y*sin(c*x*y*z)*c*z-4*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))&
     &**2*sin(c*x*y*z)*c*y*z-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*b*x+4*a**2*sin(c*x*y*z)**2*c**2*y**2*z**2*exp(-b*(x**2+y**&
     &2+z**2))**2*x+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)&
     &*c**3*y*x**2*cos(c*x*y*z)-8*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*&
     &sin(c*x*y*z)*c*y**3*cos(c*x*y*z)*b**2-4*a**2*sin(c*x*y*z)**2*c**2*&
     &x**3*z**4*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      s3 = s2-16*x*a**2*cos(c*x*y*z)**2*b**3*y**2*exp(-b*(x**2+y**2+z**2&
     &))**2*z**2-16*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**&
     &2*b**3*x+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c&
     &**2*x+2*a**2*sin(c*x*y*z)*c**3*x**4*z**3*exp(-b*(x**2+y**2+z**2))*&
     &*2*y*cos(c*x*y*z)+0.d0
      D122 = s3-4*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2&
     &*c**2*z**2*b+4*y**3*a*b*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*z-&
     &8*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z*cos(c*x*y&
     &*z)*b**2+8*x*a*cos(c*x*y*z)*b**2*y**2*exp(-b*(x**2+y**2+z**2))+4*a&
     &**2*sin(c*x*y*z)**2*c**2*x**3*z**2*exp(-b*(x**2+y**2+z**2))**2-4*z&
     &**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y**2*b*x+0.d0
      s3 = -6*a**2*sin(c*x*y*z)**2*c**2*x**3*y**3*b*exp(-b*(x**2+y**2+z*&
     &*2))**2*z-4*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*&
     &c**2*z**3*x*b-y**2*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c&
     &*cos(c*x*y*z)-8*a**2*sin(c*x*y*z)*c*x**2*z**4*exp(-b*(x**2+y**2+z*&
     &*2))**2*cos(c*x*y*z)*b**2-4*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*&
     &sin(c*x*y*z)**2*c**2*y*b*z+0.d0
      s2 = s3-z**2*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*cos(c&
     &*x*y*z)+2*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))**2*&
     &c**2*x**3*z+2*a**2*sin(c*x*y*z)*c**3*x**4*z**2*exp(-b*(x**2+y**2+z&
     &**2))**2*y**2*cos(c*x*y*z)-6*a**2*sin(c*x*y*z)**2*c**2*x**3*z**3*e&
     &xp(-b*(x**2+y**2+z**2))**2*y*b+10*x**2*a**2*cos(c*x*y*z)*b*y**2*ex&
     &p(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c+4*a**2*sin(c*x*y*z)**2*c*&
     &*2*x**3*z*exp(-b*(x**2+y**2+z**2))**2*y+0.d0
      s3 = s2+24*x*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+y**2+z**2))*&
     &*2*z-24*a**2*cos(c*x*y*z)*b**2*y**2*exp(-b*(x**2+y**2+z**2))**2*si&
     &n(c*x*y*z)*c*x**2*z**2-8*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos&
     &(c*x*y*z)*b**2*sin(c*x*y*z)*c*z**2-8*x**4*a**2*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)*b**2*y**2*sin(c*x*y*z)*c-y**2*a*exp(-b*(x**2&
     &+y**2+z**2))*sin(c*x*y*z)*c+0.d0
      s1 = s3+4*y**2*a*b*z**2*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c-y*&
     &*3*a*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)*c**2*x*z+3*z**3*a**2*ex&
     &p(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y*x+10*a**2*sin(c*x&
     &*y*z)*c*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b-2*y**&
     &5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*z*x*b+2*z*&
     &*4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*cos(c*x*y*z)*b+0.d0
      s3 = s1-z**2*a*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c-z**3*a**2*e&
     &xp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c**2*x*y+2*z**5*a**2*ex&
     &p(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c**2*x*y*b-z**3*a*exp(-b&
     &*(x**2+y**2+z**2))*cos(c*x*y*z)*c**2*x*y-8*y**4*a**2*exp(-b*(x**2+&
     &y**2+z**2))**2*cos(c*x*y*z)*b**2*x**2*sin(c*x*y*z)*c+0.d0
      s2 = s3+2*a*sin(c*x*y*z)*c*x**2*z**2*b*exp(-b*(x**2+y**2+z**2))+2*&
     &x**2*a*sin(c*x*y*z)*c*y**2*b*exp(-b*(x**2+y**2+z**2))+8*x*a*cos(c*&
     &x*y*z)*b**2*y*z*exp(-b*(x**2+y**2+z**2))-8*z**4*a**2*exp(-b*(x**2+&
     &y**2+z**2))**2*sin(c*x*y*z)*c*y**2*cos(c*x*y*z)*b**2-2*z**5*a**2*e&
     &xp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y*x*b+2*z**4*a**2*&
     &exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*y**2*x**2*cos(c*x*y*&
     &z)+0.d0
      s3 = s2-8*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z**&
     &2*cos(c*x*y*z)*b**2+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*&
     &x*y*z)**2*c**2*x*z*b-16*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(&
     &c*x*y*z)**2*b**3*x*z+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)*c*cos(c*x*y*z)*b+4*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*c&
     &os(c*x*y*z)**2*c**2*x*y**3*b+0.d0
      D123 = s3+2*a**2*cos(c*x*y*z)**2*c**2*x**3*y*z**3*exp(-b*(x**2+y*&
     &*2+z**2))**2*b-16*x*a**2*cos(c*x*y*z)**2*b**3*y*exp(-b*(x**2+y**2+&
     &z**2))**2*z**3-y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*&
     &*2*c**2*x*z+3*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**&
     &2*c**2*z*x+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*&
     &*3*z**2*x**2*cos(c*x*y*z)+12*y**2*a**2*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c*z**2*cos(c*x*y*z)*b-16*x**3*a**2*exp(-b*(x**2+y**2&
     &+z**2))**2*cos(c*x*y*z)**2*b**3*y*z+0.d0
      s2 = 4*a**2*sin(c*x*y*z)**2*c**2*x**3*y**2*exp(-b*(x**2+y**2+z**2)&
     &)**2+4*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*c**&
     &2*z**2*x+4*a**2*cos(c*x*y*z)**2*b*x**3*exp(-b*(x**2+y**2+z**2))**2&
     &*y**2*c**2*z**2-4*a*cos(c*x*y*z)*b*x*exp(-b*(x**2+y**2+z**2))-2*y*&
     &*2*a*cos(c*x*y*z)*c**2*z**2*x*exp(-b*(x**2+y**2+z**2))+4*y*a*sin(c&
     &*x*y*z)*c*x**2*z*b*exp(-b*(x**2+y**2+z**2))-8*y**2*a**2*sin(c*x*y*&
     &z)**2*c**2*x**3*z**2*exp(-b*(x**2+y**2+z**2))**2*b-4*a*sin(c*x*y*z&
     &)*c*y*z*exp(-b*(x**2+y**2+z**2))+0.d0
      s3 = s2+12*y*a**2*sin(c*x*y*z)*c*z**3*exp(-b*(x**2+y**2+z**2))**2*&
     &cos(c*x*y*z)*b+4*y**2*a**2*cos(c*x*y*z)**2*c**2*z**4*x*exp(-b*(x**&
     &2+y**2+z**2))**2*b-24*y*a**2*sin(c*x*y*z)*c*x**2*z**3*exp(-b*(x**2&
     &+y**2+z**2))**2*cos(c*x*y*z)*b**2+24*a**2*cos(c*x*y*z)**2*b**2*x*e&
     &xp(-b*(x**2+y**2+z**2))**2*z**2+0.d0
      s1 = s3-16*a**2*cos(c*x*y*z)**2*b**3*x**3*exp(-b*(x**2+y**2+z**2))&
     &**2*z**2-4*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c&
     &**2*z**2*b*x+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**&
     &2*sin(c*x*y*z)*c*z-16*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**&
     &2+z**2))**2*sin(c*x*y*z)*c*x**2*z+2*y**3*a**2*sin(c*x*y*z)*c**3*x*&
     &*2*z**3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)+0.d0
      s2 = s1+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x&
     &*y*z)*c*x**2*z*b-2*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**&
     &2*y**2*c**2*z**2*x-16*a**2*cos(c*x*y*z)*b**2*x**4*exp(-b*(x**2+y**&
     &2+z**2))**2*y*sin(c*x*y*z)*c*z-4*a**2*cos(c*x*y*z)*exp(-b*(x**2+y*&
     &*2+z**2))**2*sin(c*x*y*z)*c*y*z-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**&
     &2+y**2+z**2))**2*b*x+4*a**2*sin(c*x*y*z)**2*c**2*y**2*z**2*exp(-b*&
     &(x**2+y**2+z**2))**2*x-8*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c*y**3*cos(c*x*y*z)*b**2-16*x*a**2*cos(c*x*y*z)**2*b**3*&
     &y**2*exp(-b*(x**2+y**2+z**2))**2*z**2+0.d0
      s3 = s2+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*x+2*x**4*y**3*a**2*sin(c*x*y*z)*c**3*exp(-b*(x**2+y**2+z**2))**&
     &2*z*cos(c*x*y*z)+8*x*a*cos(c*x*y*z)*b**2*z**2*exp(-b*(x**2+y**2+z*&
     &*2))-4*x**3*y**4*a**2*sin(c*x*y*z)**2*c**2*exp(-b*(x**2+y**2+z**2)&
     &)**2*b+0.d0
      D133 = s3+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*&
     &*3*z*x**2*cos(c*x*y*z)-16*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*co&
     &s(c*x*y*z)**2*b**3*x-8*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)*c*y*cos(c*x*y*z)*b**2-4*x**5*a**2*exp(-b*(x**2+y**2+z**2))&
     &**2*sin(c*x*y*z)**2*c**2*y**2*b+4*z**3*a*b*exp(-b*(x**2+y**2+z**2)&
     &)*sin(c*x*y*z)*c*y-4*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x&
     &*y*z)**2*c**2*y**2*b*x+0.d0
      s2 = 12*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))&
     &**2*c**2*x**2*z**2+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x&
     &*y*z)*c**3*z**3*cos(c*x*y*z)*x+24*a**2*cos(c*x*y*z)**2*b**2*x**2*e&
     &xp(-b*(x**2+y**2+z**2))**2*y+2*y**2*a**2*sin(c*x*y*z)*c**3*x**3*z*&
     &*3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)-8*y**3*a**2*sin(c*x*y*&
     &z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2*b-2*a**2*cos(c*x*&
     &y*z)**2*exp(-b*(x**2+y**2+z**2))**2*y*c**2*x**2*z**2-4*a**2*cos(c*&
     &x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*b*y+0.d0
      s1 = s2+4*a*sin(c*x*y*z)*c*x**3*z*b*exp(-b*(x**2+y**2+z**2))+8*a*c&
     &os(c*x*y*z)*b**2*x**2*y*exp(-b*(x**2+y**2+z**2))-4*a*sin(c*x*y*z)*&
     &c*x*z*exp(-b*(x**2+y**2+z**2))+4*y**2*a*sin(c*x*y*z)*c*x*z*b*exp(-&
     &b*(x**2+y**2+z**2))-16*a**2*cos(c*x*y*z)**2*b**3*x**4*exp(-b*(x**2&
     &+y**2+z**2))**2*y-4*a*cos(c*x*y*z)*b*y*exp(-b*(x**2+y**2+z**2))-2*&
     &y*a*cos(c*x*y*z)*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))-4*a**2*co&
     &s(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*z+4*a**2*s&
     &in(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2*y+0.d0
      s2 = s1+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y**2*sin(&
     &c*x*y*z)*c*x*z*b-4*a**2*sin(c*x*y*z)**2*c**2*x**4*z**2*b*exp(-b*(x&
     &**2+y**2+z**2))**2*y+4*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)**2*c**2*z**2-4*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)**2*c**2*z**2*b-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x&
     &**2+y**2+z**2))**2*z**2*y+4*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**&
     &2+y**2+z**2))**2*cos(c*x*y*z)*b+2*z**4*a**2*exp(-b*(x**2+y**2+z**2&
     &))**2*sin(c*x*y*z)**2*c**2*y-4*z**4*a**2*exp(-b*(x**2+y**2+z**2))*&
     &*2*sin(c*x*y*z)**2*c**2*y**3*b+0.d0
      s3 = s2-16*a**2*cos(c*x*y*z)*b**2*y**4*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c*x*z-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2&
     &+y**2+z**2))**2*y**3+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)*c**3*y**2*cos(c*x*y*z)*x+4*y*a**2*cos(c*x*y*z)**2*c**2*x**&
     &2*z**4*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      D211 = s3-16*y**2*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z*&
     &*2))**2*cos(c*x*y*z)*b**2-4*y*a**2*sin(c*x*y*z)**2*c**2*x**2*z**4*&
     &exp(-b*(x**2+y**2+z**2))**2*b-8*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b&
     &*(x**2+y**2+z**2))**2*z**3*sin(c*x*y*z)*c+4*a**2*cos(c*x*y*z)**2*b&
     &*x**4*exp(-b*(x**2+y**2+z**2))**2*y*c**2*z**2-8*a**2*cos(c*x*y*z)*&
     &b**2*x**5*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z-24*a**2*cos&
     &(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2*y**2*sin(c*x*y*z)*&
     &c*z+0.d0
      s3 = 2*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*c**&
     &2*z**2*x+4*a**2*cos(c*x*y*z)**2*b*x**3*exp(-b*(x**2+y**2+z**2))**2&
     &*y**2*c**2*z**2-a*cos(c*x*y*z)*c**2*x**3*z**2*exp(-b*(x**2+y**2+z*&
     &*2))-4*a*cos(c*x*y*z)*b*x*exp(-b*(x**2+y**2+z**2))-y**2*a*cos(c*x*&
     &y*z)*c**2*z**2*x*exp(-b*(x**2+y**2+z**2))+0.d0
      s2 = s3+6*y*a*sin(c*x*y*z)*c*x**2*z*b*exp(-b*(x**2+y**2+z**2))-8*y&
     &**2*a**2*sin(c*x*y*z)**2*c**2*x**3*z**2*exp(-b*(x**2+y**2+z**2))**&
     &2*b-2*a*sin(c*x*y*z)*c*y*z*exp(-b*(x**2+y**2+z**2))+4*y*a**2*sin(c&
     &*x*y*z)*c*z**3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b+2*y**2*a&
     &**2*cos(c*x*y*z)**2*c**2*z**4*x*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      s3 = s2-16*y*a**2*sin(c*x*y*z)*c*x**2*z**3*exp(-b*(x**2+y**2+z**2)&
     &)**2*cos(c*x*y*z)*b**2+4*a**2*cos(c*x*y*z)**2*b**2*x*exp(-b*(x**2+&
     &y**2+z**2))**2*z**2+28*a**2*cos(c*x*y*z)**2*b**2*x*exp(-b*(x**2+y*&
     &*2+z**2))**2*y**2-16*a**2*cos(c*x*y*z)**2*b**3*x**3*exp(-b*(x**2+y&
     &**2+z**2))**2*y**2+0.d0
      s1 = s3-6*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*z**2*b*x+12*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**&
     &2*sin(c*x*y*z)*c*z-24*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**&
     &2+z**2))**2*sin(c*x*y*z)*c*x**2*z+4*a**2*cos(c*x*y*z)**2*b**2*x**3&
     &*exp(-b*(x**2+y**2+z**2))**2+2*y**3*a**2*sin(c*x*y*z)*c**3*x**2*z*&
     &*3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)+20*a**2*cos(c*x*y*z)*e&
     &xp(-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)*c*x**2*z*b+0.d0
      s3 = s1-a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*y**2*c**2&
     &*z**2*x-16*a**2*cos(c*x*y*z)*b**2*x**4*exp(-b*(x**2+y**2+z**2))**2&
     &*y*sin(c*x*y*z)*c*z-2*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**&
     &2*sin(c*x*y*z)*c*y*z-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2&
     &))**2*b*x+0.d0
      s2 = s3+4*a**2*sin(c*x*y*z)**2*c**2*y**2*z**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*x+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**&
     &3*y*x**2*cos(c*x*y*z)-8*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(&
     &c*x*y*z)*c*y**3*cos(c*x*y*z)*b**2-2*a**2*sin(c*x*y*z)**2*c**2*x**3&
     &*z**4*exp(-b*(x**2+y**2+z**2))**2*b-16*x*a**2*cos(c*x*y*z)**2*b**3&
     &*y**2*exp(-b*(x**2+y**2+z**2))**2*z**2-16*y**4*a**2*exp(-b*(x**2+y&
     &**2+z**2))**2*cos(c*x*y*z)**2*b**3*x+0.d0
      s3 = s2+z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2&
     &*x+2*a**2*sin(c*x*y*z)*c**3*x**4*z**3*exp(-b*(x**2+y**2+z**2))**2*&
     &y*cos(c*x*y*z)-2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z&
     &)**2*c**2*z**2*b+2*y**3*a*b*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*&
     &c*z-8*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*z*cos(c&
     &*x*y*z)*b**2+0.d0
      D212 = s3+8*x*a*cos(c*x*y*z)*b**2*y**2*exp(-b*(x**2+y**2+z**2))+2&
     &*a**2*sin(c*x*y*z)**2*c**2*x**3*z**2*exp(-b*(x**2+y**2+z**2))**2-a&
     &**2*cos(c*x*y*z)**2*c**2*x**3*z**2*exp(-b*(x**2+y**2+z**2))**2+2*x&
     &**5*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b*c**2*z**2+2&
     &*a**2*cos(c*x*y*z)**2*c**2*x**3*z**4*exp(-b*(x**2+y**2+z**2))**2*b&
     &-6*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y**2&
     &*b*x+0.d0
      s3 = -6*a**2*sin(c*x*y*z)**2*c**2*x**3*y**3*b*exp(-b*(x**2+y**2+z*&
     &*2))**2*z-6*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*&
     &c**2*z**3*x*b-a**2*sin(c*x*y*z)*c*x**2*exp(-b*(x**2+y**2+z**2))**2&
     &*cos(c*x*y*z)-8*a**2*sin(c*x*y*z)*c*x**2*z**4*exp(-b*(x**2+y**2+z*&
     &*2))**2*cos(c*x*y*z)*b**2-2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*&
     &sin(c*x*y*z)**2*c**2*y*b*z+0.d0
      s2 = s3-z**2*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*cos(c&
     &*x*y*z)+2*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))**2*&
     &c**2*x**3*z+2*a**2*sin(c*x*y*z)*c**3*x**4*z**2*exp(-b*(x**2+y**2+z&
     &**2))**2*y**2*cos(c*x*y*z)-4*a**2*sin(c*x*y*z)**2*c**2*x**3*z**3*e&
     &xp(-b*(x**2+y**2+z**2))**2*y*b+10*x**2*a**2*cos(c*x*y*z)*b*y**2*ex&
     &p(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c+3*a**2*sin(c*x*y*z)**2*c*&
     &*2*x**3*z*exp(-b*(x**2+y**2+z**2))**2*y+0.d0
      s3 = s2+24*x*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+y**2+z**2))*&
     &*2*z-24*a**2*cos(c*x*y*z)*b**2*y**2*exp(-b*(x**2+y**2+z**2))**2*si&
     &n(c*x*y*z)*c*x**2*z**2-8*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos&
     &(c*x*y*z)*b**2*sin(c*x*y*z)*c*z**2-8*x**4*a**2*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)*b**2*y**2*sin(c*x*y*z)*c+2*y**2*a*b*z**2*exp&
     &(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c+0.d0
      s1 = s3+3*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*y*x+12*a**2*sin(c*x*y*z)*c*x**2*z**2*exp(-b*(x**2+y**2+z**2))**&
     &2*cos(c*x*y*z)*b-4*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y&
     &*z)**2*c**2*z*x*b+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*&
     &y*z)*c*cos(c*x*y*z)*b-z**2*a*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)&
     &*c-z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c**2*x*y+0.d0
      s3 = s1+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c*&
     &*2*x*y*b-z**3*a*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)*c**2*x*y-8*y&
     &**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b**2*x**2*sin(c*&
     &x*y*z)*c+4*a*sin(c*x*y*z)*c*x**2*z**2*b*exp(-b*(x**2+y**2+z**2))+2&
     &*x**2*a*sin(c*x*y*z)*c*y**2*b*exp(-b*(x**2+y**2+z**2))+0.d0
      s2 = s3+8*x*a*cos(c*x*y*z)*b**2*y*z*exp(-b*(x**2+y**2+z**2))-8*z**&
     &4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*y**2*cos(c*x*y*z&
     &)*b**2-2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**&
     &2*y*x*b+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*&
     &y**2*x**2*cos(c*x*y*z)-8*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c*z**2*cos(c*x*y*z)*b**2-16*y**3*a**2*exp(-b*(x**2+y**2+&
     &z**2))**2*cos(c*x*y*z)**2*b**3*x*z+0.d0
      s3 = s2+2*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*c*&
     &*2*x*y**3*b+4*a**2*cos(c*x*y*z)**2*c**2*x**3*y*z**3*exp(-b*(x**2+y&
     &**2+z**2))**2*b-16*x*a**2*cos(c*x*y*z)**2*b**3*y*exp(-b*(x**2+y**2&
     &+z**2))**2*z**3+4*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*&
     &z)**2*c**2*z*x+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z&
     &)*c**3*z**2*x**2*cos(c*x*y*z)+0.d0
      D213 = s3+10*y**2*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c&
     &*z**2*cos(c*x*y*z)*b-16*x**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(&
     &c*x*y*z)**2*b**3*y*z+2*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c&
     &*x*y*z)*b*sin(c*x*y*z)*c+2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*c&
     &os(c*x*y*z)**2*b*c**2*y*z-a**2*cos(c*x*y*z)**2*c**2*x**3*y*z*exp(-&
     &b*(x**2+y**2+z**2))**2-a*cos(c*x*y*z)*c**2*x**3*y*z*exp(-b*(x**2+y&
     &**2+z**2))-a*sin(c*x*y*z)*c*x**2*exp(-b*(x**2+y**2+z**2))+0.d0
      s2 = 4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*sin(c*&
     &x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))*&
     &*2*c**2*x**2*z**2-16*a**2*cos(c*x*y*z)**2*b**3*y**3*exp(-b*(x**2+y&
     &**2+z**2))**2*z**2+8*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x**2+y&
     &**2+z**2))**2*y+2*y**2*a**2*sin(c*x*y*z)*c**3*x**3*z**3*exp(-b*(x*&
     &*2+y**2+z**2))**2*cos(c*x*y*z)-8*y**3*a**2*sin(c*x*y*z)**2*c**2*x*&
     &*2*z**2*exp(-b*(x**2+y**2+z**2))**2*b-2*a**2*cos(c*x*y*z)**2*exp(-&
     &b*(x**2+y**2+z**2))**2*y*c**2*x**2*z**2+0.d0
      s1 = s2-12*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*b*y-4*&
     &a*sin(c*x*y*z)*c*x*z*exp(-b*(x**2+y**2+z**2))+8*y**2*a*sin(c*x*y*z&
     &)*c*x*z*b*exp(-b*(x**2+y**2+z**2))-12*a*cos(c*x*y*z)*b*y*exp(-b*(x&
     &**2+y**2+z**2))-2*y*a*cos(c*x*y*z)*c**2*x**2*z**2*exp(-b*(x**2+y**&
     &2+z**2))-4*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y&
     &*z)*c*x*z+4*a**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+&
     &z**2))**2*y+28*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y**2*&
     &sin(c*x*y*z)*c*x*z*b+0.d0
      s2 = s1-8*a**2*sin(c*x*y*z)**2*c**2*x**4*z**2*b*exp(-b*(x**2+y**2+&
     &z**2))**2*y+32*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*&
     &*2*b**2+8*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+y**2+z**2))**2*&
     &z**2+4*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z**2))**2*cos(&
     &c*x*y*z)*b-24*a**2*cos(c*x*y*z)*b**2*y**4*exp(-b*(x**2+y**2+z**2))&
     &**2*sin(c*x*y*z)*c*x*z-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x&
     &**2+y**2+z**2))**2*y**3+4*y*a**2*cos(c*x*y*z)**2*c**2*x**2*z**4*ex&
     &p(-b*(x**2+y**2+z**2))**2*b+0.d0
      D222 = s2-24*y**2*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z*&
     &*2))**2*cos(c*x*y*z)*b**2-8*y*a**2*sin(c*x*y*z)**2*c**2*x**2*z**4*&
     &exp(-b*(x**2+y**2+z**2))**2*b+4*a**2*cos(c*x*y*z)**2*b*x**4*exp(-b&
     &*(x**2+y**2+z**2))**2*y*c**2*z**2-24*a**2*cos(c*x*y*z)*b**2*x**3*e&
     &xp(-b*(x**2+y**2+z**2))**2*y**2*sin(c*x*y*z)*c*z+2*z**5*a**2*exp(-&
     &b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*x**3*cos(c*x*y*z)+2*x**5*&
     &a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*z**3*cos(c*x*y*&
     &z)+8*a*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**2+z**2))-16*a**2*cos&
     &(c*x*y*z)**2*b**3*y**5*exp(-b*(x**2+y**2+z**2))**2+0.d0
      s3 = -z**3*a*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)*c**2*x**2-16*a*&
     &*2*cos(c*x*y*z)**2*b**3*y**2*exp(-b*(x**2+y**2+z**2))**2*z**3-16*a&
     &**2*cos(c*x*y*z)**2*b**3*y**4*exp(-b*(x**2+y**2+z**2))**2*z-8*a**2&
     &*cos(c*x*y*z)*b**2*y**5*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c&
     &*x+x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*z+0.d0
      s2 = s3+8*a*cos(c*x*y*z)*b**2*y**2*z*exp(-b*(x**2+y**2+z**2))+2*a*&
     &*2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*c**2*x**2*z-&
     &8*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2*y**3*sin&
     &(c*x*y*z)*c+12*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**&
     &2*sin(c*x*y*z)*c*x-24*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**&
     &2+z**2))**2*sin(c*x*y*z)*c*x*z**2+0.d0
      s3 = s2-16*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2&
     &*y*sin(c*x*y*z)*c*z**2+2*a**2*cos(c*x*y*z)**2*b*x**4*exp(-b*(x**2+&
     &y**2+z**2))**2*y**2*c**2*z-8*y**2*a**2*sin(c*x*y*z)**2*c**2*x**2*z&
     &**3*exp(-b*(x**2+y**2+z**2))**2*b-6*a**2*sin(c*x*y*z)**2*c**2*x**4&
     &*y**2*b*exp(-b*(x**2+y**2+z**2))**2*z+0.d0
      s1 = s3+4*y**2*a**2*cos(c*x*y*z)**2*c**2*x**2*z**3*exp(-b*(x**2+y*&
     &*2+z**2))**2*b+28*y**2*a**2*cos(c*x*y*z)**2*b**2*z*exp(-b*(x**2+y*&
     &*2+z**2))**2-2*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*sin(c&
     &*x*y*z)*c*x*y-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*b&
     &*z+4*a**2*sin(c*x*y*z)**2*c**2*x**2*y**2*exp(-b*(x**2+y**2+z**2))*&
     &*2*z-y**2*a*cos(c*x*y*z)*c**2*x**2*z*exp(-b*(x**2+y**2+z**2))+0.d0
      s3 = s1+6*y*a*sin(c*x*y*z)*c*x*z**2*b*exp(-b*(x**2+y**2+z**2))-2*a&
     &*sin(c*x*y*z)*c*x*y*exp(-b*(x**2+y**2+z**2))-6*a**2*sin(c*x*y*z)**&
     &2*c**2*x**2*y**4*b*exp(-b*(x**2+y**2+z**2))**2*z-4*a*cos(c*x*y*z)*&
     &b*z*exp(-b*(x**2+y**2+z**2))+0.d0
      s2 = s3+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x&
     &*y*z)*c*x*z**2*b-a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*&
     &y**2*c**2*x**2*z+4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2&
     &))**2*sin(c*x*y*z)*c*y+4*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x*&
     &*2+y**2+z**2))**2*z+2*y**3*a**2*sin(c*x*y*z)*c**3*x**3*z**2*exp(-b&
     &*(x**2+y**2+z**2))**2*cos(c*x*y*z)+4*z**3*a**2*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)**2*b**2+0.d0
      s3 = s2-16*y*a**2*sin(c*x*y*z)*c*x*z**4*exp(-b*(x**2+y**2+z**2))**&
     &2*cos(c*x*y*z)*b**2-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2&
     &+y**2+z**2))**2*y**2*z+2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*cos&
     &(c*x*y*z)**2*c**2*x**2*b-z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos&
     &(c*x*y*z)**2*c**2*x**2-2*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)**2*c**2*x**2*b+0.d0
      D223 = s3+2*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2&
     &*c**2*x**2+2*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2&
     &*c**2*z**3*b-2*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*&
     &*2*c**2*z**3*b+2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z&
     &)*c**3*z**2*cos(c*x*y*z)*y+2*a*sin(c*x*y*z)*c*x*y**3*b*exp(-b*(x**&
     &2+y**2+z**2))+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)&
     &*c**3*x**3*cos(c*x*y*z)*y+0.d0
      s2 = 4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*sin(c*&
     &x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))*&
     &*2*c**2*x**2*z**2-16*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x&
     &*y*z)**2*b**3*y-16*a**2*cos(c*x*y*z)**2*b**3*y**3*exp(-b*(x**2+y**&
     &2+z**2))**2*z**2+2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y&
     &*z)*c**3*z*y**2*cos(c*x*y*z)+2*y**2*a**2*sin(c*x*y*z)*c**3*x**3*z*&
     &*3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)-8*y**3*a**2*sin(c*x*y*&
     &z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2*b-2*a**2*cos(c*x*&
     &y*z)**2*exp(-b*(x**2+y**2+z**2))**2*y*c**2*x**2*z**2+0.d0
      s1 = s2-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*b*y-4*a&
     &*sin(c*x*y*z)*c*x*z*exp(-b*(x**2+y**2+z**2))+4*y**2*a*sin(c*x*y*z)&
     &*c*x*z*b*exp(-b*(x**2+y**2+z**2))-4*a*cos(c*x*y*z)*b*y*exp(-b*(x**&
     &2+y**2+z**2))-2*y*a*cos(c*x*y*z)*c**2*x**2*z**2*exp(-b*(x**2+y**2+&
     &z**2))-4*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z&
     &)*c*x*z+4*a**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*y+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y**2*si&
     &n(c*x*y*z)*c*x*z*b-4*a**2*sin(c*x*y*z)**2*c**2*x**4*z**2*b*exp(-b*&
     &(x**2+y**2+z**2))**2*y+0.d0
      s2 = s1+24*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+y**2+z**2))**2&
     &*z**2-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2+y**2+z**2))**&
     &2*z**2*y+12*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z**2))**2&
     &*cos(c*x*y*z)*b-16*a**2*cos(c*x*y*z)*b**2*y**4*exp(-b*(x**2+y**2+z&
     &**2))**2*sin(c*x*y*z)*c*x*z+4*y*a**2*cos(c*x*y*z)**2*c**2*x**2*z**&
     &4*exp(-b*(x**2+y**2+z**2))**2*b-24*y**2*a**2*sin(c*x*y*z)*c*x*z**3&
     &*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b**2-4*y*a**2*sin(c*x*y*&
     &z)**2*c**2*x**2*z**4*exp(-b*(x**2+y**2+z**2))**2*b-8*a**2*cos(c*x*&
     &y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2*z**3*sin(c*x*y*z)*c+0.d0
      s3 = s2+4*a**2*cos(c*x*y*z)**2*b*x**4*exp(-b*(x**2+y**2+z**2))**2*&
     &y*c**2*z**2-16*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2)&
     &)**2*y**2*sin(c*x*y*z)*c*z-8*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c*x*cos(c*x*y*z)*b**2+4*y**3*a**2*exp(-b*(x**2+y**2+&
     &z**2))**2*sin(c*x*y*z)**2*c**2*x**2+0.d0
      D233 = s3+8*y*a*b**2*z**2*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)+2&
     &*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y+4*z*&
     &*3*a*b*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*x-4*y**5*a**2*exp(-&
     &b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*x**2*b-4*x**4*a**2*exp&
     &(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y**3*b+2*y**4*a**2*e&
     &xp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*x**3*z*cos(c*x*y*z)+0.d0
      s2 = -16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2+y**2+z**2))**&
     &2*z**3+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**&
     &2*z-8*a**2*cos(c*x*y*z)*b**2*x**5*exp(-b*(x**2+y**2+z**2))**2*sin(&
     &c*x*y*z)*c*y+4*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2)&
     &)**2*c**2*x**2*z-8*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z&
     &**2))**2*y**3*sin(c*x*y*z)*c+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x*&
     &*2+y**2+z**2))**2*sin(c*x*y*z)*c*x-16*a**2*cos(c*x*y*z)*b**2*y**3*&
     &exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*z**2-4*y**4*a**2*exp(&
     &-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*z**3*b+0.d0
      s3 = s2-24*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2&
     &*y*sin(c*x*y*z)*c*z**2+4*a**2*cos(c*x*y*z)**2*b*x**4*exp(-b*(x**2+&
     &y**2+z**2))**2*y**2*c**2*z-8*y**2*a**2*sin(c*x*y*z)**2*c**2*x**2*z&
     &**3*exp(-b*(x**2+y**2+z**2))**2*b-4*a**2*sin(c*x*y*z)**2*c**2*x**4&
     &*y**2*b*exp(-b*(x**2+y**2+z**2))**2*z+0.d0
      s1 = s3+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*&
     &z**2*cos(c*x*y*z)*x+4*y**2*a**2*cos(c*x*y*z)**2*c**2*x**2*z**3*exp&
     &(-b*(x**2+y**2+z**2))**2*b+8*a*cos(c*x*y*z)*b**2*x**2*z*exp(-b*(x*&
     &*2+y**2+z**2))-4*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c*x*y-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2&
     &*b*z+0.d0
      s2 = s1+4*a**2*sin(c*x*y*z)**2*c**2*x**2*y**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*z-2*y**2*a*cos(c*x*y*z)*c**2*x**2*z*exp(-b*(x**2+y**2+z**2&
     &))+4*a*sin(c*x*y*z)*c*x**3*y*b*exp(-b*(x**2+y**2+z**2))+4*y*a*sin(&
     &c*x*y*z)*c*x*z**2*b*exp(-b*(x**2+y**2+z**2))-16*a**2*cos(c*x*y*z)*&
     &*2*b**3*x**4*exp(-b*(x**2+y**2+z**2))**2*z-4*a*sin(c*x*y*z)*c*x*y*&
     &exp(-b*(x**2+y**2+z**2))-4*a**2*sin(c*x*y*z)**2*c**2*x**2*y**4*b*e&
     &xp(-b*(x**2+y**2+z**2))**2*z-4*a*cos(c*x*y*z)*b*z*exp(-b*(x**2+y**&
     &2+z**2))+0.d0
      s3 = s2+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y*sin(c*x&
     &*y*z)*c*x*z**2*b-2*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**&
     &2*y**2*c**2*x**2*z+12*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z&
     &**2))**2*sin(c*x*y*z)*c*y+24*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b&
     &*(x**2+y**2+z**2))**2*z+0.d0
      D311 = s3+2*y**3*a**2*sin(c*x*y*z)*c**3*x**3*z**2*exp(-b*(x**2+y*&
     &*2+z**2))**2*cos(c*x*y*z)-16*y*a**2*sin(c*x*y*z)*c*x*z**4*exp(-b*(&
     &x**2+y**2+z**2))**2*cos(c*x*y*z)*b**2+4*z**3*a**2*exp(-b*(x**2+y**&
     &2+z**2))**2*sin(c*x*y*z)**2*c**2*y**2-4*z**5*a**2*exp(-b*(x**2+y**&
     &2+z**2))**2*sin(c*x*y*z)**2*c**2*y**2*b+2*z**4*a**2*exp(-b*(x**2+y&
     &**2+z**2))**2*sin(c*x*y*z)*c**3*y**3*cos(c*x*y*z)*x-16*a**2*cos(c*&
     &x*y*z)**2*b**3*x**2*exp(-b*(x**2+y**2+z**2))**2*y**2*z+0.d0
      s3 = -4*a**2*sin(c*x*y*z)**2*c**2*x**3*y**3*b*exp(-b*(x**2+y**2+z*&
     &*2))**2*z-6*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*&
     &c**2*z**3*x*b-a**2*sin(c*x*y*z)*c*x**2*exp(-b*(x**2+y**2+z**2))**2&
     &*cos(c*x*y*z)-y**2*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c&
     &*cos(c*x*y*z)-8*a**2*sin(c*x*y*z)*c*x**2*z**4*exp(-b*(x**2+y**2+z*&
     &*2))**2*cos(c*x*y*z)*b**2+0.d0
      s2 = s3-2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*y*b*z+4*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))**2&
     &*c**2*x**3*z+2*a**2*sin(c*x*y*z)*c**3*x**4*z**2*exp(-b*(x**2+y**2+&
     &z**2))**2*y**2*cos(c*x*y*z)-6*a**2*sin(c*x*y*z)**2*c**2*x**3*z**3*&
     &exp(-b*(x**2+y**2+z**2))**2*y*b+12*x**2*a**2*cos(c*x*y*z)*b*y**2*e&
     &xp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c+3*a**2*sin(c*x*y*z)**2*c&
     &**2*x**3*z*exp(-b*(x**2+y**2+z**2))**2*y+0.d0
      s3 = s2+24*x*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+y**2+z**2))*&
     &*2*z-24*a**2*cos(c*x*y*z)*b**2*y**2*exp(-b*(x**2+y**2+z**2))**2*si&
     &n(c*x*y*z)*c*x**2*z**2-8*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos&
     &(c*x*y*z)*b**2*sin(c*x*y*z)*c*z**2-8*x**4*a**2*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)*b**2*y**2*sin(c*x*y*z)*c-y**2*a*exp(-b*(x**2&
     &+y**2+z**2))*sin(c*x*y*z)*c+0.d0
      s1 = s3+2*y**2*a*b*z**2*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c-y*&
     &*3*a*exp(-b*(x**2+y**2+z**2))*cos(c*x*y*z)*c**2*x*z+4*z**3*a**2*ex&
     &p(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y*x+10*a**2*sin(c*x&
     &*y*z)*c*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b-2*y**&
     &5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*z*x*b-8*y*&
     &*4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b**2*x**2*sin(c*x&
     &*y*z)*c+0.d0
      s3 = s1+2*a*sin(c*x*y*z)*c*x**2*z**2*b*exp(-b*(x**2+y**2+z**2))+4*&
     &x**2*a*sin(c*x*y*z)*c*y**2*b*exp(-b*(x**2+y**2+z**2))+8*x*a*cos(c*&
     &x*y*z)*b**2*y*z*exp(-b*(x**2+y**2+z**2))-8*z**4*a**2*exp(-b*(x**2+&
     &y**2+z**2))**2*sin(c*x*y*z)*c*y**2*cos(c*x*y*z)*b**2-4*z**5*a**2*e&
     &xp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y*x*b+0.d0
      s2 = s3+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*&
     &y**2*x**2*cos(c*x*y*z)-8*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c*z**2*cos(c*x*y*z)*b**2+2*y**5*a**2*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)**2*c**2*x*z*b-16*y**3*a**2*exp(-b*(x**2+y**2&
     &+z**2))**2*cos(c*x*y*z)**2*b**3*x*z+2*y**4*a**2*exp(-b*(x**2+y**2+&
     &z**2))**2*sin(c*x*y*z)*c*cos(c*x*y*z)*b+2*z**3*a**2*exp(-b*(x**2+y&
     &**2+z**2))**2*cos(c*x*y*z)**2*c**2*x*y**3*b+0.d0
      s3 = s2+2*a**2*cos(c*x*y*z)**2*c**2*x**3*y*z**3*exp(-b*(x**2+y**2+&
     &z**2))**2*b-16*x*a**2*cos(c*x*y*z)**2*b**3*y*exp(-b*(x**2+y**2+z**&
     &2))**2*z**3-y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*&
     &c**2*x*z+3*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c&
     &**2*z*x+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*&
     &z**2*x**2*cos(c*x*y*z)+0.d0
      D312 = s3+10*y**2*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c&
     &*z**2*cos(c*x*y*z)*b-16*x**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(&
     &c*x*y*z)**2*b**3*y*z+2*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c&
     &*x*y*z)*b*sin(c*x*y*z)*c+2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*c&
     &os(c*x*y*z)**2*b*c**2*y*z-a**2*cos(c*x*y*z)**2*c**2*x**3*y*z*exp(-&
     &b*(x**2+y**2+z**2))**2-a*cos(c*x*y*z)*c**2*x**3*y*z*exp(-b*(x**2+y&
     &**2+z**2))-a*sin(c*x*y*z)*c*x**2*exp(-b*(x**2+y**2+z**2))+0.d0
      s3 = 2*a**2*sin(c*x*y*z)**2*c**2*x**3*y**2*exp(-b*(x**2+y**2+z**2)&
     &)**2+2*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*c**&
     &2*z**2*x+4*a**2*cos(c*x*y*z)**2*b*x**3*exp(-b*(x**2+y**2+z**2))**2&
     &*y**2*c**2*z**2-4*a*cos(c*x*y*z)*b*x*exp(-b*(x**2+y**2+z**2))-y**2&
     &*a*cos(c*x*y*z)*c**2*z**2*x*exp(-b*(x**2+y**2+z**2))+0.d0
      s2 = s3+6*y*a*sin(c*x*y*z)*c*x**2*z*b*exp(-b*(x**2+y**2+z**2))-8*y&
     &**2*a**2*sin(c*x*y*z)**2*c**2*x**3*z**2*exp(-b*(x**2+y**2+z**2))**&
     &2*b-2*a*sin(c*x*y*z)*c*y*z*exp(-b*(x**2+y**2+z**2))+12*y*a**2*sin(&
     &c*x*y*z)*c*z**3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b+2*y**2*&
     &a**2*cos(c*x*y*z)**2*c**2*z**4*x*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      s3 = s2-24*y*a**2*sin(c*x*y*z)*c*x**2*z**3*exp(-b*(x**2+y**2+z**2)&
     &)**2*cos(c*x*y*z)*b**2+28*a**2*cos(c*x*y*z)**2*b**2*x*exp(-b*(x**2&
     &+y**2+z**2))**2*z**2-16*a**2*cos(c*x*y*z)**2*b**3*x**3*exp(-b*(x**&
     &2+y**2+z**2))**2*z**2+4*a**2*cos(c*x*y*z)**2*b**2*x*exp(-b*(x**2+y&
     &**2+z**2))**2*y**2+0.d0
      s1 = s3-6*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c*&
     &*2*z**2*b*x+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c*z-16*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**2&
     &+z**2))**2*sin(c*x*y*z)*c*x**2*z+4*a**2*cos(c*x*y*z)**2*b**2*x**3*&
     &exp(-b*(x**2+y**2+z**2))**2+2*y**3*a**2*sin(c*x*y*z)*c**3*x**2*z**&
     &3*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)+20*a**2*cos(c*x*y*z)*ex&
     &p(-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)*c*x**2*z*b+0.d0
      s3 = s1-a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*y**2*c**2&
     &*z**2*x-16*a**2*cos(c*x*y*z)*b**2*x**4*exp(-b*(x**2+y**2+z**2))**2&
     &*y*sin(c*x*y*z)*c*z-2*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**&
     &2*sin(c*x*y*z)*c*y*z-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2&
     &))**2*b*x+0.d0
      s2 = s3+4*a**2*sin(c*x*y*z)**2*c**2*y**2*z**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*x-8*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*y&
     &**3*cos(c*x*y*z)*b**2-16*x*a**2*cos(c*x*y*z)**2*b**3*y**2*exp(-b*(&
     &x**2+y**2+z**2))**2*z**2+y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)**2*c**2*x+2*x**4*y**3*a**2*sin(c*x*y*z)*c**3*exp(-b*(x**&
     &2+y**2+z**2))**2*z*cos(c*x*y*z)+8*x*a*cos(c*x*y*z)*b**2*z**2*exp(-&
     &b*(x**2+y**2+z**2))+0.d0
      s3 = s2-2*x**3*y**4*a**2*sin(c*x*y*z)**2*c**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*b+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**&
     &3*z*x**2*cos(c*x*y*z)-16*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos&
     &(c*x*y*z)**2*b**3*x-8*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*&
     &x*y*z)*c*y*cos(c*x*y*z)*b**2-2*x**5*a**2*exp(-b*(x**2+y**2+z**2))*&
     &*2*sin(c*x*y*z)**2*c**2*y**2*b+0.d0
      D313 = s3+2*z**3*a*b*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*y-a*&
     &cos(c*x*y*z)*c**2*x**3*y**2*exp(-b*(x**2+y**2+z**2))-a**2*cos(c*x*&
     &y*z)**2*c**2*x**3*y**2*exp(-b*(x**2+y**2+z**2))**2+2*x**5*a**2*exp&
     &(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z)**2*b*c**2*y**2+2*x**3*y**4*a&
     &**2*cos(c*x*y*z)**2*c**2*exp(-b*(x**2+y**2+z**2))**2*b-6*z**4*a**2&
     &*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y**2*b*x+0.d0
      s2 = -16*a**2*cos(c*x*y*z)**2*b**3*y**2*exp(-b*(x**2+y**2+z**2))**&
     &2*z**3-16*a**2*cos(c*x*y*z)**2*b**3*y**4*exp(-b*(x**2+y**2+z**2))*&
     &*2*z-8*a**2*cos(c*x*y*z)*b**2*y**5*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c*x+2*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)&
     &**2*c**2*z+8*a*cos(c*x*y*z)*b**2*y**2*z*exp(-b*(x**2+y**2+z**2))+4&
     &*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*c**2*x**2&
     &*z-8*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2*y**3*&
     &sin(c*x*y*z)*c+12*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2)&
     &)**2*sin(c*x*y*z)*c*x+0.d0
      s3 = s2-24*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c*x*z**2-16*a**2*cos(c*x*y*z)*b**2*x**3*exp(-b*(x**2&
     &+y**2+z**2))**2*y*sin(c*x*y*z)*c*z**2+4*a**2*cos(c*x*y*z)**2*b*x**&
     &4*exp(-b*(x**2+y**2+z**2))**2*y**2*c**2*z-8*y**2*a**2*sin(c*x*y*z)&
     &**2*c**2*x**2*z**3*exp(-b*(x**2+y**2+z**2))**2*b+0.d0
      s1 = s3-4*a**2*sin(c*x*y*z)**2*c**2*x**4*y**2*b*exp(-b*(x**2+y**2+&
     &z**2))**2*z+4*y**2*a**2*cos(c*x*y*z)**2*c**2*x**2*z**3*exp(-b*(x**&
     &2+y**2+z**2))**2*b+24*y**2*a**2*cos(c*x*y*z)**2*b**2*z*exp(-b*(x**&
     &2+y**2+z**2))**2-4*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*s&
     &in(c*x*y*z)*c*x*y-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))*&
     &*2*b*z+0.d0
      s2 = s1+4*a**2*sin(c*x*y*z)**2*c**2*x**2*y**2*exp(-b*(x**2+y**2+z*&
     &*2))**2*z-2*y**2*a*cos(c*x*y*z)*c**2*x**2*z*exp(-b*(x**2+y**2+z**2&
     &))+4*y*a*sin(c*x*y*z)*c*x*z**2*b*exp(-b*(x**2+y**2+z**2))-4*a*sin(&
     &c*x*y*z)*c*x*y*exp(-b*(x**2+y**2+z**2))-4*a**2*sin(c*x*y*z)**2*c**&
     &2*x**2*y**4*b*exp(-b*(x**2+y**2+z**2))**2*z-4*a*cos(c*x*y*z)*b*z*e&
     &xp(-b*(x**2+y**2+z**2))+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**&
     &2))**2*y*sin(c*x*y*z)*c*x*z**2*b-2*a**2*cos(c*x*y*z)**2*exp(-b*(x*&
     &*2+y**2+z**2))**2*y**2*c**2*x**2*z+0.d0
      s3 = s2+4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*sin&
     &(c*x*y*z)*c*y+2*y**3*a**2*sin(c*x*y*z)*c**3*x**3*z**2*exp(-b*(x**2&
     &+y**2+z**2))**2*cos(c*x*y*z)-16*y*a**2*sin(c*x*y*z)*c*x*z**4*exp(-&
     &b*(x**2+y**2+z**2))**2*cos(c*x*y*z)*b**2-16*a**2*cos(c*x*y*z)**2*b&
     &**3*x**2*exp(-b*(x**2+y**2+z**2))**2*y**2*z+0.d0
      D322 = s3-4*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2&
     &*c**2*x**2*b+4*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*&
     &*2*c**2*x**2-4*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*&
     &*2*c**2*z**3*b+2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z&
     &)*c**3*z**2*cos(c*x*y*z)*y+4*a*sin(c*x*y*z)*c*x*y**3*b*exp(-b*(x**&
     &2+y**2+z**2))+2*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)&
     &*c**3*x**3*cos(c*x*y*z)*y+0.d0
      s3 = 4*a**2*cos(c*x*y*z)*b*x**3*exp(-b*(x**2+y**2+z**2))**2*sin(c*&
     &x*y*z)*c*z+4*a**2*cos(c*x*y*z)**2*b*y**3*exp(-b*(x**2+y**2+z**2))*&
     &*2*c**2*x**2*z**2-16*z**4*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x&
     &*y*z)**2*b**3*y-16*a**2*cos(c*x*y*z)**2*b**3*y**3*exp(-b*(x**2+y**&
     &2+z**2))**2*z**2-y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z&
     &)**2*c**2*x**2+0.d0
      s2 = s3+2*x**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c**3*&
     &z*y**2*cos(c*x*y*z)+4*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x**2+&
     &y**2+z**2))**2*y+2*y**2*a**2*sin(c*x*y*z)*c**3*x**3*z**3*exp(-b*(x&
     &**2+y**2+z**2))**2*cos(c*x*y*z)-8*y**3*a**2*sin(c*x*y*z)**2*c**2*x&
     &**2*z**2*exp(-b*(x**2+y**2+z**2))**2*b-a**2*cos(c*x*y*z)**2*exp(-b&
     &*(x**2+y**2+z**2))**2*y*c**2*x**2*z**2+0.d0
      s3 = s2-4*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+z**2))**2*b*y-2*a&
     &*sin(c*x*y*z)*c*x*z*exp(-b*(x**2+y**2+z**2))+6*y**2*a*sin(c*x*y*z)&
     &*c*x*z*b*exp(-b*(x**2+y**2+z**2))-4*a*cos(c*x*y*z)*b*y*exp(-b*(x**&
     &2+y**2+z**2))+0.d0
      s1 = s3-y*a*cos(c*x*y*z)*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))-2&
     &*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*z+&
     &4*a**2*sin(c*x*y*z)**2*c**2*x**2*z**2*exp(-b*(x**2+y**2+z**2))**2*&
     &y+20*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2))**2*y**2*sin(c*x*y*&
     &z)*c*x*z*b-6*a**2*sin(c*x*y*z)**2*c**2*x**4*z**2*b*exp(-b*(x**2+y*&
     &*2+z**2))**2*y+4*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*cos(c*x*y*z&
     &)**2*b**2+0.d0
      s3 = s1+28*a**2*cos(c*x*y*z)**2*b**2*y*exp(-b*(x**2+y**2+z**2))**2&
     &*z**2-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2+y**2+z**2))**&
     &2*z**2*y+12*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z**2))**2&
     &*cos(c*x*y*z)*b-16*a**2*cos(c*x*y*z)*b**2*y**4*exp(-b*(x**2+y**2+z&
     &**2))**2*sin(c*x*y*z)*c*x*z+0.d0
      s2 = s3+2*y*a**2*cos(c*x*y*z)**2*c**2*x**2*z**4*exp(-b*(x**2+y**2+&
     &z**2))**2*b-24*y**2*a**2*sin(c*x*y*z)*c*x*z**3*exp(-b*(x**2+y**2+z&
     &**2))**2*cos(c*x*y*z)*b**2-6*y*a**2*sin(c*x*y*z)**2*c**2*x**2*z**4&
     &*exp(-b*(x**2+y**2+z**2))**2*b-8*a**2*cos(c*x*y*z)*b**2*x**3*exp(-&
     &b*(x**2+y**2+z**2))**2*z**3*sin(c*x*y*z)*c+2*a**2*cos(c*x*y*z)**2*&
     &b*x**4*exp(-b*(x**2+y**2+z**2))**2*y*c**2*z**2-16*a**2*cos(c*x*y*z&
     &)*b**2*x**3*exp(-b*(x**2+y**2+z**2))**2*y**2*sin(c*x*y*z)*c*z+0.d0
      s3 = s2-8*z**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*x*c&
     &os(c*x*y*z)*b**2+2*y**3*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y&
     &*z)**2*c**2*x**2+8*y*a*b**2*z**2*exp(-b*(x**2+y**2+z**2))*cos(c*x*&
     &y*z)+x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2*c**2*y+&
     &2*z**3*a*b*exp(-b*(x**2+y**2+z**2))*sin(c*x*y*z)*c*x+0.d0
      D323 = s3-2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)**2&
     &*c**2*x**2*b-2*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*&
     &*2*c**2*y**3*b+2*y**4*a**2*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z&
     &)*c**3*x**3*z*cos(c*x*y*z)+2*x**4*a**2*exp(-b*(x**2+y**2+z**2))**2&
     &*cos(c*x*y*z)**2*b*y**3*c**2+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))*&
     &*2*cos(c*x*y*z)**2*b*c**2*x**2-y**3*a*exp(-b*(x**2+y**2+z**2))*cos&
     &(c*x*y*z)*c**2*x**2+0.d0
      s2 = -16*a**2*cos(c*x*y*z)**2*b**3*y**2*exp(-b*(x**2+y**2+z**2))**&
     &2*z**3-16*a**2*cos(c*x*y*z)**2*b**3*x**2*exp(-b*(x**2+y**2+z**2))*&
     &*2*z**3+4*a**2*cos(c*x*y*z)**2*b*y**4*exp(-b*(x**2+y**2+z**2))**2*&
     &c**2*x**2*z+4*a**2*cos(c*x*y*z)*b*y**3*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c*x-24*a**2*cos(c*x*y*z)*b**2*y**3*exp(-b*(x**2+y**2&
     &+z**2))**2*sin(c*x*y*z)*c*x*z**2-24*a**2*cos(c*x*y*z)*b**2*x**3*ex&
     &p(-b*(x**2+y**2+z**2))**2*y*sin(c*x*y*z)*c*z**2+4*a**2*cos(c*x*y*z&
     &)**2*b*x**4*exp(-b*(x**2+y**2+z**2))**2*y**2*c**2*z+0.d0
      s1 = s2-8*y**2*a**2*sin(c*x*y*z)**2*c**2*x**2*z**3*exp(-b*(x**2+y*&
     &*2+z**2))**2*b-8*a**2*sin(c*x*y*z)**2*c**2*x**4*y**2*b*exp(-b*(x**&
     &2+y**2+z**2))**2*z+4*y**2*a**2*cos(c*x*y*z)**2*c**2*x**2*z**3*exp(&
     &-b*(x**2+y**2+z**2))**2*b+8*y**2*a**2*cos(c*x*y*z)**2*b**2*z*exp(-&
     &b*(x**2+y**2+z**2))**2-4*a**2*cos(c*x*y*z)*exp(-b*(x**2+y**2+z**2)&
     &)**2*sin(c*x*y*z)*c*x*y-12*a**2*cos(c*x*y*z)**2*exp(-b*(x**2+y**2+&
     &z**2))**2*b*z+4*a**2*sin(c*x*y*z)**2*c**2*x**2*y**2*exp(-b*(x**2+y&
     &**2+z**2))**2*z-2*y**2*a*cos(c*x*y*z)*c**2*x**2*z*exp(-b*(x**2+y**&
     &2+z**2))+0.d0
      s2 = s1+8*y*a*sin(c*x*y*z)*c*x*z**2*b*exp(-b*(x**2+y**2+z**2))-4*a&
     &*sin(c*x*y*z)*c*x*y*exp(-b*(x**2+y**2+z**2))-8*a**2*sin(c*x*y*z)**&
     &2*c**2*x**2*y**4*b*exp(-b*(x**2+y**2+z**2))**2*z-12*a*cos(c*x*y*z)&
     &*b*z*exp(-b*(x**2+y**2+z**2))+28*a**2*cos(c*x*y*z)*exp(-b*(x**2+y*&
     &*2+z**2))**2*y*sin(c*x*y*z)*c*x*z**2*b-2*a**2*cos(c*x*y*z)**2*exp(&
     &-b*(x**2+y**2+z**2))**2*y**2*c**2*x**2*z+4*a**2*cos(c*x*y*z)*b*x**&
     &3*exp(-b*(x**2+y**2+z**2))**2*sin(c*x*y*z)*c*y+0.d0
      D333 = s2+8*a**2*cos(c*x*y*z)**2*b**2*x**2*exp(-b*(x**2+y**2+z**2&
     &))**2*z+2*y**3*a**2*sin(c*x*y*z)*c**3*x**3*z**2*exp(-b*(x**2+y**2+&
     &z**2))**2*cos(c*x*y*z)+32*z**3*a**2*exp(-b*(x**2+y**2+z**2))**2*co&
     &s(c*x*y*z)**2*b**2-24*y*a**2*sin(c*x*y*z)*c*x*z**4*exp(-b*(x**2+y*&
     &*2+z**2))**2*cos(c*x*y*z)*b**2+2*x**5*a**2*exp(-b*(x**2+y**2+z**2)&
     &)**2*sin(c*x*y*z)*c**3*y**3*cos(c*x*y*z)+8*a*cos(c*x*y*z)*b**2*z**&
     &3*exp(-b*(x**2+y**2+z**2))+2*y**5*a**2*exp(-b*(x**2+y**2+z**2))**2&
     &*sin(c*x*y*z)*c**3*x**3*cos(c*x*y*z)-16*a**2*cos(c*x*y*z)**2*b**3*&
     &z**5*exp(-b*(x**2+y**2+z**2))**2
