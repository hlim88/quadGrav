
      s1 = 1+2*a*cos(c*x*y)*exp(-b*(x**2+y**2))-2*y*a*sin(c*x*y)*c*x*exp&
     &(-b*(x**2+y**2))-4*a*cos(c*x*y)*b*x**2*exp(-b*(x**2+y**2))+a**2*co&
     &s(c*x*y)**2*exp(-b*(x**2+y**2))**2-2*a**2*cos(c*x*y)*exp(-b*(x**2+&
     &y**2))**2*y*sin(c*x*y)*c*x
      g11 = s1-4*a**2*cos(c*x*y)**2*b*x**2*exp(-b*(x**2+y**2))**2+y**2*a&
     &**2*sin(c*x*y)**2*c**2*x**2*exp(-b*(x**2+y**2))**2+4*a**2*cos(c*x*&
     &y)*b*x**3*exp(-b*(x**2+y**2))**2*y*sin(c*x*y)*c+4*a**2*cos(c*x*y)*&
     &*2*b**2*x**4*exp(-b*(x**2+y**2))**2+y**4*a**2*exp(-b*(x**2+y**2))*&
     &*2*sin(c*x*y)**2*c**2+4*y**3*a**2*sin(c*x*y)*c*x*exp(-b*(x**2+y**2&
     &))**2*cos(c*x*y)*b+4*a**2*cos(c*x*y)**2*b**2*x**2*exp(-b*(x**2+y**&
     &2))**2*y**2
      s1 = -x**2*a*exp(-b*(x**2+y**2))*sin(c*x*y)*c-4*x*a*exp(-b*(x**2+y&
     &**2))*cos(c*x*y)*b*y-x**2*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)*s&
     &in(c*x*y)*c-4*x*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)**2*b*y+x**3&
     &*a**2*exp(-b*(x**2+y**2))**2*y*sin(c*x*y)**2*c**2+4*x**2*a**2*exp(&
     &-b*(x**2+y**2))**2*y**2*sin(c*x*y)*c*cos(c*x*y)*b
      g12 = s1+2*x**4*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)*b*sin(c*x*y&
     &)*c+4*x**3*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)**2*b**2*y-y**2*a&
     &*exp(-b*(x**2+y**2))*sin(c*x*y)*c-y**2*a**2*exp(-b*(x**2+y**2))**2&
     &*sin(c*x*y)*c*cos(c*x*y)+y**3*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*&
     &y)**2*c**2*x+2*y**4*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c*cos(c&
     &*x*y)*b+4*y**3*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)**2*b**2*x
      g13 = 0.0
      s1 = 1+x**4*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2+4*a**2*&
     &cos(c*x*y)*b*x**3*exp(-b*(x**2+y**2))**2*y*sin(c*x*y)*c+4*a**2*cos&
     &(c*x*y)**2*b**2*x**2*exp(-b*(x**2+y**2))**2*y**2+2*a*cos(c*x*y)*ex&
     &p(-b*(x**2+y**2))-2*y*a*sin(c*x*y)*c*x*exp(-b*(x**2+y**2))
      g22 = s1-4*a*cos(c*x*y)*b*y**2*exp(-b*(x**2+y**2))+a**2*cos(c*x*y)&
     &**2*exp(-b*(x**2+y**2))**2-2*a**2*cos(c*x*y)*exp(-b*(x**2+y**2))**&
     &2*y*sin(c*x*y)*c*x-4*a**2*cos(c*x*y)**2*exp(-b*(x**2+y**2))**2*b*y&
     &**2+y**2*a**2*sin(c*x*y)**2*c**2*x**2*exp(-b*(x**2+y**2))**2+4*y**&
     &3*a**2*sin(c*x*y)*c*x*exp(-b*(x**2+y**2))**2*cos(c*x*y)*b+4*a**2*c&
     &os(c*x*y)**2*b**2*y**4*exp(-b*(x**2+y**2))**2
      g23 = 0.0
      g33 = 1.0
      s1 = -24*y**3*a**2*sin(c*x*y)*c*x**2*exp(-b*(x**2+y**2))**2*cos(c*&
     &x*y)*b**2-4*a*sin(c*x*y)*c*y*exp(-b*(x**2+y**2))-12*a*cos(c*x*y)*b&
     &*x*exp(-b*(x**2+y**2))-4*a**2*cos(c*x*y)*exp(-b*(x**2+y**2))**2*si&
     &n(c*x*y)*c*y-2*y**2*a*cos(c*x*y)*c**2*x*exp(-b*(x**2+y**2))+8*y*a*&
     &sin(c*x*y)*c*x**2*b*exp(-b*(x**2+y**2))+8*a*cos(c*x*y)*b**2*x**3*e&
     &xp(-b*(x**2+y**2))-12*a**2*cos(c*x*y)**2*exp(-b*(x**2+y**2))**2*b*&
     &x+4*a**2*sin(c*x*y)**2*c**2*y**2*exp(-b*(x**2+y**2))**2*x+28*a**2*&
     &cos(c*x*y)*exp(-b*(x**2+y**2))**2*y*sin(c*x*y)*c*x**2*b-2*a**2*cos&
     &(c*x*y)**2*exp(-b*(x**2+y**2))**2*y**2*c**2*x
      s2 = s1+32*a**2*cos(c*x*y)**2*b**2*x**3*exp(-b*(x**2+y**2))**2+2*y&
     &**3*a**2*sin(c*x*y)*c**3*x**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)-8*&
     &y**2*a**2*sin(c*x*y)**2*c**2*x**3*exp(-b*(x**2+y**2))**2*b-24*a**2&
     &*cos(c*x*y)*b**2*x**4*exp(-b*(x**2+y**2))**2*y*sin(c*x*y)*c+4*a**2&
     &*cos(c*x*y)**2*b*x**3*exp(-b*(x**2+y**2))**2*y**2*c**2
      D111 = s2-16*a**2*cos(c*x*y)**2*b**3*x**5*exp(-b*(x**2+y**2))**2-&
     &8*y**4*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2*b*x+2*y**5*a&
     &**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c**3*cos(c*x*y)+4*y**3*a**2*&
     &sin(c*x*y)*c*exp(-b*(x**2+y**2))**2*cos(c*x*y)*b+4*y**4*a**2*cos(c&
     &*x*y)**2*c**2*x*exp(-b*(x**2+y**2))**2*b+8*a**2*cos(c*x*y)**2*b**2&
     &*x*exp(-b*(x**2+y**2))**2*y**2-16*a**2*cos(c*x*y)**2*b**3*x**3*exp&
     &(-b*(x**2+y**2))**2*y**2
      s2 = 28*a**2*cos(c*x*y)**2*b**2*x**2*exp(-b*(x**2+y**2))**2*y-2*y*&
     &*5*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2*b-y**3*a*exp(-b*&
     &(x**2+y**2))*cos(c*x*y)*c**2-a**2*cos(c*x*y)**2*exp(-b*(x**2+y**2)&
     &)**2*y*c**2*x**2+4*y**3*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)**2*&
     &b**2-4*a**2*cos(c*x*y)**2*exp(-b*(x**2+y**2))**2*b*y-24*a**2*cos(c&
     &*x*y)*b**2*x**3*exp(-b*(x**2+y**2))**2*y**2*sin(c*x*y)*c
      s1 = s2+20*a**2*cos(c*x*y)*exp(-b*(x**2+y**2))**2*y**2*sin(c*x*y)*&
     &c*x*b+2*a*sin(c*x*y)*c*x**3*b*exp(-b*(x**2+y**2))-2*a*sin(c*x*y)*c&
     &*x*exp(-b*(x**2+y**2))-4*a*cos(c*x*y)*b*y*exp(-b*(x**2+y**2))+2*y*&
     &*5*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)**2*c**2*b+12*a**2*cos(c*&
     &x*y)*b*x**3*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c-16*a**2*cos(c*x*y)&
     &**2*b**3*x**2*exp(-b*(x**2+y**2))**2*y**3-y*a*cos(c*x*y)*c**2*x**2&
     &*exp(-b*(x**2+y**2))
      s2 = s1+6*y**2*a*sin(c*x*y)*c*x*b*exp(-b*(x**2+y**2))+8*a*cos(c*x*&
     &y)*b**2*x**2*y*exp(-b*(x**2+y**2))-2*a**2*cos(c*x*y)*exp(-b*(x**2+&
     &y**2))**2*sin(c*x*y)*c*x+4*a**2*sin(c*x*y)**2*c**2*x**2*exp(-b*(x*&
     &*2+y**2))**2*y+2*y**2*a**2*sin(c*x*y)*c**3*x**3*exp(-b*(x**2+y**2)&
     &)**2*cos(c*x*y)-8*y**3*a**2*sin(c*x*y)**2*c**2*x**2*exp(-b*(x**2+y&
     &**2))**2*b-6*a**2*sin(c*x*y)**2*c**2*x**4*b*exp(-b*(x**2+y**2))**2&
     &*y
      D112 = s2+2*a**2*cos(c*x*y)**2*b*x**4*exp(-b*(x**2+y**2))**2*y*c*&
     &*2-8*a**2*cos(c*x*y)*b**2*x**5*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c&
     &-16*a**2*cos(c*x*y)**2*b**3*x**4*exp(-b*(x**2+y**2))**2*y+2*y**3*a&
     &**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2+2*y**4*a**2*exp(-b*(&
     &x**2+y**2))**2*sin(c*x*y)*c**3*cos(c*x*y)*x+4*y**3*a**2*cos(c*x*y)&
     &**2*c**2*x**2*exp(-b*(x**2+y**2))**2*b-16*y**4*a**2*sin(c*x*y)*c*x&
     &*exp(-b*(x**2+y**2))**2*cos(c*x*y)*b**2-y**3*a**2*exp(-b*(x**2+y**&
     &2))**2*cos(c*x*y)**2*c**2
      D113 = 0.0
      s2 = 4*y**3*a*b*exp(-b*(x**2+y**2))*sin(c*x*y)*c-4*a**2*cos(c*x*y)&
     &*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c*y-16*y**4*a**2*exp(-b*(x**2+y&
     &**2))**2*cos(c*x*y)**2*b**3*x-16*a**2*cos(c*x*y)**2*b**3*x**3*exp(&
     &-b*(x**2+y**2))**2*y**2+8*x*a*b**2*y**2*exp(-b*(x**2+y**2))*cos(c*&
     &x*y)-24*y**3*a**2*sin(c*x*y)*c*x**2*exp(-b*(x**2+y**2))**2*cos(c*x&
     &*y)*b**2
      s1 = s2-16*a**2*cos(c*x*y)*b**2*x**4*exp(-b*(x**2+y**2))**2*y*sin(&
     &c*x*y)*c-2*y**2*a*cos(c*x*y)*c**2*x*exp(-b*(x**2+y**2))-4*x**5*a**&
     &2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2*b-4*a*sin(c*x*y)*c*y*e&
     &xp(-b*(x**2+y**2))-4*a*cos(c*x*y)*b*x*exp(-b*(x**2+y**2))+4*y*a*si&
     &n(c*x*y)*c*x**2*b*exp(-b*(x**2+y**2))-4*a**2*cos(c*x*y)**2*exp(-b*&
     &(x**2+y**2))**2*b*x
      s2 = s1+4*a**2*sin(c*x*y)**2*c**2*y**2*exp(-b*(x**2+y**2))**2*x+20&
     &*a**2*cos(c*x*y)*exp(-b*(x**2+y**2))**2*y*sin(c*x*y)*c*x**2*b-2*a*&
     &*2*cos(c*x*y)**2*exp(-b*(x**2+y**2))**2*y**2*c**2*x+2*y**3*a**2*si&
     &n(c*x*y)*c**3*x**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)-8*y**2*a**2*s&
     &in(c*x*y)**2*c**2*x**3*exp(-b*(x**2+y**2))**2*b+4*a**2*cos(c*x*y)*&
     &*2*b*x**3*exp(-b*(x**2+y**2))**2*y**2*c**2
      D122 = s2-4*y**4*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2*b&
     &*x+12*y**3*a**2*sin(c*x*y)*c*exp(-b*(x**2+y**2))**2*cos(c*x*y)*b+4&
     &*y**4*a**2*cos(c*x*y)**2*c**2*x*exp(-b*(x**2+y**2))**2*b+24*a**2*c&
     &os(c*x*y)**2*b**2*x*exp(-b*(x**2+y**2))**2*y**2+4*x**3*a**2*exp(-b&
     &*(x**2+y**2))**2*sin(c*x*y)**2*c**2+2*x**4*a**2*exp(-b*(x**2+y**2)&
     &)**2*y*sin(c*x*y)*c**3*cos(c*x*y)-8*y**5*a**2*exp(-b*(x**2+y**2))*&
     &*2*sin(c*x*y)*c*cos(c*x*y)*b**2
      D123 = 0.0
      D133 = 0.0
      s2 = 24*a**2*cos(c*x*y)**2*b**2*x**2*exp(-b*(x**2+y**2))**2*y-4*y*&
     &*5*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2*b-2*a**2*cos(c*x&
     &*y)**2*exp(-b*(x**2+y**2))**2*y*c**2*x**2-4*a**2*cos(c*x*y)**2*exp&
     &(-b*(x**2+y**2))**2*b*y-24*a**2*cos(c*x*y)*b**2*x**3*exp(-b*(x**2+&
     &y**2))**2*y**2*sin(c*x*y)*c+20*a**2*cos(c*x*y)*exp(-b*(x**2+y**2))&
     &**2*y**2*sin(c*x*y)*c*x*b
      s1 = s2+4*a*sin(c*x*y)*c*x**3*b*exp(-b*(x**2+y**2))-4*a*sin(c*x*y)&
     &*c*x*exp(-b*(x**2+y**2))-4*a*cos(c*x*y)*b*y*exp(-b*(x**2+y**2))+12&
     &*a**2*cos(c*x*y)*b*x**3*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c-16*a**&
     &2*cos(c*x*y)**2*b**3*x**2*exp(-b*(x**2+y**2))**2*y**3-2*y*a*cos(c*&
     &x*y)*c**2*x**2*exp(-b*(x**2+y**2))+4*y**2*a*sin(c*x*y)*c*x*b*exp(-&
     &b*(x**2+y**2))
      s2 = s1+8*a*cos(c*x*y)*b**2*x**2*y*exp(-b*(x**2+y**2))-4*a**2*cos(&
     &c*x*y)*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c*x+4*a**2*sin(c*x*y)**2*&
     &c**2*x**2*exp(-b*(x**2+y**2))**2*y+2*y**2*a**2*sin(c*x*y)*c**3*x**&
     &3*exp(-b*(x**2+y**2))**2*cos(c*x*y)-8*y**3*a**2*sin(c*x*y)**2*c**2&
     &*x**2*exp(-b*(x**2+y**2))**2*b-4*a**2*sin(c*x*y)**2*c**2*x**4*b*ex&
     &p(-b*(x**2+y**2))**2*y
      D211 = s2+4*a**2*cos(c*x*y)**2*b*x**4*exp(-b*(x**2+y**2))**2*y*c*&
     &*2-8*a**2*cos(c*x*y)*b**2*x**5*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c&
     &-16*a**2*cos(c*x*y)**2*b**3*x**4*exp(-b*(x**2+y**2))**2*y+4*y**3*a&
     &**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2+2*y**4*a**2*exp(-b*(&
     &x**2+y**2))**2*sin(c*x*y)*c**3*cos(c*x*y)*x+4*y**3*a**2*cos(c*x*y)&
     &**2*c**2*x**2*exp(-b*(x**2+y**2))**2*b-16*y**4*a**2*sin(c*x*y)*c*x&
     &*exp(-b*(x**2+y**2))**2*cos(c*x*y)*b**2
      s2 = 2*y**3*a*b*exp(-b*(x**2+y**2))*sin(c*x*y)*c-2*a**2*cos(c*x*y)&
     &*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c*y-16*y**4*a**2*exp(-b*(x**2+y&
     &**2))**2*cos(c*x*y)**2*b**3*x-16*a**2*cos(c*x*y)**2*b**3*x**3*exp(&
     &-b*(x**2+y**2))**2*y**2+4*a**2*cos(c*x*y)**2*b**2*x**3*exp(-b*(x**&
     &2+y**2))**2-x**3*a*exp(-b*(x**2+y**2))*cos(c*x*y)*c**2+8*x*a*b**2*&
     &y**2*exp(-b*(x**2+y**2))*cos(c*x*y)
      s1 = s2-24*y**3*a**2*sin(c*x*y)*c*x**2*exp(-b*(x**2+y**2))**2*cos(&
     &c*x*y)*b**2-16*a**2*cos(c*x*y)*b**2*x**4*exp(-b*(x**2+y**2))**2*y*&
     &sin(c*x*y)*c+2*x**5*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)**2*b*c*&
     &*2-y**2*a*cos(c*x*y)*c**2*x*exp(-b*(x**2+y**2))-2*x**5*a**2*exp(-b&
     &*(x**2+y**2))**2*sin(c*x*y)**2*c**2*b-2*a*sin(c*x*y)*c*y*exp(-b*(x&
     &**2+y**2))-4*a*cos(c*x*y)*b*x*exp(-b*(x**2+y**2))+6*y*a*sin(c*x*y)&
     &*c*x**2*b*exp(-b*(x**2+y**2))
      s2 = s1-4*a**2*cos(c*x*y)**2*exp(-b*(x**2+y**2))**2*b*x+4*a**2*sin&
     &(c*x*y)**2*c**2*y**2*exp(-b*(x**2+y**2))**2*x+20*a**2*cos(c*x*y)*e&
     &xp(-b*(x**2+y**2))**2*y*sin(c*x*y)*c*x**2*b-a**2*cos(c*x*y)**2*exp&
     &(-b*(x**2+y**2))**2*y**2*c**2*x+2*y**3*a**2*sin(c*x*y)*c**3*x**2*e&
     &xp(-b*(x**2+y**2))**2*cos(c*x*y)-8*y**2*a**2*sin(c*x*y)**2*c**2*x*&
     &*3*exp(-b*(x**2+y**2))**2*b+4*a**2*cos(c*x*y)**2*b*x**3*exp(-b*(x*&
     &*2+y**2))**2*y**2*c**2
      D212 = s2-6*y**4*a**2*exp(-b*(x**2+y**2))**2*sin(c*x*y)**2*c**2*b&
     &*x+12*y**3*a**2*sin(c*x*y)*c*exp(-b*(x**2+y**2))**2*cos(c*x*y)*b+2&
     &*y**4*a**2*cos(c*x*y)**2*c**2*x*exp(-b*(x**2+y**2))**2*b+28*a**2*c&
     &os(c*x*y)**2*b**2*x*exp(-b*(x**2+y**2))**2*y**2+2*x**3*a**2*exp(-b&
     &*(x**2+y**2))**2*sin(c*x*y)**2*c**2+2*x**4*a**2*exp(-b*(x**2+y**2)&
     &)**2*y*sin(c*x*y)*c**3*cos(c*x*y)-8*y**5*a**2*exp(-b*(x**2+y**2))*&
     &*2*sin(c*x*y)*c*cos(c*x*y)*b**2-x**3*a**2*exp(-b*(x**2+y**2))**2*c&
     &os(c*x*y)**2*c**2
      D213 = 0.0
      s1 = -4*a*sin(c*x*y)*c*x*exp(-b*(x**2+y**2))-12*a*cos(c*x*y)*b*y*e&
     &xp(-b*(x**2+y**2))-2*y*a*cos(c*x*y)*c**2*x**2*exp(-b*(x**2+y**2))+&
     &8*y**2*a*sin(c*x*y)*c*x*b*exp(-b*(x**2+y**2))-4*a**2*cos(c*x*y)*ex&
     &p(-b*(x**2+y**2))**2*sin(c*x*y)*c*x-12*a**2*cos(c*x*y)**2*exp(-b*(&
     &x**2+y**2))**2*b*y+4*a**2*sin(c*x*y)**2*c**2*x**2*exp(-b*(x**2+y**&
     &2))**2*y+28*a**2*cos(c*x*y)*exp(-b*(x**2+y**2))**2*y**2*sin(c*x*y)&
     &*c*x*b-2*a**2*cos(c*x*y)**2*exp(-b*(x**2+y**2))**2*y*c**2*x**2+4*a&
     &**2*cos(c*x*y)*b*x**3*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c+8*a**2*c&
     &os(c*x*y)**2*b**2*x**2*exp(-b*(x**2+y**2))**2*y
      s2 = s1+2*y**2*a**2*sin(c*x*y)*c**3*x**3*exp(-b*(x**2+y**2))**2*co&
     &s(c*x*y)-8*y**3*a**2*sin(c*x*y)**2*c**2*x**2*exp(-b*(x**2+y**2))**&
     &2*b-8*a**2*sin(c*x*y)**2*c**2*x**4*b*exp(-b*(x**2+y**2))**2*y-24*a&
     &**2*cos(c*x*y)*b**2*x**3*exp(-b*(x**2+y**2))**2*y**2*sin(c*x*y)*c+&
     &4*a**2*cos(c*x*y)**2*b*x**4*exp(-b*(x**2+y**2))**2*y*c**2
      D222 = s2+4*y**3*a**2*cos(c*x*y)**2*c**2*x**2*exp(-b*(x**2+y**2))&
     &**2*b-24*y**4*a**2*sin(c*x*y)*c*x*exp(-b*(x**2+y**2))**2*cos(c*x*y&
     &)*b**2-16*a**2*cos(c*x*y)**2*b**3*x**2*exp(-b*(x**2+y**2))**2*y**3&
     &+32*y**3*a**2*exp(-b*(x**2+y**2))**2*cos(c*x*y)**2*b**2+2*x**5*a**&
     &2*exp(-b*(x**2+y**2))**2*sin(c*x*y)*c**3*cos(c*x*y)+8*a*cos(c*x*y)&
     &*b**2*y**3*exp(-b*(x**2+y**2))-16*a**2*cos(c*x*y)**2*b**3*y**5*exp&
     &(-b*(x**2+y**2))**2
      D223 = 0.0
      D233 = 0.0
      D311 = 0.0
      D312 = 0.0
      D313 = 0.0
      D322 = 0.0
      D323 = 0.0
      D333 = 0.0


