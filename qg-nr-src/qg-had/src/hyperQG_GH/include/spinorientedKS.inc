
      s1 = M/2+0.d0
      s3 = z**2+y**2+x**2-a**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-&
     &2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(thet&
     &a)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(&
     &theta)*z)+0.d0
      s5 = sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x*&
     &*2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*&
     &x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*co&
     &s(theta)*x*sin(theta)*z))+0.d0
      s8 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)*&
     &*2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s9 = s8+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a**2&
     &+x**4+z**4+0.d0
      s7 = s9+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-&
     &2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(thet&
     &a)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(&
     &theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+&
     &x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a&
     &**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*co&
     &s(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqrt(&
     &2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a&
     &**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**&
     &2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s6 = 1/s7+0.d0
      s4 = s5*s6+0.d0
      s2 = s3*s4+0.d0
      H = s1*s2+0.d0
      nult = 1+0.d0
      s1 = -4*cos(theta)**2*x*a**2+2*x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y*&
     &*2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*&
     &cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta&
     &)*x*sin(theta)*z)+2*cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*&
     &sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y&
     &**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(the&
     &ta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y+2*x**3+4*s&
     &in(theta)*cos(theta)*z*a**2+2*x*a**2+2*x*y**2+2*x*z**2+0.d0
      s2 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      nulx = s1*s2+0.d0
      nuly = (sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**&
     &2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z&
     &**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**&
     &2*cos(theta)*x*sin(theta)*z))*y-2*a*cos(theta)*x+2*a*sin(theta)*z)&
     &/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z&
     &**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z*&
     &*2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z&
     &)+a**2)+0.d0
      s1 = 2*z*x**2-2*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**&
     &2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z&
     &**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**&
     &2+8*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y+2*z**3+4*cos(t&
     &heta)**2*z*a**2-2*z*a**2+4*a**2*sin(theta)*x*cos(theta)+2*z*y**2+2&
     &*z*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**&
     &2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(&
     &theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s2 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      nulz = s1*s2+0.d0
      dtH = 0+0.d0
      dtnult = 0+0.d0
      dtnulx = 0+0.d0
      dtnuly = 0+0.d0
      dtnulz = 0+0.d0
      s2 = M/2+0.d0
      s4 = 2*x+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*&
     &y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**&
     &2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*x&
     &*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2+8*sin(theta)&
     &*cos(theta)*z*a**2)/2+0.d0
      s6 = sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x*&
     &*2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*&
     &x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*co&
     &s(theta)*x*sin(theta)*z))+0.d0
      s9 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)*&
     &*2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s10 = s9+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a**&
     &2+x**4+z**4+0.d0
      s8 = s10+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2&
     &-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(the&
     &ta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin&
     &(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*&
     &a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*c&
     &os(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqrt&
     &(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+&
     &a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)*&
     &*2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s7 = 1/s8+0.d0
      s5 = s6*s7+0.d0
      s3 = s4*s5+0.d0
      s1 = s2*s3+0.d0
      s4 = M*(z**2+y**2+x**2-a**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z*&
     &*2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(t&
     &heta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*s&
     &in(theta)*z))/4+0.d0
      s6 = 1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*&
     &x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**&
     &2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*&
     &cos(theta)*x*sin(theta)*z))+0.d0
      s10 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)&
     &**2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s11 = s10+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a*&
     &*2+x**4+z**4+0.d0
      s9 = s11+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2&
     &-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(the&
     &ta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin&
     &(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*&
     &a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*c&
     &os(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqrt&
     &(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+&
     &a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)*&
     &*2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s8 = 1/s9+0.d0
      s9 = 4*x+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*&
     &y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**&
     &2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*x&
     &*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2+8*sin(theta)&
     &*cos(theta)*z*a**2)+0.d0
      s7 = s8*s9+0.d0
      s5 = s6*s7+0.d0
      s3 = s4*s5+0.d0
      s5 = -M*(z**2+y**2+x**2-a**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z&
     &**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(&
     &theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*&
     &sin(theta)*z))/2+0.d0
      s7 = sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x*&
     &*2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*&
     &x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*co&
     &s(theta)*x*sin(theta)*z))+0.d0
      s11 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)&
     &**2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s12 = s11+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a*&
     &*2+x**4+z**4+0.d0
      s10 = s12+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**&
     &2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(th&
     &eta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*si&
     &n(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z&
     &**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z*&
     &*2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z&
     &)+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2&
     &*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*&
     &cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqr&
     &t(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4&
     &+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)&
     &**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s9 = 1/s10**2+0.d0
      s11 = 4*x*y**2+4*x*a**2-8*cos(theta)**2*x*a**2+8*sin(theta)*cos(th&
     &eta)*z*a**2+4*x**3+0.d0
      s12 = s11+4*x*z**2+y**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2&
     &+8*sin(theta)*cos(theta)*z*a**2)/2+0.d0
      s13 = s12-a**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*&
     &(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2+8*sin(th&
     &eta)*cos(theta)*z*a**2)/2+0.d0
      s10 = s13+2*x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2&
     &-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*&
     &a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+x&
     &**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a*&
     &*2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos&
     &(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*x*y**2+4&
     &*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2+8*sin(theta)*cos(th&
     &eta)*z*a**2)/2+z**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**&
     &2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2&
     &*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta&
     &)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2+8*s&
     &in(theta)*cos(theta)*z*a**2)/2+0.d0
      s8 = s9*s10+0.d0
      s6 = s7*s8+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dxH = s1+s2+0.d0
      dxnult = 0+0.d0
      s4 = -4*a**2*cos(theta)**2+2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z&
     &**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(&
     &theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*&
     &sin(theta)*z)+0.d0
      s5 = s4+0.d0
      s7 = x/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2&
     &*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*&
     &cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*x*y**&
     &2+4*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2+8*sin(theta)*cos&
     &(theta)*z*a**2)+0.d0
      s8 = cos(theta)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**&
     &2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z&
     &**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**&
     &2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y*(4*x+1/sqrt(2*y**2*x**2+2&
     &*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4&
     &+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8&
     &*a**2*cos(theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**&
     &2-8*cos(theta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2))+0.d0
      s6 = s7+s8+0.d0
      s3 = s5+s6+0.d0
      s2 = s3+6*x**2+2*a**2+2*y**2+2*z**2+0.d0
      s3 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      s1 = s2*s3+0.d0
      s4 = 4*cos(theta)**2*x*a**2-2*x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**&
     &2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*c&
     &os(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)&
     &*x*sin(theta)*z)-2*cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*s&
     &qrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y*&
     &*4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(thet&
     &a)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y-2*x**3-4*si&
     &n(theta)*cos(theta)*z*a**2-2*x*a**2-2*x*y**2-2*x*z**2+0.d0
      s5 = 1/((x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)+a**2)**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*&
     &x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**&
     &4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*&
     &a**2+8*a**2*cos(theta)*x*sin(theta)*z))*(2*x+1/sqrt(2*y**2*x**2+2*&
     &a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+&
     &2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*&
     &a**2*cos(theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2&
     &-8*cos(theta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2)/2)+0.d0
      s3 = s4*s5+0.d0
      s5 = 2*cos(theta)**2*x*a**2-x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*&
     &z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos&
     &(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x&
     &*sin(theta)*z)-cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(&
     &2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a&
     &**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**&
     &2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y-x**3-2*sin(thet&
     &a)*cos(theta)*z*a**2-x*a**2-x*y**2-x*z**2+0.d0
      s6 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))**3*(4*x+1/sqrt(2*y**2*x**2+2*a*&
     &*2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*&
     &z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a*&
     &*2*cos(theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8&
     &*cos(theta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2))+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dxnulx = s1+s2+0.d0
      s1 = (1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2&
     &*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z*&
     &*2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2&
     &*cos(theta)*x*sin(theta)*z))*y*(4*x+1/sqrt(2*y**2*x**2+2*a**2*x**2&
     &+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x*&
     &*2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(&
     &theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(th&
     &eta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2))/2-2*a*cos(theta))/&
     &(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+a**2)+0.d0
      s2 = -(sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*&
     &x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**&
     &2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*&
     &cos(theta)*x*sin(theta)*z))*y-2*a*cos(theta)*x+2*a*sin(theta)*z)/(&
     &x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+&
     &a**2)**2*(2*x+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(theta)**2*x*a**2+8*sin(t&
     &heta)*cos(theta)*z*a**2)/2)+0.d0
      dxnuly = s1+s2+0.d0
      s3 = 4*x*z-1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2&
     &*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4&
     &+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8&
     &*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y*(4*x+1/sqrt(2*y**&
     &2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x&
     &**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**&
     &2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3&
     &+4*x*z**2-8*cos(theta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2))+0.d0
      s2 = s3+4*sin(theta)*cos(theta)*a**2+z/sqrt(2*y**2*x**2+2*a**2*x**&
     &2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x&
     &**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos&
     &(theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8*cos(t&
     &heta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2)+0.d0
      s3 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      s1 = s2*s3+0.d0
      s4 = -2*z*x**2+2*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x*&
     &*2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+&
     &z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a*&
     &*2+8*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y-2*z**3-4*cos(&
     &theta)**2*z*a**2+2*z*a**2-4*a**2*sin(theta)*x*cos(theta)-2*z*y**2-&
     &2*z*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a*&
     &*2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos&
     &(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s5 = 1/((x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)+a**2)**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*&
     &x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**&
     &4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*&
     &a**2+8*a**2*cos(theta)*x*sin(theta)*z))*(2*x+1/sqrt(2*y**2*x**2+2*&
     &a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+&
     &2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*&
     &a**2*cos(theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2&
     &-8*cos(theta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2)/2)+0.d0
      s3 = s4*s5+0.d0
      s5 = -z*x**2+sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2&
     &*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4&
     &+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8&
     &*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y-z**3-2*cos(theta)&
     &**2*z*a**2+z*a**2-2*a**2*sin(theta)*x*cos(theta)-z*y**2-z*sqrt(2*y&
     &**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4&
     &+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x&
     &**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s6 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))**3*(4*x+1/sqrt(2*y**2*x**2+2*a*&
     &*2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*&
     &z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a*&
     &*2*cos(theta)*x*sin(theta)*z)*(4*x*y**2+4*x*a**2+4*x**3+4*x*z**2-8&
     &*cos(theta)**2*x*a**2+8*sin(theta)*cos(theta)*z*a**2))+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dxnulz = s1+s2+0.d0
      s2 = M/2+0.d0
      s4 = 2*y+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*&
     &y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**&
     &2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*y&
     &*x**2+4*y*z**2-4*y*a**2+4*y**3)/2+0.d0
      s6 = sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x*&
     &*2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*&
     &x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*co&
     &s(theta)*x*sin(theta)*z))+0.d0
      s9 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)*&
     &*2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s10 = s9+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a**&
     &2+x**4+z**4+0.d0
      s8 = s10+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2&
     &-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(the&
     &ta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin&
     &(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*&
     &a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*c&
     &os(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqrt&
     &(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+&
     &a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)*&
     &*2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s7 = 1/s8+0.d0
      s5 = s6*s7+0.d0
      s3 = s4*s5+0.d0
      s1 = s2*s3+0.d0
      s4 = M*(z**2+y**2+x**2-a**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z*&
     &*2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(t&
     &heta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*s&
     &in(theta)*z))/4+0.d0
      s6 = 1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*&
     &x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**&
     &2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*&
     &cos(theta)*x*sin(theta)*z))+0.d0
      s10 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)&
     &**2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s11 = s10+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a*&
     &*2+x**4+z**4+0.d0
      s9 = s11+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2&
     &-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(the&
     &ta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin&
     &(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*&
     &a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*c&
     &os(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqrt&
     &(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+&
     &a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)*&
     &*2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s8 = 1/s9+0.d0
      s9 = 4*y+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*&
     &y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**&
     &2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*y&
     &*x**2+4*y*z**2-4*y*a**2+4*y**3)+0.d0
      s7 = s8*s9+0.d0
      s5 = s6*s7+0.d0
      s3 = s4*s5+0.d0
      s5 = -M*(z**2+y**2+x**2-a**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z&
     &**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(&
     &theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*&
     &sin(theta)*z))/2+0.d0
      s7 = sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x*&
     &*2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*&
     &x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*co&
     &s(theta)*x*sin(theta)*z))+0.d0
      s11 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)&
     &**2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s12 = s11+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a*&
     &*2+x**4+z**4+0.d0
      s10 = s12+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**&
     &2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(th&
     &eta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*si&
     &n(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z&
     &**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z*&
     &*2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z&
     &)+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2&
     &*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*&
     &cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqr&
     &t(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4&
     &+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)&
     &**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s9 = 1/s10**2+0.d0
      s11 = 4*y*x**2+4*y*z**2-4*y*a**2+4*y**3+0.d0
      s12 = s11+2*y*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2&
     &-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*&
     &a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+y&
     &**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a*&
     &*2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos&
     &(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+4&
     &*y*z**2-4*y*a**2+4*y**3)/2+0.d0
      s10 = s12-a**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*&
     &(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3)/2+x**2/sqrt(2*y**2*x**2+2*a**2&
     &*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z*&
     &*2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2&
     &*cos(theta)*x*sin(theta)*z)*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3)/2+&
     &z**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a&
     &**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*co&
     &s(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+&
     &4*y*z**2-4*y*a**2+4*y**3)/2+0.d0
      s8 = s9*s10+0.d0
      s6 = s7*s8+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dyH = s1+s2+0.d0
      dynult = 0+0.d0
      s3 = x/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2&
     &*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*&
     &cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*y*x**&
     &2+4*y*z**2-4*y*a**2+4*y**3)+cos(theta)/sqrt(2*z**2+2*y**2+2*x**2-2&
     &*a**2+2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**&
     &2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4&
     &*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y*(4&
     &*y+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a&
     &**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*co&
     &s(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+&
     &4*y*z**2-4*y*a**2+4*y**3))+0.d0
      s2 = s3+2*cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**&
     &2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x&
     &**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**&
     &2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a+4*x*y+0.d0
      s3 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      s1 = s2*s3+0.d0
      s4 = 4*cos(theta)**2*x*a**2-2*x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**&
     &2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*c&
     &os(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)&
     &*x*sin(theta)*z)-2*cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*s&
     &qrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y*&
     &*4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(thet&
     &a)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y-2*x**3-4*si&
     &n(theta)*cos(theta)*z*a**2-2*x*a**2-2*x*y**2-2*x*z**2+0.d0
      s5 = 1/((x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)+a**2)**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*&
     &x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**&
     &4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*&
     &a**2+8*a**2*cos(theta)*x*sin(theta)*z))*(2*y+1/sqrt(2*y**2*x**2+2*&
     &a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+&
     &2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*&
     &a**2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3&
     &)/2)+0.d0
      s3 = s4*s5+0.d0
      s5 = 2*cos(theta)**2*x*a**2-x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*&
     &z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos&
     &(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x&
     &*sin(theta)*z)-cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(&
     &2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a&
     &**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**&
     &2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y-x**3-2*sin(thet&
     &a)*cos(theta)*z*a**2-x*a**2-x*y**2-x*z**2+0.d0
      s6 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))**3*(4*y+1/sqrt(2*y**2*x**2+2*a*&
     &*2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*&
     &z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a*&
     &*2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3))+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dynulx = s1+s2+0.d0
      s2 = 1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*&
     &x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**&
     &2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*&
     &cos(theta)*x*sin(theta)*z))*y*(4*y+1/sqrt(2*y**2*x**2+2*a**2*x**2+&
     &2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**&
     &2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(t&
     &heta)*x*sin(theta)*z)*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3))/2+sqrt(&
     &2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*&
     &z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos&
     &(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x&
     &*sin(theta)*z))+0.d0
      s3 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)+0.d0
      s1 = s2*s3+0.d0
      s2 = -(sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*&
     &x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**&
     &2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*&
     &cos(theta)*x*sin(theta)*z))*y-2*a*cos(theta)*x+2*a*sin(theta)*z)/(&
     &x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+&
     &a**2)**2*(2*y+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3)/2)+0.d0
      dynuly = s1+s2+0.d0
      s3 = -1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2&
     &*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z*&
     &*2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2&
     &*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y*(4*y+1/sqrt(2*y**2*x**&
     &2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z&
     &**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**&
     &2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+4*y*z**2-4*y*a**2+4*&
     &y**3))-2*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**&
     &2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z&
     &**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**&
     &2*cos(theta)*x*sin(theta)*z))*sin(theta)*a+0.d0
      s2 = s3+4*z*y+z/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3)+0.d0
      s3 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      s1 = s2*s3+0.d0
      s4 = -2*z*x**2+2*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x*&
     &*2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+&
     &z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a*&
     &*2+8*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y-2*z**3-4*cos(&
     &theta)**2*z*a**2+2*z*a**2-4*a**2*sin(theta)*x*cos(theta)-2*z*y**2-&
     &2*z*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a*&
     &*2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos&
     &(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s5 = 1/((x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)+a**2)**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*&
     &x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**&
     &4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*&
     &a**2+8*a**2*cos(theta)*x*sin(theta)*z))*(2*y+1/sqrt(2*y**2*x**2+2*&
     &a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+&
     &2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*&
     &a**2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3&
     &)/2)+0.d0
      s3 = s4*s5+0.d0
      s5 = -z*x**2+sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2&
     &*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4&
     &+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8&
     &*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y-z**3-2*cos(theta)&
     &**2*z*a**2+z*a**2-2*a**2*sin(theta)*x*cos(theta)-z*y**2-z*sqrt(2*y&
     &**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4&
     &+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x&
     &**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s6 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))**3*(4*y+1/sqrt(2*y**2*x**2+2*a*&
     &*2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*&
     &z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a*&
     &*2*cos(theta)*x*sin(theta)*z)*(4*y*x**2+4*y*z**2-4*y*a**2+4*y**3))+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dynulz = s1+s2+0.d0
      s2 = M/2+0.d0
      s4 = 2*z+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*&
     &y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**&
     &2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*z&
     &*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*sin(t&
     &heta)*x*cos(theta))/2+0.d0
      s6 = sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x*&
     &*2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*&
     &x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*co&
     &s(theta)*x*sin(theta)*z))+0.d0
      s9 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)*&
     &*2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s10 = s9+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a**&
     &2+x**4+z**4+0.d0
      s8 = s10+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2&
     &-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(the&
     &ta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin&
     &(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*&
     &a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*c&
     &os(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqrt&
     &(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+&
     &a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)*&
     &*2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s7 = 1/s8+0.d0
      s5 = s6*s7+0.d0
      s3 = s4*s5+0.d0
      s1 = s2*s3+0.d0
      s4 = M*(z**2+y**2+x**2-a**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z*&
     &*2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(t&
     &heta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*s&
     &in(theta)*z))/4+0.d0
      s6 = 1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*&
     &x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**&
     &2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*&
     &cos(theta)*x*sin(theta)*z))+0.d0
      s10 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)&
     &**2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s11 = s10+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a*&
     &*2+x**4+z**4+0.d0
      s9 = s11+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2&
     &-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(the&
     &ta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin&
     &(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*&
     &a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*c&
     &os(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqrt&
     &(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+&
     &a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)*&
     &*2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s8 = 1/s9+0.d0
      s9 = 4*z+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*&
     &y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**&
     &2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*z&
     &*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*sin(t&
     &heta)*x*cos(theta))+0.d0
      s7 = s8*s9+0.d0
      s5 = s6*s7+0.d0
      s3 = s4*s5+0.d0
      s5 = -M*(z**2+y**2+x**2-a**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z&
     &**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(&
     &theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*&
     &sin(theta)*z))/2+0.d0
      s7 = sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*x*&
     &*2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*&
     &x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*co&
     &s(theta)*x*sin(theta)*z))+0.d0
      s11 = 2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-4*cos(theta)&
     &**2*x**2*a**2-2*y**2*a**2+y**4+a**4+0.d0
      s12 = s11+8*a**2*cos(theta)*x*sin(theta)*z+4*cos(theta)**2*z**2*a*&
     &*2+x**4+z**4+0.d0
      s10 = s12+2*z**2*x**2+y**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**&
     &2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(th&
     &eta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*si&
     &n(theta)*z)-a**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z&
     &**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z*&
     &*2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z&
     &)+x**2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2&
     &*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*&
     &cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+z**2*sqr&
     &t(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4&
     &+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)&
     &**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s9 = 1/s10**2+0.d0
      s11 = 4*z*y**2-4*z*a**2+8*a**2*sin(theta)*x*cos(theta)+8*cos(theta&
     &)**2*z*a**2+4*z**3+0.d0
      s12 = s11+4*z*x**2+y**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2&
     &+8*a**2*sin(theta)*x*cos(theta))/2+0.d0
      s13 = s12-a**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*&
     &(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*s&
     &in(theta)*x*cos(theta))/2+0.d0
      s10 = s13+x**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*&
     &(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*s&
     &in(theta)*x*cos(theta))/2+2*z*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*&
     &z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos&
     &(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x&
     &*sin(theta)*z)+z**2/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**&
     &2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2&
     &*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta&
     &)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a&
     &**2*sin(theta)*x*cos(theta))/2+0.d0
      s8 = s9*s10+0.d0
      s6 = s7*s8+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dzH = s1+s2+0.d0
      dznult = 0+0.d0
      s4 = x/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2&
     &*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*&
     &cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*z*y**&
     &2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*sin(theta&
     &)*x*cos(theta))+0.d0
      s5 = cos(theta)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**&
     &2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z&
     &**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**&
     &2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y*(4*z+1/sqrt(2*y**2*x**2+2&
     &*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4&
     &+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8&
     &*a**2*cos(theta)*x*sin(theta)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**&
     &2+8*cos(theta)**2*z*a**2+8*a**2*sin(theta)*x*cos(theta)))+0.d0
      s3 = s4+s5+0.d0
      s2 = s3+4*sin(theta)*cos(theta)*a**2+4*x*z+0.d0
      s3 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      s1 = s2*s3+0.d0
      s4 = 4*cos(theta)**2*x*a**2-2*x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**&
     &2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*c&
     &os(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)&
     &*x*sin(theta)*z)-2*cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*s&
     &qrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y*&
     &*4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(thet&
     &a)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y-2*x**3-4*si&
     &n(theta)*cos(theta)*z*a**2-2*x*a**2-2*x*y**2-2*x*z**2+0.d0
      s5 = 1/((x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)+a**2)**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*&
     &x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**&
     &4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*&
     &a**2+8*a**2*cos(theta)*x*sin(theta)*z))*(2*z+1/sqrt(2*y**2*x**2+2*&
     &a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+&
     &2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*&
     &a**2*cos(theta)*x*sin(theta)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2&
     &+8*cos(theta)**2*z*a**2+8*a**2*sin(theta)*x*cos(theta))/2)+0.d0
      s3 = s4*s5+0.d0
      s5 = 2*cos(theta)**2*x*a**2-x*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*&
     &z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos&
     &(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x&
     &*sin(theta)*z)-cos(theta)*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(&
     &2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a&
     &**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**&
     &2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z))*a*y-x**3-2*sin(thet&
     &a)*cos(theta)*z*a**2-x*a**2-x*y**2-x*z**2+0.d0
      s6 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))**3*(4*z+1/sqrt(2*y**2*x**2+2*a*&
     &*2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*&
     &z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a*&
     &*2*cos(theta)*x*sin(theta)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8&
     &*cos(theta)**2*z*a**2+8*a**2*sin(theta)*x*cos(theta)))+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dznulx = s1+s2+0.d0
      s1 = (1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2&
     &*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z*&
     &*2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2&
     &*cos(theta)*x*sin(theta)*z))*y*(4*z+1/sqrt(2*y**2*x**2+2*a**2*x**2&
     &+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x*&
     &*2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(&
     &theta)*x*sin(theta)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(th&
     &eta)**2*z*a**2+8*a**2*sin(theta)*x*cos(theta)))/2+2*a*sin(theta))/&
     &(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &+a**2)+0.d0
      s2 = -(sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2*a**2*&
     &x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**&
     &2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*&
     &cos(theta)*x*sin(theta)*z))*y-2*a*cos(theta)*x+2*a*sin(theta)*z)/(&
     &x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**&
     &2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2&
     &*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+&
     &a**2)**2*(2*z+1/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z*&
     &*2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**&
     &2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)&
     &*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*&
     &sin(theta)*x*cos(theta))/2)+0.d0
      dznuly = s1+s2+0.d0
      s3 = 2*x**2-1/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y*(4*z+1/sqrt(2*y*&
     &*2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+&
     &x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x*&
     &*2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*z*y**2-4*z*a**2+4*z**&
     &3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*sin(theta)*x*cos(theta)))&
     &+6*z**2+4*a**2*cos(theta)**2+0.d0
      s2 = s3-2*a**2+2*y**2+2*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)+z/sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*&
     &y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**&
     &2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)*(4*z&
     &*y**2-4*z*a**2+4*z**3+4*z*x**2+8*cos(theta)**2*z*a**2+8*a**2*sin(t&
     &heta)*x*cos(theta))+0.d0
      s3 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))+0.d0
      s1 = s2*s3+0.d0
      s4 = -2*z*x**2+2*sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x*&
     &*2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+&
     &z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a*&
     &*2+8*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y-2*z**3-4*cos(&
     &theta)**2*z*a**2+2*z*a**2-4*a**2*sin(theta)*x*cos(theta)-2*z*y**2-&
     &2*z*sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a*&
     &*2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos&
     &(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s5 = 1/((x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2&
     &*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta&
     &)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(t&
     &heta)*z)+a**2)**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*&
     &x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**&
     &4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*&
     &a**2+8*a**2*cos(theta)*x*sin(theta)*z))*(2*z+1/sqrt(2*y**2*x**2+2*&
     &a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+&
     &2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*&
     &a**2*cos(theta)*x*sin(theta)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2&
     &+8*cos(theta)**2*z*a**2+8*a**2*sin(theta)*x*cos(theta))/2)+0.d0
      s3 = s4*s5+0.d0
      s5 = -z*x**2+sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+2&
     &*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4&
     &+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8&
     &*a**2*cos(theta)*x*sin(theta)*z))*sin(theta)*a*y-z**3-2*cos(theta)&
     &**2*z*a**2+z*a**2-2*a**2*sin(theta)*x*cos(theta)-z*y**2-z*sqrt(2*y&
     &**2*x**2+2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4&
     &+x**4+z**4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x&
     &**2*a**2+8*a**2*cos(theta)*x*sin(theta)*z)+0.d0
      s6 = 1/(x**2+z**2+y**2+sqrt(2*y**2*x**2+2*a**2*x**2+2*y**2*z**2-2*&
     &a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*z**2*x**2+4*cos(theta)&
     &**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a**2*cos(theta)*x*sin(th&
     &eta)*z)+a**2)/sqrt(2*z**2+2*y**2+2*x**2-2*a**2+2*sqrt(2*y**2*x**2+&
     &2*a**2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**&
     &4+2*z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+&
     &8*a**2*cos(theta)*x*sin(theta)*z))**3*(4*z+1/sqrt(2*y**2*x**2+2*a*&
     &*2*x**2+2*y**2*z**2-2*a**2*z**2-2*y**2*a**2+y**4+a**4+x**4+z**4+2*&
     &z**2*x**2+4*cos(theta)**2*z**2*a**2-4*cos(theta)**2*x**2*a**2+8*a*&
     &*2*cos(theta)*x*sin(theta)*z)*(4*z*y**2-4*z*a**2+4*z**3+4*z*x**2+8&
     &*cos(theta)**2*z*a**2+8*a**2*sin(theta)*x*cos(theta)))+0.d0
      s4 = s5*s6+0.d0
      s2 = s3+s4+0.d0
      dznulz = s1+s2
