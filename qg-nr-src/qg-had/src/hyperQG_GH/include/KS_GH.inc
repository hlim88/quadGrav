      g00 = -1+2*H*nult**2
      g01 = 2*H*nulx*nult
      g02 = 2*H*nuly*nult
      g03 = 2*H*nulz*nult
      g11 = 1+2*H*nulx**2
      g12 = 2*H*nulx*nuly
      g13 = 2*H*nulx*nulz
      g22 = 1+2*H*nuly**2
      g23 = 2*H*nuly*nulz
      g33 = 1+2*H*nulz**2
      D100 = 2*dxH*nult**2+4*H*nult*dxnult
      D101 = 2*dxH*nulx*nult+2*H*dxnulx*nult+2*H*nulx*dxnult
      D102 = 2*dxH*nuly*nult+2*H*dxnuly*nult+2*H*nuly*dxnult
      D103 = 2*dxH*nulz*nult+2*H*dxnulz*nult+2*H*nulz*dxnult
      D111 = 2*dxH*nulx**2+4*H*nulx*dxnulx
      D112 = 2*dxH*nulx*nuly+2*H*dxnulx*nuly+2*H*nulx*dxnuly
      D113 = 2*dxH*nulx*nulz+2*H*dxnulx*nulz+2*H*nulx*dxnulz
      D122 = 2*dxH*nuly**2+4*H*nuly*dxnuly
      D123 = 2*dxH*nuly*nulz+2*H*dxnuly*nulz+2*H*nuly*dxnulz
      D133 = 2*dxH*nulz**2+4*H*nulz*dxnulz
      D200 = 2*dyH*nult**2+4*H*nult*dynult
      D201 = 2*dyH*nulx*nult+2*H*dynulx*nult+2*H*nulx*dynult
      D202 = 2*dyH*nuly*nult+2*H*dynuly*nult+2*H*nuly*dynult
      D203 = 2*dyH*nulz*nult+2*H*dynulz*nult+2*H*nulz*dynult
      D211 = 2*dyH*nulx**2+4*H*nulx*dynulx
      D212 = 2*dyH*nulx*nuly+2*H*dynulx*nuly+2*H*nulx*dynuly
      D213 = 2*dyH*nulx*nulz+2*H*dynulx*nulz+2*H*nulx*dynulz
      D222 = 2*dyH*nuly**2+4*H*nuly*dynuly
      D223 = 2*dyH*nuly*nulz+2*H*dynuly*nulz+2*H*nuly*dynulz
      D233 = 2*dyH*nulz**2+4*H*nulz*dynulz
      D300 = 2*dzH*nult**2+4*H*nult*dznult
      D301 = 2*dzH*nulx*nult+2*H*dznulx*nult+2*H*nulx*dznult
      D302 = 2*dzH*nuly*nult+2*H*dznuly*nult+2*H*nuly*dznult
      D303 = 2*dzH*nulz*nult+2*H*dznulz*nult+2*H*nulz*dznult
      D311 = 2*dzH*nulx**2+4*H*nulx*dznulx
      D312 = 2*dzH*nulx*nuly+2*H*dznulx*nuly+2*H*nulx*dznuly
      D313 = 2*dzH*nulx*nulz+2*H*dznulx*nulz+2*H*nulx*dznulz
      D322 = 2*dzH*nuly**2+4*H*nuly*dznuly
      D323 = 2*dzH*nuly*nulz+2*H*dznuly*nulz+2*H*nuly*dznulz
      D333 = 2*dzH*nulz**2+4*H*nulz*dznulz
      
      K11 = 1/sqrt(1+2*H*nult**2)**3*(2*dxH*nult*nulx+4*H**2*nult**3*dxn&
     &ulx+2*H*dxnult*nulx+4*H**2*nulx**3*dxnult+8*H**3*nulx*nulz*nult**3&
     &*dxnulz+8*H**3*nulx*nuly*nult**3*dxnuly+2*H*nult*dxnulx+2*H*nult*n&
     &ulx**2*dzH*nulz+4*H**2*nult**3*nulx**2*dzH*nulz+4*H**2*nulx*nulz*n&
     &ult*dxnulz+4*H*nulx*nulz**2*dxH*nult-8*H**3*nulx*nulz**2*dxnult*nu&
     &lt**2+4*H**2*nult*nulx*dznulx*nulz+8*H**3*nult**3*nulx*dznulx*nulz&
     &+2*H*nult*nulx**2*dyH*nuly+4*H**2*nult*nulx*dynulx*nuly+4*H*nulx*n&
     &uly**2*dxH*nult-8*H**3*nulx*nuly**2*dxnult*nult**2+4*H**2*nulx*nul&
     &y*nult*dxnuly+6*H*nult*nulx**3*dxH+4*H**2*nult**3*nulx**3*dxH+8*H*&
     &*2*nult*nulx**2*dxnulx+16*H**3*nult**3*nulx**2*dxnulx-8*H**3*nulx*&
     &*3*dxnult*nult**2+4*H**2*nulx*nuly**2*dxnult+4*H**2*nulx*nulz**2*d&
     &xnult-4*H**2*dxnult*nulx*nult**2+4*H**2*nult**3*nulx**2*dyH*nuly+8&
     &*H**3*nult**3*nulx*dynulx*nuly)
      s1 = -1
      s4 = -dxH*nult*nuly-H*dxnult*nuly-2*H**2*nult**3*dynulx-H*dynult*n&
     &ulx-2*H**2*nult*nulz*nulx*dznuly-2*H**2*nulx**3*dynult-2*H**2*nuly&
     &**3*dxnult-8*H**3*nult**3*nuly*nulx*dynuly-4*H*nult*dxH*nulx**2*nu&
     &ly-H*nult*dynulx-H*nult*dxnuly-dyH*nult*nulx-2*H**2*nult*nulz*dznu&
     &lx*nuly-2*H**2*nuly*nulz*nult*dxnulz-4*H**3*nuly*nulz*nult**3*dxnu&
     &lz-4*H**2*nult**3*dxH*nulx**2*nuly-4*H**2*nult*nulx*dxnulx*nuly-8*&
     &H**3*nult**3*nulx*dxnulx*nuly-4*H*nult*dyH*nulx*nuly**2-4*H**2*nul&
     &t**3*dyH*nulx*nuly**2-4*H**2*nult*nuly*nulx*dynuly-2*H*nult*nulz*d&
     &zH*nulx*nuly+4*H**3*nulx*nuly**2*dynult*nult**2-4*H**3*nult**3*nul&
     &z*dznulx*nuly-2*H*nulx*nulz**2*dyH*nult+4*H**3*nulx*nulz**2*dynult&
     &*nult**2
      s3 = s4-2*H**2*nulx*nulz*nult*dynulz-2*H**2*nult**3*dxnuly-4*H**3*&
     &nult**3*nulz*nulx*dznuly+4*H**3*nulx**2*nuly*dxnult*nult**2-4*H**2&
     &*nult**3*nulz*dzH*nulx*nuly-4*H**3*nulx*nulz*nult**3*dynulz-2*H*nu&
     &ly*nulz**2*dxH*nult+4*H**3*nuly*nulz**2*dxnult*nult**2-2*H**2*nult&
     &*nulx**2*dxnuly-2*H**2*nult*dynulx*nuly**2-4*H**3*nult**3*dynulx*n&
     &uly**2-4*H**3*nult**3*nulx**2*dxnuly+2*H**2*dxnult*nuly*nult**2+2*&
     &H**2*dynult*nulx*nult**2-2*H**2*nulx**2*nuly*dxnult-2*H*nuly**3*dx&
     &H*nult+4*H**3*nuly**3*dxnult*nult**2-2*H**2*nuly**2*nult*dxnuly-4*&
     &H**3*nuly**2*nult**3*dxnuly-2*H**2*nuly*nulz**2*dxnult-2*H*nulx**3&
     &*dyH*nult+4*H**3*nulx**3*dynult*nult**2-2*H**2*nulx**2*nult*dynulx&
     &-4*H**3*nulx**2*nult**3*dynulx-2*H**2*nulx*nuly**2*dynult-2*H**2*n&
     &ulx*nulz**2*dynult
      s4 = 1/sqrt(1+2*H*nult**2)**3
      s2 = s3*s4
      K12 = s1*s2
      s1 = 1/sqrt(1+2*H*nult**2)**3
      s3 = 2*H**2*nult**3*dxnulz+H*dznult*nulx+2*H**2*nulx**3*dznult+2*H&
     &**2*nulx**2*nult*dznulx+H*nult*dznulx+H*nult*dxnulz+2*H**2*nult**3&
     &*dznulx+dxH*nult*nulz+H*dxnult*nulz+dzH*nult*nulx+2*H**2*nulz**3*d&
     &xnult+2*H**2*nult*nuly*nulx*dynulz-4*H**3*nulx**2*nulz*dxnult*nult&
     &**2+4*H**2*nult**3*dzH*nulx*nulz**2+2*H**2*nulx*nuly*nult*dznuly+4&
     &*H**3*nulx*nuly*nult**3*dznuly+4*H**3*nuly*nulz*nult**3*dxnuly+2*H&
     &**2*nuly*nulz*nult*dxnuly+4*H**3*nult**3*nuly*nulx*dynulz+4*H*nult&
     &*dzH*nulx*nulz**2-4*H**3*nulx**3*dznult*nult**2+2*H*nult*nuly*dyH*&
     &nulx*nulz+2*H**2*nuly**2*nulz*dxnult+2*H*nulx*nuly**2*dzH*nult+8*H&
     &**3*nult**3*nulx*dxnulx*nulz+4*H**2*nult**3*dxH*nulx**2*nulz
      s2 = s3+4*H*nult*dxH*nulx**2*nulz+4*H**2*nult*nulz*nulx*dznulz+4*H&
     &**2*nult**3*nuly*dyH*nulx*nulz+2*H**2*nult*nuly*dynulx*nulz+4*H**3&
     &*nult**3*nuly*dynulx*nulz+2*H*nuly**2*nulz*dxH*nult-4*H**3*nuly**2&
     &*nulz*dxnult*nult**2+4*H**2*nult*nulx*dxnulx*nulz-4*H**3*nulx*nuly&
     &**2*dznult*nult**2+8*H**3*nult**3*nulz*nulx*dznulz-2*H**2*dxnult*n&
     &ulz*nult**2-2*H**2*dznult*nulx*nult**2-4*H**3*nulx*nulz**2*dznult*&
     &nult**2+2*H**2*nult*nulx**2*dxnulz+4*H**3*nult**3*nulx**2*dxnulz+2&
     &*H**2*nult*dznulx*nulz**2+4*H**3*nult**3*dznulx*nulz**2+2*H**2*nul&
     &x**2*nulz*dxnult-4*H**3*nulz**3*dxnult*nult**2+2*H**2*nulz**2*nult&
     &*dxnulz+4*H**3*nulz**2*nult**3*dxnulz+2*H*nulx**3*dzH*nult+2*H*nul&
     &z**3*dxH*nult+2*H**2*nulx*nuly**2*dznult+2*H**2*nulx*nulz**2*dznul&
     &t+4*H**3*nulx**2*nult**3*dznulx
      K13 = s1*s2
      s1 = -1
      s2 = (-4*H**2*nuly**3*dynult-2*H*nult*dynuly-2*dyH*nult*nuly-4*H**&
     &2*nult**3*dynuly+8*H**3*nuly*nulz**2*dynult*nult**2-2*H*dynult*nul&
     &y-4*H**2*nult**3*nuly**2*dxH*nulx-8*H**3*nulx*nuly*nult**3*dxnuly-&
     &4*H**2*nuly*nulz*nult*dynulz-4*H*nult*nulx**2*dyH*nuly-4*H**2*nult&
     &*nulx*dynulx*nuly-4*H*nuly*nulz**2*dyH*nult-2*H*nulx*nuly**2*dxH*n&
     &ult-4*H**2*nulx*nuly*nult*dxnuly-8*H**3*nuly*nulz*nult**3*dynulz-4&
     &*H**2*nult*nuly*dznuly*nulz-8*H**2*nult*nuly**2*dynuly-4*H**2*nulx&
     &**2*nuly*dynult+8*H**3*nuly**3*dynult*nult**2-4*H**2*nuly*nulz**2*&
     &dynult-16*H**3*nult**3*nuly**2*dynuly-2*H*nult*nuly**2*dzH*nulz-4*&
     &H**2*nult**3*nuly**2*dzH*nulz+4*H**2*dynult*nuly*nult**2-8*H**3*nu&
     &lt**3*nulx*dynulx*nuly-6*H*nult*nuly**3*dyH-4*H**2*nult**3*nuly**3&
     &*dyH-8*H**3*nult**3*nuly*dznuly*nulz+8*H**3*nulx**2*nuly*dynult*nu&
     &lt**2)/sqrt(1+2*H*nult**2)**3
      K22 = s1*s2
      s1 = -1
      s4 = -H*dynult*nulz-dzH*nult*nuly-H*nult*dynulz-H*dznult*nuly-4*H*&
     &*2*nult**3*dzH*nuly*nulz**2-2*H**2*nulz**3*dynult-2*H**2*nuly**3*d&
     &znult-H*nult*dznuly+4*H**3*nuly**2*nulz*dynult*nult**2-dyH*nult*nu&
     &lz-4*H*nult*dzH*nuly*nulz**2-4*H**3*nult**3*nulx*nuly*dxnulz-4*H**&
     &2*nult*nulz*nuly*dznulz+4*H**3*nulx**2*nuly*dznult*nult**2-4*H*nul&
     &t*dyH*nuly**2*nulz+4*H**3*nulx**2*nulz*dynult*nult**2-2*H**2*nult*&
     &*3*dznuly-2*H**2*nulx*nuly*nult*dznulx-4*H**3*nulx*nuly*nult**3*dz&
     &nulx+4*H**3*nuly*nulz**2*dznult*nult**2-2*H**2*nult*nuly**2*dynulz&
     &-4*H**3*nult**3*nuly**2*dynulz-2*H**2*nult*dznuly*nulz**2-4*H**3*n&
     &ult**3*dznuly*nulz**2-2*H**2*nulx**2*nulz*dynult-2*H**2*nuly**2*nu&
     &lz*dynult
      s3 = -2*H*nulz**3*dyH*nult+4*H**3*nulz**3*dynult*nult**2-2*H**2*nu&
     &lz**2*nult*dynulz-4*H**3*nulz**2*nult**3*dynulz-2*H**2*nulx**2*nul&
     &y*dznult-2*H*nuly**3*dzH*nult+4*H**3*nuly**3*dznult*nult**2-2*H**2&
     &*nuly**2*nult*dznuly-4*H**3*nuly**2*nult**3*dznuly-2*H**2*nult**3*&
     &dynulz+s4-2*H*nulx**2*nuly*dzH*nult-4*H**3*nult**3*nulx*dxnuly*nul&
     &z-2*H**2*nult*nulx*nuly*dxnulz-2*H*nulx**2*nulz*dyH*nult-4*H**2*nu&
     &lt**3*nulx*dxH*nuly*nulz-8*H**3*nult**3*nulz*nuly*dznulz-2*H**2*nu&
     &lx*nulz*nult*dynulx-4*H**3*nulx*nulz*nult**3*dynulx-2*H**2*nult*nu&
     &lx*dxnuly*nulz-2*H**2*nuly*nulz**2*dznult+2*H**2*dynult*nulz*nult*&
     &*2+2*H**2*dznult*nuly*nult**2-2*H*nult*nulx*dxH*nuly*nulz-4*H**2*n&
     &ult**3*dyH*nuly**2*nulz-4*H**2*nult*nuly*dynuly*nulz-8*H**3*nult**&
     &3*nuly*dynuly*nulz
      s4 = 1/sqrt(1+2*H*nult**2)**3
      s2 = s3*s4
      K23 = s1*s2
      K33 = 1/sqrt(1+2*H*nult**2)**3*(2*H*dznult*nulz+4*H**2*nult**3*dzn&
     &ulz+2*H*nult*dznulz+2*dzH*nult*nulz+4*H**2*nulz**3*dznult+8*H**3*n&
     &ulx*nulz*nult**3*dxnulz+4*H**2*nult**3*nulz**2*dyH*nuly-8*H**3*nul&
     &x**2*nulz*dznult*nult**2+4*H*nult*nulx**2*dzH*nulz+4*H**2*nulx*nul&
     &z*nult*dxnulz+4*H**2*nult**3*nulz**3*dzH+6*H*nult*nulz**3*dzH+8*H*&
     &*2*nult*nulz**2*dznulz+16*H**3*nult**3*nulz**2*dznulz+4*H**2*nulx*&
     &*2*nulz*dznult+4*H**2*nuly**2*nulz*dznult-8*H**3*nulz**3*dznult*nu&
     &lt**2+2*H*nulx*nulz**2*dxH*nult+4*H**2*nult*nulx*dznulx*nulz+8*H**&
     &3*nult**3*nulx*dznulx*nulz+4*H**2*nuly*nulz*nult*dynulz+2*H*nuly*n&
     &ulz**2*dyH*nult+8*H**3*nuly*nulz*nult**3*dynulz+4*H**2*nult**3*nul&
     &z**2*dxH*nulx-8*H**3*nuly**2*nulz*dznult*nult**2-4*H**2*dznult*nul&
     &z*nult**2+4*H**2*nult*nuly*dznuly*nulz+4*H*nult*nuly**2*dzH*nulz+8&
     &*H**3*nult**3*nuly*dznuly*nulz)
