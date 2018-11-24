C these routines are 1x2 refinements with a characteristic matrix 
c------------------------------------------------------------------------------
C routines in this section are
C  d3pro122c(uc,uf,d1c,d2c,d3c,charc,charf)       2x1  2 O 
C  d3pro124c(uc,uf,d1c,d2c,d3c,charc,charf)       2x1  4 O 
C  d3pro12nc(uc,uf,d1c,d2c,d3c,charc,charf,order) 2x1  p O 
c------------------------------------------------------------------------------
C 
      integer function d3pro122c( uc,uf,nxc,nyc,nzc,charmc,charmf)
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8  uc(nxc,nyc,nzc),uf(2*nxc-1,2*nyc-1,2*nzc-1)
      real*8  charmc(nxc,nyc,nzc),charmf(2*nxc-1,2*nyc-1,2*nzc-1)
      real*8  badval,bigval

      parameter ( badval = 9.9e020, bigval = 1.0e010  )

      integer pil,l

      integer   ic , jc , kc
      integer   iff , jf, kf
      integer ierr,d2lintxyzc
      integer d3pro122c_ 

      nxf = 2*nxc-1
      nyf = 2*nyc-1
      nzf = 2*nzc-1

      d3pro122c = 0

      d3pro122c =  d3pro122c_(uc,uf,nxc,nyc,nzc,charmc,charmf)

      return
      end

C second order (linear interpolation on a 2x1 mesh) with char functions
c-----------------------------------------------------------------------
      integer function  d3pro122c_(uc,uf,d1c,d2c,d3c,chc,chf)

         implicit       none

         logical        ccdim

         integer        d1c,     d2c,     d3c
         real*8         uc(d1c,d2c,d3c),
     *                  uf(2*d1c-1,2*d2c-1,2*d3c-1)
         real*8         chc(d1c,d2c,d3c),
     *                  chf(2*d1c-1,2*d2c-1,2*d3c-1)

         integer        i1,      i2,      i3,
     *                  d1f,     d2f,     d3f

      integer maxsize
      parameter ( maxsize = 2049 )
      real*8 utempf(maxsize),utempc(maxsize),cthc(maxsize)
      integer kc,jc,ic,kf,jf,iff
      integer d1122


      d3pro122c_ = 0

      d1f = 2 * d1c - 1
      d2f = 2 * d2c - 1
      d3f = 2 * d3c - 1
C at first I do all the lines that I can ......

      do kc = 1 , d3c
         kf = 2*kc - 1
         do jc = 1 , d2c
            jf = 2*jc - 1
            call d3vputv(uc,utempc,d1c,d2c,d3c,jc,kc,1,1)
            call d3vputv(chc,cthc,d1c,d2c,d3c,jc,kc,1,1)
            d3pro122c_ =  d1122(d1c,utempc,utempf,cthc)
            if (d3pro122c_.eq.1) return
            call dvput3v(utempf,uf,d1f,d2f,d3f,jf,kf,1)
         end do
      end do
C now fill all the  missing lines

      do kc = 1 , d3c
         kf =  2*kc - 1
         do iff = 1 , d1f
            call d3vputv(uf,utempc,d1f,d2f,d3f,iff,kf,2,2)
            call d3vputv(chf,cthc,d1f,d2f,d3f,iff,kf,2,2)
            d3pro122c_ =  d1122(d1c,utempc,utempf,cthc)
            call dvput3v(utempf,uf,d1f,d2f,d3f,iff,kf,2)
         end do
      end do
C now fill in the missing planes

      do jf = 1 , d2f
         do iff = 1 , d1f
            call d3vputv(uf,utempc,d1f,d2f,d3f,iff,jf,2,3)
            call d3vputv(chf,cthc,d1f,d2f,d3f,iff,jf,2,3)
            d3pro122c_ =  d1122(d1c,utempc,utempf,cthc)
            call dvput3v(utempf,uf,d1f,d2f,d3f,iff,jf,3)
         end do
      end do

      return
      end
c-----------------------------------------------------------------------
      integer function d1122(nc,uc,uf,ch)
      implicit none
      integer nc
      real*8 uc(nc),uf(2*nc-1)
      real*8 ch(nc)
      integer nf
      real*8         cc(2),      ec(2)
      data
     * cc / 0.50000e0, 0.50000e0 /,
     * ec / 1.50000e0,-0.50000e0/ 
      integer maxnh,nh
      parameter ( maxnh = 100)
      integer start(maxnh),stop(maxnh)
      real*8 val,oldval
      integer order
      parameter (order = 2)

      integer jc,jf,i,l

      d1122 = 0
      nf = 2*nc - 1

C first thing is to inject where we can inject .....

      do jc = 1 , nc
         uf(2*jc-1) = uc(jc)*ch(jc)
      end do

C next thing is to find the start and stop
      nh = 1
      start(1) = nh
      oldval = ch(1)
      stop(nh) = 0

      do i = 1 , nc
         if (ch(i).ne.oldval) then
            if (ch(i).ne.1) then
               stop(nh) = 2*(i - 1)-1
            else
               nh = nh + 1
               start(nh) = 2*i-1
               stop(nh) = 0
            end if
         end if
         oldval = ch(i)
      end do

      if (stop(nh).eq.0) stop(nh) = 2*nc-1

      do l = 1 , nh
         if ((stop(l)-start(l)).lt.order-1) then
            d1122 = 1
            return
         end if

         do jf = start(l)+1,stop(l)-1,2
            jc = int(0.5*(jf-1+1))
            uf(jf) =  uf(jf-1)*cc(1)
     .             +  uf(jf+1)*cc(2)
         end do
      end do

      if (nh.gt.1) then

         do l = 1 , nh
            if(start(l).gt.1) then
               jf = start(l) - 1
               uf(jf) = ec(1)*uf(jf+1)+ec(2)*uf(jf+3)
            end if
            if (stop(l).gt.1.and.stop(l).lt.nf) then
               jf = stop(l) + 1
               uf(jf) = ec(1)*uf(jf-1)+ec(2)*uf(jf-3)
            end if

         end do
      end if


      return
      end
C==============================================================================
      integer function d3pro124c( uc,uf,nxc,nyc,nzc,charmc,charmf)
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8  uc(nxc,nyc,nzc),uf(2*nxc-1,2*nyc-1,2*nzc-1)
      real*8  charmc(nxc,nyc,nzc),charmf(2*nxc-1,2*nyc-1,2*nzc-1)
      real*8  badval,bigval

      parameter ( badval = 9.9e007, bigval = 1.0e006  )

      integer pil,l

      integer   ic , jc , kc
      integer   iff , jf, kf
      integer ierr,d4lintxyzc
      integer icharsten3
      integer didit,d3pro124c_

      nxf = 2*nxc-1
      nyf = 2*nyc-1
      nzf = 2*nzc-1

      d3pro124c = 0
      didit = 0

      d3pro124c =  d3pro124c_(uc,uf,nxc,nyc,nzc,charmc,charmf)

      return
      end
c-----------------------------------------------------------------------
c
c     2:1 cubic interpolation ...  with a characteristic function
c
C  from a 2x1 mesh
c-----------------------------------------------------------------------

      integer function d3pro124c_(uc,uf,d1c,d2c,d3c,chc,chf)

         implicit       none

         logical        ccdim

         integer        d1c,     d2c,     d3c
         real*8         uc(d1c,d2c,d3c),
     *                  uf(2*d1c-1,2*d2c-1,2*d3c-1)
         real*8         chc(d1c,d2c,d3c),
     *                  chf(2*d1c-1,2*d2c-1,2*d3c-1)

         integer        i1,      i2,      i3,
     *                  d1f,     d2f,     d3f

      integer maxsize
      parameter ( maxsize = 2049 )
      real*8 utempf(maxsize),utempc(maxsize),cthc(maxsize)
      integer kc,jc,ic,kf,jf,iff
      integer d1124


      d1f = 2 * d1c - 1
      d2f = 2 * d2c - 1
      d3f = 2 * d3c - 1
C at first I do all the lines that I can ......

      do kc = 1 , d3c
         kf = 2*kc - 1
         do jc = 1 , d2c
            jf = 2*jc - 1
            call d3vputv(uc,utempc,d1c,d2c,d3c,jc,kc,1,1)
            call d3vputv(chc,cthc,d1c,d2c,d3c,jc,kc,1,1)
            d3pro124c_ =  d1124(d1c,utempc,utempf,cthc)
            if (d3pro124c_.eq.1) return
            call dvput3v(utempf,uf,d1f,d2f,d3f,jf,kf,1)
         end do
      end do
C now fill all the  missing lines

      do kc = 1 , d3c
         kf =  2*kc - 1
         do iff = 1 , d1f
            call d3vputv(uf,utempc,d1f,d2f,d3f,iff,kf,2,2)
            call d3vputv(chf,cthc,d1f,d2f,d3f,iff,kf,2,2)
            d3pro124c_ =  d1124(d1c,utempc,utempf,cthc)
            call dvput3v(utempf,uf,d1f,d2f,d3f,iff,kf,2)
         end do
      end do
C now fill in the missing planes

      do jf = 1 , d2f
         do iff = 1 , d1f
            call d3vputv(uf,utempc,d1f,d2f,d3f,iff,jf,2,3)
            call d3vputv(chf,cthc,d1f,d2f,d3f,iff,jf,2,3)
            d3pro124c_ =  d1124(d1c,utempc,utempf,cthc)
            call dvput3v(utempf,uf,d1f,d2f,d3f,iff,jf,3)
         end do
      end do

      return
      end
c-----------------------------------------------------------------------
      integer function  d1124(nc,uc,uf,ch)
      implicit none
      integer nc
      real*8 uc(nc),uf(2*nc-1)
      real*8 ch(nc)
      integer nf
      real*8         lc(4),      cc(4),      rc(4), ec(4)
      data
     * lc / 0.31250e0, 0.93750e0,-0.31250e0, 0.06250e0 /,
     * cc /-0.06250e0, 0.56250e0, 0.56250e0,-0.06250e0 /,
     * rc / 0.06250e0,-0.31250e0, 0.93750e0, 0.31250e0 /,
     * ec / 4.00000e0,-6.00000e0, 4.00000e0,-1.00000e0/
      integer maxnh,nh
      parameter ( maxnh = 100)
      integer start(maxnh),stop(maxnh)
      real*8 val,oldval
      integer order
      parameter ( order = 4)

      integer jc,jf,i,l

      nf = 2*nc - 1
      d1124 = 0

C first thing is to inject where we can inject .....

      do jc = 1 , nc
         uf(2*jc-1) = uc(jc)*ch(jc)
      end do

C next thing is to find the start and stop
      nh = 1
      start(1) = nh
      oldval = ch(1)
      stop(nh) = 0

      do i = 1 , nc
         if (ch(i).ne.oldval) then
            if (ch(i).ne.1) then
               stop(nh) = 2*(i - 1)-1
            else
               nh = nh + 1
               start(nh) = 2*i-1
               stop(nh) = 0
            end if
         end if
         oldval = ch(i)
      end do

      if (stop(nh).eq.0) stop(nh) = 2*nc-1

      do l = 1 , nh
         if ((stop(l)-start(l)).lt.order-1) then
            d1124 = 1
            return
         end if

         jf = start(l) + 1
         uf(jf) = uf(jf-1)*lc(1)+uf(jf+1)*lc(2)+uf(jf+3)*lc(3)+
     .            uf(jf+5)*lc(4)
   
         do jf = start(l)+3,stop(l)-3,2
            jc = int(0.5*(jf-3+1))
            uf(jf) =  uf(jf-3)*cc(1)*ch(jc)+uf(jf-1)*cc(2)*ch(jc+1)
     .             +  uf(jf+1)*cc(3)*ch(jc+2)+uf(jf+3)*cc(4)*ch(jc+3)
         end do
   
         jf = stop(l) - 1
         jc = int(0.5*(jf-3+1))
   
         uf(jf) = uf(jf+1)*rc(4)+uf(jf-1)*rc(3)+
     .            uf(jf-3)*rc(2)+uf(jf-5)*rc(1)

      end do

      if (nh.gt.1) then

         do l = 1 , nh
            if(start(l).gt.1) then
               jf = start(l) - 1
               uf(jf) = ec(1)*uf(jf+1)+ec(2)*uf(jf+2)+
     .                  ec(3)*uf(jf+3)+ec(4)*uf(jf+4)
            end if
            if (stop(l).gt.1.and.stop(l).lt.nf) then
               jf = stop(l) + 1
               uf(jf) = ec(1)*uf(jf-1)+ec(2)*uf(jf-2)+
     .                  ec(3)*uf(jf-3)+ec(4)*uf(jf-4)
            end if

         end do
      end if


      return
      end


c
c-----------------------------------------------------------------------
      logical function ccdim(n)

         implicit      none

         integer       n

         ccdim =  mod(n,2) .eq. 1  .and.  n .ge. 5

         return

      end



C----------------------------------------------------------------------------
C
C     Vector pronlogation.
C
C----------------------------------------------------------------------------

      SUBROUTINE dvprln_sc(V1,V2,INC,N,chc)

         IMPLICIT     LOGICAL (A-Z)

         INTEGER      N,      INC
         real*8       V1(N),  V2(2*N-1)
         real*8       chc(n)

         INTEGER      J

         DO 10 J = 0 , N - 1
            V2(INC * J + 1) = V1(J+1)*chc(j+1)
 10      CONTINUE

         RETURN

      END
C----------------------------------------------------------------------------

      integer function d3pro12nc( uc,uf,nxc,nyc,nzc,charmc,charmf,order)
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8  uc(nxc,nyc,nzc),uf(2*nxc-1,2*nyc-1,2*nzc-1)
      real*8  charmc(nxc,nyc,nzc),charmf(2*nxc-1,2*nyc-1,2*nzc-1)
      real*8  badval,bigval
      integer order

      parameter ( badval = 9.9e007, bigval = 1.0e006  )

      integer pil,l

      integer   ic , jc , kc
      integer   iff , jf, kf
      integer ierr,dnlintxyzc
      integer icharsten3,d3pro12nc_

      nxf = 2*nxc-1
      nyf = 2*nyc-1
      nzf = 2*nzc-1

      d3pro12nc = 0

      d3pro12nc =  d3pro12nc_(uc,uf,nxc,nyc,nzc,order,charmc,charmf)

      return
      end
c-----------------------------------------------------------------------
      integer function d3pro12nc_(uc,uf,nxc,nyc,nzc,order,chc,chf) 
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8    uc(nxc,nyc,nzc),uf(2*nxc-1,2*nyc-1,2*nzc-1) 
      real*8    chc(nxc,nyc,nzc),chf(2*nxc-1,2*nyc-1,2*nzc-1)

      integer   ic , jc , kc
      integer   iff , jf, kf

      integer order

      integer d1nint2c,ierr
      integer wksp_size
      parameter ( wksp_size = 8000 )
      real*8 wksp(wksp_size)
      integer maxsize,d1nint2ec
      parameter ( maxsize = 2049 )
      real*8 utempf(maxsize),utempc(maxsize),cthc(maxsize)
      real*8 badv
      parameter ( badv = 1.0d20 )

      integer i,j,k

C the first thing that I must do is set up the coefficients ....
        d3pro12nc_ = 0

      nxf = 2*nxc-1
      nyf = 2*nyc-1
      nzf = 2*nzc-1

      d3pro12nc_  = d1nint2c(order,wksp,wksp_size)

      if (d3pro12nc_.ne.0) then
            return
      end if

C at first I do all the lines that I can ......

      do kc = 1 , nzc
         kf = 2*kc - 1
         do jc = 1 , nyc 
            jf = 2*jc - 1
            call d3vputv(uc,utempc,nxc,nyc,nzc,jc,kc,1,1)
            call d3vputv(chc,cthc,nxc,nyc,nzc,jc,kc,1,1)
                d3pro12nc_ = d1nint2ec(nxc,utempc,utempf,
     .                       wksp,wksp_size,order,cthc)
                if (d3pro12nc_.eq.1) return
            call dvput3v(utempf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill all the  missing lines

      do kc = 1 , nzc
         kf =  2*kc - 1 
         do iff = 1 , nxf
            call d3vputv(uf,utempc,nxf,nyf,nzf,iff,kf,2,2)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,kf,2,2)
            d3pro12nc_ = 
     .         d1nint2ec(nxc,utempc,utempf,wksp,wksp_size,order,cthc)
                if (d3pro12nc_.eq.1) return
            call dvput3v(utempf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

C now fill in the missing planes

      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,utempc,nxf,nyf,nzf,iff,jf,2,3)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,jf,2,3)
            d3pro12nc_ = 
     .         d1nint2ec(nxc,utempc,utempf,wksp,wksp_size,order,cthc)
                if (d3pro12nc_.eq.1) return
            call dvput3v(utempf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do

      return
      end
C 1x2 refinements using nth order interpolation
C evaluates the stencils via lud decomposition 
C The size of the work array( mat ) = order * (order+1) *(order -1)
C this is for the lud matrix + the pivit matrix stored for
C the order -1 possibilities
C
C
C
      integer function d1nint2c(order,mat,mat_size) 
      implicit none

      integer maxorder,order
      parameter ( maxorder = 20 )
      integer mat_size
      real*8 mat(mat_size)
      integer need_size
      integer construct_matrixc

      d1nint2c = 0

C      need_size = order*(order+1)*(order-1)
      need_size = order*(order+1)*(order+1)

      if (mat_size.lt.(need_size)) then
         write(6,*) 'Need to make a larger work array in d1nint2c'
         d1nint2c = 1
         return
      end if

      d1nint2c =  construct_matrixc(order,mat)

      return
      end

C======================================================================
      integer function d1nint2(order,mat,mat_size)
      implicit none

      integer maxorder,order
      parameter ( maxorder = 20 )
      integer mat_size
      real*8 mat(mat_size)
      integer need_size
      integer construct_matrix

      d1nint2 = 0

      need_size = order*(order+1)*(order-1)

      if (mat_size.lt.(need_size)) then
         write(6,*) 'Need to make a larger work array in d1nint2'
         d1nint2 = 1
         return
      end if

      d1nint2 =  construct_matrix(order,mat)

      return
      end

C======================================================================
      integer function d1nint2ec(n,uc,uf,mat,mat_size,order,chc)
      implicit none
      integer n,mat_size
      real*8 uc(n),uf(2*n-1),chc(n)
      real*8 mat(mat_size)
      integer order
        integer d1nint2ec_

      integer maxorder
      parameter ( maxorder = 20 )

      d1nint2ec =  d1nint2ec_(n,uc,uf,mat,order,chc)

      return
      end
C======================================================================
      integer function d1nint2ec_(n,uc,uf,mat,order,ch)
      implicit none
      integer n,order
      real*8 uc(n),uf(2*n-1),ch(n)
      real*8 mat(order,order+1)
C                coef , stencil

      integer maxorder
      parameter ( maxorder = 20 )
      real*8 input_vec(maxorder),coeff_vec(maxorder)

      integer maxnh,nh,l
      parameter ( maxnh = 100)
      integer start(maxnh),stop(maxnh)
      real*8 oldval,badval
      parameter ( badval = 1.0d20)


      integer ic,iff
      integer i,j,ip,iv
      integer length
      integer k,nf
      integer istop,istart
      integer ioff

      length = nint(order*0.5)
      ioff = mod(order,2)

      nf = 2*n - 1 
        d1nint2ec_ = 0

      do iff = 1 , nf
         uf(iff) = 0.0e0
      end do

      do ic = 1 , n
         iff = 2*ic - 1
         uf(iff) = ch(ic)*uc(ic)
      end do


C next thing is to find the start and stop
      nh = 1
      start(1) = nh
      oldval = ch(1)
      stop(nh) = 0

      do i = 1 , n
         if (ch(i).ne.oldval) then
            if (ch(i).ne.1) then
               stop(nh) = 2*(i - 1)-1
            else
               nh = nh + 1
               start(nh) = 2*i-1
               stop(nh) = 0
            end if
         end if
         oldval = ch(i)
      end do
      if (stop(nh).eq.0) stop(nh) = 2*n-1


      do l = 1 , nh

C first do the interior points ......... k = 1 

         k = 1

         istart = 2*(length-ioff) + start(l) - 1
         istop  = stop(l)  - length*2 + 1

         do i = 1 , length

            do iff =istart,istop,2 
               iv = iff+2*(i-1)+1
               ic = int((iv+1)*0.5)
               uf(iff) = uf(iff) + mat(i,k)*uf(iv) *ch(ic)
            end do
         end do

         j = i

         ip = 0
         do i = j , order
            ip = ip + 1

            do iff =istart,istop,2 
               iv = iff-2*(ip-1)-1
               ic = int((iv+1)*0.5)
               uf(iff) = uf(iff) + mat(i,k)*uf(iv)*ch(ic) 
            end do
         end do

C now lets get the end points ............

C first do all the positive terms that we must modify
         istop = 1
         istart = length 

            if ((stop(l) - (length-1)*2 + 1.gt.2*n-1).or.
     .        stop(l)-1.lt.1) then
                d1nint2ec_ = 1
                return
            end if

         do iff = stop(l) - 1 , stop(l) - (length-1)*2 + 1  , -2 
            uf(iff) = 0.0e0
            k = k + 1
   
            i = 1
   
            do i = 1 , length
               ip = iff + 2*i - 1
               ic = int((ip+1)*0.5)
   
               if (ip.le.stop(l)) then
                  uf(iff) = uf(iff) + mat(i,k)*uf(ip)*ch(ic)
               end if
            end do
      
            j = i - 1 
            do i = 1 , length -ioff
               j = j + 1
               ip = iff - 2*i + 1
               ic = int((ip+1)*0.5)
               if (ip.gt.start(l)) then
                  uf(iff) = uf(iff) + mat(j,k)*uf(ip) *ch(ic)
               end if
            end do
      
            istop = istop + 1
            ip = 0
            ic = iff - 2 *(length-ioff) - 1 + 2 
                if (ic-2*(istart-istop+1).lt.1) then
                    d1nint2ec_ = 1
                    return
                end if
            do i = istart, istop , -1 
               ic = ic - 2 
               iv = int((ic+1)*0.5)
               uf(iff) = uf(iff) + mat(i,k)*uf(ic) *ch(iv)
            end do
         end do
   
C now let get the other end points
         istop = 1
         do iff = start(l)+1, start(l)+1+(length-ioff-2)*2,2 
            uf(iff) = 0.0e0
            k = k + 1
   
            do i = 1 , length
               ip = iff + 2*i - 1
               if (ip.le.stop(l)) then
                  ic = int((ip+1)*0.5)
                  uf(iff) = uf(iff) + mat(i,k)*uf(ip)*ch(ic)
               end if
            end do
   
            j = i - 1 
            do i = 1 , length
               j = j + 1
               ip = iff - 2*i + 1
               if (ip.ge.start(l)) then
                  ic  = int((ip+1)*0.5)
                  uf(iff) = uf(iff) + mat(j,k)*uf(ip)*ch(ic)
               end if
            end do

            istop = istop + 1
            j = iff  + 2*(length-1) + 1 
   
            do i = order , order -  (length-istop)+ioff , -1
               j = j + 2 
               ic = int((j+1)*0.5)
               uf(iff) = uf(iff) + mat(i,k)*uf(j) *ch(ic)
            end do
         end do

         if (start(l).ne.1) then
            k = k + 1
            iff = start(l) - 1
            uf(iff) = 0.0e0

            do i = 1 , order
               uf(iff) = uf(iff) + mat(i,k)*uf(iff+i)
            end do
         end if

         if (stop(l).ne.nf) then
            k = k + 1
            iff = stop(l) + 1
            uf(iff) = 0.0e0
            do i = 1 , order
               uf(iff) = uf(iff) + mat(i,k)*uf(iff-i)
            end do
         end if

      end do


      return
      end
C======================================================================
      integer function construct_matrix(order,mat)
      implicit none
      integer order
      real*8 mat(order,order-1)
      integer maxorder
      parameter ( maxorder = 20 )
      real*8 mat2d(maxorder,maxorder),mat2ds(maxorder,maxorder)
      real*8 work(maxorder)
      integer ipvt(maxorder)
      integer construct_matrix_

      construct_matrix= construct_matrix_(order,mat2d,mat2ds,
     .                                    ipvt,mat,work)

      return
      end
C======================================================================
      integer function construct_matrix_(order,mat2d,mat2ds,
     .        ipvt,mat,work)
      implicit none

      integer order
C mat is stored as the coefficients and for each stencil
C coef# , stencil#
C 
      real*8 mat(order,order-1)
      real*8 mat2d(order,order),mat2ds(order,order)
      real*8 factorial
      real*8 work(order),det
      integer fa
      integer ipvt(order)
      integer length
      integer ifail
      integer ioff

      integer i,j,k,l
      integer iendp1

      construct_matrix_ = 0


      length = nint(order*0.5)
      k = 1
      ioff = mod(order,2)

C do the positive direction first ...
      fa = -1
      do i = 1 , length 
         fa = fa + 2
         do j = 1 , order
            mat2d(i,j) = fa**(j-1)/factorial(j-1) 
         end do
      end do

      fa = 1
      iendp1 = i
C nod do the negative direction ....
      do i = i , order
         fa = fa - 2   
         do j = 1 , order
            mat2d(i,j) = fa**(j-1)/factorial(j-1) 
         end do
      end do

      call copy(order*order,mat2d,mat2ds)

C now do a lud decomposition on this and then an iverse and then
C store the coefficients

      call sgetrf(mat2d,order,order,ipvt,ifail)
      call sgetri(mat2d,order,order,ipvt,det,work,01)

C now get out the coefficients .....

      construct_matrix_ = ifail
      if (ifail.ne.0) return

C now copy into mat array
      call copy_mat2d_mat(order,mat2d,k,mat)

C now lets get the other stencils starting from the second point
C and then going forward ......
C assume 2x1 refinements ...

C first modify all the positive  terms that we must modify ...

      do l = length  - 1  , 1 , -1

C l is the number of terms to modify
         fa = -order + 1 + ioff
         k = k + 1
         call copy(order*order,mat2ds,mat2d)
C i is which term in the matrix ....
         do i = length , length - l + 1  , -1
            fa = fa - 2
            do j = 1 , order
               mat2d(i,j) = fa**(j-1)/factorial(j-1)
            end do
         end do

C now do a lud decomposition on this ........

         call sgetrf(mat2d,order,order,ipvt,ifail)
         call sgetri(mat2d,order,order,ipvt,det,work,01)

         construct_matrix_ = ifail
         if (ifail.ne.0) return

C now copy into mat array
         call copy_mat2d_mat(order,mat2d,k,mat)

      end do
C now modify the negative terms .....

      do l = length  - 1  , 1 + ioff , -1

C l is the number of terms to modify
C         fa = order - 1
         fa = -1 + length*2
         k = k + 1
         call copy(order*order,mat2ds,mat2d)
C i is which term in the matrix ....
         do i = order , order - l + 1 +ioff  , -1
            fa = fa + 2
            do j = 1 , order
               mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
            end do
         end do

C now do a lud decomposition on this ........

         call sgetrf(mat2d,order,order,ipvt,ifail)
         call sgetri(mat2d,order,order,ipvt,det,work,01)

         construct_matrix_ = ifail
         if (ifail.ne.0) return

C now copy into mat array
         call copy_mat2d_mat(order,mat2d,k,mat)

      end do

      return
      end
C======================================================================
      subroutine copy_mat2d_mat(order,mat2d,k,mat)
      implicit none

      integer order
      real*8 mat(order,order-1)
      real*8 mat2d(order,order)
      integer k

      integer i,j

      do i = 1 , order
         mat(i,k) = mat2d(1,i)
      end do

      return
      end
C======================================================================
      subroutine copy_mat2d_matc(order,mat2d,k,mat)
      implicit none

      integer order
      real*8 mat(order,order+1)
      real*8 mat2d(order,order)
      integer k

      integer i,j

      do i = 1 , order
         mat(i,k) = mat2d(1,i)
      end do

      return
      end
C======================================================================
      real*8 function factorial(ival)
      implicit none
      integer ival

      integer i

      factorial = 1.0e0

      do i = 1 , ival
         factorial = factorial * i
      end do

      return
      end
C======================================================================
      subroutine copy(n,ain,aout)
      implicit none
      integer n
      real*8 ain(n),aout(n)

      integer i

      do i = 1 , n
         aout(i) = ain(i)
      end do

      return
      end
C======================================================================
c-----------------------------------------------------------------------
C routines in this section are
C       d3pro122(uc,uf,d1c,d2c,d3c)          2x1       second order
C       d3pro124(uc,uf,d1c,d2c,d3c)          2x1       fourth order
C       d3pro12n(uc,uf,d1c,d2c,d2c,order)    2x1       nth    order
c-----------------------------------------------------------------------

C second order (linear interpolation on a 2x1 mesh)
c-----------------------------------------------------------------------
      integer function d3pro122( uc,uf,nxc, nyc, nzc) 
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8    uc(nxc,nyc,nzc),uf(2*nxc-1,2*nyc-1,2*nzc-1)

      integer   ic , jc , kc
      integer   if , jf, kf

      nxf = 2*nxc-1
      nyf = 2*nyc-1
      nzf = 2*nzc-1

      d3pro122 = 0

c first injection

      do kc = 1, nzc 

         kf = 2*kc -1

         do jc = 1, nyc

            jf = 2*jc -1


            do ic = 1, nxc 

               if = 2*ic -1
               uf(if,jf,kf) = uc(ic,jc,kc)

            end do

            do if = 2 , nxf - 1 , 2

               uf(if,jf,kf) = 0.5e0*( uf(if-1,jf,kf) + uf(if+1,jf,kf) )

            end do

            if (jc.gt.1) then

c fill in missing lines

               do if = 1 , nxf

                  uf(if,jf-1,kf) = 0.5e0*( uf(if,jf-2,kf) +
     .                                     uf(if,jf,kf) )
               end do

            end if


         end do

         if (kc.gt.1) then

c fill in missing planes

            do jf = 1 , nyf

               do if = 1 , nxf

                  uf(if,jf,kf-1) = 0.5e0*( uf(if,jf,kf-2) +
     .                                     uf(if,jf,kf) )

               end do
            end do

         end if

      end do


      return
      end
c-----------------------------------------------------------------------
c
c     2:1 cubic interpolation ...
c
C  from a 2x1 mesh
c-----------------------------------------------------------------------

      integer function d3pro124(uc,uf,d1c,d2c,d3c)

         implicit       none

         logical        ccdim

         integer        d1c,     d2c,     d3c
         real*8         uc(d1c,d2c,d3c),
     *                  uf(2*d1c-1,2*d2c-1,2*d3c-1)

         integer        i1,      i2,      i3,
     *                  d1f,     d2f,     d3f

         real*8         lc(4),      cc(4),      rc(4)
         data
     *     lc / 0.31250e0,  0.93750e0, -0.31250e0,  0.06250e0 /,
     *     cc /-0.06250e0,  0.56250e0,  0.56250e0, -0.06250e0 /,
     *     rc / 0.06250e0, -0.31250e0,  0.93750e0,  0.31250e0 /


         d3pro124 = 0

         if( ccdim(d1c) .and. ccdim(d2c) .and. ccdim(d3c) ) then
            d1f = 2 * d1c - 1
            d2f = 2 * d2c - 1
            d3f = 2 * d3c - 1
            do 20 i3 = 1 , d3c
               do 10 i2 = 1  , d2c
                  call dv2i4_s(uc(1,i2,i3),uf(1,2*i2-1,2*i3-1),d1c)
 10            continue
 20         continue
            do 40 i3 = 1 , d3f , 2
               call dvsma4_s(uf(1,1,i3),uf(1,3,i3),
     *                     uf(1,5,i3),uf(1,7,i3),
     *                     uf(1,2,i3),lc,d1f)
                 do 30 i2 = 4 , d2f - 3 , 2
                  call dvsma4_s(uf(1,i2-3,i3),uf(1,i2-1,i3),
     *                        uf(1,i2+1,i3),uf(1,i2+3,i3),
     *                        uf(1,i2,i3),cc,d1f)
 30            continue
               call dvsma4_s(uf(1,d2f-6,i3),uf(1,d2f-4,i3),
     *                     uf(1,d2f-2,i3),uf(1,d2f,i3),
     *                     uf(1,d2f-1,i3),rc,d1f)
 40         continue
            do 50 i2 = 1 , d2f
               call dvsma4_s(uf(1,i2,1),uf(1,i2,3),
     *                     uf(1,i2,5),uf(1,i2,7),
     *                     uf(1,i2,2),lc,d1f)
 50         continue
            do 70 i3 = 4 , d3f - 3 , 2
               do 60 i2 = 1 , d2f
                  call dvsma4_s(uf(1,i2,i3-3),uf(1,i2,i3-1),
     *                        uf(1,i2,i3+1),uf(1,i2,i3+3),
     *                        uf(1,i2,i3),cc,d1f)
 60            continue
 70         continue
            do 80 i2 = 1 , d2f
               call dvsma4_s(uf(1,i2,d3f-6),uf(1,i2,d3f-4),
     *                     uf(1,i2,d3f-2),uf(1,i2,d3f),
     *                     uf(1,i2,d3f-1),rc,d1f)
 80         continue

         else
            write(*,*) '>>> d3pro124: Invalid dimension(s): ',
     *                 d1c, d2c, d3c
            d3pro124 = 1
         end if

         return

      end
c----------------------------------------------------------------------
c
c     V5 :=   Sum  [ S(i) * Vi ] ... used by multidimensional 2:1 cubic
c           i=1..4                   interpolation ...
c
c-----------------------------------------------------------------------

      subroutine dvsma4_s(v1,v2,v3,v4,v5,s,n)

         implicit       none

         integer        n
         real*8         v1(n),       v2(n),      v3(n),      v4(n),
     *                  v5(n)
         real*8         s(4)

         integer        j

         do 10 j = 1 , n
            v5(j) = s(1) * v1(j) + s(2) * v2(j) +
     *              s(3) * v3(j) + s(4) * v4(j)
 10      continue

         return

      end

c-----------------------------------------------------------------------
c
c     2:1 cubic interpolation of VC to VF ...
c
c-----------------------------------------------------------------------

      subroutine dv2i4_s(vc,vf,nc)

         implicit      none

         integer       nc
         real*8        vc(nc),     vf(*)

         integer       i,          nf

         real*8        lc(4),      cc(4),      rc(4)
         data
     *     lc / 0.31250e0,  0.93750e0, -0.31250e0,  0.06250e0 /,
     *     cc /-0.06250e0,  0.56250e0,  0.56250e0, -0.06250e0 /,
     *     rc / 0.06250e0, -0.31250e0,  0.93750e0,  0.31250e0 /

         if( nc .ge. 4 ) then
            nf = 2 * nc - 1
            call dvprln_s(vc,vf,2,nc)
            vf(2) = lc(1) * vf(1) + lc(2) * vf(3) +
     *              lc(3) * vf(5) + lc(4) * vf(7)
            do 10 i = 4 , nf - 3 , 2
               vf(i) = cc(1) * vf(i-3) + cc(2) * vf(i-1) +
     *                 cc(3) * vf(i+1) + cc(4) * vf(i+3)
 10         continue
            vf(nf-1) = rc(1) * vf(nf-6) + rc(2) * vf(nf-4) +
     *                 rc(3) * vf(nf-2) + rc(4) * vf(nf)
         else
            write(*,*) '>>> dv2i4_s:: Too few points for cubic '//
     *                 'interpolation ...'

         end if

         return

      end

C----------------------------------------------------------------------------
C
C     Vector pronlogation.
C
C----------------------------------------------------------------------------

      SUBROUTINE dvprln_s(V1,V2,INC,N)

         IMPLICIT     LOGICAL (A-Z)

         INTEGER      N,      INC
         REAL*8       V1(N),  V2(2*N-1)

         INTEGER      J

         DO 10 J = 0 , N - 1
            V2(INC * J + 1) = V1(J+1)
 10      CONTINUE

         RETURN

      END
C this will do nth order interpolation on a 2x1 mesh .....
c-----------------------------------------------------------------------
      integer function d3pro12n(uc,uf,nxc,nyc,nzc,order) 
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8    uc(nxc,nyc,nzc),uf(2*nxc-1,2*nyc-1,2*nzc-1) 

      integer   ic , jc , kc
      integer   iff , jf, kf

      integer order

      integer d1nint2,ierr
      integer wksp_size
      parameter ( wksp_size = 8000 )
      real*8 wksp(wksp_size)
      integer maxsize
      parameter ( maxsize = 2049 )
      real*8 utempf(maxsize),utempc(maxsize)


      d3pro12n = 0

C the first thing that I must do is set up the coefficients ....
      nxf = 2*nxc-1
      nyf = 2*nyc-1
      nzf = 2*nzc-1

      ierr = d1nint2(order,wksp,wksp_size)

      if (ierr.ne.0) then
         write(6,*) 'Error in d1nint2'
         d3pro12n = 1
      end if

C at first I do all the lines that I can ......

      do kc = 1 , nzc
         kf = 2*kc - 1
         do jc = 1 , nyc 
            jf = 2*jc - 1
            call d3vputv(uc,utempc,nxc,nyc,nzc,jc,kc,1,1)
            call d1nint2e(nxc,utempc,utempf,wksp,wksp_size,order)
            call dvput3v(utempf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill all the  missing lines

      do kc = 1 , nzc
         kf =  2*kc - 1 
         do iff = 1 , nxf
            call d3vputv(uf,utempc,nxf,nyf,nzf,iff,kf,2,2)
            call d1nint2e(nxc,utempc,utempf,wksp,wksp_size,order)
            call dvput3v(utempf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

C now fill in the missing planes

      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,utempc,nxf,nyf,nzf,iff,jf,2,3)
            call d1nint2e(nxc,utempc,utempf,wksp,wksp_size,order)
            call dvput3v(utempf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do

      return
      end
C 1x2 refinements using nth order interpolation
C evaluates the stencils via lud decomposition 
C The size of the work array( mat ) = order * (order+1) *(order -1)
C this is for the lud matrix + the pivit matrix stored for
C the order -1 possibilities
C
C
C
C======================================================================
      subroutine d1nint2e(n,uc,uf,mat,mat_size,order)
      implicit none
      integer n,mat_size
      real*8 uc(n),uf(2*n-1)
      real*8 mat(mat_size)
      integer order

      integer maxorder
      parameter ( maxorder = 20 )

      call d1nint2e_(n,uc,uf,mat,order)

      return
      end
C======================================================================
      subroutine d1nint2e_(n,uc,uf,mat,order)
      implicit none
      integer n,order
      real*8 uc(n),uf(2*n-1)
      real*8 mat(order,order-1)
C                coef , stencil

      integer maxorder
      parameter ( maxorder = 20 )
      real*8 input_vec(maxorder),coeff_vec(maxorder)

      integer ic,iff
      integer i,j,ip
      integer length
      integer k,nf
      integer istop,istart
      integer ioff

      length = nint(order*0.5)
      ioff = mod(order,2)

      nf = 2*n - 1 

      do ic = 1, nf
         uf(ic) = 0.0e0
      end do

      do ic = 1 , n
         iff = 2*ic - 1
         uf(iff) = uc(ic)
      end do

C first do the interior points ......... k = 1 

      k = 1

      istart = 2*(length-ioff)
      istop  = nf - length*2 + 1

      do i = 1 , length

         do iff =istart,istop,2 
            uf(iff) = uf(iff) + mat(i,k)*uf(iff+2*(i-1)+1) 
         end do
      end do

      j = i

      ip = 0
      do i = j , order
         ip = ip + 1

         do iff =istart,istop,2 
            uf(iff) = uf(iff) + mat(i,k)*uf(iff-2*(ip-1)-1)
         end do

      end do


C now lets get the end points ............

C first do all the positive terms that we must modify
      istop = 1
      istart = length 

      do iff = nf - 1 , nf - (length-1)*2 + 1  , -2 
         uf(iff) = 0.0e0
         k = k + 1

         i = 1

         do i = 1 , length
            ip = iff + 2*i - 1
            if (ip.le.nf) then
               uf(iff) = uf(iff) + mat(i,k)*uf(ip) 
            end if
         end do
   
         j = i - 1 
         do i = 1 , length -ioff
            j = j + 1
            ip = iff - 2*i + 1
            if (ip.gt.1) then
               uf(iff) = uf(iff) + mat(j,k)*uf(ip) 
            end if
         end do
   
         istop = istop + 1
         ip = 0
         ic = iff - 2 *(length-ioff) - 1 + 2 
         do i = istart, istop , -1 
            ic = ic - 2 
            uf(iff) = uf(iff) + mat(i,k)*uf(ic) 
         end do
      end do

C now let get the other end points
      istop = 1
      do iff = 2, (length-ioff-2)*2+2,2 
         uf(iff) = 0.0e0
         k = k + 1

C         i = 1
C         uf(iff)  = uf(iff) +  mat(i,k)*uf(iff+1)
         do i = 1 , length
            ip = iff + 2*i - 1
            if (ip.le.nf) then
               uf(iff) = uf(iff) + mat(i,k)*uf(ip)
            end if
         end do

         j = i - 1 
         do i = 1 , length
            j = j + 1
            ip = iff - 2*i + 1
            if (ip.ge.1) then
               uf(iff) = uf(iff) + mat(j,k)*uf(ip)
            end if
         end do

         istop = istop + 1
         j = iff  + 2*(length-1) + 1 

         do i = order , order -  (length-istop)+ioff , -1
            j = j + 2 
            uf(iff) = uf(iff) + mat(i,k)*uf(j) 
         end do

      end do


      return
      end
C======================================================================
C these routines are 1xn refinements with a characteristic matrix
c------------------------------------------------------------------------------
C routines in this section are
C  d3pro1n2c(uc,uf,d1c,d2c,d3c,d1f,d2f,d3f,
C            charc,charf,xf,yf,zf)                          nx1  2 O
C  d3pro1n4c(uc,uf,d1c,d2c,d3c,d1f,d2f,d3f,
C            charc,charf,xf,yf,zf)                          nx1  4 O
C  d3pro1nnc(uc,uf,d1c,d2c,d3c,d1f,d2f,d3f,
C            charc,charf,xf,yf,zf,order)                    nx1  p O
c------------------------------------------------------------------------------
      integer function d3pro1n2c( uc,uf,nxc,nyc,nzc,nxf,nyf,nzf,
     .                      charmc,charmf)
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8  uc(nxc,nyc,nzc),uf(nxf,nyf,nzf) 
      real*8  charmc(nxc,nyc,nzc),charmf(nxf,nyf,nzf) 
      real*8  badval,bigval


      parameter ( badval = 9.9e020, bigval = 1.0e010  )

      integer pil,l
      integer d3pro1n2c_

      integer   ic , jc , kc
      integer   iff , jf, kf
      integer ierr,d2lintxyzc
      integer icharsten3

      d3pro1n2c = 0

      d3pro1n2c= d3pro1n2c_(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,charmc,charmf)

      return
      end
C-------------------------------------------------------------------------
C prolongs from a 1 to n grid using second order techinques
C i.e. cubic interpolation
C
      integer function d3pro1n2c_(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,chc,chf)
      implicit none
      integer  nxf,nyf,nzf
      integer  nxc,nyc,nzc

      real*8 uf(nxf,nyf,nzf), uc(nxc,nyc,nzc)
      real*8 chf(nxf,nyf,nzf), chc(nxc,nyc,nzc)

      integer i,j,k
      integer iff,jf,kf,kff
      integer ic,jc,kc
      integer nxfact,nyfact,nzfact
      integer nxskip,nyskip,nzskip
      integer iskip,jskip,kskip 
      integer maxsize
      parameter ( maxsize = 2049)

      real*8 vtc(maxsize),vtf(maxsize),cthc(maxsize)
      integer dvi1q1c

C figure out the refinement factors in each dimension 

      d3pro1n2c_ = 0

      nxfact = (nxf-1)/(nxc-1)
      nyfact = (nyf-1)/(nyc-1)
      nzfact = (nzf-1)/(nzc-1)

      nxskip = nxfact - 1
      nyskip = nyfact - 1
      nzskip = nzfact - 1
C

      do kc = 1 , nzc
         kf = nzfact*kc-nzskip
         do jc = 1 , nyc
            jf = nyfact*jc-nyskip
            call d3vputv(uc,vtc,nxc,nyc,nzc,jc,kc,1,1)
            call d3vputv(chc,cthc,nxc,nyc,nzc,jc,kc,1,1)
            d3pro1n2c_ =  dvi1q1c(vtc,vtf,nxfact,nxc,cthc)
            if (d3pro1n2c_.eq.1) return
            call dvput3v(vtf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill in the missing lines

      do kc = 1 , nzc
         kf =  nzfact*kc-nzskip
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,kf,nyfact,2)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,kf,nyfact,2)
            d3pro1n2c_ =  dvi1q1c(vtc,vtf,nyfact,nyc,cthc) 
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

C now fill in the missing planes
      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,jf,nzfact,3)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,jf,nzfact,3)
            d3pro1n2c_ =  dvi1q1c(vtc,vtf,nzfact,nzc,cthc) 
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do
            
      return
      end
C-----------------------------------------------------------------------------C
      integer function  DVI1Q1c(V1,V2,IRHO,n1,charm)
C
      IMPLICIT    none
C
      integer     IRHO, N1,RHOM1
      REAL*8       V1(N1),     V2(*), charm(n1)
C
      REAL*8       C0, C1, C2, C3, SIGMA
      INTEGER      I, iskip
      integer iff,nxf
      integer irhom1
      integer max,jf,i1,i2
      parameter ( max = 32 )
      real*8 f1(max),f2(max),oldval
      real*8 f1e(max),f2e(max)
      real*8 f11e(max),f22e(max)
      integer maxnh
      parameter ( maxnh = 100 )
      integer nh,start(maxnh),stop(maxnh),l
      integer order
      parameter ( order = 2)

      dvi1q1c = 0
      irhom1 = irho - 1 
      nxf = irho * n1 -irhom1

C first thing is to inject where we can inject
      do i = 1 , n1
         iff = irho*i - irhom1 
         v2(iff) = v1(i)*charm(i)
      end do

      do iskip = 1 , irhom1
         f2(iskip) = 1.0e0*(-iskip+irho)/irho
         f1(iskip) =  1.0e0*iskip/irho
         f2e(iskip) = 1.0e0*(iskip -irho)/irho
         f1e(iskip) = 1.0e0*(2*irho-iskip)/irho
         f22e(iskip)= 1.0e0*(-iskip)/irho
         f11e(iskip)= 1.0e0*(irho+iskip)/irho
      end do

C now find the start and stop
      nh = 1
      start(1) = nh
      oldval = charm(1)
      stop(nh) = 0

      do i = 1 , n1
         if (charm(i).ne.oldval) then
            if (charm(i).ne.1) then
               stop(nh) = i - 1
            else
               nh = nh + 1
               start(nh) = i
               stop(nh) = 0
            end if
         end if
         oldval = charm(i)
      end do

      if (stop(nh).eq.0) stop(nh) = n1
C now the start and the stop are equal to their coarse grid values ......


      do l = 1 , nh
         if ((stop(l)-start(l)).lt.order-1) then
            dvi1q1c = 1
            return
         end if

         do iskip = 1 , irhom1
            do i = start(l) , stop(l) -1
               iff = irho*(i-1) + iskip+1
               v2(iff) = f2(iskip)*v1(i)*charm(i) + 
     .                   f1(iskip)*v1(i+1)*charm(i+1)
            end do
C now do the extrapolation points

            if (nh.gt.1) then
               if (start(l).gt.1) then
                  i = irho*(start(l)-1) - irhom1 + iskip
                  i1= i + (irho-iskip)
                  i2= i1 + irho 
                  v2(i) = f1e(iskip)*v2(i1) + f2e(iskip)*v2(i2)
               end if

               if (stop(l).gt.1.and.stop(l).lt.nxf) then
                  i = irho*(stop(l)) - irhom1 + iskip
                  i1= i - iskip
                  i2= i1 - irho
                  v2(i) = f22e(iskip)*v2(i2)+f11e(iskip)*v2(i1)
               end if


            end if
         end do
      end do

      return
      end
C-----------------------------------------------------------------------------C
      integer function d3pro1n4c( uc,uf,nxc,nyc,nzc,nxf,nyf,nzf,
     .                      charmc,charmf)
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8  uc(nxc,nyc,nzc),uf(nxf,nyf,nzf) 
      real*8  charmc(nxc,nyc,nzc),charmf(nxf,nyf,nzf) 
      real*8  badval,bigval

      parameter ( badval = 9.9e020, bigval = 1.0e010  )

      integer pil,l

      integer   ic , jc , kc
      integer   iff , jf, kf
      integer ierr,d4lintxyzc
      integer d3pro1n4c_ 

      d3pro1n4c = 0

      d3pro1n4c= d3pro1n4c_(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,charmc,charmf)

      return
      end
C-------------------------------------------------------------------------
      integer function d3pro1n4c_(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,chc,chf)
      implicit none
      integer  nxf,nyf,nzf
      integer  nxc,nyc,nzc

      real*8 uf(nxf,nyf,nzf), uc(nxc,nyc,nzc)
      real*8 chf(nxf,nyf,nzf), chc(nxc,nyc,nzc)

      integer i,j,k
      integer iff,jf,kf,kff
      integer ic,jc,kc
      integer nxfact,nyfact,nzfact
      integer nxskip,nyskip,nzskip
      integer iskip,jskip,kskip 

      integer max
      parameter ( max = 32 )
      real*8 ci(4,max),cbl(4,max),cbr(4,max)
      real*8 cel(4,max),cer(4,max)
      real*8 f1,f2

      integer maxsize
      parameter ( maxsize = 2049 )
      real*8 vtc(maxsize),vtf(maxsize),cthc(maxsize)
      integer dvcuq1_sc

C figure out the refinement factors in each dimension 

      d3pro1n4c_ = 0

      nxfact = (nxf-1)/(nxc-1)
      nyfact = (nyf-1)/(nyc-1)
      nzfact = (nzf-1)/(nzc-1)

      nxskip = nxfact - 1
      nyskip = nyfact - 1
      nzskip = nzfact - 1

      call get3co4(nxfact,ci,cbl,cbr,cel,cer)
C

      do kc = 1 , nzc
         kf = nzfact*kc-nzskip
         do jc = 1 , nyc
            jf = nyfact*jc-nyskip
            call d3vputv(uc,vtc,nxc,nyc,nzc,jc,kc,1,1)
            call d3vputv(chc,cthc,nxc,nyc,nzc,jc,kc,1,1)
            d3pro1n4c_ =  dvcuq1_sc(vtc,vtf,nxfact,nxc,
     .                              cthc,ci,cbl,cbr,cel,cer)
            if (d3pro1n4c_.eq.1) return
            call dvput3v(vtf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill in the missing lines
      call get3co4(nyfact,ci,cbl,cbr,cel,cer)

      do kc = 1 , nzc
         kf =  nzfact*kc-nzskip
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,kf,nyfact,2)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,kf,nyfact,2)
            d3pro1n4c_ =  dvcuq1_sc(vtc,vtf,nyfact,nyc,
     .                              cthc,ci,cbl,cbr,cel,cer)
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

      call get3co4(nzfact,ci,cbl,cbr,cel,cer)
C now fill in the missing planes
      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,jf,nzfact,3)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,jf,nzfact,3)
            d3pro1n4c_ =  dvcuq1_sc(vtc,vtf,nzfact,nzc,
     .                              cthc,ci,cbl,cbr,cel,cer)
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do
            
      return
      end
C---------------------------------------------------------------------
      integer function dvcuq1_sc(vc,uf,nfact,nc,ch,ci,cbl,cbr,cel,cer)
      implicit none
      integer nc,nf,nfact
      real*8 vc(nc),uf(nfact*nc-nfact+1),ch(nc)
      integer nfm1
      integer i,j,l,iff,jf,ic,jc
      integer max,iskip
      parameter ( max = 32 )
      integer maxnh
      parameter ( maxnh = 40 )
      integer nh,start(maxnh),stop(maxnh)

      real*8 ci(4,max),cbl(4,max),cbr(4,max)
      real*8 cel(4,max),cer(4,max)
      real*8 oldval
      integer order
      parameter ( order = 4 )

      dvcuq1_sc = 0

      nfm1 = nfact - 1 
      nf = nfact*nc - nfm1

C first thing is to inject where we can inject
      do i = 1 , nc
         iff = nfact*i - nfm1
         uf(iff) = vc(i)*ch(i)
      end do

C now find the start and stop
      nh = 1
      start(1) = nh
      oldval = ch(1)
      stop(nh) = 0

      do i = 1 , nc
         if (ch(i).ne.oldval) then
            if (ch(i).ne.1) then
               stop(nh) = i - 1
            else
               nh = nh + 1
               start(nh) = i
               stop(nh) = 0
            end if
         end if
         oldval = ch(i)
      end do

      if (stop(nh).eq.0) stop(nh) = nc
C now the start and the stop are equal to their coarse grid values ......
      do l = 1 , nh
         if ((stop(l)-start(l)).lt.order-1) then
            dvcuq1_sc = 1
            return
         end if
         do iskip = 1 , nfm1
            i = start(l)
               iff = nfact*(i-1) + iskip+1
               uf(iff)=cbl(1,iskip)*uf(iff-iskip)+
     .                 cbl(2,iskip)*uf(iff-iskip+nfact)+
     .                 cbl(3,iskip)*uf(iff-iskip+2*nfact)+
     .                 cbl(4,iskip)*uf(iff-iskip+3*nfact)

            do i = start(l)+1,stop(l)-2
               iff = nfact*(i-1) + iskip+1
               uf(iff)=ci(1,iskip)*uf(iff-iskip-nfact)+
     .                 ci(2,iskip)*uf(iff-iskip)+
     .                 ci(3,iskip)*uf(iff-iskip+nfact)+
     .                 ci(4,iskip)*uf(iff-iskip+2*nfact)
            end do

            i = stop(l) - 1 
               iff = nfact*(i-1) + iskip+1
               uf(iff) =cbr(1,iskip)*uf(iff-iskip-2*nfact)+
     .                  cbr(2,iskip)*uf(iff-iskip-nfact)+
     .                  cbr(3,iskip)*uf(iff-iskip)+
     .                  cbr(4,iskip)*uf(iff-iskip+nfact)


C now do the extrapolation points

            if (nh.gt.1) then
               if (start(l).gt.1) then
                  iff = nfact*(start(l)-1) - nfm1 + iskip
                  uf(iff)=cel(4,iskip)*uf(iff-iskip+nfact)+
     .                    cel(3,iskip)*uf(iff-iskip+2*nfact)+
     .                    cel(2,iskip)*uf(iff-iskip+3*nfact)+
     .                    cel(1,iskip)*uf(iff-iskip+4*nfact)
               end if

               if (stop(l).gt.1.and.stop(l).lt.nc) then
                  iff = nfact*(stop(l)) - nfm1 + iskip
                  uf(iff)=cer(4,iskip)*uf(iff-iskip)+
     .                    cer(3,iskip)*uf(iff-iskip-nfact)+
     .                    cer(2,iskip)*uf(iff-iskip-2*nfact)+
     .                    cer(1,iskip)*uf(iff-iskip-3*nfact)
               end if
            end if
         end do
      end do

      return
      end
C---------------------------------------------------------------------
      subroutine get3co4(nfact,coeffi,coeffbl,coeffbr,coeffel,coeffer)
      implicit none
      integer max
      parameter ( max = 32 )
      real*8 coeffi(4,max),coeffbl(4,max),coeffbr(4,max)
      real*8 coeffel(4,max),coeffer(4,max)
      integer nfm1,nfact,l

      nfm1 = nfact - 1
      do l = 1 , nfm1
C interior
         coeffi(1,l)=-(-3.0e0*l*l*nfact+l**3+2.0e0*l*nfact**2)/
     .                  6.0e0/nfact**3
         coeffi(2,l)=-(-6.0e0*nfact**3+6.0e0*l*l*nfact+3.0e0*l*nfact**2
     .                  -3.0e0*l**3)/6.0e0/nfact**3
         coeffi(3,l)=-(-6.0e0*l*nfact**2-3.0e0*l*l*nfact+3.0e0*l**3)/
     .                  6.0e0/nfact**3
         coeffi(4,l)=-(l*nfact**2-l**3)/6.0e0/nfact**3
C right edge
         coeffbr(1,l)= -1.0e0*(l**3 - l*nfact**2)/6.0e0/nfact**3
         coeffbr(2,l)=-(-3.0e0*l*l*nfact+6.0e0*l*nfact**2-3.0e0*l**3 )
     .                /6.0e0/nfact**3 
         coeffbr(3,l)= -(-6.0e0*nfact**3-3.0e0*l*nfact**2+
     .                    6.0e0*l*l*nfact+3.0e0*l**3)/6.0e0/nfact**3
         coeffbr(4,l)= -(-3.0e0*l*l*nfact-2.0e0*l*nfact**2-l**3)/
     .                 6.0e0/nfact**3
C left edge
         coeffbl(1,l)=-(-6.0e0*nfact**3+l**3-6.0e0*l*l*nfact+
     .                 11.0e0*l*nfact**2)/6.0e0/nfact**3
         coeffbl(2,l)=-(-18.0e0*l*nfact**2-3.0e0*l*l*l+15.0e0*l*l*nfact)
     .                /6.0e0/nfact**3
         coeffbl(3,l)=-(3.0e0*l*l*l-12.0e0*l*l*nfact+9.0e0*l*nfact**2)
     .                /6.0e0/nfact**3
         coeffbl(4,l)=-(-2.0e0*l*nfact**2+3.0e0*l*l*nfact-l**3)
     .                /6.0e0/nfact**3
C extrapolation left edge
         coeffel(1,l)=(-6.0e0*nfact**3+l**3-6.0e0*l*l*nfact+
     .                 11.0e0*l*nfact**2)/6.0e0/nfact**3
         coeffel(4,l)=(-26.0e0*l*nfact**2+9.0e0*l*l*nfact-l**3
     .                 +24.0e0*nfact**3)/6.0e0/nfact**3
         coeffel(3,l)=(-36.0e0*nfact**3+3.0e0*l**3-24.0e0*l*l*nfact+
     .                  57.0e0*l*nfact**2)/6.0e0/nfact**3
         coeffel(2,l)=(-42.0e0*l*nfact**2+24.0e0*nfact**3-
     .                  3.0e0*l**3+21.0e0*l*l*nfact)/6.0e0/nfact**3
C extrapolation left edge
         coeffer(1,l)=-(l**3+3.0e0*l*l*nfact+2.0e0*l*nfact**2)
     .                  /6.0e0/nfact**3
         coeffer(2,l)=-(-9.0e0*l*nfact**2-3.0e0*l**3-12.0e0*l*l*nfact)
     .                  /6.0e0/nfact**3
         coeffer(3,l)=-(3.0e0*l**3+15.0e0*l*l*nfact+18.0e0*l*nfact**2)
     .                  /6.0e0/nfact**3
         coeffer(4,l)=-(-11.0e0*l*nfact**2-6.0e0*l*l*nfact-l**3-
     .                  6.0e0*nfact**3)/6.0e0/nfact**3
      end do

      return
      end
C---------------------------------------------------------------------
C
C     Returns cubically interpolated value at X = X0 + U * H.
C
C---------------------------------------------------------------------
C
      REAL*8 FUNCTION ngfpq3_sc(F0,F1,F2,F3,ch0,ch1,ch2,ch3,U)
      implicit none
C
         real*8      C1, C2, F0, F1, F2, F3, U
         real*8      ch0,ch1,ch2,ch3
         real*8      HALF,             THIRD
         PARAMETER ( HALF = 0.5e0,
     *               THIRD = 0.3333333333333333e0 )
C
         C1 = HALF * U * (U - 1.0e0)
         C2 = THIRD * C1 * (U - 2.0e0)

         ngfpq3_sc= F0*ch0*(1.0e0 - U + C1 - C2) +
     *            F1*ch1*(U - 2.0e0 * C1 + 3.0e0 * C2) +
     *            F2*ch2*(C1 - 3.0e0 * C2) +
     *            F3*ch3*C2
C
         RETURN
C
      END
C
C---------------------------------------------------------------------
      integer function d3pro1nnc( uc,uf,nxc,nyc,nzc,nxf,nyf,nzf,
     .                      charmc,charmf,order)
      implicit none
      integer nxc, nyc, nzc
      integer nxf, nyf, nzf
      real*8  uc(nxc,nyc,nzc),uf(nxf,nyf,nzf) 
      real*8  charmc(nxc,nyc,nzc),charmf(nxf,nyf,nzf) 
      real*8  badval,bigval
      integer order

      integer pil,l,d3pro1nnc_

      integer   ic , jc , kc
      integer   iff , jf, kf
      integer ierr,dnlintxyzc
      integer icharsten3

      d3pro1nnc = 0

      d3pro1nnc =  d3pro1nnc_(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,order,
     .                ierr,charmc,charmf)

      return
      end
C-------------------------------------------------------------------------
      integer function d3pro1nnc_(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,order,
     .                            ierr,chc,chf)
      implicit none
      integer  nxf,nyf,nzf
      integer  nxc,nyc,nzc
      integer order
      logical ccdimn

      real*8 uf(nxf,nyf,nzf), uc(nxc,nyc,nzc)
      real*8 chf(nxf,nyf,nzf), chc(nxc,nyc,nzc)

      integer wksp_size
      parameter ( wksp_size = 8000 )
C this size will suit must needs we can have a 20x1 refinement with
C 20th order 
      real*8 wksp(wksp_size)

      integer i,j,k
      integer iff,jf,kf,kff
      integer ic,jc,kc
      integer nxfact,nyfact,nzfact
      integer nxskip,nyskip,nzskip
      integer iskip,jskip,kskip 
      integer ierr
      integer d1nintnc,d1nintnec

      real*8 f1,f2

      integer maxsize
      parameter ( maxsize = 2049 )
      real*8 vtc(maxsize),vtf(maxsize),cthc(maxsize)

      ierr = 0
      d3pro1nnc_ = 0

      if( (.not.ccdimn(nxc,order)).or.(.not.ccdimn(nyc,order)).or.
     .    (.not.ccdimn(nzc,order))) then
         write(6,*) 'not enough points for d3pro1nn'
         d3pro1nnc_ = 1
         return
      end if


C figure out the refinement factors in each dimension 

      nxfact = (nxf-1)/(nxc-1)
      nyfact = (nyf-1)/(nyc-1)
      nzfact = (nzf-1)/(nzc-1)

      nxskip = nxfact - 1
      nyskip = nyfact - 1
      nzskip = nzfact - 1

      d3pro1nnc_  = d1nintnc(order,wksp,wksp_size,nxfact)
      if (d3pro1nnc_.eq.1) return
C

      do kc = 1 , nzc
         kf = nzfact*kc-nzskip
         do jc = 1 , nyc
            jf = nyfact*jc-nyskip
            call d3vputv(uc,vtc,nxc,nyc,nzc,jc,kc,1,1)
            call d3vputv(chc,cthc,nxc,nyc,nzc,jc,kc,1,1)
            d3pro1nnc_ =  d1nintnec( nxc,nxf,vtc,vtf,wksp,
     .                     wksp_size,order,nxfact,cthc )
            if (d3pro1nnc_.eq.1) return
            call dvput3v(vtf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill in the missing lines
      if (nyfact.ne.nxfact) then
         d3pro1nnc_ = d1nintnc(order,wksp,wksp_size,nyfact)
      end if

      do kc = 1 , nzc
         kf =  nzfact*kc-nzskip
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,kf,nyfact,2)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,kf,nyfact,2)
            d3pro1nnc_ =  d1nintnec( nxc,nxf,vtc,vtf,wksp,
     .                     wksp_size,order,nyfact,cthc )
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

      if (nzfact.ne.nyfact) then
         d3pro1nnc_ = d1nintnc(order,wksp,wksp_size,nzfact)
      end if
C now fill in the missing planes
      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,jf,nzfact,3)
            call d3vputv(chf,cthc,nxf,nyf,nzf,iff,jf,nzfact,3)
            d3pro1nnc_ =  d1nintnec( nxc,nxf,vtc,vtf,wksp,
     .                     wksp_size,order,nzfact,cthc )
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do
            
      return
      end
c-----------------------------------------------------------------------
c
c     Dimension checker ...
c
c-----------------------------------------------------------------------

      logical function ccdimn(n,order)

         implicit      none

         integer       n,order

         ccdimn=  mod(n,2) .eq. 1  .and.  n .gt. order

         return

      end
C 1xn refinements using nth order interpolation
C evaluates the stencils via lud decomposition 
C and then inversion
C for nx1 we had ....
C
C The size of the work array( mat ) = order * (order+1) *(order -1)
C this is for the lud matrix + the pivit matrix stored for
C the order -1 possibilities
C
C I have not yet figured out the 1xn nth order ...
C
c-----------------------------------------------------------------------
      integer function d1nintnc(order,mat,mat_size,nfact) 
      implicit none

      integer maxorder,order
      parameter ( maxorder = 20 )
      integer mat_size
      real*8 mat(mat_size)
      integer need_size
      integer construct_matrixnc
      integer nfact

      d1nintnc= 0

C      need_size = order*(order+1)*(order-1)
      need_size = order*(order+1)*(order+1)

      if (mat_size.lt.(need_size)) then
         write(6,*) 'Need to make a larger work array in d1nintnc'
         d1nintnc= 1
         return
      end if

      d1nintnc=  construct_matrixnc(order,mat,nfact)

      return
      end
C======================================================================
      integer function d1nintnec(n,nf,uc,uf,mat,mat_size,order,nfact,ch)
      implicit none
      integer n,nf,mat_size
      real*8 uc(n),uf(nf),ch(n)
      real*8 mat(mat_size)
      integer order
      integer nfact,d1nintnec_

      integer maxorder
      parameter ( maxorder = 20 )

      d1nintnec =  d1nintnec_(n,nf,uc,uf,mat,order,nfact,ch)

      return
      end
C======================================================================
      integer function d1nintnec_(n,nf,uc,uf,mat,order,nfact,ch)
      implicit none
      integer n,nf,order,nfact
      real*8 uc(n),uf(nf),ch(n)
      real*8 mat(order,order+1,nfact+1)

      integer maxorder
      parameter ( maxorder = 20 )
      real*8 input_vec(maxorder),coeff_vec(maxorder)

      integer maxnh,nh,l
      parameter ( maxnh = 100)
      integer start(maxnh),stop(maxnh)
      real*8 oldval,badval
      parameter ( badval = 1.0d20)

      integer ic,iff
      integer i,j,ip
      integer length
      integer k,ioff
      integer istop,istart
      integer iv
      integer lskip,nskip
      real*8 nfacti
      

      length = nint(order*0.5)
      ioff = mod(order,2)
      nskip = nfact - 1
      nfacti = 1.0e0/nfact
      d1nintnec_ = 0

      do iff = 1 , nf
         uf(iff) = 0.0e0
      end do

      do ic = 1 , n
         iff = nfact*ic - nskip
         uf(iff) = uc(ic)*ch(ic)
      end do

C now find the start and stop
      nh = 1
      start(1) = nh
      oldval = ch(1)
      stop(nh) = 0

      do i = 1 , n
         if (ch(i).ne.oldval) then
            if (ch(i).ne.1) then
               stop(nh) = i - 1
            else
               nh = nh + 1
               start(nh) = i
               stop(nh) = 0
            end if
         end if
         oldval = ch(i)
      end do

      if (stop(nh).eq.0) stop(nh) = n


      do l = 1 , nh
         if ((stop(l)-start(l)).lt.order-1) then
            d1nintnec_ = 1
            return
         end if

         do lskip = 1 , nskip

            k = 1

            istart = nfact*(length-ioff)+(lskip-nskip) + 
     .               nfact*start(l)-nfact
            istop  = nfact*stop(l)-nskip-(length-1)*nfact-(nskip-lskip)

            do i = 1 , length
               do iff = istart,istop , nfact
                  iv= iff+nfact*(i-1)+(nfact-lskip)
                  ic = nint((iv+nskip)*nfacti)
                  uf(iff) = uf(iff) + mat(i,k,lskip)*
     .                      uf(iv)*ch(ic)
               end do
            end do

            j = i
            ip = 0

            do i = j , order
               ip = ip + 1

               do iff =istart,istop , nfact
                  iv = iff-nfact*(ip-1)-lskip
                  ic = nint((iv+nskip)*nfacti)
                  uf(iff) = uf(iff) + mat(i,k,lskip)*
     .                      uf(iv) * ch(ic) 

               end do
            end do

C                            BOUNDARIES

C        now lets get the end points ............
C        first do all the positive terms that we must modify

            istop = 1
            istart = length 
   
            do iff=stop(l)*nfact-nskip-(nfact-lskip),
     .             stop(l)*nfact-nskip-(length-1)*nfact+lskip,-nfact
               uf(iff) = 0.0e0
               k = k + 1
               i = 1
   
               do i = 1 , length
                  ip = iff + nfact*i - lskip 
                  if (ip.le.stop(l)*nfact-nskip) then
                     ic = nint((ip+nskip)*nfacti)
                     uf(iff) = uf(iff) + mat(i,k,lskip)*uf(ip)*ch(ic)
                  end if
               end do
   
               j = i - 1 
               do i = 1 , length -ioff
                  j = j + 1
                  ip = iff - nfact*i + (nfact-lskip )
                  if (ip.gt.start(l)*nfact-nskip) then
                     ic = nint((ip+nskip)*nfacti)
                     uf(iff) = uf(iff) + mat(j,k,lskip)*uf(ip)*ch(ic)
                  end if
               end do
      
               istop = istop + 1
               ip = 0
   
               ic = iff - nfact*(length-ioff)  - lskip + nfact
               do i = istart, istop , -1 
                  ic = ic - nfact 
                  iv = nint( (nskip+ic)*nfacti)
                  uf(iff) = uf(iff) + mat(i,k,lskip)*uf(ic)*ch(iv)
               end do
            end do
   
C        now let get the other end points
            istop = 1
            do iff = start(l)*nfact+lskip-nskip,
     .         (length-ioff-2)*nfact+start(l)*nfact+lskip-nskip,nfact 
               uf(iff) = 0.0e0
               k = k + 1
   
               do i = 1 , length
                  ip = iff + nfact*i - lskip 
                  if (ip.le.nfact*stop(l)-nskip) then
                     ic = nint( (nskip+ip)*nfacti)
                     uf(iff) = uf(iff) + mat(i,k,lskip)*uf(ip)*ch(ic)
                  end if
               end do

               j = i - 1 
               do i = 1 , length
                  j = j + 1
                  ip = iff - nfact*i + (nfact-lskip )
                  if (ip.ge.nfact*start(l)-nskip) then
                     ic = nint( (nskip+ip)*nfacti)
                     uf(iff) = uf(iff) + mat(j,k,lskip)*uf(ip)*ch(ic)
                  end if
               end do

               istop = istop + 1
               j = iff  + nfact*(length-1) +(nfact-lskip)

               do i = order , order -  (length-istop) +ioff, -1
                  j = j + nfact 
                  ic =  nint( (nskip+j)*nfacti)
                  uf(iff) = uf(iff) + mat(i,k,lskip)*uf(j)*ch(ic)
               end do
            end do
         end do
      end do


      do l = 1 , nh
         do lskip = 1 , nskip
C now do the extrapolation points

            k = order - 1
            if (nh.gt.1) then
               k = k + 1
               if (start(l).gt.1) then
                  iff = nfact*(start(l)-1) - nskip + lskip

                  uf(iff) = 0.0e0
                  do i = 1 , order
                     uf(iff) = uf(iff) + mat(i,k,lskip)*
     .                         uf(iff+i+nfact-lskip-1)
                  end do
               end if

               k = k + 1
               if (stop(l).gt.1.and.stop(l).lt.n) then
                  iff = nfact*(stop(l)) - nskip + lskip
                  do i = 1 , order
                     uf(iff) = uf(iff) + mat(i,k,lskip)*
     .                         uf(iff-i-lskip+1)
                  end do
               end if
            end if


         end do
      end do

      return
      end
C======================================================================
      integer function construct_matrixnc(order,mat,nfact)
      implicit none
      integer order
      integer nfact
      real*8 mat(order,order+1,nfact-1)
      integer maxorder
      parameter ( maxorder = 20 )
      real*8 mat2d(maxorder,maxorder),mat2ds(maxorder,maxorder) 
      real*8 work(maxorder)
      integer ipvt(maxorder)
      integer construct_matrixnc_

      construct_matrixnc= construct_matrixnc_(order,mat2d,mat2ds,
     .                                      ipvt,mat,work,nfact)

      return
      end
C======================================================================
      integer function construct_matrixnc_(order,mat2d,mat2ds,
     .        ipvt,mat,work,nfact)
      implicit none

      integer order
      integer nfact,nskip
C mat is stored as the coefficients and for each stencil
C 
      real*8 mat(order,order+1,nfact-1)
      real*8 mat2d(order,order),mat2ds(order,order)
      real*8 factorial
      real*8 work(order),det
      integer fa
      integer ipvt(order)
      integer length
      integer ifail
      integer ioff

      integer i,j,k,l,jj,lskip
      integer iendp1

      construct_matrixnc_ = 0
      nskip = nfact - 1 

      length = nint(order*0.5)
      ioff = mod(order,2)

      do lskip = 1 , nskip                        ! loop over each type
         k = 1

C        do the positive direction first ...
         fa = -lskip
         do i = 1 , length 
            fa = fa + nfact 
            do j = 1 , order
               mat2d(i,j) = (1.0e0*fa**(j-1))/factorial(j-1) 
            end do
         end do
   
         fa = nfact - lskip 
         iendp1 = i
         jj = i

C        now do the negative direction ....

         do i = jj , order
            fa = fa - nfact   
            do j = 1 , order
               mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1) 
            end do
         end do

         call copy1d(order*order,mat2d,mat2ds)

C        now do a lud decomposition on this and then an iverse and then
C        store the coefficients

         call sgetrf(mat2d,order,order,ipvt,ifail)
         construct_matrixnc_  = ifail
         call sgetri(mat2d,order,order,ipvt,det,work,01)

C        now get out the coefficients .....


         construct_matrixnc_ = ifail
         if (ifail.ne.0) return

C        now copy1d into mat array

         call copy1d_mat2d_mat(order,lskip,mat2d,k,mat,nskip)

C        now lets get the other stencils starting from the second point
C        and then going forward ......
C        assume nx1 refinements ...

C        first modify all the positive  terms that we must modify ...

         do l = length  - 1  , 1 , -1    !(

C        l is the number of terms to modify

            fa = -lskip - nfact*(length-1-ioff) 
            k = k + 1
            call copy1d(order*order,mat2ds,mat2d)

C        i is which term in the matrix ....

            do i = length , length - l + 1  , -1
               fa = fa - nfact 
               do j = 1 , order
                  mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
               end do
            end do

C now do a lud decomposition on this ........

            call sgetrf(mat2d,order,order,ipvt,ifail)
            construct_matrixnc_  = ifail
            call sgetri(mat2d,order,order,ipvt,det,work,01)

            construct_matrixnc_ = ifail
            if (ifail.ne.0) return

C now copy1d into mat array
            call copy1d_mat2d_mat(order,lskip,mat2d,k,mat,nskip)

         end do                         !)
C now modify the negative terms .....

         do l = length  - 1  , 1 + ioff , -1

C l is the number of terms to modify
C            fa = -(-(order-2)*nfact + nfact - (nfact-lskip) )

            fa = -lskip + length*nfact
            k = k + 1
            call copy1d(order*order,mat2ds,mat2d)
C i is which term in the matrix ....
            do i = order , order - l + 1  +ioff , -1
               fa = fa + nfact 
               do j = 1 , order
                  mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
               end do
            end do

C now do a lud decomposition on this ........

            call sgetrf(mat2d,order,order,ipvt,ifail)
            construct_matrixnc_  = ifail
            call sgetri(mat2d,order,order,ipvt,det,work,01)

            construct_matrixnc_ = ifail
            if (ifail.ne.0) return

C now copy1d into mat array
            call copy1d_mat2d_mat(order,lskip,mat2d,k,mat,nskip)

         end do

C now get the extrapolation points

            k = k + 1
            fa = 1.0e0*(nfact-lskip)
            do i = 1 , order
               do j = 1 , order
                  mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
               end do
               fa = fa + 1.0e0
            end do
C now do a lud decomposition on this ........
            call sgetrf(mat2d,order,order,ipvt,ifail)
            call sgetri(mat2d,order,order,ipvt,det,work,01)

            construct_matrixnc_ = ifail
            if (ifail.ne.0) return
C now copy into mat array
            call copy1d_mat2d_mat(order,lskip,mat2d,k,mat,nskip)

            k = k + 1
            fa = -1.0e0*lskip
            do i = 1 , order
               do j = 1 , order
                  mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
               end do
               fa = fa - 1.0e0
            end do
C now do a lud decomposition on this ........
            call sgetrf(mat2d,order,order,ipvt,ifail)
            call sgetri(mat2d,order,order,ipvt,det,work,01)

            construct_matrixnc_ = ifail
            if (ifail.ne.0) return
C now copy into mat array
            call copy1d_mat2d_mat(order,lskip,mat2d,k,mat,nskip)


      end do ! for the skip factor

      return
      end
C======================================================================
      subroutine copy1d(n,ain,aout)
      implicit none
      integer n
      real*8 ain(n),aout(n)

      integer i,j

      do i = 1 , n
         aout(i) = ain(i)
      end do

      return
      end
C------------------------------------------------------------------------------
C the routines to call are:
C  d3pro1n2(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc)                  nx1    second order
C  d3pro1n4(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc)                  nx1    fourth order
C  d3pro1nn(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,order)            nx1    nth    order
C------------------------------------------------------------------------------


C prolongs from a 1 to n grid using second order techinques
C i.e. cubic interpolation
C

      integer function  d3pro1n2(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc)
      implicit none
      integer  nxf,nyf,nzf
      integer  nxc,nyc,nzc

      real*8 uf(nxf,nyf,nzf), uc(nxc,nyc,nzc)

      integer i,j,k
      integer iff,jf,kf,kff
      integer ic,jc,kc
      integer nxfact,nyfact,nzfact
      integer nxskip,nyskip,nzskip
      integer iskip,jskip,kskip 

      real*8 vtc(1025),vtf(1025)

C figure out the refinement factors in each dimension 

      d3pro1n2 = 0
      if (nxf.gt.1025.or.nyf.gt.1025.or.nzf.gt.1025) then
         d3pro1n2 = 1
         return
      end if

      nxfact = (nxf-1)/(nxc-1)
      nyfact = (nyf-1)/(nyc-1)
      nzfact = (nzf-1)/(nzc-1)

      nxskip = nxfact - 1
      nyskip = nyfact - 1
      nzskip = nzfact - 1
C

      do kc = 1 , nzc
         kf = nzfact*kc-nzskip
         do jc = 1 , nyc
            jf = nyfact*jc-nyskip
            call d3vputv(uc,vtc,nxc,nyc,nzc,jc,kc,1,1)
            call dvi1q1(vtc,vtf,nxfact,nxc)
            call dvput3v(vtf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill in the missing lines

      do kc = 1 , nzc
         kf =  nzfact*kc-nzskip
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,kf,nyfact,2)
            call dvi1q1(vtc,vtf,nyfact,nyc) 
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

C now fill in the missing planes
      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,jf,nzfact,3)
            call dvi1q1(vtc,vtf,nzfact,nzc) 
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do
            
      return
      end
C-----------------------------------------------------------------------------C
      subroutine  DVI1Q1(V1,V2,IRHO,n1)
C
      IMPLICIT    none
C
      INTEGER      IRHO, N1,RHOM1
      REAL*8       V1(N1),     V2(1)
C
      REAL*8       C0, C1, C2, C3, SIGMA
      INTEGER      I, iskip
      integer iff,nxf
      integer irhom1
      integer max
      parameter ( max = 32 )
      real*8 f1(max),f2(max)

      nxf = iskip * n1
      irhom1 = irho - 1 

C first thing is to inject where we can inject
      do i = 1 , n1
         iff = irho*i - irhom1 
         v2(iff) = v1(i)
      end do

      do iskip = 1 , irhom1
         f2(iskip) = 1.0e0*(-iskip+irho)/irho
         f1(iskip) =  1.0e0*iskip/irho
      end do

      do iskip = 1 , irhom1
         do i = 1 , n1-1
            iff = irho*(i-1) + iskip+1
            v2(iff) = f2(iskip)*v1(i) + f1(iskip)*v1(i+1) 
         end do
      end do

      return
      end

C prolongs from a 1 to n grid using fourth order techinques
C i.e. cubic interpolation
C

      integer function  d3pro1n4(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc)
      implicit none
      integer  nxf,nyf,nzf
      integer  nxc,nyc,nzc

      real*8 uf(nxf,nyf,nzf), uc(nxc,nyc,nzc)

      integer i,j,k
      integer iff,jf,kf,kff
      integer ic,jc,kc
      integer nxfact,nyfact,nzfact
      integer nxskip,nyskip,nzskip
      integer iskip,jskip,kskip 

      real*8 f1,f2

      real*8 vtc(1025),vtf(1025)

      d3pro1n4 = 0

      if (nxf.gt.1025.or.nyf.gt.1025.or.nzf.gt.1025) then
         d3pro1n4 = 1
         return
      end if
C figure out the refinement factors in each dimension 

      nxfact = (nxf-1)/(nxc-1)
      nyfact = (nyf-1)/(nyc-1)
      nzfact = (nzf-1)/(nzc-1)

      nxskip = nxfact - 1
      nyskip = nyfact - 1
      nzskip = nzfact - 1
C

      do kc = 1 , nzc
         kf = nzfact*kc-nzskip
         do jc = 1 , nyc
            jf = nyfact*jc-nyskip
            call d3vputv(uc,vtc,nxc,nyc,nzc,jc,kc,1,1)
            call dvcuq1_s(vtc,vtf,nxfact,nxc)
            call dvput3v(vtf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill in the missing lines

      do kc = 1 , nzc
         kf =  nzfact*kc-nzskip
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,kf,nyfact,2)
            call dvcuq1_s(vtc,vtf,nyfact,nyc) 
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

C now fill in the missing planes
      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,jf,nzfact,3)
            call dvcuq1_s(vtc,vtf,nzfact,nzc) 
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do
            
      return
      end
C
C---------------------------------------------------------------------
C
C     Returns FTOC : 1 cubic refinement of V1 in V2 (user's respons-
C     ibility to ensure that V1 contains a minimum of 4 elements).
C
C---------------------------------------------------------------------
C
      SUBROUTINE dvcuq1_s(V1,V2,FTOC,N1)
      implicit none
C
         real*8      ngfpq3_s
         integer   n1,ftoc
         real*8      V1(n1), V2(n1*ftoc)
         real*8      FTOCR
         INTEGER   FTOCP1, I1, I2, I, J, N1M2
C
         N1M2 = N1 - 2
         FTOCP1 = FTOC + 1
         I2 = 1
         FTOCR = 1.0e0 / FTOC


         DO 10 J = 1 , FTOC
            V2(I2) = ngfpq3_s(V1(1),V1(2),V1(3),V1(4),(J-1)*FTOCR)
            I2 = I2 + 1
 10      CONTINUE
         DO 30 I1 = 2 , N1M2
            DO 20 J = 1 , FTOC
               V2(I2) = ngfpq3_s(V1(I1-1),V1(I1),V1(I1+1),V1(I1+2),
     *                         1.0e0 + (J-1)*FTOCR)
               I2 = I2 + 1
 20         CONTINUE
 30      CONTINUE
         DO 40 J = 1 , FTOCP1
            V2(I2) = ngfpq3_s(V1(N1-3),V1(N1-2),V1(N1-1),V1(N1),
     *                      2.0e0 + (J-1)*FTOCR)
            I2 = I2 + 1
 40      CONTINUE
C
         RETURN
C
      END
C---------------------------------------------------------------------
C
C     Returns cubically interpolated value at X = X0 + U * H.
C
C---------------------------------------------------------------------
C
      REAL*8 FUNCTION ngfpq3_s(F0,F1,F2,F3,U)
      implicit none
C
         real*8      C1, C2, F0, F1, F2, F3, U
         real*8      HALF,             THIRD
         PARAMETER ( HALF = 0.5e0,
     *               THIRD = 0.3333333333333333e0 )
C
         C1 = HALF * U * (U - 1.0e0)
         C2 = THIRD * C1 * (U - 2.0e0)

         ngfpq3_s = F0 * (1.0e0 - U + C1 - C2) +
     *            F1 * (U - 2.0e0 * C1 + 3.0e0 * C2) +
     *            F2 * (C1 - 3.0e0 * C2) +
     *            F3 * C2
C
         RETURN
C
      END
C
C prolongs from a 1 to n grid using nth order techinques
C
C
C need a 1D nth order interpolation
C this routine will make sure that its a nx1 where n=nx=ny=nz ....
C
      integer function  d3pro1nn(nxf,nyf,nzf,nxc,nyc,nzc,uf,uc,order) 
      implicit none
      integer  nxf,nyf,nzf
      integer  nxc,nyc,nzc
      integer order
      logical ccdimn

      real*8 uf(nxf,nyf,nzf), uc(nxc,nyc,nzc)

      integer wksp_size
      parameter ( wksp_size = 8000 )
C this size will suit must needs we can have a 20x1 refinement with
C 20th order 
      real*8 wksp(wksp_size)

      integer i,j,k
      integer iff,jf,kf,kff
      integer ic,jc,kc
      integer nxfact,nyfact,nzfact
      integer nxskip,nyskip,nzskip
      integer iskip,jskip,kskip 
      integer d1nintn

      real*8 f1,f2

      real*8 vtc(1025),vtf(1025)

      d3pro1nn = 0

      if( (.not.ccdimn(nxc,order)).or.(.not.ccdimn(nyc,order)).or.
     .    (.not.ccdimn(nzc,order))) then
         write(6,*) 'not enough points for d3pro1nn'
         d3pro1nn = 1
         return
      end if
      if (nxf.gt.1025.or.nyf.gt.1025.or.nzf.gt.1025) then
         d3pro1nn = 1
         return
      end if


C figure out the refinement factors in each dimension 

      nxfact = (nxf-1)/(nxc-1)
      nyfact = (nyf-1)/(nyc-1)
      nzfact = (nzf-1)/(nzc-1)

      nxskip = nxfact - 1
      nyskip = nyfact - 1
      nzskip = nzfact - 1

      d3pro1nn = d1nintn(order,wksp,wksp_size,nxfact)
C

      do kc = 1 , nzc
         kf = nzfact*kc-nzskip
         do jc = 1 , nyc
            jf = nyfact*jc-nyskip
            call d3vputv(uc,vtc,nxc,nyc,nzc,jc,kc,1,1)
            call d1nintne( nxc,nxf,vtc,vtf,wksp,
     .                     wksp_size,order,nxfact )
            call dvput3v(vtf,uf,nxf,nyf,nzf,jf,kf,1)
         end do
      end do

C now fill in the missing lines
      if (nyfact.ne.nxfact) then
         d3pro1nn = d1nintn(order,wksp,wksp_size,nyfact)
      end if

      do kc = 1 , nzc
         kf =  nzfact*kc-nzskip
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,kf,nyfact,2)
            call d1nintne( nxc,nxf,vtc,vtf,wksp,
     .                     wksp_size,order,nyfact )
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,kf,2)
         end do
      end do

      if (nzfact.ne.nyfact) then
         d3pro1nn = d1nintn(order,wksp,wksp_size,nzfact)
      end if
C now fill in the missing planes
      do jf = 1 , nyf
         do iff = 1 , nxf
            call d3vputv(uf,vtc,nxf,nyf,nzf,iff,jf,nzfact,3)
            call d1nintne( nxc,nxf,vtc,vtf,wksp,
     .                     wksp_size,order,nzfact )
            call dvput3v(vtf,uf,nxf,nyf,nzf,iff,jf,3)
         end do
      end do
            
      return
      end
C 1xn refinements using nth order interpolation
C evaluates the stencils via lud decomposition 
C and then inversion
C for nx1 we had ....
C
C The size of the work array( mat ) = order * (order+1) *(order -1)
C this is for the lud matrix + the pivit matrix stored for
C the order -1 possibilities
C
C I have not yet figured out the 1xn nth order ...
C
      integer function d1nintn(order,mat,mat_size,nfact) 
      implicit none

      integer maxorder,order
      parameter ( maxorder = 20 )
      integer mat_size
      real*8 mat(mat_size)
      integer need_size
      integer construct_matrixn
      integer nfact

      d1nintn = 0

      need_size = order*(order+1)*(order-1)

      if (mat_size.lt.(need_size)) then
         write(6,*) 'Need to make a larger work array in d1nintn'
         d1nintn = 1
         return
      end if

      d1nintn =  construct_matrixn(order,mat,nfact)

      return
      end
C======================================================================
      subroutine d1nintne(n,nf,uc,uf,mat,mat_size,order,nfact)
      implicit none
      integer n,nf,mat_size
      real*8 uc(n),uf(nf)
      real*8 mat(mat_size)
      integer order
      integer nfact

      integer maxorder
      parameter ( maxorder = 20 )

      call d1nintne_(n,nf,uc,uf,mat,order,nfact)

      return
      end
C======================================================================
      subroutine d1nintne_(n,nf,uc,uf,mat,order,nfact)
      implicit none
      integer n,nf,order,nfact
      real*8 uc(n),uf(nf)
      real*8 mat(order,order-1,nfact-1)

      integer maxorder
      parameter ( maxorder = 20 )
      real*8 input_vec(maxorder),coeff_vec(maxorder)

      integer ic,iff
      integer i,j,ip
      integer length
      integer k,ioff
      integer istop,istart
      integer lskip,nskip
      

      length = nint(order*0.5)
      ioff = mod(order,2)
      nskip = nfact - 1

      do iff = 1 , nf
         uf(iff) = 0.0e0
      end do

      do ic = 1 , n
         iff = nfact*ic - nskip
         uf(iff) = uc(ic)
      end do


      do lskip = 1 , nskip
C        first do the interior points ......... k = 1 

         k = 1

         istart = nfact*(length-ioff)+(lskip-nskip)
         istop  = nf - length*nfact + lskip 

         do i = 1 , length
            do iff = istart,istop , nfact
               uf(iff) = uf(iff) + mat(i,k,lskip)*
     .                   uf(iff+nfact*(i-1)+(nfact-lskip)) 
            end do
         end do

         j = i
         ip = 0

         do i = j , order
            ip = ip + 1

            do iff =istart,istop , nfact
               uf(iff) = uf(iff) + mat(i,k,lskip)*
     .                   uf(iff-nfact*(ip-1)-lskip) 

            end do
         end do

C         goto 1000

C                            BOUNDARIES

C        now lets get the end points ............
C        first do all the positive terms that we must modify

         istop = 1
         istart = length 

         do iff=nf-(nfact-lskip),nf-(length-1)*nfact+lskip,-nfact
            uf(iff) = 0.0e0
            k = k + 1
            i = 1

            do i = 1 , length
               ip = iff + nfact*i - lskip 
               if (ip.le.nf) then
                  uf(iff) = uf(iff) + mat(i,k,lskip)*uf(ip) 
               end if
            end do
   
            j = i - 1 
            do i = 1 , length -ioff
               j = j + 1
               ip = iff - nfact*i + (nfact-lskip )
               if (ip.gt.1) then
                  uf(iff) = uf(iff) + mat(j,k,lskip)*uf(ip) 
               end if
            end do
   
            istop = istop + 1
            ip = 0

            ic = iff - nfact*(length-ioff)  - lskip + nfact
            do i = istart, istop , -1 
               ic = ic - nfact 
               uf(iff) = uf(iff) + mat(i,k,lskip)*uf(ic) 
            end do
         end do

C        now let get the other end points

         istop = 1
         do iff = 1+lskip, (length-ioff-2)*nfact+1+lskip,nfact 
            uf(iff) = 0.0e0
            k = k + 1
C            i = 1
C            uf(iff)  = uf(iff) +  mat(i,k,lskip)*uf(iff+1)

            do i = 1 , length
               ip = iff + nfact*i - lskip 
               if (ip.le.nf) then
                  uf(iff) = uf(iff) + mat(i,k,lskip)*uf(ip)
               end if
            end do

            j = i - 1 
            do i = 1 , length
               j = j + 1
               ip = iff - nfact*i + (nfact-lskip )
               if (ip.ge.1) then
                  uf(iff) = uf(iff) + mat(j,k,lskip)*uf(ip)
               end if
            end do

            istop = istop + 1
            j = iff  + nfact*(length-1) +(nfact-lskip)

            do i = order , order -  (length-istop) +ioff, -1
               j = j + nfact 
               uf(iff) = uf(iff) + mat(i,k,lskip)*uf(j) 
            end do

         end do

1000   continue
      end do


      return
      end
C======================================================================
      integer function construct_matrixn(order,mat,nfact)
      implicit none
      integer order
      integer nfact
      real*8 mat(order,order-1,nfact-1)
      integer maxorder
      parameter ( maxorder = 20 )
      real*8 mat2d(maxorder,maxorder),mat2ds(maxorder,maxorder) 
      real*8 work(maxorder)
      integer ipvt(maxorder)
      integer construct_matrixn_

      construct_matrixn= construct_matrixn_(order,mat2d,mat2ds,
     .                                      ipvt,mat,work,nfact)

      return
      end
C======================================================================
      integer function construct_matrixn_(order,mat2d,mat2ds,
     .        ipvt,mat,work,nfact)
      implicit none

      integer order
      integer nfact,nskip
C mat is stored as the coefficients and for each stencil
C 
      real*8 mat(order,order-1,nfact-1)
      real*8 mat2d(order,order),mat2ds(order,order)
      real*8 factorial
      real*8 work(order),det
      integer fa
      integer ipvt(order)
      integer length
      integer ifail
      integer ioff

      integer i,j,k,l,jj,lskip
      integer iendp1

      construct_matrixn_ = 0
      nskip = nfact - 1 

      length = nint(order*0.5)
      ioff = mod(order,2)

      do lskip = 1 , nskip                        ! loop over each type
         k = 1

C        do the positive direction first ...
         fa = -lskip
         do i = 1 , length 
            fa = fa + nfact 
            do j = 1 , order
               mat2d(i,j) = (1.0e0*fa**(j-1))/factorial(j-1) 
            end do
         end do
   
         fa = nfact - lskip 
         iendp1 = i
         jj = i

C        now do the negative direction ....

         do i = jj , order
            fa = fa - nfact   
            do j = 1 , order
               mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1) 
            end do
         end do

         call copy1d(order*order,mat2d,mat2ds)

C        now do a lud decomposition on this and then an iverse and then
C        store the coefficients

         call sgetrf(mat2d,order,order,ipvt,ifail)
         construct_matrixn_  = ifail
         call sgetri(mat2d,order,order,ipvt,det,work,01)

C        now get out the coefficients .....


         construct_matrixn_ = ifail
         if (ifail.ne.0) return

C        now copy1d into mat array

         call copy1d_mat2d_matnc(order,lskip,mat2d,k,mat,nskip)

C        now lets get the other stencils starting from the second point
C        and then going forward ......
C        assume nx1 refinements ...

C        first modify all the positive  terms that we must modify ...

         do l = length  - 1  , 1 , -1    !(

C        l is the number of terms to modify

            fa = -lskip - nfact*(length-1-ioff) 
            k = k + 1
            call copy1d(order*order,mat2ds,mat2d)

C        i is which term in the matrix ....

            do i = length , length - l + 1  , -1
               fa = fa - nfact 
               do j = 1 , order
                  mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
               end do
            end do

C now do a lud decomposition on this ........

            call sgetrf(mat2d,order,order,ipvt,ifail)
            construct_matrixn_  = ifail
            call sgetri(mat2d,order,order,ipvt,det,work,01)

            construct_matrixn_ = ifail
            if (ifail.ne.0) return

C now copy1d into mat array
            call copy1d_mat2d_matnc(order,lskip,mat2d,k,mat,nskip)

         end do                         !)
C now modify the negative terms .....

         do l = length  - 1  , 1 + ioff , -1

C l is the number of terms to modify
C            fa = -(-(order-2)*nfact + nfact - (nfact-lskip) )

            fa = -lskip + length*nfact
            k = k + 1
            call copy1d(order*order,mat2ds,mat2d)
C i is which term in the matrix ....
            do i = order , order - l + 1  +ioff , -1
               fa = fa + nfact 
               do j = 1 , order
                  mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
               end do
            end do

C now do a lud decomposition on this ........

            call sgetrf(mat2d,order,order,ipvt,ifail)
            construct_matrixn_  = ifail
            call sgetri(mat2d,order,order,ipvt,det,work,01)

            construct_matrixn_ = ifail
            if (ifail.ne.0) return

C now copy1d into mat array
            call copy1d_mat2d_matnc(order,lskip,mat2d,k,mat,nskip)

         end do

      end do ! for the skip factor

      return
      end
C======================================================================
      subroutine copy1d_mat2d_mat(order,lskip,mat2d,k,mat,nskip)
      implicit none

      integer order
      integer lskip,nskip
      real*8 mat(order,order+1,nskip)
      real*8 mat2d(order,order)
      integer k

      integer i,j

C      write(80,'(2i5)') k,lskip

      do i = 1 , order
         mat(i,k,lskip) = mat2d(1,i) 
      end do
C      write(80,'(1P,8e10.2)')  (mat(i,k,lskip),i=1,order)

      return
      end
C======================================================================
      subroutine copy1d_mat2d_matnc(order,lskip,mat2d,k,mat,nskip)
      implicit none

      integer order
      integer lskip,nskip
      real*8 mat(order,order-1,nskip)
      real*8 mat2d(order,order)
      integer k

      integer i,j

C      write(80,'(2i5)') k,lskip

      do i = 1 , order
         mat(i,k,lskip) = mat2d(1,i)
      end do
C      write(80,'(1P,8e10.2)')  (mat(i,k,lskip),i=1,order)

      return
      end
C======================================================================
      subroutine err(name)
      implicit none
      character*(*) name

      write(6,*) name
      write(6,*) 'Please call Scott Klasky @471-4700 about problem'

      return
      end
C-----------------------------------------------------------------------------
C  given a list of points these routines interpolate to a solution vector
C
C  d2lintxyz(nx,ny,nz,x,y,z,pil,listx,listy,listz,u,soln)       list    order 2
C  d4lintxyz(nx,ny,nz,x,y,z,pil,listx,listy,listz,u,soln)       list    order 4
C  dnlintxyz(nx,ny,nz,x,y,z,pil,listx,listy,listz,u,soln,order) list    order n
C-----------------------------------------------------------------------------
C This routine takes a list of points and interpolates to this to second
C order
      integer function d2lintxyz(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln)
CFPP$ EXPAND(getijk)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 u(nx,ny,nz)
      real*8 soln(pil)

      integer order
      parameter ( order = 2 )

C local variables
      real*8 dx,dy,dz
      real*8 dxi,dyi,dzi
      integer l
      integer i,j,k
      real*8 dx0,dy0,dz0
      integer ierr, getijk
      real*8 vx(4),vx2(2)


      d2lintxyz = 0

      dx = abs( x(2)-x(1) )
      dy = abs( y(2)-y(1) )
      dz = abs( z(2)-z(1) )
      dxi= 1.0e0/dx
      dyi= 1.0e0/dy
      dzi= 1.0e0/dz

C this will make a bounding box around each point and then interpolate
C in a series of 1d interpolations 


      do l = 1 , pil

C step 1 is to get the 2 points that surround this point in the x,y,z directions
C it returns the i-1,j-1,k-1 point and the dx0ff ...

         ierr =  getijk(nx,x,listx(l),i,dx0)
         if (ierr.eq.1) goto 1000
         ierr =  getijk(ny,y,listy(l),j,dy0)
         if (ierr.eq.1) goto 1000
         ierr =  getijk(nz,z,listz(l),k,dz0)
         if (ierr.eq.1) goto 1000

C now form a cube around this point ......
C interpolation formula is
C
C f = (h-a)/h * f(x-a) + a/h * f(x+h-a) 
C
C and we are given the a as in dx0 = a

         vx(1) =  (dx-dx0)*dxi*u(i  ,j  ,k  )+dx0*dxi*u(i+1,j  ,k  )
         vx(2) =  (dx-dx0)*dxi*u(i  ,j+1,k  )+dx0*dxi*u(i+1,j+1,k  )
         vx(3) =  (dx-dx0)*dxi*u(i  ,j  ,k+1)+dx0*dxi*u(i+1,j  ,k+1)
         vx(4) =  (dx-dx0)*dxi*u(i  ,j+1,k+1)+dx0*dxi*u(i+1,j+1,k+1)

         vx2(1) = (dy-dy0)*dyi*vx(1) + dy0*dyi*vx(2)
         vx2(2) = (dy-dy0)*dyi*vx(3) + dy0*dyi*vx(4)

         soln(l) = (dz-dz0)*dzi*vx2(1) + dz0*dzi*vx2(2)

      end do

      return
1000   continue

      d2lintxyz    = 1
      return
      end
C------------------------------------------------------------------------------
      integer function getijk(n,v,lv,l,del)
      implicit none
      integer n,l
      real*8 v(n)
      real*8 del,lv

      integer i

      getijk = 0

      do i = 1 , n
         if (v(i).ge.lv) goto 100
      end do

      getijk = 1
      return

100   continue

      l = i - 1
      if (l.lt.1) l = 1
      del = lv - v(l)

      if (del.lt.(0.0e0)) getijk = 1

      return
      end
C This routine takes a list of points and interpolates to this to fourth
C order
      integer function d4lintxyz(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln) 
CFPP$ EXPAND(getijk2)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 u(nx,ny,nz)
      real*8 soln(pil)

      integer order
      parameter ( order = 4 )

C local variables
      real*8 dx,dy,dz
      real*8 dxi,dyi,dzi
      integer l,lj,lk
      integer i,j,k
      real*8 dx0,dy0,dz0
      integer ierr, getijk2
      real*8 vx(4,4),vx2(4)
      real*8 den

      real*8 dx2,dx3,dx4,dx5
      real*8 dy2,dy3,dy4,dy5
      real*8 dz2,dz3,dz4,dz5

      real*8 dx02,dx03,dx04,dx05
      real*8 dy02,dy03,dy04,dy05
      real*8 dz02,dz03,dz04,dz05

      d4lintxyz = 0

      dx = abs( x(2)-x(1) )
      dy = abs( y(2)-y(1) )
      dz = abs( z(2)-z(1) )
      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz

      dx3 = dx**3
      dy3 = dy**3
      dz3 = dz**3

      dxi= 1.0e0/dx
      dyi= 1.0e0/dy
      dzi= 1.0e0/dz

C this will make a bounding box around each point and then interpolate
C in a series of 1d interpolations 


      do l = 1 , pil

C step 1 is to get the 2 points that surround this point in the x,y,z directions
C it returns the i-1,j-1,k-1 point and the dx0ff ...

         ierr =  getijk2(nx,x,listx(l),i,dx0)
         if (ierr.eq.1) goto 1000
         ierr =  getijk2(ny,y,listy(l),j,dy0)
         if (ierr.eq.1) goto 1000
         ierr =  getijk2(nz,z,listz(l),k,dz0)
         if (ierr.eq.1) goto 1000

         dx02 = dx0**2
         dx03 = dx0**3

         dy02 = dy0**2
         dy03 = dy0**3

         dz02 = dz0**2
         dz03 = dz0**3

C now form a cube around this point ......
C interpolation formula is
C
C f = (h-a)/h * f(x-a) + a/h * f(x+h-a) 
C
C and we are given the a as in dx0 = a

         den = 1.0e0/(6.0e0*dx3) 

         do lk = 1 , 4
            do lj = 1 , 4

               vx(lj,lk) =  den*(
     .            u(i+3,j+lj-1,k+lk-1)*(dx03+2.0e0*dx2*dx0-
     .                                  3.0e0*dx*dx02)+
     .            u(i+2,j+lj-1,k+lk-1)*(-3.0e0*dx03-9.0e0*dx2*dx0
     .                                  +12.0e0*dx*dx02) +
     .            u(i+1,j+lj-1,k+lk-1)*(3.0e0*dx03+18.0e0*dx2*dx0-
     .                                  15.0e0*dx*dx02)+
     .            u(i  ,j+lj-1,k+lk-1)*(-dx03+6.0e0*dx3-11.0e0*dx2*dx0+
     .                                  6.0e0*dx*dx02) )
            end do
         end do

         den = 1.0e0/(6.0e0*dy3) 
         do lk = 1 , 4

            vx2(lk) = den*(
     .            vx(4,lk)*(dy03+2.0e0*dy2*dy0-3.0e0*dy*dy02)+
     .            vx(3,lk)*(-3.0e0*dy03-9.0e0*dy2*dy0+12.0e0*dy*dy02)+
     .            vx(2,lk)*(3.0e0*dy03+18.0e0*dy2*dy0-15.0e0*dy*dy02)+
     .            vx(1,lk)*(-dy03+6.0e0*dy3-11.0e0*dy2*dy0 + 
     .                      6.0e0*dy*dy02) )
         end do

         den = 1.0e0/(6.0e0*dz3) 

         soln(l) = den*(
     .            vx2(4)*(dz03+2.0e0*dz2*dz0-3.0e0*dz*dz02) +
     .            vx2(3)*(-3.0e0*dz03-9.0e0*dz2*dz0+12.0e0*dz*dz02) +
     .            vx2(2)*(3.0e0*dz03+18.0e0*dz2*dz0-15.0e0*dz*dz02) +
     .            vx2(1)*(-dz03+6.0e0*dz3-11.0e0*dz2*dz0 +
     .                    6.0e0*dz*dz02) )


      end do

      return
1000   continue

      d4lintxyz    = 1
      return
      end
C------------------------------------------------------------------------------
      integer function getijk2(n,v,lv,l,del)
      implicit none
      integer n,l
      real*8 v(n)
      real*8 del,lv

      integer i

      getijk2 = 0

      do i = 1 , n
         if (v(i).ge.lv) goto 100
      end do

      getijk2 = 1
      return

100   continue

C now push this off two points ...
      l = i - 2
      if (l.lt.1) l = 1

      if (l+3.gt.n) l = n - 3

      del = lv - v(l)

      if (del.lt.(0.0e0)) getijk2 = 1

      return
      end
C This routine takes a list of points and interpolates to this to fourth
C order
      integer function dnlintxyz(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln,order) 
CFPP$ EXPAND(getijkn,getco)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 u(nx,ny,nz)
      real*8 soln(pil)

      integer order

C local variables
      real*8 dx,dy,dz
      real*8 dxi,dyi,dzi
      integer l,lj,lk,li
      integer i,j,k
      real*8 dx0,dy0,dz0
      integer ierr, getijkn
      integer maxorder
      parameter ( maxorder = 64 )
      real*8 mat2d(maxorder,maxorder)
      real*8 vx(maxorder,maxorder),vx2(maxorder)
      real*8 den,work(maxorder),sten(maxorder)
      integer getco

      real*8 dx2,dx3,dx4,dx5
      real*8 dy2,dy3,dy4,dy5
      real*8 dz2,dz3,dz4,dz5

      real*8 dx02,dx03,dx04,dx05
      real*8 dy02,dy03,dy04,dy05
      real*8 dz02,dz03,dz04,dz05

      dnlintxyz = 0

      dx = abs( x(2)-x(1) )
      dy = abs( y(2)-y(1) )
      dz = abs( z(2)-z(1) )
      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz

      dx3 = dx**3
      dy3 = dy**3
      dz3 = dz**3

      dxi= 1.0e0/dx
      dyi= 1.0e0/dy
      dzi= 1.0e0/dz

C this will make a bounding box around each point and then interpolate
C in a series of 1d interpolations 


      do l = 1 , pil

C step 1 is to get the 2 points that surround this point in the x,y,z directions
C it returns the i-1,j-1,k-1 point and the dx0ff ...

         ierr =  getijkn(nx,x,listx(l),i,dx0,order)
         if (ierr.eq.1) goto 1000
         ierr =  getijkn(ny,y,listy(l),j,dy0,order)
         if (ierr.eq.1) goto 1000
         ierr =  getijkn(nz,z,listz(l),k,dz0,order)
         if (ierr.eq.1) goto 1000

C now we must get a 1d matrix of coefficients .......
c do first for the x direction

         ierr = getco(mat2d,order,dx,dx0,sten)
         if (ierr.eq.1) goto 1000

         do lk = 1 , order
            do lj = 1 , order
               vx(lj,lk) = 0.0e0
            end do
               vx2(lk) = 0.0e0
         end do

         do lk = 1 , order 
            do lj = 1 , order 
               do li = 1 , order
                  vx(lj,lk)=sten(li)*u(i+li-1,j+lj-1,k+lk-1)+vx(lj,lk)
               end do
            end do
         end do

         ierr = getco(mat2d,order,dy,dy0,sten)
         if (ierr.eq.1) goto 1000

         do lk = 1 , order
            do li = 1 , order
               vx2(lk) = sten(li)*vx(li,lk) + vx2(lk)
            end do
         end do

         ierr = getco(mat2d,order,dz,dz0,sten)
         if (ierr.eq.1) goto 1000
         soln(l) = 0.0e0

         do li = 1 , order
            soln(l) = soln(l) + sten(li)*vx2(li)
         end do

      end do

      return

1000   continue

      dnlintxyz    = 1
      return
      end
C------------------------------------------------------------------------------
      integer function getijkn(n,v,lv,l,del,order)
      implicit none
      integer n,l
      real*8 v(n)
      real*8 del,lv

      integer i
      integer order

      getijkn = 0

      do i = 1 , n
         if (v(i).ge.lv) goto 100
      end do

      getijkn = 1
      return

100   continue

C now push this off order/2 points ...
      l = i - nint(order*0.50) 
      if (l.lt.1) l = 1

      if (l+(order-1).gt.n) l = n - (order-1) 

      del = lv - v(l)

      if (del.lt.(0.0e0)) getijkn = 1

      return
      end
C------------------------------------------------------------------------------
      integer function getco(mat2d,order,dx,dx0,sten)
      implicit none
      integer order
      real*8 dx,dx0
      real*8 mat2d(order,order)
      real*8 sten(order)
      integer lj,lk,ljs
      real*8 factorial

      integer maxorder
      parameter ( maxorder = 64 )
      integer ipvt(maxorder)
      real*8 work(maxorder)
      real*8 det
      integer ifail

      ljs = 0
      do lj = 1 , order 
         if ( abs((lj-1)*dx-dx0).lt.(1.0e-10)) then
            ljs = lj
         end if
      end do

      if (ljs.eq.0) then
         do lk = 1 , order
            do lj = 1 , order
               mat2d(lk,lj)=((lk-1)*dx-dx0)**(lj-1)/factorial(lj-1)
            end do
         end do

         call sgetrf(mat2d,order,order,ipvt,ifail)
         call sgetri(mat2d,order,order,ipvt,det,work,01)

         do lj = 1 , order
            sten(lj) = mat2d(1,lj)
         end do
      else
         do lj = 1 , order
            sten(lj) = 0.0e0
         end do

         sten(ljs) = 1.0e0
      end if

      getco = ifail

      return
      end
C======================================================================
c-----------------------------------------------------------------------
C routines in this section are
C       d2lintxyzc(nx,ny,nz,x,y,z,pil,lx,ly,lz,u,v,char,iwk,didit) 
C                                                           list  char 2nd O 
C       d4lintxyzc(nx,ny,nz,x,y,z,pil,lx,ly,lz,u,v,char,iwk,didit)
C                                                           list  char 4th O 
C       dnlintxyzc(nx,ny,nz,x,y,z,pil,lx,ly,lz,u,v,order,char,iwk,stwk,didit)
C                                                           list  char pth O 
c-----------------------------------------------------------------------
      integer function d2lintxyzc(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln,charm,iwk,didit)
CFPP$ EXPAND(getijknew,reset_offsets4,movepoint)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 u(nx,ny,nz),charm(nx,ny,nz)
      real*8 soln(pil)
      integer iwk(4*pil),moveit

      integer order,didit

C local variables
      real*8 dx,dy,dz
      integer l,lj,lk,li,ll
      integer i,j,k,iout
      real*8 dx0,dy0,dz0,dxout
      integer ierr
      real*8 mat2d(2,2)
      real*8 vx(4),vx2(2)
      real*8 den,work(2),sten(2)
      integer reset_offsets4
      integer movepoint,getijknew
      parameter ( order = 2)
      integer iw

      real*8 dxi,dyi,dzi

      d2lintxyzc= 0

      dx = abs( x(2)-x(1) )
      dy = abs( y(2)-y(1) )
      dz = abs( z(2)-z(1) )

      dxi = 1.0e0/dx
      dyi = 1.0e0/dy
      dzi = 1.0e0/dz

C this will make a bounding box around each point and then interpolate
C in a series of 1d interpolations


      iw = nint(0.5*order)

      do l = 1 , pil

         if(didit.eq.0) then
            iwk(3*pil+1) = 1
            moveit = 1
            d2lintxyzc =
     .      getijknew(nx,ny,nz,listx(l),listy(l),listz(l),x,y,z,
     .                     i,j,k,dx0,dy0,dz0,charm,order)

            if (d2lintxyzc.eq.1) return

            if (charm(i+iw,j+iw,k+iw).eq.1) then
               ierr = movepoint(nx,ny,nz,x,y,z,listx(l),listy(l),
     .                   listz(l),i,j,k,dx0,dy0,dz0,charm,order)

               if (ierr.eq.1) goto 1000
               moveit = 0
               iwk(3*pil+l) = 0
            end if

            iwk(l) = i
            iwk(1*pil+l) = j
            iwk(2*pil+l) = k
         else
            i = iwk(l)
            j = iwk(1*pil+l)
            k = iwk(2*pil+l)
            moveit = iwk(3*pil+l)
            dx0 = listx(l) - x(i)
            dy0 = listy(l) - y(j)
            dz0 = listz(l) - z(k)
         end if

         iout = i
         dxout = dx0

C
         if (moveit.eq.1) then
            ll = 0
            do lk = 1 , order
               do lj = 1 , order
                  ierr = reset_offsets4(nx,ny,nz,x,y,z,
     .             listx(l),charm,dxout,iout,j+lj-1,k+lk-1,i,dx0,order)
                   if (ierr.eq.1) goto 1000
                   ll = ll + 1
                   vx(ll) = (dx-dx0)*dxi*u(i,j+lj-1,k+lk-1) +
     .                      dx0*dxi*u(i+1,j+lj-1,k+lk-1)
               end do
            end do
         else

            ll = 0
            do lk = 1 , order
               do lj = 1 , order
                   ll = ll + 1
                   vx(ll) = (dx-dx0)*dxi*u(i,j+lj-1,k+lk-1) +
     .                      dx0*dxi*u(i+1,j+lj-1,k+lk-1)
               end do
            end do
         end if

         vx2(1) = (dy-dy0)*dyi*vx(1) + dy0*dyi*vx(2)
         vx2(2) = (dy-dy0)*dyi*vx(3) + dy0*dyi*vx(4)
         soln(l) = (dz-dz0)*dzi*vx2(1) + dz0*dzi*vx2(2)

      end do

      didit = didit + 1
      return

1000   continue

      didit = didit + 1
      d2lintxyzc   = 1
      return
      end
C======================================================================
      integer function d2lintxyzcd(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln,charm,iwk,swk,didit,wd) 
CFPP$ EXPAND(getsten_num)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 soln(pil)
      real*8 u(nx,ny,nz),charm(nx,ny,nz)
      integer iwk(*),swk(*)
      integer order
      parameter ( order = 2)
      character*1 wd

      integer d2lintxyzc
      integer didit
      integer l,m,lm
      integer maxs
      parameter ( maxs = 257 * 257 * 3)
      real*8 inlx(maxs),inly(maxs),inlz(maxs),ins(maxs)
      integer sten_num(maxs),getsten_num
      real*8 dx,dy,dz
      real*8 dxi,dyi,dzi
      real*8 coeff(5,3)
      integer ip

      d2lintxyzcd = 0

      if (pil.gt.257*257) then
         d2lintxyzcd = 1
         return
      end if

      if (wd.ne.' ') then
         call getcoeff2(coeff,2*order+1,order+1)
      end if


      if (wd.eq.' ') then

         d2lintxyzcd = d2lintxyzc(nx,ny,nz,x,y,z,pil,listx,listy,
     .                            listz,u,soln,charm,iwk,didit) 
         if (d2lintxyzcd.eq.1) return
      else if (wd.eq.'x') then
         dx = x(2) - x(1) 
         dxi= 1.0e0/dx 
         ip = 0
         do m = 1 , pil

            if (didit.eq.0) then
               sten_num(m)=getsten_num(nx,ny,nz,x,y,z,listx(m),
     .                         listy(m),listz(m),charm,order,1)
               swk(m) = sten_num(m)
            else
               sten_num(m) = swk(m)
            end if

            if (sten_num(m).eq.2) then
               l = -1
               do lm = 1,2
                  ip = ip + 1
                  inlx(ip) = listx(m) + l*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
                  l = l + 2
               end do
            else if (sten_num(m).eq.1) then
               do lm = 0,2
                  ip = ip + 1
                  inlx(ip) = listx(m) + lm*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
               end do
            else
               do lm = -2,0
                  ip = ip + 1
                  inlx(ip) = listx(m) + lm*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
               end do
            end if
         end do
         d2lintxyzcd = d2lintxyzc(nx,ny,nz,x,y,z,ip,
     .                    inlx,inly,inlz,u,ins,charm,iwk,didit) 
         if (d2lintxyzcd.eq.1) return
         ip = 0
         do m = 1 , pil
            ip = ip + 3
            if (sten_num(m).eq.2) then
               ip = ip - 1
               soln(m) = dxi*(coeff(4,2)*ins(ip)+coeff(2,2)*ins(ip-1)) 
            else if (sten_num(m).eq.1) then
               soln(m)=dxi*(coeff(3,1)*ins(ip-2)+coeff(4,1)*ins(ip-1)
     .                     +coeff(5,1)*ins(ip))
            else
               soln(m)=dxi*(coeff(1,3)*ins(ip-2)+coeff(2,3)*ins(ip-1)
     .                     +coeff(3,3)*ins(ip))
            end if
         end do

      else if (wd.eq.'y') then

         dy = y(2) - y(1)
         dyi= 1.0e0/dy
         ip = 0
         do m = 1 , pil

            if (didit.eq.0) then
               sten_num(m)=getsten_num(nx,ny,nz,x,y,z,listx(m),
     .                         listy(m),listz(m),charm,order,1)
               swk(m) = sten_num(m)
            else
               sten_num(m) = swk(m)
            end if

            if (sten_num(m).eq.2) then
               l = -1
               do lm = 1,2
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + l*dy
                  inlz(ip) = listz(m)
                  l = l + 2
               end do
            else if (sten_num(m).eq.1) then
               do lm = 0,2
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + lm*dy
                  inlz(ip) = listz(m)
               end do
            else
               do lm = -2,0
                  ip = ip + 1
                  inlx(ip) = listx(m)
                  inly(ip) = listy(m) + lm*dy
                  inlz(ip) = listz(m)
               end do
            end if
         end do
         d2lintxyzcd = d2lintxyzc(nx,ny,nz,x,y,z,ip,
     .                    inlx,inly,inlz,u,ins,charm,iwk,didit)
         if (d2lintxyzcd.eq.1) return
         ip = 0
         do m = 1 , pil
            ip = ip + 3
            if (sten_num(m).eq.2) then
               ip = ip - 1
               soln(m) = dyi*(coeff(4,2)*ins(ip)+coeff(2,2)*ins(ip-1))
            else if (sten_num(m).eq.1) then
               soln(m)=dyi*(coeff(3,1)*ins(ip-2)+coeff(4,1)*ins(ip-1)
     .                     +coeff(5,1)*ins(ip))
            else
               soln(m)=dyi*(coeff(1,3)*ins(ip-2)+coeff(2,3)*ins(ip-1)
     .                     +coeff(3,3)*ins(ip))
            end if
         end do

      else if (wd.eq.'z') then

         dz = z(2) - z(1)
         dzi= 1.0e0/dz
         ip = 0
         do m = 1 , pil

            if (didit.eq.0) then
               sten_num(m)=getsten_num(nx,ny,nz,x,y,z,listx(m),
     .                         listy(m),listz(m),charm,order,1)
               swk(m) = sten_num(m)
            else
               sten_num(m) = swk(m)
            end if

            if (sten_num(m).eq.2) then
               l = -1
               do lm = 1,2
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m) + l*dz
                  l = l + 2
               end do
            else if (sten_num(m).eq.1) then
               do lm = 0,2
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m) +lm*dz
               end do
            else
               do lm = -2,0
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m) + lm*dz
               end do
            end if
         end do
         d2lintxyzcd = d2lintxyzc(nx,ny,nz,x,y,z,ip,
     .                    inlx,inly,inlz,u,ins,charm,iwk,didit)
         if (d2lintxyzcd.eq.1) return
         ip = 0
         do m = 1 , pil
            ip = ip + 3
            if (sten_num(m).eq.2) then
               ip = ip - 1
               soln(m) = dzi*(coeff(4,2)*ins(ip)+coeff(2,2)*ins(ip-1))
            else if (sten_num(m).eq.1) then
               soln(m)=dzi*(coeff(3,1)*ins(ip-2)+coeff(4,1)*ins(ip-1)
     .                     +coeff(5,1)*ins(ip))
            else
               soln(m)=dzi*(coeff(1,3)*ins(ip-2)+coeff(2,3)*ins(ip-1)
     .                     +coeff(3,3)*ins(ip))
            end if
         end do



      else
         write(6,*) '*** This derivative is not implemented ***'
         d2lintxyzcd = 1
      end if

      return
      end
C======================================================================
      integer function d4lintxyzc(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln,charm,iwk,didit) 
CFPP$ EXPAND(getijknew,reset_offsets4,movepoint) 
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 u(nx,ny,nz),charm(nx,ny,nz)
      real*8 soln(pil)
      integer iwk(4*pil),moveit
      integer didit

      integer order

C local variables
      real*8 dx,dy,dz
      integer l,lj,lk,li,ll
      integer i,j,k
      real*8 dx0,dy0,dz0
      integer ierr, getijknc4
      real*8 vx(4,4),vx2(4),den
      integer reset_offsets4,getijknew
      integer movepoint
      parameter ( order = 4)

C      real*8 sum
      real*8 dxi,dyi,dzi
      real*8 dx02,dy02,dz02
      real*8 dx03,dy03,dz03
      real*8 dx20,dxx02
      real*8 dy20,dyy02
      real*8 dz20,dzz02
      real*8 dx2,dy2,dz2
      real*8 dx3,dy3,dz3
      integer iout,iw
      real*8 dxout

      d4lintxyzc= 0

      dx = abs( x(2)-x(1) )
      dy = abs( y(2)-y(1) )
      dz = abs( z(2)-z(1) )

      dxi = 1.0e0/dx
      dyi = 1.0e0/dy
      dzi = 1.0e0/dz

      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz

      dx3 = dx**3
      dy3 = dy**3
      dz3 = dz**3

      den = 1.0e0/(6.0e0*dx3)
      iw = nint(0.5*order)   

C this will make a bounding box around each point and then interpolate
C in a series of 1d interpolations 
      iw = nint(0.5*order)   
      
      do l = 1 , pil

         if (didit.eq.0) then
            iwk(3*pil+l) = 1
            moveit = 1
            d4lintxyzc =
     .      getijknew(nx,ny,nz,listx(l),listy(l),listz(l),
     .                x,y,z,i,j,k,dx0,dy0,dz0,charm,order)


            if (d4lintxyzc.eq.1) return

            if (charm(i+iw,j+iw,k+iw).eq.1) then
               ierr = movepoint(nx,ny,nz,x,y,z,listx(l),listy(l),
     .                     listz(l),i,j,k,dx0,dy0,dz0,charm,order)

               if (ierr.eq.1) goto 1000
               moveit = 0
               iwk(3*pil+l) = 0 
            end if

            iwk(l) = i
            iwk(1*pil+l) = j
            iwk(2*pil+l) = k
         else
            i = iwk(l)
            j = iwk(1*pil+l)
            k = iwk(2*pil+l)
            moveit = iwk(3*pil+l)
            dx0 = listx(l) - x(i)
            dy0 = listy(l) - y(j)
            dz0 = listz(l) - z(k)
         end if

         iout = i
         dxout = dx0

         dy02 = dy0**2
         dy03 = dy0**3
         dz02 = dz0**2
         dz03 = dz0**3
         dy20 = dy2*dy0
         dyy02= dy*dy02
         dz20 = dz2*dz0
         dzz02= dz*dz02

C

         if (moveit.eq.1) then
            do lk = 1 , order
               do lj = 1 , order
   
                  ierr = reset_offsets4(nx,ny,nz,x,y,z,listx(l),
     .                   charm,dxout,iout,j+lj-1,k+lk-1,i,dx0,order)
                  if (ierr.eq.1) goto 1000
   
                  dx02 = dx0**2
                  dx03 = dx0**3
                  dx20 = dx2*dx0
                  dxx02= dx*dx02
   
                  vx(lj,lk) =  den*(
     .            u(i+3,j+lj-1,k+lk-1)*(dx03+2.0e0*dx20-
     .                                  3.0e0*dxx02)+
     .            u(i+2,j+lj-1,k+lk-1)*(-3.0e0*dx03-9.0e0*dx20
     .                                  +12.0e0*dxx02) +
     .            u(i+1,j+lj-1,k+lk-1)*(3.0e0*dx03+18.0e0*dx20-
     .                                  15.0e0*dxx02)+
     .            u(i  ,j+lj-1,k+lk-1)*(-dx03+6.0e0*dx3-11.0e0*dx20+
     .                                  6.0e0*dxx02) )

               end do
            end do
         else

         do lk = 1 , order
            do lj = 1 , order

                  dx02 = dx0**2
                  dx03 = dx0**3
                  dx20 = dx2*dx0
                  dxx02= dx*dx02

                  vx(lj,lk) =  den*(
     .               u(i+3,j+lj-1,k+lk-1)*(dx03+2.0e0*dx20-
     .                                     3.0e0*dxx02)+
     .               u(i+2,j+lj-1,k+lk-1)*(-3.0e0*dx03-9.0e0*dx20
     .                                     +12.0e0*dxx02) +
     .               u(i+1,j+lj-1,k+lk-1)*(3.0e0*dx03+18.0e0*dx20-
     .                                     15.0e0*dxx02)+
     .               u(i  ,j+lj-1,k+lk-1)*(-dx03+6.0e0*dx3-11.0e0*dx20+
     .                                     6.0e0*dxx02) )

               end do
            end do

         end if

         den = 1.0e0/(6.0e0*dy3)

         do lk = 1 , 4

            vx2(lk) = den*(
     .            vx(4,lk)*(dy03+2.0e0*dy20-3.0e0*dyy02)+
     .            vx(3,lk)*(-3.0e0*dy03-9.0e0*dy20+12.0e0*dyy02)+
     .            vx(2,lk)*(3.0e0*dy03+18.0e0*dy20-15.0e0*dyy02)+
     .            vx(1,lk)*(-dy03+6.0e0*dy3-11.0e0*dy20 +
     .                      6.0e0*dyy02) )
         end do

         den = 1.0e0/(6.0e0*dz3)

         soln(l) = den*(
     .            vx2(4)*(dz03+2.0e0*dz20-3.0e0*dzz02) +
     .            vx2(3)*(-3.0e0*dz03-9.0e0*dz20+12.0e0*dzz02) +
     .            vx2(2)*(3.0e0*dz03+18.0e0*dz20-15.0e0*dzz02) +
     .            vx2(1)*(-dz03+6.0e0*dz3-11.0e0*dz20 +
     .                    6.0e0*dzz02) )

      end do
      didit = didit + 1

      return

1000   continue

      d4lintxyzc   = 1
      return
      end
C------------------------------------------------------------------------------
      integer function getijknc4(n,v,lv,l,del,order)
      implicit none
      integer n,l
      real*8 v(n)
      real*8 del,lv
      integer order,length

      integer i,ip1,im1

      getijknc4= 0

      do i = 1 , n
         if (v(i).ge.lv) goto 100
      end do

      getijknc4= 1
      return

100   continue
      length = nint(0.5*order)

      ip1 = i + 1
      im1 = i - 1

      if (ip1.gt.n.or.im1.lt.1) then
         getijknc4 = 1
         return
      end if

      if (abs(v(i-1)-lv).lt.abs(v(i)-lv)) then
            i = i - 1
      else if (abs(v(i+1)-lv).lt.abs(v(i)-lv)) then
            i = i + 1
      end if
      
C now push this off 2 points ...

      l = i - length 

      if (l.lt.1) l = 1
      if (l+length+1.gt.n) l = n - length - 1 

      del = lv - v(l)

      if (del.lt.(0.0e0)) getijknc4= 1

      return
      end
C======================================================================
      integer function reset_offsets4(nx,ny,nz,x,y,z,listx,
     .                               charm,dx,i,j,k,iout,dxout,order) 
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      real*8 listx
      real*8 charm(nx,ny,nz)
      real*8 dx
      integer i,j,k
      integer order
      integer im,ip,ims,ips
      integer istart,istop
      integer iout
      real*8    dxout


C local variables .......

      integer li,lj,lk
      logical xhit
      integer mini
      real*8  minx

      reset_offsets4= 0

      iout = i
      dxout = dx


C now scan to see if we hit into charm=1
      xhit = .false.
      
      do li = 1 , order 
         if(int(charm(i+li-1,j,k)).ne.1) xhit = .true.
      end do


      if (.not.xhit) return 
C move the i value back to the center of the stencil ........

      minx = 9.90d99
      do li = max(1,i-order),min(nx,i+order)
         if (abs(x(li)-listx).lt.minx) then
            iout = li
            minx = abs(x(li)-listx)   
         end if
      end do

C so the questions are : which way do I move the start point
C we try both directions 

      ip = 0
      im = 0
      ims = 0
      ips = 0
      
      istart= 1
      istop = nx 

C first move it in the negative direction
      do li = iout,istart, -1
         if (int(charm(li,j,k)).eq.1) then
            ims = li
            goto 100
         end if
         im = im + 1
      end do
100   continue

      do li = iout,istop
         if (int(charm(li,j,k)).eq.1) then
            ips = li
            goto 200
         end if
         ip = ip + 1
      end do
200   continue

      if (ips.eq.0.and.ims.eq.0)  then
         reset_offsets4= 1   
         return
      end if

      if (im.lt.ip) then
         iout = ims - nint(0.5*order) + 1
         if (iout.lt.1) then
            iout = 1
            reset_offsets4= 1   
         end if
      else
         iout = ips
         if (iout+order.gt.nx) then
            reset_offsets4= 1   
         end if
      end if

C one more check ....
5000  continue
      ip = 0
      do li = iout,iout+order-1
         ip = ip + int(charm(li,j,k))
      end do
      if (ip.lt.order) then
         iout = max(1,iout-(order-ip))
         goto 5000
      end if


      dxout = listx - x(iout)

      return
      end
C----------------------------------------------------------------------------
      integer function d4lintxyzcd(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln,charm,iwk,swk,didit,wd) 
CFPP$ EXPAND(getsten_num)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 soln(pil)
      real*8 u(nx,ny,nz),charm(nx,ny,nz)
      integer order
      parameter ( order = 4)
      character*1 wd

      integer d4lintxyzc,d4lintzyxc,d4lintyzxc
      integer l,m,lm,ip
      integer maxs,ierr
      integer iwk(*),swk(*)
      integer didit
      parameter ( maxs = 257 * 257 * 5)
      real*8 inlx(maxs),inly(maxs),inlz(maxs),ins(maxs)
      real*8 dx,dy,dz
      real*8 dxi,dyi,dzi
      integer sten_num(maxs),getsten_num 
      real*8 twel
      parameter ( twel = 1.0e0 / 12.0e0)
      real*8 coeff(9,5)

      d4lintxyzcd = 0
      if (pil.gt.257*257) then
         d4lintxyzcd = 1
         return
      end if
      if (wd.ne.' ') then

         call getcoeff4(coeff,2*order+1,order+1)

      end if

      if (wd.eq.' ') then

         d4lintxyzcd = d4lintxyzc(nx,ny,nz,x,y,z,pil,listx,listy,
     .                 listz,u,soln,charm,iwk,didit) 
         if (d4lintxyzcd.eq.1) return
      else if (wd.eq.'x') then
         dx = x(2) - x(1) 
         dxi= 1.0e0/dx 
         ip = 0
         do m = 1 , pil

            if (didit.eq.0) then
            sten_num(m)=getsten_num(nx,ny,nz,x,y,z,listx(m),
     .                      listy(m),listz(m),charm,order,1)

               swk(m) = sten_num(m)
            else
               sten_num(m) = swk(m)
            end if

            if (sten_num(m).eq.3) then
               l = -2
               do lm = 1,2
                  ip = ip + 1
                  inlx(ip) = listx(m) + l*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
               l = 1
               do lm = 3,4
                  ip = ip + 1
                  inlx(ip) = listx(m) + l*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
                  l = l + 1
               end do

            else if (sten_num(m).eq.1) then
               l = 0
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) + l*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else if (sten_num(m).eq.2) then
               l = -1 
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) + l*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else if (sten_num(m).eq.4) then
               l = -3
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) + l*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else if (sten_num(m).eq.5) then
               l = -4
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) + l*dx
                  inly(ip) = listy(m)
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else 
               d4lintxyzcd = 1
               return
            end if
            
         end do

         d4lintxyzcd = d4lintxyzc(nx,ny,nz,x,y,z,ip,
     .               inlx,inly,inlz,u,ins,charm,iwk,didit)
         if (d4lintxyzcd.eq.1) return


         ip = 0

         do m = 1 , pil
            ip = ip + 5
            if (sten_num(m).eq.3) then
               ip = ip - 1
               soln(m) =dxi*( coeff(3,3)*ins(ip-3) +
     .                        coeff(4,3)*ins(ip-2) +
     .                        coeff(6,3)*ins(ip-1) +
     .                        coeff(7,3)*ins(ip  ) )

            else if (sten_num(m).eq.1) then
               soln(m) =dxi*( coeff(5,1)*ins(ip-4) +
     .                        coeff(6,1)*ins(ip-3) +
     .                        coeff(7,1)*ins(ip-2) +
     .                        coeff(8,1)*ins(ip-1) +
     .                        coeff(9,1)*ins(ip  ) )

            else if (sten_num(m).eq.2) then
               soln(m) =dxi*( coeff(4,2)*ins(ip-4) +
     .                        coeff(5,2)*ins(ip-3) +
     .                        coeff(6,2)*ins(ip-2) +
     .                        coeff(7,2)*ins(ip-1) +
     .                        coeff(8,2)*ins(ip  ) )

            else if (sten_num(m).eq.4) then
               soln(m) =dxi*( coeff(2,4)*ins(ip-4) +
     .                        coeff(3,4)*ins(ip-3) +
     .                        coeff(4,4)*ins(ip-2) +
     .                        coeff(5,4)*ins(ip-1) +
     .                        coeff(6,4)*ins(ip  ) )

            else if (sten_num(m).eq.5) then
               soln(m) =dxi*( coeff(1,5)*ins(ip-4) +
     .                        coeff(2,5)*ins(ip-3) +
     .                        coeff(3,5)*ins(ip-2) +
     .                        coeff(4,5)*ins(ip-1) +
     .                        coeff(5,5)*ins(ip  ) )

            end if

         end do

      else if (wd.eq.'y') then
         dy = y(2) - y(1)
         dyi= 1.0e0/dy

         ip = 0
         do m = 1 , pil
            if (didit.eq.0) then
            sten_num(m)=getsten_num(nx,ny,nz,x,y,z,listx(m),
     .                      listy(m),listz(m),charm,order,2)

               swk(m) = sten_num(m)
            else
               sten_num(m) = swk(m) 
            end if

            if (sten_num(m).eq.3) then
               l = -2
               do lm = 1,2
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + l*dy
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
               l = 1
               do lm = 3,4
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + l*dy
                  inlz(ip) = listz(m)
                  l = l + 1
               end do

            else if (sten_num(m).eq.1) then
               l = 0
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + l*dy
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else if (sten_num(m).eq.2) then
               l = -1 
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + l*dy
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else if (sten_num(m).eq.4) then
               l = -3
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + l*dy
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else if (sten_num(m).eq.5) then
               l = -4
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) + l*dy
                  inlz(ip) = listz(m)
                  l = l + 1
               end do
            else 
               d4lintxyzcd = 1
               return
            end if
            
         end do
   
         d4lintxyzcd = d4lintxyzc(nx,ny,nz,x,y,z,ip,
     .             inlx,inly,inlz,u,ins,charm,iwk,didit)
         if (d4lintxyzcd.eq.1) return


         ip = 0

         do m = 1 , pil
            ip = ip + 5
            if (sten_num(m).eq.3) then
               ip = ip - 1
               soln(m) =dyi*( coeff(3,sten_num(m))*ins(ip-3) +
     .                        coeff(4,sten_num(m))*ins(ip-2) +
     .                        coeff(6,sten_num(m))*ins(ip-1) +
     .                        coeff(7,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.1) then
               soln(m) =dyi*( coeff(5,sten_num(m))*ins(ip-4) +
     .                        coeff(6,sten_num(m))*ins(ip-3) +
     .                        coeff(7,sten_num(m))*ins(ip-2) +
     .                        coeff(8,sten_num(m))*ins(ip-1) +
     .                        coeff(9,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.2) then
               soln(m) =dyi*( coeff(4,sten_num(m))*ins(ip-4) +
     .                        coeff(5,sten_num(m))*ins(ip-3) +
     .                        coeff(6,sten_num(m))*ins(ip-2) +
     .                        coeff(7,sten_num(m))*ins(ip-1) +
     .                        coeff(8,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.4) then
               soln(m) =dyi*( coeff(2,sten_num(m))*ins(ip-4) +
     .                        coeff(3,sten_num(m))*ins(ip-3) +
     .                        coeff(4,sten_num(m))*ins(ip-2) +
     .                        coeff(5,sten_num(m))*ins(ip-1) +
     .                        coeff(6,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.5) then
               soln(m) =dyi*( coeff(1,sten_num(m))*ins(ip-4) +
     .                        coeff(2,sten_num(m))*ins(ip-3) +
     .                        coeff(3,sten_num(m))*ins(ip-2) +
     .                        coeff(4,sten_num(m))*ins(ip-1) +
     .                        coeff(5,sten_num(m))*ins(ip  ) )

            end if

         end do

      else if (wd.eq.'z') then
         dz = z(2) - z(1)
         dzi= 1.0e0/dz
         ip = 0
         do m = 1 , pil

            if (didit.eq.0) then 
            sten_num(m)=getsten_num(nx,ny,nz,x,y,z,listx(m),
     .                      listy(m),listz(m),charm,order,3)

               swk(m) = sten_num(m)
            else
               sten_num(m) = swk(m)
            end if


            if (sten_num(m).eq.3) then
               l = -2
               do lm = 1,2
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) 
                  inlz(ip) = listz(m) + l*dz
                  l = l + 1
               end do
               l = 1
               do lm = 3,4
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) 
                  inlz(ip) = listz(m) + l*dz
                  l = l + 1
               end do

            else if (sten_num(m).eq.1) then
               l = 0
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) 
                  inlz(ip) = listz(m) + l*dz
                  l = l + 1
               end do
            else if (sten_num(m).eq.2) then
               l = -1 
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) 
                  inlz(ip) = listz(m) + l*dz
                  l = l + 1
               end do
            else if (sten_num(m).eq.4) then
               l = -3
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) 
                  inlz(ip) = listz(m) + l*dz
                  l = l + 1
               end do
            else if (sten_num(m).eq.5) then
               l = -4
               do lm = 1,order+1
                  ip = ip + 1
                  inlx(ip) = listx(m) 
                  inly(ip) = listy(m) 
                  inlz(ip) = listz(m) + l*dz
                  l = l + 1
               end do
            else 
               d4lintxyzcd = 1
               return
            end if
            
         end do
   
         d4lintxyzcd = d4lintxyzc(nx,ny,nz,x,y,z,ip,
     .            inlx,inly,inlz,u,ins,charm,iwk,didit) 
         if (d4lintxyzcd.eq.1) return


         ip = 0

         do m = 1 , pil
            ip = ip + 5
            if (sten_num(m).eq.3) then
               ip = ip - 1
               soln(m) =dzi*( coeff(3,sten_num(m))*ins(ip-3) +
     .                        coeff(4,sten_num(m))*ins(ip-2) +
     .                        coeff(6,sten_num(m))*ins(ip-1) +
     .                        coeff(7,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.1) then
               soln(m) =dzi*( coeff(5,sten_num(m))*ins(ip-4) +
     .                        coeff(6,sten_num(m))*ins(ip-3) +
     .                        coeff(7,sten_num(m))*ins(ip-2) +
     .                        coeff(8,sten_num(m))*ins(ip-1) +
     .                        coeff(9,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.2) then
               soln(m) =dzi*( coeff(4,sten_num(m))*ins(ip-4) +
     .                        coeff(5,sten_num(m))*ins(ip-3) +
     .                        coeff(6,sten_num(m))*ins(ip-2) +
     .                        coeff(7,sten_num(m))*ins(ip-1) +
     .                        coeff(8,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.4) then
               soln(m) =dzi*( coeff(2,sten_num(m))*ins(ip-4) +
     .                        coeff(3,sten_num(m))*ins(ip-3) +
     .                        coeff(4,sten_num(m))*ins(ip-2) +
     .                        coeff(5,sten_num(m))*ins(ip-1) +
     .                        coeff(6,sten_num(m))*ins(ip  ) )

            else if (sten_num(m).eq.5) then
               soln(m) =dzi*( coeff(1,sten_num(m))*ins(ip-4) +
     .                        coeff(2,sten_num(m))*ins(ip-3) +
     .                        coeff(3,sten_num(m))*ins(ip-2) +
     .                        coeff(4,sten_num(m))*ins(ip-1) +
     .                        coeff(5,sten_num(m))*ins(ip  ) )

            end if

         end do




      else

         write(6,*) '*** This derivative is not implemented ***'
         d4lintxyzcd = 1 
      end if

      return
      end

      integer function movepoint(nx,ny,nz,x,y,z,listx,listy,
     .                     listz,i,j,k,dx0,dy0,dz0,charm,order)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      real*8 charm(nx,ny,nz)
      real*8 listx,listy,listz
      integer i,j,k
      real*8 dx0,dy0,dz0
      real*8 count
      integer order
      integer li,lj,lk,l
      integer is,js,ks
      integer dirp1,dirp2,dirp3
      integer dirm1,dirm2,dirm3
      integer dir1,dir2,dir3
      integer sdir1,sdir2,sdir3
      integer ip
      real*8 dist,xip

      movepoint =  0

      count = 0.0e0

      if (i+order-1.gt.nx.or.j+order-1.gt.ny.
     .    or.k+order-1.gt.nz) then
         movepoint = 1
         return
      end if

      do li = i , i+order -1
         do lj = j , j+order-1
            do lk = k , k+order-1
               count = count + charm(li,lj,lk)
            end do
         end do
      end do

      if (int(count).eq.order**3) return

C put the point to the center of the stencil
      i = i + nint(order*0.5)
      j = j + nint(order*0.5)
      k = k + nint(order*0.5)

      dir1 = -1
      dir2 = -1
      dir3 = -1


      do l = 1 , order*4+1
         do lk = 0 , l
            do lj = 0 , l
               do li = 0 , l
                  count = 0.0e0
                  xip = 0.0e0
                  if (k+lk*dir3+order-1.gt.nz.or.
     .                j+lj*dir2+order-1.gt.ny.or.
     .                i+li*dir1+order-1.gt.nx) then
                     movepoint = 1
                     return
                  end if
                  do ks = k+lk*dir3,k+lk*dir3+order-1
                     do js = j+lj*dir2,j+lj*dir2+order-1
                        do is = i+li*dir1,i+li*dir1+order-1
                           count = count + charm(is,js,ks)
                           xip = xip + 1.0e0
                        end do
                     end do
                  end do
                  if (int(count).eq.xip) goto 1000
                 end do
            end do
         end do
      end do



      movepoint = 1
      i = i - nint(order*0.5)
      j = j - nint(order*0.5)
      k = k - nint(order*0.5)

      return

1000   continue
      i = i +  li*dir1
      j = j + lj*dir2
      k = k + lk*dir3
      dx0 = listx - x(i)
      dy0 = listy - y(j)
      dz0 = listz - z(k)

      return
      end

C-----------------------------------------------------------------------
      integer function getijknew(nx,ny,nz,listx,listy,listz,x,y,z,
     .                     i,j,k,dx0,dy0,dz0,charm,order)
CFPP$ EXPAND(getijknc4)
      implicit none
      integer nx,ny,nz
      real*8 listx,listy,listz
      real*8 x(nx),y(ny),z(nz)
      integer i,j,k
      real*8 dx0,dy0,dz0
      real*8 charm(nx,ny,nz)
      integer order

      integer li,lj,lk
      integer length
      integer ierr,getijknc4
      real*8 dist,mindist
      real*8 is,js,ks
      integer offset

      getijknew =  getijknc4(nx,x,listx,i,dx0,order)
      if (getijknew.eq.1) return
      getijknew =  getijknc4(ny,y,listy,j,dy0,order)
      if (getijknew.eq.1) return
      getijknew =  getijknc4(nz,z,listz,k,dz0,order)
      if (getijknew.eq.1) return


C make sure its the closest point

      length = nint(0.5*order)
      i = i + length
      j = j + length 
      k = k + length 

      getijknew = 0

      mindist = 9.90e09
      offset = 0
100   continue
      do li = i-offset , offset+i
         do lj = j-offset , offset+j
            do lk = k-offset , offset+k
               dist = (x(li)-listx)**2+(y(lj)-listy)**2+(z(lk)-listz)**2
               if (dist.lt.charm(li,lj,lk)*mindist) then
                  mindist = dist
                  is = li
                  js = lj
                  ks = lk
               end if
            end do
         end do
      end do


      if (mindist.lt.(100.0e0)) then
         i = is - length 
         j = js - length
         k = ks - length 
   
         dx0  = listx - x(i)
         dy0  = listy - y(j)
         dz0  = listz - z(k)
      else if (offset.lt.nx) then
         offset = offset + 1
         goto 100
        else
        getijknew = 1
        write(6,*) 'offset,nx=',offset,nx
      end if

      return
      end
C==========================================================================
      integer function getsten_num(nx,ny,nz,x,y,z,listx,listy,listz,
     .                             charm,order,dir)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      real*8 listx,listy,listz
      real*8 charm(nx,ny,nz)
      integer dir
      integer order

      integer ierr,getijknew
      integer length
      integer ip,im,l,lp,lm
      integer i,j,k
      real*8 dx0,dy0,dz0

      length = nint(0.5*order)

      getsten_num =  getijknew(nx,ny,nz,listx,listy,listz,x,y,z,
     .               i,j,k,dx0,dy0,dz0,charm,order)

      if (getsten_num.eq.1) return
      i = i + length
      j = j + length
      k = k + length

C now I have i,j,k

      if (dir.eq.1) then

         ip = 0
         im = 0

         do l = i+1 , i+order
            ip = ip + int(charm(l,j,k))
         end do

         do l = i-order,i-1
            im = im + int(charm(l,j,k))
         end do


         getsten_num = 1 + length

         lm = 0
         lp = order
         do l = 1 , order+1
            if (ip.ge.lp.and.im.ge.lm) then
               getsten_num = l
            end if
            lm = lm + 1
            lp = lp - 1
         end do

         if (ip.eq.im) getsten_num = 1 + length


      else if (dir.eq.2) then

         ip = 0
         im = 0

         do l = j+1 , j+order
            ip = ip + int(charm(i,l,k))
         end do

         do l = j-order,j-1
            im = im + int(charm(i,l,k))
         end do


         getsten_num = 1 + length

         lm = 0
         lp = order
         do l = 1 , order+1
            if (ip.ge.lp.and.im.ge.lm) then
               getsten_num = l
            end if
            lm = lm + 1
            lp = lp - 1
         end do

         if (ip.eq.im) getsten_num = 1 + length

      else

         ip = 0
         im = 0

         do l = k+1 , k+order
            ip = ip + int(charm(i,j,l))
         end do

         do l = k-order,k-1
            im = im + int(charm(i,j,l))
         end do


         getsten_num = 1 + length

         lm = 0
         lp = order
         do l = 1 , order+1
            if (ip.ge.lp.and.im.ge.lm) then
               getsten_num = l
            end if
            lm = lm + 1
            lp = lp - 1
         end do

         if (ip.eq.im.and.ip.gt.0) getsten_num = 1 + length

      end if

1000   continue

      return
      end


      subroutine getcoeff4(co,nc,ns) 
      implicit none
      integer nc,ns
      real*8 co(nc,ns)

      integer l

      co(1,3) = 0.0e0
      co(2,3) = 0.0e0
      co(3,3) = 1.0e0/12.0e0 
      co(4,3) = -2.0e0/3.0e0 
      co(5,3) = 0.0e0
      co(6,3) =  2.0e0/3.0e0 
      co(7,3) = -1.0e0/12.0e0 
      co(8,3) = 0.0e0
      co(9,3) = 0.0e0


      co(1,1) = 0.0e0
      co(2,1) = 0.0e0
      co(3,1) = 0.0e0 
      co(4,1) = 0.0e0 
      co(5,1) = -25.0e0/12.0e0 
      co(6,1) =  4.0e0
      co(7,1) = -3.0e0
      co(8,1) = 4.0e0/3.0e0
      co(9,1) = -1.0e0/4.0e0

      co(1,2) = 0.0e0
      co(2,2) = 0.0e0
      co(3,2) = 0.0e0
      co(4,2) = -0.25e0
      co(5,2) = -5.0e0/6.0e0
      co(6,2) = 1.5e0
      co(7,2) = -0.5e0
      co(8,2) = 1.0e0/12.0e0
      co(9,2) = 0.0e0

      co(1,4) = 0.0e0
      co(2,4) = -1.0e0/12.0e0 
      co(3,4) = 0.5e0 
      co(4,4) = -1.5e0 
      co(5,4) = 5.0e0/6.0e0 
      co(6,4) = 0.25e0 
      co(7,4) = 0.0e0 
      co(8,4) = 0.0e0 
      co(9,4) = 0.0e0

      co(1,5) = 0.25e0 
      co(2,5) = -4.0e0/3.0e0 
      co(3,5) = 3.0e0 
      co(4,5) = -4.0e0 
      co(5,5) = 25.0e0/12.0e0 
      co(6,5) = 0.0e0
      co(7,5) = 0.0e0
      co(8,5) = 0.0e0
      co(9,5) = 0.0e0

      return
      end
C==============================================================================
      subroutine getcoeff2(co,nc,ns) 
      implicit none
      integer nc,ns
      real*8 co(nc,ns)

      integer l

      co(1,2) = 0.0e0
      co(2,2) =-0.5e0
      co(3,2) = 0.0e0
      co(4,2) = 0.5e0
      co(5,2) = 0.0e0

      co(1,1) = 0.0e0
      co(2,1) = 0.0e0
      co(3,1) =-1.5e0
      co(4,1) = 2.0e0
      co(5,1) =-0.5e0

      co(1,3) = 0.5e0
      co(2,3) =-2.0e0
      co(3,3) = 1.5e0
      co(4,3) = 0.0e0
      co(5,3) = 0.0e0

      return
      end
C==============================================================================
C This routine takes a list of points and interpolates to this to nth
C order
      integer function dnlintxyzc(nx,ny,nz,x,y,z,pil,listx,listy,listz,
     .                           u,soln,order,charm,iwk,stwk,didit)
CFPP$ EXPAND(getijknew,reset_offsets4,movepoint)
      implicit none
      integer nx,ny,nz
      real*8 x(nx),y(ny),z(nz)
      integer pil
      real*8 listx(pil),listy(pil),listz(pil)
      real*8 u(nx,ny,nz),charm(nx,ny,nz)
      real*8 soln(pil)

      integer order
      integer didit
      integer iwk(4*pil),movepoint,moveit
      real*8  stwk((order*order*order+2*order)*pil)
      integer istwk

C local variables
      real*8 dx,dy,dz
      integer l,lj,lk,li,lli
      integer i,j,k
      real*8 dx0,dy0,dz0
      integer ierr, getijknew
      integer maxorder
      parameter ( maxorder = 64 )
      real*8 mat2d(maxorder,maxorder)
      real*8 vx(maxorder,maxorder),vx2(maxorder)
      real*8 den,work(maxorder),sten(maxorder)
      integer getcoc,reset_offsets4
      integer iw,iout,isoff
      real*8 dxout

      real*8 dy2,dy3,dy4,dy5
      real*8 dz2,dz3,dz4,dz5

      dnlintxyzc= 0

      dx = abs( x(2)-x(1) )
      dy = abs( y(2)-y(1) )
      dz = abs( z(2)-z(1) )
      iw = nint(0.5*order)

      isoff = order*order*pil
      istwk = 0

      do l = 1 , pil

         if (didit.eq.0) then
            iwk(3*pil+l) = 1
            moveit = 1
            dnlintxyzc =
     .      getijknew(nx,ny,nz,listx(l),listy(l),listz(l),
     .                x,y,z,i,j,k,dx0,dy0,dz0,charm,order)


            if (dnlintxyzc.eq.1) return
            if (charm(i+iw,j+iw,k+iw).eq.1) then
               ierr = movepoint(nx,ny,nz,x,y,z,listx(l),listy(l),
     .                     listz(l),i,j,k,dx0,dy0,dz0,charm,order)

               if (ierr.eq.1) goto 1000
               moveit = 0
               iwk(3*pil+l) = 0
            end if
            iwk(l) = i
            iwk(1*pil+l) = j
            iwk(2*pil+l) = k
         else
            i = iwk(l)
            j = iwk(1*pil+l)
            k = iwk(2*pil+l)
            moveit = iwk(3*pil+l)
            dx0 = listx(l) - x(i)
            dy0 = listy(l) - y(j)
            dz0 = listz(l) - z(k)
         end if
         iout = i
         dxout = dx0

C
         do lk = 1 , order
            do lj = 1 , order
               vx(lj,lk) = 0.0e0
            end do
               vx2(lk) = 0.0e0
         end do

         if (moveit.eq.1) then

            do lk = 1 , order
               do lj = 1 , order
                  ierr = reset_offsets4(nx,ny,nz,x,y,z,
     .             listx(l),charm,dxout,iout,j+lj-1,k+lk-1,i,dx0,order)
                  if (ierr.eq.1) goto 1000

                  if (didit.eq.0) then
                     ierr = getcoc(mat2d,order,dx,dx0,sten)
                     if (ierr.eq.1) goto 1000
                     do lli = 1 , order
                        istwk = istwk + 1 
                        stwk(istwk) = sten(lli) 
                     end do
                  else
                     do lli = 1 , order
                        istwk = istwk + 1
                        sten(lli)=stwk(istwk) 
                     end do
                  end if

                  if (ierr.eq.1) goto 1000
                  do li = 1 , order
                     vx(lj,lk)=sten(li)*u(i+li-1,j+lj-1,k+lk-1)+
     .                         vx(lj,lk)
                  end do
               end do
            end do

         else

            do lk = 1 , order
               do lj = 1 , order

                  if (didit.eq.0) then
                     ierr = getcoc(mat2d,order,dx,dx0,sten)
                     if (ierr.eq.1) goto 1000
                     do lli = 1 , order
                        istwk = istwk + 1
                        stwk(istwk) = sten(lli) 
                     end do
                  else
                     do lli = 1 , order
                        istwk = istwk + 1
                        sten(lli)=stwk(istwk) 
                     end do
                  end if

                  do li = 1 , order
                     vx(lj,lk)=sten(li)*u(i+li-1,j+lj-1,k+lk-1)+
     .                         vx(lj,lk)
                  end do
               end do
            end do
         end if

         if (didit.eq.0) then
            ierr = getcoc(mat2d,order,dy,dy0,sten)
            if (ierr.eq.1) goto 1000
            do lli = 1 , order
               istwk = istwk + 1
               stwk(istwk) = sten(lli) 
            end do
         else
            do lli = 1 , order
               istwk = istwk + 1
               sten(lli) = stwk(istwk) 
            end do
         end if

         do lk = 1 , order
            do li = 1 , order
               vx2(lk) = sten(li)*vx(li,lk) + vx2(lk)
            end do
         end do

         if (didit.eq.0) then
            ierr = getcoc(mat2d,order,dz,dz0,sten)
            if (ierr.eq.1) goto 1000
            do lli = 1 , order
               istwk = istwk + 1
               stwk(istwk) = sten(lli) 
            end do
         else
            do lli = 1 , order
               istwk = istwk + 1
               sten(lli) = stwk(istwk) 
            end do
         end if

         if (ierr.eq.1) goto 1000
         soln(l) = 0.0e0

         do li = 1 , order
            soln(l) = soln(l) + sten(li)*vx2(li)
         end do

      end do

      didit = didit + 1
      return

1000   continue

      dnlintxyzc   = 1
      return
      end
C------------------------------------------------------------------------------
      integer function getcoc(mat2d,order,dx,dx0,sten)
      implicit none
      integer order
      real*8 dx,dx0
      real*8 mat2d(order,order)
      real*8 sten(order)
      integer lj,lk,ljs
      real*8 factorial

      integer maxorder
      parameter ( maxorder = 64 )
      integer ipvt(maxorder)
      real*8 work(maxorder)
      real*8 det
      integer ifail

      ljs = 0
      do lj = 1 , order
         if ( abs((lj-1)*dx-dx0).lt.(1.0e-10)) then
            ljs = lj
         end if
      end do

      if (ljs.eq.0) then
         do lk = 1 , order
            do lj = 1 , order
               mat2d(lk,lj)=((lk-1)*dx-dx0)**(lj-1)/factorial(lj-1)
            end do
         end do

         call sgetrf(mat2d,order,order,ipvt,ifail)
         call sgetri(mat2d,order,order,ipvt,det,work,01)

         do lj = 1 , order
            sten(lj) = mat2d(1,lj)
         end do
      else
         do lj = 1 , order
            sten(lj) = 0.0e0
         end do

         sten(ljs) = 1.0e0
      end if
      getcoc = ifail

      return
      end

C======================================================================
      integer function deriv(n,f,dir,del,sout,mat) 
      implicit none
      integer n
      real*8 f(n)
      real*8 sout
      real*8 del,mat(n,n),ipvt(40)
      real*8 factorial
      real*8 fact
      integer dir
      integer i,j,k,l,ifail

      if (dir.ne.0) then
         do k = 1 , n
            mat(k,1) = 1.0e0
            do j = 2 , n
               mat(k,j)=((k-1)*del*dir)**(j-1)/factorial(j-1)
            end do
         end do
      else
         k = 0
         do l = -int(n*0.4999),int(n*0.5) 
            k = k + 1
            mat(k,1) = 1.0e0
            do j = 2 , n
               mat(k,j)=(l*del)**(j-1)/factorial(j-1)
            end do
         end do
      end if

      call sgetrf(mat,n,n,ipvt,ifail)
      deriv = ifail
      if (deriv.ne.0) return
      call sgetrs(mat,n,n,ipvt,f,0)

      sout = f(2)

      return
      end


C-------------------------------------------------------------------------------C
C  General 3 dimensional interpolation subroutine base on matt's 1D
C    interpolation routine DVINQN.
C
C---------------------------------------------------------------------
      subroutine d3vinqn(a3,a3tmp1,a3tmp2,a3bar,x,y,z,
     *                   xbar,ybar,zbar,nx,ny,nz,
     *                   nxbar,nybar,nzbar,
     *                   vout,nintrp)

      implicit none

      integer nx,ny,nz,nxbar,nybar,nzbar,nintrp, maxdim, i,j,k
      parameter (maxdim = 1025)
      real*8 a3(nx,ny,nz), a3tmp1(nxbar,ny,nz),a3tmp2(nxbar,nybar,nz)
      real*8 a3bar(nxbar,nybar,nzbar), vout

      real*8 a(maxdim), x(*), y(*), z(*)
      real*8 abar(maxdim), xbar(*), ybar(*), zbar(*)

      if (nx.gt.maxdim.or.ny.gt.maxdim.or.nz.gt.maxdim.or.
     *    nxbar.gt.maxdim.or.nybar.gt.maxdim.or.nzbar.gt.maxdim) then
         write(0,*)'ERROR! d2vinqn: maxdim=',maxdim,' too small'
         goto 900
      end if

      do k = 1,nz
         do j = 1,ny
            call d3vputv(a3,a,nx,ny,nz,j,k,1,1)
            call dvinqn(a,x,abar,xbar,nx,nxbar,vout,vout,nintrp)
            call dvput3v(abar,a3tmp1,nxbar,ny,nz,j,k,1)
         end do
      end do

      do k = 1,nz
         do j = 1,ny
            call d3vputv(a3,a,nx,ny,nz,j,k,1,1)
            call dvinqn(a,x,abar,xbar,nx,nxbar,vout,vout,nintrp)
            call dvput3v(abar,a3tmp1,nxbar,ny,nz,j,k,1)
         end do
      end do

      do k = 1,nz
         do i=1,nxbar
            call d3vputv(a3tmp1,a,nxbar,ny,nz,i,k,1,2)
            call dvinqn(a,y,abar,ybar,ny,nybar,vout,vout,nintrp)
            call dvput3v(abar,a3tmp2,nxbar,nybar,nz,i,k,2)
         end do
      end do

      do i=1,nxbar
         do j = 1,nybar
            call d3vputv(a3tmp2,a,nxbar,nybar,nz,i,j,1,3)
            call dvinqn(a,z,abar,zbar,nz,nzbar,vout,vout,nintrp)
            call dvput3v(abar,a3bar,nxbar,nybar,nzbar,i,j,3)
         end do
      end do

      return
 900  continue
      end

C-------------------------------------------------------------------------------C
C  subroutine that inserts part of a 3D vector into a 1D vector.
C
C-----------------------------------------------------------------------------
      subroutine d3vputv(a3,a,nx,ny,nz,l,m,stride,case)

      implicit none

      integer nx,ny,nz,l,m,case,n
      integer stride

      real*8 a3(nx,ny,nz), a(1025)
      integer i


      i = 0

      if (case.eq.1) then
         do n = 1,nx,stride
            i = i + 1
            a(i) = a3(n,l,m)
         end do
      else if (case.eq.2) then
         do n = 1,ny,stride
            i = i + 1
            a(i) = a3(l,n,m)
         end do
      else if (case.eq.3) then
         do n = 1,nz,stride
            i = i + 1
            a(i) = a3(l,m,n)
         end do
      else
         write(0,*)'ERROR: d3vputv - case=',case,' not correct'
         goto 900
      end if

      return
 900  continue
      stop
      end
C-------------------------------------------------------------------------------C
C  subroutine that inserts a 1D vector of a given stride into a 3D vector in
C   a specified place.
C
C-------------------------------------------------------------------------------
      subroutine dvput3v(a,a3,nx,ny,nz,l,m,case)

      implicit none

      integer nx,ny,nz,l,m,case,n

      real*8 a(1025), a3(nx,ny,nz)

      if (case.eq.1) then
         do n = 1,nx
            a3(n,l,m) = a(n)
         end do
      else if (case.eq.2) then
         do n = 1,ny
            a3(l,n,m) = a(n)
         end do
      else if (case.eq.3) then
         do n = 1,nz
            a3(l,m,n) = a(n)
         end do
      else
         write(0,*)'ERROR: dvput3v - case=',case,' not correct'
         goto 900
      end if

      return
 900  continue
      stop
      end

      integer function icharsten3(val,nx,ny,nz,f,order)
      implicit none
      integer nx,ny,nz
      real*8 val
      real*8 f(nx,ny,nz)

      integer i,j,k
      integer order

      icharsten3 = 0

      do  k = 1 , nz
         do j= 1 , ny


            do i = 1 , order 
               if (abs(f(i,j,k)).gt.val) then
                  icharsten3 = 1
                  return
               end if
            end do

            do i = nx,nx-order+1,-1
               if (abs(f(i,j,k)).gt.val) then
                  icharsten3 = 1
                  return
               end if
            end do

         end do
      end do

      return
      end
C----------------------------------------------------------------------------

c-----------------------------------------------------------------------
c
c     2:1 cubic interpolation of VC to VF ... with a char function
c
c-----------------------------------------------------------------------

      subroutine dv2i4_sc(vc,chc,vf,chf,nc)

         implicit      none

         integer       nc
         real*8        vc(nc),     vf(*)
         real*8        chc(nc),    chf(*)

         integer       i,          nf

         real*8        lc(4),      cc(4),      rc(4)
         data
     *     lc / 0.31250e0,  0.93750e0, -0.31250e0,  0.06250e0 /,
     *     cc /-0.06250e0,  0.56250e0,  0.56250e0, -0.06250e0 /,
     *     rc / 0.06250e0, -0.31250e0,  0.93750e0,  0.31250e0 /

         if( nc .ge. 4 ) then
            nf = 2 * nc - 1
C injection:
            call dvprln_sc(vc,vf,2,nc,chc)
C
            vf(2) = lc(1) * vf(1)* chf(1) + lc(2) * vf(3) * chf(3) +
     *              lc(3) * vf(5)* chf(5) + lc(4) * vf(7) * chf(7)
            do 10 i = 4 , nf - 3 , 2
               vf(i)=cc(1)*vf(i-3)*chf(i-3)+cc(2)*vf(i-1)*chf(i-1) +
     *                 cc(3)*vf(i+1)*chf(i+1)+cc(4)*vf(i+3)*chf(i+3)
 10         continue
            vf(nf-1)=rc(1)*vf(nf-6)*chf(nf-6)+rc(2)*vf(nf-4)*chf(nf-4)+
     *               rc(3)*vf(nf-2)*chf(nf-2)+rc(4)*vf(nf)*chf(nf)
         else
            write(*,*) '>>> dv2i4_s:: Too few points for cubic '//
     *                 'interpolation ...'

         end if

         return


      end
c----------------------------------------------------------------------
c
c     V5 :=   Sum  [ S(i) * Vi ] ... used by multidimensional 2:1 cubic
c           i=1..4                   interpolation ...
c
c-----------------------------------------------------------------------

      subroutine dvsma4_sc(v1,v2,v3,v4,c1,c2,c3,c4,v5,s,n)

         implicit       none

         integer        n
         real*8         v1(n),       v2(n),      v3(n),      v4(n),
     *                  v5(n)
         real*8         c1(n),       c2(n),      c3(n),      c4(n)
         real*8         s(4)

         integer        j

         do 10 j = 1 , n
            v5(j) = s(1) * v1(j) *c1(j) + s(2) * v2(j)*c2(j)  +
     *              s(3) * v3(j)*c3(j) + s(4) * v4(j)*c4(j)
 10      continue
         return

      end


c-----------------------------------------------------------------------
      subroutine d1lv(v,n,val)
      implicit none
      integer n
      real*8 v(n)
      real*8 val
      integer i

      do i = 1 , n
         v(i) = val
      end do

      return
      end
C======================================================================
      integer function construct_matrixc(order,mat)
      implicit none
      integer order
      real*8 mat(order,order+1)
      integer maxorder
      parameter ( maxorder = 20 )
      real*8 mat2d(maxorder,maxorder),mat2ds(maxorder,maxorder)
      real*8 work(maxorder)
      integer ipvt(maxorder)
      integer construct_matrixc_

      construct_matrixc= construct_matrixc_(order,mat2d,mat2ds,
     .                                    ipvt,mat,work)

      return
      end
C======================================================================
      integer function construct_matrixc_(order,mat2d,mat2ds,
     .        ipvt,mat,work)
      implicit none

      integer order
C mat is stored as the coefficients and for each stencil
C coef# , stencil#
C 
      real*8 mat(order,order+1)
      real*8 mat2d(order,order),mat2ds(order,order)
      real*8 factorial
      real*8 work(order),det
      integer fa
      integer ipvt(order)
      integer length
      integer ifail
      integer ioff

      integer i,j,k,l
      integer iendp1

      construct_matrixc_ = 0


      length = nint(order*0.5)
      k = 1
      ioff = mod(order,2)

C do the positive direction first ...
      fa = -1
      do i = 1 , length 
         fa = fa + 2
         do j = 1 , order
            mat2d(i,j) = fa**(j-1)/factorial(j-1) 
         end do
      end do

      fa = 1
      iendp1 = i
C nod do the negative direction ....
      do i = i , order
         fa = fa - 2   
         do j = 1 , order
            mat2d(i,j) = fa**(j-1)/factorial(j-1) 
         end do
      end do

      call copy(order*order,mat2d,mat2ds)

C now do a lud decomposition on this and then an iverse and then
C store the coefficients

      call sgetrf(mat2d,order,order,ipvt,ifail)
      call sgetri(mat2d,order,order,ipvt,det,work,01)

C now get out the coefficients .....

      construct_matrixc_ = ifail
      if (ifail.ne.0) return

C now copy into mat array
      call copy_mat2d_matc(order,mat2d,k,mat)

C now lets get the other stencils starting from the second point
C and then going forward ......
C assume 2x1 refinements ...

C first modify all the positive  terms that we must modify ...

      do l = length  - 1  , 1 , -1

C l is the number of terms to modify
         fa = -order + 1 + ioff
         k = k + 1
         call copy(order*order,mat2ds,mat2d)
C i is which term in the matrix ....
         do i = length , length - l + 1  , -1
            fa = fa - 2
            do j = 1 , order
               mat2d(i,j) = fa**(j-1)/factorial(j-1)
            end do
         end do

C now do a lud decomposition on this ........

         call sgetrf(mat2d,order,order,ipvt,ifail)
         call sgetri(mat2d,order,order,ipvt,det,work,01)

         construct_matrixc_ = ifail
         if (ifail.ne.0) return

C now copy into mat array
         call copy_mat2d_matc(order,mat2d,k,mat)

      end do
C now modify the negative terms .....

      do l = length  - 1  , 1 + ioff , -1

C l is the number of terms to modify
C         fa = order - 1
         fa = -1 + length*2
         k = k + 1
         call copy(order*order,mat2ds,mat2d)
C i is which term in the matrix ....
         do i = order , order - l + 1 +ioff  , -1
            fa = fa + 2
            do j = 1 , order
               mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
            end do
         end do

C now do a lud decomposition on this ........
         call sgetrf(mat2d,order,order,ipvt,ifail)
         call sgetri(mat2d,order,order,ipvt,det,work,01)

         construct_matrixc_ = ifail
         if (ifail.ne.0) return

C now copy into mat array
         call copy_mat2d_matc(order,mat2d,k,mat)

      end do

C now get the extrapolation points

      k = k + 1
      fa = 1.0e0
      do i = 1 , order
         do j = 1 , order
            mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
         end do
         fa = fa + 1.0e0
      end do
C now do a lud decomposition on this ........
      call sgetrf(mat2d,order,order,ipvt,ifail)
      call sgetri(mat2d,order,order,ipvt,det,work,01)

      construct_matrixc_ = ifail
      if (ifail.ne.0) return

C now copy into mat array
      call copy_mat2d_matc(order,mat2d,k,mat)

C 
      k = k + 1
      fa = -1.0e0
      do i = 1 , order
         do j = 1 , order
            mat2d(i,j) = (1.0e0*fa)**(j-1)/factorial(j-1)
         end do
         fa = fa - 1.0e0
      end do
C now do a lud decomposition on this ........
      call sgetrf(mat2d,order,order,ipvt,ifail)
      call sgetri(mat2d,order,order,ipvt,det,work,01)

      construct_matrixc_ = ifail
      if (ifail.ne.0) return

C now copy into mat array
      call copy_mat2d_matc(order,mat2d,k,mat)

      return
      end
