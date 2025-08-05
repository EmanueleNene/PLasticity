      subroutine v_gauss (y,n,seed)
c --- primary function:    generate a vector of normal numbers
c
c  input arguments:
c     n        (int,sc,in)       - number of normals desired
c     seed     (dp,sc,inout)     - random number generator seed
c
c  output arguments:
c     y        (dp,ar(n),out)    - vector of numbers from gaussian dist
c     seed     (dp,sc,inout)     - random number generator seed
c
#include "impcom.inc"
c
      integer  nval
      parameter  (nval=29)
      integer  n,i,j
      double precision  y(n),seed,xval(nval),yval(nval),x,con
                                                                  
      data  xval /
     & 0.000000d0, 0.023922d0, 0.047758d0, 0.071424d0, 0.094835d0,
     & 0.117911d0, 0.148027d0, 0.198468d0, 0.219043d0, 0.245373d0,
     & 0.270350d0, 0.299546d0, 0.321214d0, 0.346136d0, 0.384930d0,
     & 0.419243d0, 0.449497d0, 0.469946d0, 0.480301d0, 0.489830d0,
     & 0.495060d0, 0.497445d0, 0.499032d0, 0.499517d0, 0.499807d0,
     & 0.499912d0, 0.499952d0, 0.499983d0, 0.500000d0/
      data  yval /
     & 0.0000d0,   0.0600d0,   0.1200d0,   0.1800d0,   0.2400d0,
     & 0.3000d0,   0.3800d0,   0.5200d0,   0.5800d0,   0.6600d0,
     & 0.7400d0,   0.8400d0,   0.9200d0,   1.0200d0,   1.2000d0,
     & 1.4000d0,   1.6400d0,   1.8800d0,   2.0600d0,   2.3200d0,
     & 2.5800d0,   2.8000d0,   3.1000d0,   3.3000d0,   3.5500d0,
     & 3.7500d0,   3.9000d0,   4.1500d0,   5.0000d0/

c --- get random numbers in vector
      call v_rand (y(1),n,seed)
c --- convert to normal numbers
      do 20  i = 1,n
         x = abs(y(i)-0.5d0)
         do 15  j = 2,nval
            if (x .le. xval(j)) then
               con = (x-xval(j-1))/(xval(j)-xval(j-1))
               con = yval(j-1) + con*(yval(j)-yval(j-1))
               if (y(i) .lt. 0.5d0) con = -con
               y(i) = con
               go to 20
             endif
  15     continue
  20  continue

      return
      end
