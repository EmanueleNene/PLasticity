      subroutine TransBeam (d12,tr)
c --- primary function: geometric transform for beam elements w or w/o a 3rd node
c
c  input arguments:
c     d12      (dp,ar(3),in)     - shape function derivatives
c
c  output arguments:
c     tr       (dp,ar(3,3),in)   - transformation matrix
c
#include "impcom.inc"
c
      double precision
     & tr(3,3),aleng,d12(3),con,aleng2,alen1,
     & salp1,calp1,salp2,calp2,salp3,calp3

c --- initialization
c --- compute current length and length components

      con = d12(1)*d12(1) + d12(2)*d12(2)
      aleng2 = con + d12(3)*d12(3)
      aleng = sqrt(aleng2)
      if (aleng.eq.0.0d0) aleng = 1.0d0

c --- calculate length in the x-y plane
      alen1 = 0.0d0
      if (con .gt. 0.0d0) alen1 = sqrt(con)

c --- compute sin and cos of third angle
      salp3 = 0.0d0
      calp3 = 1.0d0

c --- calculate sin and cos of first angle
      if (alen1 .le. 0.0001d0*aleng) then
          salp1 = 0.0d0
          calp1 = 1.0d0
      else
          salp1 = d12(2)/alen1
          calp1 = d12(1)/alen1
      endif

c --- calculate sin and cos of second angle
      salp2 = d12(3)/aleng
      calp2 = alen1/aleng

c --- calculate transformation matrix
      tr(1,1) = calp1*calp2
      tr(2,1) = - calp1*salp2*salp3 - salp1*calp3
      tr(3,1) = - calp1*salp2*calp3 + salp1*salp3
      tr(1,2) = salp1*calp2
      tr(2,2) = - salp1*salp2*salp3 + calp1*calp3
      tr(3,2) = - salp1*salp2*calp3 - calp1*salp3
      tr(1,3) = salp2
      tr(2,3) = salp3*calp2
      tr(3,3) = calp3*calp2

      return
      end
