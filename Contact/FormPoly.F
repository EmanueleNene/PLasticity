      function FormPoly (coef,ksi,n)
c --- primary function: get intersection for parabola target
c
#include "impcom.inc"
c
      integer n,i
      double precision FormPoly,coef(n),ksi
c
      FormPoly = coef(n)
      do i = 1,n-1
         FormPoly = FormPoly + coef(n-i)*ksi**i
      enddo
c
      return
      end
