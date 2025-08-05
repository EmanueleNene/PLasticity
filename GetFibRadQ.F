      function GetFibRadQ(elem)
c --- primary function: return contact radius outside FibNet
c
#include "impcom.inc"
c
      integer elmiqr, i, rlgetn, elem, ireal, cstype
      
      double precision GetFibRadQ, rvr(3)

      GetFibRadQ = 0.0d0

      if (elem.le.0) goto 999
      ireal = elmiqr(elem,-3)

      if (ireal.le.0) goto 999
     
      i = rlgetn(ireal,1,3,rvr(1))
      if (i.lt.2) goto 999
      if (rvr(2).le.0.0d0) goto 999
      cstype = nint(rvr(1))
      if (cstype.le.2) then
c --- rectangular cross section - output height and not the radius
          GetFibRadQ = 0.25d0*(rvr(2)+rvr(3))
      else
          GetFibRadQ = rvr(2)
      endif      

  999 continue      
      return
      end

