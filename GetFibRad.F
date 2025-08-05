      function GetFibRad(elem)
c --- primary function: return contact radius inside FibNet
c
#include "impcom.inc"
#include "solu_info.inc"
c
      integer elmiqr, i, rlgetn, elem, ireal, cstype
      
      double precision GetFibRad, rvr(3)

      GetFibRad = 0.0d0

      if (elem.le.0) goto 999
      ireal = elmiqr(elem,-3)

      if (ireal.le.0) goto 999
     
      i = rlgetn(ireal,1,3,rvr(1))
      if (i.lt.2) goto 999
      if (rvr(2).le.0.0d0) goto 999
      cstype = nint(rvr(1))
      if (cstype.le.2) then
c --- rectangular cross section - output height and not the radius
          if (A_SpecVer.ne.2) then
              GetFibRad = 0.25d0*(rvr(2)+rvr(3))
          else
              GetFibRad = 0.5d0*min(rvr(2),rvr(3))
          endif
      else
          GetFibRad = rvr(2)
      endif      

  999 continue      
      return
      end

