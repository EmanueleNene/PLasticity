      subroutine GetIsoYield (prop, option, pleq, sigy, dsigdep, sigy0)
c --- primary function: calculate material yield stress and plastic slope
c
c  input arguments:
c      elemId   (int,sc,in)               - element # 
c      matId    (int,sc,in)               - material #
c      temp     (dp,sc,in)                - current moisture
c      pleq     (dp,sc,in)                - equivalent plastic strain
c
c  output arguments:
c      sigy     (dp,sc,out)               - yield stress at 't+dt'
c      dsigdep  (dp,sc,out)               - plastic slope at 't+dt'
c
#include "impcom.inc"
c
      double precision ONE
      parameter       (ONE = 1.0d0)

      double precision pleq, sigy, dsigdep, sigy0

      double precision prop(*), bkh, r0, rinf
      
      integer          option

      if (prop(1).gt.0.0d0) then
          sigy = prop(1)
      else
          sigy = sigy0 
      endif
      r0 = prop(2)
      rinf = prop(3)
      bkh = prop(4)
      sigy = sigy + r0*pleq + rinf*(ONE-exp(-bkh*pleq))
      dsigdep = r0 + rinf*bkh*exp(-bkh*pleq)
     
      return
      end

