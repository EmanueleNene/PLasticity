      function IfEqvDP (a,b)
c --- primary function: check if two dp numbers are equal (within machine tolerance)
c
c  input arguments:
c     a        (dp,sc,in)          - first dp number
c     b        (dp,sc,in)          - second dp number
c
c  output arguments:
c     IfEqvDP   (log,sc,out)        - .true. if the numbers are equal
c
c
#include "impcom.inc"
c
      double precision   a,b,machep
      logical IfEqvDP
      parameter (machep=1.0d-12)

      IfEqvDP = .false.
      if (abs(a-b).le.machep*abs(b)) IfEqvDP = .true.

      return
      end
