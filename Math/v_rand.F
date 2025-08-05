      subroutine v_rand (y,n,seed)
c --- primary function: fill a vector with random numbers
c
c  input arguments:
c     y        (dp,ar(n),inout)  - vector y
c     n        (int,sc,in)       - length of vector y
c     seed     (dp,sc,inout)     - control key and/or random seed
c
c     seed codes:
c              0.0d0    - initialize the default sequence
c              1.0d0    - initialize the custom sequence
c
c  output arguments:
c     y        (dp,ar(n),inout)  - vector y
c     seed     (dp,sc,inout)     - restored point
c
#include "impcom.inc"
c
      integer  n,i
      double precision
     &         seed,y(n),am1,am2,am3,am4

c --- initialize series
      data  am1,am2,am3 /0.32007779d0,0.23717810d0,0.52636370d0/

      if (seed .eq. 0.0d0) then
          seed = am1
      endif

c --- calculate random vector
      am1 = seed
      do i = 1,n
         am4 = am1 + am2 + am3
         if (am2 .lt. 0.5d0) am4 = am4 + 0.00001357d0
         if (am4 .ge. 1.0d0) am4 = am4 - 1.0d0
         if (am4 .ge. 1.0d0) am4 = am4 - 1.0d0
c --- truncate to 8 significant digits
         am4 = 1.0d-8*nint(am4*1.0d8)
         am1 = am2
         am2 = am3
         am3 = am4
         y(i) = am3
      enddo
      seed = am1

 999  return
      end
