      function GetSecChars (ireal, key)
c --- primary function: retireve the integer section data
c
c  input variables:
c   ireal           (int,sc,in)       - section number.
c   key             (int,sc,in)       - which integer value to return
c
c  outputs variables:
c   GetSecChars     (int,sc,out)      - integer data to return
c                                       NumSectPts      (key = 1)
c                                       NumSectRstNodes (key = 2)
#include "impcom.inc"
#include "constants.inc"
c
      integer  SecIdNum, key, GetSecChars, ireal, rlgetn
      integer  NumElems, numIntPts, numCirElems, numThruElems, j
      double precision rvr
      parameter (numCirElems = 8, 
     &           numThruElems = 2, 
     &           numIntPts = 4)

#if defined(DEBUG)
      call SubBeg('GetSecChars')
#endif

      j = rlgetn(ireal,1,1,rvr)

      SecIdNum = abs(nint(rvr))

      if (key.eq.1) then
          numElems = 0
          if (SecIdNum.eq.1) then
              numElems = numThruElems*numThruElems
          elseif (SecIdNum.eq.2) then
              numElems = 12 ! v13
          elseif (SecIdNum.eq.3) then
              numElems = numCirElems*numThruElems
          elseif (SecIdNum.eq.4) then
              numElems = numCirElems
          endif
          GetSecChars = numElems*numIntPts
      elseif (key.eq.2) then
          if (SecIdNum.eq.1) then
              GetSecChars = 9
          elseif (SecIdNum.eq.2) then
              GetSecChars = 24 ! v13
          elseif (SecIdNum.eq.3) then
              GetSecChars = 17
          elseif (SecIdNum.eq.4) then
              GetSecChars = 16
          endif
      else
          GetSecChars = 0
      endif

#if defined(DEBUG)
      call SubEnd('GetSecChars')
#endif

      return
      end
