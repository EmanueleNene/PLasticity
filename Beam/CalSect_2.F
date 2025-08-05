      subroutine CalSect_2(numElems,matID,Summary,MatWt,ErrorFlag)
c --- primary function: calculate section data
c
#include "impcom.inc"
#include "constants.inc"
c
      integer numElems, kel, ErrorFlag, first

      double precision
     &     matID(*), Summary(*),MatWt(*)


      first = 0
      ErrorFlag = 0
      do kel=1,numElems
         if (int(matID(kel)).gt.0 )then
             first=int(matID(kel))
             goto 5
          endif
      enddo
 5    continue
      if (first.eq.0) then
          ErrorFlag = 1
          goto 800
      endif

      Summary(20) = 0.0d0
      do kel=1,numElems
         if (int(matID(kel)).gt.0) then
             MatWt(kel)=1.0d0
             MatWt(kel+numElems)=1.0d0
             MatWt(kel+2*numElems)=1.0d0
         else
             MatWt(kel)=0.0d0
             MatWt(kel+numElems)=0.0d0
             MatWt(kel+2*numElems)=0.0d0
         endif
      enddo

      goto 900
 800  continue
      ErrorFlag = 1
 900  continue


      return
      end


