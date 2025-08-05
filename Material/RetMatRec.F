      subroutine RetMatRec (matRecordIndex, kMatRecPt, nTens, vect)
c --- primary function: retrieve current material record data
c
c  input arguments:
c     matRecordIndex(in,sc,in)        - material record index
c     kMatRecPt     (in,sc,in)        - integration point where data to be inquired
c     nTens         (in,sc,in)        - number of components to be inquired
c     numWordsPt    (in,ar(*),in)     - number of material records
c
c  output arguments:
c     vect          (dp,ar(*),ou)     - inquired data
c
#include "impcom.inc"
#include "constants.inc"
#include "elem_wrk.inc"
#include "in_mem.inc"
c
      integer          iProcSMP
      integer          matRecordIndex, kMatRecPt, nTens
      double precision vect(*)
      integer          numWordsPt, i, k1, kk, iptr, dptr
c

      iptr = 1 + E_LenMatRecIndx*iProcSMP()

      if (M_MatRecI(iptr+matRecordIndex) .gt. 0) then
          numWordsPt = M_MatRecI(iptr)
          dptr = E_PtrMatRecDp + E_LenMatRec*iProcSMP()
          dptr = dptr - 1
          kk = dptr + (kMatRecPt-1) * numWordsPt
          k1 = kk + M_MatRecI(iptr+matRecordIndex) - 1
      else
          nTens = 0
          goto 999
      endif

      do i = 1,nTens
         vect(i) = M_MatRecDp(k1+i)
      enddo



  999 continue

      return
      end

