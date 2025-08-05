      subroutine PutMatRec (matRecordIndex, kMatRecPt0, nTens, vect)
c --- primary function: put material record data
c
c  input arguments:
c     matRecordIndex (in,sc,in)        - material record index
c     kMatRecPt      (in,sc,in)        - integration point where data to be inquired
c     elemId         (in,sc,in)        - element number
c     nTens          (in,sc,in)        - number of components to be inquired
c     numWordsPt     (in,ar(*),in)     - number of material records
c
c  output arguments:
c     vect           (dp,ar(*),ou)     - inquired data
c
#include "impcom.inc"
#include "constants.inc"
#include "elem_wrk.inc"
#include "in_mem.inc"
c
c
      integer          iProcSMP
c
      integer          matRecordIndex, kMatRecPt0, nTens
      double precision vect(*)
c
      integer          numWordsPt, i, k1, kk, iptr, dptr,
     &                 iproc, kMatRecPt


      iproc = iProcSMP()
      iptr = 1 + E_LenMatRecIndx*iproc

      kMatRecPt = kMatRecPt0
      k1 = 0
c
      if (M_MatRecI(iptr+matRecordIndex) .gt. 0) then
          numWordsPt = M_MatRecI(iptr)
          dptr = E_PtrMatRecDp + E_LenMatRec*iproc
          dptr = dptr - 1
          kk = dptr + (kMatRecPt-1) * numWordsPt
          k1 = kk + M_MatRecI(iptr+matRecordIndex) - 1
      elseif (nTens .eq. 0) then
          goto 999
      endif

      do i = 1, nTens
         M_MatRecDp(k1+i) = vect(i)
      enddo


  999 continue

      return
      end
