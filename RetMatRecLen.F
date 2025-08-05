      subroutine RetMatRecLen (elChar,elData,nlkword,
     &                         lenAdd,lenmx,numWordsPt,nmatMisc)
c --- primary function: determine the maximum length of material record required
c                       for all materials used by this element
c  input arguments:
c     elID          (int,sc,in)          - element label (number)
c     elChar        (int,ar(*),in)       - array of element type characteristics
c     elData        (int,ar(10),in)      - array of element data
c
c  output arguments:
c     nlkword       (intL,ar(2,*),inout) - material flags
c     lenAdd        (int,sc,out)         - length of add. record at int. points
c     numWordsPt    (int,sc,out)         - total double precision words of
c                                          material record at integration point
c     lenmx         (int,sc,out)         - numWordsPt * number of integration points
c
#include "impcom.inc"
#include "constants.inc"
#include "elem_wrk.inc"
c

      integer          elChar(*), lenmx,  numWordsPt, lenAdd,
     &                 nmatMisc, mRecIntPidx(W_LenMatRecIndx),
     &                 elData(EL_DIM)
      integer*8        nlkword(W_nMatCol,W_nMatRow)



#if defined(DEBUG)
      call SubBeg('GetMatRecLen')
#endif

      call CalMatRecLen (nlkword(1,1), elChar(1), elData(1),
     &                   lenmx, numWordsPt, nmatMisc, mRecIntPidx)

      if (lenmx.gt.0) lenAdd = lenmx  + E_MatRecAdd

#if defined(DEBUG)
      call SubEnd('GetMatRecLen')
#endif

      return
      end
