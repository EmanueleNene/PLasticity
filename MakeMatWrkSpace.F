      subroutine MakeMatWrkSpace (leneMatRec)
c --- primary function: create elements and material workspaces
c
c  output arguments:
c     leneMatRec   (int,sc,out)  - maximum length of material record
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
c
      integer          elmget, etyget

      integer          leneMatRec


      integer          i, ii, ityp, lenWrk, elemId, matId,
     &                 lenAdd, numWordsPt, nmatMisc,
     &                 elChar(IELCSZ), elData(EL_DIM), nodes(W_NdMax)
      integer*8        nlkword(W_nMatCol,W_nMatRow)


#if defined(DEBUG)
      call SubBeg('MakeMatWrkSpace')
#endif

      leneMatRec = 0

      if (D_NumEl .le. 0)  goto 999

      do elemId = 1,D_NumEl
c --- retrieve element attributes
         i = elmget (elemId,elData(1),nodes(1))
         if (i .le. 0) cycle
         ityp  = elData(EL_TYPE)
         i = etyget (ityp,elChar(1))
         if (i.le.0) cycle
c --- skip if elements do not use this workspace
         if (elChar(NMPSVR).eq.0) cycle
         matId = elData(EL_MAT)
         call GetMatKeys (matId,nlkword(1,1))
         call RetMatRecLen (elChar(1), elData(1), nlkword(1,1),
     &                      lenAdd, lenWrk, numWordsPt, nmatMisc)

         if (leneMatRec.lt.lenWrk)  leneMatRec = lenWrk

      enddo

      ii = leneMatRec/8 + 1
      leneMatRec = ii*8


  999 continue

#if defined(DEBUG)
      call SubEnd('MakeMatWrkSpace')
#endif

      return
      end
