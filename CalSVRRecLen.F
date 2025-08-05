      subroutine CalSVRRecLen (elChar,elData,mat,lSVR)
c --- primary function: determine record lengths for the SVR file
c
c  input arguments:
c     elChar   (int,ar(*),in)    - the element type descriptor
c     mat      (int,sc,in)       - material number for this element
c
c  output arguments:
c     lSVR     (int,sc,out)      - the integer length on SVR
c
#include "impcom.inc"
#include "constants.inc"
#include "sec_info.inc"
#include "SVR_file.inc"
c
      integer elChar(IELCSZ),mat,lSVR,lenAdd,lenmx,jtype,ireal,
     &        numwd,nmatMisc,elData(EL_DIM),GetSecChars,numSectPts

      integer*8  nlkword(W_nMatCol,W_nMatRow)
#if defined(DEBUG)
      call SubBeg('CalSVRRecLen')
#endif

      lSVR = 0
c --- skip if dead or no SVRs
      if (elData(EL_DEAD) .eq. 1) goto 999
      if (elChar(NMPSVR).le.0 .and. elChar(NMSSVR).le.0) goto 999

      jtype = elChar(JETYP)
      if (jtype.le.0) goto 999
      
c --- define SVR file size
c --- general structural SVRs
      if (elChar(NMSSVR).gt.0) then
          lSVR = lSVR + W_IntPerDp*elChar(NMSSVR) + 3
          ireal = elData(EL_REAL)
          if (jtype.eq.W_eBeam2nd .or. jtype.eq.W_eBeam3nd) then
              if (elChar(KYOP7).eq.1 .and.
     &        .not.(B_PtrRealN.ne.PTRFTNNULL .and. B_RealN(ireal).gt.0)
     &           ) then
                  numSectPts = GetSecChars(ireal,1) 
                  lSVR = lSVR + 20 + W_IntPerDp*B_NumResPts*NumSectPts
              endif
          endif
      endif
c --- inelasticity SVRs
      if (elChar(NMPSVR).gt.0) then
          lenAdd = 0
c --- get nonlinear material keys
          call GetMatKeys (mat,nlkword(1,1))
          call RetMatRecLen (elChar(1),elData(1),nlkword(1,1),
     &                       lenAdd,lenmx,numwd,nmatMisc)
          if (lenmx.gt.0) lSVR = lSVR + W_IntPerDp*lenAdd + 3
      endif
c --- index (if not null)
c      if (lSVR.gt.0) lSVR = lSVR + 3 + S_nIndx

  999 continue

#if defined(DEBUG)
      call SubEnd('CalSVRRecLen')
#endif

      return
      end
