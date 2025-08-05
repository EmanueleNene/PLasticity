      subroutine CalResuLen (elData,nlkword,eomask,lerst)
c --- primary function: compute the length of RST file size for individual elems
c
c  input arguments:
c     nlkword    (intL,ar(2,*),in)         - plasticity types/options
c     eomask     (int,sc,in)               - element output mask
c
c  output arguments:
c     lerst      (int,sc,in/out)           - length of rst file
c
#include "impcom.inc"
#include "constants.inc"
#include "RST_file.inc"
c
      integer          GetSecChars, etyiqr

      integer          eomask, lerst, elData(EL_DIM)
      integer*8        nlkword(W_nMatCol,W_nMatRow)

      integer          intrst, ncomp, nlocs,
     &                 NumRstNodes, NumSectPts, numElemNodes,
     &                 beamRstNodesLen, numSectPtsLen
      integer*8        keyplL, keycrL, keynlL

#if defined(DEBUG)
      call SubBeg('CalResuLen')
#endif

      intrst = 2

      keyplL = nlkword(1,1)
      keycrL = nlkword(1,2)
      keynlL = keyplL + keycrL

c --- define the number of stress components
      ncomp = 3
c --- define number of components and number of locations
      NumRstNodes = GetSecChars(elData(EL_REAL),2)
      NumSectPts = GetSecChars(elData(EL_REAL),1)
c --- was checked with IsQuadBeam before
      if (etyiqr(elData(EL_TYPE),-JETYP) .eq. W_eBeam3nd) then
          numElemNodes = 3
          nlocs = NumSectPts * 2
          numSectPtsLen = NumSectPts * ncomp * 2
      else
          numElemNodes = 2      
          nlocs = NumSectPts
          numSectPtsLen = NumSectPts * ncomp
      endif
      beamRstNodesLen = NumRstNodes * numElemNodes * ncomp

c --- nodal stresses
      if (btest(eomask,W_OutStress)) then
          lerst = lerst + intrst*beamRstNodesLen + 3
      endif

c --- nodal elastic strains
      if (btest(eomask,W_OutElaStra)) then
          lerst = lerst + intrst*beamRstNodesLen + 3
      endif

c --- nodal plastic strains
      if (keyplL.gt.0) then
          if (btest(eomask,W_OutPlStra)) then
              lerst = lerst + (intrst*ncomp*nlocs + 3)
          endif
      endif

c --- nodal creep strains
      if (keycrL.ne.0) then
          if (btest(eomask,W_OutCrStra)) then
              lerst = lerst + intrst*ncomp*nlocs + 3
          endif
      endif

c --- nodal hygro strains
      if (btest(eomask,W_OutGyStra)) then
          lerst = lerst + intrst*NumRstNodes * numElemNodes + 3
      endif

c --- nonlinear data
      if (keynlL.ne.0) then
          lerst = lerst + intrst*  6 * nlocs + 3
      endif


#if defined(DEBUG)
      call SubEnd('CalResuLen')
#endif

      return
      end
