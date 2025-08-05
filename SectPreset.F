      subroutine SectPreset (NumEl)
c --- primary function: allocate space for contact
c
c  input arguments:
c     NumEl (int,sc,in)          - number of elements
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "sec_info.inc"
c
      integer NumEl
      PTRFTN  MemAlloc, MemRealloc
      integer GetSecChars,elmget,etyget,rlinqr,elem,nnod,ireal,jtype,
     &        elData(EL_DIM),nodes(W_NdMax),elChar(IELCSZ),NumRstNodes,
     &        ityp,itypold,maxr,MemLen,PtrNext,numSectPts
      logical IsBeam
      integer  secDataI(B_IntDatLen + 1000) 
      double precision secDataR(B_MaxSecNds*3  + B_MaxSecPts/4)
      double precision sectSummary(20)
      double precision secResults(B_NumResPts*B_MaxSecPts)
      
#if defined(DEBUG)
      call SubBeg('SectPreset')
#endif

      B_PtrRealN  = PTRFTNNULL
      B_PtrSectData = PTRFTNNULL
c --- set up the section data
      maxr = rlinqr(0,DB_MAXDEFINED)
      if (maxr.le.0) then
          call erhandler ('SectPreset', 5000, 3, 'No real constants 
     &                     defined',0.0d0,' ')
          call StopAll (-5)          
      endif
          
      B_PtrRealN = MemAlloc(maxr,W_LenInt,'B_PtrRealN      ')
      call vi_zero(B_RealN(1),maxr)
      itypold = 0
      PtrNext = 1
      MemLen = 0
      do elem = 1,NumEl
         nnod = elmget(elem,elData(1),nodes(1))
         if (nnod.le.0) cycle
         if (elData(EL_DEAD).eq.1) cycle
         ireal = elData(EL_REAL)
         if (ireal.le.0) cycle
         If (B_RealN(ireal).ne.0) cycle
         ityp = elData(EL_TYPE)
         if (ityp .ne. itypold) then
             itypold = ityp
             ityp = etyget(ityp,elChar(1))
         endif
c --- only beams are processed
         jtype = elChar(JETYP)
         IsBeam = jtype.eq.W_eBeam2nd .or. jtype.eq.W_eBeam3nd
         if (.not.IsBeam) cycle
         
         numSectPts = GetSecChars(ireal,1)
         MemLen = MemLen + B_NumResPts*numSectPts + 20
         
         if (B_PtrSectData.eq.PTRFTNNULL) then
             B_PtrSectData = MemAlloc(MemLen,
     &                                W_LenDouble,'B_PtrSectData   ')
         else
             B_PtrSectData = MemRealloc(B_PtrSectData,MemLen,
     &                                  W_LenDouble,'B_PtrSectData   ')
         endif
         call CalBeamSect (ireal,elChar(KYOP1),secDataI(1),secDataR(1),
     &                     sectSummary(1),secResults(1),NumRstNodes,
     &                     B_SectData(PtrNext))
         B_RealN(ireal) = PtrNext         
         PtrNext = MemLen + 1
      enddo

#if defined(DEBUG)
      call SubEnd('SectPreset')
#endif
      return
      end