      subroutine CalSect_10 (modelR,NNODI,numElems,numNodes,elemConn,
     &                       len_stiff,address,numSectPts, Results,
     &                       Summary, ErrorFlag)
c --- primary function: compute section data
c
#include "impcom.inc"
#include "constants.inc"
#include "sec_info.inc"
c
      integer NNODI,ixx,gProc
      integer
     &     numSectPts,len_stiff,address(*),
     &     numIntPts,numElems,numNodes,elemConn(9,*)
      integer
     &     rXcoord,
     &     rYcoord,
     &     rGxy
      integer
     &     pArea,
     &     pCgx,
     &     pCgy,
     &     pIx,
     &     pIy,
     &     pIxx,
     &     pIxy,
     &     pIyy,
     &     pJ,
     &     pIpp,
     &     pSw0,
     &     pSwx,
     &     pSwy,
     &     pSww,
     &     pSx,
     &     pSy,
     &     pLast
      parameter (
     &     pArea = 1,
     &     pCgx  = 2,
     &     pCgy  = 3,
     &     pIx   = 4,
     &     pIy   = 5,
     &     pIxx  = 6,
     &     pIxy  = 7,
     &     pIyy  = 8,
     &     pJ    = 9,
     &     pIpp  = 10,
     &     pSw0  = 11,
     &     pSwx  = 12,
     &     pSwy  = 13,
     &     pSww  = 14,
     &     pSx   = 15,
     &     pSy   = 16,
     &     pLast = 16)

      parameter (numIntPts = 4)
      integer ErrorFlag,iModWt

      double precision modelR(*),Results(B_NumResPts,*),Summary(20),
     &                 work, scL(3,3),scLI(3,3),scR(3),sc(3),avWarp

      double precision MW1(B_MaxLenStiff), MW2(B_MaxSecNds), 
     &                 ModWt(B_MaxElems), 
     &                 tmp(B_MaxSecNds*(B_MaxSecNds+1))

      rXcoord = 1
      rYcoord = numNodes+1
      rGxy = 2*numNodes+1

      call v_zero(Summary(1),pLast)
      call CalSect_2 (numElems,modelR(rGxy),Summary(1),ModWt(1),
     &                ErrorFlag)
      if (ErrorFlag .eq. 1) goto 999
      call v_zero(Results(1,1),B_NumResPts*numElems*numIntPts)
      gProc = 1
      call CalSect_5 (NNODI,len_stiff,numNodes,numElems,address,
     &                elemConn,modelR(rXcoord),modelR(rYcoord),ModWt(1),
     &                MW1(1),MW2(1),Summary(1),ErrorFlag,gProc)
      if (ErrorFlag .eq. 1) goto 999
          call CalSect_1 (MW1(1), address, len_stiff, MW2(1), 
     &                    numNodes, tmp(1), ErrorFlag)
      if (ErrorFlag .eq. 1) goto 999
      if (numNodes .gt. 0) then
          call CalSect_9(NNODI,gProc,numNodes,numElems,elemConn,
     &                   modelR(rXcoord),modelR(rYcoord),
     &                   MW2(1),ModWt(1),Results(1,1))
          numSectPts = numElems * numIntPts
      endif
c --- calculate centroid and modify coordinates
      call v_zero(Summary(1), pLast)
      do ixx=1,numSectPts
         if (Results(7,ixx).gt.0.0d0) then
             Summary(pArea) = Summary(pArea) + Results(7,ixx)
             Summary(pIy)   = Summary(pIy) + Results(7,ixx)*
     &            Results(6,ixx)
             Summary(pIx)   = Summary(pIx) + Results(7,ixx)*
     &            Results(5,ixx)
         endif
      enddo
      
      Summary(pCgx) = Summary(pIx)/Summary(pArea)
      Summary(pCgy) = Summary(pIy)/Summary(pArea)
      do ixx = 1,numSectPts
         Results(5,ixx) = Results(5,ixx) - Summary(pCgx)
         Results(6,ixx) = Results(6,ixx) - Summary(pCgy)
      enddo
c --- calculate shear center
      Summary(pIx)  = 0.0d0
      Summary(pIy)  = 0.0d0
      Summary(pIxx) = 0.0d0
      Summary(pIyy) = 0.0d0
      Summary(pIxy) = 0.0d0
      Summary(pSw0) = 0.0d0
      Summary(pSwx) = 0.0d0
      Summary(pSwy) = 0.0d0
      do ixx=1,numSectPts
         if (Results(7,ixx) .gt. 0.0d0) then
             Summary(pIxx) = Summary(pIxx) +
     &            Results(6,ixx)**2 * Results(7,ixx)
             Summary(pIxy) = Summary(pIxy) +
     &            Results(6,ixx)*
     &            Results(5,ixx)* Results(7,ixx)
             Summary(pIyy) = Summary(pIyy) +
     &            Results(5,ixx)**2 * Results(7,ixx)
             Summary(pSw0) = Summary(pSw0)+
     &            Results(2,ixx)*Results(7,ixx)
             Summary(pSwx) = Summary(pSwx)+
     &            Results(5,ixx)*Results(2,ixx)*Results(7,ixx)
             Summary(pSwy) = Summary(pSwy)+
     &            Results(6,ixx)*Results(2,ixx)*Results(7,ixx)
         endif
      enddo
      call v_zero(scL(1,1),9)
      scL(1,1) = Summary(pArea)
      scL(1,2) = 0.0d0
      scL(1,3) = 0.0d0
      scL(2,1) = scL(1,2)
      scL(2,2) = Summary(pIxx)
      scL(2,3) =-Summary(pIxy)
      scL(3,1) = scL(1,3)
      scL(3,2) = scL(2,3)
      scL(3,3) = Summary(pIyy)
      scR(1) = Summary(pSw0)
      scR(2) =-Summary(pSwy)
      scR(3) = Summary(pSwx)
      call invA3x3(scL(1,1),scLI(1,1),work)
      
      do ixx=1,3
         sc(ixx)=
     &        scLI(ixx,1)*scR(1)+scLI(ixx,2)*scR(2)+scLI(ixx,3)*scR(3)
      enddo
      
      Summary(pSx) = sc(2)
      Summary(pSy) = sc(3)
      Summary(pJ)   = 0.0d0
      Summary(pSww) = 0.0d0
      do ixx=1,numSectPts
         if (Results(7,ixx) .gt. 0.0d0) then
            Results(2,ixx)=Results(2,ixx)-Summary(pSw0)/Summary(pArea)
            work=(Summary(pIxx)*Summary(pIyy)-Summary(pIxy)**2)
            sc(2)=(Summary(pSwy)*Summary(pIxx)-
     &             Summary(pSwx)*Summary(pIxy))*Results(5,ixx)/work
            sc(3)=(Summary(pSwx)*Summary(pIyy)-
     &             Summary(pSwy)*Summary(pIxy))*Results(6,ixx)/work
            work=((Results(2,ixx)-sc(2)-sc(3))**2)*Results(7,ixx)
            Summary(pSww) = Summary(pSww) + work
            Summary(pJ) = Summary(pJ) +
     &          (Results(5,ixx)**2+Results(6,ixx)**2)*Results(7,ixx)
     &          -Results(6,ixx)*Results(3,ixx)*Results(7,ixx)+
     &           Results(5,ixx)*Results(4,ixx)*Results(7,ixx)
            Summary(pIpp) = Summary(pIpp) + Results(7,ixx)*
     &           (Results(5,ixx)**2+Results(6,ixx)**2)
         endif
      enddo
c --- resolve for warping function w.r.t shear center as pole
      do ixx=1,numNodes
         modelR(rXcoord+ixx-1)=modelR(rXcoord+ixx-1)-Summary(pSx)
         modelR(rYcoord+ixx-1)=modelR(rYcoord+ixx-1)-Summary(pSy)
      enddo
      gProc = 2
      call CalSect_5(NNODI,len_stiff,numNodes,numElems,address,
     &               elemConn,modelR(rXcoord),modelR(rYcoord),ModWt(1),
     &               MW1(1),MW2(1),Summary(1),ErrorFlag,gProc)
      if (ErrorFlag .eq. 1) goto 999
      call CalSect_1 (MW1(1),address,len_stiff,MW2(1),
     &                numNodes,tmp(1),ErrorFlag)
      if (ErrorFlag .eq. 1) goto 999
      if (numNodes .gt. 0) then
          call CalSect_9(NNODI,gProc,numNodes,numElems,elemConn,
     &                   modelR(rXcoord),modelR(rYcoord),
     &                   MW2(1),ModWt(1),Results(1,1))
      endif
      Summary(pSw0) = 0.0d0
      Summary(pSwx) = 0.0d0
      Summary(pSwy) = 0.0d0
      Summary(pSww) = 0.0d0
      Summary(pJ)   = 0.0d0
      do ixx=1,numSectPts
         if (Results(7,ixx) .gt. 0.0d0) then
            Summary(pSw0) = Summary(pSw0)+
     &           Results(2,ixx)*Results(7,ixx)
         endif
      enddo
      avWarp = Summary(pSw0)/Summary(pArea)
      Summary(pSw0) = 0.0d0
      do ixx=1,numSectPts
         if (Results(7,ixx).gt.0.0d0) then
            Results(2,ixx)=Results(2,ixx)- avWarp
            Summary(pSww)=Summary(pSww)+Results(2,ixx)**2*Results(7,ixx)
            Summary(pJ) = Summary(pJ) +
     &           (Results(5,ixx)**2+Results(6,ixx)**2)*Results(7,ixx)
     &           -Results(6,ixx)*Results(3,ixx)*Results(7,ixx)+
     &            Results(5,ixx)*Results(4,ixx)*Results(7,ixx)
            Summary(pSw0) = Summary(pSw0)+
     &           Results(2,ixx)*Results(7,ixx)
            Summary(pSwx) = Summary(pSwx)+
     &           Results(5,ixx)*Results(2,ixx)*Results(7,ixx)
            Summary(pSwy) = Summary(pSwy)+
     &           Results(6,ixx)*Results(2,ixx)*Results(7,ixx)
         endif
      enddo
      if (ErrorFlag .eq. -99 ) then
          Summary(pSw0) = 0.0d0
          Summary(pSwx) = 0.0d0
          Summary(pSwy) = 0.0d0
          Summary(pSww) = 0.0d0
          Summary(pSx) = Summary(pCgx)
          Summary(pSy) = Summary(pCgy)
      endif
c
      do ixx = 1,numNodes
         modelR(rXcoord+ixx-1) = modelR(rXcoord+ixx-1)
     &                         + Summary(pSx)-Summary(pCgx)
         modelR(rYcoord+ixx-1) = modelR(rYcoord+ixx-1)
     &                         + Summary(pSy)-Summary(pCgy)
      enddo
      gProc = 3
      iModWt = numElems+1
      call CalSect_5(NNODI, len_stiff,numNodes,numElems,address,
     &               elemConn, modelR(rXcoord),modelR(rYcoord),
     &               ModWt(iModWt),MW1(1),MW2(1),Summary(1),
     &               ErrorFlag,gProc)
      if (ErrorFlag .eq. 1) goto 999
      call CalSect_1 (MW1(1),address,len_stiff,MW2(1),
     &                numNodes,tmp(1),ErrorFlag)
      if (ErrorFlag .eq. 1) goto 999
      if (numNodes .gt. 0) then
          call CalSect_9(NNODI,gProc,numNodes,numElems,elemConn,
     &                   modelR(rXcoord),modelR(rYcoord),
     &                   MW2(1),ModWt(iModWt),Results(1,1))
      endif
      gProc = 4
      iModWt = numElems+1
      call CalSect_5(NNODI, len_stiff,numNodes,numElems,address,
     &               elemConn, modelR(rXcoord),modelR(rYcoord),
     &               ModWt(iModWt),MW1(1),MW2(1),Summary(1),
     &               ErrorFlag,gProc)
      if (ErrorFlag.eq.1) goto 999
      call CalSect_1 (MW1(1),address,len_stiff,MW2(1),
     &                numNodes,tmp(1),ErrorFlag)
      if (ErrorFlag .eq. 1) goto 999
      if (numNodes .gt. 0) then
          call CalSect_9(NNODI,gProc,numNodes,numElems,elemConn,
     &                   modelR(rXcoord), modelR(rYcoord),
     &                   MW2(1),ModWt(iModWt),Results(1,1))
      endif
      Summary(17)=0.0d0
      Summary(18)=0.0d0
      Summary(19)=0.0d0
      do ixx=1,numSectPts
         if (Results(7,ixx) .gt. 0.0d0) then
             Summary(17) = Summary(17)+Results(7,ixx)*
     &                    (Results(8,ixx)**2+Results(9,ixx)**2)
             Summary(18) = Summary(18)+Results(7,ixx)*
     &                    (Results(8,ixx)*Results(10,ixx)+
     &                     Results(9,ixx)*Results(11,ixx))
             Summary(19) = Summary(19)+Results(7,ixx)*
     &                    (Results(10,ixx)**2+Results(11,ixx)**2)
         endif
      enddo
      work = 1.0d0/
     &     (Summary(pArea)*(Summary(17)*Summary(19)-Summary(18)**2))
      Summary(17) = Summary(17)*work
      Summary(18) =-Summary(18)*work
      Summary(19) = Summary(19)*work

 999  continue
      
      return
      end
