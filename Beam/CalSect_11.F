      subroutine CalSect_11(secDataR,NNODI,numElems,numNodes,elemConn,
     &                      numSectPts,Results,Summary,ErrorFlag, 
     &                      nSecRstNodes)
c --- primary function: compute section data
c
#include "impcom.inc"
#include "constants.inc"
#include "sec_info.inc"
c
      integer numSectPts,NNODI,numElems,numNodes,
     &        elemConn(9,numElems),ErrorFlag,len_stiff,nSecRstNodes
      double precision secDataR(*),Results(*),Summary(*)

      integer Add(B_MaxSecNds), ColHt(B_MaxSecNds)
      double precision MR(B_DpDatLen)

      call v_move(secDataR(1),MR(1),B_DpDatLen)

      call CalSect_8 (NNODI,numNodes,numElems,elemConn,
     &                len_stiff,Add(1),ColHt(1),nSecRstNodes)
      call CalSect_10 (MR(1),NNODI,numElems,numNodes,elemConn,
     &                 len_stiff,Add(1),numSectPts,
     &                 Results,Summary,ErrorFlag)

      return
      end




