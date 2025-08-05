      subroutine PutSVRIndx (elem,locsvr,svindx)
c --- primary function:    write the svr element index onto file
c --- secondary functions: update the locsvr pointer to next element
c
c  input arguments:
c     locsvr   (int*8,sc,inout)      - pointer to start of SVRs for element
c     svindx   (int,ar(S_nIndx),in)  - index to SVRs for this element
c                                      low and high parts of 64 bit address
c
c  output arguments:
c     locsvr   (int*8,sc,inout)      - pointer to start of SVRs for next element
c
#include "impcom.inc"
#include "constants.inc"
#include "SVR_file.inc"
#include "solu_info.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer*8   MakeLongInt

      integer*8  locsvr, loc8
      integer    elem, svindx(2+S_nIndx), j, m
      
#if defined(DEBUG)
c      call SubBeg('PutSVRIndx')
#endif

      if (S_BlockNum.ne.-1 .and. A_KeyFirstPass.eq.1) then
          m = S_nIndx*(elem-1) + D_NumEl*2
          do j = 1,S_nIndx
             M_SvrIndx(m+j) = svindx(j+2)
          enddo
          locsvr = MakeLongInt(svindx(1),svindx(2)) + svindx(3)
      endif

#if defined(DEBUG)
c      call SubEnd('PutSVRIndx')
#endif

      return
      end
