      subroutine GetSVRIndx (elem, locsvr, svindx)
c --- primary function: get the index for SVRs
c
c  input arguments:
c     locsvr   (int*8,sc,in)    - pointer to location of index
c
c  output arguments:
c     svindx   (int,ar(S_nIndx),out) - the word index of SVR variables
c
#include "impcom.inc"
#include "constants.inc"
#include "SVR_file.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer*8 Buf_Inq, locsvr, loc8, krw
      integer   elem, svindx(2+S_nIndx), j, m, kbf

#if defined(DEBUG)
      call SubBeg('GetSVRIndx')
#endif

c --- set to zero if file not being read
      if (S_BlockNum.eq.-1)  then
          call vi_zero (svindx(1),S_nIndx)
          goto 999
      endif

      krw = Buf_Inq(S_BlockNum)
c --- is this file being written?
      if (krw.ne.2) then
c --- this index has been written once before: read it in
          call MakeShortInt (locsvr+S_RecStart,svindx(1),svindx(2))
          m = S_nIndx*(elem-1) + D_NumEl*2
          do j = 1,S_nIndx
             svindx(j+2) = M_SvrIndx(m+j)
          enddo
      else
c --- this index has not been written before: initialize it 
c     and mark that the onset points are not yet initialized
          call vi_init (svindx(4),S_nIndx-1,-1)
          call MakeShortInt (locsvr+S_RecStart,svindx(1),svindx(2))
          svindx(3) = 0
      endif
  999 continue

#if defined(DEBUG)
      call SubEnd('GetSVRIndx')
#endif

      return
      end
