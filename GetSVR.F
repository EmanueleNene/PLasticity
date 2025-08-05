      subroutine GetSVR (svindx,nset,nsvr,svr)
c --- primary function: get svr data set for an element
c
c  input arguments:
c     svindx   (int,ar(20),in)   - index for svr for this element
c     nset     (int,sc,in)       - the set number in this index
c
c  output arguments:
c     nsvr     (int,sc,inout)    - number of dp words in this set
c     svr      (dp,ar(nsvr),in)  - data in this set
c
#include "impcom.inc"
#include "constants.inc"
#include "SVR_file.inc"
#include "in_mem.inc"
c
      integer*8   MakeLongInt

      integer*8   loc8
      integer  svindx(20), nset, nsvr, i,kbf,rloc
      double precision  svr(*)

#if defined(DEBUG)
      call SubBeg('GetSVR')
#endif

c --- get the data
c --- retrieve the pointer to the requested svr set
      rloc = svindx(nset+3)
      if (rloc.ge.0) then
c --- this svr set is defined: retrieve it
          loc8 = MakeLongInt(svindx(1),svindx(2)) + rloc
          if (S_InCore) then
              call v_move (M_SvrD(loc8+1),svr(1),nsvr)
          else
              i = 10
              call lockSMP (W_LockSvr)
              call Buf_Read (S_BlockNum,loc8,i,svr(1),kbf)
              call unlockSMP (W_LockSvr)
              nsvr = i/W_IntPerDp
          endif
      else
c --- this svr set was not previously defined: initialize it
          call v_zero (svr(1),nsvr)
      endif

#if defined(DEBUG)
      call SubEnd('GetSVR')
#endif

      return
      end
