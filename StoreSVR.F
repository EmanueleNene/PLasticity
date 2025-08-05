      subroutine StoreSVR (svindx,nset,nsvr,svr)
c --- primary function: write out a svr data set for an element
c
c  input arguments:
c     nset     (int,sc,in)          - the set number in this index (same as GetSVR)
c     nsvr     (int,sc,in)          - number of dp words in this set
c     svr      (dp,ar(nsvr),in)     - data in this set
c
c  output arguments:
c     svindx   (int,ar(10,2),inout) - updated index
c
#include "impcom.inc"
#include "constants.inc"
#include "SVR_file.inc"
#include "solu_info.inc"
#include "in_mem.inc"
c
      integer*8   MakeLongInt

      integer*8  loc8
      integer    svindx(20), nset, nsvr, rloc
      double precision  svr(*)

#if defined(DEBUG)
      call SubBeg('StoreSVR')
#endif

c --- only write if requested!
      if (A_KeySvrUpd.eq.0) goto 999

c --- retrieve the file location of this set
      rloc = svindx(nset+3)
      if (rloc.lt.0) then
c --- location not defined: set the location and write (also update loc)
          loc8 = MakeLongInt(svindx(1),svindx(2))
          loc8 = loc8 + svindx(3)
          svindx(nset+3) = svindx(3)
          if (S_InCore) then
              call v_move (svr(1),M_SvrD(loc8+1),nsvr)
          else
              call lockSMP (W_LockSvr)
              call Buf_Write (S_BlockNum,loc8,nsvr*W_IntPerDp,svr(1),0)
              call unlockSMP (W_LockSvr)
          endif
          svindx(3) = svindx(3) + nsvr*W_IntPerDp + 3
      else
c --- location defined (written once before): overwrite the data
          loc8 = MakeLongInt(svindx(1),svindx(2)) + rloc
          if (S_InCore) then
              call v_move (svr(1),M_SvrD(loc8+1),nsvr)
          else
              call lockSMP (W_LockSvr)
              call Buf_Write (S_BlockNum,loc8,nsvr*W_IntPerDp,svr(1),0)
              call unlockSMP (W_LockSvr)
          endif
      endif

  999 continue

#if defined(DEBUG)
      call SubEnd('StoreSVR')
#endif

      return
      end
