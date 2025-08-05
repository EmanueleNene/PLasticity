      subroutine OvrToSvr
c --- primary function: restore OVR file to SVR File (or to memory)
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "SVR_file.inc"
#include "OVR_file.inc"
#include "dof_info.inc"
c
      integer   Bin_Set
      integer*8   MakeLongInt
      PTRFTN    MemAlloc, MemAllocL

      integer*8   maxrecL, iL
      integer   i, k, ivect(100), kbf, nSVR
      double precision   dvect(20)
      character*(B_PathLen) pname

      pointer (pdWorkL,dWork)
      double precision  dWork(*)

#if defined(DEBUG)
      call SubBeg('OvrToSvr')
#endif

c --- check for memory resident (not used any more)
      if (M_PtrOvr.ne.PTRFTNNULL) then
          call v_move (M_Ovr(1),dvect(1),20)
c --- we don't need to restore A_TimeCur and A_TimeInc, it is
c --- already in memory in this case, and they may be overwritten
          A_TimeBeg = dvect(4)
          A_TimeEnd = dvect(5)
          i = nint(dvect(19))
          i = i*A_NumSolVects
          call v_move (M_Ovr(21),M_Disp(1),i)
          i = i + 21
          call v_move (M_Ovr(i),M_ERot(1,1),3*D_BacLen)
          goto 999
      endif

c --- open OVR file
      call GetPath (O_Name,pname)
      O_BlockNum = 2
      O_Unit = 17
c --- force file type of OVR file to be same as SVR file
      i = Bin_Set(O_BlockNum,O_Unit,W_Read,B_NumOfPages,pname)

c --- restore OVR header and write to SVR file
      O_FilePos = 0
      S_FilePos = 0
      i = 100
      call Buf_Read  (O_BlockNum,O_FilePos,i,ivect(1),kbf)
      call Buf_Write (S_BlockNum,S_FilePos,i,ivect(1),i)

c --- restore OVR integer header
      i = 80
      call Buf_Read (O_BlockNum,O_FilePos,i,ivect(1),kbf)
      if (i.lt.80) call vi_zero (ivect(41),40)
      nSVR = i

c --- get the maximum record written and write integer header to SVR file
      maxrecL = MakeLongInt(ivect(40),ivect(39))
      call Buf_Write (S_BlockNum,S_FilePos,nSVR,ivect(1),nSVR)

c --- restore OVR dp header and write to SVR file
      i = 20
      call Buf_Read  (O_BlockNum,O_FilePos,i,dvect(1),kbf)
      call Buf_Write (S_BlockNum,S_FilePos,i,dvect(1),0)

c --- NOTE: other data such as the D_DofList/M_BackList/M_ElemOrd/M_DofBits/Indexes
c ---       arrays do not need to be copied (they are already on the SVR record)
c --- restore OVR u vectors and write to SVR file
      O_FilePos = MakeLongInt(ivect(36),ivect(46))
      S_FilePos = O_FilePos
      k = D_ULen*W_IntPerDp

      call v_zero (M_Disp(1),D_ULen*A_NumSolVects)
      call Buf_Read  (O_BlockNum,O_FilePos,k,M_Disp(1),kbf)
      call Buf_Write (S_BlockNum,S_FilePos,k,M_Disp(1),0)
      call Buf_Read  (O_BlockNum,O_FilePos,k,M_IncI(1),kbf)
      call Buf_Write (S_BlockNum,S_FilePos,k,M_IncI(1),0)
      call Buf_Read  (O_BlockNum,O_FilePos,k,M_IncS(1),kbf)
      call Buf_Write (S_BlockNum,S_FilePos,k,M_IncS(1),0)
      call Buf_Read  (O_BlockNum,O_FilePos,k,M_Velo(1),kbf)
      call Buf_Write (S_BlockNum,S_FilePos,k,M_Velo(1),0)
      if (M_PtrAcce.ne.PTRFTNNULL .and. A_NumSolVects.eq.5) then
          call Buf_Read  (O_BlockNum,O_FilePos,k,M_Acce(1),kbf)
          call Buf_Write (S_BlockNum,S_FilePos,k,M_Acce(1),0)
      endif
c --- restore OVR global rotations and write to SVR file (if necessary)
      if (ivect(10).gt.0 .and. M_PtrERot.ne.PTRFTNNULL) then
          call v_zero (M_ERot(1,1),3*D_BacLen)
          k = 3*D_BacLen*W_IntPerDp
          call Buf_Read  (O_BlockNum,O_FilePos,k,M_ERot(1,1),kbf)
          call Buf_Write (S_BlockNum,S_FilePos,k,M_ERot(1,1),0)
      endif

c --- loop thru all element saved data on OVR file and write to SVR file
      O_FilePos = MakeLongInt(ivect(37),ivect(47))
      if (S_InCore) then
          iL = S_FileLengL
          call Buf_ReadL (O_BlockNum,O_FilePos,iL,M_SvrD(1),kbf)
      else
          i = 10000 ! same in SvrToOvr
          pdWorkL = MemAlloc(i,W_LenDouble,'pdWorkL         ')
  60      i = 10
              if (O_FilePos.ge.maxrecL) goto 70
              call Buf_Read (O_BlockNum,O_FilePos,i,dWork(1),kbf)
              if (i.le.0) goto 70
              call Buf_Write (S_BlockNum,S_FilePos,i,dWork(1),kbf)
          goto 60
c --- deallocate the temporary work memory
          call MemFree (pdWorkL)
      endif

c --- deallocate the temporary work memory
  70  call Bin_Clo (O_BlockNum,'KEEP')


  999 continue

#if defined(DEBUG)
      call SubEnd('OvrToSvr')
#endif

      return
      end