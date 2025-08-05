      subroutine SvrToOvr
c --- primary function: copy the SVR file to OVR
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "dof_info.inc"
#include "SVR_file.inc"
#include "OVR_file.inc"
c
      integer  Bin_Set
      integer*8  MakeLongInt
      PTRFTN   MemAlloc, MemAllocL

      integer*8  maxrecL, ihedstL, iL
      integer  kbf, ivect(100), i, j, nSVR
      double precision   dvect(20)
      character*(B_PathLen) pname


      pointer (piWorkL,iWork)
      pointer (pdWorkL,dWork)
      integer  iWork(*)
      double precision  dWork(*)

#if defined(DEBUG)
      call SubBeg('SvrToOvr')
#endif
c --- put in-memory resident structure for no SVR case (not used any more)
      if (S_BlockNum.lt.0) then
          if (M_PtrOvr .eq. PTRFTNNULL) then
              i = A_NumSolVects*D_ULen + 20
              i = i + 3*D_BacLen
              M_PtrOvr = MemAlloc(i,W_LenDouble,'M_PtrOvr        ')
          endif
          call v_zero (dvect(1),20)
          dvect(1) = A_TimeCur
          dvect(2) = A_TimeInc
          dvect(4) = A_TimeBeg
          dvect(5) = A_TimeEnd
          dvect(19) = D_ULen
          dvect(20) = A_NumSolVects
          call v_move (dvect(1),M_Ovr(1),20)
          call v_move (M_SoluVects(1),M_Ovr(21),A_NumSolVects*D_ULen)
          i = A_NumSolVects*D_ULen + 21
          call v_move (M_ERot(1,1),M_Ovr(i),3*D_BacLen)
          goto 999
      endif

c --- open OVR file
      call GetPath (O_Name,pname)
      O_BlockNum = 2
      O_Unit = 17
c --- force file type of OVR file to be same as SVR file
      i = Bin_Set(O_BlockNum,O_Unit,W_Write,B_NumOfPages,pname)

      if (A_TimeCur.eq.A_TimeEnd) then
          call Bin_Clo (O_BlockNum,'DELE')
          goto 999
      endif
      
c --- set aside memory space for element index vector and element SVRs
      piWorkL = MemAlloc((2+S_nIndx)*D_NumEl,
     &                   W_LenInt,'piWorkL         ')

c --- copy SVR to OVR
c --- file header
      S_FilePos = 0
      O_FilePos = 0
      i = 100
      call Buf_Read  (S_BlockNum,S_FilePos,i,ivect(1),kbf)
      call Buf_Write (O_BlockNum,O_FilePos,i,ivect(1),i)
c --- integer header
      i = 80
      call Buf_Read  (S_BlockNum,S_FilePos,i,ivect(1),kbf)
      if (i.lt.80) call vi_zero (ivect(41),40)
      nSVR = i
c --- get the maximum record written
      maxrecL = MakeLongInt(ivect(40),ivect(39))
c --- set the counters to indicate the current solution status
      ivect(11) = A_LoadStepNum
      ivect(17) = A_SubStepNum
      ivect(18) = 1
c --- mark this file as an OVR file
      ivect(20) = 1
c --- mark file to identify that writing is in progress
      ihedstL = O_FilePos
      call Buf_Write (O_BlockNum,O_FilePos,nSVR,ivect(1),nSVR)

c --- dp header
      i = 20
      call Buf_Read (S_BlockNum,S_FilePos,i,dvect(1),kbf)
      dvect(1) = A_TimeCur
      dvect(2) = A_TimeInc
      dvect(4) = A_TimeBeg
      dvect(5) = A_TimeEnd

      call Buf_Write (O_BlockNum,O_FilePos,i,dvect(1),0)

c --- D_NumDof
      call MakeShortInt (O_FilePos,ivect(31),ivect(41))
      call Buf_Write (O_BlockNum,O_FilePos,D_NumDof,
     &                D_DofList(1),D_NumDof)

c --- M_BackList
      call MakeShortInt (O_FilePos,ivect(32),ivect(42))
      call Buf_Write (O_BlockNum,O_FilePos,D_BacLen,
     &                M_BackList(1),D_BacLen)

c --- M_ElemOrd
      call MakeShortInt (O_FilePos,ivect(33),ivect(43))
      call Buf_Write (O_BlockNum,O_FilePos,D_NumEl,M_ElemOrd(1),D_NumEl)

c --- M_DofBits
      call MakeShortInt (O_FilePos,ivect(49),ivect(50))
      call Buf_Write (O_BlockNum,O_FilePos,D_ULen,M_DofBits(1),D_ULen)

c --- record index vector
      S_FilePos = MakeLongInt(ivect(38),ivect(48))
      i = (2+S_nIndx)*D_NumEl
      call Buf_Read  (S_BlockNum,S_FilePos,i,iWork(1),kbf)
      call MakeShortInt (O_FilePos,ivect(38),ivect(48))
      call Buf_Write (O_BlockNum,O_FilePos,i,iWork(1),i)

c --- u vectors: write current values
      call MakeShortInt (O_FilePos,ivect(36),ivect(46))
      j = 1
      do i = 1,A_NumSolVects
         call Buf_Write (O_BlockNum,O_FilePos,
     &                   D_ULen*W_IntPerDp,M_SoluVects(j),0)
         j = j + D_ULen
      enddo
      call MakeShortInt (O_FilePos,ivect(51),ivect(52))
      if (M_PtrERot.ne.PTRFTNNULL) then
          j = 3*D_BacLen*W_IntPerDp
          call Buf_Write (O_BlockNum,O_FilePos,j,M_ERot(1,1),0)
      endif
c --- loop thru all element saved data
      S_FilePos = MakeLongInt(ivect(37),ivect(47))
      call MakeShortInt (O_FilePos,ivect(37),ivect(47))
      if (S_InCore) then
          iL = S_FileLengL
          call Buf_WriteL (O_BlockNum,O_FilePos,iL,M_SvrD(1),1)
      else
          i = 10000 ! same in OvrToSvr
          pdWorkL = MemAlloc(i,W_LenDouble,'pdWorkL         ')
  60      i = 10
               if (S_FilePos.ge.maxrecL) goto 70
               call Buf_Read (S_BlockNum,S_FilePos,i,dWork(1),kbf)
               if (i.le.0) goto 70
               call Buf_Write (O_BlockNum,O_FilePos,i,dWork(1),kbf)
          goto 60
          call MemFree (pdWorkL)
      endif
c --- mark the file to indicate that writing is finished
  70  call MakeShortInt (O_FilePos,ivect(40),ivect(39))
      call MemFree (piWorkL)

c --- close the files
      call Bin_Clo (O_BlockNum,'KEEP')
      
  999 continue

#if defined(DEBUG)
      call SubEnd('SvrToOvr')
#endif

      return
      end