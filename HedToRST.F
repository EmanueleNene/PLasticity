      subroutine HedToRST
c --- primary function: put heading information on results file
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "RST_file.inc"
#include "dof_info.inc"

      PTRFTN   MemAlloc
      integer*8  Bin_Get

      integer*8  leng8
      integer    ivect(40)

      pointer (piWorkL,iWork)
      pointer (pdWorkL,dWork)
      integer  iWork(*)
      double precision  dWork(*)

#if defined(DEBUG)
      call SubBeg('HedToRST')
#endif

c --- see if file is being written
      if (R_BlockNum .le. 0) goto 999

c --- see if already done
      leng8 = Bin_Get(R_BlockNum) + R_FilePos
      if (leng8 .gt. 0) goto 999

c --- write standard file header
      call Bin_Head (R_BlockNum,R_Unit,R_FilePos)

c --- write results file header
      call vi_zero (ivect(1),40)
      ivect(1) = R_Unit
      ivect(2) = D_NumNd
      ivect(3) = D_BacLen
      ivect(4) = R_MaxSet
      ivect(5) = D_NumDof
      ivect(6) = D_NumEl
      ivect(7) = D_NumEl
      ivect(8) = 0
      ivect(34) = 0
      ivect(35) = 1

c --- header
      R_FilePos = R_FilePos + 40 + 3
c --- dof index
      call Buf_Write (R_BlockNum,R_FilePos,D_NumDof,
     &                D_DofList(1),D_NumDof)
c --- node index
      ivect(15) = R_FilePos
      call Buf_Write (R_BlockNum,R_FilePos,D_BacLen,
     &                M_BackList(1),D_BacLen)
c --- eorder vector
      if (D_NumEl .gt. 0) then
          ivect(14) = R_FilePos
          call Buf_Write (R_BlockNum,R_FilePos,D_NumEl,
     &                    M_ElemOrd(1),D_NumEl)
      endif

c --- write out integer and dp index records
      piWorkL = MemAlloc(3*R_MaxSet,W_LenInt,'piWorkL         ')
      pdWorkL = MemAlloc(R_MaxSet  ,W_LenDouble, 'pdWorkL         ')
      call vi_zero (iWork(1),2*R_MaxSet)
c --- index vector
      R_ResIndxLoc = R_FilePos
      ivect(11) = R_FilePos
      call Buf_Write (R_BlockNum,R_FilePos,2*R_MaxSet,
     &                iWork(1),2*R_MaxSet)
c --- time vector
      call v_zero (dWork(1),R_MaxSet)
      R_TimIndxLoc = R_FilePos
      ivect(12) = R_FilePos
      call Buf_Write (R_BlockNum,R_FilePos,R_MaxSet*W_IntPerDp,
     &                dWork(1),0)
c --- loadstep/iteration and cumulative iteration
      call vi_zero (iWork(1),3*R_MaxSet)
      R_LdStepLoc = R_FilePos
      ivect(13) = R_FilePos
      call Buf_Write (R_BlockNum,R_FilePos,R_MaxSet*3,
     &                iWork(1),R_MaxSet*3)

      call MemFree (pdWorkL)
      call MemFree (piWorkL)

c --- update the header record
      ivect(10) = R_FilePos
      leng8 = R_FilePos
      call MakeShortInt (leng8,ivect(23),ivect(24))
      leng8 = W_HedLen
      call Buf_Write (R_BlockNum,leng8,40,ivect(1),40)

      call GeomToRST (M_ElemOrd(1))

  999 continue

#if defined(DEBUG)
      call SubEnd('HedToRST')
#endif

      return
      end
