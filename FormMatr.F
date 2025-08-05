      subroutine FormMatr (NumEl,eorder,disp,FirstRound)
c --- primary function: compute element matrices and load vectors as requested
c
c  input arguments:
c     NumEl       (int,sc,in)           - the number of elements
c     eorder      (int,ar(NumEl),in)    - the element order vector
c     disp        (dp,ar(*),in)         - past global disp vectors
c     FirstRound  (log,sc,in)           - indicator of elimination round
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "dof_info.inc"
#include "SVR_file.inc"
#include "cont_info.inc"
#include "solver_info.inc"
#if defined(MOIST)
#include "beam_hyg.inc"
#endif
c
     
      external  CompDriver
      integer   GetErr,InqSMP,GetWriteUnit
      integer*8 MakeLongInt
      PTRFTN    MemAlloc, MemAllocL
      
      integer*8 loc8, iL
      integer   NumEl,eorder(NumEl),i,j,ivect(80),kbf,
     &          emsize,nprocs,nSVR,ndel,iott
      integer   parallelOff
      double precision  disp(*),dvect(20)
      logical   jump,FirstRound

      pointer (piElptch,Elptch)
      pointer (pdWork,Work)
      pointer (pdEmatrix,Ematrix)
      pointer (pdSaveForcNR,SaveForcNR)
      double precision Work(*),SaveForcNR(*),Ematrix(*)
      integer Elptch(*)

#if defined(DEBUG)
      call SubBeg('FormMatr')
#endif

c --- reset error status
      C_CutMat = 0
      C_CutEle = 0
      pdEmatrix = PTRFTNNULL
      iott = GetWriteUnit()            
      
c --- contact status key
      if (K_nMast.gt.0) then
          call SetCell (K_nPerSeg,K_nSeg,K_nPerFounI,K_nPerFounDp,
     &                  K_nFoun,K_NdCoord(1),M_BackList(1),
     &                  M_FrwdList(1),K_Seg(1),K_FounI(1),K_FounDp(1),
     &                  M_SoluVects(1),K_ContInfo(1),K_PiloNd(1),
     &                  K_nPerMast,K_nMast,NumEl)
      endif
      
c --- zero out the force vector
      call v_zero (M_ForcTot(1),D_ULen)

      pdSaveForcNR = PTRFTNNULL
      if (A_KeyFirstIter.eq.1 .and. A_BisecNum.gt.0 .and. 
     &    A_KeyFirstPass.eq.0) then
          pdSaveForcNR = MemAlloc(D_ULen,W_LenDouble,'pdSaveForcNR    ')
          call v_move (M_ForcNR(1),SaveForcNR(1),D_ULen)
      endif    
      call v_zero (M_ForcNR(1),D_ULen)

      if (A_KeyStab.gt.0) call v_zero (M_ForcStab(1),D_ULen)

c --- set up solution integer arrays
      call vi_zero (M_DofBits(1),D_ULen)

c --- define the element displacement pointers from the global pointers
      call SetDispPointers

      if (GetErr().eq.1) goto 999
      if (FirstRound) goto 100

c --- SVR file
      if (A_KeyFirstPass.eq.1) then
c --- set aside space for file indices
          M_PtrSvrIndx = MemAlloc((2+S_nIndx)*NumEl,
     &                            W_LenInt,'M_PtrSvrIndx    ')
c --- first pass
          call Buf_Inf (S_BlockNum,W_Write)
          call vi_zero (ivect(1),80)
          ivect(1) = 3
          ivect(2) = NumEl
          ivect(3) = D_NumDof
          ivect(4) = D_ULen
          ivect(5) = D_BacLen
          ivect(6) = D_FwdLen
          S_FilePos = W_HedLen
c --- space for ivect(80)
          S_FilePos = S_FilePos + 80 + 3
c --- space for dvect(20)
          call Bin_Inc (S_BlockNum,S_FilePos,20*W_IntPerDp,0)
          call MakeShortInt (S_FilePos,ivect(31),ivect(41))
          call Buf_Write (S_BlockNum,S_FilePos,D_NumDof,
     &                    D_DofList(1),D_NumDof)
c --- space for M_BackList
          call MakeShortInt (S_FilePos,ivect(32),ivect(42))
          call Buf_Write (S_BlockNum,S_FilePos,
     &                    D_BacLen,M_BackList(1),D_BacLen)
c --- element order
          call MakeShortInt (S_FilePos,ivect(33),ivect(43))
          call Buf_Write (S_BlockNum,S_FilePos,NumEl,eorder(1),NumEl)
c --- space for M_DofBits
          call MakeShortInt (S_FilePos,ivect(49),ivect(50))
          S_FilePos = S_FilePos + D_ULen + 3
c --- space for record index vector
          call MakeShortInt (S_FilePos,ivect(38),ivect(48))
          S_FilePos = S_FilePos + (2+S_nIndx)*NumEl + 3
c --- space for u vectors
          pdWork = MemAlloc(D_ULen,W_LenDouble,'pdWork          ')
          call v_zero (Work(1),D_ULen)
          call MakeShortInt (S_FilePos,ivect(36),ivect(46))
          do i = 1,A_NumSolVects
             j = D_ULen*W_IntPerDp
             call Buf_Write (S_BlockNum,S_FilePos,j,Work(1),0)
          enddo
c --- space for global rotations
          call MakeShortInt (S_FilePos,ivect(51),ivect(52))
          i = 3*D_BacLen*W_IntPerDp
          call Buf_Write (S_BlockNum,S_FilePos,i,Work(1),0)
          call MemFree(pdWork)
c --- start of saved data
          call MakeShortInt (S_FilePos,ivect(37),ivect(47))
c --- write dummy record to extend file
          i = 40
          call Buf_Write (S_BlockNum,S_FilePos,i,ivect(1),i)
c --- write out index vector
          S_FilePos = W_HedLen
          call Buf_Write (S_BlockNum,S_FilePos,80,ivect(1),80)
c --- position for SVR data
          S_FilePos = MakeLongInt(ivect(37),ivect(47))
      else
c --- subsequent passes
          call Buf_Inf (S_BlockNum,W_ReadAndWrite)
          S_FilePos = W_HedLen
          i = 80
          call Buf_Read (S_BlockNum,S_FilePos,i,ivect(1),kbf)
          if (i.lt.80) call vi_zero (ivect(41),40)
          nSVR = i
c --- get the M_ElemOrd vector
          NumEl = ivect(2)
          S_FilePos = MakeLongInt(ivect(33),ivect(43))
          call Buf_Read (S_BlockNum,S_FilePos,NumEl,eorder(1),kbf)
c --- set aside space for file indices
          if (.not.S_InCore) then
              M_PtrSvrIndx = MemAlloc((2+S_nIndx)*NumEl,
     &                                W_LenInt,'M_PtrSvrIndx    ')
c --- get the file index
              S_FilePos = MakeLongInt(ivect(38),ivect(48))
              if (S_FilePos.gt.0) then
                  i = (2+S_nIndx)*NumEl
                  call Buf_Read(S_BlockNum,S_FilePos,i,M_SvrIndx(1),kbf)
              else
                  call vi_zero (M_SvrIndx(1),(2+S_nIndx)*NumEl)
              endif
          endif
c --- start of the element SVR data
          S_FilePos = MakeLongInt(ivect(37),ivect(47))
      endif

c --- pre-define SVR file index (for parallel writing)
c     once its written, you can't change it, that's why
c     we do not update it every iteration
      
      if (A_KeyFirstPass.eq.1) then
          P_Reorder = 0          
          M_PtrSvrD = PTRFTNNULL
          call MakeSvrIndex (NumEl,M_SvrIndx(1),S_FileLengL)
          if (S_InCore) then
              M_PtrSvrD = MemAllocL(S_FileLengL,W_LenInt,
     &                              'M_PtrSvrD       ')
          endif
          if (M_PtrSvrD.gt.PTRFTNNULL) then
              S_InCore = .true.
              S_FilePos = 0
          else
              M_PtrSvrD = PTRFTNNULL
              S_InCore = .false.
          endif
          S_RecStart = S_FilePos
      endif
      
  100 continue
c --- initialize element calculated data (bisection criteria)
      A_CrMax = TINY
      A_StraMax = TINY
      A_MaxEqStra = TINY
      
      call ModiDofBits

      emsize = W_NrMax*W_NrMax
c --- get space for all the element matrices
      emsize = 4*emsize
      nprocs = InqSMP(0)
      pdEmatrix = MemAlloc(emsize*nprocs,W_LenDouble,'pdEmatrix       ')

c --- update connectivity of contact elements
      if (A_KeyFirstPass.eq.1) K_ContStat = 0
c --- for the update on bisection
      if (A_BisecNum.gt.0 .and. A_IterNum.eq.1) K_ContStat = 1
      
      if (K_nMast.gt.0) then
          C_ConnUpd = 1
#if defined(MOIST)
c          if (A_InLineSrch.eq.0) then
          if (A_IterNum.eq.1 .and. A_InLineSrch.eq.0) then
              call v_zero(M_HygStrains(1),H_NumBeams*H_Col*H_Row)
              call InitHygMem()
          endif
#endif
          call UpdContConn (NumEl,disp(1),M_SvrIndx(1),
     &                      pdSaveForcNR,K_FounI(1),K_ContInfo(1),
     &                      K_Seg(1),K_FounDp(1),nprocs,
     &                      FirstRound)
#if defined(MOIST)
c          if (A_InLineSrch.eq.0) then
          if (A_IterNum.eq.1 .and. A_InLineSrch.eq.0) then
              call MemFree (H_pDpVArs)
              call MemFree (H_pIntVArs)
          endif          
#endif          
      endif

      C_ConnUpd = 0
      
c --- remove free end of the fibers for Björn
      if (A_SpecVer.eq.1 .and. FirstRound) then
          piElptch = MemAlloc(NumEl,W_LenInt,'piElptch        ')
          call RemoveEnds (Elptch(1),ndel)
          if (ndel.gt.0) then
              write (iott,2002) ndel
 2002         format('   ->  Database was reduced by',i9,
     &               ' free-end beams')
          endif
          call MemFree(piElptch)
      endif
      
      
      call MemFree (pdSaveForcNR)
      if (FirstRound) goto 999

c    _______________  this part is to be run in parallel  ______________
      parallelOff = 0
      if (NumEl.lt.InqSMP(2) .or. nprocs.eq.1) parallelOff = 1
c --- get maximum matrix size that will be encountered
c --- initialize based on assembly method
c     when A_InLineSrch=1, the call comes from line search
c     and we only need internal forces vector
      if (A_InLineSrch.eq.0) then
          if (M_PtrSparseMtrx.eq.PTRFTNNULL .or. K_ContStat.ne.0) then
              write (iott,*) ' --> START SYMBOLIC ASSEMBLY'
c --- assembly
              call InfoSMP (1,parallelOff)
              call SpDealloc (M_PtrSparseMtrx)
              call ResumeSMP
              call SparseInterface (M_PtrSparseMtrx)
              call SuspendSMP
              P_Reorder = 1
              write (iott,*) ' --> END SYMBOLIC ASSEMBLY'
              if (K_nMast.gt.0) K_ContStat = 0              
          endif
c --- the load vector will be copied, so we zero the matrix only
          call SpZero(M_PtrSparseMtrx,0)
      endif

      call InfoSMP (1,parallelOff)
      if (parallelOff.eq.0) then
          call ResumeSMP
          call ForkSMP5 (CompDriver,
     &                   NumEl,disp(1),M_SvrIndx(1),Ematrix(1),emsize)
          call SuspendSMP
      else
          call CompDriver(NumEl,disp(1),M_SvrIndx(1),Ematrix(1),emsize)
      endif
c    ___________________  end of parallel portion  _____________________

      
      jump = .false.

c --- check for element error
      if (GetErr().eq.1) then
          jump = .true.
c --- print warnings if element jacobians is extremely small
          if (A_IterNum.eq.1 .and. A_SubStepNum.eq.1 .and.
     &        A_LoadStepNum.eq.1) then
              call DistortMessage(3)
          else
c --- set the error flag here in case no error has yet been output
c --- need to force error so we exit cleanly upon abort
             call PutErr(1)
          endif
      endif

      if (jump) then
          if (.not.S_InCore) call MemFree (M_PtrSvrIndx)
          goto 999
      endif

c --- file header double precision data
      call v_zero (dvect(1),20)
      dvect(1) = A_TimeCur
      dvect(2) = A_TimeInc
      dvect(3) = 0.0d0
      dvect(4) = A_TimeBeg
      dvect(5) = A_TimeEnd

c --- finish the SVR file header
      loc8 = S_FilePos
      S_FilePos = W_HedLen
      i = 80
      call Buf_Read (S_BlockNum,S_FilePos,i,ivect(1),kbf)
      if (i.lt.80) call vi_zero (ivect(41),40)
      nSVR = i
      ivect(10) = 1
      ivect(11) = A_LoadStepNum
      ivect(12) = A_CumIterNum
      ivect(13) = A_Trans
      ivect(16) = A_NumSolVects
      ivect(17) = A_SubStepNum
      ivect(18) = A_IterNum
      ivect(21) = A_RowDisp
      ivect(23) = A_RowIncI
      ivect(24) = A_RowIncS
      ivect(26) = A_RowVelo
      ivect(27) = A_RowAcce

      if (A_KeyFirstPass.eq.1) then
          call MakeShortInt (loc8,ivect(40),ivect(39))
      endif
      S_FilePos = W_HedLen
      call Buf_Write (S_BlockNum,S_FilePos,nSVR,ivect(1),nSVR)
      call Buf_Write (S_BlockNum,S_FilePos,20*W_IntPerDp,dvect(1),0)

c --- M_DofBits
      S_FilePos = MakeLongInt(ivect(49),ivect(50))
      if (S_FilePos.gt.0) then
          call Buf_Write (S_BlockNum,S_FilePos,D_ULen,
     &                    M_DofBits(1),D_ULen)
      endif
c --- file record index
      if (A_KeyFirstPass.eq.1) then
          S_FilePos = MakeLongInt(ivect(38),ivect(48))
          i = (2+S_nIndx)*NumEl
          call Buf_Write (S_BlockNum,S_FilePos,i,M_SvrIndx(1),i)
      endif
      if (.not.S_InCore) call MemFree (M_PtrSvrIndx)

  999 continue

      call MemFree (pdEmatrix)

#if defined(DEBUG)
      call SubEnd('FormMatr')
#endif

      return
      end

