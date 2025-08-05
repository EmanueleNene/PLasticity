      subroutine MakeElemOutput (numel,eorder,disp,rforce)
c --- primary function: produce element output
c
c  input arguments:
c     numel    (int,sc,in)           - the number of elements
c     eorder   (int,ar(numel),in)    - the element order vector
c     disp     (dp,ar(*),in)         - global disp vectors
c
c  output arguments:
c     rforce   (dp,ar(*),inout)      - reaction force vector
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "dof_info.inc"
#include "SVR_file.inc"
#include "RST_file.inc"
#include "cont_info.inc"
#if defined(MOIST)
#include "beam_hyg.inc"
#endif
c
      external   OutputDriver
      integer    vi_find,InqSMP,cpget
      PTRFTN     MemAlloc,MemAllocL

      integer*8  Bin_Get, MakeLongInt

      integer*8  loc8,curloc8,leng8,idxlocL,iL
      integer    numel,eorder(numel),i,j,parallelOff,
     &           ivect(80),kbf,idf,ii,jj,kk,nsm,node,
     &           nprocs,kReformSave,nSVR

      double precision  disp(*),rforce(*),condp,dvalue(4)

      pointer (plResuIdxL,ResuIdxL)
      pointer (piCoupL,Coup)

      integer*8  ResuIdxL(*)
      integer    Coup(*)

#if defined(DEBUG)
      call SubBeg('MakeElemOutput')
#endif

      kReformSave = A_KeyReform
c --- set aside space for file indices
      plResuIdxL = MemAlloc(numel,W_LenLongInt,'plResuIdxL      ')

c --- initialize results file element index
      do i = 1,numel
         ResuIdxL(i) = 0
      enddo
      idxlocL = 0
 
      if (btest(R_OutMask,W_OutElIndx)) then
          idxlocL = R_FilePos
          if (R_FilePos .gt. 0) then
c --- otherwise skip index space for later
              R_FilePos = R_FilePos + numel*W_IntPerLong + 3
          endif
      endif

      if (.not.S_InCore) then
          M_PtrSvrIndx = MemAlloc((2+S_nIndx)*numel,
     &                            W_LenInt,'M_PtrSvrIndx    ')
      endif
      
c --- read SVR file
      if (S_BlockNum.gt.0) then
c --- change to read/write
          call Buf_Inf (S_BlockNum,W_ReadAndWrite)
c --- header
          S_FilePos = W_HedLen
          i = 80
          call Buf_Read (S_BlockNum,S_FilePos,i,ivect(1),kbf)
          if (i.lt.80) call vi_zero (ivect(41),40)
c --- get the M_ElemOrd vector
          numel  = ivect(2)
          S_FilePos = MakeLongInt(ivect(33),ivect(43))
          call Buf_Read (S_BlockNum,S_FilePos,numel,eorder(1),kbf)
c --- get the SVR index
          if (.not.S_InCore) then
              S_FilePos = MakeLongInt(ivect(38),ivect(48))
              if (S_FilePos.gt.0) then
                  j = (2+S_nIndx)*numel
                  call Buf_Read(S_BlockNum,S_FilePos,j,M_SvrIndx(1),kbf)
              else
                  call vi_zero (M_SvrIndx(1),(2+S_nIndx)*numel)
              endif
          endif
c --- point to the data records
          S_FilePos = MakeLongInt(ivect(37),ivect(47))
      else
          call vi_zero (M_SvrIndx(1),(2+S_nIndx)*D_NumEl)
      endif

c --- initialize element calculated data
      call v_zero (rforce(1),D_ULen)

      A_CrMax = TINY
      A_StraMax = TINY
      A_MaxEqStra = TINY
      A_DelRealBeam = 0

#if defined(MOIST)
      call InitHygMem()
#endif      
      
c --- get space for the element matrices
      nprocs = InqSMP(0)
c    _______________  this part is to be run in parallel  ______________
      parallelOff = 0
      if (numel.lt.InqSMP(2) .or. nprocs.eq.1) parallelOff = 1
      call InfoSMP (1,parallelOff)
      if (parallelOff.eq.0) then
          call ResumeSMP
          call ForkSMP6 (OutputDriver,numel,disp(1),D_ULen,rforce(1),
     &                   ResuIdxL(1),M_SvrIndx(1))
          call SuspendSMP
      else
          call OutputDriver(numel,disp(1),D_ULen,rforce(1),
     &                      ResuIdxL(1),M_SvrIndx(1))
      endif
c    _______________  end of parallel portion  _______________

#if defined(MOIST)
      call MemFree (H_pDpVArs)
      call MemFree (H_pIntVArs)
#endif          

c --- write out the element index on the results file
      if (idxlocL.gt.0) then
          do i = 1,numel
             if (ResuIdxL(i).gt.0) ResuIdxL(i) = ResuIdxL(i) - idxlocL
          enddo
          i = numel*W_IntPerLong
          call MakeLongIntV (numel,ResuIdxL(1))
          call Buf_Write (R_BlockNum,idxlocL,i,ResuIdxL(1),i)
c --- update the max length written so far
          curloc8 = Bin_Get(R_BlockNum)
          leng8 = curloc8 + R_FilePos
          loc8 = 0
          call Bin_Put (R_BlockNum,loc8)
          idxlocL = W_HedLen
          j = 40
          call Buf_Read (R_BlockNum,idxlocL,j,ivect(1),kbf)
          ivect(10) = leng8
          call MakeShortInt (leng8,ivect(23),ivect(24))
          idxlocL = W_HedLen
          call Buf_Write (R_BlockNum,idxlocL,40,ivect(1),40)
          call Bin_Put (R_BlockNum,curloc8)
      endif
      call MemFree (plResuIdxL)
      if (.not.S_InCore) call MemFree (M_PtrSvrIndx)

c --- calculate the response frequency or eigenvalue
      call CalcTransFreq
      
      A_KeyReform = kReformSave

      do i = 1,D_ULen
c --- adjust reaction forces for applied forces
         if (btest(M_DofBits(i),W_ForcBit)) then
             call InqDof (i,node,idf)
             call GetF (node,idf,dvalue(1))
             rforce(i) = rforce(i) - dvalue(1)
c --- calculate total reaction force
         elseif (btest(M_DofBits(i),W_DispBit)) then
             call InqDof (i,node,idf)
             if (idf.le.3)
     &           A_TotReactCur = A_TotReactCur + abs(rforce(i))
            endif
      enddo
    
c --- sum reaction forces at coupled nodes
      if (D_NumCp.gt.0) then
          piCoupL = MemAlloc(D_LenCp,W_LenInt,'piCoupL         ')
          do 70 i = 1,D_NumCp
             nsm = cpget(i,Coup(1))
             j = Coup(nsm+1)
             idf = vi_find(D_DofList(1),j,D_NumDof)
             if (idf.le.0) goto 70
c --- search for used master node
             do 100 j = 1,nsm
                node = Coup(j)
                node = M_FrwdList(node)
c --- check for used node
                if (node.gt.0) then
                    jj = (node-1)*D_NumDof + idf
                    condp = disp(jj)
                    if (condp.ne.HUGE) goto 110
                endif
 100         continue
c --- none of the nodes are used
             goto 70
c --- node j is a valid master
 110         do 171 ii = 1,nsm
c --- skip master node
                if (ii.eq.j) goto 171
                node = Coup(ii)
                node = M_FrwdList(node)
c --- bypass unused nodes
                if (node.le.0) goto 171
                kk = (node-1)*D_NumDof + idf
                rforce(jj) = rforce(jj) + rforce(kk)
                rforce(kk) = 0.0d0
 171         continue
  70      continue
          call MemFree (piCoupL)
      endif
 
#if defined(DEBUG)
      call SubEnd('MakeElemOutput')
#endif

      return
      end
