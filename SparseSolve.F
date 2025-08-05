      subroutine SparseSolve (plSpL,uLeng,u)
c --- primary function: finish equation assemble and solve equations
c
c  input arguments:
c   plSpL    (PTR,sc,in)    - pointer of sparse matrix
c   uLeng    (int,sc,in)    - length of u vector
c
c  output arguments:
c   u        (dp,ar(uLeng),out) - the solution
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "dof_info.inc"
#include "in_mem.inc"
#include "cont_info.inc"
#include "DMP_file.inc"
#include "file_sys.inc" 
c
      integer   conchk, foriqr, disiqr, SpSolve,
     &          Bin_Set
      PTRFTN    MemAlloc

      PTRFTN    plSpL
      integer   uLeng, nNod, nForce, nDisp, nEqn, dofgl,
     &          NodeInt, DofInt, NodeExt, i, kerr, kbf,
     &          j, sz1, sz2, sz3
      integer*8 pL2, pL3
      character*(B_PathLen) pname 
      
      double precision  dValue(4), u(uLeng,*)

      pointer (piWorkD,iWorkD)
      pointer (pdWorkD,dWorkD)
      pointer (piWorkF,iWorkF)
      pointer (pdWorkF,dWorkF)
      integer iWorkD(*), iWorkF(*)
      double precision dWorkD(*), dWorkF(*)

      pointer (pdLVL,dLoadVect)
      pointer (pxSpL,SpInfo)
      
      double precision  dLoadVect(*)
      integer*8 SpInfo(*)      

#if defined(DEBUG)
      call SubBeg('SparseSolve')
#endif

#if !defined(ANSYS_MEM)
c --- dump space for contact database to file
      if (K_nMast.gt.0 .or. K_nFoun.gt.0) then
c --- open DMP file
          call GetPath (Z_Name,pname)
          Z_BlockNum = 2
          Z_Unit = 2
          Z_FilePos = 0
          i = Bin_Set(Z_BlockNum,Z_Unit,W_Write,2,pname)
          sz1 = K_nFoun*K_nPerFounDp + D_BacLen*D_NumDof
          j = sz1*W_IntPerDp
          call Buf_Write (Z_BlockNum,Z_FilePos,j,K_FounDp(1),0)
          pL2 = Z_FilePos
          call MemFree (K_PtrFounDp)
          sz2 = K_nSeg*K_nPerSeg + K_nFoun*K_nPerFounI
          call Buf_Write (Z_BlockNum,Z_FilePos,sz2,K_Seg(1),sz2)
          pL3 = Z_FilePos
          call MemFree (K_PtrSeg)
          sz3 = K_nMast*K_nPerMast
          call Buf_Write (Z_BlockNum,Z_FilePos,sz3,K_ContInfo(1),sz3)
          call MemFree (K_PtrContInfo)
      endif
#endif
      
      pxSpL = plSpL
      pdLVL = SpInfo(W_SpPtrLv)
      
      piWorkD = PTRFTNNULL
      pdWorkD = PTRFTNNULL
      piWorkF = PTRFTNNULL
      pdWorkF = PTRFTNNULL

      nNod = D_BacLen
      nEqn = nNod*D_NumDof
      nForce = foriqr(0,DB_MAXDEFINED)
      nDisp  = disiqr(0,DB_MAXDEFINED)
      i = 1
      if (nDisp .gt. 0) i = nDisp
      piWorkD = MemAlloc(i,W_LenInt,'piWorkD         ')
      pdWorkD = MemAlloc(i,W_LenDouble,'pdWorkD         ')
      if (nForce.gt.0) then
          i = nForce
          piWorkF = MemAlloc(i,W_LenInt,'piWorkF         ')
          pdWorkF = MemAlloc(i,W_LenDouble,'pdWorkF         ')
      endif
c --- copy the accumulated NR vector to the sparse storage before adjusting
c      call cv1_v2 (M_ForcNR(1),dLoadVect(1),nEqn,-1.0d0)
      call v1mv2_v3 (M_ForcTot(1),M_ForcNR(1),dLoadVect(1),nEqn)
      
c --- add forces to force vector(made redundant)
c      if (nForce.gt.0) then
c          nForce = 0
c          do NodeInt = 1,nNod
c             NodeExt = M_BackList(NodeInt)
c             do DofInt = 1,D_NumDof
c                call GetF (NodeExt,D_DofList(DofInt),dValue(1))
c                if (dValue(1).ne.0.0d0) then
c                    nForce = nForce + 1
c                    iWorkF(nForce) = (NodeInt-1)*D_NumDof+DofInt
c                    dWorkF(nForce) = dValue(1)
c                endif
c             enddo
c          enddo
c          if (D_NumCp .gt. 0) 
c     &        call AddDofCpI (nForce,iWorkF(1))
c          call SpAddVector(plSpL,nForce,dWorkF(1),iWorkF(1))
c      endif

c --- add constraints
      if (nDisp .gt. 0) then
          nDisp = 0
          do NodeInt = 1,nNod
             NodeExt = M_BackList(NodeInt)
             do DofInt = 1,D_NumDof
                if (conchk(NodeExt,D_DofList(DofInt)).ne.2) then
                    call GetDofNumber (M_BackList(NodeInt),
     &                                 D_DofList(DofInt),dofgl)
                    call GetD(NodeExt,D_DofList(DofInt),dValue(1))
                    dValue(1) = dValue(1) - M_Disp(dofgl)
                    nDisp = nDisp + 1
                    iWorkD(nDisp) = (NodeInt-1)*D_NumDof + DofInt
                    dWorkD(nDisp) = dValue(1)
                endif
             enddo
          enddo
          call SpAddConstrain (plSpL,1,nDisp,iWorkD(1),dWorkD(1))
      endif

      call MemFree (pdWorkF)
      call MemFree (piWorkF)

c --- for dummy DOFs, put diagonal equal 1.0d0, put force = 0.0d0
      call SpDeadDof (plSpL,1,nEqn)

      kerr = SpSolve (plSpL,u(1,1))
      if (kerr .ne. 0) then
          call erhandler ('SparseSolve',5000,3,
     &       ' Error in Sparse solver leg= %i',dble(kerr),' ')
          goto 999
      endif

      if (nDisp .gt. 0) then
c --- recover the solution due to any constraints
          call SpBCAdjust(nEqn,nDisp,iWorkD(1),dWorkD(1),u(1,1))
      endif


      if (D_NumCp .gt. 0) then
c --- enforce CPs
          call AddDofCpDp (uLeng,u(1,1))
      endif

 999  call MemFree (piWorkD)
      call MemFree (pdWorkD)

#if !defined(ANSYS_MEM)
c --- retrive space for contact data
      if (K_nMast.gt.0 .or. K_nFoun.gt.0) then
c --- open DMP file
         Z_FilePos = 0
c --- get dp-memory space for foundation and contact databases
         K_PtrFounDp = MemAlloc(sz1,W_LenDouble,'K_PtrFounDp     ')
         K_PtrNdCoord =K_PtrFounDp + K_nFoun*K_nPerFounDp*W_AddrStrideDp
         j = sz1*W_IntPerDp
         call Buf_Read (Z_BlockNum,Z_FilePos,j,K_FounDp(1),kbf)
c --- get in-memory space for foundation and contact databases
         K_PtrSeg = MemAlloc(sz2,W_LenInt,'K_PtrSeg        ')
         K_PtrFounI  = K_PtrSeg + K_nSeg*K_nPerSeg*W_AddrStrideI
         Z_FilePos = pL2
         call Buf_Read (Z_BlockNum,Z_FilePos,sz2,K_Seg(1),kbf)
         K_PtrContInfo = MemAlloc(sz3,W_LenInt,'K_PtrContInfo   ')
         Z_FilePos = pL3
         call Buf_Read (Z_BlockNum,Z_FilePos,sz3,K_ContInfo(1),kbf)
c --- close and delete the file
         call Bin_Clo (Z_BlockNum,'DELE')
      endif
#endif

#if defined(DEBUG)
      call SubEnd('SparseSolve')
#endif

      return
      end




