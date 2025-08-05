      subroutine CompDriver (numel,disp,SVRIndex,MtrxLoc,emsize)
c --- primary function: control the computing of the element matrices and load vectors
c
c  input arguments:
c     numel        (int,sc,in)             - number of elements (D_NumEl)
c     disp         (dp,ar(D_ULen,*),in)    - global disp vectors in nodal cs
c     SVRIndex     (int,ar(D_NumEl,2),out) - SVR file index vector
c     MtrxLoc      (dp,ar(*),in)           - element matrix array
c     emsize       (int,sc,in)             - matrix array space per processor
c
c  referenced variables:
c     ls       (int,ar(*),local)           - dof list for the element
c     elChar   (int,ar(IELCSZ),loc)        - array of element characteristics
c     elData   (int,ar(EL_DIM),local)      - array of element data from the database
c
c  internal variables:
c     kelfil   (int,ar(10))        - keys indicating matrices read/written
c     elvol    (dp,sc)             - element mass
c     center   (dp,ar(3))          - element centroid location
c     eomask   (int,sc,in)         - bit pattern for element output
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "SVR_file.inc"
#include "conv_stat.inc"
#include "in_mem.inc"
#include "cont_info.inc"
c
      integer   GetErr,NextSMP,iProcSMP,InqSMP,
     &          elmget,etyget
      integer*8 MakeLongInt

      integer numel,elord,elem,ityp,jtyp,nnod,nr,kelfil(10),i,
     &        eomask,kback,emsize,mloc,mat,erflag,kelreqLoc(10),
     &        kReform,SVRIndex(numel,2),lensav,elData(EL_DIM)
      double precision  disp(*),elvol,center(3), dperr(4)
      integer*8  lcSVRL
      logical    lParallel, IfFixed

      integer          ls(W_NrMax),nodes(W_NdMax)
      double precision xyzang(6,W_NdMax),MtrxLoc(*),Fnr(2*W_NrMax),
     &                 uelm(5*W_NrMax),uelmSave(2*W_NrMax)

      pointer (pdzsL,Kmtrx)

      double precision Kmtrx(emsize)

      integer elChar(IELCSZ), ndUpdate


#if defined(DEBUG)
      call SubBeg('CompDriver')
#endif

c --- for line search, only update the LV if requested
      if (A_InLineSrch.gt.0) kReform = 0
c --- otherwise always treat as reforming
      if (A_InLineSrch.eq.0) kReform = 1

c --- get location for matrices
      mloc = emsize*iProcSMP() + 1

c --- default is to use local arrays (not dynamically allocated ones)
c --- this is one exeption though
      pdzsL = loc(MtrxLoc(mloc))

c --- check for parallel
      lParallel = .false.
      if (InqSMP(4).gt.0) lParallel = .true.

c --- loop on elements
      elord = 0

 500  continue

         if (lParallel) then
             elord = NextSMP()
         else
             elord = elord + 1
         endif

         kback = 999

         if (elord.gt.D_NumEl) goto 499
         kback = 500

         elem = M_ElemOrd(elord)
c --- get the element
         nnod = elmget(elem,elData(1),nodes(1))
         if (nnod.le.0) goto 499
         if (nnod.lt.W_NdMax)
     &       call vi_zero (nodes(nnod+1),W_NdMax-nnod)
                          
c --- check for dead elements
         if (elData(EL_DEAD).eq.1) goto 499
             
         ityp = elData(EL_TYPE)
c --- get the element characteristics
         erflag = etyget(ityp,elChar(1))
         jtyp = elChar(JETYP)
         if (jtyp.eq.W_eSlave) goto 499

         mat = elData(EL_MAT)         
         call vi_move (A_MatrReq(1),kelreqLoc(1),10)         
c --- get the SVR pointers
         lcSVRL = MakeLongInt(SVRIndex(elord,1),SVRIndex(elord,2))
         if (A_KeyFirstPass.eq.1) then
             if (lParallel) then
                 call CalSVRRecLen (elChar(1),elData(1),mat,lensav)
                 call lockSMP (W_LockSvr)
                 S_FilePos = S_FilePos + lensav
                 call unlockSMP (W_LockSvr)
             endif
         endif

         kback = 0

 499     continue

         if (kback.eq.500) goto 500
         if (kback.eq.999) goto 999

c --- get the DOFs for this matrix
         call DefLsVect (elChar(1),nodes(1),ls(1),nr)

c --- skip forming stiffness matrix during line search
         if (A_InLineSrch.eq.1) then
             kelreqLoc(1) = 0
             kelreqLoc(2) = 0
             kelreqLoc(3) = 0
             if (A_Trans.gt.0) then
c --- transient always needs mass
                 kelreqLoc(2) = 1
                 kelreqLoc(3) = 1
             endif
             if (A_Gravity.eq.1) then
c --- gravity needs mass
                 kelreqLoc(2) = 1
             endif

c --- check for any imposed disps on this element
             if (IfFixed(nr,ls)) kelreqLoc(1) = 1
         endif

         call vi_zero (kelfil(1),10)
c --- initialize element mass and centroid
         elvol = 0.0d0
         call v_zero (center(1),3)

c --- form the appropriate matrices and load vectors
c --- get element nodal coordinates
         call RetXYZ (nnod,nodes(1),xyzang(1,1))
c --- get the element displacements
         call DispInGlobal (nr,elChar(1),nodes(1),ls(1),D_ULen,
     &                      disp(1),A_NumSolVects,uelm(1))

c --- save unmodified velocities and acelerations for DynMatAdj
         if (A_Trans.gt.0 .and. A_IterNum.gt.1) then
             call v_move (uelm(nr*(A_RowVelo-1)+1),uelmSave(1   ),nr)
             call v_move (uelm(nr*(A_RowAcce-1)+1),uelmSave(nr+1),nr)
             call SmoothAcce (nr,ls(1),uelm(1))
         endif         
     
c --- call for element formation
         eomask = 0
         call ElemMatForm (elem,elChar(1),elData(1),eomask,nodes(1),nr,
     &                     kelreqLoc(1),kelfil(1),xyzang(1,1),uelm(1),
     &                     Kmtrx(1),Fnr(1),elvol,center(1),lcSVRL,
     &                     ndUpdate)

         if (GetErr().eq.1 .or. C_CutMat.gt.0 .or. C_CutEle.gt.0) 
     &       goto 999

c --- check if the nodes were updated (for contact elements)
         if (ndUpdate.gt.0) then
             if (kReform.eq.1) then
c --- stop the program since it should not happen here if the matrix is formed
                 call erhandler('CompDriver',5000,3,
     &                          'DOF vector is corrupted',
     &                           dperr(1),' ')
                 call StopAll (-1)
             else
c --- update LS vector if in line search
                 call DefLsVect (elChar(1),nodes(1),ls(1),nr)
             endif
         endif
c --- restore unmodified velocities and acelerations for DynMatAdj
         if (A_Trans.gt.0 .and. A_IterNum.gt.1) then
             call v_move (uelmSave(1   ),uelm(nr*(A_RowVelo-1)+1),nr)
             call v_move (uelmSave(nr+1),uelm(nr*(A_RowAcce-1)+1),nr)
         endif
c --- add NR force, some other minor jobs
         call PrepElem (elem,kelfil(1),nr,elData(1),elChar(1),
     &                  ls(1),Kmtrx(1),Fnr(1),uelm(1),D_ULen,
     &                  M_ForcTot(1),M_ForcNR(1),M_ForcStab(1),kReform)
         
c --- assemble matrices for solution
c     when A_InLineSrch=1, the call comes from line search
c     and we only need internal forces vector
         if (kReform.eq.1) then
c --- define mat and element number: passed down to solver
c --- perform the required assembly of the matrices and load vectors
             call AssemKCL (M_PtrSparseMtrx,nr,kelfil(1),ls(1),Kmtrx(1))
         endif
c --- store the new file pointers (in non-parallel mode)
         if (.not.lParallel) then
             if (A_KeyFirstPass.eq.1) then
                 call CalSVRRecLen (elChar(1),elData(1),mat,lensav)
                 if (lcSVRL.ne.0 .and. A_KeySvrUpd.gt.0 .and.
     &               lcSVRL-S_FilePos.ne.lensav) then
                     if (InqSMP(0).gt.1) then
c --- kill if parallel
                         call erhandler('CompDriver',5000,3,
     &                                  'SVR storage error',
     &                                   dperr(1),' ')
                     endif
                 endif
             endif
             if (lcSVRL.gt.0) S_FilePos = lcSVRL
         endif
c --- process next element
      goto 500

 999  continue

#if defined(DEBUG)
      call SubEnd('CompDriver')
#endif

      return
      end

