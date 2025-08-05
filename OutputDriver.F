      subroutine OutputDriver (numel,disp,lendisp,rforce,
     &                         ElemIndexL,SVRIndex)
c --- primary function: control the computing of the element output quantities
c
c  input arguments:
c     numel    (int,sc,in)             - number of elements
c     disp     (dp,ar(lendisp,*),in)   - global disp vectors in nodal cs
c     lendisp  (int,sc,in)             - length of disp array
c     SVRIndex(int,ar(numel,2),in)     - index for SVR file
c
c  output arguments:
c     ElemIndexL(intL,ar(numel),inout) - results file location for this element
c
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "RST_file.inc"
#include "dof_info.inc"
#include "cont_info.inc"
#include "solu_info.inc"
c
      integer  GetErr,NextSMP,InqSMP,elmget,etyget
      integer*8  MakeLongInt

      integer elord,elem,lendisp,ityp,mat,nnod,nr,kback,
     & kelfil(10),i,eomask,lenrst,erflag,numel,
     & EleResultIdx(25),elData(EL_DIM),SVRIndex(numel,2)
      
      double precision disp(lendisp,*),rforce(lendisp),
     & elener(10),elvol,dperr(4)
      
      integer*8 lcSVRL,lcerstL,locedxL,EleResultIdxL(25),
     & ElemIndexL(numel)
      
      logical lParallel

      integer  elChar(IELCSZ), ls(W_NrMax),nodes(W_NdMax)
      double precision xyzang(6,W_NdMax),Fnr(2*W_NrMax),uelm(5*W_NrMax)

      double precision Kmtrx(W_NrMax*W_NrMax*4)


#if defined(DEBUG)
      call SubBeg('OutputDriver')
#endif

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

         elem = M_ElemOrd(elord)
         kback = 500
         if (elem.le.0) goto 499

c --- get the element
         nnod = elmget (elem,elData(1),nodes(1))
         if (nnod.le.0) goto 499

         if (nnod.lt.W_NdMax)
     &       call vi_zero (nodes(nnod+1),W_NdMax-nnod)

c --- set element output mask
         eomask = 0
         if (elData(EL_DEAD) .eq. 0) then
             eomask = M_ElemOutMask(elem)
         else
             goto 499
         endif    
         ityp = elData(EL_TYPE)
c --- get the element characteristics
         erflag = etyget(ityp,elChar(1))
         if (elChar(JETYP).eq.W_eSlave) goto 499
         mat = elData(EL_MAT)
c --- get file pointers
         lcSVRL = MakeLongInt(SVRIndex(elord,1),SVRIndex(elord,2))
         call lockSMP (W_LockElemInq)
         lcerstL = R_FilePos
         if (btest(eomask,W_OutElIndx)) then
             if (lParallel) then
                 call CalRSTRecLen(elChar(1),mat,elem,elData(1),
     &                             eomask,lenrst)
                 R_FilePos = R_FilePos + lenrst
            endif
         endif
         call unlockSMP (W_LockElemInq)

c --- leave space for element results header record
         do i = 1,25
            EleResultIdxL(i) = 0
         enddo

         if (btest(eomask,W_OutElIndx)) then
             ElemIndexL(elord) = lcerstL
             locedxL = lcerstL
             lcerstL = lcerstL + 25 + 3
         endif

         kback = 0

 499     continue

         if (kback.eq.999) goto 999
         if (kback.eq.500) goto 500
c --- get the dofs for this matrix
         call DefLsVect (elChar(1),nodes(1),ls(1),nr)

         call vi_zero (kelfil(1),10)
c --- initialize element volume and energies
         elvol = 0.0d0
         call v_init (elener(1),10,TINY)
c --- calculate element results
c --- check if element stresses needed
c --- get element nodal coordinates
         call RetXYZ (nnod,nodes(1),xyzang(1,1))
c --- get the element displacements
         call DispInGlobal (nr,elChar(1),nodes(1),ls(1),
     &                      lendisp,disp(1,1),A_NumSolVects,uelm(1))

         if (A_Trans.gt.0 .and. A_IterNum.gt.1) 
     &       call SmoothAcce (nr,ls(1),uelm(1))

c --- call for element output
         call ElemOutForm (elem,elChar(1),elData(1),eomask,nodes(1),nr,
     &                     ls(1),rforce(1),kelfil(1),xyzang(1,1),
     &                     uelm(1),Kmtrx(1),Fnr(1),elvol,elener(1),
     &                     lcSVRL,EleResultIdxL(1),lcerstL)

         if (GetErr().eq.1) goto 999

         if (A_Trans.gt.0) then
             if (elChar(JETYP).ne.W_eSlave .and. 
     &           elChar(JETYP).ne.W_eMaster) then
c --- accumulate terms for CalcRayleigh quotient frequency estimate
c --- skip if contact elements
                 call CalcRayleigh (kelfil(1),nr,Kmtrx(1),uelm(1),
     &                              elener(3))
             endif
         endif

c --- write out the element index
         if (btest(eomask,W_OutElIndx)) then
             do i = 1,25
                if (EleResultIdxL(i).gt.0) then
                    EleResultIdxL(i) = EleResultIdxL(i) - locedxL
                endif
             enddo
             do i = 1,25
                EleResultIdx(i) = EleResultIdxL(i)
             enddo
             call lockSMP (W_LockRst)
             call Buf_Write (R_BlockNum,locedxL,25,EleResultIdx(1),25)
             call unlockSMP (W_LockRst)
          endif
c --- store the new file pointers (in non-parallel mode)
 400      if (.not.lParallel .and. lcerstL.gt.R_FilePos) then
              if (A_KeyFirstPass.eq.1) then
                  call CalRSTRecLen (elChar(1),mat,elem,elData(1),
     &                               eomask,lenrst)
                  if (lcerstL-R_FilePos.gt.lenrst) then
                      if (InqSMP(0).gt.1) then
c --- kill if is parallel
                          call erhandler ('OutputDriver',5000,3,
     &                                    'RST storage error',
     &                                     dperr(1),' ')
                      endif
                  endif
              endif
              R_FilePos = lcerstL
          endif
c --- process next element
      goto 500

 999  continue


#if defined(DEBUG)
      call SubEnd('OutputDriver')
#endif

      return
      end

