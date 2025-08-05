      subroutine UpdElemConn (numel,disp,SVRIndex,ififo)
c --- primary function: update connectivity of contact elems
c
c  input arguments:
c     numel         (int,sc,in)             - number of elements (D_NumEl)
c     disp          (dp,ar(D_ULen,*),in)    - global disp vectors
c     SVRIndex      (int,ar(D_NumEl,2),out) - SVR file index vector
c     Kmtrx         (dp,ar(*),in)           - element matrix array
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "in_mem.inc"
#include "cont_info.inc"
c
      integer*8 MakeLongInt
      integer*8 lcSVRL
      integer numel,elord,elem,SVRIndex(numel,2),imast,elData(EL_DIM),
     &        nnod,nr,kelfil(10),eomask,NextSMP,elmget,etyget,ival(3),
     &        kelreqLoc(10),elChar(IELCSZ),nodes(W_NdMax),ifoun,mloc,
     &        ififo(K_nPerFounI,K_nFoun),ityp,i,InqSMP,idummy
      double precision disp(D_ULen),elmass,center(3)
      integer ls(W_NrMax)
      double precision xyzang(6,W_NdMax),uelm(5*W_NrMax),ddummy
      logical lParallel
      
#if defined(DEBUG)
      call SubBeg('UpdElemConn')
#endif

c --- check for parallel
      lParallel = .false.
      if (InqSMP(4).gt.0) lParallel = .true.

c --- loop on contact elements
      imast = 0

100   continue
         if (lParallel) then
             imast = NextSMP()
         else
             imast = imast + 1
         endif 
         if (imast.gt.K_nMast) goto 999
         
         elem = M_ElSet(imast)
         nnod = elmget(elem,elData(1),nodes(1))
         if (nnod.le.0) goto 100
         if (elData(EL_DEAD).eq.1) goto 100
         call conget (elem,ival(1))
         ifoun = abs(ival(3))
         if (ififo(1,ifoun).le.0) goto 100
         if (ififo(2,ifoun).gt.0.and.ififo(1,ifoun).eq.1) goto 100
         if (ififo(13,ifoun).ne.3) goto 100
         if (ififo(15,ifoun).le.0) goto 100
         elord = M_ElemOrd(elem)
         if (elord.le.0) goto 100
c --- update contact element connectivity
         ityp = elData(EL_TYPE)
         i = etyget(ityp,elChar(1))
         if (nnod.lt.W_NdMax)
     &       call vi_zero (nodes(nnod+1),W_NdMax-nnod)

c --- get the SVR file pointer
         lcSVRL = MakeLongInt(SVRIndex(elord,1),SVRIndex(elord,2))

         call DefLsVect (elChar(1),nodes(1),ls(1),nr)

c --- no manipulation of matrices!
         call vi_zero (kelfil(1),10)
         call vi_zero (kelreqLoc(1),10)

c --- initialize element mass and centroid
         elmass = 0.0d0
         call v_zero (center(1),3)
c --- form the appropriate matrices and load vectors
c --- get element nodal coordinates
         call RetXYZ (nnod,nodes(1),xyzang(1,1))
c --- get the element displacements
         call DispInGlobal (nr,elChar(1),nodes(1),ls(1),D_ULen,
     &                      disp(1),A_NumSolVects,uelm(1))
         if (A_Trans.gt.0 .and. A_IterNum.gt.1) 
     &       call SmoothAcce (nr,ls(1),uelm(1))
c --- call for element formation
         eomask = 0
         call ElemMatForm (elem,elChar(1),elData(1),eomask,nodes(1),nr,
     &                     kelreqLoc(1),kelfil(1),xyzang(1,1),uelm(1),
     &                     ddummy,ddummy,elmass,center(1),lcSVRL,
     &                     idummy)
         goto 100
 999  continue

#if defined(DEBUG)
      call SubEnd('UpdElemConn')
#endif

      return
      end
