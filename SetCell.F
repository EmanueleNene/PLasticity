      subroutine SetCell (nPerSeg,nSeg,nPerFounI,
     &                    nPerFounDp,nFoun,icoor,nodbac,nodfwd,isegda,
     &                    ififo,founda,dispn,icnda,nopilo,
     &                    nPerMast,nMast,numel)
c --- primary function: determine the current cell of foundations
c
c  input arguments:
c     D_BacLen   (int,sc,in)         - number of active nodes
c     D_NumDof   (int,sc,in)         - number of active dof per node
c     nPerSeg    (int,sc,in)         - number of items per segment database entry
c     nSeg    (int,sc,in)            - number of segments on foundation
c     nPerFounI   (int,sc,in)        - number of integer items per foundation
c                                      database entry
c     nPerFounDp   (int,sc,in)       - number of double precision items per
c                                      foundation database entry
c     nFoun    (int,sc,in)           - number of foundations
c     icoor    (dp,ar(*),in)         - initial nodal coordinate
c     nodbac   (int,ar(D_BacLen),out)- internal to exter list
c
c     isegda   (int,ar(nPerSeg,nSeg),in)
c                                    - segment information
c          isegda(1,iseg)        = code of segment type (icode)
c                                = 1   : straight line
c                                = 2   : parabolic arc
c          isegda(2,iseg)        = fondation number
c          isegda(j,iseg)        = 4 nodes (j=3-10) on seg. iseg
c          isegda(j,iseg)        = 4 adjacent segments (j=11-14)
c                                = 2 adjacent segments (j=11,12)
c                                     for seg. iseg
c                                = 0  (icode > 10)
c                                  --> none adjacent segment
c          isegda(15,iseg)       = segment element number
c                                  < slave element death
c          isegda(16,iseg)       = attached beam element number, if it is
c          isegda(17,iseg)       = attached contact element number, if it is
c          isegda(18,iseg)       = remark
c
c  output arguments:
c     dispn     (dp,ar(*),inout)   - corrent nodal displacement
c     ififo    (dp,ar(nPerFounI,nFoun),inout)
c                                  - foundation information (integer)
c          ififo(1,ifoun)       = number of segments
c          ififo(4,ifoun)       = constrain condition in x-dir.
c          ififo(5,ifoun)       = constrain condition in y-dir.
c          ififo(6,ifoun)       = constrain condition in z-dir.
c          ififo(7,ifoun)       = dimensions
c          ififo(8,ifoun)       = 1st segment numbering on this foundation
c          ififo(9,ifoun)       = 2nd segment numbering on this foundation
c          ififo(10,ifoun)      = lowest nodal numbering on this foundation
c          ififo(12,ifoun)      = real constant number (ireal)
c          ififo(13,ifoun)      = contact pair key
c                                 0 - no contact pair
c                                 1 - only slave element is defined
c                                 2 - only contact element is defined
c                                 3 - contact pair is defined
c          ififo(14,ifoun)      = type of associated contant element type (ityp)
c          ififo(15,ifoun)      = type of deformable slave surface
c                                -2 - mixed type surface (may wrong)
c                                -1 - rigid slave surface
c                                 0 - none
c                                 1 - deformable slave surface
c                                 2 - deformable slave surface (two pass)
c          ififo(16,ifoun)      = contact element number who is colsest to current slave
c          ififo(17,ifoun)      = slave element number who is pair of closest contact element
c          ififo(18,ifoun)      = iteration number of surface adjust
c     founda   (dp,ar(nperf,nFoun),inout)
c                                  - foundation  cell region
c          founda(1,ifoun)       = minimum value in x-dir.
c          founda(2,ifoun)       = maximum value in x-dir.
c          founda(3,ifoun)       = minimum value in y-dir.
c          founda(4,ifoun)       = maximum value in y-dir.
c          founda(5,ifoun)       = minimum value in z-dir.
c          founda(6,ifoun)       = maximum value in z-dir.
c          founda(7,ifoun)       = first real constant
c          founda(8,ifoun)       = second real constant
c          founda(9,ifoun)       = gap or penetration of closest contact element
c          founda(10,ifoun)      = x component of founda(9,ifoun)
c          founda(11,ifoun)      = y component of founda(9,ifoun)
c          founda(12,ifoun)      = z component of founda(9,ifoun)
c          founda(13,ifoun)      = x component of founda(9,ifoun) (total)
c          founda(14,ifoun)      = y component of founda(9,ifoun) (total)
c          founda(15,ifoun)      = z component of founda(9,ifoun) (total)
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "in_mem.inc"
c
      integer nSeg,nPerSeg,nPerFounI,nPerFounDp,nFoun,
     &  nodbac(*),ififo(nPerFounI,nFoun),isegda(nPerSeg,nSeg),
     &  nPerMast,nMast,numel,nodfwd(*), icnda(nPerMast,nMast),
     &  nopilo(*)
      double precision icoor(*),dispn(*),founda(nPerFounDp,nFoun)
      integer iseg1,inopil,iropil,iseg2,iddl,inode,knode,
     &  icode,nnod,kfoun,
     &  ln1,ln2,ln3,ln4,j,iseg,ifoun,ln10,ln20,ln30,ln40

      double precision ZERO,x1,x2,x3,x4,dx

      pointer (pdULocui0L,ULocui0)
      double precision  ULocui0(*)

      parameter (ZERO=0.0d0)

#if defined(DEBUG)
      call SubBeg('SetCell')
#endif

c --- get space for global rotations
      pdULocui0L = M_PtrIncI

c --- update contact info
      do 180 kfoun = 1,nFoun
         ifoun = ififo(24,kfoun)
         if (ififo(1,ifoun).le.0) goto 180
         ififo(23,ifoun) = 0
         founda(10,ifoun) = founda(9,ifoun)
         founda(12,ifoun) = founda(11,ifoun)
         founda(29,ifoun) = 1.d0
         if (A_KeyFirstIter.eq.1.and.A_SubStepNum.gt.1) then
c --- first far field
            if (ififo(16,ifoun).eq.0.and.ififo(25,ifoun).ge.0)
     &          ififo(25,ifoun) = -1
c --- first near field
            if (ififo(16,ifoun).gt.0.and.ififo(25,ifoun).le.-1) then
               if (founda(9,ifoun).ge.0.d0) then
                   ififo(25,ifoun) = -2
               endif
            endif
         endif
         if (A_KeyMatUpd.eq.1) then
             if (ififo(16,ifoun).gt.0.and.ififo(25,ifoun).le.-2)
     &           ififo(25,ifoun) = 0
         endif
         founda(11,ifoun) = ZERO
         founda(30,ifoun) = ZERO
         founda(31,ifoun) = ZERO
         founda(32,ifoun) = ZERO
         founda(33,ifoun) = ZERO
         founda(36,ifoun) = 1.d30
         founda(37,ifoun) = ZERO
         founda(28,ifoun) = 1.d30
         founda(35,ifoun) = ZERO !!!
 180  continue
      
c --- check the current cell of foundations
      do 250 kfoun = 1,nFoun
         ifoun = ififo(24,kfoun)
         if (ififo(1,ifoun).le.0) goto 250
         iseg1 = ififo(8,ifoun)
         iseg2 = ififo(9,ifoun)
         inopil = ififo(2,ifoun)
         iropil = ififo(3,ifoun)
         if (iropil.gt.0) iropil = (iropil-1)/D_NumDof + 1
         if (((inopil-1)/D_NumDof+1).gt.0.and.
     &         A_IterNum.ge.2.and.iropil.le.0) then
c --- for rigid translation motion
             do 190 j  = 1,ififo(7,ifoun)
                if (ififo(3+j,ifoun).eq.0) goto 190
                if (ififo(3+j,ifoun).eq.1 .and. A_IterNum.gt.2) goto 190
                dx = ULocui0(inopil+j-1)
                if (abs(dx).ge.HUGE) dx = 0.d0
                founda(2*j-1,ifoun) = founda(2*j-1,ifoun) + dx
                founda(2*j  ,ifoun) = founda(2*j  ,ifoun) + dx
  190        continue
             goto 250
         endif
         do 230 iseg = iseg1,iseg2
            icode = isegda(1,iseg)
            ln10 = isegda(3,iseg)
            ln20 = isegda(4,iseg)
            ln30 = isegda(5,iseg)
            ln40 = isegda(6,iseg)
            do 220 j  = 1,ififo(7,ifoun)
               if (ififo(3+j,ifoun) .eq. 0) goto 220
               if (ififo(3+j,ifoun).eq.1 .and. A_IterNum.gt.2 ) goto 220
               if (iseg .eq. iseg1) then
                   founda(2*j-1,ifoun) =  HUGE
                   founda(2*j  ,ifoun) = -HUGE
               endif
               if (icode.eq.1) then
                   ln1 = ln10 + j -1
                   ln3 = ln30 + j -1
                   x1 = icoor(ln1) + dispn(ln1) + founda(12+j,ifoun)
                   x2 = icoor(ln3) + dispn(ln3) + founda(12+j,ifoun)
                   founda(2*j-1,ifoun) =
     &                   min(founda(2*j-1,ifoun),x1,x2)
                   founda(2*j,ifoun) =
     &                   max(founda(2*j,ifoun),x1,x2)
                   goto 220
               endif
               if (icode.eq.2) then
                   ln1 = ln10 + j -1
                   ln2 = ln20 + j -1
                   ln3 = ln30 + j -1
                   x1 = icoor(ln1) + dispn(ln1) + founda(12+j,ifoun)
                   x2 = icoor(ln2) + dispn(ln2) + founda(12+j,ifoun)
                   x3 = icoor(ln3) + dispn(ln3) + founda(12+j,ifoun)
                   founda(2*j-1,ifoun) =
     &                    min(founda(2*j-1,ifoun),x1,x2,x3)
                   founda(2*j,ifoun) =
     &                    max(founda(2*j,ifoun),x1,x2,x3)
                   goto 220
               endif
               if (icode.eq.7) then
                   ln1 = ln10 + j - 1
                   ln2 = ln20 + j - 1
                   ln3 = ln30 + j - 1
                   ln4 = ln40 + j - 1
                   x1 = icoor(ln1) + dispn(ln1) + founda(12+j,ifoun)
                   x2 = icoor(ln2) + dispn(ln2) + founda(12+j,ifoun)
                   x3 = icoor(ln3) + dispn(ln3) + founda(12+j,ifoun)
                   x4 = icoor(ln4) + dispn(ln4) + founda(12+j,ifoun)
                   founda(2*j-1,ifoun) =
     x                   min(founda(2*j-1,ifoun),x1,x2,x3,x4)
                   founda(2*j,ifoun) =
     x                   max(founda(2*j,ifoun),x1,x2,x3,x4)
                   goto 220
               endif
 220        continue
 230     continue
 250  continue

      call UpdPiloNd (nPerFounI,nFoun,icoor(1),nodbac(1),nopilo(1),
     &                ififo(1,1),dispn(1))
     
#if defined(DEBUG)
      call SubEnd('SetCell')
#endif

      return
      end