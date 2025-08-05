      subroutine Cont_2 (BacLen,NumDof,nodfwd,nPerSeg,nSeg,nPerFounI,
     &                   nMast,nPerFounDp,nFoun,icoor,
     &                   dispn,isegda,ififo,founda)
c --- primary function: form foundation and segment database
c
c  input arguments:
c     BacLen   (int,sc,in)         - number of active nodes
c     NumDof   (int,sc,in)         - number of active dof per node
c     nodfwd   (int,ar(BacLen),in) - forward list for compressed nodes
c     nPerSeg  (int,sc,in)         - number of items per segment database entry
c     nSeg     (int,sc,in)         - number of segments on foundation
c     nPerMast (int,sc,in)         - number of items per contact element
c     nPerFounI(int,sc,in)         - number of integer items per foundation
c                                    database entry
c     nPerFounDp   (int,sc,in)     - number of double precision items per
c                                     foundation database entry
c     nMast    (int,sc,in)         - number of contact elements
c     nFoun    (int,sc,in)         - number of foundations
c     icoor    (dp,ar(*),in)       - initial nodal coordinate
c     dispn    (dp,ar(*),in)       - corrent nodal displacement
c
c  output arguments:
c     isegda   (int,ar(nPerSeg,nSeg),out)
c                                  - segment information
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

c     ififo    (dp,ar(nPerFounI,nFoun),out)
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
c          ififo(14,ifoun)     = type of associated contant element type(ityp)
c          ififo(15,ifoun)      = type of deformable slave surface
c                                -2 - mixed type surface (may wrong)
c                                -1 - rigid slave surface
c                                 0 - none
c                                 1 - deformable slave surface
c          ififo(16,ifoun)      = contact element number who is closest to
c                                 current slave
c          ififo(17,ifoun)      = slave element number who is pair of
c                                 closest contact element
c          ififo(18,ifoun)      = iteration number of surface adjust
c          ififo(19,ifoun)      = max. node per slave element on each
c                                 slave surface
c
c     founda   (dp,ar(nperf,nFoun),out)
c                                  - foundation  cell region
c          founda(1,ifoun)       = minimum value in x-dir.
c          founda(2,ifoun)       = maximum value in x-dir.
c          founda(3,ifoun)       = minimum value in y-dir.
c          founda(4,ifoun)       = maximum value in y-dir.
c          founda(5,ifoun)       = minimum value in z-dir.
c          founda(6,ifoun)       = maximum value in z-dir.
c          founda(7,ifoun)       = first real constant
c          founda(8,ifoun)       = second real constant
c          founda(9,ifoun)     = gap or penetration of closest contact element
c          founda(10,ifoun)      = x component of founda(9,ifoun) (increment)
c          founda(11,ifoun)      = y component of founda(9,ifoun) (increment)
c          founda(12,ifoun)      = z component of founda(9,ifoun) (increment)
c          founda(13,ifoun)      = x component of founda(9,ifoun) (total)
c          founda(14,ifoun)      = y component of founda(9,ifoun) (total)
c          founda(15,ifoun)      = z component of founda(9,ifoun) (total)
c          founda(16,ifoun)      = thickness of founda(9,ifoun)
c          founda(17,ifoun)      = length of founda(9,ifoun)
c     icnda    (int,ar(nPerMast,nMast),out) - contact element information
c                               = 1 element type (ktyp)
c                               = 2 attached beam element number (ibeam)
c                               = 3 element in contact
c                               = 4 - 7 nodes number (nodes)
c                               = 8 contact element number (ielm)
c                               = 9 slave surface number (ifoun)
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
c
      integer BacLen,NumDof,nPerSeg,nSeg,nPerFounI,
     &  nPerFounDp,nFoun,nMast,etyiqr,rlget,elmget
      integer nodfwd(BacLen),isegda(nPerSeg,nSeg),
     &  ififo(nPerFounI,nFoun)
      double precision founda(nPerFounDp,nFoun),
     &  icoor(BacLen*NumDof),dispn(BacLen*NumDof)
      integer ifoun,i,iseg,iel,ielm,nnod,ityp,icode,inode,iseg2,
     &  ietypz,
     &  ireal,elData(EL_DIM),nodes(W_NdMax),ival(3), itypold
      double precision x1,x2,x3,x4,rvr(W_NumContParam),r00


#if defined(DEBUG)
      call SubBeg('Cont_2')
#endif

c --- initiation of foundation database
      do 100 ifoun = 1,nFoun
         ififo(2,ifoun) = -NumDof + 1
         ififo(3,ifoun) = 0

         ireal = ififo(12,ifoun)

         if (ififo(13,ifoun).eq.2.and.ififo(11,ifoun).lt.0)
     &       ififo(13,ifoun) = 3
         if (ififo(13,ifoun).eq.0) then
             call erhandler ('Cont_2',5000,2,
     &       'Contact pair identified by (set %I) is
     &        missing slave and master elements.',dble(ireal),' ')
         endif

         do i = 1,3
            founda(2*i-1,ifoun) =  HUGE
            founda(2*i  ,ifoun) = -HUGE
         enddo

         if (ifoun.eq.1) then
             ififo(9,ifoun) = 0
         else
             ififo(9,ifoun) = ififo(9,ifoun-1) + ififo(1,ifoun-1)
         endif
         ififo(8,ifoun) = ififo(9,ifoun) + 1
 100  continue


c --- rewrite contact element and segments database
      itypold = 0
      do 400 iel = nMast+1,nMast+nSeg
         ielm = M_ElSet(iel)
         nnod = elmget(ielm,elData(1),nodes(1))
         ityp = elData(EL_TYPE)
         if (ityp .ne. itypold) then
             ietypz = etyiqr(ityp,-JETYP)
             itypold = ityp
         endif
c --- foundation and segment information
         if (ietypz.eq.W_eSlave) then
             call conget (ielm,ival(1))
             ifoun = ival(1)
             icode = elData(EL_SHAPE)
             ififo(7,ifoun) = 3
             nnod = 3

             if (icode.eq.99) nnod = 2 
             if (icode.eq.7) nnod = 4

             iseg = ififo(9,ifoun) + 1
             isegda(15,iseg) = ielm
             isegda(16,iseg) = ival(2)
             if (icode.eq.99) isegda(16,iseg) = 0
             ival(2) = iseg
             call conput (ielm,ival(1))            
             if (ifoun .eq. nFoun) then
                 iseg2 = nSeg + 1
             else
                 iseg2 = ififo(8,ifoun+1)
             endif
             if (iseg .ge. iseg2) then
                 call erhandler ('Cont_2',5000,3,
     &          'Segment numbering is incorrect.',0.d0,' ')
             endif
             ififo(9,ifoun) = iseg
             do 340 inode = 1,nnod
                if (nodes(inode).le.0) then
c                    call erhandler ('Cont_2',5000,3,
c     &              'Node numbering is incorrect.',0.d0,' ')
                    nodes(inode) = -NumDof + 1
                else
                    nodes(inode) = (nodfwd(nodes(inode))-1)*NumDof + 1
                endif
                isegda(2+inode,iseg) = nodes(inode)
 340         continue


             isegda(1,iseg) = icode
             isegda(2,iseg) = ifoun
c
             ififo(19,ifoun) = max(ififo(19,ifoun),nnod)

             if (icode.eq.1) then
                 nnod = 3
                 isegda(3,iseg) = nodes(2)
                 isegda(4,iseg) = -NumDof + 1
                 isegda(5,iseg) = nodes(1)
                 call v_zero (rvr(1),W_NumContParam)
                 i = rlget (ififo(12,ifoun),rvr(1))
                 founda(7,ifoun) = rvr(1)
                 founda(8,ifoun) = rvr(2)
             endif
             if (icode.eq.2) then
                 nnod = 3
                 isegda(3,iseg) = nodes(2)
                 isegda(4,iseg) = nodes(3)
                 isegda(5,iseg) = nodes(1)
                 call v_zero (rvr(1),W_NumContParam)
                 i = rlget (ififo(12,ifoun),rvr(1))
                 founda(7,ifoun) = rvr(1)
                 founda(8,ifoun) = rvr(2)
             endif

             do i = 1,nnod
                if (isegda(2+i,iseg).le.0)
     &              isegda(2+i,iseg)= -NumDof + 1
             enddo

             if (icode.eq.99) then
                 isegda(4,iseg) = nodes(2)
                 isegda(5,iseg) = -numdof + 1
                 ififo(10,ifoun) = iseg
                 ififo(2,ifoun) = nodes(1)
                 ififo(3,ifoun) = nodes(2)
                 if (ififo(3,ifoun).le.0) then
                     ififo(3,ifoun) = 0
                     call erhandler ('Cont_2.F',5000,2,'Missing internal
     &                                node %i.', dble(ielm),' ')
                endif
             endif
            
c --- initialize foundation cell
            if (icode.eq.1) then
                x1 = icoor(nodes(1)) + dispn(nodes(1))
                x2 = icoor(nodes(2)) + dispn(nodes(2))
                r00 = rvr(2)
                founda(1,ifoun) = min(founda(1,ifoun),x1-r00,x2-r00)
                founda(2,ifoun) = max(founda(2,ifoun),x1+r00,x2+r00)
                x1 = icoor(nodes(1)+1) + dispn(nodes(1)+1)
                x2 = icoor(nodes(2)+1) + dispn(nodes(2)+1)
                founda(3,ifoun) = min(founda(3,ifoun),x1-r00,x2-r00)
                founda(4,ifoun) = max(founda(4,ifoun),x1+r00,x2+r00)
                x1 = icoor(nodes(1)+2) + dispn(nodes(1)+2)
                x2 = icoor(nodes(2)+2) + dispn(nodes(2)+2)
                founda(5,ifoun) = min(founda(5,ifoun),x1,x2)
                founda(5,ifoun) = min(founda(5,ifoun),x1-r00,x2-r00)
                founda(6,ifoun) = max(founda(6,ifoun),x1+r00,x2+r00)
                i = rlget (ififo(12,ifoun),rvr(1))
                founda(7,ifoun) = rvr(1)
                founda(8,ifoun) = rvr(2)
            endif
            if (icode.eq.2) then
                r00 = rvr(2)
                x1 = icoor(nodes(1)) + dispn(nodes(1))
                x2 = icoor(nodes(2)) + dispn(nodes(2))
                x3 = icoor(nodes(3)) + dispn(nodes(3))
                founda(1,ifoun) = min(founda(1,ifoun),x1-r00,
     &                                         x2-r00,x3-r00)
                founda(2,ifoun) = max(founda(2,ifoun),x1+r00,
     &                                         x2+r00,x3+r00)
                x1 = icoor(nodes(1)+1) + dispn(nodes(1)+1)
                x2 = icoor(nodes(2)+1) + dispn(nodes(2)+1)
                x3 = icoor(nodes(3)+1) + dispn(nodes(3)+1)
                founda(3,ifoun) = min(founda(3,ifoun),x1-r00,
     &                                         x2-r00,x3-r00)
                founda(4,ifoun) = max(founda(4,ifoun),x1+r00,
     &                                         x2+r00,x3+r00)
                x1 = icoor(nodes(1)+2) + dispn(nodes(1)+2)
                x2 = icoor(nodes(2)+2) + dispn(nodes(2)+2)
                x3 = icoor(nodes(3)+2) + dispn(nodes(3)+2)
                founda(5,ifoun) = min(founda(5,ifoun),x1-r00,
     &                                                x2-r00,x3-r00)
                founda(6,ifoun) = max(founda(6,ifoun),x1+r00,
     &                                                x2+r00,x3+r00)
                call v_zero (rvr(1),W_NumContParam)
                i = rlget (ififo(12,ifoun),rvr(1))
                founda(7,ifoun) = rvr(1)
                founda(8,ifoun) = rvr(2)
            elseif (icode.eq.7) then
                x1 = icoor(nodes(1)) + dispn(nodes(1))
                x2 = icoor(nodes(2)) + dispn(nodes(2))
                x3 = icoor(nodes(3)) + dispn(nodes(3))
                x4 = icoor(nodes(4)) + dispn(nodes(4))
                founda(1,ifoun) = min(founda(1,ifoun),x1,x2,x3,x4)
                founda(2,ifoun) = max(founda(2,ifoun),x1,x2,x3,x4)
                x1 = icoor(nodes(1)+1) + dispn(nodes(1)+1)
                x2 = icoor(nodes(2)+1) + dispn(nodes(2)+1)
                x3 = icoor(nodes(3)+1) + dispn(nodes(3)+1)
                x4 = icoor(nodes(4)+1) + dispn(nodes(4)+1)
                founda(3,ifoun) = min(founda(3,ifoun),x1,x2,x3,x4)
                founda(4,ifoun) = max(founda(4,ifoun),x1,x2,x3,x4)
                x1 = icoor(nodes(1)+2) + dispn(nodes(1)+2)
                x2 = icoor(nodes(2)+2) + dispn(nodes(2)+2)
                x3 = icoor(nodes(3)+2) + dispn(nodes(3)+2)
                x4 = icoor(nodes(4)+2) + dispn(nodes(4)+2)
                founda(5,ifoun) = min(founda(5,ifoun),x1,x2,x3,x4)
                founda(6,ifoun) = max(founda(6,ifoun),x1,x2,x3,x4)
            endif
         endif
 400  continue

#if defined(DEBUG)
      call SubEnd('Cont_2')
#endif
      return
      end
      