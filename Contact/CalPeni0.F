      subroutine CalPeni0 (cnodes,x0,y0,z0,xa,ya,za,
     &                     r,dmax,nSeg,nPerSeg,nPerFounI,nPerFounDp,
     &                     ififo,coor0,dispn,isegda,founda,nocon,iseg,
     &                     epsil,delxy,tols,gapshft,cnnode)
c --- primary function: compute the relations between an integration point
c                       of the contact element and the corresponding foundation
c
c  input arguments:
c      cnodes (int,ar(2),in) - contact node list
c      xa,ya,za (dp,ar(3),in) - coordinates of contact nodes
c      dmax     (dp,sc,in)  - maximum distance for a valid contact
c                                (if dabs(d)>dmax, false contact is assumed)
c      nSeg    (int,sc,in) - total number of segments defining the foundations
c      nPerSeg    (int,sc,in) - number of items per segment database entry
c      ififo    (dp,ar(6),in)
c                                  - foundation information (integer)
c          ififo(1)       = number of segments
c          ififo(4)       = constrain condition in x-dir.
c          ififo(5)       = constrain condition in y-dir.
c          ififo(6)       = constrain condition in z-dir.
c          ififo(7)       = dimensions
c          ififo(8)       = 1st segment numbering on this foundation
c          ififo(9)       = 2nd segment numbering on this foundation
c          ififo(10)      = lowest nodal numbering on this foundation
c          ififo(12)      = real constant number (ireal)
c          ififo(13)      = contact pair key
c                            0 - no contact pair
c                            1 - only slave element is defined
c                            2 - only contact element is defined
c                            3 - contact pair is defined
c          ififo(14)      = type of associated contant element type
c          ififo(15)      = type of deformable slave surface
c                           -2 - mixed type surface (may wrong)
c                           -1 - rigid slave surface
c                            0 - none
c                            1 - deformable slave surface
c                            2 - deformable slave surface (two pass)
c          ififo(16) = contact element number who is closest to current slave
c          ififo(17) = slave element number who is pair of closest
c                       contact element
c          ififo(18) = iteration number of surface adjust

c      coor0       (dp,ar(*),in)- initial nodal coordinate coordinate
c      dispn       (dp,ar(*),in)     - corrent nodal displacement
c      isegda   (int,ar(nPerSeg,nSeg),in)
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
c
c      founda   (dp,ar(nPerFounDp),in)
c                                  - foundation  cell region
c          founda(1)        = minimum value in x-dir.
c          founda(2)        = maximum value in x-dir.
c          founda(3)        = minimum value in y-dir.
c          founda(4)        = maximum value in y-dir.
c          founda(5)        = minimum value in z-dir.
c          founda(6)        = maximum value in z-dir.
c          founda(7)        = first real constant
c          founda(8)        = second real constant
c          founda(9)        = gap or penetration of closest contact element
c          founda(10)       = x component of founda(9)
c          founda(11)       = y component of founda(9)
c          founda(12)       = z component of founda(9)
c          founda(13)       = x component of founda(9) (total)
c          founda(14)       = y component of founda(9) (total)
c          founda(15)       = z component of founda(9) (total)

c       fcell   (dp,ar(4),in) - foundation  cell region
c          fcell(1)         = minimum value in x-dir.
c          fcell(2)         = maximum value in x-dir.
c          fcell(3)         = minimum value in y-dir.
c          fcell(4)         = maximum value in y-dir.
c
c  output arguments:
c      x0 ,y0 ,z0  (dp,sc,out)  - global coordinates of the integration point
c      r(3,3)   (dp,ar(3,3),out)
c                         - components of the normal vector to the structure
c                                only for 3-d case
c       epsil    (dp,ar(7),out)
c        epsil(3)
c                d  = distance between the integration point
c                     and the foundation
c                     ( d < 0 : contact    )
c                     ( d > 0 : no contact )
c       nocon   (int,sc,out)
c                   = 0  : contact
c                   = 2  : no contact
c       iseg  (int,sc,inout)    - number of the foundation segment
c                                 with which contact is established
c                                 (=0 if no contact with the foundation)
c       delxy    (dp,ar(3),out) - distance between integration point
c                                 to slave node
c         delxy(1)          - x component
c         delxy(2)          - y component
c         delxy(3)          - z component
c
#include "impcom.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer nSeg,nPerSeg,nPerFounI,nPerFounDp,nocon,iseg,
     &        isegda(nPerSeg,nSeg),ififo(nPerFounI),cnodes(7)
      double precision x0,y0,z0,xa(7),ya(7),za(7),dmax,
     &                 epsil(7),coor0(*),dispn(*),r(3,3),delxy(3),
     &                  founda(nPerFounDp),ptols,ro(3,3)

      integer is,i1,kddlx,kddly,kddlz,iseg1,iseg2,is0,
     &        inode,nnode,n1,n2,cnnode(2),contactNode(2),tnodes(7),
     &        iseg0,iseg00,ipass1,iconnect(2),locu,locv,i2,numRoot,
     &        debug
      double precision ZERO,HALF,d0,xb(7),yb(7),zb(7),t(7),tols,
     &                 gapshft,fifoa(7),fifob(7),ddum,
     &                 ta(3),tb(3),ac(3),bc(3)

      parameter (ZERO=0.0d0, HALF=0.5d0)

#if defined(DEBUG)
      call SubBeg('CalPeni0')
#endif

      debug = 0
      
      call v_zero (epsil(1),7)
      call v_zero (ro(1,1),9)
      epsil(3) = dmax
      is0 = 0
c --- degenerated case
      locu = 0
      nnode = 7
      do inode = 1,nnode
         if (cnodes(inode).le.0) then
             n1 = 0
             if (inode.eq.3) then
                 n1 = 1
                 n2 = 2
             endif
             if (inode.eq.4) then
                 n1 = 1
                 n2 = 1
                 locu = -1
             endif
             if (inode.eq.5) then
                 n1 = 2
                 n2 = 2
                 if (locu.eq.-1) then
                     locu = 2
                 else
                     locu = 1
                 endif
             endif
             if (inode.eq.6) then
                 n1 = 1
                 n2 = 4
             endif
             if (inode.eq.7) then
                 n1 = 2
                 n2 = 5
             endif

             if (n1.gt.0) then
                 xa(inode) = HALF*(xa(n1)+xa(n2))
                 ya(inode) = HALF*(ya(n1)+ya(n2))
                 za(inode) = HALF*(za(n1)+za(n2))
             else
                 xa(inode) = ZERO
                 ya(inode) = ZERO
                 za(inode) = ZERO
             endif
          endif
      enddo

      d0 = dmax

c --- contact node
      call vi_zero (contactNode(1),2)
      do i2 =1,2
         if (cnnode(i2).gt.0) then
             contactNode(i2) = cnnode(i2)
             cnnode(i2) = (cnnode(i2)-1)/D_NumDof + 1
             cnnode(i2) = M_BackList(cnnode(i2))
         endif
      enddo
      ipass1 = 0
      iseg0 = 0
      iseg00 = 0
      iseg1 = ififo(8)
      iseg2 = ififo(9)
      do is = iseg1,iseg2
         isegda(18,is) = 0
         if (isegda(15,is).lt.0) isegda(18,is) = 1
      enddo
      is = iseg1 -1
      is0 = iseg

      if (iseg.gt.0) then
          if (isegda(18,iseg).eq.0) goto 100
      endif

      is0 = iseg1

  050 continue
      iseg00 = 0
      ipass1 = 0
      do is = iseg1,iseg2
         if (isegda(18,is).ne.0) cycle
         iseg = is
         goto 100
      enddo

      iseg = is0
      nocon = 2
      epsil(3) = dmax
      goto 999

 100  continue
      iseg00 = iseg0
      iseg0 = iseg
      isegda(18,iseg) = 1
      nnode = 7

      call v_zero (t(1),nnode)
      call vi_zero (iconnect(1),2)
      locv = 0

      tnodes(1) = isegda(3,iseg)
      tnodes(2) = isegda(5,iseg)
      tnodes(3) = isegda(4,iseg)
      tnodes(4) = isegda(6,iseg)
      tnodes(5) = isegda(7,iseg)
      tnodes(6) = isegda(8,iseg)
      tnodes(7) = isegda(9,iseg)

      do inode = 1,nnode
         i1 = tnodes(inode)
         if (i1.le.0) then
            n1 = 0
            if (inode.eq.3) then
                n1 = 1
                n2 = 2
            endif
            if (inode.eq.4) then
                n1 = 1
                n2 = 1
                locv = -1
            endif
            if (inode.eq.5) then
               n1 = 2
               n2 = 2
               if (locv.eq.-1) then
                   locv = 2
               else
                   locv = 1
               endif
            endif
            if (inode.eq.6) then
                n1 = 1
                n2 = 4
            endif
            if (inode.eq.7) then
                n1 = 2
                n2 = 5
            endif
            if (n1.gt.0) then
                xb(inode) = HALF*(xb(n1)+xb(n2))
                yb(inode) = HALF*(yb(n1)+yb(n2))
                zb(inode) = HALF*(zb(n1)+zb(n2))
           else
                xb(inode) = ZERO
                yb(inode) = ZERO
                zb(inode) = ZERO
            endif
            cycle
         endif

c --- contact node can not contact its own elements
         if (inode.le.2) then
            do i2  = 1,2
               if (contactNode(i2).gt.0) then
                  if (i1.eq.contactNode(i2)) then
                     if (ififo(15).eq.3) then
                        nocon = 2
                        epsil(3) = d0
                        goto 200
                     endif
                     iconnect(i2) = 1
                  endif
               endif
            enddo
         endif
         kddlx = i1
         kddly = kddlx + 1
         kddlz = kddlx + 2
         xb(inode) = dispn(kddlx) + coor0(kddlx) + founda(13)
         yb(inode) = dispn(kddly) + coor0(kddly) + founda(14)
         zb(inode) = dispn(kddlz) + coor0(kddlz) + founda(15)
      enddo

      d0 = epsil(3)
      ddum = HUGE
c --- reduce slip off tolerance for the initial search
      ptols = 0.8d0*tols
c --- beam segment
      call BeamToBeam (founda(7),x0,y0,z0,xa(1),ya(1),za(1),
     &                 t(1),t(1),t(1),xb(1),yb(1),zb(1),t(1),t(1),t(1),
     &                 r(1,1),ro(1,1),epsil(1),fifoa(1),fifob(1),nocon,
     &                 dmax,delxy(1),isegda(11,iseg),ptols,gapshft,
     &                 ipass1,t(1),t(1),t(1),locu,locv,numRoot,
     &                 ddum,ddum,ddum,ta(1),tb(1),ac(1),bc(1),
     &                 debug)
 200  continue

c --- contact state
      if (abs(epsil(3)).gt.dmax) then
          epsil(3) = dmax
          if (nocon.eq.0) nocon = 2
      endif
      if (nocon.eq.0) goto 999
      if (epsil(3).lt.d0) then
          is0 = iseg
      else
          epsil(3) = d0
      endif
      iseg1 = is + 1
      if (nocon.eq.2) goto 050
      if (nocon.eq.-1) iseg = isegda(11,iseg)
      if (nocon.eq. 1) iseg = isegda(12,iseg)

      if (iseg.eq.0) goto 050
      if (isegda(18,iseg).eq.0) then
          ipass1 = 0
          goto 100
      endif
      if (iseg.eq.iseg00) then
          ipass1 = ipass1 + 1
          if (ipass1.lt.5) goto 100
      endif
      goto 050

 999  continue
      if (nocon.eq.0) then
          do i2  = 1,2
             if (iconnect(i2).eq.1) cnnode(i2) = -cnnode(i2)
          enddo
      endif

#if defined(DEBUG)
      call SubEnd('CalPeni0')
#endif

      return
      end
