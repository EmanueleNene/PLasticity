      subroutine CalIniPeni (elem,nnode,ityp,rvr,ival,nSeg,
     &  nPerSeg,nPerFounI,nPerFounDp,nPerMast,nMast,nodes,ififo,coor0,
     &  dispn,isegda,founda,icnda,icon,da)
c --- primary function: calculate initial gap and penetration for 3d beam contact
c
c  input arguments:
c       elem   (int,sc,in)        - element number
c       nnode  (int,sc,in)        - number of node in element
c       nodes  (int,ar(nnod),in)  - array of element node numbers
c       ival   (int,ar(10),in)    - element options
c       ityp   (int,sc,in)        - element type
c       nSeg  (int,sc,in)         - total number of segments
c       nPerSeg  (int,sc,in)      - number of items per segment database entry
c       nPerFounI (int,sc,in)     - number of integer items per foundation
c       fpin   (dp,sc,in)         - pinb factor
c       rvr    (dp,ar(18),in)     - real constant
c       isegda (int,ar(nPerSeg,nSeg),in)
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
c                                  --> none adjacent segmentc
c          isegda(15,iseg)       = segment element number
c                                  < slave element death
c          isegda(16,iseg)       = attached beam element number, if it is
c          isegda(17,iseg)       = attached contact element number, if it is
c          isegda(18,iseg)        = remark
c
c       ififo    (int,ar(nPerFounI),inout)
c                                  - foundation information (integer)
c          ififo(1)             = number of segments
c          ififo(4)             = constrain condition in x-dir.
c          ififo(5)             = constrain condition in y-dir.
c          ififo(6)             = constrain condition in z-dir.
c          ififo(7)             = dimensions
c          ififo(8)             = 1st segment numbering on this foundation
c          ififo(9)             = 2nd segment numbering on this foundation
c          ififo(10)            = lowest nodal numbering on this foundation
c          ififo(12)            = real constant number (ireal)
c          ififo(13)            = contact pair key
c                                 0 - no contact pair
c                                 1 - only slave element is defined
c                                 2 - only contact element is defined
c                                 3 - contact pair is defined
c          ififo(14)            = type of associated contant element type (ityp)
c          ififo(15)            = type of deformable slave surface
c                                -2 - mixed type surface (may wrong)
c                                -1 - rigid slave surface
c                                 0 - none
c                                 1 - deformable slave surface
c                                 2 - deformable slave surface (two pass)
c          ififo(16)            = contact element number who is closest to current slave
c          ififo(17)            = slave element number who is pair of closest contact element
c          ififo(18)            = iteration number of surface adjust
c
c       founda   (dp,ar(nPerFounDp),inout)
c                               - foundation  cell region
c          founda(1)            = minimum value in x-dir.
c          founda(2)            = maximum value in x-dir.
c          founda(3)            = minimum value in y-dir.
c          founda(4)            = maximum value in y-dir.
c          founda(5)            = minimum value in z-dir.
c          founda(6)            = maximum value in z-dir.
c          founda(7)            = first real constant
c          founda(8)            = second real constant
c          founda(9)            = gap or penetration of closest contact element
c          founda(10)           = x component of founda(9)
c          founda(11)           = y component of founda(9)
c          founda(12)           = z component of founda(9)
c          founda(13)           = x component of founda(9,ifoun) (total)
c          founda(14)           = y component of founda(9,ifoun) (total)
c          founda(15)           = z component of founda(9,ifoun) (total)
c          founda(16)           = thickness of founda(9,ifoun)
c          founda(17)           = length of founda(9,ifoun)
c
c       coor0  (dp,ar(*),in)       - initial nodal coordinate coordinate
c       dispn  (dp,ar(*),in)       - corrent nodal displacement
c       iseg   (int,sc,inout)      - contact pair remark
c
#include "impcom.inc"
#include "constants.inc"
c

      integer
     &  elem,nnode,nPerSeg,nSeg,nPerFounI,nPerFounDp,npi,icon,
     &  nodes(nnode),isegda(nPerSeg,nSeg),ififo(nPerFounI),iseg,
     &  elmgct,etyiqr,ityp,ival(3),nPerMast,nMast,icnda(nPerMast,nMast)
c
      double precision CalDist2P,coor0(*),dispn(*),founda(nPerFounDp),
     &  fpin,cradi,tols,da
c
      integer
     &  i,ipi,nocon,ickcel,k,ibeam,ivolu,methi,cnnode(2)

c
      double precision
     &  ZERO,ONE,length,pinb,x,y,z,dajc,pinb2,
     &  xmax,ymax,zmax,xmin,ymin,zmin,
     &  r(3,3),elvol !,gapcri
c
      double precision
     &  xyt(4,7),epsil(7),gapshft,
     &  fcell(6),delxy(3),rvr(W_NumContParam),dvet(11),ebeam,
     &  xa(7),ya(7),za(7)
c
      parameter (ZERO=0.0d0, ONE=1.0d0)
c

#if defined(DEBUG)
      call SubBeg('CalIniPeni')
#endif
      npi = 1
      methi = etyiqr(ityp,-KYOP8)
      ibeam = icnda(2,abs(ival(2)))
c --- zero stress for underlying element
      if (ibeam.gt.0.and.ibeam.ne.elem) then
          ivolu = elmgct (ibeam,dvet(1))
          if (ivolu.ne.0) then
              dvet(6) = 0.d0
              call elmpct (ibeam,dvet(1))
          endif
      endif

c --- slave edge tolerence to prevent slipping off
      if (rvr(25).le.TINY) then
          tols = ONE  ! in % of the segment length
      else
          tols = rvr(25)
      endif

      if (rvr(6).gt.TINY) then
          fpin = rvr(6)
      else
          fpin = 3.d0
          if (ififo(15).ge.1) fpin = 2.d0
          if (abs(rvr(7)).gt.TINY.or.abs(rvr(8)).gt.TINY)
     &        fpin = 4.d0
          if (methi.eq.2.or.methi.eq.4) fpin = 4.d0
          if (etyiqr(ityp,-KYOP1).eq.5.or.etyiqr(ityp,-KYOP1).eq.4)
     &        fpin = fpin/4.d0
      endif

c --- nodal coordinates
c
      do 10 i = 1,nnode
         if (nodes(i).eq.0) goto 10
         xyt(1,i) = dispn(nodes(i)) + coor0(nodes(i))
         xyt(2,i) = dispn(nodes(i)+1) + coor0(nodes(i)+1)
         xyt(3,i) = dispn(nodes(i)+2) + coor0(nodes(i)+2)
  10  continue

c --- calculate element length
      length = CalDist2P(xyt(1,1),xyt(1,2))

c --- get thickness of attached element
c
      cradi = min(abs(rvr(1)),abs(rvr(2)))
      if (cradi.lt.TINY) cradi = 0.125d0*length
c
      ivolu = elmgct (elem,dvet(1))
      ivolu = ivolu

c --- find stiffness
      ebeam = dvet(5)
      dvet(6) = cradi

      call elmpct (elem,dvet(1))

c --- calculate contact pinball range inside which contact may occur
      pinb = fpin*cradi
      pinb = max(pinb,1.1d0*abs(rvr(10)))
      if (rvr(6).lt.ZERO) pinb = -rvr(6)
      pinb = pinb + abs(rvr(1)+rvr(2))
      pinb2 = pinb
      if (pinb2.le.2.d0*length) pinb2 = 2.d0*length

c
      if (ififo(11).lt.0.and.ififo(2).gt.0) goto 03

c --- check whether the current element is in the foundation
      xmax = xyt(1,1)
      ymax = xyt(2,1)
      zmax = xyt(3,1)
      xmin = xyt(1,1)
      ymin = xyt(2,1)
      zmin = xyt(3,1)
      do 02 i = 2,nnode
         if (nodes(i).eq.0) goto 02
         xmax = max (xyt(1,i),xmax)
         xmin = min (xyt(1,i),xmin)
         ymax = max (xyt(2,i),ymax)
         ymin = min (xyt(2,i),ymin)
         zmax = max (xyt(3,i),zmax)
         zmin = min (xyt(3,i),zmin)
   02  continue
      fcell(1) = founda(1) - pinb2
      fcell(2) = founda(2) + pinb2
      fcell(3) = founda(3) - pinb2
      fcell(4) = founda(4) + pinb2
      fcell(5) = founda(5) - pinb2
      fcell(6) = founda(6) + pinb2
      ickcel = 0
      if (xmax .lt. fcell(1)) ickcel = 1
      if (xmin .gt. fcell(2)) ickcel = 1
      if (ymax .lt. fcell(3)) ickcel = 1
      if (ymin .gt. fcell(4)) ickcel = 1
      if (zmax .lt. fcell(5)) ickcel = 1
      if (zmin .gt. fcell(6)) ickcel = 1

 03   continue
c --- initial contact detection
c
c --- get contact normal

      cnnode(1) = nodes(1)
      cnnode(2) = nodes(2)

      do 20 k = 1,nnode
         if (nodes(k).eq.0) goto 20
c
         xa(k) = xyt(1,k)
         ya(k) = xyt(2,k)
         za(k) = xyt(3,k)

  20  continue

      elvol = ZERO
      iseg = 0
      nocon = 2
      do 30 ipi=1,npi

         if (nodes(ipi).le.0) goto 30
c --- coordinates and velocities in initial configuration
         dajc  = length
         if (rvr(2).gt.TINY) then
             elvol = dajc*abs(rvr(2))
         else
             elvol = dajc*0.125d0*length
         endif
         if (ickcel.eq.1) goto 30

         if (ififo(11).lt.0.and.ififo(2).gt.0) goto 30

         if (methi.le.2) then
             gapshft = rvr(10)
         else
             gapshft = 0.0d0
         endif
         call v_zero(r(1,1),9)
c --- compute geometric relations between structure and foundation
         call CalPeni0 (nodes(1),x,y,z,xa(1),ya(1),za(1),
     &        r(1,1),pinb,nSeg,nPerSeg,nPerFounI,nPerFounDp,
     &        ififo(1),coor0(1),dispn(1),isegda(1,1),
     &        founda(1),nocon,iseg,epsil(1),delxy(1),
     &        tols,gapshft,cnnode(1))
         
         if (nocon.eq.2) iseg = 0

c --- prevent spurious contact
         if (nocon.eq.0.and.iseg.gt.0.and.-epsil(3).gt.0.1d0*pinb.and.
     &       ififo(15).ge.1) then
             if (ibeam.eq.abs(isegda(16,iseg))) then
                 epsil(3) = pinb
             endif
         endif

         if ((methi.eq.3.or.methi.eq.4.or.etyiqr(ityp,-KYOP1).eq.6)
     &       .and.abs(epsil(3)).lt.pinb) then
c             gapcri = rvr(5)
c             if (gapcri.eq.0.0d0) then
c                 gapcri = 0.01d0
c             endif
c             if (gapcri.gt.0.0d0) then
c                 gapcri = min(gapcri*length,gapcri*cradi)
c             else
c                 gapcri = - gapcri
c             endif
             if (methi.eq.3.or.methi.eq.4) then
c                 if (epsil(3).le.gapcri) then
c                     epsil(3) = -rvr(10)
c                 else
                     epsil(3) = epsil(3)-rvr(10)
c                 endif
             endif
         endif

         if (nocon.eq.0.and.abs(epsil(3)).gt.pinb) nocon = 2
c --- get the closest contact point information

         if (abs(epsil(3)).lt.pinb.and.iseg.gt.0) then
             if (epsil(3).lt.founda(9)) then
                 founda(9) = epsil(3)
                 ififo(16) = elem
                 ififo(17) = isegda(15,iseg)
                 founda(10) = delxy(1)
                 founda(11) = delxy(2)
                 founda(12) = delxy(3)
                 founda(25) = cradi
                 founda(27) = length
             endif
             if (abs(epsil(3)).gt.abs(founda(32)).and.nocon.eq.0) then
                 if (epsil(3).lt.0.or.etyiqr(ityp,-KYOP1).eq.5.or.
     &               etyiqr(ityp,-KYOP1).eq.4) then
                     founda(32) = epsil(3)
                     founda(33) = elem
                     founda(34) = isegda(15,iseg)
                 endif
             endif
         endif
  30  continue
      founda(16) = founda(16) + length
      founda(17) = founda(17) + cradi
      founda(18) = founda(18) + ebeam
      founda(19) = founda(19) + 1.d0
      founda(20) = founda(20) + pinb
      founda(21) = founda(21) + elvol

      da = HUGE
      if (nocon.eq.0) then
          icnda(3,icon) = iseg
          isegda(19,iseg) = icon
          da = epsil(3)
      endif
         
#if defined(DEBUG)
      call SubEnd('CalIniPeni')
#endif

      return
      end
