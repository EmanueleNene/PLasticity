      subroutine CalIniPeni4 (elem,ityp,rvr,ival,nSeg,ifoun,
     &                        nPerSeg,nPerFounI,nPerFounDp,nPerMast,
     &                        nMast,nodes,ififo,coor0,
     &                        dispn,isegda,founda,icnda,icon,da)
c --- primary function: calculate initial gap and penetration for beam-surface
c
c  input arguments:
c       elem   (int,sc,in)         - element number
c       nodes  (int,ar(nnod),in)   - array of element node numbers
c       ival   (int,ar(10),in)     - element options
c       ityp   (int,sc,in)         - element type
c       nSeg  (int,sc,in)          - total number of segments
c       ifoun  (int,sc,in)         - foundation number
c       nPerSeg  (int,sc,in)       - number of items per segment database entry
c       nPerFounI (int,sc,in)      - number of integer items per foundation
c       rvr    (dp,ar(18),in)      - real constant
c       isegda (int,ar(nPerSeg,nSeg),in)
c                                  - segment information 
c          isegda(1,iseg)        = code of segment type (icode)
c                                = 1   : straight line
c                                = 2   : parabolic arc
c                                = 7   : quadrilateral
c
c          isegda(2,iseg)        = fondation number           
c          isegda(j,iseg)        = 4 nodes (j=3-10) on seg. iseg
c          isegda(j,iseg)        = 4 adjacent segments (j=11-14)
c                                    for seg. iseg    (3-d triangle)
c                                = 2 adjacent segments (j=11,12)
c                                     for seg. iseg    (2-d)
c                                = 0  (icode > 10)
c                                  --> none adjacent segment
c          isegda(15,iseg)       = segment element number
c                                  < slave element death
c          isegda(16,iseg)       = attached beam element number, if it is
c          isegda(17,iseg)       = attached contact element number, if it is
c          isegda(18,iseg)        = remark
c
c       ififo    (int,ar(nPerFounI),inout)
c                                  - foundation information (integer)
c          ififo(1)             = number of segments
c          ififo(2)             = pilot nodal number
c          ififo(3)             = pilot nodal number on rotation
c          ififo(4)             = constrain condition in x-dir.
c          ififo(5)             = constrain condition in y-dir.
c          ififo(6)             = constrain condition in z-dir.
c          ififo(7)             = dimensions
c          ififo(8)             = 1st segment numbering on this foundation
c          ififo(9)             = 2nd segment numbering on this foundation
c          ififo(10)            = lowest nodal numbering on this foundation
c          ififo(11)            = 0 pilot node does not connect to other element
c                               =   element number that pilot node connects to
c          ififo(12)            = real constant number
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
c                                 2 - deformable slave surface (TWO pass)
c          ififo(16)            = contact element number who is colsest to current slave
c          ififo(17)            = slave element number who is pair of closest contact element
c          ififo(18)            = iteration number of surface adjust

c       founda   (dp,ar(nPerFounDp),inout)
c                                  - foundation  cell region
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
c
c  output arguments:
c       iseg   (int,sc,inout)      - contact pair remark
c       thick  (dp,sc,out)         - thickness of attached element (dummy)
c
#include "impcom.inc"
#include "constants.inc"
c
      integer etyiqr,elmgct,rlget,elmget

      integer
     &  elem,nPerSeg,nSeg,nPerFounI,nPerFounDp,ityp,ifoun,icon,
     &  nodes(*),isegda(nPerSeg,nSeg),ififo(nPerFounI),iseg,ival(3),
     &  nPerMast,nMast,icnda(nPerMast,nMast)

      double precision CalDist2P,coor0(*),dispn(*),founda(nPerFounDp),
     &  rvr(W_NumContParam),rvrb(W_NumContParam),thick,crad,da

      integer
     &  i,ipi,nocon,ickcel,ibeam,npi,ivolu,methi,hi,
     &  nnodr,cnnode,irealb,elDatab(EL_DIM),nodesb(W_NdMax)
    
      double precision
     &  ZERO,length,pinb,x,y,dajc,elvol,pinb2,
     &  xmax,ymax,xmin,ymin,fpin,dvet(11)

      double precision
     &  xyt(3,3),epsil(8),weight,fcell(6),delxy(3),tols,esolid,
     &  r(3,3),gapshft

      double precision
     &  z,zmax,zmin

      parameter (ZERO=0.0d0)
      
      nnodr = 3
      if (nodes(3).eq.0.or.nodes(2).eq.nodes(3)) nnodr = 2
      npi  = abs(ival(1))
      methi = etyiqr(ityp,-KYOP8)
      hi = 1
      if (etyiqr(ityp,-KYO12).eq.1) hi = 0
      
      ibeam = icnda(2,abs(ival(2)))

c --- shift to account for beam radius
      crad = ZERO

c --- zero stress for underlying element
      if (ibeam.gt.0.and.ibeam.ne.elem) then
          ivolu = elmgct (ibeam,dvet(1))
          if (ivolu.ne.0) then
              dvet(6) = 0.d0
              call elmpct (ibeam,dvet(1))
          endif
          if (etyiqr(ityp,-KYO11).eq.1) then
c --- get attached beam element radius
              i = elmget(ibeam,elDatab(1),nodesb(1))
              irealb = elDatab(EL_REAL)
              i = rlget(irealb,rvrb(1))
              crad = rvrb(2)
          endif
      endif

c --- slave edge tolerence to prevent slipping off  
      if (rvr(25).lt.TINY) then
          tols = 0.2d0
      else
          tols = rvr(25)
      endif
      if (rvr(6).gt.TINY) then
          fpin = rvr(6)
      else
          fpin = 3.d0
          if (ififo(15).ge.1) fpin = 2.0d0
          if (abs(rvr(7)).gt.TINY.or.abs(rvr(8)).gt.TINY)
     &        fpin = 4.d0
          if (methi.eq.2.or.methi.eq.4) fpin = 4.d0
          if (etyiqr(ityp,-KYOP1).eq.5.or.etyiqr(ityp,-KYOP1).eq.4)
     &        fpin = fpin/4.d0
      endif

c --- nodal coordinates
      do 10 i = 1,nnodr
         if (nodes(i).eq.0) goto 10
         xyt(1,i) = dispn(nodes(i)) + coor0(nodes(i))
         xyt(2,i) = dispn(nodes(i)+1) + coor0(nodes(i)+1)
         xyt(3,i) = dispn(nodes(i)+2) + coor0(nodes(i)+2)
  10  continue

c --- calculate element length and centroid location
      length = CalDist2P(xyt(1,1),xyt(1,2))

c --- get thickness of attached beam element
      thick = min(abs(rvr(1)),abs(rvr(2)))
      if (thick.lt.TINY) thick = length*0.125d0
      ivolu = elmgct (elem,dvet(1))

c --- find stiffness
      esolid = dvet(5)
      dvet(6) = thick

      call elmpct (elem,dvet(1))

c --- calculate contact pinball range inside which contact may occur
      pinb = fpin*thick
      pinb = max(pinb,1.1*abs(rvr(10)))
      if (rvr(6).lt.ZERO) pinb = -rvr(6)
      pinb2 = pinb
      if (pinb2.le.2.d0*length) pinb2 = 2.d0*length

      if (ififo(11).lt.0.and.ififo(2).gt.0) goto 03

c --- check whether the current element is in the foundation cell
      xmax = xyt(1,1)
      ymax = xyt(2,1)
      zmax = xyt(3,1)
      xmin = xyt(1,1)
      ymin = xyt(2,1) 
      zmin = xyt(3,1)
      do 02 i = 2,3
         if (nodes(i).le.0) goto 02
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

      elvol = ZERO
      iseg = 0
      nocon = 2
      da = HUGE
      do 30 ipi=1,npi

         cnnode = nodes(ipi)
         if (cnnode.le.0) goto 30         
         x = xyt(1,ipi)
         y = xyt(2,ipi)
         z = xyt(3,ipi)
         dajc= 1.d0
         weight = 1.d0

c --- compute element volume and energy

         if (abs(rvr(2)).gt.TINY) then
             elvol = elvol + weight*dajc*abs(rvr(2))
         else
             elvol = elvol + weight*dajc*0.125*length
         endif
         if (ickcel.eq.1) goto 30

         if (methi.le.2) then
             gapshft = rvr(10) + crad
         else
             gapshft = crad
         endif

         call v_zero(r(1,1),9)
c --- compute geometric relations between structure and foundation
         call CalPeni04 (x,y,z,r(1,1),pinb,nSeg,ifoun,
     &                   nPerSeg,nPerFounI,nPerFounDp,ififo(1),
     &                   coor0(1),dispn(1),isegda(1,1),fcell(1),
     &                   founda(1),nocon,iseg,epsil(1),
     &                   delxy(1),tols,gapshft,elem,cnnode,
     &                   ibeam,hi)
         
         if (nocon.eq.0.and.abs(epsil(3)).gt.pinb) nocon = 2

         if (nocon.eq.2) then
             iseg = 0
             goto 30
         endif

c --- get the closest contact point information
         if (abs(epsil(3)).lt.pinb.and.iseg.gt.0) then
             if (epsil(3).lt.founda(9)) then
                 founda(9) = epsil(3)
                 ififo(16) = elem
                 ififo(17) = isegda(15,iseg)
                 founda(10) = delxy(1)
                 founda(11) = delxy(2)
                 founda(25) = thick
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
         if (nocon.eq.0 .and. abs(epsil(3)).lt.da) then
             icnda(3,icon) = iseg
             isegda(19,iseg) = icon
             da = epsil(3)
         endif
  30  continue
      founda(16) = founda(16) + length
      founda(17) = founda(17) + thick
      founda(18) = founda(18) + esolid
      founda(19) = founda(19) + 1.d0
      founda(20) = founda(20) + pinb
      founda(21) = founda(21) + elvol

      return
      end


