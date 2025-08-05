      subroutine CalPeni (cnodes,ickcee,x0,y0,z0,xa,ya,za,vax,vay,vaz,
     &                    r,pinb,nSeg,ifoun,nPerSeg,nPerFounI,
     &                    nPerFounDp,iendp,ififo,coor0,isegda,founda,
     &                    nocon,isega,iseg,epsil,fifoa,fifob,tols,
     &                    gapshft,cnnode,krot,qab,ta,tb,ac,bc,elem)
c --- primary function: compute the relations between an integration point
c                       of the contact element and the corresponding foundation
c  input arguments:
c      cnodes (int,ar(2),in) - contact node list
c      xa,ya,za (dp,ar(3),in) - coordinates of contact nodes
c      vax,vya,vaz(dp,ar(3),in)  - components of the velocity of the integration
c                                point in global cartesian coordinates
c      dx,dy       (dp,sc,in)  - components of a unit vector tangent to
c                                the element at the integration point,
c                                positive in the direction (node1 to node3))
c      r(3,3)      (dp,ar(3,3),in)
c                            - components of the normal vector to the structure
c                                only for 3-d case
c      dmax        (dp,sc,in)  - maximum distance for a valid contact
c                                (if dabs(d)>dmax, false contact is assumed)
c      nSeg       (int,sc,in) - total number of segments defining the
c                                 foundations
c      nPerSeg       (int,sc,in) - number of items per segment database entry
c      nPerFounI      (int,sc,in) - number of integer items per foundation
c      iendp       (int,sc,in) = 0 beginning of sub-step
c                              = 1 end of sub-step
c      ififo    (dp,ar(6),in)
c                                  - foundation information (integer)
c          ififo(1)       = number of segments
c          ififo(4)       = constrain condition in x-dir.c          
c          ififo(5)       = constrain condition in y-dir.c 
c          ififo(6)       = constrain condition in z-dir.c
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
c
c
c      coor0     (dp,ar(*),in)     - initial nodal coordinate coordinate
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
c
c  output arguments:
c      x0 ,y0 ,z0  (dp,sc,out)  - global coordinates of the integration point
c      r(3,3)   (dp,ar(3,3),out)
c                         - components of the normal vector to the structure
c                                only for 3-d case
c       epsil    (dp,ar(4),out)
c        epsil(1)   = tangent velocity of structure in first dir.
c        epsil(2)   = tangent velocity of structure in second dir.
c        epsil(3)
c                d  = distance between the integration point
c                     and the foundation
c                     ( d < 0 : contact    )
c                     ( d > 0 : no contact )
c        epsil(7)   = angle between beams
c
c       nocon   (int,sc,out)
c                   = 0  : contact
c       iseg  (int,sc,inout)    - number of the foundation segment at current
c       isega (int,sc,inout)    - number of the foundation segment at initial
c                                 with which contact is established
c                                 (=0 if no contact with the foundation)
c       ta,tb (dp,ar(3),out)    - tangents to two beams at contact  
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "cont_info.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
c
      integer nSeg,ifoun,nPerSeg,nPerFounI,nPerFounDp,iendp,nocon,
     &  iseg,ififo(nPerFounI),isegda(nPerSeg,nSeg),isega,
     &  cnnode(2),locu,locv,cnodes(7),krot,numRoot,try,
     &  iproc,iProcSMP,debug,elem
      double precision x0,y0,z0,xa(7),ya(7),za(7),vax(7),vay(7),vaz(7),
     &  dmax,founda(nPerFounDp),fifoa(7),fifob(7),r(3,3),epsil(7),
     &  coor0(*),delxy(3),qab(34),ro(3,3),ta(3),tb(3),ac(3),bc(3)

      integer is,i1,kddlx,kddly,kddlz,piseg,pnocon,
     &  ickcel,iseg1,iseg2,ickcee,is0,i2,tnodes(7),
     &  ipass,nnode,inode,locdis,locdin,iconnect(2),
     &  n1,n2,contactNode(2),iseg0,iseg00,iseg000,ipass1,sseg
      double precision ZERO,HALF,ONE,gapshft,
     &  xb(7),yb(7),zb(7),vbx(7),vby(7),vbz(7),
     &  d0,pinb,tols,ptols,rx(7),dum,
     &  ry(7),rz(7),vrx(7),vry(7),vrz(7),
     &  BestU,BestV,BestD,
     &  pBestU,pBestV,pBestD     

      pointer (pdDinL,Din)
      pointer (pdDisL,Dis)
      pointer (piCalPene,CalPene)
      double precision  Din(*), Dis(*)
      byte CalPene(*)
      
      parameter (ZERO=0.0d0, HALF=0.5d0, ONE=1.0d0)

#if defined(DEBUG)
      call SubBeg('CalPeni')
#endif

      debug = 0

c      if (ififo(12).eq.8812) then
c          debug = 1
c      endif

c      if (ififo(12).eq.14842) then
c          debug = 1
c      endif

c      if (elem.eq.91648) then
c          debug = 1
c      endif
   
    
c      if (debug.gt.0) then
c          write (6,*) ' <<<<<<< START ---- ',elem,ififo(12)
c      endif


      call v_zero (fifoa(1),7)
      call v_zero (fifob(1),7)
      call v_zero (epsil(1),7)
      call v_move (r(1,1),ro(1,1),9)
      
      dmax = pinb
      epsil(3) = dmax
      is0 = 0
      try = 0

      if (A_KeyFirstPass.eq.0) then
          pBestU = qab(19)
          pBestV = qab(20)
          pBestD = qab(21)
          sseg = int(qab(22))
      else
          pBestU = HUGE
          pBestV = HUGE
          pBestD = HUGE
          sseg = 0
      endif
          
      piseg = iseg
c      pnocon = 2
      
      call vi_zero (iconnect(1),2)

      iproc = iProcSMP()
      piCalPene = K_PtrConSrch + iproc*nSeg

      if (ifoun.le.0) goto 070

      ickcel = ickcee
      if (ickcee.ge.2) ickcel = 0
      if (ickcel.eq.1) goto 070
      ickcel = ickcee

c --- degenerated case
c
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
                 vax(inode) = HALF*(vax(n1)+vax(n2))
                 vay(inode) = HALF*(vay(n1)+vay(n2))
                 vaz(inode) = HALF*(vaz(n1)+vaz(n2))
             else
                 xa(inode) = ZERO
                 ya(inode) = ZERO
                 za(inode) = ZERO
                 vax(inode) = ZERO
                 vay(inode) = ZERO
                 vaz(inode) = ZERO
             endif
         endif
      enddo

      d0 = dmax

c --- contact node
      call vi_zero (contactNode(1),2)
      do i2 =1,2
         if (cnnode(i2).gt.0) then
             contactNode(i2) = (M_FrwdList(cnnode(i2))-1)*D_NumDof + 1
         endif
      enddo
      ipass = 0
      ipass1 = 0
      iseg000 = iseg
      iseg00 = 0
      iseg0 = 0

      if (iseg.eq.isega.or.iseg.eq.0) then
          ipass = 1
          iseg = isega
      elseif (iseg.gt.0) then
          ipass = 1
          iseg0 = iseg
      endif

      iseg1 = ififo(8)
      iseg2 = ififo(9)
      
      do is = iseg1,iseg2
         CalPene(is) = 0
         if (isegda(15,is).lt.0) CalPene(is) = 1
      enddo
      is = iseg1 - 1
  020 continue
      is0 = iseg
      if (iseg.gt.0) then
          if (CalPene(iseg).eq.0)  goto 100
      endif
      
      if (ipass.eq.0) then
          ipass = 1
          iseg = isega
          goto 020
      endif
      is0 = iseg1
  050 continue
      iseg0 = 0
      ipass1 = 0
      do 060 is = iseg1,iseg2
         if (CalPene(is).ne.0) goto 060
         if (elem.eq.isegda(17,is)) then
             CalPene(is) = 1
             goto 060
         endif
         iseg = is
         goto 100
  060 continue
  070 iseg = is0

      if (abs(epsil(3)).gt.dmax) epsil(3) = dmax
      nocon = 2
      goto 999

  100 continue
      iseg00 = iseg0
      iseg0 = iseg
      CalPene(iseg) = 1

c --- geometrical constants at the end of the step

      nnode = 7

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
             tnodes(inode) = 0
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
                 vbx(inode) = HALF*(vbx(n1)+vbx(n2))
                 vby(inode) = HALF*(vby(n1)+vby(n2))
                 vbz(inode) = HALF*(vbz(n1)+vbz(n2))
                 rx(inode) = HALF*(rx(n1)+rx(n2))
                 ry(inode) = HALF*(ry(n1)+ry(n2))
                 rz(inode) = HALF*(rz(n1)+rz(n2))
                 vrx(inode) = HALF*(vrx(n1)+vrx(n2))
                 vry(inode) = HALF*(vry(n1)+vry(n2))
                 vrz(inode) = HALF*(vrz(n1)+vrz(n2))
             else
                 xb(inode) = ZERO
                 yb(inode) = ZERO
                 zb(inode) = ZERO
                 vbx(inode) = ZERO
                 vby(inode) = ZERO
                 vbz(inode) = ZERO
                 rx(inode) = ZERO
                 ry(inode) = ZERO
                 rz(inode) = ZERO
                 vrx(inode) = ZERO
                 vry(inode) = ZERO
                 vrz(inode) = ZERO
             endif
             cycle
         endif

c --- contact node cannot contact its own elements
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
         pdDisL = M_PtrDisp
         locdis = kddlx
         pdDinL = M_PtrIncS
         locdin = kddlx

         if (abs(Dis(locdis)).eq.HUGE) then
             xb(inode) = coor0(kddlx) + founda(13)
         else
             xb(inode) = Dis(locdis) + coor0(kddlx) + founda(13)
         endif
         if (abs(Dis(locdis+1)).eq.HUGE) then
             yb(inode) = coor0(kddly) + founda(14)
         else
             yb(inode) = Dis(locdis+1) + coor0(kddly) + founda(14)
         endif
         vbx(inode) = Din(locdin)
         vby(inode) = Din(locdin+1)

c --- initial configuration
         if (iendp.eq.-1) then
             xb(inode) = coor0(kddlx) + founda(13)
             yb(inode) = coor0(kddly) + founda(14)
             vbx(inode) = ZERO
             vby(inode) = ZERO
         endif

         if (abs(Dis(locdis+2)).eq.HUGE) then
             zb(inode) = coor0(kddlz) + founda(15)
         else
             zb(inode) = Dis(locdis+2) + coor0(kddlz) + founda(15)
         endif
         vbz(inode) = Din(locdin+2)

         if (iendp.eq.-1) then
             zb(inode) = coor0(kddlz) + founda(15)
             vbz(inode) = ZERO
         endif

         if (krot.eq.1) then
             pdDisL = M_PtrSoluVects
             locdis = kddlx
             pdDinL = M_PtrIncS
             locdin = kddlx
             if (abs(Dis(locdis+3)).eq.HUGE) then
                 rx(inode) = ZERO
             else
                 rx(inode) = Dis(locdis+3)
             endif
             if (abs(Din(locdin+3)).eq.HUGE) then
                 vrx(inode) = ZERO
             else
                 vrx(inode) = Din(locdin+3)
             endif
             if (abs(Dis(locdis+4)).eq.HUGE) then
                 ry(inode) = ZERO
             else
                 ry(inode) = Dis(locdis+4)
             endif
             if (abs(Din(locdin+4)).eq.HUGE) then
                 vry(inode) = ZERO
             else
                 vry(inode) = Din(locdin+4)
             endif
             if (abs(Dis(locdis+5)).eq.HUGE) then
                 rz(inode) = ZERO
             else
                 rz(inode) = Dis(locdis+5)
             endif
             if (abs(Din(locdin+5)).eq.HUGE) then
                 vrz(inode) = ZERO
             else
                 vrz(inode) = Din(locdin+5)
             endif
         else
             rx(inode) = ZERO
             vrx(inode) = ZERO
             ry(inode) = ZERO
             vry(inode) = ZERO
             rz(inode) = ZERO
             vrz(inode) = ZERO
         endif

      enddo

c --- at the beginning of the step

      if (iendp.eq.0) then
          do inode = 1,nnode
             xb(inode) = xb(inode) - vbx(inode)
             yb(inode) = yb(inode) - vby(inode)
             zb(inode) = zb(inode) - vbz(inode)
             rx(inode) = rx(inode) - vrx(inode)
             ry(inode) = ry(inode) - vry(inode)
             rz(inode) = rz(inode) - vrz(inode)
          enddo
      endif

c --- case of a straight line
      d0 = epsil(3)

c --- fix slave segment
      if (ickcel.eq.3) ipass1 = -1

c --- beam segment
      numRoot = 0

      if (A_KeyFirstPass.eq.1 .or. iseg.ne.sseg) then
          BestU = HUGE
          BestV = HUGE
          if (A_KeyFirstPass.ne.1) then
              ptols = tols*0.8d0
          else
              ptols = tols
          endif
      else
          BestU = qab(19)
          BestV = qab(20)
          BestD = qab(21)
          ptols = tols
      endif
      
      call BeamToBeam (founda(7),x0,y0,z0,xa(1),ya(1),za(1),
     &                 vax(1),vay(1),vaz(1),xb(1),yb(1),zb(1),
     &                 vbx(1),vby(1),vbz(1),r(1,1),ro(1,1),epsil(1),
     &                 fifoa(1),fifob(1),nocon,dmax,delxy(1),
     &                 isegda(11,iseg),ptols,gapshft,ipass1,
     &                 rx(1),ry(1),rz(1),locu,locv,numRoot,
     &                 BestU,BestV,BestD,ta(1),tb(1),ac(1),bc(1),
     &                 debug)
     
      if (nocon.lt.2 .or. abs(epsil(3)).lt.dmax) then
          if (abs(BestV).le.(ONE+0.01d0*tols)) then
c              if (debug.gt.0) then
c                  write (6,*) 'U 1, nocon,epsil(3),BestV,pBestV',
c     &                     nocon,epsil(3),BestV,pBestV,(ONE+0.01d0*tols)
c                  write (6,*) 'cnodes', cnodes(1:14)
c                  write (6,*) 'tnodes', tnodes(1:7)
c              endif
              pBestU = BestU
              pBestV = BestV
              pBestD = BestD
              piseg = iseg
c              pnocon = nocon
          elseif ((pBestU**2+pBestV**2).gt.(BestU**2+BestV**2)) then
c              if (debug.gt.0) then
c                 write (6,*) 'U 2, nocon,epsil(3),BestV,pBestV',
c     &                     nocon,epsil(3),BestV,pBestV,(ONE+0.01d0*tols)
c                  write (6,*) 'cnodes', cnodes(1:14)
c                  write (6,*) 'tnodes', tnodes(1:7)
c              endif
              pBestU = BestU
              pBestV = BestV
              pBestD = BestD
              piseg = iseg
c              pnocon = nocon
          endif
      endif

 200  continue

c --- contact state
      if (abs(epsil(3)).gt.dmax) then
          epsil(3) = dmax
          if (nocon.eq.0) nocon = 2
      endif

      if (nocon.eq.0 .and. ipass1.gt.0 .and. iseg000.gt.0) then
          if (iseg.ne.iseg000 .and. try.eq.0) then
              try = 1
              nocon = 2
              iseg = iseg000
              goto 100
          endif    
      elseif (try.eq.1 .and. nocon.ne.0) then
          try = 2
          nocon = 2
          iseg = iseg00
          goto 100
      elseif (try.eq.2 .and. nocon.ne.0) then
          call erhandler ('CalPeni',5000,2,'Error in CalPeni.
     &                     This should never happen.',ZERO,' ')
      endif
      
c      if (nocon.eq.0 .and. ipass1.gt.0) then
c          epsil(3) = dmax
c          nocon = 2
c          goto 999
c      endif
      
      if (nocon.eq.0) goto 999
      if (ickcel.eq.3) then
c --- contact element loses contact with slave element
          ickcel = 0
          ipass1 = 0
      endif

      if (epsil(3).lt.d0) then
          is0 = iseg
      else
          epsil(3) = d0
      endif

      if (ipass.eq.0) then
          ipass = 1
          iseg = isega
          goto 020
      endif
      iseg1 = is + 1
      if (nocon.eq.2) goto 050

      if (nocon.eq.-1) iseg = isegda(11,iseg)
      if (nocon.eq. 1) iseg = isegda(12,iseg)

      if (iseg.eq.0) goto 050

      if (CalPene(iseg).eq.0) then
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
          if (cnodes(3).eq.0) then
              fifoa(1) = fifoa(1)+HALF*fifoa(3)
              fifoa(2) = fifoa(2)+HALF*fifoa(3)
              fifoa(3) = ZERO
          endif
          if (cnodes(6).eq.0.and.cnodes(4).gt.0) then
              fifoa(1) = fifoa(1)+HALF*fifoa(6)
              fifoa(4) = fifoa(4)+HALF*fifoa(6)
              fifoa(6) = ZERO
          endif
          if (cnodes(7).eq.0.and.cnodes(5).gt.0) then
              fifoa(2) = fifoa(2)+HALF*fifoa(7)
              fifoa(5) = fifoa(5)+HALF*fifoa(7)
              fifoa(7) = ZERO
          endif
          if (tnodes(3).eq.0) then
              fifob(1) = fifob(1)+HALF*fifob(3)
              fifob(2) = fifob(2)+HALF*fifob(3)
              fifob(3) = ZERO
          endif
          if (tnodes(6).eq.0.and.tnodes(4).gt.0) then
              fifob(1) = fifob(1)+HALF*fifob(6)
              fifob(4) = fifob(4)+HALF*fifob(6)
              fifob(6) = ZERO
          endif
          if (tnodes(7).eq.0.and.tnodes(5).gt.0) then
              fifob(2) = fifob(2)+HALF*fifob(7)
              fifob(5) = fifob(5)+HALF*fifob(7)
              fifob(7) = ZERO
          endif
          dum = fifob(2)           
          fifob(2) = fifob(3)
          fifob(3) = dum
      endif


      if (nocon.eq.0 .and. iseg.ne.isega .and. gapshft.eq.ZERO .and. 
     &    A_KeyFirstPass.eq.0. and. abs(qab(21)).lt.dmax) then
          if ((epsil(3)-qab(21)).le.ZERO) then
              gapshft = qab(21)
              epsil(3) = epsil(3) - gapshft
          elseif (epsil(3).lt.ZERO) then
              gapshft = epsil(3)
              epsil(3) = ZERO
          endif    
      endif

c      if (pnocon.lt.2 .or. isega.gt.0 .or. iseg000.gt.0) then
      if (piseg.gt.0) then
          qab(19) = pBestU
          qab(20) = pBestV
          qab(22) = dble(piseg)
      else
          qab(19) = HUGE
          qab(20) = HUGE
          qab(22) = ZERO
      endif

      qab(21) = pBestD
      iseg = piseg
      
c      if (debug.gt.0) then
c          write (6,*) 'isegs', iseg,isega,iseg000
c          write (6,*) 'cnodes', cnodes(1:14)
c          write (6,*) 'tnodes', tnodes(1:7)
c          write (6,*) 'DP qab:', qab(19:21)
c          write (6,*) 'DP gapshft:', gapshft,epsil(3)
c          write (6,*) 'INT   :', ipass1,numRoot
c          write (6,*) ' <<<<<<< FINISH ---- ',nocon
c          write (6,*) ''
c      endif
      

#if defined(DEBUG)
      call SubEnd('CalPeni')
#endif

      return
      end
