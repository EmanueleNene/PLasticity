      subroutine CalPeni4 (ickcee,x0,y0,z0,vx0,vy0,vz0,dx,r,
     &                     pinb,nSeg,ifoun,nPerSeg,nPerFounI,nPerFounDp,
     &                     iendp,ififo,coor0,isegda,fcell,founda,
     &                     nocon,isega,iseg,epsil,fifon,tols,gapshft,
     &                     elem,cnnode,krot,natcor,df,ibeam,
     &                     rad,hi)
c --- primary function: compute the relations between an integration point
c                       of the contact element and the corresponding foundation
c
c  input arguments:
c
c      x0 ,y0 ,z0  (dp,sc,in)  - global coordinates of the integration point
c      vx0,vy0,vz0 (dp,sc,in)  - components of the velocity of the integration 
c                                point in global cartesian coordinates
c      r(3,3)      (dp,ar(3,3),in)
c                              -  components of the normal vector to the structure
c      dmax        (dp,sc,in)  - maximum distance for a valid contact
c                                (if dabs(d)>dmax, false contact is assumed)
c      nSeg        (int,sc,in) - total number of segments defining the
c                                 foundations
c      nPerSeg      (int,sc,in) - number of items per segment database entry
c      nPerFounI    (int,sc,in) - number of integer items per foundation
c      iendp       (int,sc,in) = 0 beginning of sub-step
c                              = 1 end of sub-step
c      ififo       (dp,ar(6),in)
c                              - foundation information (integer)
c          ififo(1)       = number of segments
c          ififo(2)       = pilot nodal number
c          ififo(3)       = pilot nodal number on rotation
c          ififo(4)       = constrain condition in x-dir.
c          ififo(5)       = constrain condition in y-dir.
c          ififo(6)       = constrain condition in z-dir.
c          ififo(7)       = dimensions
c          ififo(8)       = 1st segment numbering on this foundation
c          ififo(9)       = 2nd segment numbering on this foundation
c          ififo(10)      = lowest nodal numbering on this foundation
c          ififo(11)      = 0 pilot node does not connect to other element
c                         =   element number that pilot node connects to
c          ififo(12)      = real constant number (ireal)
c          ififo(13)      = contact pair key
c                            0 - no contact pair
c                            1 - only slave element is defined
c                            2 - only contact element is defined
c                            3 - contact pair is defined
c          ififo(14)      = type of associated contant element type (ityp)
c          ififo(15)      = type of deformable slave surface
c                           -2 - mixed type surface (may wrong)
c                           -1 - rigid slave surface
c                            0 - none
c                            1 - deformable slave surface
c                            2 - deformable slave surface (TWO pass)
c          ififo(16)      = contact element number who is closest to
c                            current slave
c          ififo(17)      = slave element number who is pair of closest
c                            contact element
c          ififo(18)      = iteration number of surface adjust
c          ififo(19,ifoun)      = max. node per slave element on each
c                                 slave surface
c
c      coor0     (dp,ar(*),in)     - initial nodal coordinate coordinate
c      isegda   (int,ar(nPerSeg,nSeg),in)
c                                  - segment information 
c          isegda(1,iseg)        = code of segment type (icode)
c                                = 1   : straight line                
c                                = 2   : parabolic arc                
c                                = 7   : quadrilateral                
c            
c          isegda(2,iseg)        = fondation number           
c          isegda(j,iseg)        = 4 nodes (j=3-10) on seg. iseg
c          isegda(j,iseg)        = 4 adjacent segments (j=11-14)
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
c       fcell   (dp,ar(4),in) - foundation  cell region
c          fcell(1)         = minimum value in x-dir.
c          fcell(2)         = maximum value in x-dir.
c          fcell(3)         = minimum value in y-dir.
c          fcell(4)         = maximum value in y-dir.
c
c  output arguments:
c
c       epsil    (dp,ar(4),out)
c       epsil(1)   = tangent velocity of structure in first dir.
c       epsil(2)   = tangent velocity of structure in second dir.
c       epsil(3)    
c               d  = distance between the integration point
c                    and the foundation
c                     ( d < 0 : contact    )
c                     ( d > 0 : no contact )
c        epsil(4)   = foundation moisture
c
c       fifon   (dp,ar(8),out) - shape functions of the foundation
c       nocon   (int,sc,out)
c                   = 0  : contact
c                   = 11,12 or 13 | one adjacent segment
c       iseg  (int,sc,inout)    - number of the foundation segment at current 
c       isega (int,sc,inout)    - number of the foundation segment at initial
c                                 with which contact is established
c                                 (=0 if no contact with the foundation)
c
#include "impcom.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "cont_info.inc"
#include "solu_info.inc"
c#include "debug_tmp.inc"
c
      integer nSeg,ifoun,nPerSeg,nPerFounI,nPerFounDp,iendp,nocon,
     &        iseg,ififo(nPerFounI),isegda(nPerSeg,nSeg),isega,ibeam,
     &        elem,cnnode,iProcSMP,ipi,i,j,k,hi
c     &        ,sseg,try
      double precision x0,y0,z0,vx0,vy0,vz0,dmax,founda(nPerFounDp),
     &        fifon(7),r(3,3),epsil(8),coor0(*),fcell(*),
     &        delxy(3)

      integer icrot,is,icode,inopil,iropil,i1,kddlx,kddly,kddlz,
     &        ickcel,iseg1,iseg2,ickcee,is0,iproc,n1,n2,
     &        ipass,nnode,inode,locdis,locdin,iconnect,tnodes(7),
     &        contactNode,iseg0,iseg00,iseg000,ipass1,krot,locp,
     &        trg3dbeam,jbeam,debug
c     &        ,piseg
      double precision ZERO,HALF,ONE,xp,yp,zp,gapshft,
     &        x(7),y(7),z(7),vx(7),vy(7),vz(7),vxp,vyp,vzp,xpn,ypn,
     &        zpn,d0,xpn1,ypn1,zpn1,omeg1,omeg2,omeg3,pinb,
     &        tols,nx(7),ny(7),nz(7),rx(7),dum,dx,
     &        ry(7),rz(7),vrx(7),vry(7),vrz(7),natcor(5),
     &        dfi(2,9,4),df(10),rad(2),GetFibRad,
     &        ro(3,3),dxdksi,dxdeta,dydksi,dydeta,dzdksi,dzdeta,
     &        dtdx,dtdy,dtdz,dajc,ptols
c     &        ,roots(2),pBestU,pBestV,BestU,BestV

      parameter (ZERO=0.0d0,HALF=0.5d0,ONE=1.0d0)
      
      data  (((dfi(j,i,k),j=1,2),i=1,4),k=1,4) /
     & -0.5d0,-0.5d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0,
     & -0.5d0, 0.0d0, 0.5d0,-0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0,
     &  0.0d0, 0.0d0, 0.0d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0, 0.0d0,
     &  0.0d0,-0.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0,-0.5d0, 0.5d0/

      pointer (pdDinL,Din)
      pointer (pdDisL,Dis)
      pointer (piCalPene,CalPene)
      double precision  Din(*), Dis(*)
      byte CalPene(*)
      
#if defined(DEBUG)
      call SubBeg('CalPeni4')
#endif

c --- initialize
      call v_zero (fifon(1),7)
      call v_zero (epsil(1),8)
      call v_zero (df(1),10)
      call v_move (r(1,1),ro(1,1),9)
      dmax = pinb
      epsil(3) = dmax
      is0 = 0
      iconnect = 0

c      iseg = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
c      piseg = iseg
c      pBestV = HUGE
c      try = 0
c      if (A_KeyFirstPass.eq.0) then
c          pBestU = roots(1)
c          sseg = int(roots(2))
c      else
c          pBestU = HUGE
c          sseg = 0
c      endif      
      
c --- create space for each processor to mark the need for calculations
      iproc = iProcSMP()
      piCalPene = K_PtrConSrch + iproc*nSeg

      rad(1) = GetFibRad(ibeam)
      if (founda(7).lt.ZERO) rad(1) = -rad(1)
      if (abs(rad(1)).le.TINY) rad(1) = founda(7)
      
      if (ifoun.le.0) goto 070
    
      ickcel = ickcee

c --- check the cell region
      if (ickcel.ne.1) then
          if (x0.lt.fcell(1) .or. x0.gt.fcell(2)) ickcel = 1
          if (y0.lt.fcell(3) .or. y0.gt.fcell(4)) ickcel = 1
          if (z0.lt.fcell(5) .or. z0.gt.fcell(6)) ickcel = 1
       endif
      if (ickcee.ge.2) ickcel = 0
      if (ickcel.eq.1) goto 070

      ickcel = ickcee
c --- indicate that the calculation for rotation has done (icrot.ne.0)
c                                                or no    (icrot=0)
      icrot = 0
      d0 = dmax

c --- contact node
      contactNode = 0
      if (cnnode.gt.0) contactNode = (M_FrwdList(cnnode)-1)*D_NumDof + 1
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
          if (CalPene(iseg).eq.0 .and. isegda(1,iseg).ne.99) goto 100
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
         if (isegda(1,is).eq.99) goto 060 
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
      icode = isegda(1,iseg)
      
      CalPene(iseg) = 1
      inopil = ififo(2)
      iropil = ififo(3)

c --- calculate the effect for the pilot foundation with rotation
      if (iropil.gt.0 .and. icrot.eq.0) then
          kddlx = inopil
          kddly = kddlx + 1
          kddlz = kddlx + 2
         
          pdDisL = M_PtrDisp
          locdis = kddlx
          pdDinL = M_PtrIncS
          locdin = kddlx

          if (abs(Dis(locdis)).eq.HUGE) then
              xpn = coor0(kddlx) + founda(13)
          else
              xpn = Dis(locdis) + coor0(kddlx) + founda(13)
          endif
          if (abs(Dis(locdis+1)).eq.HUGE) then
              ypn = coor0(kddly) + founda(14)
          else
              ypn = Dis(locdis+1) + coor0(kddly) + founda(14)
          endif
          if (abs(Dis(locdis+2)).eq.HUGE) then
              zpn = coor0(kddlz) + founda(15)
          else
              zpn = Dis(locdis+2) + coor0(kddlz) + founda(15)
          endif

          xp = xpn - Din(locdin)
          yp = ypn - Din(locdin+1)
          zp = zpn - Din(locdin+2)

          vxp = Din(locdin)
          vyp = Din(locdin+1)
          vzp = Din(locdin+2)

          omeg1 = Din(locdin+3)
          omeg2 = Din(locdin+4)
          omeg3 = Din(locdin+5)

          icrot = 1
      endif

c --- geometrical constants 

      nnode = 3
      if (icode.eq.7) nnode = 4
      trg3dbeam = 0
      if (icode.eq.1.or.icode.eq.2) then
c --- 3d beam target
          tnodes(1) = isegda(3,iseg)
          tnodes(2) = isegda(5,iseg)
          tnodes(3) = isegda(4,iseg)
          tnodes(4) = isegda(6,iseg)
          tnodes(5) = isegda(7,iseg)
          tnodes(6) = isegda(8,iseg)
          tnodes(7) = isegda(9,iseg)
          if (tnodes(3).le.0.and.tnodes(6).le.0.and.tnodes(7).le.0) then
              icode = 1
              nnode = 5
          else
              icode = 2
              trg3dbeam = 1
              nnode = 7
          endif
      endif
      iconnect = 0
      do inode = 1,nnode
         if (trg3dbeam.eq.1) then
             i1 = tnodes(inode)
         else
             i1 = isegda(2+inode,iseg)
         endif
         if (i1.le.0) then
             n1 = 0
             if (trg3dbeam.eq.1) then
                 tnodes(inode) = 0
                 if (inode.eq.3) then
                     n1 = 1
                     n2 = 2
                 endif
                 if (inode.eq.4) then
                     n1 = 1
                     n2 = 1
                 endif
                 if (inode.eq.5) then
                     n1 = 2
                     n2 = 2
                 endif
                 if (inode.eq.6) then
                     n1 = 1
                     n2 = 4
                 endif
                 if (inode.eq.7) then
                     n1 = 2
                     n2 = 5
                 endif
             endif
             if (n1.gt.0) then
                 x(inode) = HALF*(x(n1)+x(n2))
                 y(inode) = HALF*(y(n1)+y(n2))
                 z(inode) = HALF*(z(n1)+z(n2))
                 vx(inode) = HALF*(vx(n1)+vx(n2))
                 vy(inode) = HALF*(vy(n1)+vy(n2))
                 vz(inode) = HALF*(vz(n1)+vz(n2))
                 rx(inode) = HALF*(rx(n1)+rx(n2))
                 ry(inode) = HALF*(ry(n1)+ry(n2))
                 rz(inode) = HALF*(rz(n1)+rz(n2))
                 vrx(inode) = HALF*(vrx(n1)+vrx(n2))
                 vry(inode) = HALF*(vry(n1)+vry(n2))
                 vrz(inode) = HALF*(vrz(n1)+vrz(n2))
             else
                 x(inode) = ZERO
                 y(inode) = ZERO
                 z(inode) = ZERO
                 vx(inode) = ZERO
                 vy(inode) = ZERO
                 vz(inode) = ZERO
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
         if (contactNode.gt.0) then
             if (i1.eq.contactNode) then
                 if (ififo(15).eq.3) then
                     nocon = 2
                     epsil(3) = d0
                     goto 200
                 endif
                 iconnect = 1
             endif
         endif
         kddlx = i1
         kddly = kddlx + 1
         kddlz = kddlx + 2
         pdDisL = M_PtrDisp
         locdis = kddlx
         pdDinL = M_PtrIncS
         locdin = kddlx

         if (abs(Dis(locdis)).eq.HUGE) then
             x(inode) = coor0(kddlx) + founda(13)
         else
             x(inode) = Dis(locdis) + coor0(kddlx) + founda(13)
         endif
         if (abs(Dis(locdis+1)).eq.HUGE) then
             y(inode) = coor0(kddly) + founda(14)
         else
             y(inode) = Dis(locdis+1) + coor0(kddly) + founda(14)
         endif
         if (abs(Dis(locdis+2)).eq.HUGE) then
             z(inode) = coor0(kddlz) + founda(15)
         else
             z(inode) = Dis(locdis+2) + coor0(kddlz) + founda(15)
         endif

         vx(inode) = Din(locdin)
         vy(inode) = Din(locdin+1)
         vz(inode) = Din(locdin+2)
c --- initial configuration
         if (iendp.eq.-1) then
             x(inode) = coor0(kddlx) + founda(13)
             y(inode) = coor0(kddly) + founda(14)
             z(inode) = coor0(kddlz) + founda(15)
             vx(inode) = ZERO
             vy(inode) = ZERO
             vz(inode) = ZERO
         endif

         if (krot.eq.1) then
             pdDisL = M_PtrDisp
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
             x(inode) = x(inode) - vx(inode)
             y(inode) = y(inode) - vy(inode)
             z(inode) = z(inode) - vz(inode)
             rx(inode) = rx(inode) - vrx(inode)
             ry(inode) = ry(inode) - vry(inode)
             rz(inode) = rz(inode) - vrz(inode)
          enddo
      endif

c --- calculate the velocity on foundation due to rotation on pilot node
      if (iropil.gt.0) then
          if (iendp.eq.0) then
              xpn1 = xp
              ypn1 = yp
              zpn1 = zp
          else
              xpn1 = xpn
              ypn1 = ypn
              zpn1 = zpn
          endif
          if (iendp.eq.-1) then
              xpn1 = ZERO
              ypn1 = ZERO
              zpn1 = ZERO
              vxp = ZERO
              vyp = ZERO
              vzp = ZERO
              omeg1 = ZERO
              omeg2 = ZERO
              omeg3 = ZERO
          endif
      endif

      d0 = epsil(3)

c --- fix slave segment
      if (ickcel.eq.3) ipass1 = -1
      if (ickcel.eq.4) ipass1 = -2

      if (icode.eq.1 .or. icode.eq.2) then
          jbeam = abs(isegda(16,iseg))
          rad(2) = GetFibRad(jbeam)
          if (rad(2).le.TINY) rad(2) = founda(8)
c          if (A_KeyFirstPass.eq.1 .or. iseg.ne.sseg) then
c              BestU = HUGE
c          else
c              BestU = roots(1)
c          endif      
c          pBestV = ZERO
      endif

      if (A_KeyFirstPass.eq.1 .or. iseg.ne.iseg000) then
          if (A_KeyFirstPass.ne.1) then
              ptols = tols*0.8d0
          else
              ptols = tols
          endif
      else
          ptols = tols
      endif 
      
      if (icode.eq.1) then
c --- 3D line segment
          debug = 0
          dmax = dmax + abs(rad(1)) + abs(rad(2))
          call ParBeam2nd (rad(1),x(1),y(1),z(1),vx(1),vy(1),vz(1),
     &                     x0,y0,z0,vx0,vy0,vz0,r(1,1),ro(1,1),
     &                     epsil(1),fifon(1),df(1),nocon,dmax,iropil,
     &                     vxp,vyp,vzp,xpn1,ypn1,zpn1,omeg1,omeg2,omeg3,
     &                     delxy(1),isegda(11,iseg),ptols,gapshft,
     &                     ipass1,rx(1),ry(1),rz(1),locp,natcor(1),HUGE,
     &                     hi,debug)
c          if (debug.eq.1.and.nocon.ne.2) write (6,*) 'iseg', iseg
          if (nocon.eq.0) then
              dx = locp
              fifon(2) = ZERO
              fifon(6) = ZERO
              fifon(7) = ZERO
              df(2) = ZERO
              df(6) = ZERO
              df(7) = ZERO
c --- no curvature effect
              df(9) = ZERO
              df(10) = ZERO          
          endif
      elseif (icode.eq.2) then
          debug = 0
c --- 3D parabolic segment
c          X_DummyI = X_DummyI + 1
          
c          if (elem.eq.65998) debug = 1
         
          dmax = dmax + abs(rad(1)) + abs(rad(2))
          call ParBeam3nd (rad(1),x(1),y(1),z(1),vx(1),vy(1),vz(1),
     &                     x0,y0,z0,vx0,vy0,vz0,r(1,1),ro(1,1),
     &                     epsil(1),fifon(1),df(1),nocon,dmax,iropil,
     &                     vxp,vyp,vzp,xpn1,ypn1,zpn1,omeg1,omeg2,omeg3,
     &                     delxy(1),isegda(11,iseg),ptols,gapshft,
     &                     ipass1,rx(1),ry(1),rz(1),locp,natcor(1),HUGE,
     &                     hi,debug)

c          if (debug.eq.1) write (6,*)'elem,iseg,iseg0',elem,iseg,iseg0
c          if (debug.eq.1) write (6,*)'nocon,ipass,dmax',nocon,ipass,dmax
c          if (debug.eq.1) write (6,*) 'd,ks', epsil(3:4)     
             
          if (nocon.eq.0) then
              dx = locp
              if (tnodes(3).eq.0) then
                  fifon(1) = fifon(1)+0.5d0*fifon(3)
                  fifon(2) = fifon(2)+0.5d0*fifon(3)
                  fifon(3) =ZERO
                  df(1) = df(1)+0.5d0*df(3)
                  df(2) = df(2)+0.5d0*df(3)
                  df(3) =ZERO
              endif
              if (tnodes(6).eq.0.and.locp.eq.-1) then
                  if (tnodes(4).gt.0) then
                      fifon(1) = fifon(1)+0.5d0*fifon(6)
                      fifon(4) = fifon(4)+0.5d0*fifon(6)
                      df(1) = df(1)+0.5d0*df(6)
                      df(4) = df(4)+0.5d0*df(6)
                  else
                      fifon(4) =ZERO
                      df(4) =ZERO
                  endif
                  fifon(6) =ZERO
                  df(6) =ZERO
              endif
              if (tnodes(7).eq.0.and.locp.eq.1) then
                  if (tnodes(5).gt.0) then
                      fifon(2) = fifon(2)+0.5d0*fifon(7)
                      fifon(5) = fifon(5)+0.5d0*fifon(7)
                      df(2) = df(2)+0.5d0*df(7)
                      df(5) = df(5)+0.5d0*df(7)
                  else
                      fifon(5) = ZERO
                      df(5) = ZERO
                  endif
                  fifon(7) = ZERO
                  df(7) = ZERO
              endif
c --- no curvature effect
              df(9) = ZERO
              df(10) = ZERO
          endif
          dum = fifon(2)           
          fifon(2) = fifon(3)
          fifon(3) = dum
          dum = df(2)           
          df(2) = df(3)
          df(3) = dum
c --- case of a quadrilateral segment     
      elseif (icode.eq.7) then
          do 240 ipi = 1,4
             dxdksi = ZERO
             dxdeta = ZERO
             dydksi = ZERO
             dydeta = ZERO
             dzdksi = ZERO
             dzdeta = ZERO
             do 230 inode = 1,4
                dxdksi = dxdksi+dfi(1,inode,ipi)*x(inode)
                dxdeta = dxdeta+dfi(2,inode,ipi)*x(inode)
                dydksi = dydksi+dfi(1,inode,ipi)*y(inode)
                dydeta = dydeta+dfi(2,inode,ipi)*y(inode)
                dzdksi = dzdksi+dfi(1,inode,ipi)*z(inode)
                dzdeta = dzdeta+dfi(2,inode,ipi)*z(inode)
 230         continue
             dtdx = (dydksi*dzdeta -dydeta*dzdksi)
             dtdy = (dzdksi*dxdeta -dzdeta*dxdksi)
             dtdz = (dxdksi*dydeta -dxdeta*dydksi)
             dajc = sqrt(dtdx*dtdx+dtdy*dtdy+dtdz*dtdz)
             nx(ipi) = dtdx/dajc
             ny(ipi) = dtdy/dajc
             nz(ipi) = dtdz/dajc
 240      continue
          
          call BeamToSurf (4,x(1),y(1),z(1),vx(1),vy(1),vz(1),
     &                     x0,y0,z0,vx0,vy0,vz0,nx(1),ny(1),nz(1),
     &                     r(1,1),epsil(1),fifon(1),nocon,dmax,iropil,
     &                     vxp,vyp,vzp,xpn1,ypn1,zpn1,omeg1,omeg2,omeg3,
     &                     delxy(1),isegda(11,iseg),ptols,gapshft,elem,
     &                     isegda(15,iseg),cnnode,ififo(12),ipass1,
     &                     rx(1),ry(1),rz(1),natcor(1))
      endif

c      BestU = epsil(4)
c      BestV = epsil(5)
c      if (nocon.lt.2 .or. abs(epsil(3)).lt.dmax) then
c          if (abs(BestU).le.(ONE+0.01d0*tols)) then
c              pBestU = BestU
c              pBestV = BestV
c              piseg = iseg
c          elseif ((pBestU**2+pBestV**2).gt.(BestU**2+BestV**2)) then
c              pBestU = BestU
c              pBestV = BestV
c              piseg = iseg
c          endif
c      endif
      
200   continue

c --- contact state
      if (abs(epsil(3)).gt.dmax) then
          epsil(3) = dmax
          if (nocon.eq.0) nocon = 2
      endif

c      if (nocon.eq.0 .and. ipass1.gt.0 .and. iseg000.gt.0) then
c          if (iseg.ne.iseg000 .and. try.eq.0) then
c              try = 1
c              nocon = 2
c              iseg = iseg000
c              goto 100
c          endif    
c      elseif (try.eq.1 .and. nocon.ne.0) then
c          try = 2
c          nocon = 2
c          iseg = iseg00
c          goto 100
c      elseif (try.eq.2 .and. nocon.ne.0) then
c          call erhandler ('CalPeni4',5000,2,'Error in CalPeni4.
c     &                     This should never happen.',ZERO,' ')
c      endif      
      
      
      if (ipass1.ne.-2.and.ickcee.eq.4) ickcee = 3
      if (nocon.eq.0) goto 999
      ickcee = 0
      if (ickcel.ge.3) then
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
      if (iseg.eq.iseg2) iseg2 = iseg - 1
      if (nocon.eq.2) goto 050
      if (icode.le.2) then
c --- beams
          if (nocon.eq.-1) iseg = isegda(11,iseg)
          if (nocon.eq. 1) iseg = isegda(12,iseg)
      else
c --- surface
          if (nocon.eq.11) iseg = isegda(11,iseg)
          if (nocon.eq.12) iseg = isegda(12,iseg)
          if (nocon.eq.13) iseg = isegda(13,iseg)
          if (nocon.eq.14) iseg = isegda(14,iseg)
      endif
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
      
      if (nocon.eq.0.and.iconnect.eq.1) cnnode = -cnnode

c      if (piseg.gt.0) then
c          roots(1) = pBestU
c          roots(2) = dble(piseg)
c      else
c          roots(1) = HUGE
c          roots(2) = ZERO
c      endif
c      iseg = piseg
c      write (8,*) X_DummyI
      
#if defined(DEBUG)
      call SubEnd('CalPeni4')
#endif

      return
      end