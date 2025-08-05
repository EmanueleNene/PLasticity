      subroutine CalPeni04 (x0,y0,z0,r,dmax,nsegt,ifoun,npers,
     &                      nPerFounI,nperfd,ififo,coor0,dispn,isegda,
     &                      fcell,founda,nocon,iseg,epsil,delxy,tols,
     &                      gapshft,elem,cnnode,ibeam,hi)
c --- primary function: compute the initial relations between an integration point
c                       of the contact element and the corresponding foundation
c                     
#include "impcom.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer nsegt,ifoun,npers,nPerFounI,nperfd,nocon,iseg,
     &  isegda(npers,nsegt),ififo(nPerFounI),ibeam,
     &  ipi,i,j,k,hi
      double precision x0,y0,z0,dmax,fifon(8),epsil(8),coor0(*),
     &  dispn(*),r(3,3),fcell(*),delxy(3),founda(nperfd)

      integer is,icode,i1,kddlx,kddly,kddlz,ickcel,iseg1,iseg2,is0,
     &  inode,nnode,elem,cnnode,contactNode,tnodes(7),n1,n2,
     &  iseg0,iseg00,ipass1,iconnect,locp,trg3dbeam,jbeam
     
      double precision ZERO,HALF,d0,x(7),y(7),z(7),t(7),tols,natcor(5),
     &  nx(7),ny(7),nz(7),gapshft,ro(3,3),ptols,
     &  dfi(2,9,4),df(10),rad(2),GetFibRad,
     &  dxdksi,dxdeta,dydksi,dydeta,dzdksi,dzdeta,dtdx,dtdy,dtdz,dajc

      parameter (ZERO=0.0d0,HALF=0.5d0)
      
      data  (((dfi(j,i,k),j=1,2),i=1,4),k=1,4) /
     & -0.5d0,-0.5d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0,
     & -0.5d0, 0.0d0, 0.5d0,-0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0,
     &  0.0d0, 0.0d0, 0.0d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0, 0.0d0,
     &  0.0d0,-0.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0,-0.5d0, 0.5d0/

      
c --- initialize
      call v_zero (epsil(1),8)
      call v_zero (ro(1,1),9)
      d0 = dmax
      epsil(3) = dmax
      is0 = 0

      rad(1) = GetFibRad(ibeam)
      if (rad(1).le.TINY) rad(1) = founda(7)
      
      if (ifoun.le.0) goto 070

      ickcel = 0
      iconnect = 0
c --- check the cell region
      if (x0.lt.fcell(1) .or. x0.gt.fcell(2)) ickcel = 1
      if (y0.lt.fcell(3) .or. y0.gt.fcell(4)) ickcel = 1
      if (z0.lt.fcell(5) .or. z0.gt.fcell(6)) ickcel = 1
      if (ickcel .eq. 1) goto 070

c --- contact node
      contactNode = 0
      if (cnnode.gt.0) then
          contactNode = cnnode
          cnnode = (cnnode-1)/D_NumDof + 1
          cnnode = M_BackList(cnnode)
      endif

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
          if (isegda(18,iseg).eq.0.and.isegda(1,iseg).ne.99)
     &        goto 100
      endif
      is0 = iseg1      

  050 continue
      iseg00 = 0
      ipass1 = 0
      do 060 is = iseg1,iseg2
         if (isegda(18,is).ne.0) goto 060
         if (isegda(1,is).eq.99) goto 060 
         if (elem.eq.isegda(17,is)) then
             isegda(18,is) = 1
             goto 060
         endif
         iseg = is
         goto 100
  060 continue

  070 iseg = is0
      nocon = 2
      epsil(3) = dmax
      goto 999

 100  continue
      iseg00 = iseg0
      iseg0 = iseg
      icode = isegda(1,iseg)
      isegda(18,iseg) = 1
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

      call v_zero (t(1),7)

      iconnect = 0
      do 110 inode = 1,nnode
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
             else
                 x(inode) = ZERO
                 y(inode) = ZERO
                 z(inode) = ZERO
             endif
             goto 110 
         endif
c --- contact node can not contact its own elements
         if (contactNode.gt.0) then
             if (i1.eq.contactNode) then
                 nocon = 2
                 epsil(3) = d0
                 iconnect = 1
                 goto 200
             endif
         endif
         kddlx = i1
         kddly = kddlx + 1
         kddlz = kddlx + 2
         x(inode) = dispn(kddlx) + coor0(kddlx) + founda(13)
         y(inode) = dispn(kddly) + coor0(kddly) + founda(14)
         z(inode) = dispn(kddlz) + coor0(kddlz) + founda(15)
  110 continue

      d0 = epsil(3)

      if (icode.eq.1 .or. icode.eq.2) then
          jbeam = abs(isegda(16,iseg))
          rad(2) = GetFibRad(jbeam)
          if (rad(2).le.TINY) rad(2) = founda(8)
      endif
      ptols = tols*0.8d0
      if (icode.eq.1) then
c --- 3D line segment
          dmax = dmax + abs(rad(1)) + abs(rad(2))
          call ParBeam2nd (rad(1),x(1),y(1),z(1),t(1),t(1),t(1),
     &                     x0,y0,z0,ZERO,ZERO,ZERO,r(1,1),ro(1,1),
     &                     epsil(1),fifon(1),df(1),nocon,dmax,0,
     &                     ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     &                     delxy(1),isegda(11,iseg),ptols,gapshft,
     &                     ipass1,t(1),t(1),t(1),locp,natcor(1),HUGE,
     &                     hi,0)
      elseif (icode.eq.2) then
c --- 3D parabolic segment
          dmax = dmax + abs(rad(1)) + abs(rad(2))
          call ParBeam3nd (rad(1),x(1),y(1),z(1),t(1),t(1),t(1),
     &                     x0,y0,z0,ZERO,ZERO,ZERO,r(1,1),ro(1,1),
     &                     epsil(1),fifon(1),df(1),nocon,dmax,0,
     &                     ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     &                     delxy(1),isegda(11,iseg),ptols,gapshft,
     &                     ipass1,t(1),t(1),t(1),locp,natcor(1),HUGE,
     &                     hi,0)
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
  240     continue
          call BeamToSurf (4,x(1),y(1),z(1),t(1),t(1),t(1),
     &                     x0,y0,z0,ZERO,ZERO,ZERO,nx(1),ny(1),nz(1),
     &                     r(1,1),epsil(1),fifon(1),nocon,dmax,0,
     &                     ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     &                     delxy(1),isegda(11,iseg),ptols,gapshft,elem,
     &                     isegda(15,iseg),cnnode,ififo(12),ipass1,
     &                     t(1),t(1),t(1),natcor(1))
      endif
 
 200  continue

c --- contact state
      if (epsil(3).gt.dmax) then 
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
      if (nocon.eq.0.and.iconnect.eq.1) cnnode = -cnnode
      return
      end

