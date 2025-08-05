      subroutine BeamToSurf (nfacn,x,y,z,vx,vy,vz,x0,y0,z0,
     &  vx0,vy0,vz0,nx,ny,nz,r,epsil,fifon,nocon,dmax0,iropil,
     &  vxp,vyp,vzp,xpn,ypn,zpn,omegx,omegy,omegz,delxy,isegda,
     &  tols,gapshft,elContact,elTarget,contactNode,kerr,ipass,
     &  rx,ry,rz,natcor)
c --- primary function: compute the relations between an integration point
c                       and a 4-node quadrilateral surface element
c
c  input arguments:
c
c      x0 ,y0 ,z0  (dp,sc,in) - global coordinates of the integration point
c      vx0,vy0vz0  (dp,sc,in) - components of the velocity of the integration 
c                               point in global cartesian coordinates
c      r (dp,ar(3,3),in)      - components of the normal vector to the 
c                               structure
c      dmax        (dp,sc,in) - the maximum distance for false contact
c                               (if abs(d).gt.dmax, false contact is assumed)
c      iropil   (int,sc,in)   - rotaion option
c                                = 0 no rotation
c                                = 1 with rotation
c      xpn      (int,sc,in)   - pilot node loaction x componemt
c      ypn      (int,sc,in)   - pilot node loaction y componemt
c      zpn      (int,sc,in)   - pilot node loaction z componemt
c      vxp      (int,sc,in)   - pilot node velocity x componemt
c      vyp      (int,sc,in)   - pilot node velocity y componemt
c      vzp      (int,sc,in)   - pilot node velocity z componemt
c      omegx    (int,sc,in)   - pilot node angle velocity x componemt
c      omegy    (int,sc,in)   - pilot node angle velocity y componemt
c      omegz    (int,sc,in)   - pilot node angle velocity z componemt
c
c  output arguments:
c
c        epsil    (dp,ar(4),inout)
c
c        epsil(1)   = tangent velocity of structure in first dir.
c        epsil(2)   = tangent velocity of structure in second dir.
c        epsil(3)    
c                d  = distance between the integration point
c                     and the foundation
c                     ( d < 0 : contact    )
c                     ( d > 0 : no contact )
c
c        fifon (dp,ar(8),out) - shape functions of the foundation
c        nocon = 0  : if there is probably contact 
c                     ( normal coherence but no test on the distance)
c              = 11 : if there is no contact but : 
c                     eta < -1   and 1 > ksi > -1
c              = 12 : if there is no contact but : 
c                     ksi > 1    and 1 > eta > -1
c              = 13 : if there is no contact but :    
c                     eta > 1    and 1 > ksi > -1
c              = 14 : if there is no contact but :    
c                     ksi < -1   and 1 > eta > -1
c
c       delxy    (dp,ar(3),out) - distance between integartion point to target node
c         delxy(1)          - x component
c         delxy(2)          - y component
c         delxy(3)          - z component
c
#include "impcom.inc"
c
      integer nfacn,nocon,iropil,inors,isegda(*),
     &        elContact,elTarget,contactNode,kerr
c
      double precision x0,y0,z0,vx0,vy0,vz0,dmax,epsil(8),fifon(8),
     &  vxp,vyp,vzp,xpn,ypn,zpn,omegx,omegy,omegz,delxy(3),gapshft,
     &  ContCheck4,r(3,3),x(8),y(8),z(8),vx(8),vy(8),vz(8),off(3),
     &  propen,tols,nx(8),ny(8),nz(8),dmax0,
     &  rx(8),ry(8),rz(8),natcor(5)
c
      double precision ZERO,ONE,d,vpx,vpy,vpz,QUART,norm,
     &   xmin,xmax,ymin,ymax,zmin,zmax,dd,dmin,dpres,
     &   root1(65),root2(65),tmpSolutions(65),bestU,
     &   bestV,bestD,xx,yy,zz,xs,ys,zs,d0(8),con1,vn(3),
     &   vn1,vn2,vn3,prec,dperr(5)
      integer numRoot,i,j,tmpNumSolutions,found,iroot,itry,ipass
c
      parameter (QUART=0.25d0, ZERO=0.0d0, ONE=1.0d0)

      itry = 0
      dmax = dmax0
      nocon = 2
      d = dmax
      bestU = HUGE
      bestV = HUGE

c --- check out of pinball region
c 
      xmin = x(1)
      ymin = y(1)
      zmin = z(1)
      xmax = x(1)
      ymax = y(1)
      zmax = z(1)
      do i =2,nfacn
         if (x(i).lt.xmin) xmin = x(i)
         if (x(i).gt.xmax) xmax = x(i)
         if (y(i).lt.ymin) ymin = y(i)
         if (y(i).gt.ymax) ymax = y(i)
         if (z(i).lt.zmin) zmin = z(i)
         if (z(i).gt.zmax) zmax = z(i)
      enddo
      xmin = xmin - dmax
      if (x0.le.xmin) goto 999
      xmax = xmax + dmax
      if (x0.ge.xmax) goto 999
      ymin = ymin - dmax
      if (y0.le.ymin) goto 999
      ymax = ymax + dmax
      if (y0.ge.ymax) goto 999
      zmin = zmin - dmax
      if (z0.le.zmin) goto 999
      zmax = zmax + dmax
      if (z0.ge.zmax) goto 999

 10   continue

      off(1) = x0
      off(2) = y0
      off(3) = z0
      vn(1) = 0
      vn(2) = 0
      vn(3) = 0
      do i = 1,nfacn
         norm = sqrt(nx(i)*nx(i)+ny(i)*ny(i)+nz(i)*nz(i))
         nx(i) = nx(i)/norm
         ny(i) = ny(i)/norm
         nz(i) = nz(i)/norm
         vn(1) = vn(1) + nx(i)
         vn(2) = vn(2) + ny(i)
         vn(3) = vn(3) + nz(i)
      enddo

c --- select non-zero normal
      vn1 = vn(2)*vn(2)+vn(3)*vn(3)
      vn2 = vn(1)*vn(1)+vn(3)*vn(3)
      vn3 = vn(1)*vn(1)+vn(2)*vn(2)
      inors = 1
      if (vn2.gt.norm) then
          inors  = 2
          norm = vn2
      endif
      if (vn3.gt.norm) then
          inors = 3
      endif

 100  continue
      d = dmax
      nocon = 2
      call BeamToSurf_1 (nfacn,x(1),y(1),z(1),off(1),nx(1),ny(1),nz(1),
     &                   inors,itry,root1(1),root2(1),numRoot,ipass)

      tmpNumSolutions = 0

      do 550 j=1,numRoot
         if (ContCheck4 (nfacn,root1(j),root2(j),nx(1),ny(1),nz(1),
     &       x(1),y(1),z(1),dmax,x0,y0,z0,dd,ipass).ge.ZERO) goto 550
         tmpSolutions(tmpNumSolutions*3+1) = root1(j)
         tmpSolutions(tmpNumSolutions*3+2) = root2(j)
         tmpSolutions(tmpNumSolutions*3+3) = dd
         tmpNumSolutions = tmpNumSolutions + 1
 550  continue     
      
      if (tmpNumSolutions.le.0) goto 999

c --- find the closest solution
      bestD = dmax
      iroot = 0
      con1 = ONE + 1.d-8
      do 700 i=0, tmpNumSolutions-1
         if (abs(tmpSolutions(i*3+1)).gt.con1 .or. 
     &       abs(tmpSolutions(i*3+2)).gt.con1) goto 700
         dd = tmpSolutions(i*3+3)
         found = 0
         if (dd.le.ZERO.and.bestD.gt.ZERO) found = 1
         if ((dd.le.ZERO.and.bestD.le.ZERO).or.
     &       (dd.ge.ZERO.and.bestD.ge.ZERO)) then
            if (abs(dd).lt. bestD) found = 1
         endif                            
         if (found.eq.1) then
             bestU = tmpSolutions(i*3+1) 
             bestV = tmpSolutions(i*3+2) 
             bestD = dd
         endif
         if (found.eq.1) iroot = 1
 700  continue
 701  continue
      if (itry.eq.2) then
          call SelectRoot4 (x(1),y(1),z(1),x0,y0,z0,vn(1),nx(1),ny(1),
     &                     nz(1),dmax,bestU,bestV,bestD,iroot,ipass)
          if (iroot.eq.1) then
              tmpSolutions(1) = bestU
              tmpSolutions(2) = bestV
              tmpSolutions(3) = bestD
              tmpNumSolutions = 1
          else
              d = dmax
              nocon = 2
              goto 999
          endif
      endif
      if (iroot.eq.0) goto 800
      if (abs(bestD).ge.dmax) goto 800

 702  continue
      nocon = 0
      d = bestD

c --- interpolation function of contact
c
      fifon(1)= QUART*(ONE-bestU)*(ONE-bestV)
      fifon(2)= QUART*(ONE+bestU)*(ONE-bestV)
      fifon(3)= QUART*(ONE+bestU)*(ONE+bestV)
      fifon(4)= QUART*(ONE-bestU)*(ONE+bestV)
      natcor(1) = bestU
      natcor(2) = bestV

c --- verify the accuracy of the intersection
      xx  = ZERO
      yy  = ZERO
      zz  = ZERO
      vn(1) = ZERO
      vn(2) = ZERO
      vn(3) = ZERO
      do 710 i = 1,nfacn
         xx = xx + fifon(i)*x(i)
         yy = yy + fifon(i)*y(i)
         zz = zz + fifon(i)*z(i)
         vn(1) = vn(1) + fifon(i)*nx(i)
         vn(2) = vn(2) + fifon(i)*ny(i)
         vn(3) = vn(3) + fifon(i)*nz(i)
 710  continue

      natcor(3) = xx
      natcor(4) = yy
      natcor(5) = zz
      vpx = sqrt(vn(1)*vn(1)+vn(2)*vn(2)+vn(3)*vn(3))     
      vn(1) = vn(1)/vpx
      vn(2) = vn(2)/vpx
      vn(3) = vn(3)/vpx
      xx  = xx + bestD*vn(1)
      yy  = yy + bestD*vn(2)
      zz  = zz + bestD*vn(3)
      vpx = (xx-x0)*(xx-x0)+(yy-y0)*(yy-y0)+(zz-z0)*(zz-z0)
      prec = 1.e-6
      if (bestD.gt.0) prec = 10.d0*prec 
      if (vpx.gt.prec*dmax*dmax) then
         if (itry.eq.0) then
             itry = 1
             goto 100
         endif
         if (itry.eq.1) then
             itry = 2
             goto 701
         endif
         dperr(1) = bestD
         dperr(2) = contactNode
         dperr(3) = elContact
         dperr(4) = elTarget
         dperr(5) = kerr
         call erhandler ('BeamToSurf',5000,2,
     &      'Gap(+)/penetration(-) %g  between contact node %i of 
     &      element %i and target element %i (real ID %i) is
     &      inaccurate',dperr(1),' ')
      endif


c --- calculate contact normal and tangent trasformation
      call CalInvTr (nfacn,x(1),y(1),z(1),bestU,bestV,vn(1),r(1,1))

c --- relative velocity
c
      if (iropil.le.0) then
          vpx = ZERO
          vpy = ZERO
          vpz = ZERO
          do 720 i = 1,nfacn
             vpx = vpx + fifon(i)*vx(i)
             vpy = vpy + fifon(i)*vy(i)
             vpz = vpz + fifon(i)*vz(i)
  720     continue
      else
          xs = x0 + d*r(3,1)
          ys = y0 + d*r(3,2)
          zs = z0 + d*r(3,3)
          vpx = vxp - omegz*(ys-ypn) + omegy*(zs-zpn)
          vpy = vyp - omegx*(zs-zpn) + omegz*(xs-xpn)
          vpz = vzp - omegy*(xs-xpn) + omegx*(ys-ypn)
      endif
      vpx  = vx0 - vpx
      vpy  = vy0 - vpy
      vpz  = vz0 - vpz

      epsil(1)= vpx*r(1,1) + vpy*r(1,2) + vpz*r(1,3)
      epsil(2)= vpx*r(2,1) + vpy*r(2,2) + vpz*r(2,3)

c --- rotations
c
      epsil(6) = ZERO
      epsil(7) = ZERO
      epsil(8) = ZERO
      do 730 i = 1,nfacn
         epsil(6) = epsil(6) + rx(i)*fifon(i)
         epsil(7) = epsil(7) + ry(i)*fifon(i)
         epsil(8) = epsil(8) + rz(i)*fifon(i)
  730 continue
      goto 900
  800 continue

c --- if the search in the patch fails but there are some intersections
c     found, we search the extended patch.
      bestD = dmax
      do 750 i=0, tmpNumSolutions-1
         dd = tmpSolutions(i*3+3)
         if (abs(dd).lt.bestD) then
             found = 1
             bestU = tmpSolutions(i*3+1) 
             bestV = tmpSolutions(i*3+2) 
             bestD = dd
         endif
 750  continue
      if (found.eq.0) goto 999
      if (abs(bestD).ge.dmax) goto 999
      d = bestD
      if (ipass.eq.-1) then
          nocon = 0
          goto 702
      endif

c --- select contact status
      con1 = ONE + 1.e-8
      if (abs(bestD).lt.0.1d0*dmax.and.bestD.lt.sqrt(TINY)) 
     &    con1 = ONE + 1.d-5*tols
      if (ipass.ge.1) then
          con1 = ONE + tols*10.d0**(ipass-6)
      endif
      propen = ONE + 0.01d0*tols
      if (bestV.lt.-con1.and.abs(bestV).ge.abs(bestU)) then
          nocon = 11
          if (isegda(1).ne.0) goto 900
          if (abs(bestV).gt.propen.or.abs(bestU).gt.propen) goto 900
      endif
      if (bestU.gt.con1.and.abs(bestU).ge.abs(bestV)) then
          nocon = 12
          if (isegda(2).ne.0) goto 900
          if (abs(bestV).gt.propen.or.abs(bestU).gt.propen) goto 900
      endif
      if (bestV.gt.con1.and.abs(bestV).ge.abs(bestU)) then
          nocon = 13
          if (isegda(3).ne.0) goto 900
          if (abs(bestV).gt.propen.or.abs(bestU).gt.propen) goto 900
      endif
      if (bestU.lt.-con1.and.abs(bestU).ge.abs(bestV)) then
          nocon = 14
          if (isegda(4).ne.0) goto 900
          if (abs(bestV).gt.propen.or.abs(bestU).gt.propen) goto 900
      endif

      goto 702
c
 900  continue
      delxy(1)= vn(1)
      delxy(2)= vn(2)
      delxy(3)= vn(3)
      d = d - gapshft

      if (nocon.eq.0) goto 999
      if (abs(d).lt.0.1*dmax) goto 999
      if (abs(bestV).lt.2.0d0.and.abs(bestU).lt.2.0d0) goto 999

      do 910 i = 1,nfacn
         d0(i) = (x(i)-x0)*(x(i)-x0) + 
     &           (y(i)-y0)*(y(i)-y0) +
     &           (z(i)-z0)*(z(i)-z0)
 910  continue

      dmin = d0(1)
      do i = 2,nfacn
         dmin = min(dmin,d0(i))
      enddo

      dpres = dmax*1.d-9

      do 920 i = 1,nfacn
         if (abs(d0(i)-dmin).le.dpres) then
             d = sqrt(d0(i))
             if (d.ge.dmax) then
                 d = dmax
                 goto 999
             endif
             if (d.lt.TINY*dmax) goto 999
             delxy(1)= (x0-x(i))/d
             delxy(2)= (y0-y(i))/d
             delxy(3)= (z0-z(i))/d
             goto 999
         endif
 920  continue
 
 999  continue

      epsil(3) = d
      epsil(4) = bestU
      epsil(5) = bestV

      return
      end
