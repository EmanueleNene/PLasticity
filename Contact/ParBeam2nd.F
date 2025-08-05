      subroutine ParBeam2nd (rad,x,y,z,vx,vy,vz,x0,y0,z0,
     &                       vx0,vy0,vz0,r,ro,epsil,fifon,df,nocon,
     &                       dmax0,iropil,vxp,vyp,vzp,xpn,ypn,zpn,
     &                       omegx,omegy,omegz,delxy,isegda,tols,
     &                       gapshft,ipass,rx,ry,rz,locp,natcor,
     &                       BestU,hi,debug)
c --- primary function: compute the relations between an integration point
c                       of the contact element and a 3d line segment
#include "impcom.inc"
c 
      integer nocon,iropil,isegda(*),debug,key2d,hi
c
      double precision x0,y0,z0,vx0,vy0,vz0,dmax,epsil(8),fifon(5),
     &  vxp,vyp,vzp,xpn,ypn,zpn,omegx,omegy,omegz,delxy(3),gapshft,
     &  r(3,3),ro(3,3),x(5),y(5),z(5),vx(5),vy(5),vz(5),propen,BestU,
     &  tols,dmax0,rx(5),ry(5),rz(5),rad(2),df(10),df2(5),natcor(5),
     &  v_sum

      double precision ZERO,ONE,d,vpx,vpy,vpz,xmin,xmax,ymin,ymax,
     &   zmin,zmax,xs,ys,zs,con1,dis1,a(4),ks,xc,yc,zc,dis,d1,d2,adis,
     &   b(3),c,d0,a1x,a2x,a3x,a1y,a2y,a3y,a1z,a2z,a3z,roots(3),ks0,
     &   alpha,dir(3),norm,vect(6)
      integer iroot,ipass,numRoots,locp
c
      parameter (ZERO = 0.0d0, ONE = 1.0d0)
c
      key2d = 0
      if ((abs(z0)+abs(z(1))+abs(z(3))).lt.TINY) then
c --- 2-D logid
          key2d = 1
      endif
      dmax = dmax0
      ks = BestU
      if (ipass.eq.-2) then
          alpha = 0.2d0
          ks = natcor(1)
          locp = natcor(2)
          goto 702
      endif 
      if (ipass.eq.-1) then
          dmax = 5.0d0*dmax0
      endif
      nocon = 2
      d = dmax
      locp = 0
      delxy(1) = ZERO
      delxy(2) = ZERO
      delxy(3) = ZERO

c --- check out of pinball region
c 
      if (ipass.eq.-1) then
          nocon = 0
          locp = 0
          goto 10
      endif
      xmin = min(x(1),x(3)) - dmax
      if (x0.le.xmin) goto 999
      xmax = max(x(1),x(3)) + dmax
      if (x0.ge.xmax) goto 999
      ymin = min(y(1),y(3)) - dmax
      if (y0.le.ymin) goto 999
      ymax = max(y(1),y(3)) + dmax
      if (y0.ge.ymax) goto 999
      zmin = min(z(1),z(3)) - dmax
      if (z0.le.zmin) goto 999
      zmax = max(z(1),z(3)) + dmax
      if (z0.ge.zmax) goto 999
 10   continue

      a(1) = x(3) - x(1) 
      a(2) = y(3) - y(1) 
      a(3) = z(3) - z(1) 
      b(1) = x0 - 0.5d0*(x(3)+x(1))
      b(2) = y0 - 0.5d0*(y(3)+y(1))
      b(3) = z0 - 0.5d0*(z(3)+z(1))
      c = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
      d0 = 2.0d0*(a(1)*b(1) + a(2)*b(2) + a(3)*b(3))
      ks = d0/c
      if (ipass.eq.-1) then
          nocon = 0
          locp = 0
          alpha = 0.2d0
          if (abs(ks).le.(ONE-2.d0*alpha).or.abs(alpha).le.TINY)goto 702
          if (ks.lt.ZERO) then
              if (isegda(1).eq.0) goto 702
              locp = -1
          else
              if (isegda(2).eq.0) goto 702
              locp = 1
          endif
          goto 702
      endif
      if (hi.eq.1) then
          alpha = 0.01d0*tols
      else
          alpha = ZERO
      endif
      con1 = ONE + alpha + 1.d-8
      if (ks.gt.con1) then
          nocon = 1
          goto 900
      endif
      if (ks.lt.-con1) then
          nocon = -1
          goto 900
      endif
      locp = 0
      if (abs(ks).le.(ONE-2.0d0*alpha).or.abs(alpha).le.TINY) then
          goto 702
      else
          ks0 = ks
          if (ks.lt.ZERO) then
              if (isegda(1).eq.0) goto 702
c --- left side
              a1x = x(4) + x(3) - 2.d0*x(1)
              a2x = x(3) - x(4) 
              a3x = x(1) - x0
              a1y = y(4) + y(3) - 2.d0*y(1)
              a2y = y(3) - y(4) 
              a3y = y(1) - y0
              a1z = z(4) + z(3) - 2.d0*z(1)
              a2z = z(3) - z(4) 
              a3z = z(1) - z0
              locp = -1
          else
              if (isegda(2).eq.0) goto 702
c --- right side
              a1x = x(1) + x(5) - 2.d0*x(3)
              a2x = x(5) - x(1) 
              a3x = x(3) - x0
              a1y = y(1) + y(5) - 2.d0*y(3)
              a2y = y(5) - y(1) 
              a3y = y(3) - y0
              a1z = z(1) + z(5) - 2.d0*z(3)
              a2z = z(5) - z(1) 
              a3z = z(3) - z0
              locp = 1
          endif
          a(1) = alpha*(a1x*a1x+a1y*a1y+a1z*a1z)
          a(2) = 3.d0*alpha*(a1x*a2x+a1y*a2y+a1z*a2z)
          a(3) = 2.d0*alpha*(a2x*a2x+a2y*a2y+a2z*a2z)+a(1)+
     &           4.d0*(a3x*a1x+a3y*a1y+a3z*a1z)
          a(4) = a(2)/3.d0+4.d0*(a3x*a2x+a3y*a2y+a3z*a2z)
          call FindRoot3K (a(1), roots(1), numRoots)
          ks = roots(1)

c          if (debug.eq.1) write (8,*) 'ru',roots(1:numRoots)
              
          do iroot = 2,numRoots
             if (abs(ks).gt.abs(roots(iroot))) ks = roots(iroot)
          enddo

c          if (debug.eq.1) write (8,*) 'ks',ks,ks0

          con1 = ONE + 1.d-8
          if (abs(ks).gt.con1) then
c --- if we go over the line, reduce number of found contacts
              d = dmax
              goto 900
          endif
          
c          if (locp.eq.-1.and.ks.le.-1.d0) then
c              locp = 0
c              ks = ks0
c              goto 702
cc              nocon = -1
cc              goto 900
c          endif
c          if (locp.eq.1.and.ks.gt.1.d0) then
cc              locp = 0
cc              ks = ks0
cc              goto 702
c              nocon = 1
c              goto 900
c          endif
      endif
 702  continue
      nocon = 0
      call v_zero (df2(1),5)
      if (locp.eq.0) then
          fifon(1) = 0.5d0*(ONE-ks)
          fifon(2) = ZERO
          fifon(3) = 0.5d0*(ONE+ks)
          fifon(4) = ZERO
          fifon(5) = ZERO
          df(1) = -0.5d0
          df(2) = ZERO
          df(3) = 0.5d0
          df(4) = ZERO
          df(5) = ZERO
      elseif (locp.eq.-1) then
          fifon(1) = ONE - 0.5d0*alpha*(ks*ks+ONE)
          fifon(2) = ZERO
          fifon(3) = 0.25d0*alpha*(ks*ks+2.d0*ks+ONE)
          fifon(4) = 0.25d0*alpha*(ks*ks-2.d0*ks+ONE)
          fifon(5) = ZERO
          df(1) = -alpha*ks
          df(2) = ZERO
          df(3) = 0.5d0*alpha*(ks+ONE)
          df(4) = 0.5d0*alpha*(ks-ONE)
          df(5) = ZERO
          df2(1) = -alpha
          df2(3) = 0.5d0*alpha
          df2(4) = 0.5d0*alpha
      else
          fifon(1) = 0.25d0*alpha*(ks*ks-2.d0*ks+ONE)
          fifon(2) = ZERO
          fifon(3) = ONE - 0.5d0*alpha*(ks*ks+ONE)
          fifon(4) = ZERO
          fifon(5) = 0.25d0*alpha*(ks*ks+2.d0*ks+ONE)
          df(1) = 0.5d0*alpha*(ks-ONE)
          df(2) = ZERO
          df(3) = -alpha*ks
          df(4) = ZERO
          df(5) = 0.5d0*alpha*(ks+ONE)
          df2(1) = 0.5d0*alpha
          df2(3) = -alpha
          df2(5) = 0.5d0*alpha
      endif
      natcor(1) = ks
      natcor(2) = locp
 
      xc = fifon(1)*x(1)+fifon(3)*x(3)+fifon(4)*x(4)+fifon(5)*x(5)
      yc = fifon(1)*y(1)+fifon(3)*y(3)+fifon(4)*y(4)+fifon(5)*y(5)
      zc = fifon(1)*z(1)+fifon(3)*z(3)+fifon(4)*z(4)+fifon(5)*z(5)
      natcor(3) = xc
      natcor(4) = yc
      natcor(5) = zc
      dis = (x0-xc)*(x0-xc)+(y0-yc)*(y0-yc)+(z0-zc)*(z0-zc)
      dis = sqrt(dis)

      if (rad(1).lt.ZERO) then
          call v_zero(dir(1),3)
          dis1 = abs(abs(rad(1))-abs(rad(2))) - dis
      else
c --- save the current normal for later checking
          if (key2d.eq.1) then
              dir(1) = y(1) - y(3)
              dir(2) = -(x(1) - x(3))
              dir(3) = ZERO
          else
              dir(1) = ro(3,1)
              dir(2) = ro(3,2)
              dir(3) = ro(3,3) 
          endif
      endif

      a(1) = df(1)*x(1)+df(3)*x(3)+df(4)*x(4)+df(5)*x(5)
      a(2) = df(1)*y(1)+df(3)*y(3)+df(4)*y(4)+df(5)*y(5)
      a(3) = df(1)*z(1)+df(3)*z(3)+df(4)*z(4)+df(5)*z(5)
c --- calculate contact normal and tangent trasformatin
      if (dis.gt.1.0d-12*dmax) then 
          if (rad(1).lt.ZERO) then 
              r(3,1) = (x0-xc)/dis
              r(3,2) = (y0-yc)/dis
              r(3,3) = (z0-zc)/dis
          else
              r(3,1) = -(x0-xc)/dis
              r(3,2) = -(y0-yc)/dis
              r(3,3) = -(z0-zc)/dis
c --- check on the normal change
              norm = dir(1)*r(3,1)+dir(2)*r(3,2)+dir(3)*r(3,3)
              if (norm.lt.ZERO) then
c --- reverse normal
                  r(3,1) = -r(3,1)
                  r(3,2) = -r(3,2)
                  r(3,3) = -r(3,3)
                  dis = -dis
              endif
          endif
      elseif ((dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3)).gt.
     &         sqrt(TINY)) then
c --- the normal will be retained
c --- check on the normal change
          if (rad(1).lt.ZERO) then 
              r(3,1) = (x0-xc)
              r(3,2) = (y0-yc)
              r(3,3) = (z0-zc)
          else
              r(3,1) = -(x0-xc)
              r(3,2) = -(y0-yc)
              r(3,3) = -(z0-zc)
          endif
          norm = dir(1)*r(3,1)+dir(2)*r(3,2)+dir(3)*r(3,3)
          if (norm.lt.ZERO) then
              dis = -dis
          endif
          r(3,1) = dir(1)
          r(3,2) = dir(2)
          r(3,3) = dir(3)
      else
          a(1) = df(1)*x(1)+df(3)*x(3)+df(4)*x(4)+df(5)*x(5)
          a(2) = df(1)*y(1)+df(3)*y(3)+df(4)*y(4)+df(5)*y(5)
          a(3) = df(1)*z(1)+df(3)*z(3)+df(4)*z(4)+df(5)*z(5)
          vect(1) = x(1)
          vect(2) = y(1)
          vect(3) = z(1)
          vect(4) = x(1)+a(1)
          vect(5) = y(1)+a(2)
          vect(6) = z(1)+a(3)
          call LineGeomTrans (vect(1),3,r(1,1))
          goto 200     
      endif
      adis = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      r(1,1) = a(1)/adis
      r(1,2) = a(2)/adis
      r(1,3) = a(3)/adis
      r(2,1) = r(3,2)*r(1,3) - r(3,3)*r(1,2)
      r(2,2) = r(3,3)*r(1,1) - r(3,1)*r(1,3)
      r(2,3) = r(3,1)*r(1,2) - r(3,2)*r(1,1)
 200  continue 

      if (rad(1).ge.ZERO) then
          dis1 = dis - rad(1) - abs(rad(2))
      endif

      if (abs(dis1).lt.dmax) then
          nocon = 0
          d = dis1
      else
          nocon = 2
          d = dmax
          goto 999
      endif

c --- segment curvature
c     
      a(1) = df2(1)*x(1)+df2(3)*x(3)+df2(4)*x(4)+df2(5)*x(5)
      a(2) = df2(1)*y(1)+df2(3)*y(3)+df2(4)*y(4)+df2(5)*y(5)
      a(3) = df2(1)*z(1)+df2(3)*z(3)+df2(4)*z(4)+df2(5)*z(5)
      adis = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
      if (adis.gt.TINY) then
          df(9) = -(r(3,1)*a(1)+r(3,2)*a(2)+r(3,3)*a(3))/adis
          df(10) = -(r(1,1)*a(1)+r(1,2)*a(2)+r(1,3)*a(3))/adis
      else
          df(9) = ZERO
          df(10) = ZERO
      endif
c --- relative velocity
c
      if (iropil.le.0) then
          vpx = fifon(1)*vx(1)+fifon(3)*vx(3)+
     &          fifon(4)*vx(4)+fifon(5)*vx(5)
          vpy = fifon(1)*vy(1)+ fifon(3)*vy(3)+
     &          fifon(4)*vy(4)+fifon(5)*vy(5)
          vpz = fifon(1)*vz(1) + fifon(3)*vz(3)+
     &          fifon(4)*vz(4)+fifon(5)*vz(5)
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
c
      epsil(1) = vpx*r(1,1) + vpy*r(1,2) + vpz*r(1,3)
      epsil(2) = vpx*r(2,1) + vpy*r(2,2) + vpz*r(2,3)
      epsil(6) = rx(1)*fifon(1) + rx(3)*fifon(3) +
     &           rx(4)*fifon(4) + rx(4)*fifon(4)
      epsil(7) = ry(1)*fifon(1) + ry(3)*fifon(3) + 
     &           ry(4)*fifon(4) + ry(4)*fifon(4)
      epsil(8) = rz(1)*fifon(1) + rz(3)*fifon(3) +
     &           rz(4)*fifon(4) + rz(4)*fifon(4)
      delxy(1) = r(3,1)
      delxy(2) = r(3,2)
      delxy(3) = r(3,3)
      d = d - gapshft

      if (rad(1).ge.ZERO) then
          dis1 = gapshft -rad(1)-abs(rad(2))
      else
          dis1 = gapshft + abs(abs(rad(1))-abs(rad(2)))
      endif

      if (ipass.eq.-2) then
          natcor(3) = natcor(3) - r(3,1)*dis1
          natcor(4) = natcor(4) - r(3,2)*dis1
          natcor(5) = natcor(5) - r(3,3)*dis1
      endif
      goto 999
c
 900  continue
c
c --- select contact status
      con1 = ONE + 1.e-8
      if (abs(d).lt.0.1d0*dmax.and.d.lt.sqrt(TINY)) 
     &    con1 = ONE + 1.0d-5*tols
      if (ipass.ge.1) then
          con1 = ONE + tols*10.d0**(ipass-6)
      endif
      propen = ONE + 0.01d0*tols

c      if (abs(ks).gt.propen) then
c          nocon = 2
c          goto 999
c      endif 

      if (ks.gt.con1) then
          nocon = 1
          if (isegda(2).ne.0.or.d.gt.ZERO) goto 950
          if (abs(ks).gt.propen) goto 950
      elseif (ks.lt.-con1) then
          nocon = -1
          if (isegda(1).ne.0.or.d.gt.ZERO) goto 950
          if (abs(ks).gt.propen) goto 950
      endif
      goto 702
 950  continue
 
      if (nocon.eq.0) goto 999
      if (abs(d).lt.0.1d0*dmax) goto 999
      if (abs(ks).lt.2.0d0) goto 999
      d1 = (x(1)-x0)*(x(1)-x0)+(y(1)-y0)*(y(1)-y0)+(z(1)-z0)*(z(1)-z0)
      d2 = (x(3)-x0)*(x(3)-x0)+(y(3)-y0)*(y(3)-y0)+(z(3)-z0)*(z(3)-z0)
      if (d1.ge.d2) then
          d = sqrt(d2)
          if (d.ge.dmax) then
              d = dmax
              goto 999
          endif
          if (d.gt.sqrt(TINY)*dmax) then
              delxy(1) = (x0-x(3))/d
              delxy(2) = (y0-y(3))/d
              delxy(3) = (z0-z(3))/d
          endif
      else
          d = sqrt(d1)
          if (d.ge.dmax) then
              d = dmax
              goto 999
          endif
          if (d.gt.sqrt(TINY)*dmax) then
              delxy(1) = (x0-x(1))/d
              delxy(2) = (y0-y(1))/d
              delxy(3) = (z0-z(1))/d
         endif
      endif
 
 999  continue
      epsil(3) = d
      epsil(4) = ks
      epsil(5) = ZERO
c      if (debug.eq.1 .and. nocon.ne.2) write (6,*) 'ks_fin',ks, d, nocon
          
      return
      end

