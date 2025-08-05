      subroutine BeamToBeam (rad,x0,y0,z0,xa,ya,za,
     &                       vax,vay,vaz,xb,yb,zb,vbx,vby,vbz,
     &                       r,ro,epsil,fifoa,fifob,nocon,dmax0,
     &                       delxy,isegda,tols,gapshft,ipass,
     &                       rx,ry,rz,locu,locv,numRoot,
     &                       pBestU,pBestV,pBestD,a,b,ac,bc,
     &                       debug)
c --- primary function: compute the relations between beam-beam
c
c  input arguments:
c
c      x0 ,y0 ,z0  (dp,sc,in) - global coordinates of the integration point
c      vx0,vy0,vz0 (dp,sc,in) - components of the velocity of the integration
c                               point in global cartesian coordinates
c      r (dp,ar(3,3),in)      - components of the normal vector to the
c                               structure
c      dmax        (dp,sc,in) - the maximum distance for false contact
c                               (if abs(d).gt.dmax, false contact is assumed)
c  output arguments:
c
c       epsil    (dp,ar(7),inout)
c
c        epsil(1)   = tangent velocity of structure in first dir.
c        epsil(2)   = tangent velocity of structure in second dir.
c        epsil(3)
c                d  = distance between the integration point
c                     and the foundation
c                     ( d < 0 : contact    )
c                     ( d > 0 : no contact )
c        epsil(4-6) = rotations
c        epsil(7)   = angle between the beams
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
c              = 2  : if there is no contact because of normal consistency
c                                      or because abs(d) > dmax
c       delxy    (dp,ar(2),out) - distance between integartion point to slave node
c         delxy(1)          - x component
c         delxy(2)          - y component
c         delxy(3)          - z component
c
#include "impcom.inc"
c
      integer nocon,isegda(*),locu,locv,debug
      double precision v_dot
c
      double precision rad(2),x0,y0,z0,xa(7),ya(7),za(7),vx0,vy0,vz0,
     &  dmax,epsil(7),fifoa(7),fifob(7),vax(7),vay(7),vaz(7),
     &  delxy(3),gapshft,dir(3),ro(3,3),norm,
     &  r(3,3),xb(7),yb(7),zb(7),vbx(7),vby(7),vbz(7),
     &  propen,tols,dmax0,rx(7),ry(7),rz(7)

      double precision ZERO,ONE,TWO,d,vpx,vpy,vpz,d1,d2,
     &   ra,rb,dd,xac,yac,zac,xbc,ybc,zbc,r1,r2,r3,ac(3),bc(3),
     &   root1(65),root2(65),tmpSolutions(65),dfa(7),dfb(7),
     &   BestU,BestV,BestD,con1,alpha,adis,a(3),b(3),ab(3),vect(6),
     &   amax,amin,bmax,bmin,pBestU,pBestV,pBestD,dist,disto,
     &   angtols
      integer numRoot,i,j,tmpNumSolutions,found,iroot,ipass
c
      parameter (ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0,angtols=0.02d0)

#if defined(DEBUG)
      call SubBeg('BeamToBeam')
#endif
      BestU = HUGE
      BestV = HUGE
      if (ipass.eq.-1) then
          dmax = 5.d0*dmax0
          goto 10
      else
          dmax = dmax0
      endif
      nocon = 2
      d = dmax
      BestD = dmax
      x0 = (xa(1)+xa(2)+xa(3))/3.d0
      y0 = (ya(1)+ya(2)+ya(3))/3.d0
      z0 = (za(1)+za(2)+za(3))/3.d0
      amin = min(xa(1),xa(2),xa(3))
      bmax = max(xb(1),xb(2),xb(3))
      if (amin.gt.bmax+dmax) goto 999
      amax = max(xa(1),xa(2),xa(3))
      bmin = min(xb(1),xb(2),xb(3))
      if (bmin.gt.amax+dmax) goto 999

      amin = min(ya(1),ya(2),ya(3))
      bmax = max(yb(1),yb(2),yb(3))
      if (amin.gt.bmax+dmax) goto 999
      amax = max(ya(1),ya(2),ya(3))
      bmin = min(yb(1),yb(2),yb(3))
      if (bmin.gt.amax+dmax) goto 999

      amin = min(za(1),za(2),za(3))
      bmax = max(zb(1),zb(2),zb(3))
      if (amin.gt.bmax+dmax) goto 999
      amax = max(za(1),za(2),za(3))
      bmin = min(zb(1),zb(2),zb(3))
      if (bmin.gt.amax+dmax) goto 999

c --- check of pinball region
c
      xac = x0
      yac = y0
      zac = z0
      r1 = (xac-xa(1))*(xac-xa(1))+(yac-ya(1))*(yac-ya(1))+
     &     (zac-za(1))*(zac-za(1))
      r2 = (xac-xa(2))*(xac-xa(2))+(yac-ya(2))*(yac-ya(2))+
     &     (zac-za(2))*(zac-za(2))
      r3 = (xac-xa(3))*(xac-xa(3))+(yac-ya(3))*(yac-ya(3))+
     &     (zac-za(3))*(zac-za(3))
      ra = max (r1,r2,r3)
      xbc = (xb(1)+xb(2)+xb(3))/3.d0
      ybc = (yb(1)+yb(2)+yb(3))/3.d0
      zbc = (zb(1)+zb(2)+zb(3))/3.d0
      r1 = (xbc-xb(1))*(xbc-xb(1))+(ybc-yb(1))*(ybc-yb(1))+
     &     (zbc-zb(1))*(zbc-zb(1))
      r2 = (xbc-xb(2))*(xbc-xb(2))+(ybc-yb(2))*(ybc-yb(2))+
     &     (zbc-zb(2))*(zbc-zb(2))
      r3 = (xbc-xb(3))*(xbc-xb(3))+(ybc-yb(3))*(ybc-yb(3))+
     &     (zbc-zb(3))*(zbc-zb(3))
      rb = max (r1,r2,r3)
      ra = sqrt(ra)+sqrt(rb)+dmax
      rb = (xac-xbc)*(xac-xbc)+(yac-ybc)*(yac-ybc)+(zac-zbc)*(zac-zbc)
      rb = sqrt(rb)

      if (rb.ge.ra) goto 999


c      if (debug.gt.0) then
c          write (6,*) 'Enter BeamToBeam'
c          write (6,*) 'INI pBestU,pBestV,pBestD', pBestU,pBestV,pBestD
c      endif    


 10   continue
      d = dmax
      nocon = 2
      alpha = 0.2d0
      root1(1) = pBestU
      root2(1) = pBestV
      call BeamToBeam_1 (xa(1),ya(1),za(1),xb(1),yb(1),zb(1),alpha,
     &                   locu,locv,root1(1),root2(1),numRoot,ipass)

      
c      if (debug.gt.0) then
c          write (6,*) 'numRoot',numRoot
c          write (6,*) 'root1(1:numRoot)',root1(1:numRoot)
c          write (6,*) 'root2(1:numRoot)',root2(1:numRoot)          
c          write (6,*) 'xa(1:3)=[', xa(1),xa(3),xa(2),'];'
c          write (6,*) 'ya(1:3)=[', ya(1),ya(3),ya(2),'];'
c          write (6,*) 'xb(1:3)=[', xb(1),xb(3),xb(2),'];'
c          write (6,*) 'yb(1:3)=[', yb(1),yb(3),yb(2),'];'
c      endif
      
      if (numRoot.eq.0 .and. pBestU.ne.HUGE) then
          numRoot = 1
          root1(1) = pBestU
          root2(1) = pBestV
c          if (debug.gt.0) then
c              write (6,*) 'numRoot CORRECTED',numRoot,locu,locv      
c              write (6,*) 'root1(1:numRoot)',root1(1:numRoot)
c              write (6,*) 'root2(1:numRoot)',root2(1:numRoot)          
c          endif
      endif
      
      if (pBestU.eq.HUGE) then
          if (numRoot.gt.2) then
c --- !!! do not handle number of contacts higher than 2
              goto 999
c              write (6,*) 'WTF 1'
          elseif (numRoot.eq.2) then
c --- !!! do not handle number of contacts higher than 1 per segment length
              if ((abs(root1(1)).le.ONE.and.abs(root2(1)).le.ONE) .and.
     &            (abs(root1(2)).le.ONE.and.abs(root2(2)).le.ONE)) then
                   goto 999
c              write (6,*) 'WTF 2'
              endif
          endif    
      endif

      tmpNumSolutions = 0
      do 550 j=1,numRoot
         call BeamToBeam_2 (xa(1),ya(1),za(1),xb(1),yb(1),zb(1),alpha,
     &                      locu,locv,root1(j),root2(j),
     &                      fifoa(1),fifob(1),dd,
     &                      ac(1),bc(1),dfa(1),dfb(1),0)
         dd = dd - abs(rad(1)) - abs(rad(2))
         if (abs(dd).gt.dmax) goto 550
         tmpSolutions(tmpNumSolutions*3+1) = root1(j)
         tmpSolutions(tmpNumSolutions*3+2) = root2(j)
         tmpSolutions(tmpNumSolutions*3+3) = dd
         tmpNumSolutions = tmpNumSolutions + 1
 550  continue


c      if (debug.gt.0) then
c          write (6,*) 'tmpNumSolutions',tmpNumSolutions,dd,dmax,numRoot
c      endif

      if (tmpNumSolutions.le.0) goto 999

c --- find the closest solution
      disto = HUGE
      iroot = 0
      con1 = ONE + 1.d-8

      do 700 i=0,tmpNumSolutions-1
         vect(1) = tmpSolutions(i*3+1)
         vect(2) = tmpSolutions(i*3+2)
         dd = tmpSolutions(i*3+3)
         found = 0
         if (pBestU.eq.HUGE) then
             if (abs(vect(1)).gt.con1.or.abs(vect(2)).gt.con1) goto 700
             if (dd.le.ZERO.and.BestD.gt.ZERO) found = 1
             if ((dd.le.ZERO.and.BestD.le.ZERO).or.
     &           (dd.ge.ZERO.and.BestD.ge.ZERO)) then
                  if (abs(dd).lt.BestD) found = 1
             endif             
         else
c --- logic for existing solution
             dist = (pBestU-vect(1))**2 +
     &              (pBestV-vect(2))**2 +
     &              ((pbestD-dd)/dmax)**2
             if (dist.lt.disto) then 
                 found = 1
                 disto = dist
             endif
         endif
         if (found.eq.1) then
             BestU = vect(1)
             BestV = vect(2)
             BestD = dd
             iroot = 1
         endif
 700  continue

      if (iroot.eq.1 .and. 
     &   (abs(BestU).gt.con1 .or. abs(BestV).gt.con1)) then
c --- if we go over the line, reduce number of found contacts
          iroot = 0
      endif

c      if (debug.gt.0) then
c          write (6,*) 'iroot, found', iroot, found
c      endif

      if (iroot.eq.0) goto 800
      if (abs(BestD).ge.dmax) goto 800

 702  continue
c --- because of logic the position of nocon is different from that in the fast version 
      nocon = 0
      
      d = BestD

c --- interpolation function of contact
c
      call BeamToBeam_2 (xa(1),ya(1),za(1),xb(1),yb(1),zb(1),alpha,
     &                   locu,locv,BestU,BestV,fifoa(1),fifob(1),dd,
     &                   ac(1),bc(1),dfa(1),dfb(1),1)

      dir(1) = ro(3,1)
      dir(2) = ro(3,2)
      dir(3) = ro(3,3) 
      
c --- calculate contact normal and tangent trasformatin
      if (dd.gt.1.0d-12*dmax) then 

          do i =1,3
             r(3,i) = (bc(i)-ac(i))/dd
          enddo

c      if (debug.gt.0) then
c          write (6,*) 'large distanse', dd, r(3,1:3)
c      endif

c --- check on the normal change
          norm = dir(1)*r(3,1)+dir(2)*r(3,2)+dir(3)*r(3,3)
          if (norm.lt.ZERO) then

c      if (debug.gt.0) then
c          write (6,*) 'normal reverse', norm, dir(1:3)
c      endif

c --- reverse normal
              r(3,1) = -r(3,1)
              r(3,2) = -r(3,2)
              r(3,3) = -r(3,3)
              d = -dd - abs(rad(1)) - abs(rad(2))
          endif
      elseif ((dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3)).gt.
     &         sqrt(TINY)) then
c      if (debug.gt.0) then
c          write (6,*) 'normal retained', dir(1:3)
c      endif

c --- the normal will be retained
c --- check on the normal change
          do i =1,3
             r(3,i) = (bc(i)-ac(i))
          enddo
          norm = dir(1)*r(3,1)+dir(2)*r(3,2)+dir(3)*r(3,3)
          if (norm.lt.ZERO) then
              d = -dd - abs(rad(1)) - abs(rad(2))
          endif
          r(3,1) = dir(1)
          r(3,2) = dir(2)
          r(3,3) = dir(3)
      else

c --- this fibers cross !!!
          a(1) = dfa(1)*xa(1)
          a(2) = dfa(1)*ya(1)
          a(3) = dfa(1)*za(1)
          b(1) = dfb(1)*xb(1)
          b(2) = dfb(1)*yb(1)
          b(3) = dfb(1)*zb(1)
          do i = 2,7
             a(1) = a(1) + dfa(i)*xa(i)
             a(2) = a(2) + dfa(i)*ya(i)
             a(3) = a(3) + dfa(i)*za(i)
             b(1) = b(1) + dfb(i)*xb(i)
             b(2) = b(2) + dfb(i)*yb(i)
             b(3) = b(3) + dfb(i)*zb(i)
          enddo
          adis = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
          a(1) = a(1)/adis
          a(2) = a(2)/adis
          a(3) = a(3)/adis
          adis = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
          b(1) = b(1)/adis
          b(2) = b(2)/adis
          b(3) = b(3)/adis
          
          ab(1) = a(1) + b(1)
          ab(2) = a(2) + b(2)
          ab(3) = a(3) + b(3)

          adis = sqrt(ab(1)*ab(1)+ab(2)*ab(2)+ab(3)*ab(3))

c      if (debug.gt.0) then
c          write (6,*) 'fibers cross',adis,pBestU
c      endif

          if (pBestU.eq.HUGE) then
c --- these fibers share locations!
              if (adis.lt.angtols .or. (TWO-adis).lt.angtols) then
                  d = dmax
                  nocon = 2
c                  write (6,*) 'WTF 3'
                  goto 999
              endif
          else
c              if (adis.lt.sqrt(TINY) .or. (TWO-adis).lt.sqrt(TINY)) then
c --- these fibers share locations!
c                  d = dmax
c                  nocon = 2
c                  goto 999
c              endif
          endif
          
          vect(1) = xa(1)
          vect(2) = ya(1)
          vect(3) = za(1)
          vect(4) = xa(1)+ab(1)
          vect(5) = ya(1)+ab(2)
          vect(6) = za(1)+ab(3)

          call LineGeomTrans (vect(1),3,r(1,1))          
                    
          goto 200
      endif
      
      a(1) = dfa(1)*xa(1)
      a(2) = dfa(1)*ya(1)
      a(3) = dfa(1)*za(1)
      b(1) = dfb(1)*xb(1)
      b(2) = dfb(1)*yb(1)
      b(3) = dfb(1)*zb(1)
      do i = 2,7
         a(1) = a(1) + dfa(i)*xa(i)
         a(2) = a(2) + dfa(i)*ya(i)
         a(3) = a(3) + dfa(i)*za(i)
         b(1) = b(1) + dfb(i)*xb(i)
         b(2) = b(2) + dfb(i)*yb(i)
         b(3) = b(3) + dfb(i)*zb(i)
      enddo
      adis = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      if (adis.lt.sqrt(TINY)) then
          d = dmax
          nocon = 2
c          write (6,*) 'WTF 4'
          goto 999
      endif

      a(1) = a(1)/adis
      a(2) = a(2)/adis
      a(3) = a(3)/adis
      adis = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
      if (adis.lt.sqrt(TINY)) then
          d = dmax
          nocon = 2
c          write (6,*) 'WTF 5'
          goto 999
      endif

      b(1) = b(1)/adis
      b(2) = b(2)/adis
      b(3) = b(3)/adis
      
      epsil(7) = v_dot(a(1),b(1),3)
      if (epsil(7).lt.-ONE) epsil(7) = -ONE
      if (epsil(7).gt.ONE) epsil(7) = ONE
      
      epsil(7) = acos(epsil(7))
      
      ab(1) = a(1) + b(1)
      ab(2) = a(2) + b(2)
      ab(3) = a(3) + b(3)
      adis = sqrt(ab(1)*ab(1)+ab(2)*ab(2)+ab(3)*ab(3))
      if (pBestU.eq.HUGE) then
c --- these fibers are nearly parallel !!!
          if (abs(sin(epsil(7))) .lt. angtols) then
              d = dmax
c              write (6,*) 'WTF 2'
              nocon = 2
              goto 999
          endif
      else
          if (adis.lt.sqrt(TINY) .or. (TWO-adis).lt.sqrt(TINY)) then
c --- these fibers are parallel !!!
              vect(1) = xa(1)
              vect(2) = ya(1)
              vect(3) = za(1)
              vect(4) = xa(1)+ab(1)
              vect(5) = ya(1)+ab(2)
              vect(6) = za(1)+ab(3)

              call LineGeomTrans (vect(1),3,r(1,1))          
                       
              goto 200
c              d = dmax
c              nocon = 2
c              goto 999
          endif
      endif
c      if (debug.gt.0) then
c          write (6,*) 'adis, dd', adis, dd
c      endif

      r(1,1) = ab(1)/adis
      r(1,2) = ab(2)/adis
      r(1,3) = ab(3)/adis
      r(2,1) = -r(3,2)*r(1,3) + r(3,3)*r(1,2)
      r(2,2) = -r(3,3)*r(1,1) + r(3,1)*r(1,3)
      r(2,3) = -r(3,1)*r(1,2) + r(3,2)*r(1,1)
 200  continue

c --- finding intersection coords
c
      x0 = ac(1)
      y0 = ac(2)
      z0 = ac(3)
      vx0 = ZERO
      vy0 = ZERO
      vz0 = ZERO

      do i=1,7
         vx0 = vx0 + vax(i)*fifoa(i)
         vy0 = vy0 + vay(i)*fifoa(i)
         vz0 = vz0 + vaz(i)*fifoa(i)
      enddo

c --- relative velocity
c
      vpx = ZERO
      vpy = ZERO
      vpz = ZERO
      do 720 i = 1,7
         vpx = vpx + fifob(i)*vbx(i)
         vpy = vpy + fifob(i)*vby(i)
         vpz = vpz + fifob(i)*vbz(i)
  720 continue
      vpx = vx0 - vpx
      vpy = vy0 - vpy
      vpz = vz0 - vpz
c
      epsil(1)= vpx*r(1,1) + vpy*r(1,2) + vpz*r(1,3)
      epsil(2)= vpx*r(2,1) + vpy*r(2,2) + vpz*r(2,3)

c --- moisture and rotations on fondation
c
      epsil(4) = ZERO
      epsil(5) = ZERO
      epsil(6) = ZERO
      do 730 i = 1,7
         epsil(4) = epsil(4) + rx(i)*fifob(i)
         epsil(5) = epsil(5) + ry(i)*fifob(i)
         epsil(6) = epsil(6) + rz(i)*fifob(i)
  730 continue
      goto 900
  800 continue

c --- if the search in the patch fails but there are some intersections
c --- found, we take the solution from the extended path
      bestD = dmax
      disto = HUGE
      do 750 i=0,tmpNumSolutions-1
         vect(1) = tmpSolutions(i*3+1)
         vect(2) = tmpSolutions(i*3+2)
         dist = vect(1)**2+vect(2)**2
         if (dist.lt.disto) then
             BestU = vect(1)
             BestV = vect(2)
             disto = dist
             BestD = tmpSolutions(i*3+3)
         endif
 750  continue

      if (abs(bestD).lt.dmax) found = 1
      if (found.eq.0) goto 999
      if (BestU.eq.HUGE) goto 999

      d = BestD
      if (ipass.eq.-1) goto 702

c --- select contact status
      con1 = ONE + 1.d-8
      if (abs(BestD).lt.0.1d0*dmax.and.BestD.lt.sqrt(TINY))
     &    con1 = ONE + 1.d-5*tols
      if (ipass.ge.1) then
          con1 = ONE + tols*10.d0**(ipass-6)
      endif

     
      propen = ONE  + 0.01d0*tols

c      if (debug.gt.0) then
c          write (6,*) 'found, BestU, BestV',found,BestU,BestV
c          write (6,*) 'propen,con1',ipass,propen,con1
c      endif

      if (abs(BestU).gt.propen) then
          nocon = 2
          goto 999
      endif    

      if (BestV.gt.con1) then
          nocon = 1
c          if (isegda(2).ne.0.or.BestD.gt.ZERO) goto 900
          if (isegda(2).ne.0) goto 900
          if (abs(BestV).gt.propen) goto 900
      elseif (BestV.lt.-con1) then
          nocon = -1
c          if (isegda(1).ne.0.or.BestD.gt.ZERO) goto 900
          if (isegda(1).ne.0) goto 900
          if (abs(BestV).gt.propen) goto 900
      endif

c      if (debug.gt.0) then
c          write (6,*) '702'
c      endif

      goto 702

 900  continue

      delxy(1) = r(3,1)
      delxy(2) = r(3,2)
      delxy(3) = r(3,3)
      d = d - gapshft
      if (nocon.eq.0) goto 999
      if (abs(d).lt.0.1d0*dmax) goto 999
      if (abs(BestV).lt.2.d0.and.abs(BestU).lt.2.d0) goto 999
      if (iroot.eq.0) goto 999
      d1 = (xb(1)-x0)*(xb(1)-x0)+(yb(1)-y0)*(yb(1)-y0)+
     &     (zb(1)-z0)*(zb(1)-z0)
      d2 = (xb(3)-x0)*(xb(3)-x0)+(yb(3)-y0)*(yb(3)-y0)+
     &     (zb(3)-z0)*(zb(3)-z0)
      if (d1.ge.d2) then
          d = sqrt(d2)
          if (d.ge.dmax) then
              d = dmax
              goto 999
          endif
          if (d.gt.sqrt(TINY)*dmax) then
              delxy(1) = (x0-xb(3))/d
              delxy(2) = (y0-yb(3))/d
              delxy(3) = (z0-zb(3))/d
          endif
      else
          d = sqrt(d1)
          if (d.ge.dmax) then
              d = dmax
              goto 999
          endif
          if (d.gt.sqrt(TINY)*dmax) then
              delxy(1) = (x0-xb(1))/d
              delxy(2) = (y0-yb(1))/d
              delxy(3) = (z0-zb(1))/d
          endif
      endif

 999  continue
      epsil(3) = d
      pBestU = BestU
      pBestV = BestV
      pBestD = BestD

c --- the exact location of the contact point
      if (nocon.eq.0) then
          x0 = x0 + rad(1)*r(3,1)
          y0 = y0 + rad(1)*r(3,2)
          z0 = z0 + rad(1)*r(3,3)
      endif
      
c      if (debug.gt.0) then
c          write (6,*) 'BestU,BestV,BestD,nocon', BestU,BestV,BestD,nocon
c          write (6,*) 'isegda(1:2)',isegda(1:2)
c          write (6,*) 'Exit BeamToBeam'
c      endif

#if defined(DEBUG)
      call SubEnd('BeamToBeam')
#endif

      return
      end 
