      subroutine BeamToBeamFast (rad,xa,ya,za,xb,yb,zb,dmax,epsil,
     &                           debug,cont,tols,angtols)
c --- primary function: compute the relations between beam-beam
c
c  input arguments:
c
c      rad         (dp,ar(2),in)    - radii of two beams
c      xa,ya,za    (dp,ar(3),in)    - global coordinates of the the first beam
c      xb,yb,zb    (dp,ar(3),in)    - global coordinates of the the second beam
c      dmax        (dp,sc,in)       - the maximum distance for false contact
c                                    (if abs(d).gt.dmax, false contact is assumed)
c
c  output arguments:
c
c      epsil       (dp,sc,out)      - distance (gap > 0, overlap < 0)
c
c
#include "impcom.inc"
c
      integer nocon
c
      double precision v_dot
      
      double precision rad(2),x0,y0,z0,
     &  xa(3),ya(3),za(3),xb(3),yb(3),zb(3),
     &  dmax,epsil,fifoa(3),fifob(3)

      double precision ZERO,ONE,TWO,d,d1,d2,propen,tols,angtols,
     &   ra,rb,dd,xac,yac,zac,xbc,ybc,zbc,r1,r2,r3,ac(3),bc(3),
     &   root1(65),root2(65),tmpSolutions(65),dfa(3),dfb(3),
     &   BestU,BestV,BestD,adis,a(3),b(3),vect(6),ang,
     &   amax,amin,bmax,bmin,con1
      integer numRoot,i,j,tmpNumSolutions,found,iroot,debug,cont
c
      parameter (ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0)
     
      BestU = HUGE
      BestV = HUGE
      numRoot = 0
      nocon = 2
      BestD = dmax
      d = dmax
      
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
c      if (debug.eq.1) then
c          write (8,*) 'rb, ra', numRoot
c      endif

      if (rb.ge.ra) goto 999

 10   continue
      d = dmax
      nocon = 2
      root1(1) = HUGE
      root2(1) = HUGE
      call BeamToBeam_1Fast (xa(1),ya(1),za(1),xb(1),yb(1),zb(1),
     &                       root1(1),root2(1),numRoot,0)
      
      propen = ONE + 0.01d0*tols
c      if (debug.eq.1) then
c          write (8,*) 'numRoot', numRoot
c          write (8,*) 'root1(1:numRoot)', root1(1:numRoot)
c          write (8,*) 'root2(1:numRoot)', root2(1:numRoot)          
c      endif

      if (numRoot.gt.2) then
c --- !!! do not handle number of contacts higher than 2
          d = -HUGE
          nocon = 2
c          write (6,*) 'wtf 1'
          goto 999
      elseif (numRoot.eq.2) then
c --- !!! do not handle number of contacts higher than 1 per segment length
          if((abs(root1(1)).le.propen.and.abs(root2(1)).le.propen) .and.
     &       (abs(root1(2)).le.propen.and.abs(root2(2)).le.propen)) then
              d = -HUGE
              nocon = 2
c              write (6,*) 'wtf 2'
              goto 999
          endif   
      endif

      tmpNumSolutions = 0

      do 550 j=1,numRoot
         call BeamToBeam_2Fast (xa(1),ya(1),za(1),xb(1),yb(1),zb(1),
     &                          root1(j),root2(j),fifoa(1),fifob(1),dd,
     &                          ac(1),bc(1),dfa(1),dfb(1),0)
         dd = dd - abs(rad(1)) - abs(rad(2))

         if (abs(dd).gt.dmax) goto 550
         tmpSolutions(tmpNumSolutions*3+1) = root1(j)
         tmpSolutions(tmpNumSolutions*3+2) = root2(j)
         tmpSolutions(tmpNumSolutions*3+3) = dd
         tmpNumSolutions = tmpNumSolutions + 1
 550  continue
c
      if (tmpNumSolutions.le.0) goto 999
c --- find the closest solution
      iroot = 0
      con1 = ONE + 1.d-8

      do 700 i=0,tmpNumSolutions-1
         vect(1) = tmpSolutions(i*3+1)
         vect(2) = tmpSolutions(i*3+2)
         dd = tmpSolutions(i*3+3)
         found = 0
c --- logic for first pass or no previous solution
         if (abs(vect(1)).gt.con1.or.abs(vect(2)).gt.con1) goto 700
         if (dd.le.ZERO.and.BestD.gt.ZERO) found = 1
         if ((dd.le.ZERO.and.BestD.le.ZERO).or.
     &       (dd.ge.ZERO.and.BestD.ge.ZERO)) then
              if (abs(dd).lt.BestD) found = 1
         endif             
         if (found.eq.1) then
             BestU = vect(1)
             BestV = vect(2)
             BestD = dd
             iroot = 1
         endif
 700  continue

      if (iroot.eq.0) goto 800
      if (abs(BestD).ge.dmax) goto 800
c --- because of logic the position of nocon is different from that in the full version 
      nocon = 0
 702  continue

      d = BestD

c --- interpolation function of contact
c
      if (cont.ne.1 .and. 
     &    nocon.eq.0 .and. -d.ge.min(rad(1),rad(2))) then
c --- too large interpenetration
          nocon = 2
c          write (6,*) 'wtf 3'          
          goto 999
      endif

      call BeamToBeam_2Fast (xa(1),ya(1),za(1),xb(1),yb(1),zb(1),
     &                       BestU,BestV,fifoa(1),fifob(1),dd,
     &                       ac(1),bc(1),dfa(1),dfb(1),1)


      a(1) = dfa(1)*xa(1)
      a(2) = dfa(1)*ya(1)
      a(3) = dfa(1)*za(1)
      b(1) = dfb(1)*xb(1)
      b(2) = dfb(1)*yb(1)
      b(3) = dfb(1)*zb(1)
      do i = 2,3
         a(1) = a(1) + dfa(i)*xa(i)
         a(2) = a(2) + dfa(i)*ya(i)
         a(3) = a(3) + dfa(i)*za(i)
         b(1) = b(1) + dfb(i)*xb(i)
         b(2) = b(2) + dfb(i)*yb(i)
         b(3) = b(3) + dfb(i)*zb(i)
      enddo
      adis = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      if (adis.lt.sqrt(TINY)) then
          d = -HUGE
          nocon = 2
c          write (6,*) 'wtf 4'          
          goto 999
      endif
      a(1) = a(1)/adis
      a(2) = a(2)/adis
      a(3) = a(3)/adis
      adis = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
      if (adis.lt.sqrt(TINY)) then
          d = -HUGE
          nocon = 2
c          write (6,*) 'wtf 5'          
          goto 999
      endif
      b(1) = b(1)/adis
      b(2) = b(2)/adis
      b(3) = b(3)/adis

      ang = v_dot(a(1),b(1),3)
      if (ang.lt.-ONE) ang = -ONE
      if (ang.gt.ONE) ang = ONE
      
      ang = acos(ang)
      
      a(1) = a(1) + b(1)
      a(2) = a(2) + b(2)
      a(3) = a(3) + b(3)
      adis = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

c      if (debug.eq.1) then
c          write (8,*) 'adis, nocon', adis, nocon
c      endif

c --- these fibers are nearly parallel !!!
      if (abs(sin(ang)) .lt. angtols) then
          d = -HUGE
          nocon = 2
c          write (6,*) 'wtf 6'          
          goto 999
      endif

c --- finding intersection coords
c
      x0 = ac(1)
      y0 = ac(2)
      z0 = ac(3)

      goto 900
  800 continue

c --- if the search in the patch fails but there are some intersections
c --- found, we take the solution from the extended path
      bestD = dmax
      d = BestD
c      nocon = 2
      do 750 i=0, tmpNumSolutions-1
         dd = tmpSolutions(i*3+3)
         if (abs(dd).lt.BestD) then
             found = 1
             BestU = tmpSolutions(i*3+1) 
             BestV = tmpSolutions(i*3+2) 
             BestD = dd
         endif
 750  continue

c      if (debug.eq.1) then
c          write (8,*) 'BestD, dmax, nocon', BestD, dmax, nocon
c      endif
      
c --- if the search in the patch fails but there are some intersections
c --- found, we search the extended patch.
      if (abs(BestD).ge.dmax) goto 999
      d = BestD

c --- select contact status for "slide pass" contacts (% of the length)
      if (abs(BestU).gt.propen) then
          nocon = 2
          if (abs(d).ge.dmax) goto 999
          goto 702
      endif    
      if (BestV.gt.propen) then
          nocon = 1
      elseif (BestV.lt.-propen) then
          nocon = -1
      else
          nocon = 0
      endif
c --- force to go through 702 at any case
      goto 702
c
 900  continue
c      if (debug.eq.1) then
c          write (8,*) 'nocon, BestV, BestU', nocon, BestV, BestU
c          write (8,*) 'iroot, d', iroot, d
c      endif
    
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
      else
          d = sqrt(d1)
          if (d.ge.dmax) then
              d = dmax
              goto 999
          endif
      endif

 999  continue

      if (abs(d).lt.HUGE.and.nocon.eq.2) d = dmax
      if (nocon.ne.0) d = dmax
      
      epsil = d

c      if (debug.eq.1) then
c          write (8,*) 'nocon, epsil', nocon, epsil
c      endif

      return
      end