      subroutine BeamToBeamXY (rad,xa,ya,za,xb,yb,zb,dmax,epsil,
     &                         debug,tols)
c --- primary function: compute the relations between beam-beam crossing in the plane
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
      double precision rad(2),x0,y0,
     &  xa(3),ya(3),za(3),xb(3),yb(3),zb(3),
     &  dmax,epsil,fifoa(3),fifob(3)

      double precision ONE,d,d1,d2,propen,tols,
     &   ra,rb,dd,xac,yac,xbc,ybc,r1,r2,r3,ac(3),bc(3),
     &   root1(65),root2(65),tmpSolutions(65),dfa(3),dfb(3),
     &   BestD,a(3),b(3),zv0(3),amax,amin,bmax,bmin,ltol
      integer numRoot,i,j,tmpNumSolutions,debug
c
      parameter (ONE=1.0d0)
     
      numRoot = 0
      nocon = 2
      BestD = dmax
      d = dmax
      
c --- check only in-plane for the stats
      ltol = 
     &      sqrt((xa(1)-xa(2))**2+(ya(1)-ya(2))**2+(za(1)-za(2))**2)
      ltol = 0.01d0*tols*max(ltol,
     &      sqrt((xb(1)-xb(2))**2+(yb(1)-yb(2))**2+(zb(1)-zb(2))**2))
      ltol = 0.0d0
          
      amin = min(xa(1),xa(2),xa(3))
      bmax = max(xb(1),xb(2),xb(3))
      if (amin.gt.bmax+ltol) goto 999
      amax = max(xa(1),xa(2),xa(3))
      bmin = min(xb(1),xb(2),xb(3))
      if (bmin.gt.amax+ltol) goto 999
          
      amin = min(ya(1),ya(2),ya(3))
      bmax = max(yb(1),yb(2),yb(3))
      if (amin.gt.bmax+ltol) goto 999

      amax = max(ya(1),ya(2),ya(3))
      bmin = min(yb(1),yb(2),yb(3))
      if (bmin.gt.amax+ltol) goto 999

     
      x0 = (xa(1)+xa(2)+xa(3))/3.d0
      y0 = (ya(1)+ya(2)+ya(3))/3.d0

c --- check of pinball region
c
      xac = x0
      yac = y0
      r1 = (xac-xa(1))*(xac-xa(1))+(yac-ya(1))*(yac-ya(1))
      r2 = (xac-xa(2))*(xac-xa(2))+(yac-ya(2))*(yac-ya(2))
      r3 = (xac-xa(3))*(xac-xa(3))+(yac-ya(3))*(yac-ya(3))
      ra = max (r1,r2,r3)
      xbc = (xb(1)+xb(2)+xb(3))/3.d0
      ybc = (yb(1)+yb(2)+yb(3))/3.d0
      r1 = (xbc-xb(1))*(xbc-xb(1))+(ybc-yb(1))*(ybc-yb(1))
      r2 = (xbc-xb(2))*(xbc-xb(2))+(ybc-yb(2))*(ybc-yb(2))
      r3 = (xbc-xb(3))*(xbc-xb(3))+(ybc-yb(3))*(ybc-yb(3))
      rb = max (r1,r2,r3)
      ra = sqrt(ra)+sqrt(rb)+ltol
      rb = (xac-xbc)*(xac-xbc)+(yac-ybc)*(yac-ybc)
      rb = sqrt(rb)

      if (rb.ge.ra) goto 999

      call v_zero (zv0(1),3)
      d = dmax
      nocon = 2
      root1(1) = HUGE
      root2(1) = HUGE
      call BeamToBeam_1Fast (xa(1),ya(1),zv0(1),xb(1),yb(1),zv0(1),
     &                       root1(1),root2(1),numRoot,0)
      
c      if (debug.eq.1) then
c          write (8,*) 'numRoot', numRoot
c          write (8,*) 'root1(1:numRoot)', root1(1:numRoot)
c          write (8,*) 'root2(1:numRoot)', root2(1:numRoot)          
c      endif
      
      propen = ONE + 0.01d0*tols

      if (numRoot.gt.2) then
c --- !!! do not handle number of contacts higher than 2
          d = -HUGE
          nocon = 2
          goto 999
      elseif (numRoot.eq.2) then
c --- !!! do not handle number of contacts higher than 1 per segment length
          if((abs(root1(1)).le.propen.and.abs(root2(1)).le.propen) .and.
     &       (abs(root1(2)).le.propen.and.abs(root2(2)).le.propen)) then
              d = -HUGE
              nocon = 2
              goto 999
          endif   
      endif

      tmpNumSolutions = 0

      do 550 j=1,numRoot
          if (abs(root1(j)).gt.propen.or.abs(root2(j)).gt.propen) 
     &        goto 550
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

c --- if the search in the patch fails but there are some intersections
c --- found, we take the solution from the extended path
      bestD = dmax
      do 750 i=0, tmpNumSolutions-1
         dd = tmpSolutions(i*3+3)
         if (abs(dd).lt.BestD) then
             BestD = dd
             nocon = 0
         endif
 750  continue
      d = BestD
      
 999  continue

      if (abs(d).lt.HUGE.and.nocon.eq.2) d = dmax
      if (nocon.ne.0) d = dmax
      
      epsil = d

      return
      end