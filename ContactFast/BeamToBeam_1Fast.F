      subroutine BeamToBeam_1Fast (xa,ya,za,xb,yb,zb,u,v,numInter,ipass)
c --- primary function: beam-beam contact calculations (see parent function)
c
#include "impcom.inc"
c
      double precision xa(3),ya(3),za(3),xb(3),yb(3),zb(3)
      double precision u(20),v(20)
      integer numInter,ipass
c --- the coefficients of the two equations:
      double precision am(4,3),bm(4,3)
      double precision ZERO,ONE,TWO,delta,s(4,4),t(4,4),
     &                 st(2,100),lb,hb,tol
      integer i,j,k,numU,numV,numIt,numTmpInter,inFlag
      equivalence (numU,numV)
      parameter (ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0)

      call v_zero (st(1,1),200)
      
      if (u(1).eq.HUGE) then
          numU = 4
          delta=TWO/(dble(numU)-ONE)
      else
          numU = 1
          delta = ZERO
      endif
      
      if (ipass.eq.-1) then
         lb = -900.d0
         hb = 900.d0
      else
         lb = -10.d0
         hb = 10.d0
      endif

      numInter = 0
      call GetEqCoefFast (xa(1),ya(1),za(1),xb(1),yb(1),zb(1),
     &                am(1,1),bm(1,1))

      if (delta.gt.0.0d0) then
          do 140 i=1,numU
             do 150 j=1,numV
                s(j,i) = -ONE + (i-1)*delta
                t(j,i) = -ONE + (j-1)*delta
 150         continue
 140      continue
      else
         s(1,1) = u(1)
         t(1,1) = v(1)
      endif

      numIt=3
      tol = 0.007d0
      do 240 i=1,numU
        do 250 j=1,numV
           call SearchRootsFast (am(1,1),bm(1,1),
     &                       s(j,i),t(j,i),numIt,tol)
 250    continue
 240  continue

      numTmpInter=0
      do 340 i=1,numU
         do 350 j=1,numV
            if ((s(j,i).ge.lb .and. s(j,i).le.hb) .and.
     &          (t(j,i).ge.lb .and. t(j,i).le.hb)) then
               if (numTmpInter .gt. 0)then
                  inFlag = 0
                  do 360 k=1,numTmpInter
                    if (abs(s(j,i)-st(1,k)).lt.0.001d0 .and.
     &                  abs(t(j,i)-st(2,k)).lt.0.001d0) then
                       inFlag = 1
                       goto 380
                    endif
 360              continue
 380              continue
                  if (inFlag .eq. 0)then
                     st(1,numTmpInter+1)=s(j,i)
                     st(2,numTmpInter+1)=t(j,i)
                     numTmpInter = numTmpInter+1
                  endif
               else
                  st(1,numTmpInter+1)=s(j,i)
                  st(2,numTmpInter+1)=t(j,i)
                  numTmpInter = numTmpInter+1
               endif
            endif
 350    continue
 340  continue
      
      tol = 0.0000000001d0
      numIt=20
      do 440 i=1,numTmpInter
          call SearchRootsFast (am(1,1),bm(1,1),
     &                      st(1,i),st(2,i),numIt,tol)
 440  continue

      numInter=0
      do 540 i=1,numTmpInter
        if ((st(1,i).ge.lb .and. st(1,i).le.hb) .and.
     &      (st(2,i).ge.lb .and. st(2,i).le.hb))then
           inFlag = 0
           do 460 k=1,numInter
             if (abs(u(k)-st(1,i)).lt.0.00000001d0 .and.
     &           abs(v(k)-st(2,i)).lt.0.00000001d0)then
                 inFlag = 1
                 goto 480
             endif
 460       continue
 480       continue
           if (inFlag .eq. 0)then
               u(numInter+1)=st(1,i)
               v(numInter+1)=st(2,i)
               numInter = numInter+1
           endif
        endif
 540  continue

      return
      end
