      subroutine BeamToSurf_1 (nfacn,x,y,z,off,nx,ny,nz,inors,itry,
     &                         u,v,numInter,ipass)
c --- primary function: solve coef(1)*x^2 + coef(2)*x + coef(3) = 0
c
#include "impcom.inc"
c
      double precision x(8),y(8),z(8),off(3),nx(8),ny(8),nz(8)
      double precision u(20),v(20)
      integer numInter,ipass
      double precision a(16,2)
      double precision delta,s(4,4),t(4,4),st(2,100),lb,hb,tol
      integer nfacn,i,j,k,numU,numV,numIt,numTmpInter,inFlag,inors,itry
      equivalence (numU,numV)

      call v_zero (st(1,1),200)

c --- checking one segment should be OK here
      numU = 1

      if (numU.gt.1) then
          delta=2.0d0/(dble(numU)-1.0d0)
          do 140 i=1,numU
             do 150 j=1,numV
                s(j,i) = -1.0d0 + (i-1)*delta
                t(j,i) = -1.0d0 + (j-1)*delta
 150         continue
 140      continue

      else
          s(1,1) = 0.0d0
          t(1,1) = 0.0d0
      endif    

      if (ipass.eq.-1) then
          lb = -900
          hb = 900
      else
          lb = -10
          hb = 10
      endif

      numInter = 0
      call GetEqCoef4 (x(1),y(1),z(1),nx(1),ny(1),nz(1),off(1),
     &                 inors,itry,a(1,1))


      numIt=3
      tol = 0.007d0
      do 240 i=1,numU
        do 250 j=1,numV
           call SearchRoots4 (itry,a(1,1),s(j,i),t(j,i),numIt,tol)
 250    continue
 240  continue

      numTmpInter=0
      do 340 i=1,numU
         do 350 j=1,numV
            if ((s(j,i).ge.lb .and. s(j,i).le.hb) .and.
     &         (t(j,i).ge.lb .and. t(j,i).le.hb))then
               if (numTmpInter .gt. 0)then
                   inFlag = 0
                   do 360 k=1,numTmpInter
                      if (abs(s(j,i)-st(1,k)).lt.0.001d0 .and.
     &                    abs(t(j,i)-st(2,k)).lt.0.001d0)then
                         inFlag = 1
                         goto 380
                      endif
 360               continue
 380               continue
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
         call SearchRoots4 (itry,a(1,1),st(1,i),st(2,i),numIt,tol)
 440  continue

      numInter=0
      do 540 i=1,numTmpInter
         if((st(1,i).ge.lb .and. st(1,i).le.hb) .and.
     &      (st(2,i).ge.lb .and. st(2,i).le.hb))then
            inFlag = 0
            do 460 k=1,numInter
               if(abs(u(k)-st(1,i)).lt.0.00000001d0 .and.
     &            abs(v(k)-st(2,i)).lt.0.00000001d0)then
                  inFlag = 1
                  goto 480
               endif
 460        continue
 480        continue
            if (inFlag .eq. 0)then
                u(numInter+1)=st(1,i)
                v(numInter+1)=st(2,i)
                numInter = numInter+1
            endif
         endif
 540  continue

      return
      end
