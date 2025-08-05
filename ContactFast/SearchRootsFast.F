      subroutine SearchRootsFast (a,b,s,t,numIt,error)
c --- primary function: solve system of quadratic equations
c                       coef(1)*x^2 + coef(2)*x + coef(3) = 0.
c
c  input arguments:
c     a,b   (dp,ar(4,3),in)    - the coeficiects
c     numIt (int,sc,in)      - number of interval for root search
c
c  output arguments:
c     s,t   (dp, sc, out)    - the roots.
c
#include "impcom.inc"
c
      double precision a(4,3),b(4,3)
      double precision s,t
      integer numIt
      double precision oldS,oldT,delta,f1,pf1_s,pf1_t,f2,pf2_s,pf2_t,
     &                 newS,newT,error,rs,rt
      integer i

      oldS = s
      oldT = t

      call SolveQuadEqFast (a(1,1),b(1,1),oldS,oldT,
     &                  f1,pf1_s,pf1_t,f2,pf2_s,pf2_t)
      if (abs(oldS).gt.1.d20 .or. abs(oldT).gt.1.d20) then
         s=-1000.0d0
         t=-1000.0d0
         goto 999
      endif
      do 100 i=1,numIt
            rs = pf2_t * f1 - pf1_t * f2
            rt = pf1_s * f2 - pf2_s * f1
            if (abs(rs).lt.TINY.and.abs(rt).lt.TINY) then
                s = oldS
                t = oldT
                goto 999
            endif
            delta = pf1_s * pf2_t - pf1_t * pf2_s
            if(abs(delta) .lt. 100.0d0*TINY)then
                  s=-1000.0d0
                  t=-1000.0d0
              goto 999
            endif
            delta = 1.0d0 / delta

            newS = oldS - rs * delta
            newT = oldT - rt * delta

            if(abs(newS-oldS).lt.error .and.
     &         abs(newT-oldT).lt.error)then
             s = newS
             t = newT
             goto 999
          endif

            oldS = newS
            oldT = newT

          if (abs(oldS).gt.1.d20 .or. abs(oldT).gt.1.d20) then
             s=-1000.0d0
             t=-1000.0d0
             goto 999
          endif
          call SolveQuadEqFast(a(1,1),b(1,1),oldS,oldT,
     &                     f1,pf1_s,pf1_t,f2,pf2_s,pf2_t)
 100  continue

      s=-1000.d0
      t=-1000.d0
c
 999  continue

      return
      end
