      subroutine EstimIterNeeded(c,nfit,IterNum,IterMax,nIterEsti)
c --- primary function: estimate the number of iteration to attain convergence
c
c     Estimate is obtained by linear extrapolation of
c     log of normalized convergence values to zero using
c     least squares linear curve fit to last nfit points
c
c                            |  1   2  ....  n  | t
c                       a =  |                  |
c                            |  1   1  ....  1  |
c
c                           |  n*(n+1)*(2n+1)/6       n*(n+1)/2  |
c     [(a-transpose)*(a)] = |                                    |
c                           |     n*(n+1)/2               1      |
c
c             [(a-transpose)*(a)]{x} = [(a-transpose)]{c}
c
c             nIterEsti = - x2 / x1 + IterNum - n
c
c  input arguments:
c     c        (dp,ar(nfit),in) - array of logs of normalized
c                                 convergence values (ordering is cyclic)
c     n        (int,sc,in)      - number of values used in estimate
c     IterNum   (int,sc,in)     - current iteration number
c     IterMax   (int,sc,in)     - maximum iterations this substep
c
c  output arguments:
c     nIterEsti    (int,sc,out) - estimated iterations at convergence
c
#include "impcom.inc"
c
      integer i,j,joff,n,nfit,IterNum,IterMax,nIterEsti
      double precision c(nfit),a(2,2),y(2),con,con1
c

#if defined(DEBUG)
      call SubBeg('EstimIterNeeded')
#endif

      if (IterNum.lt.nfit) then
          n = IterNum
      else
          n = nfit
      endif

c --- check for enough points to draw a line
c
      if (n.le.1) then
          nIterEsti = IterMax*2
          goto 999
      endif

c --- calculate location of latest value
      joff = mod( IterNum-1,nfit ) + 1
      if (c(joff).le.0.0d0) then
          nIterEsti = IterNum
          goto 999
      endif

c --- calculate r.h.s.
c
      y(1) = 0.0d0
      y(2) = 0.0d0
      do 100 i = 1, n
         j = mod( (i+IterNum-1),n ) + 1
         y(1) = y(1) + dble(i)*c(j)
         y(2) = y(2) + c(j)
  100 continue

c --- calculate reciprocal of determinant of (a-transpose)*(a)
c
      con = -1.0d0 / dble( n*n*(n+1)*(n-1)/12 )

c --- form inverse of (a-transpose)*(a)
c
      a(1,1) = con*dble(n)
      a(1,2) = -con*dble( n*(n+1)/2 )
      a(2,1) = a(1,2)
      a(2,2) = con*dble( n*(n+1)*(2*n+1)/6 )

c --- estimate equilibrium iteration at convergence
c
      con  = a(1,1)*y(1) + a(1,2)*y(2)
      con1 = a(2,1)*y(1) + a(2,2)*y(2)
      if (abs(con) .gt. abs(con1)/200.d0) then
          nIterEsti = -nint(  con1 / con ) + IterNum - n
          if (nIterEsti .lt. (IterNum-n) .or.
     &        nIterEsti .gt. 2*IterMax ) then
              nIterEsti = IterMax*2
          elseif ( nIterEsti .le. IterNum ) then
              nIterEsti = IterNum + n
          endif
      else
          nIterEsti = IterMax*2
      endif
c
  999 continue
c

#if defined(DEBUG)
      call SubEnd('EstimIterNeeded')
#endif

      return
      end
