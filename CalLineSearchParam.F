      subroutine CalLineSearchParam (lsitr,lspar,lsval,enrm0,lsp,
     &                               curBracket,oldBracket,lspch)
c --- primary function: estimate the optimal line search parameter
c                       using Brent's algorithm; 9.3 of Numerical Recipes
c
c  input arguments:
c     lsitr  (dp,sc,in)     - line search iteration counter
c     enrm0  (dp,sc,in)     - initial energy norm (at lsp = 0)
c
c  input/output arguments:
c     lspar  (dp,ar(3),inout)  - array of 3 previous line search parameters
c     lsval  (dp,ar(3),inout)  - array of 3 previous energy norms
c
c  output arguments:
c     lsp    (dp,sc,out)    - estimated parameter which zeroes residual
c     lspch  (int,sc,out)   - flag indicates whether line search should
c                             continue
c  referenced arguments:
c     A_MinLineSrchVal  (dp,sc)        - smallest allowable line search parameter
c     A_LineSrchTol     (dp,sc)        - line search tolerance for convergence
c     A_MaxLineSrchIter (int,sc)       - maximum number of iterations
c
#include "impcom.inc"
#include "solu_info.inc"
c
      double precision lspar(3),lsval(3),enrm0,lsp,
     & curBracket,oldBracket,
     & dlspar,lsparo,bracket,PP,QQ,RR,SS,con
      integer lsitr,lspch

#if defined(DEBUG)
      call SubBeg('CalLineSearchParam')
#endif

      lsparo = lsp
      lspch = 1

c --- 1st iteration
      if (lsitr.eq.1) then
c --- check if line search should be performed
          if (lsval(2)*lsval(1).ge.0.0d0) then
c --- positive values - no line search needed
              lspch = 0
              goto 999
          endif
c --- interpolate between s=0 and s=1
          con = (lspar(1)*lsval(2) - lspar(2)*lsval(1))
          lsp = con/(lsval(2) - lsval(1))
          if (lsp.lt.A_MinLineSrchVal) lsp = A_MinLineSrchVal
          lspar(3) = lspar(2)
          lsval(3) = lsval(2)
          lspar(2) = lsp
          curBracket = lsp
          oldBracket = 1.0d0
c --- check if energy change is small - no line search needed
          if (abs(lsval(2)).lt.A_LineSrchTol*abs(enrm0)) then
c --- tolerance satisfied
              lspch = 0
              goto 999
          endif
          goto 999
      endif

c --- check if we are still bracketed; if not force bisection
      if (lsval(2)*lsval(3).gt.0.0d0) then
          lspar(3) = lspar(1)
          lsval(3) = lsval(1)
          curBracket = lspar(2) - lspar(1)
          oldBracket = curBracket
      endif

c --- check if lsval(2) is still our best estimate; if not swap
      if (abs(lsval(3)).lt.abs(lsval(2))) then
          con = lspar(2)
          lspar(2) = lspar(3)
          lspar(3) = con
          con = lsval(2)
          lsval(2) = lsval(3)
          lsval(3) = con
          lsp = lspar(2)
      endif

c --- define current half bracket
      bracket = 0.5d0*(lspar(3) - lspar(2))

c --- check for convergence
      if (abs(bracket).lt.A_MinLineSrchVal) then
c --- small bracket
          if (lsp.lt.A_MinLineSrchVal) lsp = A_MinLineSrchVal
          if (lsp.gt.1.0d0) lsp = 1.0d0
          lspch = 0
          goto 999
      endif

c --- check if we do quadratic interpolation or simple bisection
      if (abs(oldBracket).ge.A_MinLineSrchVal .and.
     &    abs(lsval(1)).gt.abs(lsval(2))) then
          SS = lsval(2)/lsval(1)
          if (lspar(1).eq.lspar(3)) then
c --- have only 2 points: use regula falsi
              PP = 2.0d0*bracket*SS
              QQ = 1.0d0 - SS
          else
c --- have 3 points, try quadratic interpolation
              QQ = lsval(1)/lsval(3)
              RR = lsval(2)/lsval(3)
              PP = SS*( 2.0d0*bracket*QQ*(QQ-RR) -
     &            (lspar(2)-lspar(1))*(RR-1.0d0) )
              QQ = (QQ-1.0d0)*(RR-1.0d0)*(SS-1.0d0)
          endif
          if (PP.gt.0.0d0) then
              QQ = -QQ
          else
              PP = -PP
          endif
          con = 3.0d0*bracket*QQ - 1.0d-6*QQ
          con = min(con,abs(oldBracket*QQ))
          if (2.0d0*PP.lt.con) then
c --- quadratic interpolation can be safely used
              oldBracket = curBracket
              curBracket = PP/QQ
          else
c --- can't use quadratic interp so use bisection
              curBracket = bracket
              oldBracket = curBracket
          endif
      else
c --- can't use quadratic interp
          if (abs(lsval(3)).lt.abs(10.0d0*enrm0)) then
c --- close enough, so use regula-falsi
              SS = lsval(2)/lsval(3)
              PP = 2.0d0*bracket*SS
              QQ = 1.0d0 - SS
              oldBracket = curBracket
              curBracket = -PP/QQ
          else
c --- still far away; use bisection
              curBracket = bracket
              oldBracket = curBracket
          endif
      endif

c --- move last best guess to the 1st slot
      lspar(1) = lspar(2)
      lsval(1) = lsval(2)

c --- update the new search parameter
      lspar(2) = lspar(2) + curBracket
      lsp = lspar(2)

c --- check limits   A_MinLineSrchVal <= lsp <= 1.0d0
      if (lsp.lt.A_MinLineSrchVal) lsp = A_MinLineSrchVal
      if (lsp.gt.1.0d0) lsp = 1.0d0
c
c --- check if value changed. done if not.
      dlspar = abs (lsp - lsparo)
      if (dlspar.lt.0.5d0*A_LineSrchTol*(lsp+lsparo)) then
c --- lsp the same
          lspch = 0
      endif

c --- check if tolerance satisfied
      if (abs(lsval(2)).lt.A_LineSrchTol*abs(enrm0)) then
c --- tolerance satisfied
          lspch = 0
      endif

c --- check if maximum number of iterations exceeded
      if (lsitr.ge.A_MaxLineSrchIter) then
c --- max iter exceeded
          lspch = 0
      endif
c
  999 continue

#if defined(DEBUG)
      call SubEnd('CalLineSearchParam')
#endif

      return
      end
