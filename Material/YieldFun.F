      subroutine YieldFun(nprop,prop,proppl,qEl,pleq,dpleq,sigy,dsigdep)
c --- primary function: calculate material yield stress and tangent slope
c
c input arguments
c      nprop    (int,sc,in)              -  number of properties defined
c      prop     (dp, ar,in)              -  elastic material property table
c      proppl   (dp, ar,in)              -  plastic material property table
c      qEl      (dp, ar,in)              -  trial stress
c      pleq     (dp, sc,in)              -  equivalent plastic strain
c      dpleq    (dp, sc,in)              -  increment of equivalent plastic strain
c
c output arguments
c      sigy     (dp,sc,out)              -  yield stress at 't+dt'
c      dsigdep  (dp,sc,out)              -  plastic slope at 't+dt'
c
#include "impcom.inc"
c
      integer          nprop
      double precision prop(*),proppl(nprop),pleq,sigy,dsigdep,dpleq,qEl
      double precision Ex,sigy0,Et,q1,q2,b,con1,con2,sigy_0,ds0dep
      double precision pleq_t

      double precision ONE, SMALL, ZERO
      parameter       (ZERO       = 0.d0,
     &                 ONE        = 1.d0,
     &                 SMALL      = 1.d-10)
      
      Ex = prop(1)
      sigy0 = proppl(1)
      Et = proppl(2)
      q1 = Ex * Et / (Ex - Et)
      if (nprop.eq.4) then
c --- sigy = sigy0 + q1 *ep + q2 * (1-exp(-b * ep))
          q2 = proppl(3)
          b = proppl(4)
          con1 = exp (-b * pleq)
          sigy = sigy0 + q1 * pleq + q2 * ( ONE - con1)
          dsigdep = q1 + q2 * b * con1
          if (qEl .gt. sigy) then
              pleq_t = pleq - dpleq
              if (q1.lt.ZERO .and. pleq_t.gt.ZERO) then
c --- account for softening and failure
                  con2 = sigy0 + q1 * pleq_t 
     &                 + q2 * (ONE-exp(-b * pleq_t))
                  if (sigy.lt.0.03d0*sigy0) then
                      sigy = 0.005d0*sigy0
                      dsigdep = SMALL*Ex
                  endif
              endif
          endif
      else
c --- sigy = sigy0 + q1 *ep
          sigy = sigy0 + pleq * q1
          dsigdep = q1
      endif

      return
      end