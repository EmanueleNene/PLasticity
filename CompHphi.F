      subroutine CompHphi(phi, Hmatx, InvFlag)
c --- primary function: compute H(phi) matrix or its inverse
c
c     NOTES:
c     
c     H(phi) = n*transpose(n) + (sin(theta)/theta)*(1-n*transpose(n)) +
c                              ((1-cos(theta))/theta^2)*skew(theta)
c
c --- Inverse of H(phi)
c     Hinv(phi) = n*transpose(n) + 
c                 [(thetam/2)/tan(thetam/2)]*(1 - n*transpose(n)) -
c                 skew(theta)
c --- n = theta/thetam, thetam = sqrt(theta*theta)
c --- theta - incremental rotation vector for substep
c
c  input arguments:
c     phi(dp,ar(3),in)       - incremental rotation vector
c     InvFlag(int,sc,in)     - flag to indicate:
c                               = 0 : return H(phi)
c                               = 1 : return inverse of H(phi)
c
c  output arguments:
c     Hmatx(dp,ar(3,3),out)  - H(phi) matrix or its inverse
c                              depending on the InvFlag
c
#include "impcom.inc"
c
      integer InvFlag
      double precision phi(*), Hmatx(3,*)
      integer i, j
      double precision phiMag, phimag2, phimag4, phimag6, phimag8
      double precision rVec(3), fac1, fac2, phiSkew(3,3)
      double precision SMALL1, SMALL2, C1, C2, C3, C4, D1, D2, D3
      double precision ZERO, ONE, TWO, HALF

      parameter(ZERO=0.D0, ONE=1.D0, TWO=2.D0, HALF=0.5D0)
      parameter (SMALL1=1.0d-2, SMALL2=1.0d-16)
      parameter (C1=0.33333333d0, C2=0.13333333d0, C3=0.053968254d0,
     x           C4=0.0218694885d0)
      parameter (D1=0.16666667D0, D2=0.00833333D0, D3=0.00019841D0)

#if defined(DEBUG)
      call SubBeg('CompHphi')
#endif

      call v_zero(Hmatx,9)
c --- compute magnitude of incremental rotation
      phiMag = sqrt(phi(1)*phi(1)+phi(2)*phi(2)+phi(3)*phi(3))
      if (phiMag.gt.SMALL2) then
c --- compute rotation axis for incremental rotation
          rVec(1) = phi(1)/phiMag
          rVec(2) = phi(2)/phiMag
          rVec(3) = phi(3)/phiMag
c --- compute skew symmetric matrix for incremental rotation
          call skewA3x3(phi(1), phiSkew(1,1))
      else
          Hmatx(1,1) = ONE
          Hmatx(2,2) = ONE
          Hmatx(3,3) = ONE
          goto 999
      endif

      if (InvFlag.eq.0) then
c --- compute H(phi) matrix
          if (phiMag.gt.SMALL1) then
              fac1 = sin(phiMag)/phiMag
              fac2 = (ONE-cos(phiMag))/(phiMag*phiMag)
          else
c --- evaluate using series expansion of sin(x)/x
c     fac1 = sin(x)/x = 1 - (1/6)x^2 + (1/120)x^4 - (1/5040)x^6+(1/362880)*x^8
c     fac2 = (1/x)(1-cos(x)) = 0.5* (1/y)sin(y) * (1/y)sin(y)... y = x/2
              phimag2 = phiMag*phiMag
              phimag4 = phimag2*phimag2
              phimag6 = phimag4*phimag2
              fac1 = ONE - D1*phimag2 + D2*phimag4 - D3*phimag6

              phiMag = phiMag/TWO
              phimag2 = phiMag*phiMag
              phimag4 = phimag2*phimag2
              phimag6 = phimag4*phimag2
              fac2 = ONE - D1*phimag2 + D2*phimag4 - D3*phimag6
              fac2 = HALF*fac2*fac2
         endif
         do i = 1, 3
            do j = 1, 3
               if (i.eq.j) then
                   Hmatx(i,j) = fac1 + rVec(i)*rVec(j)*(ONE-fac1)
               else
                   Hmatx(i,j) = rVec(i)*rVec(j)*(ONE-fac1) +
     $                          fac2*phiSkew(i,j)
               endif
            enddo
         enddo
      elseif (InvFlag.eq.1) then
c --- compute inverse of H(phi) matrix
         phiMag = phiMag/TWO
         fac2 = ZERO
         if (phiMag.gt.SMALL1) then
             fac2 = tan(phiMag)/phiMag
         else
c --- evaluate tan(thetam/2)/(thetam/2) using Taylor series:
c --- tan(x)/x = 1 + (1/3)x^2 + (2/15)*x^4 + (17/315)*x^6 + (62/2835)*x^8 ... for x < pi/2
             phimag2 = phiMag*phiMag
             phimag4 = phimag2*phimag2
             phimag6 = phimag4*phimag2
             phimag8 = phimag6*phimag2
             fac2 = ONE + C1*phimag2 + C2*phimag4 + C3*phimag6 +
     &                                              C4*phimag8
         endif
c --- assemble terms for inverse of H(phi)
         if (fac2.gt.ZERO) then
             fac2 = ONE/fac2
             do i = 1,3
                do j = 1,3
                   if (i.eq.j) then
                       Hmatx(i,j) = fac2 + 
     &                              rVec(i)*rVec(j)*(ONE-fac2)
                   else
                       Hmatx(i,j) = rVec(i)*rVec(j)*(ONE-fac2) -
     &                              HALF*phiSkew(i,j)
                   endif
                enddo
             enddo
         endif
      endif

 999  continue

#if defined(DEBUG)
      call SubEnd('CompHphi')
#endif

      return
      end
