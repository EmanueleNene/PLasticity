      subroutine UpdRot(rVec,rCor)
c --- primary function: given a new correction vector, update the incremental rotations
c
c  input arguments:
c     rCor (dp, ar(3), in )      - correction rotation vector
c
c  output arguments:
c     rVec (dp, ar(3), inout)    - compounded rotation vector
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision rVec(3),rCor(3),wInc(3),wCor(3),
     & w0, dw, aSmallNumber, ZERO, ONE, HALF, TWO,
     & cw0, sw0, cdw, sdw, f0, df, Q(3),
     & DQ(3), Q0(3), w1, sw1, cw1, f1
c
      parameter (ZERO=0.0d0, ONE=1.0d0, TWO=2.0d0, HALF=0.5d0)
      parameter (aSmallNumber=1.0d-24)

      wInc(1) = rVec(1)
      wInc(2) = rVec(2)
      wInc(3) = rVec(3)
      wCor(1) = rCor(1)
      wCor(2) = rCor(2)
      wCor(3) = rCor(3)
      rVec(1) = rVec(1)+rCor(1)
      rVec(2) = rVec(2)+rCor(2)
      rVec(3) = rVec(3)+rCor(3)
      w0 = wInc(1)*wInc(1)+wInc(2)*wInc(2)+wInc(3)*wInc(3)
      dw = wCor(1)*wCor(1)+wCor(2)*wCor(2)+wCor(3)*wCor(3)
      if (w0.gt.aSmallNumber.and.dw.gt.aSmallNumber) then
          w0 = HALF*sqrt(w0)
          dw = HALF*sqrt(dw)
          cw0 = cos(w0)
          cdw = cos(dw)
          sw0 = sin(w0)
          sdw = sin(dw)
          f0 = HALF*sw0/w0
          df = HALF*sdw/dw
          Q0(1) = wInc(1)*f0
          Q0(2) = wInc(2)*f0
          Q0(3) = wInc(3)*f0
          DQ(1) = wCor(1)*df
          DQ(2) = wCor(2)*df
          DQ(3) = wCor(3)*df
          Q(1) = cdw*Q0(1)+cw0*DQ(1)-Q0(2)*DQ(3)+Q0(3)*DQ(2)
          Q(2) = cdw*Q0(2)+cw0*DQ(2)-Q0(3)*DQ(1)+Q0(1)*DQ(3)
          Q(3) = cdw*Q0(3)+cw0*DQ(3)-Q0(1)*DQ(2)+Q0(2)*DQ(1)
          sw1 = Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3)
          if (sw1.ne.ZERO) then
              sw1 = sqrt(sw1)
              if (sw1.gt.aSmallNumber) then
                  cw1 = cw0*cdw-(Q0(1)*DQ(1)+Q0(2)*DQ(2)+Q0(3)*DQ(3))
                  w1 = atan2(sw1,cw1)
              else
                  w1 = sw1
              endif
             dw = w1*(ONE-(Q(1)*Q0(1)+Q(2)*Q0(2)+Q(3)*Q0(3))/(sw0*sw1))
             if (dw.gt.HALF*W_Pi) w1 = w1-W_Pi
             f1 = TWO*w1/sw1
             rVec(1) = f1*Q(1)
             rVec(2) = f1*Q(2)
             rVec(3) = f1*Q(3)
          else
             rVec(1) = ZERO
             rVec(2) = ZERO
             rVec(3) = ZERO
          endif
      endif

      return
      end
