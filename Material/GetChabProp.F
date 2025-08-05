      subroutine GetChabProp (keycpL, moist, numpt, kkh, ckh, gkh, dckh)
c *** primary function: obtain material parameters for Chaboche model
c
c  input arguments:
c     keycpL   (int*8,sc,in)       - key indicating the material and property number
c     moist    (dp,sc,in)          - current moisture
c     e        (dp,sc,in)          - young's modulus
c     posn     (dp,sc,in)          - Poisson's ratio
c
c  output arguments:
c     numpt    (int,sc,out)        - NPTS defined by TB command
c     ckh      (dp,ar(*),out)      - plastic strains
c     gkh      (dp,ar(*),out)      - input stress
c
c  table layout:
c     t1, k, (C1,g1), (C2,g2), (C3,g3), (C4,g4), (C5,g5)
c     t2, k, (C1,g1), (C2,g2), (C3,g3), (C4,g4), (C5,g5)
c     ...................
c
#include "impcom.inc"
c
      double precision ZERO, Temptol, ONE
      parameter       (ZERO = 0.0d0, 
     &                 ONE = 1.0d0,
     &                 Temptol = 1.d-3)

      integer          nlgetPtrs

      integer*8        keycpL
      
      integer          numpt
      double precision moist, ckh(*), gkh(*), dckh(*), kkh

      integer          mat, i, k1, k2, numTemps, propNum,
     &                 TempFirst, TempLast, TempStride,
     &                 tbStride
      double precision propnl(*), Header(*), tratio, dtemp1
      pointer         (dPtrHeader,Header)
      pointer         (dPtrTable,propnl)

c --- get the pointer to the start of the table of property data
      call MakeShortInt(keycpL,mat,propNum)
      i = nlgetPtrs(mat,propNum,dPtrHeader,dPtrTable)
      
      numTemps = nint(Header(2))
      numpt = nint(Header(15))

      call v_zero(dckh(1),numpt)

c --- get the table for a given moisture
      TempFirst = 1
      TempStride = 1
      tbStride = numTemps
      TempLast = TempFirst+TempStride*(numTemps-1)
      if (numTemps.gt.1) then
          if (moist.le.propnl(TempFirst)) then
              k1 = TempFirst+tbStride
              kkh = propnl(k1)
              do i = 1,numpt
                 k1 = k1+tbStride
                 ckh(i) = propnl(k1)
                 k1 = k1+tbStride
                 gkh(i) = propnl(k1)
              enddo
          elseif (moist.ge.propnl(TempLast)) then
              k1 = TempLast+tbStride
              kkh = propnl(k1)
              do i = 1,numpt
                 k1 = k1+tbStride
                 ckh(i) = propnl(k1)
                 k1 = k1+tbStride
                 gkh(i) = propnl(k1)
              enddo
          else
              k1 = TempFirst -TempStride
              do i = 1,numTemps-1
                 k1 = k1 + TempStride
                 k2 = k1 + TempStride
                 if (moist.gt.propnl(k1).and.moist.le.propnl(k2))goto 10
              enddo
 10           continue
              tratio = propnl(k2)-propnl(k1)
              dtemp1 = ZERO
              if (tratio.gt.Temptol) then
                  dtemp1 = ONE/tratio
                  tratio = (moist-propnl(k1))/tratio
              endif
              k1 = k1+tbStride
              k2 = k1+1
              kkh = propnl(k1) + tratio*(propnl(k2)-propnl(k1))
              k1 = k1+tbStride
              k2 = k1+1
              do i = 1,numpt
                 ckh(i) = propnl(k1)+tratio*(propnl(k2)-propnl(k1))
                 dckh(i) = (propnl(k2)-propnl(k1))*dtemp1/ckh(i)
                 k1 = k1+tbStride
                 k2 = k1+1
                 gkh(i) = propnl(k1)+tratio*(propnl(k2)-propnl(k1))
                 k1 = k1+tbStride
                 k2 = k1+1
              enddo
          endif
      else
          k1 = TempFirst+1
          kkh = propnl(k1)
          k1 = k1+1
          do i = 1,numpt
             ckh(i) = propnl(k1)
             gkh(i) = propnl(k1+1)
             k1 = k1+2
          enddo
      endif

      return
      end
