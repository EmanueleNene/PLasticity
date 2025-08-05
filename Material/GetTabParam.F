      subroutine GetTabParam (keyL, moist, numpm, prop, option)
c --- primary function: obtain material parameters from a table
c
c  input arguments:
c     keyL     (int*8,sc,in)      - key indicating the material and property number
c     moist    (dp,sc,in)         - current temperature
c
c  output arguments:
c     numpm   (int,sc,out)         - NPTS defined by TB command
c     prop    (dp,ar(*),out)       - material properties
c     option  (in,sc, out)         - option from the TB command
c
#include "impcom.inc"
c
      double precision MoistTol
      parameter       (MoistTol=1.0d-3)

      integer          nlgetPtrs

      integer          numpm, option
      integer*8        keyL
      double precision moist, prop(*)

      integer          i, k1, k2, numMoists, mat, propNum,
     &                 MoistFirst, MoistLast, MoistStride, numPoints,
     &                 tbStride, locMoist, plopt
      double precision propnl(*), Header(*), mratio
      pointer         (dPtrHeader,Header)
      pointer         (dPtrTable,propnl)

#if defined(DEBUG)
      call SubBeg('GetTabParam')
#endif

c --- get the pointer to the start of the table of property data
      call MakeShortInt(keyL,mat,propNum)
      i = nlgetPtrs(mat,propNum,dPtrHeader,dPtrTable)
      
      if (i.eq.0) then
          option = 0
          numpm = 0
          goto 999
      endif

      plopt = nint(Header(1))
      numMoists = nint(Header(2))
      locMoist = nint(Header(3))
      numPoints = nint(Header(4))
      option = nint(Header(13))

      numpm = numPoints

c --- N.B. this is a patch for version 13 +
      if (numpm.eq.0 .and. option.gt.0 .and. plopt.eq.9) then
          numpm = 6
          MoistFirst = 42
      else
          MoistFirst = 1
      endif

c --- get the table for a given temperature      
      MoistStride = 1
      tbStride = numMoists
      MoistLast = MoistFirst+MoistStride*(numMoists-1)

      if (numMoists.gt.1) then
          if (moist.le.propnl(MoistFirst)) then
              k1  = MoistFirst + tbStride
              do i = 1,numpm
                 prop(i) = propnl(k1)
                 k1 = k1+tbStride
              enddo
          elseif (moist.ge.propnl(MoistLast)) then
              k1 = MoistLast + tbStride
              do i = 1,numpm
                 prop(i) = propnl(k1)
                 k1 = k1 + tbStride
              enddo
          else
              k1 = MoistFirst - MoistStride
              do i = 1,numMoists-1
                 k1 = k1 + MoistStride
                 k2 = k1 + MoistStride
                 if (moist.gt.propnl(k1).and.moist.le.propnl(k2))goto 10
              enddo
 10           continue
              mratio = propnl(k2)-propnl(k1)
              if (mratio.gt.MoistTol) mratio = (moist-propnl(k1))/mratio
              k1 = k1 + tbStride
              k2 = k1 + MoistStride
              do i = 1,numpm
                 prop(i) = propnl(k1) + mratio*(propnl(k2)-propnl(k1))
                 k1 = k1 + tbStride
                 k2 = k1 + MoistStride
              enddo
          endif
      else
          k1 = MoistFirst + locMoist
          do i = 1,numpm
             prop(i) = propnl(k1)
             k1 = k1+1
          enddo
      endif
      
  999 continue

#if defined(DEBUG)
      call SubEnd('GetTabParam')
#endif

      return
      end
