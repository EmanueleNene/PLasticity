      subroutine GetMLIN (keyplL, moist, peeq, sigy, slope)
c --- primary function: get the yield stress for a multilinear curve table
c
c  input arguments:
c     keyplL   (int*8,sc,in)       - key indicating the material and property number
c     moist    (dp,sc,in)          - current moisture
c     peeq     (dp,sc,in)          - accumulated effective strain
c
c  output arguments:
c     sigy     (dp,sc,out)         - current yield stress
c     slope    (dp,sc,out)         - current plastic slope
c
c  table layout:
c     t1 e1,s1 e2,s2 ... enpts,snpts      (not all npts may be used)
c     t2 e1,s1 e2,s2 ... enpts,snpts
c      ...
c     tnmoist e1,s1 e2,s2 ... enpts,snpts
c
#include "impcom.inc"
c
      integer          k, kt1, kt2, nmoist, dtpmoist, npts, 
     &                 nStride, tloc
      integer*8        keyplL
      double precision moist, peeq, sigy, sigy2, mratio, slope,
     &                 con1, sigy1, sigy3, Header(*)
      pointer         (dPtrHeader,Header)
c --- from ansysdef.inc in v15
      integer    NLTABHEADLENG
      parameter (NLTABHEADLENG=15)

      slope = 0.0d0
      
c --- get the pointer to the start of the table of property data
      call nlgetPtrFromKey(keyplL,dPtrHeader)
     
c --- retrieve the number of moistures from the header record
      nmoist = nint(Header(2))
c --- retrieve the number of points from the header record
      dtpmoist = nint(Header(4))
      npts = dtpmoist/2
c --- get the moisture location
      nStride = dtpmoist + 1
c --- skip past the header and the user-input data to the temps
      k = NLTABHEADLENG
      call GetTabRow (Header(k+1),nmoist,nStride,moist,tloc,mratio)
      nStride = 2
c --- at first moisture
      kt1 = k + (tloc-1)*(dtpmoist+1) + 1
      call GetTabCol (Header(kt1+1),npts,nStride,peeq,sigy,slope)
c --- at second moisture if needed
      if (mratio .ne. 0.0d0) then
          kt2 = kt1 + dtpmoist + 1
          call GetTabCol (Header(kt2+1),npts,nStride,peeq,sigy2,slope)
          sigy = sigy  + mratio*(sigy2 - sigy)
          con1 = peeq+1.d-6
          call GetTabCol(Header(kt1+1),npts,nStride,con1,sigy1,slope)
          call GetTabCol(Header(kt2+1),npts,nStride,con1,sigy2,slope)
          sigy3= sigy1 + mratio*(sigy2 - sigy1)
          con1 = peeq-1.d-6
          if (con1 .gt. 0.d0) then
              call GetTabCol(Header(kt1+1),npts,nStride,con1,
     &                       sigy1,slope)
              call GetTabCol(Header(kt2+1),npts,nStride,con1,
     &                       sigy2,slope)
              sigy1 = sigy1 + mratio*(sigy2-sigy1)
              slope = (sigy3-sigy1)*0.5d+6
          else
              slope = (sigy3-sigy)*1.0d+6
          endif
      endif
      
      return
      end
