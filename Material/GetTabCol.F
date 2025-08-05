      subroutine GetTabCol (vect,nmax,nStride,peeq, sigy,slope)
c --- primary function: find the position of a value in a virtual table
c                       (assumes the table may contain less than nmax values)
c
c  input arguments:
c     vect     (dp,ar(nStride,*),in) - data table (vector)
c     nmax     (int,sc,in)           - maximum number of values in vector
c     nstrid   (int,sc,in)           - stride in the vector (in integers)
c     peeq     (dp,sc,in)            - equivalent plastic strain 
c
c  output arguments:
c     sigy     (dp,sc,out)           - yield stress
c     slope    (dp,sc,out)           - slope
c
#include "impcom.inc"
c
      integer          nmax, nStride, i, nsearch
      double precision vect(nStride,*), peeq, sigy, slope, peeq1, sigy1
      double precision TOL
      parameter       (TOL=1.0d-12)

      sigy = vect(2,1)
      slope= 0.0d0
      if (nmax .le. 1) goto 999

      peeq1 = vect(1,1)
      if (peeq1 .ne .0.0d0) then
          call erhandler ('GetTabCol', 5000, 3, 'Wrong input data.
     &                    First strain value in MISO should be zero.',
     &                    0.0d0,' ')
          call StopAll (-5)           
      endif
      
      sigy1 = sigy
      
      if (peeq .le. peeq1) then 
          slope = (vect(2,2) - sigy1) / (vect(1,2) - peeq1)
          goto 999
      endif
      
      nsearch = nmax
      do i = 2,nsearch
c --- check for end of the table
         if (vect(1,i) .le. peeq1) goto 10
c --- check if we are within this segment 
c --- use a tolerance for case when right on a breakpoint!
         if (peeq-vect(1,i) .lt. -TOL) goto 20
         sigy1 = vect(2,i)
         peeq1 = vect(1,i)
      enddo 

  10  sigy = sigy1
      slope= 0.0d0
      
      goto 999

  20  slope = (vect(2,i) - sigy1) / (vect(1,i) - peeq1)
      sigy  = sigy1 + slope*(peeq-peeq1)

 999  continue
      
      return
      end
