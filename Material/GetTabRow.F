      subroutine GetTabRow (Table,nmax,nStride,value,iloc,ratio)
c --- primary function: find the position of a value in a table
c                       (assumes the table may contain less than nmax values)
c      
c  input arguments:
c     Table    (dp,ar(nmax),in)    - the table data vector
c     nmax     (int,sc,in)         - maximum number of values in vector
c     nStride  (int,sc,in)         - stride in the vector
c     value    (dp,sc,in)          - value for which the position is required
c
c  output arguments:
c     iloc     (int,sc,out)        - index below the value location
c     ratio    (dp,sc,out)         - the ratio between the two values
c
#include "impcom.inc"
c
      integer           nmax, nStride, iloc, i
      double precision  Table(nStride,nmax), value, ratio, vb, con 

c --- check for only 1 value in table
      iloc = 1
      ratio = 0.0d0
      if (nmax .le. 1) goto 999 

c --- check for below the table
      vb = Table(1,1)
      if (value .le. vb) goto 999

c --- loop thru the table and find the location
      do i = 2,nmax
c --- retrieve the table value
         con = Table(1,i)
c --- check for end of the table: value is off the end
         if (con .le. vb) then
            iloc = i - 1
            ratio = 0.0d0
            goto 999
         endif
c --- check if we found the location
         if (value .lt. con) goto 20
         vb = con
      enddo

c --- value is past end of table
      iloc = nmax
      ratio = 0.0d0
      goto 999

c --- location found within table
  20  iloc = i - 1
      ratio = (value - vb)/(con - vb)

 999  continue
      
      return
      end
