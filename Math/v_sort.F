      subroutine v_sort(vect,n,indx)
c --- primary function: integer vector sort routine
c
c  input arguments:
c     n     (int,sc,in)       - number of elements to sort
c
c  output arguments:
c     vect  (dp,ar(*),inout)  - vector to sort
c     indx  (int,ar(*),inout) - sort index
c
#include "impcom.inc"
c      
      double precision vect(n), hold
      integer  n, indx(n), i, j, k, l, m, ihold
      
      do i = 1,n
         indx(i) = i
      enddo
      
      m = n
 10   m = m/2
         if (m .eq. 0) goto 999
         k = n - m
         j = 1
 20      i = j
 30         l = i + m
            if (vect(i) .le. vect(l)) goto 40

            ihold  = indx(i)
            indx(i) = indx(l)
            indx(l) = ihold

            hold  = vect(i)
            vect(i) = vect(l)
            vect(l) = hold

            i = i - m
            if (i .ge. 1) goto 30
 40      j = j + 1
         if (j .gt. k) goto 10
      goto 20

 999  return
      end