      subroutine Axv1_v2 (a,v1,v2, nr,nc)
c --- primary function: multiply a vector by a matrix
c                       if  nc.gt.0 the result is summed to v2
c
c  input arguments:
c     a        (dp,ar(nr,*),in)  - matrix a
c     v1       (dp,ar(*),in)     - vector v1
c     nr       (int,sc,in)       - number of rows in matrix a
c     nc       (int,sc,in)       - number of columns to multiply in matrix a
c
c  output arguments:
c     v2       (dp,ar(*),out)    - product vector v2
c
#include "impcom.inc"

      integer  nr, nc, n, i, j, irem
      double precision   a(nr,*),v1(*),v2(nr)
 

      n = - nc
      if (n .le. 0) then
c --- result is not accumulated with input values
         n = nc

         irem = mod(n,8)
         goto (80,10,20,30,40,50,60,70), irem+1

 10      do i=1,nr
            v2(i) = a(i,1) * v1(1)
         enddo
         goto 99

 20      do i=1,nr
            v2(i) = a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
         enddo
         goto 99

 30      do i=1,nr
            v2(i) = a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
         enddo
         goto 99

 40      do i=1,nr
            v2(i) = a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
         enddo
         goto 99

 50      do i=1,nr
            v2(i) = a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
     &            + a(i,5) * v1(5)
         enddo
         goto 99


 60      do i=1,nr
            v2(i) = a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
     &            + a(i,5) * v1(5)
     &            + a(i,6) * v1(6)
         enddo
         goto 99

 70      do i=1,nr
            v2(i) = a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
     &            + a(i,5) * v1(5)
     &            + a(i,6) * v1(6)
     &            + a(i,7) * v1(7)
         enddo
         goto 99

 80      do i=1,nr
            v2(i) = a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
     &            + a(i,5) * v1(5)
     &            + a(i,6) * v1(6)
     &            + a(i,7) * v1(7)
     &            + a(i,8) * v1(8)
         enddo
         irem = 8
         goto 99

      else
c --- result is accumulated with input values
         irem = mod(n,8)
         goto (99,110,120,130,140,150,160,170), irem+1

 110     do i=1,nr
            v2(i) = v2(i)
     &            + a(i,1) * v1(1)
         enddo
         goto 99

 120     do i=1,nr
            v2(i) = v2(i)
     &            + a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
         enddo
         goto 99

 130     do i=1,nr
            v2(i) = v2(i)
     &            + a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
         enddo
         goto 99

 140     do i=1,nr
            v2(i) = v2(i)
     &            + a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
         enddo
         goto 99

 150     do i=1,nr
            v2(i) = v2(i)
     &            + a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
     &            + a(i,5) * v1(5)
         enddo
         goto 99

 160     do i=1,nr
            v2(i) = v2(i)
     &            + a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
     &            + a(i,5) * v1(5)
     &            + a(i,6) * v1(6)
         enddo
         goto 99

 170     do i=1,nr
            v2(i) = v2(i)
     &            + a(i,1) * v1(1)
     &            + a(i,2) * v1(2)
     &            + a(i,3) * v1(3)
     &            + a(i,4) * v1(4)
     &            + a(i,5) * v1(5)
     &            + a(i,6) * v1(6)
     &            + a(i,7) * v1(7)
         enddo

      endif

c --- finish up the remaining number of columns
 99      do j=irem+1,n,8
            do i=1,nr
               v2(i) = v2(i)
     &               + a(i,j  ) * v1(j  )
     &               + a(i,j+1) * v1(j+1)
     &               + a(i,j+2) * v1(j+2)
     &               + a(i,j+3) * v1(j+3)
     &               + a(i,j+4) * v1(j+4)
     &               + a(i,j+5) * v1(j+5)
     &               + a(i,j+6) * v1(j+6)
     &               + a(i,j+7) * v1(j+7)
            enddo
         enddo


      return
      end
