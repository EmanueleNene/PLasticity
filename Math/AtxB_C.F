      subroutine AtxB_C (a,b,c, na,nb,nc, n1,n2,n3)
c --- primary function:  (a)t * (b) = (c)     t means transpose
c                        if  nc.gt.0 the result is summed to c
c
c  input arguments:
c     a        (dp,ar(na,*),in)   - matrix a
c     b        (dp,ar(nb,*),in)   - matrix b
c     na       (int,sc,in)        - number of rows in matrix a
c     nb       (int,sc,in)        - number of rows in matrix b
c     nc       (int,sc,in)        - number of rows in matrix c
c     n1       (int,sc,in)        - number of rows in matrix c to fill
c     n2       (int,sc,in)        - number of columns in matrix c to fill
c     n3       (int,sc,in)        - number of columns in matrix a and
c                                   number of columns of matrix b
c                                   if n3 is negative, accumulate results in c
c
c  output arguments:
c     c        (dp,ar(nc,*),out)  - product matrix c
c
#include "impcom.inc"
c
      integer  na, nb, nc, n1, n2, n3, n, i, j, k, irem
      double precision   a(na,*),b(nb,*),c(nc,*)

      n = - n3
      if (n .le. 0) then

c --- new target matrix
c
         n = n3
         irem = mod(n,8)
         goto (80,10,20,30,40,50,60,70) irem + 1
        
 10      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
            enddo
         enddo
         goto 99


 20      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
     &                  + a(2,i)*b(2,j)
            enddo
         enddo
         goto 99


 30      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
     &                  + a(2,i)*b(2,j)
     &                  + a(3,i)*b(3,j)
            enddo
         enddo
         goto 99


 40      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
     &                  + a(2,i)*b(2,j)
     &                  + a(3,i)*b(3,j)
     &                  + a(4,i)*b(4,j)
            enddo
         enddo
         goto 99


 50      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
     &                  + a(2,i)*b(2,j)
     &                  + a(3,i)*b(3,j)
     &                  + a(4,i)*b(4,j)
     &                  + a(5,i)*b(5,j)
            enddo
         enddo
         goto 99


 60      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
     &                  + a(2,i)*b(2,j)
     &                  + a(3,i)*b(3,j)
     &                  + a(4,i)*b(4,j)
     &                  + a(5,i)*b(5,j)
     &                  + a(6,i)*b(6,j)
            enddo
         enddo
         goto 99


 70      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
     &                  + a(2,i)*b(2,j)
     &                  + a(3,i)*b(3,j)
     &                  + a(4,i)*b(4,j)
     &                  + a(5,i)*b(5,j)
     &                  + a(6,i)*b(6,j)
     &                  + a(7,i)*b(7,j)
            enddo
         enddo
         goto 99

 80      irem=8
         do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(1,i)*b(1,j)
     &                  + a(2,i)*b(2,j)
     &                  + a(3,i)*b(3,j)
     &                  + a(4,i)*b(4,j)
     &                  + a(5,i)*b(5,j)
     &                  + a(6,i)*b(6,j)
     &                  + a(7,i)*b(7,j)
     &                  + a(8,i)*b(8,j)
            enddo
         enddo
         goto 99

      else

         irem = mod(n,8)
         goto (99,110,120,130,140,150,160,170) irem + 1

 110     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(1,i)*b(1,j)
            enddo
         enddo
         goto 99

 120     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(1,i)*b(1,j)
     &                         + a(2,i)*b(2,j)
            enddo
         enddo
         goto 99

 130     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(1,i)*b(1,j)
     &                         + a(2,i)*b(2,j)
     &                         + a(3,i)*b(3,j)
            enddo
         enddo
         goto 99

 140     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(1,i)*b(1,j)
     &                         + a(2,i)*b(2,j)
     &                         + a(3,i)*b(3,j)
     &                         + a(4,i)*b(4,j)
            enddo
         enddo
         goto 99


 150     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(1,i)*b(1,j)
     &                         + a(2,i)*b(2,j)
     &                         + a(3,i)*b(3,j)
     &                         + a(4,i)*b(4,j)
     &                         + a(5,i)*b(5,j)
            enddo
         enddo
         goto 99

 160     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(1,i)*b(1,j)
     &                         + a(2,i)*b(2,j)
     &                         + a(3,i)*b(3,j)
     &                         + a(4,i)*b(4,j)
     &                         + a(5,i)*b(5,j)
     &                         + a(6,i)*b(6,j)
               enddo
         enddo
         goto 99

 170     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(1,i)*b(1,j)
     &                         + a(2,i)*b(2,j)
     &                         + a(3,i)*b(3,j)
     &                         + a(4,i)*b(4,j)
     &                         + a(5,i)*b(5,j)
     &                         + a(6,i)*b(6,j)
     &                         + a(7,i)*b(7,j)
            enddo
         enddo

      endif

 99   do k = irem+1,n,8
         do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(k  ,i)*b(k  ,j)
     &                         + a(k+1,i)*b(k+1,j)
     &                         + a(k+2,i)*b(k+2,j)
     &                         + a(k+3,i)*b(k+3,j)
     &                         + a(k+4,i)*b(k+4,j)
     &                         + a(k+5,i)*b(k+5,j)
     &                         + a(k+6,i)*b(k+6,j)
     &                         + a(k+7,i)*b(k+7,j)
            enddo
         enddo
      enddo

      return
      end
