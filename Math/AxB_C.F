      subroutine AxB_C (a,b,c, na,nb,nc, n1,n2,n3)
c --- primary function: multiply two double precision matrices
c                       and stores the results to matrix c if n3.gt.0
c                       or add to matrix c otherwise 
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
c                                   number of rows in matrix b
c                                   to work with (the two need
c                                   to be the same for the inner product)
c                                   if n3 is negative, accumulate in c
c
c  output arguments:
c     c        (dp,ar(nc,*),out)  - product matrix c
c
#include "impcom.inc"
c
      integer  na,nb,nc, n1,n2,n3, n, i,j,k, irem
      double precision   a(na,*),b(nb,*),c(nc,*)
 

      if (n3 .gt. 0) then

c --- new matrix
         n = n3
         irem = mod(n,8)
         goto (80,10,20,30,40,50,60,70) irem + 1
        
 10      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
            enddo
         enddo
         goto 99

 20      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
     &                  + a(i,2)*b(2,j)
            enddo
         enddo
         goto 99

 30      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
     &                  + a(i,2)*b(2,j)
     &                  + a(i,3)*b(3,j)
            enddo
         enddo
         goto 99

 40      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
     &                  + a(i,2)*b(2,j)
     &                  + a(i,3)*b(3,j)
     &                  + a(i,4)*b(4,j)
            enddo
         enddo
         goto 99

 50      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
     &                  + a(i,2)*b(2,j)
     &                  + a(i,3)*b(3,j)
     &                  + a(i,4)*b(4,j)
     &                  + a(i,5)*b(5,j)
            enddo
         enddo
         goto 99

 60      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
     &                  + a(i,2)*b(2,j)
     &                  + a(i,3)*b(3,j)
     &                  + a(i,4)*b(4,j)
     &                  + a(i,5)*b(5,j)
     &                  + a(i,6)*b(6,j)
            enddo
         enddo
         goto 99

 70      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
     &                  + a(i,2)*b(2,j)
     &                  + a(i,3)*b(3,j)
     &                  + a(i,4)*b(4,j)
     &                  + a(i,5)*b(5,j)
     &                  + a(i,6)*b(6,j)
     &                  + a(i,7)*b(7,j)
            enddo
         enddo
         goto 99

 80      irem = 8
         do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(1,j)
     &                  + a(i,2)*b(2,j)
     &                  + a(i,3)*b(3,j)
     &                  + a(i,4)*b(4,j)
     &                  + a(i,5)*b(5,j)
     &                  + a(i,6)*b(6,j)
     &                  + a(i,7)*b(7,j)
     &                  + a(i,8)*b(8,j)
            enddo
         enddo
         goto 99
         
      else


c --- result is accumulated with input values
         n = -n3
         irem = mod(n,8)
         goto (99,110,120,130,140,150,160,170) irem + 1

 110     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(1,j)
            enddo
         enddo
         goto 99

 120     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(1,j)
     &                         + a(i,2)*b(2,j)
            enddo
         enddo
         goto 99

 130     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(1,j)
     &                         + a(i,2)*b(2,j)
     &                         + a(i,3)*b(3,j)
            enddo
         enddo
         goto 99

 140     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(1,j)
     &                         + a(i,2)*b(2,j)
     &                         + a(i,3)*b(3,j)
     &                         + a(i,4)*b(4,j)
            enddo
         enddo
         goto 99

 150     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(1,j)
     &                         + a(i,2)*b(2,j)
     &                         + a(i,3)*b(3,j)
     &                         + a(i,4)*b(4,j)
     &                         + a(i,5)*b(5,j)
            enddo
         enddo
         goto 99

 160     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(1,j)
     &                         + a(i,2)*b(2,j)
     &                         + a(i,3)*b(3,j)
     &                         + a(i,4)*b(4,j)
     &                         + a(i,5)*b(5,j)
     &                         + a(i,6)*b(6,j)
            enddo
         enddo
         goto 99

 170     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(1,j)
     &                         + a(i,2)*b(2,j)
     &                         + a(i,3)*b(3,j)
     &                         + a(i,4)*b(4,j)
     &                         + a(i,5)*b(5,j)
     &                         + a(i,6)*b(6,j)
     &                         + a(i,7)*b(7,j)
            enddo
         enddo
         goto 99

      endif

  99  do k = irem+1,n,8
         do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,k  )*b(k  ,j)
     &                         + a(i,k+1)*b(k+1,j)
     &                         + a(i,k+2)*b(k+2,j)
     &                         + a(i,k+3)*b(k+3,j)
     &                         + a(i,k+4)*b(k+4,j)
     &                         + a(i,k+5)*b(k+5,j)
     &                         + a(i,k+6)*b(k+6,j)
     &                         + a(i,k+7)*b(k+7,j)
            enddo
         enddo
      enddo

      return
      end
