      subroutine Atxv1_v2 (a,v1,v2, nr,nc)
c --- primary function: multiply vector by full transposed matrix
c                       if  nc.gt.0 the result is summed to v2
c
c  input arguments:
c     a        (dp,ar(nr,*),in)   - matrix a (first dimension must = nr)
c     v1       (dp,ar(nv),in)     - vector v1 (nv must be greater or equal to nr)
c     nr       (int,sc,in)        - first dimension and number of active
c                                   rows of the untransposed matrix a
c                                   (also the number of active rows of vector v1)
c     nc       (int,sc,in)        - number of columns of the untransposed matrix a
c                                   (also the number of computed items in the product vector v2)
c                                   if negative, accumulate
c
c  output arguments:
c     v2       (dp,ar(na,*),out)  - product vector v2
c
#include "impcom.inc"
c
      integer  nr, nc, n, i, j, irem
      double precision  a(nr,*), v1(nr), v2(*)
 
      n = - nc
      if (n .le. 0) then
c --- result is not accumulated with input values
         n = nc

         irem = mod(nr,8)
         goto (80,10,20,30,40,50,60,70) irem+1

 10      do i=1,n
            v2(i) = a(1,i) * v1(1)
         enddo
         goto 99

 20      do i=1,n
            v2(i) = a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
         enddo
         goto 99

 30      do i=1,n
            v2(i) = a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
         enddo
         goto 99

 40      do i=1,n
            v2(i) = a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
         enddo
         goto 99

 50      do i=1,n
            v2(i) = a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
         enddo
         goto 99


 60      do i=1,n
            v2(i) = a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
     &            + a(6,i) * v1(6)
         enddo
         goto 99

 70      do i=1,n
            v2(i) = a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
     &            + a(6,i) * v1(6)
     &            + a(7,i) * v1(7)
         enddo
         goto 99

 80      do i=1,n
            v2(i) = a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
     &            + a(6,i) * v1(6)
     &            + a(7,i) * v1(7)
     &            + a(8,i) * v1(8)
         enddo
         irem=8
         goto 99

      else
c --- result is accumulated with input values
         irem = mod(nr,8)
         goto (99,110,120,130,140,150,160,170,180) irem+1

 110     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
         enddo
         goto 99

 120     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
         enddo
         goto 99

 130     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
         enddo
         goto 99

 140     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
         enddo
         goto 99

 150     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
         enddo
         goto 99

 160     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
     &            + a(6,i) * v1(6)
         enddo
         goto 99

 170     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
     &            + a(6,i) * v1(6)
     &            + a(7,i) * v1(7)
         enddo
         goto 99

 180     do i=1,n
            v2(i) = v2(i)
     &            + a(1,i) * v1(1)
     &            + a(2,i) * v1(2)
     &            + a(3,i) * v1(3)
     &            + a(4,i) * v1(4)
     &            + a(5,i) * v1(5)
     &            + a(6,i) * v1(6)
     &            + a(7,i) * v1(7)
     &            + a(8,i) * v1(8)
         enddo
         goto 99

      endif

c --- finish up the remaining number of columns
 99      do j=irem+1,nr,8
            do i=1,n
               v2(i) = v2(i)
     &               + a(j  ,i) * v1(j  )
     &               + a(j+1,i) * v1(j+1)
     &               + a(j+2,i) * v1(j+2)
     &               + a(j+3,i) * v1(j+3)
     &               + a(j+4,i) * v1(j+4)
     &               + a(j+5,i) * v1(j+5)
     &               + a(j+6,i) * v1(j+6)
     &               + a(j+7,i) * v1(j+7)
          enddo
         enddo

      return
      end
