      recursive subroutine vi_sort(ivect,n)
c --- primary function: integer vector sort routine
c
c  input arguments:
c     n     (int,sc,in)       - number of elements to sort
c
c  output arguments:
c     ivect (int,ar(*),inout) - integer vector to sort
c
#include "impcom.inc"
c      
      integer  n, ivect(n)

      call r_sort(ivect(1),1,n)
      
      return

      contains

            recursive subroutine r_sort(a, first, last)
c --- primary function: recursive sort
c
c
            integer a(*), x, t
            integer first, last
            integer i, j

            x = a( (first+last) / 2 )
            i = first
            j = last
            do
                do while (a(i).lt.x)
                    i = i + 1
                enddo
                do while (x.lt.a(j))
                    j = j - 1
                enddo
                if (i.ge.j) exit
                t = a(i);  a(i) = a(j);  a(j) = t
                i = i + 1
                j = j - 1
            enddo
            if (first.lt.i-1) call r_sort(a(1), first, i-1)
            if (j+1.lt.last)  call r_sort(a(1), j+1, last)
      
            return
            end

      end

      