      subroutine vi_intersect (na,ia, nb,ib, nr,ir)
c --- primary function: form the intersection of two vectors 
c                                                                 
c  input arguments:
c     na       (int,sc,in)       - number of terms in vector ia
c     ia       (int,ar(na),in)   - the first vector ia (sorted)
c     nb       (int,sc,in)       - number of terms in vector ib
c     ib       (int,ar(nb),in)   - the second vector ib (sorted)
c                                             
c  output arguments:
c     nr       (int,sc,out)      - number of terms in union vector
c     ir       (int,ar(nr),out)  - the intersection vector (sorted)
c                                  (must allocate max(na,nb))
c                                             
#include "impcom.inc"
c                     
      integer  na, ia(na), nb, ib(nb), nr, ir(*), ja, jb


      ja = 1
      jb = 1
      nr = 0

  10     if (ja .gt. na) goto 999
  15     if (jb .gt. nb) goto 999

         if (ia(ja) - ib(jb)) 20,30,40

  20     ja  = ja + 1 
         goto 10

  30     nr = nr + 1
         ir(nr) = ia(ja)
         ja = ja + 1
         jb = jb + 1
         goto 10

  40     jb = jb + 1
         goto 15
                                 
 999  return
      end
