      subroutine GetMatProps (elem,mat,lp,moist,prop,n)
c --- primary function: to evaluate a group of material properties
c
c  input arguments:
c     elem  (int,sc,in)       - element number
c     mat   (int,sc,in)       - material number(input quantity mat, mat comma
c     lp    (int,ar(n),in)    - keys for which specific value is requested
c     moist (dp,sc,in)        - moisture at which to evaluate material
c     n     (int,sc,in)       - number of properties to be evaluated
c
c  output arguments:
c     prop     (dp,ar(n),out) - values of material property
c
c     GetMatProps is used to several material property defined by lp array
c     Can be called as as:
c        integer lp(3)
c        data lp /1,10,13/
c        call GetMatProps (elem,mat,lp(1),moist,prop(1),3)
c
c --- MP command labels
c        EX  = 1, EY  = 2, EZ  = 3, NUXY= 4, NUYZ= 5, NUXZ= 6, GXY = 7, GYZ = 8,
c        GXZ = 9, ALPX=10, ALPY=11, ALPZ=12, DENS=13, MU  =14, MUX =16, MUY = 14
c
#include "impcom.inc"
#include "constants.inc"
c
      integer  elem,mat,n,lp(n),i
      double precision  prop(n),moist
      
      if (n.gt.1) then
          do i = 1,n
             call prope1 (elem,mat,lp(i),moist,prop(i))
          enddo
      else
          call prope1 (elem,mat,lp(1),moist,prop(1))
      endif      

      return
      end

