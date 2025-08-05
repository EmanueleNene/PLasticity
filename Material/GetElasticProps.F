      subroutine GetElasticProps (elem,matin,moist,prop)
c --- primary function: evaluate material properites
c --- secondary function: adjust for orthotropic input
c
c  input arguments:
c     elem     (int,sc,in)         - element number
c     matin    (int,sc,in)         - material reference number
c                                    if negative, no required props
c     moist      (dp,sc,in)        - moisture for evaluation
c
c  output arguments:
c     prop     (dp,ar(13),out)     - material properties following order:
c                                    ex,ey,ez,gxy,gyz,gxz,
c                                    uxy,uyz,uxz,alpx,alpy,alpz,dens
c
#include "impcom.inc"
c
      integer           elem, matin, mat, lp(13), i
      double precision  moist,prop(13)

      data lp /1,2,3, 7,8,9, 4,5,6, 10,11,12, 13/

c --- evaluate properties

      mat = abs(matin)
      call GetMatProps (elem,mat,lp(1),moist,prop(1),13)

c --- define shear moduli defaults
c --- check for undefined mat props
      if (prop(1).gt.0.0d0) then
          if (prop(4).eq.0.0d0) prop(4) = (prop(1)*prop(2))/
     &       (prop(1) + prop(2) + 2.0d0*prop(7)*prop(1))
          if (prop(5).eq.0.0d0) prop(5) = prop(4)
          if (prop(6).eq.0.0d0) prop(6) = prop(4)
      endif


      return
      end
