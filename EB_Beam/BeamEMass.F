      subroutine BeamEMass (dens,ex,gxy,area,aiz,aiy,shearz,sheary,
     &                      aleng,epinit,krm,nr,zass,lumpm)
c --- primary function: compute mass matrix for an euler beam
c
c     ref 1: eqn. 11.32, przemieniecki, page 294
c     ref 2: "vibrations of a hanging timoshenko beam under gravity"
c             t. yokoyama, journal of sound and vibration, 1990
c             141 (2), 245-258
c     note:  set zass matrix to zero before call to this routine
c
c  input arguments:
c     dens     (dp,sc,in)        - density
c     ex       (dp,sc,in)        - young's modulus
c     gxy      (dp,sc,in)        - shear modulus
c     area     (dp,sc,in)        - cross-sectional area
c     aiz      (dp,sc,in)        - z-axis area moment of inertia
c     aiy      (dp,sc,in)        - y-axis area moment of inertia
c     shearz   (dp,sc,in)        - z-axis shear constant
c     sheary   (dp,sc,in)        - y-axis shear constant
c     aleng    (dp,sc,in)        - length of beam
c     epinit   (dp,sc,in)        - initial strain
c     krm      (int,sc,in)       - reduced mass key (0 - consistent)
c                                                   (1 - reduced)
c     nr       (int,sc,in)       - dimension of zass matrix
c     lumpm    (int,sc,in)       - key to form the lumped matrix
c
c  output arguments:
c     zass     (dp,ar(nr,nr),out) - mass matrix
c
#include "impcom.inc"
c
      integer  krm,nr,lumpm
      double precision
     & dens,ex,gxy,area,aiz,aiy,shearz,sheary,aleng,epinit,
     & zass(nr,nr), aleng2,cy,cy2,ucy,cz,cz2,ucz,emass,con,
     & emassy,emassz,r2y,r2z
c
c
      emass = dens*area*aleng*(1.0d0 - epinit)
      if (lumpm .eq. 0) then
c --- consistent mass matrix
          zass(1,1) = emass/3.0d0
          zass(7,1) = emass/6.0d0
          zass(1,7) = zass(7,1)
          zass(7,7) = zass(1,1)
          if (krm .eq. 1) then
             zass(2,2) = zass(1,1)
             zass(3,3) = zass(1,1)
             zass(8,8) = zass(7,7)
             zass(9,9) = zass(7,7)
             zass(2,8) = zass(1,7)
             zass(3,9) = zass(1,7)
             zass(8,2) = zass(2,8)
             zass(9,3) = zass(3,9)
          else
             aleng2 = aleng*aleng
             con = 12.0d0*ex/(gxy*area*aleng2)
             cy = con*sheary*aiy
             cy2 = cy**2
             ucy = (1.0d0 + cy)**2
             emassy = emass/ucy
             r2y = aiy/area
             cz = con*shearz*aiz
             cz2 = cz**2
             ucz = (1.0d0 + cz)**2
             emassz = emass/ucz
             r2z = aiz/area
             zass(2,2) = (13.0d0/35.0d0 +  0.7d0*cz + 1.0d0/3.0d0*cz2 +
     &                    1.2d0*r2z/aleng2)*emassz
             zass(8,8) = zass(2,2)
             zass(3,3) = (13.0d0/35.0d0 +  0.7d0*cy + 1.0d0/3.0d0*cy2 +
     &                    1.2d0*r2y/aleng2)*emassy
             zass(9,9) = zass(3,3)
             zass(4,4) = (aiz + aiy)*emass/(3.0d0*area)
             zass(10,4) = zass(4,4)*0.5d0
             zass(4,10) = zass(10,4)
             zass(10,10) = zass(4,4)
             zass(5,5) = ((1.0d0/105.0d0+cy/60.0d0+cy2/120.0d0)*aleng2 +
     &                (2.0d0/15.0d0 + cy/6.0d0 + cy2/3.0d0)*r2y)*emassy
             zass(6,6) = ((1.0d0/105.0d0+cz/60.0d0+cz2/120.0d0)*aleng2 +
     &                (2.0d0/15.0d0 + cz/6.0d0 + cz2/3.0d0)*r2z)*emassz
             zass(11,11) = zass(5,5)
             zass(12,12) = zass(6,6)
             zass(6,2) = ((11.0d0/210.0d0+cz*11.0d0/120.0d0+cz2/24.0d0)*
     &                   aleng + (0.1d0 - cz*0.5d0)*r2z/aleng)*emassz
             zass(2,6) = zass(6,2)
             zass(5,3) =-((11.0d0/210.0d0+cy*11.0d0/120.0d0+cy2/24.0d0)*
     &                    aleng + (0.1d0 - cy*0.5d0)*r2y/aleng)*emassy
             zass(3,5) = zass(5,3)
             zass(11,9) = - zass(5,3)
             zass(9,11) = zass(11,9)
             zass(12,8) = - zass(6,2)
             zass(8,12) = zass(12,8)
             zass(8,2) = (9.0d0/70.0d0 + cz*0.3d0  + cz2/6.0d0 -
     &                    1.2d0*r2z/aleng2)*emassz
             zass(2,8) = zass(8,2)
             zass(9,3) = (9.0d0/70.0d0 + cy*0.3d0  + cy2/6.0d0 -
     &                    1.2d0*r2y/aleng2)*emassy
             zass(3,9) = zass(9,3)
             zass(11,3) = ((13.0d0/420.0d0+cy*3.0d0/40.0d0+cy2/24.0d0)*
     &                    aleng - (0.1d0 - cy*0.5d0)*r2y/aleng)*emassy
             zass(3,11) = zass(11,3)
             zass(12,2) =-((13.0d0/420.0d0+cz*3.0d0/40.0d0+cz2/24.0d0)*
     &                    aleng - (0.1d0 - cz*0.5d0)*r2z/aleng)*emassz
             zass(2,12) = zass(12,2)
             zass(9,5) = - zass(11,3)
             zass(5,9) = zass(9,5)
             zass(8,6) = - zass(12,2)
             zass(6,8) = zass(8,6)
             zass(11,5) =-((1.0d0/140.0d0+cy/60.0d0+cy2/120.0d0)*aleng2+
     &                 (1.0d0/30.0d0 + cy/6.0d0 - cy2/6.0d0)*r2y)*emassy
             zass(5,11) = zass(11,5)
             zass(12,6) =-((1.0d0/140.0d0+cz/60.0d0+cz2/120.0d0)*aleng2+
     &                 (1.0d0/30.0d0 + cz/6.0d0 - cz2/6.0d0)*r2z)*emassz
             zass(6,12) = zass(12,6)
         endif
      else
c --- lumped mass matrix
          con = emass*0.5d0
          call v_initi (zass(1,1),13,3,con)
          call v_initi (zass(7,7),13,3,con)
c --- no rotational effects in this case because subsequent local to global
c     transformation puts terms in off-diagonal
          if (krm .eq. 0) then
              zass(4,4) = (aiz + aiy)*emass/(2.0d0*area)
              zass(5,5) = zass(4,4)
              zass(6,6) = zass(4,4)
              zass(10,10) = zass(4,4)
              zass(11,11) = zass(4,4)
              zass(12,12) = zass(4,4)
          endif
      endif
      return
      end
