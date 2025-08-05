      subroutine BeamEMass2D (dens,ex,gxy,area,aiz,shearz,aleng,
     &                        epinit,krm,nr,zass,lumpm)
c --- primary function:  compute mass matrix for 2-d beams
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
c     shearz   (dp,sc,in)        - z-axis shear constant
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
     & dens,ex,gxy,area,aiz,shearz,aleng,epinit,
     & zass(nr,nr),aleng2,cz,cz2,ucz,emass,con,emassz,r2z
c
c
      emass = dens*area*aleng*(1.0d0 - epinit)
      if (lumpm .eq. 0) then
          zass(1,1) = emass/3.0d0
          zass(4,1) = emass/6.0d0
          zass(1,4) = zass(4,1)
          zass(4,4) = zass(1,1)
          if (krm .eq. 1) then
              zass(2,2) = zass(1,1)
              zass(5,5) = zass(4,4)
              zass(2,5) = zass(1,4)
              zass(5,2) = zass(2,5)
          else
              aleng2 = aleng*aleng
              con = 12.0d0*ex/(gxy*area*aleng2)
              cz = con*shearz*aiz
              cz2 = cz**2
              ucz = (1.0d0 + cz)**2
              emassz = emass/ucz
              r2z = aiz/area
              zass(2,2) = (13.0d0/35.0d0 +  0.7d0*cz + 1.0d0/3.0d0*cz2 +
     &                     1.2d0*r2z/aleng2)*emassz
              zass(5,5) = zass(2,2)
              zass(3,3) = ((1.0d0/105.0d0+cz/60.0d0+cz2/120.0d0)*aleng2+
     &                 (2.0d0/15.0d0 + cz/6.0d0 + cz2/3.0d0)*r2z)*emassz
              zass(6,6) = zass(3,3)
              zass(3,2) = ((11.0d0/210.0d0+cz*11.0d0/120.0d0+cz2/24.0d0)
     &                  * aleng + (0.1d0 - cz*0.5d0)*r2z/aleng)*emassz
               zass(2,3) = zass(3,2)
               zass(6,5) = - zass(3,2)
               zass(5,6) = zass(6,5)
               zass(5,2) = (9.0d0/70.0d0 + cz*0.3d0  + cz2/6.0d0 -
     &                      1.2d0*r2z/aleng2)*emassz
               zass(2,5) = zass(5,2)
               zass(6,2) =-((13.0d0/420.0d0+cz*3.0d0/40.0d0+cz2/24.0d0)*
     &                      aleng - (0.1d0 - cz*0.5d0)*r2z/aleng)*emassz
               zass(2,6) = zass(6,2)
               zass(5,3) = - zass(6,2)
               zass(3,5) = zass(5,3)
               zass(6,3) =-((1.0d0/140.0d0+cz/60.0d0+cz2/120.0d0)*aleng2
     &                +(1.0d0/30.0d0 + cz/6.0d0 - cz2/6.0d0)*r2z)*emassz
            zass(3,6) = zass(6,3)
         endif
      else
c --- lumped mass matrix
          con = emass*0.5d0
          call v_initi (zass(1,1),7,2,con)
          call v_initi (zass(4,4),7,2,con)
      endif

      return
      end
