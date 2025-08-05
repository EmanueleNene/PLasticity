      subroutine UpdPiloNd (nPerFounI,nFoun,icoor,nodbac,nopilo,
     &                      ififo,dispn)
c --- primary function:   determine/update the current cell of foundations
c --- secondary function: update u-vector for foundations with pilot node
c
c  input arguments:
c     nPerFounI   (int,sc,in)         - number of integer items per foundation
c                                       database entry
c     nFoun    (int,sc,in)            - number of foundations
c     icoor    (dp,ar(*),in)          - initial nodal coordinate
c     nodbac   (int,ar(D_BacLen),out) - internal to external list
c
c     dispn     (dp,ar(*),inout)      - corrent nodal displacement
c
c     nopilo   (int,ar(D_BacLen),out) - pilot nodes lists
c
#include "impcom.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "solu_info.inc"
c
      integer nPerFounI,nFoun,nodbac(*),ififo(nPerFounI,nFoun),nopilo(*)
      double precision icoor(*),dispn(*)

      integer inopil,iropil,kddlx,kddly,kddlz,
     &  inopi,ino,inol,i,ifoun,kfoun,jfoun,err

      double precision ZERO,xp,yp,zp,x0,y0,xp0,yp0,z0,zp0,
     &  excent1,excent2,excent3,x,y,z,IM(9)

      double precision tgrotat(3,3),tet(3),dis0,con,duiOrig

      pointer (pdULocui0L,ULocui0)
      pointer (pdULocun0L,ULocun0)
      double precision  ULocui0(*), ULocun0(*)

      parameter (ZERO=0.0d0)

#if defined(DEBUG)
      call SubBeg('UpdPiloNd')
#endif

      pdULocui0L = M_PtrIncI 
      pdULocun0L = M_PtrIncS

c --- unit vector for rotation
      call v_zero (IM(1),9)
      call v_initi (IM(1),4,3,1.0d0)

c --- for slave surface with pilot node, update u vector

      do 180 kfoun = 1,nFoun
         ifoun = ififo(24,kfoun)
         if (ififo(1,ifoun).le.0) goto 180
         inopil = ififo(2,ifoun)
         if (((inopil-1)/D_NumDof + 1).le.0) goto 180
         iropil = abs(ififo(3,ifoun))
c --- correct position for pilot nodes
         if (iropil.gt.0) then
             tet(1) = ULocun0(inopil+3)
             tet(2) = ULocun0(inopil+4)
             tet(3) = ULocun0(inopil+5)
c --- calculate the rotation
             call BasRot (IM(1),tgrotat(1,1),tet(1),3)
         endif

         kddlx = inopil
         kddly = kddlx + 1
         kddlz = kddlx + 2
         
         if (abs(dispn(kddlx)).eq.HUGE) dispn(kddlx) = ZERO
         xp = icoor(kddlx) + dispn(kddlx)
         if (abs(dispn(kddly)).eq.HUGE) dispn(kddly) = ZERO
         yp = icoor(kddly) + dispn(kddly)
         if (abs(dispn(kddlz)).eq.HUGE) dispn(kddlz) = ZERO
         zp = icoor(kddlz) + dispn(kddlz)

         xp0 = xp - ULocun0(kddlx)
         yp0 = yp - ULocun0(kddly)
         zp0 = zp - ULocun0(kddlz)
         
         if (iropil.gt.0) then
             err = 0
             if (abs(ULocui0(inopil+3)).ge.HUGE) then
                 ULocui0(inopil+3) = ZERO
                 err = 1
             endif
             if (abs(ULocui0(inopil+4)).ge.HUGE) then
                 ULocui0(inopil+4) = ZERO
                 err = 1
             endif
             if (abs(ULocui0(inopil+5)).ge.HUGE) then
                 ULocui0(inopil+5) = ZERO
                 err = 1
             endif
             if (err.eq.1.and.A_IterNum.eq.2.and.A_SubStepNum.eq.1) then
                 err = (ififo(2,ifoun)-1)/D_NumDof + 1
                 err = nodbac(err)
                 call erhandler ('UpdPiloNd',5000,
     &                 2,'The rotational DOF of pilot node %i is not 
     &                 constrained which can cause rigid body motion.',
     &                 dble(err),' ')
             endif
         endif

c --- calculate increments of nodes on foundation with pilot node
         inopi = (inopil-1)/D_NumDof + 1
         do 175 ino = 1,D_BacLen
            if (ino .eq. inopi) goto 175
            jfoun = nopilo(ino)
            if (jfoun.le.0) goto 175
            if (jfoun.ne.ifoun) goto 175
            
            inol = (ino-1)*D_NumDof + 1
            if (iropil.le.0) then
c --- no rotation (same as pilot node)
                do 150 i = 0,2
                   dispn(inol+i) = dispn(inol+i)-ULocun0(inol+i)
                   if (abs(ULocui0(inopil+i)).ge.HUGE) then
                       if (A_IterNum.eq.2.and.A_SubStepNum.eq.1) then
                           err = (ififo(2,ifoun)-1)/D_NumDof + 1
                           err = nodbac(err)
                           call erhandler ('UpdPiloNd',5000,2,
     &                     'The translational DOF of pilot node %i is
     &                     not constrained which can cause rigid 
     &                     body motion.',dble(err),' ')
                        endif
                        ULocui0(inopil+i) = ZERO
                   endif                  
                   ULocun0(inol+i) = ULocun0(inopil+i)
                   dispn(inol+i) =  dispn(inol+i) + ULocun0(inopil+i)
 150            continue
                goto 175
            endif

c --- rotation case
            kddlx = inol
            kddly = kddlx + 1
            kddlz = kddlx + 2

            x0 = icoor(kddlx) + dispn(kddlx) - ULocun0(kddlx)
            y0 = icoor(kddly) + dispn(kddly) - ULocun0(kddly)
            z0 = icoor(kddlz) + dispn(kddlz) - ULocun0(kddlz)
            
            excent1 = x0 - xp0
            excent2 = y0 - yp0
            excent3 = z0 - zp0

            x = xp + tgrotat(1,1)*excent1
     &             + tgrotat(1,2)*excent2
     &             + tgrotat(1,3)*excent3
            y = yp + tgrotat(2,1)*excent1
     &             + tgrotat(2,2)*excent2
     &             + tgrotat(2,3)*excent3
            z = zp + tgrotat(3,1)*excent1
     &             + tgrotat(3,2)*excent2
     &             + tgrotat(3,3)*excent3
            
            if (ULocui0(kddlx).eq.HUGE) ULocui0(kddlx) = ZERO
            dis0 = dispn(kddlx) - ULocui0(kddlx)
            dispn(kddlx) = x - icoor(kddlx)
            ULocui0(kddlx) = dispn(kddlx) - dis0
            ULocun0(kddlx) = x - x0
 
            if (ULocui0(kddly).eq.HUGE) ULocui0(kddly) = ZERO
            dis0 = dispn(kddly) - ULocui0(kddly)
            dispn(kddly) = y - icoor(kddly)
            ULocui0(kddly) = dispn(kddly) - dis0
            ULocun0(kddly) = y - y0
 
            if (ULocui0(kddlz).eq.HUGE) ULocui0(kddlz) = ZERO
            dis0 = dispn(kddlz) - ULocui0(kddlz)
            dispn(kddlz) = z - icoor(kddlz)
            ULocui0(kddlz) = dispn(kddlz) - dis0
            ULocun0(kddlz) = z - z0
 
 175     continue
 180  continue

 999  continue

#if defined(DEBUG)
      call SubEnd('UpdPiloNd')
#endif

      return
      end