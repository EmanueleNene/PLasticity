      subroutine Cont_5 (NumEl,BacLen,NumDof,nodfwd,nodbac,nPerSeg,
     &                   nSeg,nPerFounI,nMast,nPerMast,nPerFounDp,nFoun,
     &                   icoor,dispn,isegda,ififo,founda,icnda,ncon)
c --- primary function: get initial gap and penetration information
c                       for each entire slave surface
c
c  input arguments:
c     BacLen   (int,sc,in)         - number of active nodes
c     NumDof   (int,sc,in)         - number of active dof per node
c     nodfwd   (int,ar(BacLen),in) - forward list for compressed nodes
c     nPerSeg  (int,sc,in)         - number of items per segment database entry
c     nSeg     (int,sc,in)         - number of segments on foundation
c     nPerMast (int,sc,in)         - number of items per contact element
c     nPerFounI(int,sc,in)         - number of integer items per foundation
c                                    database entry
c     nPerFounDp   (int,sc,in)     - number of double precision items per
c                                     foundation database entry
c     nMast   (int,sc,in)          - number of contact elements
c     nFoun    (int,sc,in)         - number of foundations
c     icoor    (dp,ar(*),in)       - initial nodal coordinate
c     dispn    (dp,ar(*),in)       - corrent nodal displacement
c
c  output arguments:
c     isegda   (int,ar(nPerSeg,nSeg),out)
c                                  - segment information
c          isegda(1,iseg)        = code of segment type (icode)
c                                = 1   : straight line
c                                = 2   : parabolic arc
c          isegda(2,iseg)        = fondation number
c          isegda(j,iseg)        = 4 nodes (j=3-10) on seg. iseg
c          isegda(j,iseg)        = 4 adjacent segments (j=11-14)
c                                = 2 adjacent segments (j=11,12)
c                                     for seg. iseg
c                                = 0  (icode > 10)
c                                  --> none adjacent segment
c          isegda(15,iseg)       = segment element number
c                                  < slave element death
c          isegda(16,iseg)       = attached beam element number, if it is
c          isegda(17,iseg)       = attached contact element number, if it is
c          isegda(18,iseg)        = remark
c
c     ififo    (dp,ar(nPerFounI,nFoun),out)
c                                  - foundation information (integer)
c          ififo(1,ifoun)       = number of segments
c          ififo(4,ifoun)       = constrain condition in x-dir.
c          ififo(5,ifoun)       = constrain condition in y-dir.
c          ififo(6,ifoun)       = constrain condition in z-dir.
c          ififo(7,ifoun)       = dimensions
c          ififo(8,ifoun)       = 1st segment numbering on this foundation
c          ififo(9,ifoun)       = 2nd segment numbering on this foundation
c          ififo(10,ifoun)      = lowest nodal numbering on this foundation
c          ififo(12,ifoun)      = real constant number (ireal)
c          ififo(13,ifoun)      = contact pair key
c                                 0 - no contact pair
c                                 1 - only slave element is defined
c                                 2 - only contact element is defined
c                                 3 - contact pair is defined
c          ififo(14,ifoun)      = type of associated contant element type (ityp)
c          ififo(15,ifoun)      = type of deformable slave surface
c                                -2 - mixed type surface (may wrong)
c                                -1 - rigid slave surface
c                                 0 - none
c                                 1 - deformable slave surface
c                                 2 - deformable slave surface (two pass)
c          ififo(16,ifoun)      = contact element number who is closest to current slave
c          ififo(17,ifoun)      = slave element number who is pair of closest contact
c                                 element
c          ififo(18,ifoun)      = iteration number of surface adjust
c     founda   (dp,ar(nperf,nFoun),out)
c                                  - foundation  cell region
c          founda(1,ifoun)       = minimum value in x-dir.
c          founda(2,ifoun)       = maximum value in x-dir.
c          founda(3,ifoun)       = minimum value in y-dir.
c          founda(4,ifoun)       = maximum value in y-dir.
c          founda(5,ifoun)       = minimum value in z-dir.
c          founda(6,ifoun)       = maximum value in z-dir.
c          founda(7,ifoun)       = first real constant
c          founda(8,ifoun)       = second real constant
c          founda(9,ifoun)       = gap or penetration of closest contact element
c          founda(10,ifoun)      = x component of founda(9,ifoun) (increment)
c          founda(11,ifoun)      = y component of founda(9,ifoun) (increment)
c          founda(12,ifoun)      = z component of founda(9,ifoun) (increment)
c          founda(13,ifoun)      = x component of founda(9,ifoun) (total)
c          founda(14,ifoun)      = y component of founda(9,ifoun) (total)
c          founda(15,ifoun)      = z component of founda(9,ifoun) (total)
c          founda(16,ifoun)      = length of contact surface
c          founda(17,ifoun)      = depth of contact surface
c          founda(18,ifoun)      = stiffness of founda(9,ifoun)
c          founda(19,ifoun)      = number of founda(9,ifoun)
c          founda(20,ifoun)      = area of founda(9,ifoun)
c     icnda    (int,ar(nPerMast,nMast),out) - contact element information
c                               = 1 element type (ktyp)
c                               = 2 attached beam element number (ibeam)
c                               = 3 integration type (intyp)
c                               = 4 - 7 nodes number (nodes)
c                               = 8 contact element number (ielm)
c                               = 9 slave surface number (ifoun)
c
#include "impcom.inc"
#include "constants.inc"
#include "conv_stat.inc"
#include "solu_info.inc"
c
      PTRFTN   MemAlloc
      integer BacLen,NumDof,nPerSeg,nSeg,nPerFounI,
     &  nPerFounDp,nFoun,elmget,NumEl,nPerMast,nMast,
     &  etyiqr,rlget,elmgct,GetWriteUnit
      integer nodfwd(BacLen),isegda(nPerSeg,nSeg), nContClo, nContIni,
     &  ififo(nPerFounI,nFoun),icnda(nPerMast,nMast),iott,
     &  nodbac(BacLen),inoexc,iroexc,inoexc1,iroexc1, 
     &  closi,methi,keyopt_1
      double precision founda(nPerFounDp,nFoun),icoor(BacLen*NumDof),
     &  dispn(BacLen*NumDof)
      integer elData(EL_DIM),nodes(W_NdMax),ival(3),ifoun,i,ielm,nnod,
     &  inode,ietypz,icon,ncon,ireal,nrvr,nnode,itetry,ityp,
     &  go50key,iseg,kfoun
      double precision rvr(W_NumContParam),cradi,dincre,pmax,pmin,
     &  dfact,ZERO,length,igtol,xyzc(7),dval(4),gapcri,da,dperr(1)

      parameter (ZERO=0.0d0)

      pointer (piDist,Dist)      
      double precision  Dist(*)


#if defined(DEBUG)
      call SubBeg('Cont_5')
#endif

      itetry = 20
      C_nPen = 0
      do ifoun = 1,nFoun
         ififo(18,ifoun) = 0
         ififo(22,ifoun) = 0
         founda(9,ifoun) = HUGE
         founda(10,ifoun) = ZERO
         founda(16,ifoun) = ZERO
         founda(17,ifoun) = ZERO
         founda(18,ifoun) = ZERO
         founda(19,ifoun) = ZERO
         founda(20,ifoun) = ZERO
         founda(21,ifoun) = ZERO
         founda(26,ifoun) = ZERO
      enddo
  50  continue

      piDist = MemAlloc(ncon,W_LenDouble,'piDist          ')
            
      do 100 icon = 1,ncon
         ielm = icnda(8,icon)
         nnod = elmget(ielm,elData(1),nodes(1))
         if (nnod .le. 0) goto 100
         ityp = elData(EL_TYPE)
         ietypz = etyiqr(ityp,-JETYP)
         if (ietypz.ne.W_eMaster) goto 100
         call conget (ielm,ival(1))
         ifoun = abs(ival(3))
         if (ifoun.eq.0) goto 100
         if (ififo(18,ifoun).lt.0) goto 100
         ireal = elData(EL_REAL)
         call v_zero (rvr(1),W_NumContParam)
         nrvr = rlget(ireal,rvr(1))
         if (ififo(15,ifoun).eq.-1.or.ififo(15,ifoun).le.-3) then
c --- rigid-deformable contact
             nnode = 3
         elseif (ififo(15,ifoun).gt.0) then
c --- deformable-deformable contact
             nnode = 7
         else
             dperr(1) = dble(ireal)
             call erhandler ('Cont_5', 5000, 3, 'The real constant set
     &                        %I is refered but undefined',
     &                        dperr(1),' ')
             call StopAll (-5)
         endif             

         do 70 inode = 1,nnode
            if (nodes(inode).gt.0)
     &          nodes(inode) = (nodfwd(nodes(inode))-1)*NumDof + 1
  70     continue
         
         if (etyiqr(ififo(14,ifoun),-KYO10).eq.0) then
             call CalIniPeni (ielm,nnode,ityp,rvr(1),ival(1),nSeg,
     &             nPerSeg,nPerFounI,nPerFounDp,nPerMast,nMast,
     &             nodes(1),ififo(1,ifoun),icoor(1),dispn(1),
     &             isegda(1,1),founda(1,ifoun),icnda(1,1),icon,
     &             Dist(icon))
         else !if  (elChar(KYO10).eq.1 .or. elChar(KYO10).eq.2) then
             call CalIniPeni4 (ielm,ityp,rvr(1),ival(1),nSeg,
     &             ifoun,nPerSeg,nPerFounI,nPerFounDp,nPerMast,nMast,
     &             nodes(1),ififo(1,ifoun),icoor(1),dispn(1),
     &             isegda(1,1),founda(1,ifoun),icnda(1,1),icon,
     &             Dist(icon))
         endif
  100 continue


      do 500 ifoun = 1,nFoun
         if (ififo(18,ifoun).lt.0) goto 500
         if (founda(19,ifoun).ge.1.0d0) then
             founda(16,ifoun) = founda(16,ifoun)/founda(19,ifoun)
             founda(17,ifoun) = founda(17,ifoun)/founda(19,ifoun)
             founda(18,ifoun) = founda(18,ifoun)/founda(19,ifoun)
             founda(20,ifoun) = founda(20,ifoun)/founda(19,ifoun)
         endif
         if (ififo(1,ifoun).le.0) then
             ififo(18,ifoun) = -2
             goto 500
         endif
c
         if (ififo(4,ifoun).ne.0) ififo(4,ifoun) = 1
         if (ififo(5,ifoun).ne.0) ififo(5,ifoun) = 1
         if (ififo(6,ifoun).ne.0) ififo(6,ifoun) = 1
         if (ififo(7,ifoun).eq.2) ififo(6,ifoun) = 0
         if (ififo(4,ifoun)+ififo(5,ifoun)+ififo(6,ifoun).eq.0)
     &       ififo(18,ifoun) = -2
         ireal = ififo(12,ifoun)
         call v_zero (rvr(1),W_NumContParam)
         nrvr = rlget(ireal,rvr(1))
         pmax = rvr(7)
         pmin = rvr(8)
         cradi = founda(17,ifoun)
         length = founda(16,ifoun)
         if (ififo(18,ifoun).eq.0) founda(26,ifoun) = founda(9,ifoun)
         if (etyiqr(ififo(14,ifoun),-KYOP3).ge.3) then
             cradi = founda(25,ifoun)
             length = founda(27,ifoun)
         endif
         if (abs(pmin).le.TINY .and. abs(pmax).le.TINY) then
             pmax = 1.d6*cradi
             ififo(18,ifoun) = -2
         else
            if (abs(pmin).le.TINY) then
                pmin = 0.005d0*cradi
            else
               if (pmin.lt.ZERO) then
                   pmin = -pmin
               else
                   pmin = pmin*cradi
               endif
            endif
            if (abs(pmax).le.TINY) then
                if (rvr(4).gt.ZERO) then
                    pmax = 0.5d0*rvr(4)*cradi
                else
                    pmax = 0.5d0*rvr(4)
                endif
            else
                if (pmax.gt.ZERO) pmax = pmax*cradi
            endif
            if (abs(pmax).le.TINY) pmax = 0.05d0*cradi
            if (pmax.lt.ZERO) pmax = -pmax
            pmax = min(pmax,0.5d0*length,0.5d0*cradi)
            if (abs(rvr(7)).le.TINY) then
                pmax = max(2.0d0*pmin,pmax)
            else
                pmax = max(1.1111d0*pmin,pmax)
            endif
            if (abs(rvr(8)).le.TINY) then
                pmin = min(pmin,0.5d0*pmax)
            else
                pmin = min(pmin,0.9d0*pmax)
            endif
         endif
         if (ififo(16,ifoun).eq.0) then
            ififo(18,ifoun) = -1
         else
            dfact = max(abs(founda(9,ifoun)),1.d-6*cradi)
            if (founda(9,ifoun).le.ZERO) then
               if (abs(founda(9,ifoun)).ge.pmax) then
                  itetry = 50
                  ififo(18,ifoun) = ififo(18,ifoun)+1
                  dfact = 0.99d0*pmax-dfact
               else
                  itetry = 100
                  if (abs(founda(9,ifoun)).ge.pmin.or.
     &                    ififo(18,ifoun).gt.itetry) then
                     if (ififo(18,ifoun).eq.-2) goto 350
                     ififo(18,ifoun) = -3
                     do 300 i = 0,ififo(7,ifoun)-1
                        founda(2*i+1,ifoun) = founda(2*i+1,ifoun) +
     &                                        founda(13+i,ifoun)
                        founda(2*i+2,ifoun) = founda(2*i+2,ifoun) +
     &                                        founda(13+i,ifoun)
  300                continue

                     iseg = isegda(15,ififo(8,ifoun))
                     nrvr = elmgct(iseg,xyzc(1))
                     xyzc(5) = founda(13,ifoun)
                     xyzc(6) = founda(14,ifoun)
                     call elmpct (iseg,xyzc(1))

  350                continue
                     if (abs(rvr(5)).le.TINY.and.
     &                   (etyiqr(ififo(14,ifoun),-KYOP5).eq.0.or.
     &                   (etyiqr(ififo(14,ifoun),-KYOP5).eq.4.and.
     &                    etyiqr(ififo(14,ifoun),-KYOP1).eq.6))) then
                        igtol = abs(founda(9,ifoun))/cradi
                        igtol = min(igtol,0.01d0)
                        igtol = max(igtol,0.001d0)
                        if (etyiqr(ififo(14,ifoun),-KYOP1).eq.6)
     &                      igtol = max(igtol,0.05d0)
                        founda(22,ifoun) = igtol
                     elseif (etyiqr(ififo(14,ifoun),-KYOP5).eq.2.or.
     &                       etyiqr(ififo(14,ifoun),-KYOP5).eq.3) then
                        igtol = 0.99d0*founda(9,ifoun)
                        founda(22,ifoun) = igtol
                     endif
                     goto 450
                  else
                     ififo(18,ifoun) = ififo(18,ifoun)+1
                     dfact = 1.01d0*pmin-dfact
                  endif
               endif
            else
               if (ififo(18,ifoun).ne.-2) goto 360

               if (abs(rvr(5)).le.TINY.and.
     &             (etyiqr(ififo(14,ifoun),-KYOP5).eq.0.or.
     &             (etyiqr(ififo(14,ifoun),-KYOP5).eq.4.and.
     &              etyiqr(ififo(14,ifoun),-KYOP1).eq.6))) then

                   igtol = 2.0d0*abs(founda(9,ifoun))/cradi
                   igtol = min(0.03d0,igtol)
                   igtol = max(0.01d0,igtol)
                   if (etyiqr(ififo(14,ifoun),-KYOP1).eq.6)
     &                 igtol = max(igtol,0.05d0)
                   founda(22,ifoun) = igtol
               elseif (etyiqr(ififo(14,ifoun),-KYOP5).eq.1.or.
     &                 etyiqr(ififo(14,ifoun),-KYOP5).eq.3) then
                   igtol = 1.01d0*abs(founda(9,ifoun))
                   founda(22,ifoun) = igtol
               endif
 360           continue
               dfact = dfact+pmin
               ififo(18,ifoun) = ififo(18,ifoun) + 1
            endif
            if (ififo(15,ifoun).ge.2) ififo(18,ifoun) = -1
            if (ififo(18,ifoun).lt.0) goto 450
            if (ififo(18,ifoun).ge.itetry) then
                founda(13,ifoun) = ZERO
                founda(14,ifoun) = ZERO
                founda(15,ifoun) = ZERO
                if (ififo(18,ifoun).eq.itetry) then
                    ififo(18,ifoun) = 1+ififo(18,ifoun)
                else
                    ififo(18,ifoun) = -ififo(18,ifoun)
                endif
                goto 450
            endif
            if (ififo(18,ifoun).le.1) founda(20,ifoun) = founda(9,ifoun)
            founda(9,ifoun) = HUGE
            do 400 i = 0,ififo(7,ifoun)-1
               dincre = 1.01d0*dfact*founda(10+i,ifoun)*ififo(4+i,ifoun)
               founda(13+i,ifoun) = founda(13+i,ifoun) + dincre
  400       continue
  450       continue
         endif
  500 continue

      go50key = 0
      do 550 ifoun = 1,nFoun
         if (ififo(18,ifoun).ge.0) then
             founda(16,ifoun) = ZERO
             founda(17,ifoun) = ZERO
             founda(18,ifoun) = ZERO
             founda(19,ifoun) = ZERO
             founda(20,ifoun) = ZERO
             founda(21,ifoun) = ZERO
             go50key = 1
         endif
  550 continue
      if (go50key.eq.1) goto 50

      
      nContIni = 0
      nContClo = 0

c --- define the status for individual contact elements
      do 200 icon = 1,ncon
         if (icnda(3,icon).eq.0) goto 200
         ielm = icnda(8,icon)
         nnod = elmget(ielm,elData(1),nodes(1))
         if (nnod .le. 0) goto 200
         ityp = elData(EL_TYPE)
         ietypz = etyiqr(ityp,-JETYP)
         if (ietypz.ne.W_eMaster) goto 200
         call conget (ielm,ival(1))
         ifoun = abs(ival(3))
         if (ifoun.eq.0) goto 200
         ireal = elData(EL_REAL)
         call v_zero (rvr(1),W_NumContParam)
         nrvr = rlget(ireal,rvr(1))
         
         keyopt_1 = etyiqr(ityp,-KYOP1)
         closi = etyiqr(ityp,-KYOP5)
         methi = etyiqr(ityp,-KYOP8)

c --- select critical gap (copy of the elements' coding)
         gapcri = rvr(5)
         if (abs(gapcri).le.TINY.and.closi.eq.0) 
     &       gapcri = founda(22,ifoun)

         if (abs(gapcri).le.TINY.and.closi.eq.4.and.keyopt_1.eq.6) 
     &       gapcri = founda(22,ifoun)

         if (gapcri.gt.ZERO) then
             gapcri = min(gapcri*founda(16,ifoun),
     &                    gapcri*founda(17,ifoun))
            if (gapcri.ge.0.9d0*founda(17,ifoun)) then
                gapcri = min(0.05d0*founda(16,ifoun),
     &                       0.05d0*founda(17,ifoun))
            endif
         else
            gapcri = - gapcri
         endif
         da = Dist(icon)
c --- check if the gap/penetration will be closed
         if ((da.gt.ZERO.and.da.le.gapcri).or.
     &      ((da.lt.ZERO.or.keyopt_1.eq.4.or.keyopt_1.eq.5)
     &        .and.(methi.eq.3.or.methi.eq.4.or.methi.eq.0))) then
c --- closed, nothing to do
         elseif (da.ne.ZERO) then
c --- cannot be closed
              icnda(3,icon) = 0
         endif

  200 continue
      call MemFree (piDist)


c --- define the status for the entire foundation and count contacts
      do 1001 kfoun = 1,nFoun
         ifoun = ififo(24,kfoun)
         ireal = ififo(12,ifoun)
c --- check if self-contact
         if (ififo(19,ifoun).eq.ireal.and.ififo(15,ifoun).eq.2) then
             ififo(15,ifoun) = 3
             founda(18,ifoun) = 0.5d0*founda(18,ifoun)
             founda(24,ifoun) = 0.5d0*founda(24,ifoun)
             founda(20,ifoun) = 0.5d0*founda(20,ifoun)
         endif

         da = founda(9,ifoun)
         ififo(26,ifoun) = 0
c --- if there is no contact at all, fast
         if (ififo(16,ifoun).eq.0) then
             if (abs(da).gt.founda(20,ifoun)) then
                 ififo(26,ifoun) = -2
             else
                 ififo(26,ifoun) = -1
             endif
             goto 1001
         endif

         ireal = ififo(12,ifoun)
         ityp = ififo(14,ifoun)
         keyopt_1 = etyiqr(ityp,-KYOP1)
         closi = etyiqr(ityp,-KYOP5)
         methi = etyiqr(ityp,-KYOP8)
         call v_zero (rvr(1),W_NumContParam)
         nrvr = rlget(ireal,rvr(1))

c --- select critical gap (copy of the elements' coding)
         gapcri = rvr(5)
        
         if (abs(gapcri).le.TINY.and.closi.eq.0) 
     &       gapcri = founda(22,ifoun)

         if (abs(gapcri).le.TINY.and.closi.eq.4.and.keyopt_1.eq.6) 
     &       gapcri = founda(22,ifoun)


         if (gapcri.gt.ZERO) then
             gapcri = min(gapcri*founda(16,ifoun),
     &                    gapcri*founda(17,ifoun))
            if (gapcri.ge.0.9d0*founda(17,ifoun)) then
                gapcri = min(0.05d0*founda(16,ifoun),
     &                       0.05d0*founda(17,ifoun))
            endif
         else
            gapcri = - gapcri
         endif

         
c --- check if the gap/penetration will be closed
         if ((da.gt.ZERO.and.da.le.gapcri).or.
     &      ((da.lt.ZERO.or.keyopt_1.eq.4.or.keyopt_1.eq.5)
     &        .and.(methi.eq.3.or.methi.eq.4.or.methi.eq.0))) then
c --- closed, nothing to do
            nContClo = nContClo + 1
         elseif (da.eq.ZERO) then
            nContIni = nContIni + 1
         else
c --- cannot be closed
            if (abs(da).gt.founda(20,ifoun)) then
                ififo(26,ifoun) = -2
            else
                ififo(26,ifoun) = -1
            endif
         endif

 1001 continue

  610 iott = GetWriteUnit()
      write (iott,2000) nContIni, nContClo, nContIni + nContClo
 2000 format(/'---> NUMBER OF CONTACTS:   INITIAL=',i7,
     &         '   CLOSED=',i7,'   TOTAL=',i7/)

c --- transfer nodal constraints and loading on rot DOF to the second pilot node
      do 700 kfoun = 1,nFoun
         ifoun = ififo(24,kfoun)
         if (ififo(1,ifoun).le.0) goto 700
         if (ififo(13,ifoun).le.0) goto 700
         inoexc1 = (ififo(2,ifoun)-1)/NumDof + 1
         ififo(3,ifoun) = abs(ififo(3,ifoun))
         if (ififo(15,ifoun).ge.0) goto 700
         if (inoexc1.gt.0) then
            inoexc = nodbac(inoexc1)
            if (ififo(3,ifoun).eq.0) then
               iroexc1 = 0
            else
               iroexc1 = (ififo(3,ifoun)-1)/NumDof + 1
            endif
            if (iroexc1.gt.0) then
                iroexc = nodbac(iroexc1)
                call v_zero (dval(1),4)

                do i = 1,3
                   call disput (iroexc,3+i,dval(1))
                enddo

                do i = 1,3
                   call disget (inoexc,3+i,dval(1))
                   call disput (iroexc,i,dval(1))
                   call forget (inoexc,3+i,dval(1))
                   call forput (iroexc,i,dval(1))
                enddo
            endif
         endif
  700 continue 

  999 continue

#if defined(DEBUG)
      call SubEnd('Cont_5')
#endif
      return
      end