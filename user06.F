      function user06()
c --- primary functions: 1) setup the contact
c                        2) report contact information
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer user06, ndgxyz, etyiqr,
     &        elmget, elmiqr, ndinqr, rlinqr, etyget
      integer kerr, i, j, i1, j1, k, node, ncount,
     &        elem, ecount, elData(EL_DIM),
     &        elChar(IELCSZ), ireal, rcount, ET1, ET2,
     &        intinfun,jreal,kreal, iFib, jFib,
     &        cRl, nRl, ni, nj, debug, iott, GetWriteUnit,
     &        ver,progress,cont,etyp,pareType

      double precision rvr(6),CalDist2P,length1,length2,
     &                 Ex1,Ex2,seed,dpinfun,fract,rnd,
     &                 xa(3),ya(3),za(3),xb(3),yb(3),zb(3),
     &                 xi(3),yi(3),zi(3),xj(3),yj(3),zj(3),
     &                 tol,epsil,tols,GetFibRadQ,minD,angtols,
     &                 radSF      

      character*28  elcdn
      logical ElFibOverlap,ElFibOverlapXY  
      PTRFTN  MemAlloc, MemRealloc, MemAllocL
      pointer (PtrElNd,ElNd)
      pointer (PtrNdXYZ,NdXYZ)
      pointer (PtrMxMn,MxMn)
      pointer (PtrRl,Rl)
      pointer (PtreReal,eReal)
      pointer (PtreR,eR)
      pointer (PtreEx,eEx)
      pointer (PtriCont,iCont)
      pointer (PtrjCont,jCont)
      pointer (PtreLength,eLength)
      pointer (Ptrndpel,ndpel)
      integer*8 LongL, iL, jL, ishift
      
      byte Rl(*)
      
      integer  ElNd(4,*),  eReal(*), iCont(*), jCont(*), ndpel(*)
      double precision  NdXYZ(3,*), MxMn(3,2,*), eEx(*), eR(*), 
     &                  eLength(*)

      integer ioerr, nott, Sys_Open, Sys_Close
      character*(B_NameLen)  jobnam
      character*(B_DirLen)   fdir  
      character*(B_PathLen)  pname
      
      ecount = elmiqr(0,DB_MAXDEFINED)

      ver = intinfun(2)
      iott = GetWriteUnit()
      write (iott,1000)
 1000 format (/'-> CONTACT SETUP STARTED...')

      tols = 1.0d-6
      if (ver.eq.1) then
          write (iott,*) 'Aggressive trimming was requested'
      elseif (ver.eq.2) then
          tols = 2.0d0
          write (iott,*) 'Moderate trimming was requested'
      elseif (ver.eq.3) then
          tols =  2.0d0
          write (iott,*) 'Reporting the contact info only'
      endif
       
      cont = intinfun(3)
      fract = min(abs(dpinfun(4)),1.0d0)
      seed = 0.0d0
      
      pareType = intinfun(5)
      angtols = 0.0d0
      if (pareType.gt.3) then
          angtols = 0.02d0
          write (iott,*) 'Will include parallel beam contact'
      else
          pareType = 0
      endif
            
      if (cont.eq.1) then
          write (iott,*) 'Will close large initial gap/penetration'
      endif

      if (fract.gt.TINY) then
          write (iott,2000) fract
 2000     format ('Fraction of active bonds is: ',g11.4)          
      endif

      radSF = dpinfun(6)
      if (radSF.gt.TINY) then
          write (iott,2001) radSF
 2001     format ('Scale the radius by ',g11.4)          
      else
          radSF = 1.0d0
      endif

      PtrElNd  = MemAlloc(ecount*4,W_LenInt,'PtrElNd         ')
      PtreReal = MemAlloc(ecount,W_LenInt,'PtreReal        ')
      PtreR = MemAlloc(ecount,W_LenDouble,'PtreR           ')
      PtreEx = MemAlloc(ecount,W_LenDouble,'PtreEx          ')
      Ptrndpel = MemAlloc(ecount,W_LenInt,'Ptrndpel        ')
      
      call vi_zero (eReal(1),ecount)
      call vi_zero (ElNd(1,1),ecount*4)
      call v_zero (eEx(1),ecount)
      call v_zero (eR(1),ecount)
      call vi_zero (ndpel(1),ecount)

c --- get the beam elements only
      rcount = 0
      ncount = 0
      do elem = 1,ecount
         i = elmget (elem,elData(1),ElNd(1,elem))
         if (ElNd(1,elem).le.0.or.i.le.0) then
             ecount = elem - 1
             goto 10
         endif
c --- store number of nodes depending on the element type
         etyp = etyiqr(elData(EL_TYPE),-JETYP)
         if (etyp.eq.W_eBeam3nd) then
             ndpel(elem) = 3
         elseif (etyp.eq.W_eBeam2nd) then
             ndpel(elem) = 2
         elseif (etyp.eq.W_eBeamEB) then
             ndpel(elem) = 2
         else
             ecount = elem - 1
             goto 10
         endif
c --- store real numbers
         eReal(elem) = elData(EL_REAL)
         rcount = max (rcount,eReal(elem))
         ncount = max (ncount,ElNd(1,elem),ElNd(2,elem),ElNd(3,elem))
         eR(elem) = GetFibRadQ(elem)*radSF
         call prope1 (elem,elData(EL_MAT),1,0.0d0,eEx(elem))
      enddo

  10  continue

      if (rcount.eq.0 .or. ncount.eq.0) then
          goto 997
      endif
        
      PtrElNd = MemRealloc(PtrElNd,ecount*4,W_LenInt,'PtrElNd         ')
      PtreReal = MemRealloc(PtreReal,ecount,W_LenInt,'PtreReal        ')
      PtreR = MemRealloc(PtreR,ecount,W_LenDouble,'PtreR           ')
      PtreEx = MemRealloc(PtreEx,ecount,W_LenDouble,'PtreEx          ')
      Ptrndpel = MemRealloc(Ptrndpel,ecount,W_LenInt,'Ptrndpel        ')
      PtrNdXYZ = MemAlloc(ncount*3,W_LenDouble,'PtrNdXYZ        ')
      LongL = rcount*rcount
      PtrRl = MemAllocL(LongL,W_LenByte,'PtrRl           ')
      PtrMxMn = MemAlloc(3*2*rcount,W_LenDouble,'PtrMxMn         ')

      PtreLength = MemAlloc(ecount,W_LenDouble,'PtreLength      ')
      PtriCont = MemAlloc(ecount,W_LenInt,'PtriCont        ')
      PtrjCont = MemAlloc(ecount,W_LenInt,'PtrjCont        ')

      call v_zero (NdXYZ(1,1),ncount*3)      
      call v_zero (eLength(1),ecount)
      do iL = 1,LongL
         Rl(iL) = 0
      enddo
      
      
c --- get nodes
      do node = 1,ncount
         i = ndgxyz (node,NdXYZ(1,node))
      enddo
c --- first, fill some values from any node in the set (need for comparision)
      do elem = 1,ecount
         ireal = eReal(elem)
         if (Rl(ireal).eq.1) cycle
         Rl(ireal) = 1
         do k = 1,3
            MxMn(k,1,ireal) = NdXYZ(k,ElNd(1,elem))
         enddo
         call v_move (MxMn(1,1,ireal),MxMn(1,2,ireal),3)
      enddo

      do i = 1,rcount
         Rl(i) = 0
      enddo
      
c --- to prevent self contact
      do jL = 1,rcount
         iL = rcount*(jL-1) + jL
         Rl(iL) = 1
      enddo

c --- fill-in boxes comparing with the existing values
      do elem = 1,ecount
         ireal = eReal(elem)
c --- loop thru coordinates
         do k = 1,3
c --- loop thru the nodes
            do j = 1,ndpel(elem)
               MxMn(k,1,ireal) = min(MxMn(k,1,ireal),
     &                               NdXYZ(k,ElNd(j,elem)))
               MxMn(k,2,ireal) = max(MxMn(k,2,ireal),
     &                               NdXYZ(k,ElNd(j,elem)))
            enddo
         enddo
         eLength(elem) = CalDist2P(NdXYZ(1,ElNd(1,elem)),
     &                             NdXYZ(1,ElNd(2,elem)))
      enddo
      
      if (ver.eq.3) goto 300

      i = etyget(2,elChar(1))
      if (i.eq.0) call vi_zero (elChar(1),IELCSZ)
      call uec101 (elcdn,elChar(1),kerr)
      call etyput (2,IELCSZ,elChar(1))

      i = etyget(3,elChar(1))
      if (i.eq.0) call vi_zero (elChar(1),IELCSZ)
c --- mark crossing beam formalation      
      elChar(KYO10) = 0
      call uec102 (elcdn,elChar(1),kerr)
      call etyput (3,IELCSZ,elChar(1))

      call vi_zero (elData(1),EL_DIM)
      elData(EL_SHAPE) = 2
      elData(EL_DIM) = 3
      elData(EL_MAT) = 0

      elem = ecount

      call v_zero(rvr(1),6)
      nRl = rcount
      progress = 10
      
      do 100 i = 1,ecount
c --- get radius for the first fiber
         rvr(1) = eR(i)
c --- get Ex for the first fiber
         ireal = eReal(i)
         ishift = (ireal-1)*nRl
         Ex1 = eEx(i)

         call MoveCoords (xa(1),ya(1),za(1),
     &                    NdXYZ(1,1),ElNd(1,i),ndpel(i))
         
         length1 = eLength(i)

         ni = 0
c --- skip checking if only one contact
         if (ver.ne.1) then
c --- form the entire fibers
             do k = 1,ecount
                kreal = eReal(k)
                if (kreal.eq.ireal) then
                    ni = ni + 1
                    iCont(ni) = k
                endif
            enddo
         endif

         if ((i*100/ecount).ge.progress) then
             write (iott,*) 'Done (crossing)',progress,'%'
             progress = progress + 10
         endif

         do 200 j = i+1, ecount
            jreal = eReal(j)
c --- do not consider self contact
            cRl = Rl(ishift+jreal)
            
            if (cRl.ne.0 .and. cRl.ne.2) goto 200
            
c --- get radius for the second fiber
            rvr(2) = eR(j)
            tol = 4.0d0*max(rvr(1),rvr(2))

c --- find out if fibers overlap (super fast)
            
            if (cRl.ne.2) then
                if (.not.ElFibOverlap(MxMn(1,1,1),ireal,jreal,tol)) then
                    Rl(ishift+jreal) = 1
                    goto 200
                else
c --- mark than no further check is done
                    Rl(ishift+jreal) = 2
                endif
            endif

            length2 = eLength(j)

            if (min(length1,length2).lt.TINY) then
c --- capture zero length here
                Rl(ishift+jreal) = 1
                goto 200
            elseif 
     &         (max(length1,length2)/min(length1,length2).gt.10.d0) then
c --- skip if the aspect ratio is large
                Rl(ishift+jreal) = 1
                goto 200
            endif            

            call MoveCoords (xb(1),yb(1),zb(1),
     &                       NdXYZ(1,1),ElNd(1,j),ndpel(j))

c --- find out if elements overlap (super fast)
            debug = 0

c            if (ireal.eq.1 .and. jreal.eq.22) then
c                write (8,*) 'i,j',i,j
c                debug = 1
c            endif
            
            call BeamToBeamFast (rvr(1),xa(1),ya(1),za(1),
     &                           xb(1),yb(1),zb(1),tol,epsil,
     &                           debug,cont,tols,angtols)

            if (epsil.eq.-HUGE) then
                Rl(ishift+jreal) = 3
                goto 200
            endif

            if (abs(epsil).lt.tol) then
                if (fract.gt.0.0d0 .and. fract.lt.1.0d0) then
                    call v_rand(rnd,1,seed)
                    if (rnd.gt.fract) then
                        Rl(ishift+jreal) = 1
                        goto 200
                    endif
                endif

                nj = 0
c --- skip checking if only one contact
                if (ver.ne.1) then
c --- form the entire fibers
                    do k = 1,ecount
                       kreal = eReal(k)
                       if (kreal.eq.jreal) then
                           nj = nj + 1
                           jCont(nj) = k
                       endif
                    enddo
c --- in ver=2 we only retain contacts which are initially closed
                    if (ver.eq.2) then
                        do i1 = 1,ni
c --- iCont is updated once per cycle and can be changed within, thus we need abs()
                           iCont(i1) = - abs(iCont(i1))
                        enddo
                        do j1 = 1,nj
                           jCont(j1) = - jCont(j1)
                        enddo
                    endif
c --- go through all the segment and insure there is no problems                
                    do i1 = 1,ni
                       iFib = abs(iCont(i1))

                       call MoveCoords (xi(1),yi(1),zi(1),
     &                                  NdXYZ(1,1),ElNd(1,iFib),
     &                                  ndpel(iFib))

                       do j1 = 1,nj
                          jFib = abs(jCont(j1))
                          call MoveCoords (xj(1),yj(1),zj(1),
     &                                     NdXYZ(1,1),ElNd(1,jFib),
     &                                     ndpel(jFib))

                          debug = 0
                          call BeamToBeamFast (rvr(1),xi(1),yi(1),zi(1),
     &                                         xj(1),yj(1),zj(1),tol,
     &                                         epsil,debug,cont,
     &                                         tols,angtols)
     
                          if (epsil.eq.-HUGE) then
                              Rl(ishift+jreal) = 3
                              goto 200
                          elseif (ver.eq.2.and.abs(epsil).lt.tol) then
c --- indicate that these contacts will be formed
                              iCont(i1) = abs(iCont(i1))
                              jCont(j1) = abs(jCont(j1))
                              if (i1.gt.1) iCont(i1-1)=abs(iCont(i1-1))
                              if (i1.lt.ni) iCont(i1+1)=abs(iCont(i1+1))
                              if (j1.gt.1) jCont(j1-1)=abs(jCont(j1-1))
                              if (j1.lt.nj) jCont(j1+1)=abs(jCont(j1+1))
                          endif
                       enddo
                    enddo
                endif
c --- if this is the first contact for this pair of beams
                rcount = rcount + 1

                cRl = rcount
c --- get Ex the for second fiber
                Ex2 = eEx(j)
c --- store real data for this foundation
                call rlput (rcount,5,rvr(1))
c --- decide which one should be slave (stiffer, thicker)
                if (ni.lt.nj) then
                    ET1 = 3
                    ET2 = 2
                elseif (ni.gt.nj) then
                    ET1 = 2
                    ET2 = 3
                else
                    if (length1.lt.length2) then
                        ET1 = 3
                        ET2 = 2
                    else
                        ET1 = 2
                        ET2 = 3
                    endif
                endif

                ET1 = 3
                ET2 = 2

c --- the desision about the contact pair is defined in the sign of cRl
                Rl(ishift+jreal) = 4 ! mark as having contact
                elData(EL_REAL) = cRl
                
                if (ver.eq.1) then
c --- one contact per pair (for bonded case)
                    elem = elem + 1
                    elData(EL_TYPE) = ET1
                    call elmput (elem,elData(1),ndpel(i),ElNd(1,i))
                    elem = elem + 1
                    elData(EL_TYPE) = ET2
                    call elmput (elem,elData(1),ndpel(j),ElNd(1,j))
                else
c --- full contact per pair (for slinding case)
                    elData(EL_TYPE) = ET1
                    do i1 = 1,ni
                       k = iCont(i1)
                       if (k.le.0) cycle
c --- store first series of segments
                       elem = elem + 1
                       call elmput (elem,elData(1),ndpel(k),ElNd(1,k))
                    enddo

                    elData(EL_TYPE) = ET2
                    do j1 = 1,nj
                       k = jCont(j1)
                       if (k.le.0) cycle
c --- store second series of segments
                       elem = elem + 1
                       call elmput (elem,elData(1),ndpel(k),ElNd(1,k))
                    enddo
                endif
            endif

 200     continue
 100  continue

      
c --- setup the parallel beam contacts      
      if (pareType.gt.0) then
          i = etyget(pareType,elChar(1))
          if (i.eq.0) call vi_zero (elChar(1),IELCSZ)
          elChar(KYO10) = 1
          call uec102 (elcdn,elChar(1),kerr)      
          call etyput (pareType,IELCSZ,elChar(1))          
      else
          goto 996
      endif

      call v_zero(rvr(1),6)
      progress = 10
      do 400 i = 1,ecount
c --- get radius for the first fiber
         rvr(1) = eR(i)
c --- get Ex for the first fiber
         ireal = eReal(i)
         ishift = (ireal-1)*nRl
         Ex1 = eEx(i)

         call MoveCoords (xa(1),ya(1),za(1),
     &                    NdXYZ(1,1),ElNd(1,i),ndpel(i))
         
         length1 = eLength(i)

         ni = 0
c --- skip checking if only one contact
c --- form the entire fibers
         do k = 1,ecount
            kreal = eReal(k)
            if (kreal.eq.ireal) then
                ni = ni + 1
                iCont(ni) = k
            endif
         enddo

         if ((i*100/ecount).ge.progress) then
             write (iott,*) 'Done (parallel)',progress,'%'
             progress = progress + 10
         endif

         do 500 j = i+1, ecount
            jreal = eReal(j)
c --- do not consider self contact
            cRl = Rl(ishift+jreal)
            
            if (cRl.eq.4 .or. cRl.eq.1) goto 500
c --- get radius for the second fiber
            rvr(2) = eR(j)
            tol = 2.0d0*max(rvr(1),rvr(2))

            length2 = eLength(j)

            if (min(length1,length2).lt.TINY) then
c --- capture zero length here
                Rl(ishift+jreal) = 1
                goto 500
            elseif 
     &         (max(length1,length2)/min(length1,length2).gt.10.d0) then
c --- skip if the aspect ratio is large
                Rl(ishift+jreal) = 1
                goto 500
            endif            

            call MoveCoords (xb(1),yb(1),zb(1),
     &                       NdXYZ(1,1),ElNd(1,j),ndpel(j))

c --- find out if elements overlap (super fast)
            debug = 0

            call ParBeam2ndFast (rvr(1),
     &                           xa(1),ya(1),za(1),
     &                           xb(1),yb(1),zb(1),
     &                           tol,epsil,cont,tols)

            if (abs(epsil).lt.tol) then
                if (fract.gt.0.0d0 .and. fract.lt.1.0d0) then
                    call v_rand(rnd,1,seed)
                    if (rnd.gt.fract) then
                        Rl(ishift+jreal) = 1
                        goto 500
                    endif
                endif
                nj = 0
c --- skip checking if only one contact
c --- form the entire fibers
                do k = 1,ecount
                   kreal = eReal(k)
                   if (kreal.eq.jreal) then
                       nj = nj + 1
                       jCont(nj) = k
                   endif
                enddo
c --- in ver=2 we only retain contacts which are initially closed
                if (ver.le.2) then
                    do i1 = 1,ni
c --- iCont is updated once per cycle and can be changed within, thus we need abs()
                       iCont(i1) = - abs(iCont(i1))
                    enddo
                    do j1 = 1,nj
                       jCont(j1) = - jCont(j1)
                    enddo
                endif
c --- go through all the segment and insure there is no problems                
                do i1 = 1,ni
                   iFib = abs(iCont(i1))

                   call MoveCoords (xi(1),yi(1),zi(1),
     &                              NdXYZ(1,1),ElNd(1,iFib),
     &                              ndpel(iFib))
                   do j1 = 1,nj
                      jFib = abs(jCont(j1))
                      call MoveCoords (xj(1),yj(1),zj(1),
     &                                 NdXYZ(1,1),ElNd(1,jFib),
     &                                 ndpel(jFib))

                      call ParBeam2ndFast (rvr(1),
     &                                     xi(1),yi(1),zi(1),
     &                                     xj(1),yj(1),zj(1),
     &                                     tol,epsil,cont,tols)
                          
                      if (ver.le.2.and.abs(epsil).lt.tol) then
c --- indicate that these contacts will be formed
                          iCont(i1) = abs(iCont(i1))
                          jCont(j1) = abs(jCont(j1))
                          if (i1.gt.1) iCont(i1-1)=abs(iCont(i1-1))
                          if (i1.lt.ni) iCont(i1+1)=abs(iCont(i1+1))
                          if (j1.gt.1) jCont(j1-1)=abs(jCont(j1-1))
                          if (j1.lt.nj) jCont(j1+1)=abs(jCont(j1+1))
                      endif
                   enddo
                enddo
c --- if this is the first contact for this pair of beams
                rcount = rcount + 1

                cRl = rcount
c --- get Ex the for second fiber
                Ex2 = eEx(j)
c --- store real data for this foundation
                call rlput (rcount,5,rvr(1))
c --- decide which one should be slave (stiffer, thicker)
                if (ni.lt.nj) then
                    ET1 = pareType
                    ET2 = 2 
                elseif (ni.gt.nj) then
                    ET1 = 2 
                    ET2 = pareType
                else
                    if (length1.lt.length2) then
                        ET1 = pareType
                        ET2 = 2 
                    else
                        ET1 = 2 
                        ET2 = pareType
                    endif
                endif
                ET1 = 3
                ET2 = 2
                
c --- the desision about the contact pair is defined in the sign of cRl
                Rl(ishift+jreal) = 4 ! mark as having contact
                elData(EL_REAL) = cRl

c --- full contact per pair (for slinding case)
                elData(EL_TYPE) = ET1
                do i1 = 1,ni
                   k = iCont(i1)
                   if (k.le.0) cycle
c --- store first series of segments
                   elem = elem + 1
                   call elmput (elem,elData(1),ndpel(k),ElNd(1,k))
                enddo

                elData(EL_TYPE) = ET2
                do j1 = 1,nj
                   k = jCont(j1)
                   if (k.le.0) cycle
c --- store second series of segments
                   elem = elem + 1
                   call elmput (elem,elData(1),ndpel(k),ElNd(1,k))
                enddo
            endif

 500     continue
 400  continue
      
      goto 996
      
c --- REPORT THE STATS >>> NEW
     
300   call Buf_Ini
      jobnam = ''
      fdir = ''
      call getjobinfo ('JOBN',jobnam,i)
      call getjobinfo ('CDIR',fdir,i)

c --- make new file for the stats
      call MakeName (jobnam,'cst     ',fdir,pname,kerr)
      
      nott = 45
      ioerr = Sys_Open (nott,pname,'UN','SERIAL','FORM','WRIT','BEGIN')

      if (ioerr .ne. 0) goto 998
      
      debug = 0
      nRl = rcount
      progress = 10
      
      do 110 i = 1,ecount
c --- get Ex for the first fiber
         ireal = eReal(i)
         ishift = (ireal-1)*nRl
         rvr(1) = 0.0d0 !eR(i)
         ni = 0
c --- skip checking if only one contact
c --- form the entire fibers
         do k = 1,ecount
            kreal = eReal(k)
            if (kreal.eq.ireal) then
                ni = ni + 1
                iCont(ni) = k
            endif
         enddo
         
         if ((i*100/ecount).ge.progress) then
             write (iott,*) 'Done',progress,'%'
             progress = progress + 10
         endif

         do 210 j = i+1, ecount
            jreal = eReal(j)
c --- do not consider self contact
            cRl = Rl(ishift+jreal)
            rvr(2) = 0.0d0 !eR(j)
            if (cRl.ne.0) goto 210
c --- find out if fibers overlap (super fast)
            if (.not.ElFibOverlapXY(MxMn(1,1,1),ireal,jreal)) then
                Rl(ishift+jreal) = 1
                goto 210
            endif

c --- form the entire fiber
            nj = 0
            do k = 1,ecount
               kreal = eReal(k)
               if (kreal.eq.jreal) then
                   nj = nj + 1
                   jCont(nj) = k
               endif
            enddo
c --- go through all the segment and insure there is no problems                
            minD = HUGE
            do i1 = 1,ni
               iFib = abs(iCont(i1))
               call MoveCoords (xi(1),yi(1),zi(1),
     &                          NdXYZ(1,1),ElNd(1,iFib),
     &                          ndpel(iFib))

               do j1 = 1,nj
                  jFib = abs(jCont(j1))
                  call MoveCoords (xj(1),yj(1),zj(1),
     &                             NdXYZ(1,1),ElNd(1,jFib),
     &                             ndpel(jFib))
                  call BeamToBeamXY (rvr(1),xi(1),yi(1),zi(1),
     &                               xj(1),yj(1),zj(1),HUGE,
     &                               epsil,debug,tols)
                  minD = min(minD,abs(epsil))
                enddo
            enddo
            if (abs(minD).lt.HUGE) write (nott,*) ireal, jreal, minD
            Rl(ishift+jreal) = 4 ! mark as having contact
 210     continue
 110  continue
      i = Sys_Close(nott,'KEEP')
      
996   call MemFree (PtrNdXYZ)
      call MemFree (PtrRl)
      call MemFree (PtrMxMn)
      call MemFree (PtriCont)
      call MemFree (PtrjCont)
      call MemFree (PtreLength)
997   call MemFree (PtrElNd)
      call MemFree (PtreEx)
      call MemFree (PtreR)
      call MemFree (PtreReal)      
      call MemFree (Ptrndpel)
      
      user06 = 0
      goto 999

998   continue

      user06 = -999

999   continue

      return
      end

      function ElFibOverlap (MxMn,i1,i2,tol)
c --- primary function: find out if fibers overlap
c
#include "impcom.inc"
c
      logical ElFibOverlap
      integer i1, i2
      double precision MxMn(3,2,*), tol, aL, aU

      ElFibOverlap = .false.
c --- try to build overlap volume
      aU = min(MxMn(1,2,i1),MxMn(1,2,i2))
      aL = max(MxMn(1,1,i1),MxMn(1,1,i2))
 
      if (aL.gt.(aU+tol)) goto 999
 
      aU = min(MxMn(2,2,i1),MxMn(2,2,i2))
      aL = max(MxMn(2,1,i1),MxMn(2,1,i2))

      if (aL.gt.(aU+tol)) goto 999
      
      aU = min(MxMn(3,2,i1),MxMn(3,2,i2))
      aL = max(MxMn(3,1,i1),MxMn(3,1,i2))

      if (aL.gt.(aU+tol)) goto 999
      
      ElFibOverlap = .true.
 
 999  continue
 
      return
      end

      function ElFibOverlapXY (MxMn,i1,i2)
c --- primary function: find out if fibers overlap
c
#include "impcom.inc"
c
      logical ElFibOverlapXY
      integer i1, i2
      double precision MxMn(3,2,*), aL, aU

      ElFibOverlapXY = .false.
c --- try to build overlap volume
      aU = min(MxMn(1,2,i1),MxMn(1,2,i2))
      aL = max(MxMn(1,1,i1),MxMn(1,1,i2))
      if (aL.gt.aU) goto 999
 
      aU = min(MxMn(2,2,i1),MxMn(2,2,i2))
      aL = max(MxMn(2,1,i1),MxMn(2,1,i2))

      if (aL.gt.aU) goto 999
      
      ElFibOverlapXY = .true.
 
 999  continue
 
      return
      end

      subroutine MoveCoords (xa,ya,za,NdXYZ,ElNd,ndpel)
c --- primary function: find out if fibers overlap
c
#include "impcom.inc"
c
      integer          i,ndpel,ElNd(3)
      double precision xa(3),ya(3),za(3),NdXYZ(3,*)
      do i = 1,ndpel
         xa(i) = NdXYZ(1,ElNd(i))
         ya(i) = NdXYZ(2,ElNd(i))
         za(i) = NdXYZ(3,ElNd(i))
      enddo
c --- put the third node in the middle
      if (ndpel.eq.2) then
          xa(3) = 0.5d0*( NdXYZ(1,ElNd(1)) + NdXYZ(1,ElNd(2)) )
          ya(3) = 0.5d0*( NdXYZ(2,ElNd(1)) + NdXYZ(2,ElNd(2)) )
          za(3) = 0.5d0*( NdXYZ(3,ElNd(1)) + NdXYZ(3,ElNd(2)) )
      endif
      return
      end

     
c      n1 = node(4987.2,1000.0,424.81)
c      n2 = node(5105.9,854.09,442.39)
c      nsel,,,,n1
c      nsel,a,,,n2
c      esln
c      esel,r,type,,1
c      elist
