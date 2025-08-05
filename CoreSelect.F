      subroutine CoreSelect(BeamStat,BeamConn,cBeams,
     &                      NumR,NumEl,nconn,maxconn)
c --- primary function: select the core of fibers
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "cont_info.inc"
#include "dof_info.inc" 
#include "solu_info.inc" 
c
      integer elmget, etyget, vi_compress
      PTRFTN  MemAlloc
      integer nconn, maxconn, NumR, NumEl,
     &        BeamStat(NumR),BeamConn(maxconn)
     
      integer i, j, nnod, select, elData(EL_DIM), ls(W_NrMax), nr,
     &        elChar(IELCSZ), nodes(W_NdMax), ityp, itypold, ii, jj, nn,
     &        iconn, elem, ireal, areal, kdf, etyp, trh,
     &        Bit(2), ipass, ntrans, npass, RotDofInc, k, numo
      
      logical  IfBit, Removed
      pointer (ptrSele,Sele)
      pointer (ptrSeleIdx,SeleIdx)
      integer  Sele(nconn,4), SeleIdx(NumR),cBeams(NumEl)


#if defined(DEBUG)
      call SubBeg('CoreSelect')
#endif

      ptrSele = MemAlloc(4*nconn,W_LenInt,'ptrSele         ')
      ptrSeleIdx = MemAlloc(NumR,W_LenInt,'ptrSeleIdx      ')
      call vi_zero (SeleIdx(1),NumR)

      iconn = 1
c --- copy selected elements to tmp arrays
      do i = 1,nconn
c --- the entry can be negative if the element touches the surface 
c     see FillConn.F for details
         j = abs(BeamConn(iconn))
         Sele(i,1) = BeamConn(iconn)
c --- place for iconn
         Sele(i,3) = iconn
c --- the original number of contacts
         numo = BeamConn(iconn+1)
c --- the compressed number of contacts
         BeamConn(iconn+1) = vi_compress (BeamConn(iconn+2),numo)
c --- place for compressed number of contacts
         Sele(i,4) = BeamConn(iconn+1)
c --- form selection index
         SeleIdx(j) = i
c --- next pointer (also available in BeamStat(*) for this matter)
         iconn = iconn + numo + 2
      enddo
      
      call vi_move(Sele(1,1),Sele(1,2),nconn)
      Bit(1) = W_DispBit
      Bit(2) = W_TranBit
      RotDofInc = 3
      ntrans = 0
      
      do 100 ipass = 1,2
c --- select fixed fibers
         itypold = 0
         etyp = 0
         do elem = 1,NumEl
c --- again, the entry can be negative if the element touches the surface 
c     in this case, the element was already selected in FillConn.F
            nnod = elmget(elem,elData(1),nodes(1))
            if (nnod.le.0) cycle
            ityp = elData(EL_TYPE)
            if (ityp .ne. itypold) then
                itypold = ityp
                ityp = etyget(ityp,elChar(1))
                etyp = elChar(JETYP)
            endif
c --- only beams are processed
            if (elChar(JBEAM).eq.0) cycle 
            ireal = elData(EL_REAL)
            if (ireal.le.0) cycle
c --- check if it has a contact for the first pass only
            if (ipass.eq.1) then
                if (cBeams(elem).eq.-1) then
                    cBeams(elem) = -ireal
                else
                    cBeams(elem) = ireal
                endif
            endif
            ii = SeleIdx(ireal)
c --- skip if it is already selected
            if (ii .le. 0) cycle
            if (Sele(ii,ipass).le.0 .and. cBeams(elem).lt.0) cycle
            call DefLsVect (elChar(1),nodes(1),ls(1),nr)
            IfBit = .true.
c --- check the bit in interest (only translational DOFs)
            do j = 1,nr,D_NumDof
               do k = 0,RotDofInc-1
                  kdf = ls(j+k)
                  if (kdf .gt. 0) then
                      if (btest(M_DofBits(kdf),Bit(ipass))) goto 110
                  endif
               enddo
            enddo
            IfBit = .false.
110         continue
c --- select the element by making the element number negative
            if (IfBit) then
c --- mark this beam as connected
                cBeams(elem) = -ireal
                Sele(ii,ipass) = -ireal
                ii = BeamStat(ireal)
c --- increase number of contacts for the fiber that touches the boundary
c     this fiber can contribute to the stiffness if it touches another fiber !!!
                BeamConn(ii+1) = BeamConn(ii+1) + 1
                if (ipass.eq.2) ntrans = ntrans + 1
            endif
         enddo
100   continue         
      
      npass = 2
c --- if there is no translational contraints, then we need only one ipass
      if (ntrans.eq.0) npass = 1
      
      do 200 ipass = 1,npass
c --- gradually select the core      
  210    continue
            select = 0
            do i = 1,nconn
               iconn = Sele(i,3)
c --- if element is selected
               if (Sele(i,ipass) .lt. 0) then
c --- mark all elements as selected which are connected to previously selected
                   nn = Sele(i,4)
                   do j = 1,nn
                      ireal = BeamConn(iconn+1+j)
                      if (ireal.le.0) cycle
                      ii = SeleIdx(ireal)
                      if (ii .le. 0) cycle
                      areal = Sele(ii,ipass)
                      if (areal .gt. 0) Sele(ii,ipass) = -areal
                   enddo
c --- mark the element as secondly passed for selection (zero means reachable here)
                   Sele(i,ipass) = 0
                   select = select + 1
               endif
            enddo
c --- no selection event occured? if so, then quit the loop
            if (select .eq. 0) goto 200
         goto 210

200   continue

c --- transfer only elements reachable from both ends
      do i = 1,nconn
         iconn = Sele(i,3)
         if (ntrans.gt.0) then
             Sele(i,1) = max(Sele(i,1),Sele(i,2))
         endif
c --- transfer beam number to Sele(*,2)
         Sele(i,2) = BeamConn(iconn)
         BeamConn(iconn) = Sele(i,1)
      enddo
      
c      write (8,*) 'do', BeamConn(1:maxconn)
c --- theshold on number of contacts per fiber
      trh = 1 !min(A_BisecNum+1,2)
c --- removed fibers with one "effective" contact per fiber
c     iterations required
 300  continue
c --- continue until this logical variable is false
            Removed = .true.
            do i = 1,nconn
c --- get the current index
               iconn = Sele(i,3)
c --- if selected and has just one contact
               if (BeamConn(iconn).eq.0 .and.
     &             BeamConn(iconn+1).le.trh) then
c --- it happened again then...
                   Removed = .false.
c --- mark for deletion with any positive number
                   BeamConn(iconn) = 999
c --- zero the number of contacts
                   BeamConn(iconn+1) = 0
c --- get the real number from the original copy
                   areal = Sele(i,2)
c --- get the connected fiber
                   ireal = BeamConn(iconn+2)
c --- remove the entry
                   BeamConn(iconn+2) = 0
c --- does this fiber exist?
                   if (ireal.le.0) goto 320
c --- if yes, get its table index
                   ii = BeamStat(ireal)
                   if (ii .le. 0) goto 320
c --- get the selection index to retrive the original number of contacts
c     it could have changed if the fiber is connected to the BCs boundary
                   jj = SeleIdx(ireal)
c --- this is the original number of contacts
                   nn = Sele(jj,4)
c --- find and remove the contact
                   do j = 1,nn
                      if (BeamConn(ii+1+j).eq.areal) then
                          BeamConn(ii+1+j) = 0
                          goto 310
                      endif
                   enddo
c --- compress the array
 310               numo = vi_compress(BeamConn(ii+2),nn)
c --- decrease the number of contacts accordingly
                   Sele(jj,4) = numo
                   BeamConn(ii+1) = numo
c --- no other contacts for this fiber? then delete it as well
                   if (BeamConn(ii+1).eq.0) BeamConn(ii) = 999
               endif
 320           continue
            enddo
      if (.not.Removed) goto 300

c      write (9,*) 'posle', BeamConn(1:maxconn)

      call MemFree(ptrSele)
      call MemFree(ptrSeleIdx)

#if defined(DEBUG)
      call SubEnd('CoreSelect')
#endif

      return
      end

