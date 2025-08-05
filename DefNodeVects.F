      subroutine DefNodeVects (NumEl,eorder,NumNd,BacLen)
c --- primary function: define the M_FrwdList and M_BackList vectors
c
c  input arguments:
c     NumNd      (int,sc,in)     - maximum node number
c
c  output arguments:
c     BacLen     (int,sc,out)    - number of nodes actually used
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
c
      integer   elmget, etyget

      integer  NumEl,eorder(NumEl),NumNd, BacLen,
     & ityold,ielm,nnod,ityp,i,j, nodes(W_NdMax),elData(EL_DIM),
     & iel,nstif, elChar(IELCSZ), jtyp

#if defined(DEBUG)
      call SubBeg('DefNodeVects')
#endif

      call vi_zero (M_FrwdList(1),NumNd)
      BacLen = 0
      ityold = 0
      
      do iel = 1,NumEl
         ielm = eorder(iel)
         if (ielm .le. 0) cycle
         call vi_zero (nodes(1),W_NdMax)
         nnod = elmget(ielm,elData(1),nodes(1))
         if (nnod .le. 0) cycle
         ityp = elData(EL_TYPE)
         if (ityp .ne. ityold) then
             i = etyget(ityp,elChar(1))
             ityold = ityp
         endif
         jtyp = elChar(JETYP)
         nstif = elChar(NMNDNE)
c --- go through the element
         do i = 1,nstif
            j = nodes(i)
c --- ignore the extra nodes appeared from the eliminations
            if (jtyp.eq.W_eMaster.and.i.ge.3) cycle
            if (j .gt. 0) then
                if (M_FrwdList(j) .le. 0) then
                    BacLen = BacLen + 1
                    M_FrwdList(j) = BacLen
                    M_BackList(BacLen) = j
                endif
            endif
         enddo
      enddo

      if (NumNd.ne.BacLen) then
          call erhandler ('DefNodeVects',5000,2,'Detected presence of %i
     &                     unconnected nodes. Issue ALLS and NDEL,ALL, 
     &                     if these are not orrientaiton nodes.',
     &                     dble(NumNd-BacLen),' ') 
      endif
      
#if defined(DEBUG)
      call SubEnd('DefNodeVects')
#endif

      return
      end
