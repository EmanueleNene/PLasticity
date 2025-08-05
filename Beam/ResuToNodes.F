      subroutine ResuToNodes(numElemNodes,numDomainIntPts,
     &                       numSecNodes,nDofNode,secDataR,offY, offZ,
     &                       nResults,eDisp, nBasis, nBasis_0, yzScale, 
     &                       pResults, nCrnr, khygro, thn, thp, 
     &                       keyplsL, warpFun,lenSecDataR)
c --- primary function: calculate element output at element nodes
c
#include "impcom.inc"
c
      integer numElemNodes, numDomainIntPts, nDofNode,
     &        numSecNodes, jj, nod, kspt, domInt,ii,offY,offZ,
     &        nCrnr(*), khygro,kk,lenSecDataR
      integer*8  keyplsL

      double precision nResults(6,2,*),secDataR(lenSecDataR,*),
     &                 work(7),pResults(9,3,*),ZERO, THREE,ONE,HALF,
     &                 yzScale(*),eDisp(nDofNode,4),yzup(2),
     &                 nBasis(3,3,3),nBasis_0(3,3,3),thn(2,*),thp(3,*),
     &                 warpFun(*)
      logical isWarp

      parameter (ZERO=0.0d0, ONE=1.0d0, THREE=3.0d0, HALF=0.5d0)

#if defined(DEBUG)
      call SubBeg('ResuToNodes')
#endif

      isWarp = .FALSE.
      if (nDofNode.eq.7) isWarp = .TRUE.
      if (numElemNodes.eq.2) then
          work(1) = ONE
          work(2) = ONE
          work(3) = ONE
          work(4) = ONE
      else
         if (keyplsL.eq.0) then
             work(1) = HALF*(ONE+sqrt(THREE))
             work(2) = HALF*(ONE-sqrt(THREE))
         else
             work(1) = 1.0d0
             work(2) = 0.0d0
         endif
         work(3) = work(2)
         work(4) = work(1)
         work(5) = HALF
         work(6) = HALF
      endif
      call v_zero(thp(1,1),numSecNodes*3)
      call v_zero(pResults(1,1,1),9*3*numSecNodes)
      do nod = 1,numElemNodes
         do kspt = 1,numSecNodes
            if (nCrnr(kspt).eq.1) then
                ii = (nod-1)*2
                work(7) = ZERO
                do domInt = 1,numDomainIntPts
                   do jj = 1,3
                      pResults(jj+3,nod,kspt) = pResults(jj+3,nod,kspt)+
     &                     work(ii+domInt)*nResults(jj,domInt,kspt)
                      pResults(jj+6,nod,kspt) = pResults(jj+6,nod,kspt)+
     &                     work(ii+domInt)*nResults(jj+3,domInt,kspt)
                   enddo
                   work(7) = work(7)+work(ii+domInt)*yzScale(domInt)
                   if (khygro.eq.1)
     &                 thp(nod,kspt) = thp(nod,kspt)+
     &                 work(ii+domInt)*thn(domInt,kspt)
                enddo
                   kk = 1
                   yzup(1) = secDataR(offY+kspt-1,kk)
                   yzup(2) = secDataR(offZ+kspt-1,kk)
                   do jj = 1,3
                      pResults(jj,nod,kspt)=eDisp(jj,nod)+
     &                     (work(7)*nBasis(jj,2,nod)-
     &                      nBasis_0(jj,2,nod))*yzup(1)+
     &                     (work(7)*nBasis(jj,3,nod)-
     &                      nBasis_0(jj,3,nod))*yzup(2)
                  enddo
            endif
         enddo
      enddo

      if (isWarp) then
          do nod=1,numElemNodes
             do kspt=1,numSecNodes
                if (nCrnr(kspt).eq.1) then
                    do jj=1,3
                       pResults(jj,nod,kspt) = pResults(jj,nod,kspt)+
     &                       nBasis(jj,1,nod)*eDisp(7,nod)*warpFun(kspt)
                    enddo
                endif
             enddo
          enddo
      endif


#if defined(DEBUG)
      call SubEnd('ResuToNodes')
#endif

      return
      end
