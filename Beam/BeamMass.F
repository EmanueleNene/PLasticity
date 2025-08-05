      subroutine BeamMass (numDomainIntPts,numElemNodes, nDofNode,
     &                     Gauss2, Gauss3, xBeg, Gaussw2,
     &                     Gaussw3, hbasis, sectSummary,
     &                     pIzz, pIyy, pIyz, pArea, nGenStr, cDensity,
     &                     density, offset, lump, stab,
     &                     nUsrDof, secStiff, gradMat, work, eMass)
c --- primary function: form mass matrix for beam
c
c  input arguments:
c
c     numDomainIntPts   (int,sc,in)    - number of domain integration points
c     xBeg        (dp,ar(nr,5),in)     - nodal coords (orig) and rot angles
c     numElemNodes(int,sc,in)          - number of element nodes
c     nUsrDof     (int,sc,in)          - matrix and LV size
c     Gauss2      (dp,ar(2),in)        - location of gaussian integration points
c     Gauss3      (dp,ar(3),in)        - location of gaussian integration points
c     Gaussw2     (dp,ar(2),in)        - weight factors of at gaussian points
c     Gaussw3     (dp,ar(3),in)        - weight factors of at gaussian points
c     hbasis      (dp,ar(3,3,3),in)    - base vectors at gaussian points
c     sectSummary (int,ar(20,8),in)    - section data
c     pIzz        (int,sc,in)          - Izz location in sectSummary
c     pIyy        (int,sc,in)          - Iyy location in sectSummary
c     pIyz        (int,sc,in)          - Iyz location in sectSummary
c     pArea       (int,sc,in)          - Area location in sectSummary
c     offset      (dp,ar(3,3,3),in)    - off set
c     lump        (int,sc,in)          - key for lumped mass
c     stab        (int,sc,in)          - key indicating the call from stabilization
c     nDofNode    (int,sc,in)          - number of DOFs per node
c     cDensity    (dp,sc,in)           - density for added mass
c     density     (dp,sc,in)           - beam density
c
c  output arguments:
c     eMass  (dp,ar(nUsrDof,nUsrDof),inout)- m matrix (A_MatrReq(2))
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer
     &     numDomainIntPts, numElemNodes, pIzz, pIyy, pIyz,
     &     pArea, nGenStr, lump, nDofNode, nUsrDof, stab

      double precision
     &     Gauss2(2), Gauss3(3), xBeg(6,4), Gaussw2(2),Gaussw3(3),
     &     hbasis(3,3,3), sectSummary(20,8), cDensity, density,
     &     offset(3,3,3),eMass(nUsrDof,nUsrDof)

      double precision ONE, ZERO, HALF, TWO
      parameter (ZERO = 0.0d0, HALF = 0.5D0, ONE = 1.0d0, TWO = 2.0d0)
      integer
     &     domInt, domIntP, j, nod, i, jj, kk, ll

      double precision
     &     shpe(3,2), dL, dxdg(3), dV, n1(3), n2(3), dxDx(3,3),
     &     secStiff(nGenStr,nGenStr), gradMat(nGenStr,nDofNode*3),
     &     work(*),Izz,Iyy,Iyz

#if defined(DEBUG)
      call SubBeg('BeamMass')
#endif

c --- integration loop
      do domInt = 1,numDomainIntPts+1
         domIntP = 1
         if (numElemNodes.eq.2) then
             shpe(1,1) = HALF*(ONE-Gauss2(domInt))
             shpe(2,1) = HALF*(ONE+Gauss2(domInt))
             shpe(1,2) =-HALF
             shpe(2,2) = HALF
         elseif(numElemNodes.eq.3) then
             shpe(1,1) = HALF*(Gauss3(domInt)-ONE)*Gauss3(domInt)
             shpe(2,1) = HALF*(Gauss3(domInt)+ONE)*Gauss3(domInt)
             shpe(3,1) = ONE-Gauss3(domInt)*Gauss3(domInt)
             shpe(1,2) = Gauss3(domInt)-HALF
             shpe(2,2) = Gauss3(domInt)+HALF
             shpe(3,2) =-TWO*Gauss3(domInt)
         endif
         dL = ZERO
         call v_zero (dxdg(1),3)
         do j = 1,3
            do nod = 1,numElemNodes
               dxdg(j) = dxdg(j)+shpe(nod,2)*xBeg(j,nod)
            enddo
            dL = dL + dxdg(j)*dxdg(j)
         enddo
         dL = sqrt(dL)
         if (numElemNodes.eq.2) then
             dV = dL*Gaussw2(domInt)
         elseif (numElemNodes.eq.3) then
             dV = dL*Gaussw3(domInt)
         endif
         n1(1) = hbasis(1,2,domInt)
         n1(2) = hbasis(2,2,domInt)
         n1(3) = hbasis(3,2,domInt)
         n2(1) = hbasis(1,3,domInt)
         n2(2) = hbasis(2,3,domInt)
         n2(3) = hbasis(3,3,domInt)
         
         Izz = sectSummary(pIzz,domIntP)
         Iyy = sectSummary(pIyy,domIntP) 
         Iyz = sectSummary(pIyz,domIntP)
         
         dxDx(1,1)=
     &        Izz*(n1(3)**2+n1(2)**2)+
     &        Iyy*(n2(3)**2+n2(2)**2)+
     &        TWO*Iyz*(n1(3)*n2(3)+n1(2)*n2(2))
         dxDx(2,2)=
     &        Izz*(n1(3)**2+n1(1)**2)+
     &        Iyy*(n2(3)**2+n2(1)**2)+
     &        TWO*Iyz*(n1(3)*n2(3)+n1(1)*n2(1))
         dxDx(3,3)=
     &        Izz*(n1(2)**2+n1(1)**2)+
     &        Iyy*(n2(2)**2+n2(1)**2)+
     &        TWO*Iyz*(n1(2)*n2(2)+n1(1)*n2(1))
         dxDx(1,2)=
     &        -Izz*(n1(1)*n1(2))
     &        -Iyy*(n2(1)*n2(2))
     &        -Iyz*(n1(1)*n2(2)+n2(1)*n1(2))
         dxDx(1,3)=
     &        -Izz*(n1(3)*n1(1))
     &        -Iyy*(n2(3)*n2(1))
     &        -Iyz*(n1(3)*n2(1)+n1(1)*n2(3))
         dxDx(2,3)=
     &        -Izz*(n1(3)*n1(2))
     &        -Iyy*(n2(3)*n2(2))
     &        -Iyz*(n1(3)*n2(2)+n2(3)*n1(2))

         dxDx(2,1) = dxDx(1,2)
         dxDx(3,1) = dxDx(1,3)
         dxDx(3,2) = dxDx(2,3)
         call v_zero (secStiff(1,1),nGenStr*nGenStr)
         if (A_StabFact.gt.ZERO .or. stab.eq.0) then
             do i = 1,3
                secStiff(i,i) = sectSummary(pArea,domIntP)*cDensity
                do j = 1,3
                   secStiff(i+3,j+3) = dxDx(i,j)*density
                enddo
             enddo
         else
c --- viscous damping in the Z-direction only
             secStiff(3,3) = sectSummary(pArea,domIntP)*cDensity
         endif
         

         call v_zero(gradMat(1,1),nGenStr*nUsrDof)
         do j = 1,numElemNodes
            jj = (j-1)*nDofNode
            do i = 1,6
               gradMat(i,jj+i) = shpe(j,1)
            enddo
            jj = jj + 3
            do kk = 1,3
               do ll = 1,3
                  gradMat(kk,jj+ll) = offset(kk,ll,j)*shpe(j,1)
               enddo
            enddo
         enddo
         call cAtxBxA_C (gradMat(1,1),secStiff(1,1),eMass(1,1),
     &                   nGenStr,nGenStr,nUsrDof,nGenStr,nUsrDof,
     &                   work(1),dV)
      enddo

      call m_symLtoU(eMass(1,1),nUsrDof,nUsrDof)

      if (lump.eq.1) then
          if (numElemNodes.eq.2) then
              call m_diag(eMass(1,1),nUsrDof,nUsrDof,0)
          else
              call m_diag(eMass(1,1),nUsrDof,nUsrDof,1)
          endif
      endif      

#if defined(DEBUG)
      call SubEnd('BeamMass')
#endif

      return
      end
