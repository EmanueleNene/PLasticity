      subroutine CalSect_4(nnodi,xy,Gxy,eStiff,eForce,gProc,
     &                     Summary,ErrorFlag)
c --- primary function: compute section data
c
#include "impcom.inc"
c
      integer MNODI,eSize,nnodi
      parameter (MNODI = 9, eSize=45)
      double precision eStiff(eSize),eForce(MNODI),xy(2,*),yz(2),
     &                 parInt(MNODI,3),eJacob(4),spaDer(MNODI,2),
     &                 Gxy,ky,kz
      double precision ZERO,ONE,deltaV,detJ,work,t1,t2,weight,
     &                 Summary(*),vy,vz
      integer kintk, knode, k1, k2, k3, nloc, ErrorFlag, gProc
      parameter (ZERO = 0.0d0, ONE=1.0d0)
      integer
     &     pIxx,
     &     pIxy,
     &     pIyy
      parameter (
     &     pIxx = 6,
     &     pIxy = 7,
     &     pIyy = 8 )

c --- torsion stress function
c
      call v_zero(eStiff(1),eSize)
      call v_zero(eForce(1),nnodi)
      if (Gxy.eq.ZERO) return
      ky = ZERO
      kz = ZERO
      if (gProc.gt.2) then
          work=2.0d0*(Summary(pIxx)*Summary(pIyy)-Summary(pIxy)**2)
          if (gProc.eq.3) then
              vy=ZERO
              vz=ONE
          else
              vy=ONE
              vz=ZERO
          endif
          ky = (vy*Summary(pIxx)-vz*Summary(pIxy))/work
          kz = (vz*Summary(pIyy)-vy*Summary(pIxy))/work
      endif
      do kintk = 1,9
         call CalSect_7(kintk, nnodi,parInt(1,1), weight)
         eJacob(1) = ZERO
         eJacob(2) = ZERO
         eJacob(3) = ZERO
         eJacob(4) = ZERO
         k3=0
         do k1=2,3
            do k2=1,2
               k3=k3+1
               do knode=1,nnodi
                  eJacob(k3)=eJacob(k3)+parInt(knode,k1)*xy(k2,knode)
               enddo
            enddo
         enddo
         deltaV=eJacob(1)*eJacob(4)-eJacob(2)*eJacob(3)
         if (deltaV .le. ZERO) then
             ErrorFlag = 1
             call erhandler ('CalSect_4',5000,3,
     &            'Section element has negative area %G. %/
     &            Check your input carefully.', deltaV, ' ')
             goto 999
         endif
         detJ = ONE/deltaV
         work = eJacob(1)
         eJacob(1) = eJacob(4)*detJ
         eJacob(2) =-eJacob(2)*detJ
         eJacob(3) =-eJacob(3)*detJ
         eJacob(4) = work*detJ
         deltaV=deltaV*weight*Gxy
         yz(1)=0.0d0
         yz(2)=0.0d0
         do knode=1,nnodi
            spaDer(knode,1)=
     &           eJacob(1)*parInt(knode,2)+eJacob(2)*parInt(knode,3)
            spaDer(knode,2)=
     &           eJacob(3)*parInt(knode,2)+eJacob(4)*parInt(knode,3)
            yz(1) = yz(1) + parInt(knode,1)*xy(1,knode)
            yz(2) = yz(2) + parInt(knode,1)*xy(2,knode)
         enddo
         nloc=0
         do k1=1,nnodi
            t1 = spaDer(k1,1)
            t2 = spaDer(k1,2)
            if (gProc.le.2) then
                eForce(k1)=eForce(k1)+(yz(2)*t1-yz(1)*t2)*deltaV
            else
                eForce(k1)=eForce(k1)+2.0d0*parInt(k1,1)*
     &                    (yz(2)*kz+yz(1)*ky)*deltaV
            endif
            do k2=1,k1
               nloc = nloc+1
               eStiff(nloc) = eStiff(nloc)+
     &              (t1*spaDer(k2,1)+t2*spaDer(k2,2))*deltaV
            enddo
         enddo
      enddo
  999 continue

      return
      end


