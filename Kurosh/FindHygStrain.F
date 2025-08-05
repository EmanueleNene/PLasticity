      subroutine FindHygStrain(BeamProps,Ang,Xi,BndEps,Eps)
#include "impcom.inc"
#include "beam_hyg.inc"

      
      double precision    BeamProps(11,2), Ang, Xi(2), BndEps(8,2)
      double precision    Eps(H_Row,H_Col)
      
      double precision    L1,L2, W1,W2, nis
      double precision    C1,C2, w
      integer             ind(2,8),i
      double precision    sgn(8)
      double precision    Hlf
      parameter           (Hlf = 0.5d0)
      
      L1      = BeamProps(1,1)
      L2      = BeamProps(1,2)
      
      W1      = BeamProps(3,1)
      W2      = BeamProps(3,2)
      
      nis = 1.0d0/sin(Ang)
      
      if (H_key(1) .eq. 2) then
          BndEps(1,1)  = BndEps(1,1)*W2*nis
          BndEps(1,2)  = BndEps(1,2)*W1*nis
          
c          Eps(1:4,1) = BndEps(1,1)/L1*(1.0d0-3.0d0*Xi(1))
c          Eps(5:8,1) = BndEps(1,1)/L1*(1.0d0+3.0d0*Xi(1))
          
c          Eps(1:4,2) = BndEps(1,2)/L2*(1.0d0-3.0d0*Xi(2))
c          Eps(5:8,2) = BndEps(1,2)/L2*(1.0d0+3.0d0*Xi(2))
          
          
          Eps(1:2,1) = BndEps(1,1)/L1*(1.0d0-3.0d0*Xi(1))
          Eps(5:6,1) = BndEps(1,1)/L1*(1.0d0-3.0d0*Xi(1))          
          Eps(3:4,1) = BndEps(1,1)/L1*(1.0d0+3.0d0*Xi(1))
          Eps(7:8,1) = BndEps(1,1)/L1*(1.0d0+3.0d0*Xi(1))          
          
          Eps(1:2,2) = BndEps(1,2)/L2*(1.0d0-3.0d0*Xi(2))
          Eps(5:6,2) = BndEps(1,2)/L2*(1.0d0-3.0d0*Xi(2))
          Eps(3:4,2) = BndEps(1,2)/L2*(1.0d0+3.0d0*Xi(2))
          Eps(7:8,2) = BndEps(1,2)/L2*(1.0d0+3.0d0*Xi(2))          
          
          
          
      else
          ind(1,1) = 1
          ind(2,1) = 4
          ind(1,2) = 2
          ind(2,2) = 3
          ind(1,3) = 2
          ind(2,3) = 3
          ind(1,4) = 1
          ind(2,4) = 4
          ind(:,5:8) = ind(:,1:4) + 4
          
          sgn(1:2) = -1.0d0
          sgn(3:4) = +1.0d0
          sgn(5:8) = sgn(1:4)
          
          w = W2 * nis
          C1 = (3.0d0*xi(1))*Hlf
          C2 = w/(2.0d0*(L1 + 1.0d0))
          do i=1,8
              Eps(i,1) = BndEps(ind(1,i),1)*sgn(i)*(C1-C2+sgn(i)*Hlf)
     &                  +BndEps(ind(2,i),1)*sgn(i)*(C1+C2+sgn(i)*Hlf)
          enddo
          Eps(:,1) = Eps(:,1) * w/L1
          
          w = W1 * nis
          C1 = (3.0d0*xi(2))*Hlf
          C2 = w/(2.0d0*(L2 + 1.0d0))
          do i=1,8
              Eps(i,2) = BndEps(ind(1,i),2)*sgn(i)*(C1-C2+sgn(i)*Hlf)
     &                  +BndEps(ind(2,i),2)*sgn(i)*(C1+C2+sgn(i)*Hlf)
          enddo
          Eps(:,2) = Eps(:,2) * w/L2
      endif
      
      return
      end