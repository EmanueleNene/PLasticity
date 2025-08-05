      subroutine RotMtx(Ang,Rot,RotT,T,TT,ElDim,ElDof)
c --- primary function: Finds the rotation matrix
c     Local  = Rot  * Global
c     Global = RotT * Local
c
c input arguments:
c   Ang     (dp,sc,in)  -Angle between two beams
c   ElDim    (int,sc,in) -Dimension of the element
c   ElDof    (int,sc,in) -Dimension of element rotation matrix (ElDim*H_nNode)
c
c output arguments:
c   Rot      (dp,arr(ElDim,ElDim),out) -2D rotation matrix for point
c   T        (dp,arr(ElDof,ElDof),out) -2D rotation matrix for 4-node element
c   TT       (dp,arr(ElDof,ElDof),out) -Transpose of T
c
#include "impcom.inc"
c

      
      integer            ElDim,ElDof
      double precision   Ang
      double precision   Rot(3,3),RotT(3,3)
      double precision   T(60,60),TT(60,60)
      
      double precision   c,s
      integer            i, loc1, loc2
      
      double precision   Zero, One
      parameter          (Zero = 0.0d0, One=1.0d0)
      
      c = cos(Ang)
      s = sin(Ang)
      
      RotT(1,1) =  c
      RotT(1,2) = -s
      RotT(2,1) =  s
      RotT(2,2) =  c
      if (ElDim .eq. 3) then
          RotT(1,3) =  Zero
          RotT(2,3) =  Zero
          RotT(3,1) =  Zero
          RotT(3,2) =  Zero
          RotT(3,3) =  One
      endif
            
      call v_zero (TT(1,1),3600)
      loc1 = 1
      do i=1,ElDof/ElDim
          loc2 = loc1+ElDim-1
          TT(loc1:loc2,loc1:loc2) = RotT
          loc1 = loc1 + ElDim
      enddo
      
      Rot = RotT
      Rot(1,2) =  s
      Rot(2,1) = -s
      
      call v_zero (T(1,1),3600)
      loc1=1
      do i=1,ElDof/ElDim
          loc2 = loc1+ElDim-1
          T(loc1:loc2,loc1:loc2) = Rot
          loc1 = loc1 + ElDim
      enddo

      return
      end