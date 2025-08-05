      subroutine DynMatAdj (nr,kelfil,delu,vel,acel,eMass,eStiff,
     &                      fExt,fInt,kReform)
c --- primary function: modify element matrices for transient effects
c
c  input arguments:
c     kelfil  (int,ar(6),in)       - element matrix key
c     delu    (dp,ar(nr),in)       - time step displacement increment
c     vel     (dp,ar(nr),in)       - velocity
c     acel    (dp,ar(nr),in)       - acceleration
c     eMass   (dp,ar(nr,nr) in)    - mass matrix
c     eStiff  (dp,ar(nr,nr),inout) - stiffness and equivalent stiffness matrix
c     fExt    (dp,ar(nr),inout)    - load and equivalent load vector
c     kReform (int,sc,in)          - stiffness matrix reform flag
c
c  output arguments:
c     eStiff (dp,ar(nr,nr),inout) - stiffness and equivalent stiffness matrix
c     fExt   (dp,ar(nr),inout)    - load and equivalent load vector
c     fInt   (dp,ar(nr),inout)    - newton-raphson force vector
c
#include "impcom.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
c
      integer  nr,kelfil(6),kReform,i
      double precision
     & delu(nr),vel(nr),acel(nr),eMass(nr,nr),
     & eStiff(nr,nr),fExt(nr),fInt(nr),con,alpha
      equivalence (alpha,I_Consts(2))

#if defined(DEBUG)
      call SubBeg('DynMatAdj')
#endif 

c --- second order integration
c --- mass matrix adjustments
      if (kelfil(2).ge.1) then
          kelfil(1) = 1
          kelfil(5) = 1
c --- set up system equations
c --- process one column at a time
          do i = 1,nr
             if (alpha.gt.1.0d-6) then
c --- implicit load vector adjustment
                 if (A_DynFirstPass.eq.1.and.A_InLineSrch.eq.0) then
                     con = I_Prm(5)*vel(i) + I_Prm(3)*acel(i)
c --- implicit direct iteration adjustment
                 else
                     con = -acel(i)
                 endif
                 call v1pcv2_v1 (fInt(1),eMass(1,i),nr,-con)
c --- implicit stiffness adjustment
                 if (kReform.eq.1) then
                     call v1pcv2_v1 (eStiff(1,i),eMass(1,i),nr,I_Prm(1))
                 endif
             else
c --- explicit load vector adjustment
                 call v1pcv2_v1 (fExt(1),eMass(1,i),nr,delu(i)*I_Prm(3))
c --- explicit stiffness adjustment
                 if (kReform.eq.1) then
                     call v1pcv2_v1 (eStiff(1,i),eMass(1,i),nr,I_Prm(1))
                 endif
             endif
          enddo
      endif

#if defined(DEBUG)
      call SubEnd('DynMatAdj')
#endif 
      
      return
      end
