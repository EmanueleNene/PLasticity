      subroutine SetDisp (nr,ls,kelout,dStiff,dLV,DofBit)
c --- primary function: adjust element load vector for prescribed disps
c
c  input arguments:
c     nr       (int,sc,in)         - matrix and lv size
c     ls       (int,ar(nr),in)     - element dof vector
c     kelout   (int,ar(10),in)     - keys indicating available matrices and load vectors
c     dStiff   (dp,ar(nr,nr),in)   - matrices
c     DofBit   (int,ar(*),in)      - BCs bit patterns
c
c  output arguments:
c     dLV      (dp,ar(nr),inout)   - load vectors, adjusted
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "in_mem.inc"
c
      integer nr,ls(nr),kelout(10),DofBit(*),node,idf,kdf,i
      double precision  dStiff(nr,nr),dLV(nr),dval(4)

c --- loop on element rows
      do 100 i = 1,nr

         kdf = ls(i)
         if (kdf.le.0) goto 100
c --- is this degree of freedom a specified displacement?
         if (.not.(btest(DofBit(kdf),W_DispBit))) goto 100
         if (.not.(btest(DofBit(kdf),W_TranBit))) goto 100

c --- get non-zero displacements
         call InqDof (kdf,node,idf)
         call GetD (node,idf,dval(1))

         dval(1) = dval(1) - M_Disp(kdf)
         
         if (dval(1).ne.0.0d0) then
             call v1pcv2_v1 (dLV(1),dStiff(1,i),nr,-dval(1))
             kelout(5) = 1
         endif
         
  100 continue

      return
      end
