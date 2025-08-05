      subroutine ElemMatForm (elem,elChar,elData,eomask,nodes,nr,
     &                        kelreq,kelfil,xyzang,uelm,Kmtrx,Fnr,
     &                        elmass,center,lcSVRL,ndUpdate)
c --- primary function: set up the call to the element routines
c
c  input arguments:
c       elem    (int,sc,in)           - element number
c       elChar  (int,ar(*),in)        - array of element type characteristics
c       elData  (int,ar(EL_DIM),in)   - array of element data
c       eomask  (int,sc,in)           - bit pattern for element output
c       nodes   (int,ar(*),in)        - array of element node numbers
c       nr      (int,sc,in)           - matrix and lv size
c       kelreq (int,ar(10),in)        - matrix and load vector form requests
c       kelfil  (int,ar(10),in)       - keys indicating incoming mats and LVs
c       xyzang  (dp,ar(6,ndim),in)    - nodal coords (orig) and rotation angles
c       uelm    (dp,ar(nr,5),in)      - element nodal solution values
c
c  output arguments:
c       Kmtrx  (dp,ar(nr,nr,4),inout) - matrices
c       Fnr    (dp,ar(nr,2),out)      - load vectors
c       elmass (dp,sc,out)            - element mass
c       center (dp,ar(3),out)         - centroid location
c       ndUpdate(int,sc,out)          - keys indicating the update nodes      
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer*8  lcSVRL,lcerstL,edindxL(25)
      integer elem,elChar(*),elData(EL_DIM),nodes(*),nr,kelreq(10),
     & kelfil(10),kelout(10),i,eomask,ndUpdate
      double precision
     & xyzang(6,*),uelm(nr,*),Kmtrx(nr,nr,4),Fnr(nr,2),elmass,center(3),
     & elener(10),elvol

#if defined(DEBUG)
      call SubBeg('ElemMatForm')
#endif

c --- initialize the output form keys
      call vi_zero (kelout(1),10)
      call v_zero (elener(1),10)

c --- ZERO out variables
      do i = 1,25
         edindxL(i) = 0
      enddo
      lcerstL = 0

      call ElemDrivers (elem,elChar(1),elData(1),eomask,nodes(1),lcSVRL,
     &                  kelreq(1),nr,xyzang(1,1),uelm(1,1),kelout(1),
     &                  Kmtrx(1,1,1),Kmtrx(1,1,2),
     &                  Kmtrx(1,1,3),Kmtrx(1,1,4),
     &                  Fnr(1,1),elmass,center(1),elvol,elener(1),
     &                  edindxL(1),lcerstL,ndUpdate)
      
c --- the first pass establishes what matrices are on the file
      if (A_KeyFirstPass.eq.1) then
          call vi_move (kelout(1),kelfil(1),10)
      else
c --- 2 load vectors are always saved so that kelfil(5,6,7) can change
          call vi_move (kelout(5),kelfil(5),3)
      endif

      if (kelfil(1).eq.0) then
          call vi_move (kelout(1),kelfil(1),10)
      endif

#if defined(DEBUG)
      call SubEnd('ElemMatForm')
#endif

      return
      end
