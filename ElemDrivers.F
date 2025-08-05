      subroutine ElemDrivers (elem,elChar,elData,eomask,nodes,SVRfp,
     &                        kelreq,nr,xyzang,uelm,kelout,
     &                        Kmtrx,Mmtrx,Cmtrx,gStif,
     &                        Fnr,elmass,center,elvol,elener,
     &                        edindxL,lcerstL,ndUpdate)
c --- primary function: call appropriate element routine
c
c  input arguments:
c       elem   (int,sc,in)         - element number
c       elChar (int,ar(*),in)      - array of element type characteristics
c       elData (int,ar(EL_DIM),in) - array of element data
c       eomask (int,sc,in)         - bit pattern for element output
c       nodes  (int,ar(*),in)      - array of element node numbers
c       SVRfp  (int*8,sc,in)       - location of this elem SVR's on file SVR
c       kelreq (int,ar(10),in)     - matrix and load vector form requests
c       nr     (int,sc,in)         - matrix and lv size
c       xyzang (dp,ar(6,*),in)     - nodal coords (orig) and rotation angles
c       uelm   (dp,ar(nr,*),in)    - element nodal solution values
c
c  output arguments:
c       kelout (int,ar(10),out)    - keys indicating created mats and LV's
c       Kmtrx  (dp,ar(nr,nr),inout)- k matrix
c       Mmtrx  (dp,ar(nr,nr),inout)- m matrix
c       gStif  (dp,ar(nr,nr),inout)- s matrix
c       Fnr    (dp,ar(nr,2),out)   - load vectors:
c                                     Fnr(1,1) = applied f vector
c                                     Fnr(1,2) = n-r restoring vector
c       elmass (dp,sc,out)         - element mass
c       center (dp,ar(3),out)      - centroid location
c       elvol  (dp,sc,out)         - element volume
c       elener (dp,ar(10),out)     - element energies
c       edindxL(int*8,ar(25),out)  - index to element result file data
c       lcerstL(int*8,sc,inout)    - result file pointer
c       ndUpdate(int,sc,out)       - keys indicating the update nodes
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer*8  SVRfp,lcerstL,edindxL(25)
      integer elem,elChar(*),elData(EL_DIM),nodes(*),kelreq(10),
     & nr,kelout(10),jtyp,eomask,ndUpdate
      double precision
     & xyzang(6,*),uelm(nr,*),Kmtrx(nr,nr),Mmtrx(nr,nr),Cmtrx(nr,nr),
     & gStif(nr,nr),Fnr(nr,2),elmass,center(3),elvol,elener(10)

#if defined(DEBUG)
      call SubBeg('ElemDrivers')
#endif

      jtyp = elChar(JETYP)
      ndUpdate = 0
      if (jtyp.eq.W_eMaster) then
          call ContactElemDriver (elem,elChar(1),elData(1),eomask,
     &                            nodes(1),SVRfp,kelreq(1),nr,
     &                            xyzang(1,1),uelm(1,1),kelout(1),
     &                            Kmtrx(1,1),Fnr(1,1),Fnr(1,2),
     &                            elvol,center(1),elener(1),
     &                            edindxL(1),lcerstL,ndUpdate)
      elseif (jtyp.eq.W_eBeam3nd .or. jtyp.eq.W_eBeam2nd) then
          call BeamElemDriver (elem,elChar(1),elData(1),eomask,nodes(1),
     &                         SVRfp,kelreq(1),nr,xyzang(1,1),
     &                         uelm(1,1),kelout(1),
     &                         Kmtrx(1,1),Mmtrx(1,1),gStif(1,1),
     &                         Fnr(1,1),Fnr(1,2),elvol,
     &                         elmass,center(1),elener(1),
     &                         edindxL(1),lcerstL)
      elseif (jtyp.eq.W_eLink) then
         call LinkElemDriver (elem,elChar(1),elData(1),eomask,
     &                        nodes(1),SVRfp,kelreq(1),nr,
     &                        xyzang(1,1),uelm(1,1),kelout(1),
     &                        Kmtrx(1,1),Mmtrx(1,1),gStif(1,1),
     &                        Fnr(1,1),Fnr(1,2),elvol,
     &                        elmass,center(1),elener(1),
     &                        edindxL(1),lcerstL)
      elseif (jtyp.eq.W_eBeamEB) then
         call EBeam3D (elem,elChar(1),elData(1),eomask,
     &                 nodes(1),SVRfp,kelreq(1),nr,
     &                 xyzang(1,1),uelm(1,1),kelout(1),
     &                 Kmtrx(1,1),Mmtrx(1,1),gStif(1,1),
     &                 Fnr(1,1),Fnr(1,2),elvol,
     &                 elmass,center(1),elener(1),
     &                 edindxL(1),lcerstL)
      endif

#if defined(DEBUG)
      call SubEnd('ElemDrivers')
#endif

      return
      end
