      subroutine LinkElemDriver (elemID,elemChar,elData,eomask,eConn,
     &                           locsvrL,kelreq,nUsrDof,xBeg,eleDisp,
     &                           kelout,eStiff,eMass,gStiff,fExt,
     &                           fInt,elemVolume,elemMass,elemCG,
     &                           eEnergy,edindxL,lcerstL)
c --- primary function: prepare input for beam driver
c     
c  input arguments:
c     elemID  (int,sc,in)             - element number
c     elemChar(int,ar(*),in)          - array of element type characteristics
c     elData  (int,ar(EL_DIM),in)     - array of element data
c     eomask  (int,sc,in)             - bit pattern of element output
c     eConn   (int,ar(2),in)          - array of element node numbers
c     locsvrL (int*8,sc,in)           - location of this eles SVRs
c     kelreq  (int,ar(10),in)         - matrix and load vector form requests
c     nUsrDof (int,sc,in)             - matrix and LV size
c     xBeg    (dp,ar(6,2),in)         - nodal coords and rot angles
c     eleDisp (dp,ar(nUsrDof,5),in)   - element nodal solution values
c     
c  output arguments:
c     kelout (int,ar(10),out)         - keys indicating created mats and LVs
c     eStiff (dp,ar(nUsrDof,nUsrDof),inout) - k matrix (kelreq(1))
c     eMass  (dp,ar(nUsrDof,nUsrDof),inout) - m matrix (kelreq(2))
c     gStiff (dp,ar(nUsrDof,nUsrDof),inout) - s matrix (kelreq(4))
c     fExt   (dp,ar(nUsrDof),out)     - applied f vector (kelreq(5))
c     fInt   (dp,ar(nUsrDof),out)     - n-r restoring f vector (kelreq(6))
c     elemVolume  (dp,sc,out)         - element volume
c     elemMass (dp,sc,out)            - element mass
c     elemCG (dp,ar(3),out)           - centroid location
c     eEnergy(dp,ar(5),out)           - element energies
c     edindxL(int*8,ar(25),out)       - element result data file indexes
c     lcerstL(int*8,sc,inout)         - position on result file
c     
#include "impcom.inc"
#include "constants.inc"
c
      integer*8  locsvrL, lcerstL, edindxL(25)
      integer
     &     elemID,elemChar(*),elData(EL_DIM),eomask,eConn(2),
     &     kelreq(10),nUsrDof,kelout(10)

      double precision
     &     xBeg(6,2),eleDisp(nUsrDof,8),elemCG(3),eEnergy(5),
     &     eStiff(nUsrDof,nUsrDof),eMass(nUsrDof,nUsrDof),
     &     gStiff(nUsrDof,nUsrDof),fExt(nUsrDof),fInt(nUsrDof),
     &     elemVolume,elemMass

#if defined(DEBUG)
      call SubBeg('LinkElemDriver')
#endif

      call LinkElement (elemID,elemChar(1),elData(1),eomask,eConn(1),
     &                  locsvrL,kelreq(1),nUsrDof,xBeg(1,1),
     &                  eleDisp(1,1),kelout(1),eStiff(1,1),eMass(1,1),
     &                  gStiff(1,1),fExt(1),fInt(1),elemVolume,elemMass,
     &                  elemCG(1),eEnergy(1),edindxL(1),lcerstL)

#if defined(DEBUG)
      call SubEnd('LinkElemDriver')
#endif

      return
      end
