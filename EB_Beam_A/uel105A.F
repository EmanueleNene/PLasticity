c      subroutine uel105 
      subroutine el04
     &                  (elem,ielc,elmdat,eomask,nodes,locsvrL,kelreq,
     &                   kelfil,nr,xyz,u,kelout,zs,zass,damp,gstif,
     &                   zsc,zscnr,elvol,elmass,center,elener,
     &                   edindxL,lcerstL)
c --- primary function:
c        1. compute element matrices, load vectors, and results
c --- secondary functions:
c        2. maintain element solution data
c
c  input arguments:
c       elem   (int,sc,in)          - element label (number)
c       ielc   (int,ar(IELCSZ),in)  - array of element type characteristics
c       elmdat (int,ar(EL_DIM),in)  - array of element data
c       eomask (int,sc,in)          - bit pattern for element output
c       nodes  (int,ar(nnod),in)    - array of element node numbers
c       locsvrL (LONGINT,sc,in)     - location of this eles svrs on file esav
c       kelreq (int,ar(10),in)      - matrix and load vector form requests
c       kelfil (int,ar(10),in)      - keys indicating incoming mats and lvs
c       nr     (int,sc,in)          - matrix and lv size
c       xyz    (dp,ar(6,nnod),in)   - nodal coords (orig) and rotation angles
c       u      (dp,ar(nr,5),in)     - element nodal solution values
c
c  output arguments:
c       kelout (int,ar(10),out)     - keys indicating created mats and lvs
c       zs     (dp,ar(nr,nr),inout) - k matrix (kelreq(1))
c       zass   (dp,ar(nr,nr),inout) - m matrix (kelreq(2))
c       damp   (dp,ar(nr,nr),inout) - c matrix (kelreq(3))
c       gstif  (dp,ar(nr,nr),inout) - s matrix (kelreq(4))
c       zsc    (dp,ar(nr),out)      - applied f vector (kelreq(5))
c       zscnr  (dp,ar(nr),out)      - n-r restoring f vector (kelreq(6))
c       elvol  (dp,sc,out)          - element volume
c       elmass (dp,sc,out)          - element mass
c       center (dp,ar(3),out)       - centroid location
c       elener (dp,ar(5),out)       - element energies
c       edindxL(LONG,ar(25),out)    - element result data file indexes
c       lcerstL(LONG,sc,inout)      - position on result file
c
#include "impcom.inc"
#include "ansysdef.inc"
#include "echprm.inc"
c
      LONGINT  locsvrL, lcerstL, edindxL(25)

      integer
     & elem,ielc(IELCSZ),elmdat(EL_DIM),eomask,nodes(3),kelreq(10),
     & kelfil(10),nr,kelout(10)

      double precision
     & xyz(6,3),u(nr,5),zs(nr,nr),zass(nr,nr),damp(nr,nr),gstif(nr,nr),
     & zsc(nr),zscnr(nr),elvol,elmass,center(3),elener(5)


      if (ielc(KYOP1).eq.1) then
         call EBeam2DA (elem,ielc(1),elmdat(1),eomask,nodes(1),locsvrL,
     &                 kelreq(1),kelfil(1),nr,xyz(1,1),u(1,1),kelout(1),
     &                 zs(1,1),zass(1,1),gstif(1,1),zsc(1),zscnr(1),
     &                 elvol,elmass,center(1),elener(1),edindxL,lcerstL)
      else
         call EBeam3DA (elem,ielc(1),elmdat(1),eomask,nodes(1),locsvrL,
     &                 kelreq(1),kelfil(1),nr,xyz(1,1),u(1,1),kelout(1),
     &                 zs(1,1),zass(1,1),gstif(1,1),zsc(1),zscnr(1),
     &                 elvol,elmass,center(1),elener(1),edindxL,lcerstL)
      endif

      return
      end
