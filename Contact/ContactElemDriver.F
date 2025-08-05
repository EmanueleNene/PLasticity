      subroutine ContactElemDriver (elem,elChar,elData,eomask,nodes,
     &                              locsvrL,kelrqq,nr,xyz,u,kelout,
     &                              Kmtrx,Fext,Fnr,elvol,center,elener,
     &                              edindxL,lcerstL,ContStat)
c --- primary function: compute element matrices, load vectors, and results
c
c  input arguments:
c       elem   (int,sc,in)         - element number
c       elChar (int,ar(IELCSZ),in) - array of element type characteristics
c       elData (int,ar(EL_DIM),in) - array of element data
c       nodes  (int,ar(nnod),in)   - array of element node numbers
c       locsvrL(integer*8,sc,in)   - location of this eles SVRs on file SVR
c       kelrqq (int,ar(10),in)     - matrix and load vector form requests
c       xyz    (dp,ar(6,nnod),in)  - nodal coords (orig) and rotation angles
c       nr     (int,sc,in)         - matrix and LV size
c       u      (dp,ar(nr,5),in)    - element nodal solution values
c
c  output arguments:
c       kelout (int,ar(10),out)    - keys indicating created mats and LVs
c       Kmtrx  (dp,ar(nr,nr),inout)- k matrix
c       Fnr    (dp,ar(nr),out)     - n-r restoring f vector
c       elvol  (dp,sc,out)         - element volume
c       center (dp,ar(3),out)      - centroid location
c       elener (dp,ar(5),out)      - element energies
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "cont_info.inc"
c
      integer*8  locsvrL,edindxL(25),lcerstL
      integer    elem,elChar(IELCSZ),elData(EL_DIM),nodes(*),
     &           kelrqq(10),nr,kelout(10),eomask

      integer    nssvr,nvmax,ntmax,svindx(20),ival(3),npi, ContStat


      double precision
     &           xyz(6,*),u(nr,5),Kmtrx(nr,nr),Fext(nr),Fnr(nr),
     &           elener(*),center(3),elvol

      double precision ssvr(372)

#if defined(DEBUG)
      call SubBeg('ContactElemDriver')
#endif

      ContStat = 0
c --- get information for this specific master element
      call conget (elem,ival(1))
      nssvr = elChar(NMSSVR)
      nvmax = elChar(NMCNSVR)
c      write (8,*) '1. nvmax, nssvr',nvmax,nssvr,elem
c --- get the SVR index vector
      call GetSVRIndx (elem,locsvrL,svindx(1))
c      write (8,*) '1.1 nvmax, nssvr',nvmax,nssvr,svindx(1)
c --- get the element SVRs
      call GetSVR (svindx(1),1,nssvr,ssvr(1))

c      write (8,*) '2. nvmax, nssvr',nvmax,nssvr

      if (elChar(KYO10).eq.0) then
c --- beam-beam contact
          npi = 1 ! number of effective int. point
          ntmax = nvmax*npi
          call ContElement (elem,elChar(1),elData(1),eomask,
     &                      nodes(1),kelrqq(1),nr,nvmax,npi,
     &                      ival(1),K_nSeg,K_nFoun,K_nPerSeg,
     &                      K_nPerFounI,K_nPerFounDp,K_nPerMast,K_nMast,
     &                      K_Seg(1),K_FounI(1),K_FounDp(1),
     &                      K_ContInfo(1),xyz(1,1),u(1,1),K_NdCoord(1),
     &                      M_BackList(1),ssvr(1),ssvr(ntmax+1),
     &                      ssvr(2*ntmax+1),kelout(1),
     &                      Kmtrx(1,1),Fext(1),Fnr(1),
     &                      elvol,center(1),elener(1),
     &                      ContStat,edindxL(1),lcerstL)
      else !if  (elChar(KYO10).eq.1 .or. elChar(KYO10).eq.2) then
c --- beam-surface contact and parallel beam contact
          npi = 3 ! number of effective int. point
          ntmax = nvmax*npi
          call ContElement4 (elem,elChar(1),elData(1),eomask,
     &                       nodes(1),kelrqq(1),nr,nvmax,npi,
     &                       ival(1),K_nSeg,K_nFoun,K_nPerSeg,
     &                       K_nPerFounI,K_nPerFounDp,K_nPerMast,
     &                       K_nMast,K_Seg(1),K_FounI(1),K_FounDp(1),
     &                       K_ContInfo(1),xyz(1,1),u(1,1),
     &                       M_SoluVects(1),K_NdCoord(1),M_BackList(1),
     &                       ssvr(1),ssvr(ntmax+1),ssvr(2*ntmax+1),
     &                       kelout(1),Kmtrx(1,1),Fnr(1),elvol,
     &                       center(1),elener(1),ContStat,
     &                       edindxL(1),lcerstL)
      endif
c      write (8,*) '3. nvmax, nssvr',nvmax, nssvr
c --- write out the element saved variables
c --- base variables
      call StoreSVR (svindx(1),1,nssvr,ssvr(1))
c --- write out the SVR index vector
      call PutSVRIndx (elem,locsvrL,svindx(1))
      
      if (ContStat.ne.0) then
          if (K_ContStat.eq.0) then
              K_ContStat = max(ContStat,K_ContStat)
          endif
      endif
c      write (8,*) '4. ntmax, nvmax, nssvr',ntmax, nvmax, nssvr
      
#if defined(DEBUG)
      call SubEnd('ContactElemDriver')
#endif

      return
      end