      subroutine Cont_3 (BacLen,NumDof,nodbac,nPerFounI,nMast,nPerMast,
     &                   nFoun,ififo,icnda,ncon)
c --- primary function: set contact element
c
c  input arguments:
c     BacLen   (int,sc,in)     - number of active nodes
c     nPerMast   (int,sc,in)   - number of items per contact element
c     nPerFounI   (int,sc,in)  - number of integer items per foundation
c     nMast   (int,sc,in)      - number of contact elements
c     nFoun    (int,sc,in)     - number of foundations
c
c  output arguments:
c     ififo    (dp,ar(nPerFounI,nFoun),out) - foundation information (integer)
c          ififo(1,ifoun)       = number of segments
c          ififo(4,ifoun)       = constrain condition in x-dir.
c          ififo(5,ifoun)       = constrain condition in y-dir.
c          ififo(6,ifoun)       = constrain condition in z-dir.
c          ififo(7,ifoun)       = dimensions
c          ififo(8,ifoun)       = 1st segment numbering on this foundation
c          ififo(9,ifoun)       = 2nd segment numbering on this foundation
c          ififo(10,ifoun)      = lowest nodal numbering on this foundation
c          ififo(12,ifoun)      = real constant number (ireal)
c          ififo(13,ifoun)      = contact pair key
c                                 0 - no contact pair
c                                 1 - only slave element is defined
c                                 2 - only contact element is defined
c                                 3 - contact pair is defined
c          ififo(14,ifoun)      = type of associated contant element type (ityp)
c          ififo(15,ifoun)      = type of deformable slave surface
c                                -3 - rigid slave surface no node fix
c                                -2 - mixed type surface (may wrong)
c                                -1 - rigid slave surface
c                                 0 - none
c                                 1 - deformable slave surface
c                                 2 - deformable slave surface (two pass)
c          ififo(16,ifoun)      = contact element number who is colsest to current slave
c          ififo(17,ifoun)      = slave element number who is pair of closest contact
c                                 element
c          ififo(18,ifoun)      = iteration number of surface adjust
c          ififo(19,ifoun)      = max. node per slave element on each slave surface
c     icnda    (int,ar(nPerMast,nMast),out) - contact element information
c                               = 1 element type (ktyp)
c                               = 2 attached beam element number (ibeam)
c                               = 3 element in contact
c                               = 4 - 7 nodes number (nodes)
c                               = 8 contact element number (ielm)
c                               = 9 slave surface number (ifoun)
c     ncon     (int,sc,in)         - active number of contact elements
c
c
#include "impcom.inc"
#include "constants.inc"
c
      integer  elmget,etyiqr,etyget

      integer nPerFounI,nFoun,nPerMast,nMast,BacLen,NumDof,
     &  nodbac(BacLen),ififo(nPerFounI,nFoun),icnda(nPerMast,nMast),
     &  nnod,elData(EL_DIM),nodes(W_NdMax),ival(3),ifoun,ielm,
     &  ietypz,icon,ncon,npi,ityp,elChar(IELCSZ),i,
     &  inopil,iropil

      double precision dperr(3)

#if defined(DEBUG)
      call SubBeg('Cont_3')
#endif

      do 100 icon = 1,ncon
         ielm = icnda(8,icon)
 
         call vi_zero (nodes(1),W_NdMax)
         nnod = elmget(ielm,elData(1),nodes(1))

         ityp = elData(EL_TYPE)
         ietypz = etyiqr(ityp,-JETYP)
c --- contact information
         i = etyget(ityp,elChar(1))
         call conget (ielm,ival(1))
         npi  = abs(ival(1))
         ifoun = icnda(9,icon)
          
         if (ival(2).ne.icon) call erhandler ('Cont_3',5000,3,
     &                       'Contact specs are corrupted.', 0.0d0,' ')
c
         ival(1) = npi
         if (icnda(1,icon).eq.2) ival(1) = -abs(npi)

         if (ififo(15,ifoun).eq.-1.or.ififo(15,ifoun).le.-3) then
c --- rigid-deformable contact
             nnod = 3
             inopil = (ififo(2,ifoun)-1)/NumDof + 1
             nnod = nnod + 1
             if (inopil.gt.0) then
                 nodes(nnod) = nodbac(inopil)
                 if (ififo(3,ifoun).eq.0) then
                     iropil = 0
                 else
                     iropil = (abs(ififo(3,ifoun))-1)/NumDof + 1
                 endif
                 if (iropil.gt.0) then
                     nnod = nnod + 1
                     nodes(nnod) = nodbac(iropil)
                 endif
             else
                 nodes(nnod) = nodbac(ififo(10,ifoun))
             endif
         elseif (ififo(15,ifoun).gt.0) then
c --- deformable-deformable contact
             nnod = 7
             if (ietypz.eq.W_eMaster) then
                 ififo(19,ifoun) = 7
             endif

             nnod = nnod + npi*ififo(19,ifoun)

             if (nnod.gt.W_NdMax) then
                 if (ififo(18,ifoun).ne.3) then
                     ififo(18,ifoun) = 3
                     dperr(1) = nnod
                     dperr(2) = ityp
                     dperr(3) = W_NdMax
                     call erhandler ('Cont_3',5000,2,'For this
     &               contact, %i nodes per element with type %i
     &               may be required, which exceeds W_NdMax %i.',
     &               dperr(1),' ')
                 endif
                 nnod = W_NdMax
             endif
         endif

         elChar(NMNDST) = nnod

c --- store the local DOFs for this element type
         call etyput (ityp,IELCSZ,elChar(1))
         call elmpnd (ielm,nnod,nodes(1))

c --- put new contact element attributes
         call conput (ielm,ival(1))

c --- verify the foundation type
         if (ifoun.le.0) goto 100
         ififo(7,ifoun) = 3
c
 100  continue

#if defined(DEBUG)
      call SubEnd('Cont_3')
#endif

      return
      end
      
