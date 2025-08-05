      subroutine ContManager (numel,lenbacl,numdofl,nodfwd,nodbac,
     &                        nPerSeg,nSeg,nPerFounI,nMast,nPerMast,
     &                        nPerFounDp,nFoun,icoor,dispn,isegda,
     &                        ififo,founda,icnda,nopilo,kerr)
c --- primary function: write foundation and segment database
c
c  input arguments:
c     numel     (int,sc,in)          - number of elements
c     lenbacl   (int,sc,in)          - number of active nodes
c     numdofl   (int,sc,in)          - number of active dof per node
c     nodfwd   (int,ar(*),in)        - forward list for compressed nodes
c     nodbac   (int,ar(*),out)       - internal to exter. list
c     nPerSeg    (int,sc,in)         - number of items per segment database entry
c     nSeg    (int,sc,in)            - number of segments on foundation
c     nPerMast   (int,sc,in)         - number of items per contact element
c     nPerFounI   (int,sc,in)        - number of integer items per foundation
c                                      database entry
c     nPerFounDp   (int,sc,in)       - number of double precision items per
c                                      foundation database entry
c     nMast   (int,sc,in)            - number of contact elements
c     nFoun    (int,sc,in)           - number of foundations
c     icoor    (dp,ar(*),in)         - initial nodal coordinate
c     dispn    (dp,ar(*),in)         - corrent nodal displacement
c
c  output arguments:
c     isegda   (int,ar(nPerSeg,nSeg),out)
c                                  - segment information
c          isegda(1,iseg)        = code of segment type (icode)
c                                = 1   : straight line
c                                = 2   : parabolic arc
c          isegda(2,iseg)        = fondation number
c          isegda(j,iseg)        = 4 nodes (j=3-10) on seg. iseg
c          isegda(j,iseg)        = 4 adjacent segments (j=11-14)
c                                = 2 adjacent segments (j=11,12)
c                                     for seg. iseg
c                                = 0  (icode > 10)
c                                  --> none adjacent segment
c          isegda(15,iseg)       = segment element number
c                                  < slave element death
c          isegda(16,iseg)       = attached beam element number, if it is
c          isegda(17,iseg)       = attached contact element number, if it is
c          isegda(18,iseg)       = remark
c                                  < slave element death
c
c     ififo    (dp,ar(nPerFounI,nFoun),out)
c                                  - foundation information (integer)
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
c                                -2 - mixed type surface (may wrong)
c                                -1 - rigid slave surface
c                                 0 - none
c                                 1 - deformable slave surface
c                                 2 - deformable slave surface (two pass)
c          ififo(16,ifoun)      = contact element number who is closest to current slave
c          ififo(17,ifoun)      = slave element number who is pair of closest contact
c                                 element
c          ififo(18,ifoun)      = iteration number of surface adjust
c          ififo(19,ifoun)      = max. node per slave element on each
c                                 slave surface

c     founda   (dp,ar(nperf,nFoun),out)
c                                  - foundation  cell region
c          founda(1,ifoun)       = minimum value in x-dir.
c          founda(2,ifoun)       = maximum value in x-dir.
c          founda(3,ifoun)       = minimum value in y-dir.
c          founda(4,ifoun)       = maximum value in y-dir.
c          founda(5,ifoun)       = minimum value in z-dir.
c          founda(6,ifoun)       = maximum value in z-dir.
c          founda(7,ifoun)       = first real constant
c          founda(8,ifoun)       = second real constant
c          founda(9,ifoun)       = gap or penetration of closest contact element
c          founda(10,ifoun)      = x component of founda(9,ifoun) (increment)
c          founda(11,ifoun)      = y component of founda(9,ifoun) (increment)
c          founda(12,ifoun)      = z component of founda(9,ifoun) (increment)
c          founda(13,ifoun)      = x component of founda(9,ifoun) (total)
c          founda(14,ifoun)      = y component of founda(9,ifoun) (total)
c          founda(15,ifoun)      = z component of founda(9,ifoun) (total)
c          founda(16,ifoun)      = thickness of founda(9,ifoun)
c          founda(17,ifoun)      = length of founda(9,ifoun)
c          founda(18,ifoun)      = stiffness of founda(9,ifoun)
c          founda(19,ifoun)      = number of founda(9,ifoun)
c          founda(20,ifoun)      = area of founda(9,ifoun)
c     icnda    (int,ar(nPerMast,nMast),out) - contact element information
c                               = 1 element type (ktyp)
c                               = 2 attached beam element number (ibeam)
c                               = 3 element in contact
c                               = 4 - 7 nodes number (nodes)
c                               = 8 contact element number (ielm)
c                               = 9 slave surface number (ifoun)
c     kerr     (int,sc,out)        - error = 1
c
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "dof_info.inc"
c
      integer numel,lenbacl,numdofl,nPerSeg,nSeg,nPerFounI,
     &        nPerFounDp,nFoun, nPerMast,nMast,GetErr
      integer nodfwd(lenbacl),isegda(nPerSeg,nSeg),
     &        ififo(nPerFounI,nFoun),icnda(nPerMast,nMast),
     &        nodbac(lenbacl),nopilo(lenbacl),kerr
      double precision founda(nPerFounDp,nFoun),
     &                 icoor(lenbacl*numdofl),dispn(lenbacl*numdofl)
      integer elData(EL_DIM),nodes(W_NdMax),ielm,nCon,ipara

#if defined(DEBUG)
      call SubBeg('ContManager')
#endif

      ipara = 0
c --- initialize contact status keys
      C_nPen = 0
      C_nConA = 0
      C_nConF = 0
      C_nDeBo = 0
      C_nFraBo = 0      
      C_nPenIni = 0
      C_ConnUpd = 0

      call vi_zero (isegda(1,1),nPerSeg*nSeg)
      call vi_zero (icnda(1,1),nPerMast*nMast)
      call vi_zero (ififo(1,1),nPerFounI*nFoun)
      call vi_zero (nopilo(1),lenbacl)
      call v_zero (founda(1,1),nPerFounDp*nFoun)

c --- find underlying elements for contact and slave elements
c
      call Cont_1 (numel,lenbacl,numdofl,nodfwd(1),nodbac(1),
     &             nSeg,nPerFounI,nMast,nPerMast,nFoun,
     &             ififo(1,1),icnda(1,1),nCon,ipara)

c --- form slave and contact surface database
c
      call Cont_2 (lenbacl,numdofl,nodfwd(1),nPerSeg,
     &             nSeg,nPerFounI,nMast,nPerFounDp,nFoun,
     &             icoor(1),dispn(1),isegda(1,1),ififo(1,1),founda(1,1))

      if (GetErr().eq.1) goto 999
      if (nCon.le.0 .or. nSeg.le.0) then
          kerr = 2
          goto 999
      endif

c --- set contact element type
c
      call Cont_3 (lenbacl,numdofl,nodbac(1),nPerFounI,nMast,nPerMast,
     &             nFoun,ififo(1,1),icnda(1,1),nCon)

      if (GetErr().eq.1) goto 999

c --- find the adjacent segments for each slave element
c
      call Cont_4 (numel,lenbacl,numdofl,nodbac(1),nPerSeg,
     &             nSeg,nPerFounI,nPerFounDp,nFoun,
     &             isegda(1,1),ififo(1,1),founda(1,1),icnda(1,1),
     &             nopilo(1),nCon,nPerMast,nMast)

      if (GetErr().eq.1) goto 999

c --- get initial gap and penetration information
c
      call Cont_5 (numel,lenbacl,numdofl,nodfwd(1),nodbac(1),nPerSeg,
     &             nSeg,nPerFounI,nMast,nPerMast,nPerFounDp,nFoun,
     &             icoor(1),dispn(1),isegda(1,1),ififo(1,1),founda(1,1),
     &             icnda(1,1),nCon)

  999 continue

#if defined(DEBUG)
      call SubEnd('ContManager')
#endif

      return
      end
   