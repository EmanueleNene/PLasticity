      subroutine MaterialMain (elemId, IntId, matId, ireal,
     &                         nTens, nlkword, numWordsPt,
     &                         mRecIntPidx, kDomIntPt, kSectPt,
     &                         kMatPoint, moisture, strainInc,
     &                         lwork, work, lengthRec,
     &                         mref, khygro, coords, basis,
     &                         iniFlag, matRec, strain, stress,
     &                         matcJac, tMatProp, siglv, siglvi, 
     &                         eEnerPt
#if defined(MOIST)
     &                         ,eh)
#else
     &                          )
#endif      
c --- primary function: interface to elements
c
c  input arguments:
c      elemId       (int,sc,in)         - element number
c      matId        (int,sc,in)         - material number
c      nlkword      (intL,ar(*,2),in)   - number of words per point
c      numWordsPt   (int,sc,in)         - material record index
c      mRecIntPidx  (int,ar(*),in)      - pointers for material record
c                                         always starts at 1
c      kDomIntPt    (int,sc,in)         - "k"th domain int pt
c      kSectPt      (int,sc,in)         - "k"th section point
c      kMatPoint    (int,sc,in)         - "k"th material point of the element
c      moisture     (dp,sc,in)          - moisture at the current mat pt
c      strainInc    (dp,ar(nTens),in)   - incremental strain
c      lwork        (int,sc,in)         - length of work space
c      work         (dp,ar(lwork),in)   - work space
c      lengthRec    (int,sc,in)         - length of material record
c      mref         (dp,sc,in)          - reference moisture
c      coords       (dp,ar(3),in)       - coordinates of this point
c      basis        (dp,ar(3,3),in)     - basis directions at this point
c      iniFlag      (log,sc,in)         - stress initialization flag
c      grad         (dp,ar(3),in)       - hygro gradients at this point
c
c  output arguments:
c      matRec   (dp,ar(lengthRec),inout) - stress material point saved variables
c      khygro       (int,sc,out)         - moisture flag
c      strain       (dp,ar(nTens),out)   - strain at 't+Dt'
c      stress       (dp,ar(nTens),out)   - Cauchy stress at 't+Dt'
c                                          (these stress values are in the
c                                           global coordinate system, as they
c                                           are for load vector calculation)
c      matcJac      (dp,ar(nTens,*),out) - consistent Jacobian material matrix
c      tMatProp     (dp,ar(5),out)       - mat prop array for passing transverse shear
c      siglv        (dp,ar(4),out)       - hygro load vector
c      siglvi       (dp,ar(*),out)       - initial stress load vector
c      eEnerPt      (dp,ar(10),out)      - energies at this point
c
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "in_mem.inc"
#if defined(MOIST)
      double precision eh
#endif
c
      integer          elemId, IntId, matId, ireal,
     &                 mRecIntPidx(*),
     &                 kDomIntPt, kSectPt, kMatPoint,
     &                 lwork, lengthRec,
     &                 khygro, numWordsPt
      integer*8        nlkword(W_nMatCol,W_nMatRow)
      double precision moisture
      double precision
     &                 strainInc(6), strain(*),
     &                 work(lwork), matRec(lengthRec),
     &                 stress(*), matcJac(*),
     &                 tMatProp(5), mref, siglv(4), eEnerPt(10),
     &                 siglvi(*), coords(3), basis(3,3)
      logical          iniFlag

      integer          nTens, keycut, firstIndex, lastWkIdx,
     &                 matMWIndex  (W_LenMatWrkIndx),
     &                 matRecIndex (W_LenMatRecIndx),
     &                 matWorkIndex(W_LenWrkIndx),
     &                 i, j, lp(2)
      data lp /13,4/

      double precision dperr(2)


#if defined(DEBUG)
      call SubBeg('MaterialMain')
#endif

      keycut   = 0

      call v_zero (siglv(1),nTens)

c --- set up index
      lastWkIdx = lwork + 1 - W_FirstMatIndx
      call vi_zero (matMWIndex,   W_LenMatWrkIndx)
      call vi_zero (matRecIndex,  W_LenMatRecIndx)
      call vi_zero (matWorkIndex, W_LenWrkIndx)

c --- set pointers for deformation gradient workspaces
c     unpack the deformation gradient and pack it properly

      firstIndex = 1

c --- set material point variables indices

      i = (kMatPoint-1)*numWordsPt
      do j = 1, W_LenMatRecIndx
         matRecIndex(j) = i + mRecIntPidx(j)
      enddo

c --- set some often used work arrays
      call DefMatWrkIndx  (nlkword(1,1), nTens, i, matMWIndex(1))

      firstIndex = firstIndex + i

c --- work areas

      call DefTmpIndx (nTens, i, matWorkIndex)

      firstIndex = firstIndex + i

      if (firstIndex.gt.lastWkIdx) then
          dperr(1) = dble(firstIndex)
          dperr(2) = dble(lastWkIdx)
          call erhandler('MaterialMain',5000,3,
     &          'Memory shortage! Array must atleast
     &           be %i, current size is %i.', dperr(1),' ')
      endif

      call MaterialDriver (ireal, IntId, keycut, matRecIndex(1), 
     &                     matRec(1), matMWIndex(1), 
     &                     work(W_FirstMatIndx), matWorkIndex(1), 
     &                     work(W_FirstWrkIndx),
     &                     strainInc, moisture, nTens, elemId, matId,
     &                     kDomIntPt, kSectPt, kMatPoint,
     &                     mref, coords(1), basis(1,1),
     &                     iniFlag, tMatProp(1), khygro,
     &                     stress(1), matcJac(1), strain(1),
     &                     siglv(1), siglvi(1), eEnerPt(1),
     &                     nlkword(1,1)
#if defined(MOIST)
     &                     ,eh)
#else
     &                     )
#endif        

      if (A_KeyMatUpd.eq.1) then
          call UpdMatRecIndx (kMatPoint, nTens, khygro, nlkword(1,1),
     &                        strain(1), stress(1),
     &                        matMWIndex(1) , work(W_FirstMatIndx))

      endif

c
      if (keycut .gt. 0) then
          call PutErr(1)
          if (C_CutMat .ne. 1) then
              call lockSMP (W_LockCut)
              C_CutMat = 1
              call unlockSMP (W_LockCut)
          endif
      endif


#if defined(DEBUG)
      call SubEnd('MaterialMain')
#endif

      return
      end

