      subroutine OrthoToIsoFric(nvmax,cinem,qa,qb,c,lTrans,keyopt)
c --- primary function: transform tangential stress and sliding to transform orthotropic friction
c                       condition to isotropic friction condition
c
c  input arguments:
c       nvmax  (int,sc,in)            - number of state variables
c
c  output arguments:
c       cinem  (dp,ar(35),inout)      - kinematic parameters used for contact
c       qa     (dp,ar(*),inout)       - stresses and state variables at initial
c                                       configuration
c       qb     (dp,ar(*),inout)       - stresses and state variables at updated
c                                       configuration for each integration point
c       c      (dp,ar(4,4),inout)     - material tensor
c
c
#include "impcom.inc"
c
      logical lTrans
      integer nvmax
      integer keyopt(12)
      double precision cinem(35),qa(nvmax),qb(nvmax),c(4,4)
      integer keyopt_1
      double precision
     &     scalefacr,scalefacs,taura,tausa,vtfra,vtfsa,lmra,
     &     lmsa,tvtfra,tvtfsa,avtfra,avtfsa,oldvtfra,oldvtfsa,
     &     taurb,tausb,vtfrb,vtfsb,lmrb,lmsb,tvtfrb,tvtfsb,
     &     avtfrb,avtfsb,cvtfr,cvtfs

#if defined(DEBUG)
      call SubBeg('OrthoToIsoFric')
#endif

c --- retrieve scale factors, tangential stress and sliding
c
      scalefacr = cinem(34)
      scalefacs = cinem(35)
      keyopt_1 = keyopt(1)
c
      taura = qa(2)
      tausa = qa(3)
      vtfra = qa(12)
      lmra = qa(22)
      lmsa = qa(23)
      tvtfra = qa(29)
      tvtfsa = qa(30)
      avtfra = qa(31)
      avtfsa = qa(32)
      oldvtfra = qa(nvmax-3)
      oldvtfsa = qa(nvmax-2)
c
      taurb = qb(2)
      tausb = qb(3)
      vtfrb = qb(12)
      vtfsb = qb(13)
      lmrb = qb(22)
      lmsb = qb(23)
      tvtfrb = qb(29)
      tvtfsb = qb(30)
      avtfrb = qb(31)
      avtfsb = qb(32)
      cvtfr = 0
      cvtfs = 0
      if (keyopt_1.ge.3) then
          cvtfr = cinem(28)
          cvtfs = cinem(29)
      endif

c --- apply scale factor
c
      if (lTrans) then
c --- scale tangential stress and sliding before friction computation
c
          taura = taura/scalefacr
          tausa = tausa/scalefacs
          vtfra = vtfra*scalefacr
          vtfsa = vtfra*scalefacs
          lmra = lmra*scalefacr  ! LM computed using tangential sliding and penalty
          lmsa = lmsa*scalefacs  ! so scale it like tangential sliding
          tvtfra = tvtfra*scalefacr
          tvtfsa = tvtfsa*scalefacs
          avtfra = avtfra*scalefacr
          avtfsa = avtfsa*scalefacs
          oldvtfra = oldvtfra*scalefacr
          oldvtfsa = oldvtfsa*scalefacs
c
          taurb = taurb/scalefacr
          tausb = tausb/scalefacs
          vtfrb = vtfrb*scalefacr
          vtfsb = vtfsb*scalefacs
          lmrb = lmrb*scalefacr  ! LM computed using tangential sliding and penalty
          lmsb = lmsb*scalefacs  ! so scale it like tangential sliding
          tvtfrb = tvtfrb*scalefacr
          tvtfsb = tvtfsb*scalefacs
          avtfrb = avtfrb*scalefacr
          avtfsb = avtfsb*scalefacs
      else

c --- scale tangential stress and sliding back after friction computation
c
          taura = taura*scalefacr
          tausa = tausa*scalefacs
          vtfra = vtfra/scalefacr
          vtfsa = vtfra/scalefacs
          lmra = lmra/scalefacr
          lmsa = lmsa/scalefacs
          tvtfra = tvtfra/scalefacr
          tvtfsa = tvtfsa/scalefacs
          avtfra = avtfra/scalefacr
          avtfsa = avtfsa/scalefacs
          oldvtfra = oldvtfra/scalefacr
          oldvtfsa = oldvtfsa/scalefacs
c
          taurb = taurb*scalefacr
          tausb = tausb*scalefacs
          vtfrb = vtfrb/scalefacr
          vtfsb = vtfsb/scalefacs
          lmrb = lmrb/scalefacr
          lmsb = lmsb/scalefacs
          tvtfrb = tvtfrb/scalefacr
          tvtfsb = tvtfsb/scalefacs
          avtfrb = avtfrb/scalefacr
          avtfsb = avtfsb/scalefacs
          if (keyopt_1.ge.3) then
              cvtfr = cvtfr/scalefacr
              cvtfs = cvtfs/scalefacs
          endif

c --- scale stiffness
c
          c(1,1) = c(1,1)*scalefacr*scalefacr
          c(2,2) = c(2,2)*scalefacs*scalefacs
          c(1,2) = c(1,2)*scalefacr*scalefacs
          c(2,1) = c(1,2)
          c(1,3) = c(1,3)*scalefacr
          c(2,3) = c(2,3)*scalefacs
      endif

c --- store tangential stress and sliding
c
      qa(2) = taura
      qa(3) = tausa
      qa(12) = vtfra
      qa(13) = vtfsa
      qa(22) = lmra
      qa(23) = lmsa
      qa(29) = tvtfra
      qa(30) = tvtfsa
      qa(31) = avtfra
      qa(32) = avtfsa
      qa(nvmax-3) = oldvtfra
      qa(nvmax-2) = oldvtfsa
c
      qb(2) = taurb
      qb(3) = tausb
      qb(12) = vtfrb
      qb(13) = vtfsb
      qb(22) = lmrb
      qb(23) = lmsb
      qb(29) = tvtfrb
      qb(30) = tvtfsb
      qb(31) = avtfrb
      qb(32) = avtfsb
      if (keyopt_1.ge.3) then
          cinem(28) = cvtfr
          cinem(29) = cvtfs
      endif
c

#if defined(DEBUG)
      call SubEnd('OrthoToIsoFric')
#endif

      return
      end
