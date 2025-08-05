      subroutine CalTransMat(xyz,qlc,tr,alen)
c --- primary function: compute the new ek matrix transpose for beam element
c
c     fast calculation of the new (updated) ek matrix transpose
c     vector e1 lies along element axis i-j, vectors e2 & e3
c     ref: eqns. (36) and (37) of 
c          crisfield, 'a consistent co-rotational formulation for non-linear, 
c          three-dimensional, beam-element', CMAME, 81 (1990) 131 - 150
c
c  input arguments:
c     tr    (dp,ar(3,3),in)  - mid-point triad matrix from nodal triads
c     xyz   (dp,ar(3,2),in)  - updated global coordinates
c
c  output arguments:
c     tr    (dp,ar(3,3),out) - ek matrix transpose; the new element transform matrix
c     alen  (dp,sc,out)      - element length in updated coordinate system
c     qlc   (dp,ar(3),out)   - reference vector which records the element
c                              rotations using mid-triad matrix
c
#include "impcom.inc"
c
      double precision xyz(3,2),qlc(3),tr(3,3),e1(3),e2(3),e3(3),
     &                 con,con2,alen,v_dot,r1(3),r2(3),r3(3),con3
c
      call v1mv2_v3 (xyz(1,2),xyz(1,1),e1(1),3)
      con = v_dot (e1(1),e1(1),3)
      alen = sqrt(con)
      con = 1.0d0/alen
      call cv1_v1 (e1(1),3,con)
      call v_movei (tr(2,1),3,r2(1),1,3)
      call v_movei (tr(3,1),3,r3(1),1,3)
      call v_movei (tr(1,1),3,r1(1),1,3)
      con2 = v_dot (r2(1),e1(1),3)/2.0d0
      con3 = v_dot (r3(1),e1(1),3)/2.0d0
      call v1pv2_v3 (e1(1),r1(1),qlc(1),3)
      call cv1_v2 (qlc(1),r1(1),3,con2)
      call v1mv2_v3 (r2(1),r1(1),e2(1),3)
      call cv1_v2 (qlc(1),r1(1),3,con3)
      call v1mv2_v3 (r3(1),r1(1),e3(1),3)
      call v_move (e2(1),qlc(1),3)
      call v_movei (e1(1),1,tr(1,1),3,3)
      call v_movei (e2(1),1,tr(2,1),3,3)
      call v_movei (e3(1),1,tr(3,1),3,3)
c
      return
      end
