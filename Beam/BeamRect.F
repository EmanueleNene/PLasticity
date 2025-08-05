      subroutine BeamRect(xySection, np_x, np_y, rData,
     &                    ElemConn, numNodes, numElems)
c
#include "impcom.inc"
#include "ansysdef.inc"
      double precision 
     &     rData(*),xySection(2,4),ZERO, HALF,g,h, ONE, FOURTH,
     &     a1,a2,a3,a4, dg,dh, dgm, dhm
      integer ElemConn(9,*),numElems, numNodes,
     &     k1, k2, np_x, np_y, knod,i,j, kel, n0,n1,n2,nelx,nely

      parameter (ZERO=0.0d0, HALF=0.5d0, ONE=1.0d0, FOURTH=0.25d0)

      double precision vg(7), vh(7)

      numNodes = np_x * np_y
      numElems = (np_x-1)/2 * (np_y-1)/2
      nelx = (np_x-1)/2
      nely = (np_y-1)/2
      dg   = ONE / dble (nelx)
      dh   = ONE / dble (nely)
      call v_zero(vg,np_x)
      call v_zero(vh,np_y)
      vg(1) = -ONE
      if (nelx.ge.3) then
          vg(2) = vg(1) + 0.125*dg
          vg(3) = vg(2) + 0.125*dg
          dgm   = (dg*nelx-0.25*dg)/dble(nelx-2)
      else
          vg(2) = vg(1) + dg
          vg(3) = vg(2) + dg
          vg(4) = vg(3) + dg
          vg(5) = ONE
      endif
      if (nelx.ge.3) then
          k2 = (nelx-1)*2+1
          do k1 = 4, k2
             vg(k1) = vg(k1-1)+dgm
          enddo
          vg(k2+1) = vg(k2) + 0.125*dg
          vg(k2+2) = ONE
      endif
      vh(1) = -ONE
      if (nely.ge.3) then
          vh(2) = vh(1) + 0.125*dh
          vh(3) = vh(2) + 0.125*dh
          dhm   = (dh*nely-0.25*dh)/dble(nely-2)
      else
          vh(2) = vh(1) + dh
          vh(3) = vh(2) + dh
          vh(4) = vh(3) + dh
          vh(5) = ONE
      endif
      if (nely.ge.3) then
          k2 = (nely-1)*2+1
          do k1 = 4, k2
             vh(k1) = vh(k1-1)+dhm
          enddo
          vh(k2+1) = vh(k2) + 0.125*dh
          vh(k2+2) = ONE
      endif
      knod = 0
      do k2=1,np_y
         h=vh(k2)
         do k1=1,np_x
            g=vg(k1)
            a1=FOURTH*(ONE-g)*(ONE-h)
            a2=FOURTH*(ONE+g)*(ONE-h)
            a3=FOURTH*(ONE+g)*(ONE+h)
            a4=FOURTH*(ONE-g)*(ONE+h)
            knod = knod+1
            rData(knod)=
     &           a1*xySection(1,1)+a2*xySection(1,2)+
     &           a3*xySection(1,3)+a4*xySection(1,4)
            rData(knod+numNodes)=
     &           a1*xySection(2,1)+a2*xySection(2,2)+
     &           a3*xySection(2,3)+a4*xySection(2,4)
         enddo
      enddo
      kel=0
      do j = 1,(np_y-1)/2
         do i = 1,(np_x-1)/2
            kel=kel+1
            n0=2*np_x*(j-1)+2*(i-1)
            n1=n0+np_x
            n2=n1+np_x
            ElemConn(1,kel)=n0+1
            ElemConn(2,kel)=n0+3
            ElemConn(3,kel)=n2+3
            ElemConn(4,kel)=n2+1
            ElemConn(5,kel)=n0+2
            ElemConn(6,kel)=n1+3
            ElemConn(7,kel)=n2+2
            ElemConn(8,kel)=n1+1
            ElemConn(9,kel)=n1+2
         enddo
      enddo
      do k1=1,numElems
         rData(2*numNodes+k1) = ONE
      enddo

      return
      END
