      subroutine Sys_Lim (nrec,leng)
c --- primary function: check operating system file limit
c
c  input arguments:
c     leng     (int,sc,in)       - record length
c     nrec     (int,sc,in)       - record number
c
#include "impcom.inc"
c
      integer  nrec,leng,nmax,i,SplitMax
      parameter (SplitMax = 65536)  ! 256 GBytes

      nmax = (SplitMax*16384*64)/leng
      i = (nrec-1)/nmax
      if (i.gt.0) then
          call erhandler('Sys_Lim',5000,3,
     &                   'The maximum file size of %I GB reached'
     &                   ,dble(SplitMax/256),' ')
          call StopAll(-1)
      endif

      return
      end
