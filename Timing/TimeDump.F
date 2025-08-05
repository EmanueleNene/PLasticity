      subroutine TimeDump
c --- primary function: dump timing info
c
#include "impcom.inc"
#include "time_info.inc"
#include "constants.inc"
#include "file_sys.inc" 
c 
      integer   Sys_Open, Sys_Close
      integer   i, j, k, nunit, n, nk, NameMax, indx(T_MaxSubNum)
      character*(B_PathLen)  pname 
      character*32      AllNames(T_MaxSubNum)
      double precision  AllTimes(T_MaxSubNum), TotalTime, fract
      integer*8         AllCount(T_MaxSubNum)   
      
      T_Stat = 0  ! stop timing
      
      call MakeName (B_JobName,'tim     ',B_Dir,pname,i)

      nunit = 37

      i = Sys_Open (nunit,pname,'UN','SERIAL','FORM','WRITE','BEGIN') 

      call v_zero (AllTimes(1),T_MaxSubNum)

      do i = 1,T_MaxSubNum
         AllCount(i) = 0
      enddo

      do i = 1,T_MaxSubNum
         AllNames(i) = ''
      enddo   
      
      NameMax = 0
      TotalTime = 0.0
      do i = 1,T_MaxLev
         n = T_NumPerLev(i)
         do j = 1,n
            do k = NameMax,1,-1       ! search in reverse order
               if (T_SubName(j,i).eq.AllNames(k)) then
                   nk = k
                   goto 20
               endif
            enddo
            NameMax = NameMax + 1
            if (NameMax.gt.T_MaxSubNum) goto 30
            nk = NameMax
            AllNames(nk) = T_SubName(j,i)
   20       AllTimes(nk) = AllTimes(nk) + T_SubTime(j,i)
            AllCount(nk) = AllCount(nk) + T_NumCalls(j,i)
            TotalTime = TotalTime + T_SubTime(j,i)
         enddo
      enddo

      if (TotalTime .lt. TINY) goto 998

  30  write (nunit,*)
      write (nunit,*) '--------------------> Performance data'
      write (nunit,*) '1 - Subroutine name'
      write (nunit,*) '2 - Total number of calls'
      write (nunit,*) '3 - Total time consumed'
      write (nunit,*) '4 - Procentage from total'
      write (nunit,*)

c --- sort upon the time consumed      
      call v_sort(AllTimes(1),NameMax,indx(1))
      fract = 100.0d0/TotalTime
      do nk = 1,NameMax
         i = indx(nk)
         write (nunit,2000) 
     &   AllNames(i),AllCount(i),AllTimes(nk),AllTimes(nk)*fract
 2000    format (1x,a32,1x,i15,g13.6,g13.3)
      enddo
      
 998  write (nunit,3000) TotalTime
3000  format (/'TOTAL CPU TIME',g13.6/)
      
      if (NameMax. gt. T_MaxSubNum)
     &    write (nunit,*) 'Note, number of routines exceeded limit of ',
     &                     T_MaxSubNum
      
      i = Sys_Close(nunit,'KEEP') 

      return
      end



     
     
