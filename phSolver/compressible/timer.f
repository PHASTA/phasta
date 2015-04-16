        subroutine timer (tcode)
c
c----------------------------------------------------------------------
c
c This routine keeps track of the CPU-time statistics.
c
c input:
c  tcode        : timer codes
c                  name of processes to be timed
c                  including 'Begin   ', 'End     ','Back    '
c                  which are the control input
c
c
c Flop count:
c
c  x + y        :  1 flop    ->  same as Cray Y-MP
c  x - y        :  1 flop    ->  same as Cray Y-MP
c  x * y        :  1 flop    ->  same as Cray Y-MP
c  x / y        :  4 flops   ->  same as Cray Y-MP
c  sqrt(x)      :  8 flops   ->  lower than Cray Y-MP
c  exp(x)       : 24 flops   ->  same as Cray Y-MP
c  log(x)       : 32 flops   ->  same as Cray Y-MP
c
c
c
c **** WARNING: this routine makes calls to the SECS function.
c
c Farzin Shakib, Summer 1985.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        character*8 tcode
c
        dimension   ratio(11)
c
c.... Find the code
c
        do i = 1, icd + 2
          ii = i
          if (tcode .eq. ccode(i)) exit
          if (ii .eq. icd + 2) call error ('timer   ', tcode, 0)
        enddo
c
c.... ------------------------>  'Begin   '  <-------------------------
c
c.... initialize the timers and the flop counter
c
        if (ii .eq. icd) then 
          comtim  = zero
          do i = 1, icd
            cpu0(i)   = zero
            cpu(i)    = zero
            nacess(i) = 0
          enddo
c
          flops  = 0
          gbytes = 0
          sbytes = 0
c
          cpu0(icd)   = secs()
          nacess(icd) = 1
          call system_clock (iclock)
          return
        endif
c
c.... ------------------------>  'Back    '  <-------------------------
c
c.... if Back, stop the clock
c
        if (ii .eq. icd+2) then
          cpu(icode) = cpu(icode) + secs() - cpu0(icode)
          icode  = icode2
          icode2 = icode3
          return
        endif
c
c.... ------------------->  Individual Processes  <--------------------
c
c.... restart the clock
c
        if (ii .ne. icd+1) then
          icode3 = icode2
          icode2 = icode
          icode  = ii
          cpu0(icode)   = secs()
          nacess(icode) = nacess(icode) + 1
          return
        endif
c
c.... ------------------------>  'End     '  <-------------------------
c
c.... print out execution time statistics
c
        if (ccode(icd) .eq. 'Total   ') return  ! redundant print
c
        call system_clock (icend, icrate, icmax)
        iclock = icend - iclock
        if (iclock .le. 0) iclock = iclock + icmax
        wclock = float(iclock) / float(icrate)
c
        cpu(icd) = secs() - cpu0(icd)
        ccode(icd) = 'Total   '
c
c.... CPU time returned by secs is in 100th of a second
c
        do i = 1, icd
          cpu(i) = cpu(i) / 100. 
        enddo
c
c.... compute the percentage of total time spent
c
        do i = 1, icd
          ratio(i) = 100.d0 * cpu(i) / cpu(icd)
        enddo
c
        if (myrank == master) then
          write (iecho,1000) ititle
          write (iecho,1100)
          write (iecho,1200) (ccode(i), nacess(i), cpu(i), ratio(i),
     &                                                   i=1,icd-4)
          write (iecho,1100)
          write (iecho,1200) (ccode(i), nacess(i), cpu(i), ratio(i),
     &                                                   i=icd-3,icd-1)
          write (iecho,1100)
          write (iecho,1300) (ccode(i),            cpu(i), ratio(i),
     &                                                   i=icd,icd)
          write (iecho,1400)  wclock
        endif
c
c.... print out performance statistics
c
        solvr1 = cpu(3) + cpu(4) + cpu(5)
        solvr2 = cpu(3) + cpu(4) + cpu(5) - cpu(9) - cpu(10)
        gather = cpu(9)
        scattr = cpu(10)
c
        if ((solvr1 .eq. zero) .or. (gather .eq. 0) .or.
     &      (scattr .eq. zero)) return
c
        if (myrank == master) then
          write (iecho,2000)
          write (iecho,2100) flops
          write (iecho,2200) flops/(1.0d6*solvr1)
          write (iecho,2300) flops/(1.0d6*solvr2)
          write (iecho,2400) 8*gbytes
          write (iecho,2500) 8*sbytes
          write (iecho,2600) 8*gbytes/(float(2**20)*gather)
          write (iecho,2700) 8*sbytes/(float(2**20)*scattr)
        endif
c
c.... return
c
        return
c
1000    format(a80,//,
     &  ' E x e c u t i o n   T i m e   S t a t i s t i c s    ',//,
     &  ' name        No. access       CPU-time       %_total')
1100    format(1x,
     &  '---------------------------------------------------')
1200    format(1x,a8,4x,i7,5x,f11.2,8x,f7.2)
1300    format(1x,a8,16x,f11.2,8x,f7.2)
1400    format(/,' Wall-clock time : ',f8.1,' seconds')
2000    format(///,
     &  ' P e r f o r m a n c e   S t a t i s t i c s    ',/)
2100    format(1x,'No. of floating-point operations : ',i20,/)
2200    format(1x,'Flop rate with gather/scatter :   ',f5.1,
     &  ' Mflops/s',/)
2300    format(1x,'Flop rate without gather/scatter :',f5.1,
     &  ' Mflops/s',//)
2400    format(1x,'No. of bytes gathered  : ',i20,/)
2500    format(1x,'No. of bytes scattered : ',i20,/)
2600    format(1x,'Gather transfer rate : ',f6.2, ' Mbytes/s',/)
2700    format(1x,'Scatter transfer rate :',f6.2, ' Mbytes/s')
c
        end
