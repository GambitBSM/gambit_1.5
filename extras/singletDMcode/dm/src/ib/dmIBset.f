***********************************************************************
*** Routine dsibset sets default settings for internal bremsstrahlung
*** (IB) calculations.
*** 
***  c - character string specifying choice to be made 
***  (currently for photon and positron yield only)
***  possible values for c:
***  'accurate'             - all channels, high prec. integration
***  'medium'               - all channels, med. prec. integration
***  'fast'                 - most important channels, med. prec. integration
***  'dynamic' or 'default' - decide model by model which channels to include,
***                           med acc of integration
***  'off'                  - no IB contribution
***
*** author: joakim edsjo, 2007-10-19
*** positron switches added: torsten bringmann, 2008-02-29
***********************************************************************

      subroutine dmibset(c)

      implicit none
      include 'dsibcom.h'
      character*(*) c
      integer i

c...all channels, high prec. integration
      if (c.eq.'accurate') then
         ibhow=1  ! fixed set of channels
         IBacc=0.001d0 ! accuracy of numerical integrations
          do i=1,12
            IBflag(i)=1
          enddo 
         IBhow_pos=2 ! calculate positron yield 
                     ! for all channels with IBflag(i)=1
                     ! NB:this can be very time-consuming!


c...all channels, med. prec. integration
       elseif (c.eq.'medium') then
         ibhow=1  ! fixed set of channels
         IBacc=0.01d0 ! accuracy of numerical integrations
          do i=1,12
            IBflag(i)=1
          enddo 
         IBhow_pos=2 ! calculate positron yield
                     ! for all channels with IBflag(i)=1
                     ! NB:this can be very time-consuming!


c...most important channels, med. prec. integration
      elseif (c.eq.'fast') then
         ibhow=1  ! fixed set of channels
         IBacc=0.01d0 ! accuracy of numerical integrations
         do i=1,12
            IBflag(i)=0
         enddo
c...Select which channels to include. Note that these choices can neglect
c...important contributions in special regions, e.g. the stop coannihilation
c...region in mSUGRA
         IBflag(1)=1 ! W+W-
         IBflag(2)=0 ! W+H- and W-H+
         IBflag(3)=0 ! H+H-
         IBflag(4)=1 ! e+ e-
         IBflag(5)=1 ! mu+ mu-
         IBflag(6)=1 ! tau+ tau-
         IBflag(7)=0 ! u u-bar
         IBflag(8)=0 ! d d-bar
         IBflag(9)=0 ! c c-bar
         IBflag(10)=0 ! s s-bar
         IBflag(11)=0 ! t t-bar
         IBflag(12)=0 ! b b-bar
         IBhow_pos=1 ! calculate positron yield only for e+e- gamma


c...decide model by model which channels to include, med acc of integration.
      elseif (c.eq.'dynamic'.or.c.eq.'default') then
         ibhow=2 ! decide channels model by model depending on mass degeneracies
         IBacc=0.01d0 ! accuracy of numerical integrations
c...ibmfr sets up to which mass fraction of the LSP we should include a given
c...final state. If the t-channel exchange sfermion is below this limit
c...that final state is included
c         ibmfr=1.25d0
         ibmfr=2.0d0 ! change to be on the safe side (slower though)
c...calculate positron yield only for e+e- gamma (if IBflag(4)=1)
         IBhow_pos=1


c...Off option
      elseif (c.eq.'off') then
         ibhow=1           
         do i=1,12
            IBflag(i)=0
         enddo
         IBhow_pos=0

      else
         write(*,*) 'ERROR in dsIBset -- unknown option: ',c
         write(*,*) 'Stopping...'
         stop
      endif

      return

      end





        

