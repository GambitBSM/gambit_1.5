      real*8 function dshrdbardiff(td,solarmod,how)

**********************************************************************
*** function dshrdbardiff calculates the differential flux of
*** antideuterons for the antideuteron kinetic energy per nucleon td 
*** as a result of neutralino annihilation in the halo.
*** compared to dshrdbdiff0, dshrdbardiff uses the rescaled local density
***   input:
***     td - antideuteron kinetic energy per nucleon in gev
***     solarmod - 0 no solar modulation
***                1 solar modulation a la perko
***     how - 1 calculate t_diff only for requested momentum
***           2 tabulate t_diff for first call and use table for
***             subsequent calls
***           3 as 2, but also write the table to disk as pbtd.dat
***             at the first call
***           4 read table from disk on first call, and use that for
***             subsequent calls
*** units: gev^-1 cm^-2 s^-1 sr^-1
*** author: 00-07-19 paolo gondolo
**********************************************************************
      implicit none
      include 'dshmcom.h'
      real*8 td,dshrdbdiff0
      integer solarmod,how
      dshrdbardiff = (rhox/rho0)**2 * dshrdbdiff0(td,solarmod,how)
      return
      end
