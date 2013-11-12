****************************************************************
***                                                          ***
*** halo local velocity distribution function DF(\vec{v})    ***
*** for the case of an isotropic distribution, i.e. for      ***
***       DF(\vec{v}) = DF(|\vec{v}|) = DF(v)                ***
***                                                          ***
*** dshmDFisotr is normalized such that                      ***
*** int d^3v DF(v) = 4 pi int_0^\infty dv v^2 DF(v) = 1      ***
***                                                          ***
*** Input:  |\vec{v}| = Speed in km/s                        ***
*** Output: DF(|\vec{v}|)  in (km/s)^(-3)                    ***
***                                                          ***
*** Calls other routines depending of choice of velocity     ***
*** distribution function (as set by isodf in the common     ***
*** blocks in dshmcom.h)                                     ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-01-30                                         ***
****************************************************************

      real*8 function dshmDFisotr(v)
      implicit none
      include 'dshmcom.h'
      real*8 v,dshmDFisotrnum
      if (isodf.eq.'num') then
        dshmDFisotr=dshmDFisotrnum(v) 
      else
        write(*,*) 'in function dshmDFisotr unrecognized isodftype flag'
        write(*,*) 'isodf = ',isodf
        write(*,*) 'program stopped'
        stop

      endif
      
      return
      end


