*
* Dummy routine to easily integrate bsir364.f with DarkSUSY
* Author: Joakim Edsjo, edsjo@physto.se
* Date: September 13, 2000
*
      SUBROUTINE MTLPRT(NAME,ERC,TEXT)
      CHARACTER*(*) NAME,ERC,TEXT

      write(*,*) 'Error in bsir364: ',name,erc,text

      RETURN

      END
