      subroutine dsanset(c)
c...set parameters for annihilation routines
c...  c - character string specifying choice to be made
c...author: joakim edsjo, 2001-09-12

      implicit none
      include 'dsandwcom.h'
      character*(*) c

c...edsjo, gondolo, schelke and ullio, 2001
      if (c.eq.'egsu01'.or.c.eq.'default') then
         slcode=1

c...ellis, falk, olive, srednicki, 1998, 2000, bug-fixed routines
      else if (c.eq.'efos00') then
         slcode=2

c...ellis, falk, olive, srednicki, 1998, 2000, published formulae
      else if (c.eq.'efos00-publ') then
         slcode=3

c...invalid choice
      else
         write (*,*) 'dsanset: unrecognized option ',c
         stop
      endif

      return
      end
