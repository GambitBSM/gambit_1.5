**********************************************************************
*** function to check if a real*8 number is NaN, returns true if it is
**********************************************************************
      logical function dsisnan(a)
      implicit none
      real*8 a

c...On most machines a comparison with NaN always give false
 
      dsisnan = .not.(a.gt.0d0.or.a.lt.0d0.or.a.eq.0d0)
      if (dsisnan) return

c...On some machine NaN is both >0 and 0 but not <0 at the same time
      dsisnan = (a.gt.0.0d0.and.a.eq.0.0d0)

      end

