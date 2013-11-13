*****************************************************************************
*** function dshayields4 calculates the yield above threshold (yieldk=1) or the
*** differential yield (yieldk=2) from a given scalar
*** boson decaying in flight, the energy of the scalar boson should be given
*** in eh.
*** scalars hno = 1-4 are supported (S10, S20, S30 and S+/S-)
*** units: 1.0e-30 m**-2 (annihilation)**-1
*****************************************************************************

      real*8 function dmhayields4(eh,emuth,hno,yieldk,istat)
      implicit none
      include 'dshacom.h'
      include 'dsidtag.h'

c------------------------ variables ------------------------------------

      real*8 eh,emuth
      integer hno,istat,yieldk
      integer dch

      real*8 yield

c------------------------ functions ------------------------------------

      real*8 dmhayieldfth,dmhayielddec

c-----------------------------------------------------------------------


c...loop through different decay channels
      yield=0.0d0

c---------- neutral scalar bosons ----------
      if (hno.le.3) then

c..."fundamental" channels
         yield=yield+dmhayielddec(eh,hno,emuth,yieldk,istat)

c---------- charged scalar bosons ----------
      else

        if (hascbr(3).gt.0.0d0) then ! u b-bar
          yield=yield+hascbr(3)*
     &    dmhayieldfth(eh,hascm,map(2),mqu,emuth,2,yieldk,istat)
        endif

        if (hascbr(4).gt.0.0d0) then ! c d-bar
          yield=yield+hascbr(4)*
     &    dmhayieldfth(eh,hascm,map(1),mqd,emuth,1,yieldk,istat)
        endif

        if (hascbr(5).gt.0.0d0) then ! c s-bar
          yield=yield+hascbr(5)*
     &    dmhayieldfth(eh,hascm,map(1),mqs,emuth,1,yieldk,istat)
        endif

        if (hascbr(6).gt.0.0d0) then ! c b-bar
          yield=yield+hascbr(6)*
     &    dmhayieldfth(eh,hascm,map(1),map(2),emuth,
     &    1,yieldk,istat)
          yield=yield+hascbr(6)*
     &    dmhayieldfth(eh,hascm,map(2),map(1),emuth,
     &    2,yieldk,istat)
        endif

        if (hascbr(7).gt.0.0d0) then ! t d-bar
          yield=yield+hascbr(7)*
     &    dmhayieldfth(eh,hascm,map(3),mqd,emuth,
     &    3,yieldk,istat)
        endif

        if (hascbr(8).gt.0.0d0) then ! t s-bar
          yield=yield+hascbr(8)*
     &    dmhayieldfth(eh,hascm,map(3),mqs,emuth,
     &    3,yieldk,istat)
        endif

        if (hascbr(9).gt.0.0d0) then ! t b-bar
          yield=yield+hascbr(9)*
     &    dmhayieldfth(eh,hascm,map(3),map(2),emuth,
     &    3,yieldk,istat)
          yield=yield+hascbr(9)*
     &    dmhayieldfth(eh,hascm,map(2),map(3),emuth,
     &    2,yieldk,istat)
        endif

        if (hascbr(12).gt.0.0d0) then ! tau nu_tau
          yield=yield+hascbr(12)*
     &    dmhayieldfth(eh,hascm,map(4),0.0d0,emuth,
     &    4,yieldk,istat)
        endif

      endif

      dmhayields4 = yield

c...check for inconsistent declarations in has0br and hascbr
      if (hno.le.3) then  ! neutral scalars
        do dch=1,11
          if (has0br(dch,hno).gt.0.0d0) then
          write(6,*) 'error in dmhayields4: inconsistent decay widths',
     +      ' declared in has0br'
          write(6,*) 'model: ',idtag
          endif
        enddo
      else                ! charged scalar
        do dch=13,15
          if (hascbr(dch).gt.0.0d0) then
          write(6,*) 'error in dmhayields4: inconsistent decay widths',
     +      ' declared in h0scbr'
          write(6,*) 'model: ',idtag
          endif
        enddo
      endif

      end




