*****************************************************************************
*** function dshayields2 calculates the yield above threshold (yieldk=1) or the
*** differential yield (yieldk=2) from a given scalar
*** boson decaying in flight, the energy of the scalar boson should be given
*** in eh.
*** scalars hno = 1-4 are supported (S10, S20, S30 and S+/S-)
*** units: 1.0e-30 m**-2 (annihilation)**-1
*****************************************************************************

      real*8 function dmhayields2(eh,emuth,hno,yieldk,istat)
      implicit none
      include 'dshacom.h'

c------------------------ variables ------------------------------------

      real*8 eh,emuth
      integer hno,istat,yieldk

      real*8 e1,e2,yield

c------------------------ functions ------------------------------------

      real*8 dmhayields3,dmhaemean,dmhayieldfth,dmhayielddec

c-----------------------------------------------------------------------

c...loop through different decay channels
      yield=0.0d0

c---------- neutral scalar bosons ----------
      if (hno.le.3) then

c..."fundamental" channels
        yield=yield+dmhayielddec(eh,hno,emuth,yieldk,istat)

c..."complex" channels
        if (has0br(1,hno).gt.0.0d0) then    ! S10 S10 channel
          e1=dmhaemean(eh,has0m(hno),has0m(1),has0m(1))
          yield=yield+2.0*has0br(1,hno)*dmhayields3(e1,emuth,
     &      1,yieldk,istat)
        endif

        if (has0br(3,hno).gt.0.0d0) then    ! S20 S20 channel
          e1=dmhaemean(eh,has0m(hno),has0m(2),has0m(2))
          yield=yield+2.0*has0br(3,hno)*dmhayields3(e1,emuth,
     &      2,yieldk,istat)
        endif

        if (has0br(4,hno).gt.0.0d0) then    ! S30 S30 channel
          e1=dmhaemean(eh,has0m(hno),has0m(3),has0m(3))
          yield=yield+2.0*has0br(4,hno)*dmhayields3(e1,emuth,
     &      3,yieldk,istat)
        endif

        if (has0br(7,hno).gt.0.0d0) then    ! S+ S- channel
          e1=dmhaemean(eh,has0m(hno),hascm,hascm)
          yield=yield+2.0*
     &      has0br(7,hno)*dmhayields3(e1,emuth,4,yieldk,istat)
        endif

        if (has0br(8,hno).gt.0.0d0) then    ! z0 S10 channel
          e2=dmhaemean(eh,has0m(hno),has0m(1),map(6))
          yield=yield+has0br(8,hno)*
     +    dmhayieldfth(eh,has0m(hno),map(6),
     &       has0m(1),emuth,6,yieldk,istat)
          yield=yield+has0br(8,hno)*dmhayields3(e2,emuth,1,yieldk,istat)
        endif

        if (has0br(9,hno).gt.0.0d0) then    ! z0 S20 channel
          e2=dmhaemean(eh,has0m(hno),has0m(2),map(6))
          yield=yield+has0br(9,hno)*
     +    dmhayieldfth(eh,has0m(hno),map(6),
     &       has0m(2),emuth,6,yieldk,istat)
          yield=yield+has0br(9,hno)*dmhayields3(e2,emuth,2,yieldk,istat)
        endif

        if (has0br(10,hno).gt.0.0d0) then    ! z0 S30 channel
          e2=dmhaemean(eh,has0m(hno),has0m(3),map(6))
          yield=yield+has0br(10,hno)*
     +    dmhayieldfth(eh,has0m(hno),map(6),
     &       has0m(3),emuth,6,yieldk,istat)
          yield=yield
     &       +has0br(10,hno)*dmhayields3(e2,emuth,3,yieldk,istat)
        endif

        if (has0br(11,hno).gt.0.0d0) then    ! w-S+ w+S- channel
          e2=dmhaemean(eh,has0m(hno),hascm,map(5))
          yield=yield+has0br(11,hno)*
     +    dmhayieldfth(eh,has0m(hno),map(5),hascm,emuth,5,yieldk,istat)
          yield=yield
     &      +has0br(11,hno)*dmhayields3(e2,emuth,4,yieldk,istat)
        endif

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

        if (hascbr(13).gt.0.0d0) then ! w+ h1
          yield=yield+hascbr(13)*
     &    dmhayieldfth(eh,hascm,map(5),has0m(1),emuth,
     &    5,yieldk,istat)
          e2=dmhaemean(eh,hascm,has0m(1),map(5))
          yield=yield+hascbr(13)*
     &      dmhayields3(e2,emuth,1,yieldk,istat)
        endif

        if (hascbr(14).gt.0.0d0) then ! w+ h2
          yield=yield+hascbr(14)*
     &    dmhayieldfth(eh,hascm,map(5),has0m(2),emuth,
     &    5,yieldk,istat)
          e2=dmhaemean(eh,hascm,has0m(2),map(5))
          yield=yield+hascbr(14)*
     &      dmhayields3(e2,emuth,2,yieldk,istat)
        endif

        if (hascbr(15).gt.0.0d0) then ! w+ h3
          yield=yield+hascbr(15)*
     &    dmhayieldfth(eh,hascm,map(5),has0m(3),emuth,
     &    5,yieldk,istat)
          e2=dmhaemean(eh,hascm,has0m(3),map(5))
          yield=yield+hascbr(15)*
     &      dmhayields3(e2,emuth,3,yieldk,istat)
        endif

      endif

      dmhayields2 = yield

      end


