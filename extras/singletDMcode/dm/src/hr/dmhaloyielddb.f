*****************************************************************************
***   function dshaloyielddb is the version of dshaloyield appropriate
***   for antideuterons
***   yieldk = 159 - antideuteron differential yield
***   author: piero ullio, ullio@sissa.it
***   date: 03-01-14
*****************************************************************************

      real*8 function dshaloyielddb(egev,yieldk,istat)
      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'



c------------------------ functions ------------------------------------

      real*8 dshayield

c------------------------ variables ------------------------------------

      real*8 egev,yield,partial,egevint
      integer ch,istat,yieldk,yieldkint
      real*8 massp,massd,ed,ep,td,tp,nucleon
      parameter (massp=0.93827231d0,massd=1.875612762d0)

c----------------------------------------------- set-up common variables

      if (.not.dshasetupcalled) then
        write(*,*) 'error in dshaloyielddb: dshasetup must be called',
     &    ' before any halo rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

      if(yieldk.ne.159) then
        write(*,*) 'error in dshaloyielddb: dshaloyielddb called',
     &    'with wrong yieldk =',yieldk
        stop
      endif  

c... shift yieldk to the value for pbar
      yieldkint=154
      td=egev
      nucleon=2.d0
      ed=nucleon*egev+massd
      ep=ed/2.d0
      tp=ep-massp
      if(tp.lt.0.d0) then
        write(*,*) 'negative tp in dshaloyielddb'
        write(*,*) 'td, tp = ',td,tp
        tp=0.d0
      endif
      egevint=tp

c...loop through different channels and calculate the yield above threshold
c...for each channel.

c      write(*,*)
c      write(*,*) 'model: ',idtag,'  eth = ',egev
      haistat=0
      yield=0.0d0
      do 100 ch=1,30
        if (habr(ch).gt.0.0d0) then
          partial=dshayield(hamwimp,egevint,ch,yieldkint,istat,mass(kt))
          haistat=or(haistat,istat)
          yield=yield+habr(ch)*0.25d0*partial**2
        endif
  100 continue
      dshaloyielddb=yield
      hristat=or(hristat,haistat)
      istat=haistat

      end



















