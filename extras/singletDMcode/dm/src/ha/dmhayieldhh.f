c********************************************************************
c hh channel from PPPC for h_mass = 125 GeV
c 
c 24+ May 2013, Christoph Weniger
c
c Returns differential gamma-ray flux [annihilation**-1 GeV**-1]
c Input: x = energy/mass, mx = mass [GeV]
c********************************************************************

      function dmhayieldhh(x,mx)

      include 'dshacom.h'

      character (LEN=100) filein
      logical first
      integer iim, iix, channel
      real*8 mList(62), xList(180)
      real*8 dmhayieldhh,x,mx,iplm,iplx
      data first/.true./
      save first, mList, xList

      channel = 14 ! hh

      ! Load PPPC data first time dmhayieldhh is called
      ! channels (3-14): e mu tau q c b t W Z g gamma h
      ! Tables are in units of dN/d(log10(x))
      if (first) then
        filein = "PPPC/AtProductionNoEW_gammas.dat"
        open(unit=13,file=filein,status='old',form='formatted')
        read(13,'(A)') ! read first comment line
        do i=1,11160
          read(13,2000) (phidiffPPPC(i, l), l=1,14)
        enddo
        close(13)
        first=.false.
        ! 62 mass values (5 GeV - 10 TeV)
        do i=1,62
          mList(i) = phidiffPPPC(1+(i-1)*180, 1)
        enddo
        ! 180 energy bins
        do i=1,180
          xList(i) = phidiffPPPC(i, 2) ! log10(energy/mDM)
        enddo
      endif

      ! find indices
      call dmhaifind(mx, mList, iplm, iim, 1, 61)
      call dmhaifind(log10(x), xList, iplx, iix, 1, 179)

      ! higgs mass is 125 GeV, extrapolate below 130 GeV
      if (mx.ge.125.and.mx.le.130) then
        iim = 18
        iplm = 0d0
      else if (mx.le.125) then
        iim = -5
        iplm = 0d0
      end if

      if (iim.ge.0.and.iix.ge.0) then
        ! interpolate
        dmhayieldhh = 
     &     iplm*iplx*phidiffPPPC(iim*180 + iix + 1, channel)
     &   + (1d0-iplm)*iplx*phidiffPPPC((iim-1)*180 + iix + 1, channel)
     &   + iplm*(1d0-iplx)*phidiffPPPC(iim*180 + (iix-1) + 1, channel)
     &   + (1d0-iplm)*(1d0-iplx)*phidiffPPPC((iim-1)*180 + (iix-1) + 1
     &     , channel)
        ! convert to correct units 
        ! PPPC units are dN/dLog10(x)
        dmhayieldhh = dmhayieldhh/mx/x*0.434294
      else
        dmhayieldhh = 0d0
      end if 

      return

      ! could have been "loadtxt" in python, but this is more fun:
2000  format(F12.0, F12.10, 3(F11.10), 4(F13.10), 2(F12.10), 3(F13.10))
      end
c**********************
