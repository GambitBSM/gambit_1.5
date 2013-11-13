*****************************************************************************
***   subroutine dmhainit initializes and loads (from disk) the common
***   block variables needed by the other halo yield routines.  yieldk is
***   the yield type (51,52 or 53 (or 151, 152, 153)) for
***   positron yields, cont. gamma or muon neutrino yields respectively.
***   yieldk is used to check that the provided data file is of the
***   correct type.  if yieldk=51,52 or 53 integrated yields are loaded and
***   if yieldk =151, 152 or 153, differential yields are loaded.
***   author: joakim edsjo
***   edsjo@physto.se date: 96-10-23 (based on dsmuinit.f version
***   3.21) 
***   modified: 98-01-26
*****************************************************************************

      subroutine dmhainit(yieldk)

      implicit none

      include 'dssusy.h'
      include 'dshacom.h'

      data map/1.35d0,5.0d0,175.0d0,1.7841d0,80.25d0,91.2d0,0.10566d0,
     &  0.0d0/
      data mqu,mqd,mqs/0.005d0,0.010d0,0.2d0/
      data lb/10.0d0,10.0d0,176.0d0,10.0d0,80.3d0,91.2d0,10.0d0,
     &  10.0d0/
      data ub/8*5000.0d0/
      data thn/90/
      data mi/10.0d0,25.0d0,50.0d0,80.3d0,91.2d0,100.0d0,150.0d0,
     &  176.0d0,200.0d0,250.0d0,350.0d0,500.0d0,750.0d0,1000.0d0,
     &  1500.0d0,2000.0d0,3000.0d0,5000.0d0/
      data milow/1,1,8,1,4,5,1,1/
      data haexhi,hasmooth,haerr/0,0,0/
      data haftype/'a'/

c...Start-up value
      data dshasetupcalled/.false./



c------------------------ variables ------------------------------------

      integer i,j,k,l,m,yieldk,yieldkfile,fi,fltype,fl
      logical first
      character*200 filein
      character*255 scr
      character*10 scr2
      character*15 filepref
      character*3 filenr
      character*10 filesuf

      data first/.true./
      save first

c...to go from new channel numbers to compressed channels
c...these compressed channel numbers are the indices for the arrays
      data chcomp/ 
     &   11*0,6,5,3*0,7,0,4,0,0,1,0,3,2,8,3*0/



c----------------------------------------- set-up common block variables

      if (first) then

c...masses of annihilation products c, b, t, tau, w and z.
        map(1)=1.35
        map(2)=5.0
        map(3)=175.0
        map(4)=1.7841
        map(5)=80.25
        map(6)=91.2
        map(7)=0.10566
        map(8)=0.0d0

c...lower and upper bounds on neutralino masses for which different
c...annihilation channels can be used
        lb(1)=10.0
        lb(2)=10.0
        lb(3)=176.0
        lb(4)=10.0
        lb(5)=80.3
        lb(6)=91.2
        lb(7)=10.0
        lb(8)=10.0d0
        do i=1,8
          ub(i)=5000.0
        enddo

c...masses for simulation corresponding to mass index i
        mi(1)=10.0
        mi(2)=25.0
        mi(3)=50.0
        mi(4)=80.3
        mi(5)=91.2
        mi(6)=100.0
        mi(7)=150.0
        mi(8)=176.0
        mi(9)=200.0
        mi(10)=250.0
        mi(11)=350.0
        mi(12)=500.0
        mi(13)=750.0
        mi(14)=1000.0
        mi(15)=1500.0
        mi(16)=2000.0
        mi(17)=3000.0
        mi(18)=5000.0

c...lowest mass index for channel j
        milow(1)=1   ! c c-bar
        milow(2)=1   ! b b-bar
        milow(3)=8   ! t t-bar
        milow(4)=1   ! tau+ tau-
        milow(5)=4   ! w+ w-
        milow(6)=5   ! z z
        milow(7)=1   ! mu+ mu-
        milow(8)=1   ! gluons

c...initialize eindex array where the energies for the bins are stored
c...integrated yields (lower end of each bin)
        do i=-1,zn
          zindex(i,1)=dble(i)/dble(zn)
        enddo

c...differential yields (center of each bin)
        do i=-1,zn
          zindex(i,2)=dble(i)/dble(zn)+0.5d0/dble(zn)
          dz(i)=1.0d0/dble(zn)
        enddo

        first=.false.
      endif

c------------------------------------------------------ load yield tables


      call dmhadec(yieldk,fltype,fi)

c...generate file name
      if (haftype.eq.'a') then
        filesuf='.dat'
      else
        filesuf='.bin'
      endif
      if (fltype.eq.1) then   ! integrated
        filepref='simint'
        write(filenr,'(i3)') yieldk
      else                    ! differential
        filepref='simdiff'
        write(filenr,'(i3)') yieldk
      endif
      filein=dmroot//'lib/'//filepref//filenr//filesuf

c...delete spaces in file name
      fl=200
      do l=1,fl
 40     if (filein(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            filein(m:m)=filein(m+1:m+1)
          enddo
          if (fl.eq.l) goto 50
          goto 40
        endif
      enddo
 50   continue 
     
      if (fltype.eq.1) then
c...integrated yields
c        write(*,*) 'enter file name for integrated yield tables:'
c        read(*,'(a)') filein
c        write(*,*) 'is this file ascii (a) or binary (b)?'
c        read(*,*) filetype

c...clear the table
        do j=1,8
          do k=1,18
            do l=0,zn
              phiint(l,k,j,fi)=0.0d0
            enddo
          enddo
        enddo

        if (prtlevel.gt.1) 
     &       write(*,*) 'loading integrated yield tables from file ',
     &       filein
        if (haftype.eq.'a'.or.haftype.eq.'a') then
          open(unit=13,file=filein,status='old',form='formatted')
          read(13,500) scr2,yieldkfile
          read(13,'(a)') scr
        else
          open(unit=13,file=filein,status='old',form='unformatted')
          read(13) yieldkfile
        endif
 500    format(1x,a10,1x,i3)

        if (yieldk.ne.yieldkfile.and.prtlevel.gt.0) then
          write(*,*) 'error in dmhainit: '
          write(*,*) 'the requested yield type is inconsistent with',
     &      ' the provided data file.'
          write(*,*) 'data file: ',filein
          write(*,*) 'requested yieldtype:   ',yieldk
          write(*,*) 'data file''s yieldtype: ',yieldkfile
          write(*,*) 'program stopped.'
          stop
        endif

        yieldtype(1,fi)=yieldkfile
        do j=1,8
          if (prtlevel.gt.1) write(*,*) '      channel number ',j
          do k=1,18
            if (k.ge.milow(j)) then
              if (haftype.eq.'a'.or.haftype.eq.'a') then
                read(13,2000) (phiint(l,k,j,fi),l=0,zn-1)
              else
                read(13) (phiint(l,k,j,fi),l=0,zn-1)
              endif
            endif
          enddo
        enddo
        close(13)

      else   ! differential yields

c...differential yields
c        write(*,*) 'enter file name for differential yield tables:'
c        read(*,'(a)') filein
c        write(*,*) 'is this file ascii (a) or binary (b)?'
c        read(*,*) filetype

c...clear yield table
        do j=1,8
          do k=1,18
            do l=-1,zn
              phidiff(l,k,j,fi)=0.0d0
            enddo
          enddo
        enddo

        if (prtlevel.gt.1) 
     &       write(*,*) 'loading differential yield tables from file ',
     &       filein
        if (haftype.eq.'a'.or.haftype.eq.'a') then
          open(unit=13,file=filein,status='old',form='formatted')
          read(13,500) scr2,yieldkfile
          read(13,'(a)') scr
        else
          open(unit=13,file=filein,status='old',form='unformatted')
          read(13) yieldkfile
        endif

        if (yieldk.ne.yieldkfile.and.prtlevel.gt.0) then
          write(*,*) 'error in dmhainit: '
          write(*,*) 'the requested yield type is inconsistent with',
     &      ' the provided data file.'
          write(*,*) 'data file: ',filein
          write(*,*) 'requested yieldtype:   ',yieldk
          write(*,*) 'data file''s yieldtype: ',yieldkfile
          write(*,*) 'program stopped.'
          stop
        endif

        yieldtype(2,fi)=yieldkfile
        do j=1,8
          if (prtlevel.gt.1) write(*,*) '      channel number ',j
          do k=1,18
            if (k.ge.milow(j)) then
              if (haftype.eq.'a'.or.haftype.eq.'a') then
                read(13,2000) (phidiff(l,k,j,fi),l=0,zn-1)
              else
                read(13) (phidiff(l,k,j,fi),l=0,zn-1)
              endif
              do l=0,zn-1 ! correct units of dyield / dz
                phidiff(l,k,j,fi)=phidiff(l,k,j,fi)/
     &            (dz(l))
              enddo
              phidiff(-1,k,j,fi)=phidiff(0,k,j,fi)
              phidiff(zn,k,j,fi)=phidiff(zn-1,k,j,fi)
            endif
          enddo
        enddo
        close(13)

      endif

      if (prtlevel.gt.1) 
     &     write(*,*) 'loading of halo yield tables finished.'

 2000 format(1000(1x,e12.6))

      return

      end


















