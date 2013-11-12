***********************************************************************
*** dshmudfearthtab returns the halo velocity distribution
*** (same as dshmudfearth.f), but reads it from a file.
*** The file should have two header lines (with arbitrary content)
*** and then lines with two columns each  with u and u*DF(u).
*** u should be in units of km/s and u*DF(u) (or f(u)/u) in units
*** of (km/s)^(-2).
*** 
*** The file loaded is given by the option type.
*** Some possible types are velocity distributions as obtained from
*** numerical simulations of WIMP propagation in the solar system
*** including solar capture. 
*** 
*** For the simulations made by Johan Lundberg, see astro-ph/0401113,
*** available options are
***   type = 1, reads file <ds-root>/dat/vdfearth-sdbest.dat :
***      best estimate of distribution at Earth from numerical sims
***   type = 2, reads file <ds-root>/dat/vdfearth-sdconserv.dat :
***      conservative estimate, only including free orbits and
***      jupiter-crossing orbits
***   type = 3, reads file <ds-root>/dat/vdfearth-sdultraconserv.dat :
***      ultraconservative estimate, only including free orbits
***   type = 4, reads file <ds-root>/dat/vdfearth-sdgauss.dat :
***      as if Earth was in free space, i.e. gaussian approx.
***   Note: tot.txt is the best estimate of the distribution at Earth
***   and should be used as a default
***
*** There are also other options, like 
***   type = 5, read a user-supplied file with file name given
***     by udfearthfile in dshmcom.h. If you change the file or
***     for any other reason want to reload it here, you have to
***     set the flag udfearthload to true, in which case it will
***     be loaded here on next call.
***
*** Input: velocity relative to earth [ km s^-1 ]
*** Output: f(u) / u [ (km/s)^(-2) ]
*** Date: January 30, 2004
***********************************************************************

      real*8 function dshmuDFearthtab(u,type)
      implicit none
      include 'dshmcom.h'
      include 'dsidtag.h'
      include 'dssusy.h'

      integer sdmax,sdopt
      parameter (sdmax=1000,  ! max number of lines in files
     &           sdopt=6)     ! number of available options in this file
      real*8 sdu(sdmax,sdopt),sdf(sdmax,sdopt)
      integer sdn(sdopt),type
      logical sdload(sdopt)
      data sdload/sdopt*.true./
      common /ntsdvel/sdu,sdf,sdn
      save /ntsdvel/,sdload
      
      integer i,fl,l,m
      real*8 u,upl

      character*200 file,scr

c...Read in velocity distribution if first call
      if (type.lt.1.or.type.gt.sdopt) then
        write(*,*) 'ERROR in dshmudfearthtab: invalid type = ',type
        stop
      endif

c...Read in the file if not already loaded
      if (sdload(type).or.(type.eq.5.and.udfearthload)) then 
        sdload(type)=.false.

        if (type.eq.1) then
          file=dmroot//'dat/'
          call dscharadd(file,'vdfearth-sdbest.dat')
        elseif (type.eq.2) then
          file=dmroot//'dat/'
          call dscharadd(file,'vdfearth-sdconserv.dat')
        elseif (type.eq.3) then
          file=dmroot//'dat/'
          call dscharadd(file,'vdfearth-sdultraconserv.dat')
        elseif (type.eq.4) then
          file=dmroot//'dat/'
          call dscharadd(file,'vdfearth-sdgauss.dat')
        elseif (type.eq.5) then
          file=udfearthfile
          udfearthload=.false.  ! remember to set to .true. if udfearthfile
                                ! is changed an a reload is needed
        endif

        write(*,*) 'dshmuDFearthtab: Opening file ',file
        open(unit=13,file=file,
     &  form='formatted',status='old')
        read(13,'(a)') scr  ! read header lines
        read(13,'(a)') scr  ! read header lines
        sdn(type)=1
 100    read(13,*,end=110) sdu(sdn(type),type),sdf(sdn(type),type)
        sdn(type)=sdn(type)+1
        if (sdn(type).gt.sdmax) then
          write(*,*) 'ERROR in dshmuDFearthtab:',
     &    '  The data file contains too many values.'
          write(*,*) '  Increase sdmax to include all values.'
          sdn(type)=sdn(type)-1
          goto 110
        endif
        goto 100
 110    continue
        sdn(type)=sdn(type)-1
        close(13)
        write(*,*) '...done, file opened'
      endif

c...Now interpolate in the table
      if (u.le.sdu(1,type).or.u.ge.sdu(sdn(type),type)) then
        dshmuDFearthtab=0.0d0
        return
      endif

      call dshunt(sdu(1,type),sdn(type),u,i)
      upl=(u-sdu(i,type))/(sdu(i+1,type)-sdu(i,type))
      dshmuDFearthtab=(1.0-upl)*sdf(i,type)+upl*sdf(i+1,type)

      goto 120

      write(*,*) 
     &  'ERROR in dshmuDFearthtab: data file corrupt. Please fix ',file
      dshmuDFearthtab=0.0d0
      return

 120  continue

      return
      end
