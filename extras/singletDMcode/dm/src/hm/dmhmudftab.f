***********************************************************************
*** dshmudftab returns the halo velocity distribution
*** (same as dshmudf.f), but reads it from a file.
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
*** Available options
***   type = 1, read a user-supplied file with file name given
***     by udffile in dshmcom.h. If you change the file or
***     for any other reason want to reload it here, you have to
***     set the flag udfload to true, in which case it will
***     be loaded here on next call.
***
*** Input: velocity relative to earth [ km s^-1 ]
*** Output: f(u) / u [ (km/s)^(-2) ]
*** Date: January 30, 2004
***********************************************************************

      real*8 function dshmudftab(u,type)
      implicit none
      include 'dshmcom.h'
      include 'dsidtag.h'
      include 'dssusy.h'

      integer dfmax,dfopt
      parameter (dfmax=1000,  ! max number of lines in files
     &           dfopt=6)     ! number of available options in this file
      real*8 dfu(dfmax,dfopt),dff(dfmax,dfopt)
      integer dfn(dfopt),type
      logical dfload(dfopt)
      data dfload/dfopt*.false./
      common /ntdfvel/dfu,dff,dfn
      save /ntdfvel/,dfload
      
      integer i,fl,l,m
      real*8 u,upl

      character*200 file,scr

c...Read in velocity distribution if first call
      if (type.lt.1.or.type.gt.dfopt) then
        write(*,*) 'ERROR in dshmudftab: invalid type = ',type
        stop
      endif

c...Read in the file if needed
      if (dfload(type).or.(type.eq.1.and.udfload)) then
        dfload(type)=.true.

        if (type.eq.1) then
          file=udffile
          udfload=.false.
        endif

        write(*,*) 'dshmudftab: Opening file ',file
        open(unit=13,file=file,
     &  form='formatted',status='old')
        read(13,'(a)') scr  ! read header lines
        read(13,'(a)') scr  ! read header lines
        dfn(type)=1
 100    read(13,*,end=110) dfu(dfn(type),type),dff(dfn(type),type)
        dfn(type)=dfn(type)+1
        if (dfn(type).gt.dfmax) then
          write(*,*) 'ERROR in dshmudftab:',
     &    '  The data file contains too many values.'
          write(*,*) '  Increase sdmax to include all values.'
          dfn(type)=dfn(type)-1
          goto 110
        endif
        goto 100
 110    continue
        dfn(type)=dfn(type)-1
        close(13)
        write(*,*) '...done, file opened'
      endif

c...Now interpolate in the table
      if (u.le.dfu(1,type).or.u.ge.dfu(dfn(type),type)) then
        dshmudftab=0.0d0
        return
      endif

      call dshunt(dfu(1,type),dfn(type),u,i)
      upl=(u-dfu(i,type))/(dfu(i+1,type)-dfu(i,type))
      dshmudftab=(1.0-upl)*dff(i,type)+upl*dff(i+1,type)

      goto 120

      write(*,*) 
     &  'ERROR in dshmudftab: data file corrupt. Please fix ',file
      dshmudftab=0.0d0
      return

 120  continue

      return
      end
