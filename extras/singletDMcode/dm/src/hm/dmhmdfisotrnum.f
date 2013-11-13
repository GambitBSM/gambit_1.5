****************************************************************
***                                                          ***
*** halo local velocity distribution function DF(\vec{v})    ***
*** for the case of an isotropic distribution, i.e. for      ***
***       DF(\vec{v}) = DF(|\vec{v}|) = DF(v)                ***
*** as loaded from table in file provided by user            ***
***                                                          ***
*** on first call the function loads from file a table of    ***
*** values and then interpolates between them.               ***
***                                                          ***
*** the file name is set by the isodfnumfile variable        ***
*** it is assumed that this file has no header and v DF(v)   ***
*** are given with the format 1000 below                     ***
***                                                          ***
*** to reload a (different) file the int. flag dfisonumset   ***
*** into the DFisosetcom common block has to be manually     ***
*** reset to 0                                               ***
***                                                          ***
*** v in km s**-1                                            ***
*** DF(v) in km**-3 s**3                                     ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-01-30                                         ***
****************************************************************

      real*8 function dshmDFisotrnum(v)
      implicit none
      include 'dshmcom.h'
      real*8 v,res
      integer ii
      character*200 scr
      real*8 xdfisonum(2000),ydfisonum(2000),ydfisonum2(2000),redfisonum
      common/dfisonumsp/xdfisonum,ydfisonum,ydfisonum2,redfisonum
ccc
      if(isodfload) then
c  load table for dshmuDFnum from file and set up splines
        open(unit=dfunit,file=isodfnumfile,status='old')
        ii=1
        read(dfunit,'(a)') scr   ! read in header line
        read(dfunit,'(a)') scr   ! read in header line
 2      read(dfunit,1000,end=1) xdfisonum(ii),ydfisonum(ii)
        ii=ii+1
        if(ii.gt.2000) then
          write(*,*) 'exceeded dimension allocated for spline'
          write(*,*) 'variables of dfisonumsp common block'
          write(*,*) 'program stopped'
          stop
        endif
        goto 2
 1      close(dfunit)
        redfisonum=dble(ii-1)
        call spline(xdfisonum,ydfisonum,int(redfisonum),1.d31,1.d31
     &           ,ydfisonum2)
        vgalesc=xdfisonum(int(redfisonum))
        isodfload=.false.  ! do not load next time
      endif
      if(v.ge.xdfisonum(1).and.v.le.xdfisonum(int(redfisonum))) then
        call splint(xdfisonum,ydfisonum,ydfisonum2,int(redfisonum)
     &  ,v,res)
        dshmDFisotrnum=res
        return
      else
        dshmDFisotrnum=0.d0
        return
      endif
 1000 format (2(1x,e14.8))
      end


