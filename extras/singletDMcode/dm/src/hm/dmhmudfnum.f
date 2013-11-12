****************************************************************
***                                                          ***
*** function which gives u*DF(u) where:                      ***
***                                                          ***
***    u is the modulus of \vec{u} = \vec{v} - \vec{v}_{ob}  ***
***    with \vec{v} the 3-d velocity of a WIMP in the        ***
***    galactic frame, and \vec{v}_{ob} the projection on    ***
***    the frame you are considering                         ***
***                                                          ***
***    DF(u) = int dOmega DF(\vec{v}), where                 ***
***    DF(\vec{v}) is the halo local velocity distribution   ***
***    function in the galactic frame                        ***
***                                                          ***
*** on first call the function loads from file a table of    ***
*** values and then interpolates between them.               ***
***                                                          ***
*** the file name is set by the udfnumfile variable          ***
*** it is assumed that this file has no header and u uDF(u)  ***
*** are given with the format 1000 below                     ***
***                                                          ***
*** to reload a (different) file the integer flag uDFnumset  ***
*** into the uDFnumsetcom common block has to be manually    ***
*** reset to 0                                               ***
***                                                          ***
*** u in km s**-1                                            ***
*** u*DF(u) in km**-2 s**2                                   ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-01-30                                         ***
****************************************************************

      real*8 function dshmuDFnum(u)
      implicit none
      include 'dshmcom.h'
      real*8 u,res
      integer ii
      character*200 scr
      real*8 xudfnum(2000),yudfnum(2000),yudfnum2(2000),reudfnum
      common/udfnumsp/xudfnum,yudfnum,yudfnum2,reudfnum
ccc
      if(udfload) then
c  load table for dshmuDFnum from file and set up splines
        open(unit=dfunit,file=udfnumfile,status='old')
        read(dfunit,'(a)') scr  ! read in header line
        read(dfunit,'(a)') scr  ! read in header line
        ii=1
 2      read(dfunit,1000,end=1) xudfnum(ii),yudfnum(ii)
        ii=ii+1
        if(ii.gt.2000) then
          write(*,*) 'exceeded dimension allocated for spline'
          write(*,*) 'variables of udfnumsp common block'
          write(*,*) 'program stopped'
          stop
        endif
        goto 2
 1      close(dfunit)
        reudfnum=dble(ii-1)
        call spline(xudfnum,yudfnum,int(reudfnum),1.d31,1.d31
     &           ,yudfnum2)
        udfload=.false.  ! do not load next time
      endif
      if(u.ge.xudfnum(1).and.u.le.xudfnum(int(reudfnum))) then
        call splint(xudfnum,yudfnum,yudfnum2,int(reudfnum),u,res)
        dshmuDFnum=res
        return
      else
        dshmuDFnum=0.d0
        return
      endif
 1000 format (2(1x,e14.8))
      end


