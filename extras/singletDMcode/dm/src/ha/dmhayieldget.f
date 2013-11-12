*****************************************************************************
*** function dmhayieldget gives the information in the differential and
*** integrated arrays phiit and phidiff for given array indices.
*** compared to getting the results directly from the array, this
*** routine performs smoothing if requested.
*** the smoothing is controlled by the parameter hasmooth in the
*** following manner.
*** hasmooth = 0 - no smoothing
***            1 - smoothing of zi-1,zi and zi+1 bins if z>.3
***            2 - smoothing of zi-2,zi-1,zi,zi+1 and zi+2 if z>0.3
*****************************************************************************

      real*8 function dmhayieldget(zi,mxi,ch,fi,ftype,istat)
      implicit none
      include 'dshacom.h'

c------------------------ variables ------------------------------------

      integer zi,mxi,ch,fi,ftype,istat,zsmstart,zimin

c-----------------------------------------------------------------------

      dmhayieldget=0.d0
      zsmstart = int(0.3*zn)
      if (ftype.eq.1) then
        zimin=0
      else
        zimin=-1
      endif

      if (zi.lt.zimin) then
        write(*,*) 'error in dmhayieldget: zi = ',zi
        write(*,*) '  fi = ',fi
        write(*,*) '  ftype = ',ftype
        haerr = 1       
        return  
!        stop
      endif
      if (zi.gt.zn) then
        write(*,*) 'error in dmhayieldget: zi = ',zi
        write(*,*) '  fi = ',fi
        write(*,*) '  ftype = ',ftype
        haerr = 1
        return  
!        stop
      endif


      if (hasmooth.eq.0.or.zi.le.zsmstart) then

        if (ftype.eq.1) then
          dmhayieldget=dble(phiint(zi,mxi,ch,fi))
        else
          dmhayieldget=dble(phidiff(zi,mxi,ch,fi))
        endif

      elseif (hasmooth.eq.1) then

        if (ftype.eq.1) then
          if (zi.ge.1.and.zi.le.(zn-1)) then
            dmhayieldget=0.25d0*dble(phiint(zi-1,mxi,ch,fi))
     &        +0.5d0*dble(phiint(zi,mxi,ch,fi))
     &        +0.25d0*dble(phiint(zi+1,mxi,ch,fi))
          elseif (zi.eq.0) then
            dmhayieldget=0.75d0*dble(phiint(zi,mxi,ch,fi))+
     &        0.25d0*dble(phiint(zi+1,mxi,ch,fi))
          elseif (zi.eq.zn) then
            dmhayieldget=0.75d0*dble(phiint(zi,mxi,ch,fi))+
     &        0.25d0*dble(phiint(zi-1,mxi,ch,fi))
          endif
        else
          if (zi.ge.1.and.zi.le.(zn-1)) then
            dmhayieldget=0.25d0*dble(phidiff(zi-1,mxi,ch,fi))
     &        +0.5d0*dble(phidiff(zi,mxi,ch,fi))
     &        +0.25d0*dble(phidiff(zi+1,mxi,ch,fi))
          elseif (zi.eq.0) then
            dmhayieldget=0.75d0*dble(phidiff(zi,mxi,ch,fi))+
     &        0.25d0*dble(phidiff(zi+1,mxi,ch,fi))
          elseif (zi.eq.zn) then
            dmhayieldget=0.75d0*dble(phidiff(zi,mxi,ch,fi))+
     &        0.25d0*dble(phidiff(zi-1,mxi,ch,fi))
          endif
        endif

      else  ! hasmooth.eq.2

        if (ftype.eq.1) then
          if (zi.le.(zn-2)) then
            dmhayieldget=0.10d0*dble(phiint(zi-2,mxi,ch,fi))
     &        +0.225d0*dble(phiint(zi-1,mxi,ch,fi))
     &        +0.350d0*dble(phiint(zi,mxi,ch,fi))
     &        +0.225d0*dble(phiint(zi+1,mxi,ch,fi))
     &        +0.100d0*dble(phiint(zi+2,mxi,ch,fi))
          elseif (zi.eq.zn-1) then
            dmhayieldget=0.10d0*dble(phiint(zi-2,mxi,ch,fi))
     &        +0.225d0*dble(phiint(zi-1,mxi,ch,fi))
     &        +0.450d0*dble(phiint(zi,mxi,ch,fi))
     &        +0.225d0*dble(phiint(zi+1,mxi,ch,fi))
          else
            dmhayieldget=0.10d0*dble(phiint(zi-2,mxi,ch,fi))
     &        +0.225d0*dble(phiint(zi-1,mxi,ch,fi))
     &        +0.675d0*dble(phiint(zi,mxi,ch,fi))
          endif
        else
          if (zi.le.(zn-2)) then
            dmhayieldget=0.10d0*dble(phidiff(zi-2,mxi,ch,fi))
     &        +0.225d0*dble(phidiff(zi-1,mxi,ch,fi))
     &        +0.350d0*dble(phidiff(zi,mxi,ch,fi))
     &        +0.225d0*dble(phidiff(zi+1,mxi,ch,fi))
     &        +0.100d0*dble(phidiff(zi+2,mxi,ch,fi))
          elseif (zi.eq.zn-1) then
            dmhayieldget=0.10d0*dble(phidiff(zi-2,mxi,ch,fi))
     &        +0.225d0*dble(phidiff(zi-1,mxi,ch,fi))
     &        +0.450d0*dble(phidiff(zi,mxi,ch,fi))
     &        +0.225d0*dble(phidiff(zi+1,mxi,ch,fi))
          else
            dmhayieldget=0.10d0*dble(phidiff(zi-2,mxi,ch,fi))
     &        +0.225d0*dble(phidiff(zi-1,mxi,ch,fi))
     &        +0.675d0*dble(phidiff(zi,mxi,ch,fi))
          endif

        endif

      endif

      end
