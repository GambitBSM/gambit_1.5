!Module just wrapping all Jim Cline's routines and expanding a bit on them.
!Pat Scott patscott@physics.mcgill.ca
!April 28 2013

module jimlib

implicit none

integer, parameter, private :: Nhw = 103
double precision :: hm(Nhw),Gh(Nhw),hwspl(Nhw)

contains

!*************************************************************************
      subroutine read_data
      integer i,Nll,Nlln,Np,Nm,Nst
      parameter (Nll=22, Nlln=27,Np=31, Nm=36, Nst=37)
      integer hmi(Nhw)
      real*8 lm(Nll),ll(Nll),llspl(Nll),&
            lmn(Nlln),lln(Nlln),llspln(Nlln),dummy,&
            lmp(Np),lsp(Np),lspspl(Np),lmm(Nm),lsm(Nm),lsmspl(Nm),&
            logmdm(Nst),sig(Nst),sigspl(Nst)
      common/splb/lm,ll,llspl,lmn,lln,llspln
      common/splb2/lmp,lsp,lspspl,lmm,lsm,lsmspl
      common/splsig/logmdm,sig,sigspl

      open(1,file='data/xenon-limit-log.dat',status='old')
      do i=1,Nll
      read(1,*) lm(i), ll(i)
      enddo
      close(1)
      call spline(lm,ll,Nll,1.d31,1.d31,llspl)

      open(1,file='data/new-xenon-p2sig-log.dat',status='old')
      do i=1,Np
      read(1,*) lmp(i), lsp(i)
      enddo
      close(1)
      call spline(lmp,lsp,Np,1.d31,1.d31,lspspl)

      open(1,file='data/new-xenon-m2sig-log.dat',status='old')
      do i=1,Nm
      read(1,*) lmm(i), lsm(i)
      enddo
      close(1)
      call spline(lmm,lsm,Nm,1.d31,1.d31,lsmspl)

!      open(1,file='data/higgs-width2.dat',status='old')
      open(1,file='data/higgs-width3.dat',status='old')
      do i=1,Nhw
!      read(1,*) hm(i), Gh(i)
      read(1,*) hmi(i), dummy,dummy,dummy,dummy,dummy,Gh(i)
      hm(i) = hmi(i)
      enddo
      close(1)
      call spline(hm,Gh,Nhw,1.d31,1.d31,hwspl)

      open(1,file='data/new-xenon-limit-log.dat',status='old')
      do i=1,Nlln
      read(1,*) lmn(i), lln(i)
      enddo
      close(1)
      call spline(lmn,lln,Nlln,1.d31,1.d31,llspln)

      open(1,file='data/steigman.dat',status='old')
      do i=1,Nst
      read(1,*) logmdm(i), sig(i)
      enddo
      close(1)
      call spline(logmdm,sig,Nst,1.d31,1.d31,sigspl)      
      end subroutine

!*************************************************************************

      double precision function higgswid(sqrts)
      use types
      double precision, intent(IN) :: sqrts
      
      if (sqrts .le. 2.d0*mt) then
        call splint(hm,Gh,hwspl,Nhw,sqrts,higgswid)
      else      
        higgswid = GammaGG(sqrts,mW/sqrts) + GammaGG(sqrts,mZ/sqrts)*0.5d0 + Gammaff(sqrts,mt,3,as) 
      endif

      end function

!*************************************************************************

      double precision function GammaGG(m,x)
      use types
      double precision, intent(IN) :: m, x
      GammaGG = m**3/(16.d0*pi*v0*v0) * (1.d0 - 4.d0*x*x + 12*x**4 )*sqrt(1.d0-4.d0*x*x)
      end function GammaGG

!*************************************************************************

      double precision function Gammaff(mh_in,mf,Nc,alphaS)
      use types
      double precision, intent(IN) :: mh_in, mf, alphaS
      integer, intent(IN) :: Nc
      double precision :: vf
      double precision, parameter :: CF = 4.d0/3.d0
      vf = 1.d0-4.d0*mf*mf/(mh_in*mh_in)
      Gammaff = dble(Nc)/(8.d0*pi) * (mf/v0)**2 * mh_in * vf**1.5d0
      Gammaff = Gammaff*(1.d0 + ( 1.5d0 * log(mf*mf/(mh_in*mh_in)) + 2.25d0) * (CF*alphaS/pi) )
      end function Gammaff

!*************************************************************************
      function sigma_SI(l_m,ms)
!      compute spin-independent cross-section in cm^2 
      real*8 yN,l_m,ms,mN,mh,pi,cm2,sigma_SI
      parameter (yN=0.345d0, mN=0.931d0, mh=125.d0, pi=3.1415926536d0, &
        cm2 = 1.d26/0.197**2)
      sigma_SI = (yN*l_m/ms)**2*(mN/mh)**4/(4*pi)
      sigma_SI = sigma_SI/cm2
      end function

!*************************************************************************
      function xenon_limit(l_m,ms)
!      compute the limiting fraction f of the DM component to 
!      saturate XENON100 bound
      logical old
      real*8 yN,l_m,ms,mN,mh,pi,xenon_limit,x,y,logf,cm2,sigma,yp,ym, &
        dy
      parameter (yN=0.345d0, mN=0.931d0, mh=125.d0, pi=3.1415926536d0, &
        cm2 = 1.d26/0.197**2)
      integer i,Nll,Nlln,Np,Nm,XeCL,Xe_astro
      parameter (Nll=22, Nlln=27,Np=31, Nm=36)
      real*8 lm(Nll),ll(Nll),llspl(Nll), &
             lmn(Nlln),lln(Nlln),llspln(Nlln), &
        lmp(Np),lsp(Np),lspspl(Np),lmm(Nm),lsm(Nm),lsmspl(Nm)
      common/splb/lm,ll,llspl,lmn,lln,llspln
      common/splb2/lmp,lsp,lspspl,lmm,lsm,lsmspl
      common/lfb/logf,y
      common/oldb/old
      common/XeCLb/XeCL,Xe_astro


      x = log10(ms)
      if (old) then
         if (x .le. lm(Nll)) then
           call splint(lm,ll,llspl,Nll,x,y)
         else
           y = ll(Nll) + x - lm(Nll)
         endif
      else
         if (x .le. lmn(Nlln)) then
           call splint(lmn,lln,llspln,Nlln,x,y)
         else
           y = lln(Nlln) + x - lmn(Nlln)
         endif
      endif

      if (XeCL.eq.95) then
         call splint(lmp,lsp,lspspl,Np,x,yp)
         call splint(lmm,lsm,lsmspl,Nm,x,ym)
         dy = (yp-ym)/4*1.14 ! increase in log(sigma) from 90 to 95 c.l.
      else
         dy = 0.d0
      endif      

      if (Xe_astro.eq.1) dy = dy -0.29128 + 0.4557*x - 0.081349*x**2
      if (Xe_astro.eq.-1) dy = dy -1.3329 + 1.4365*x - 0.59337*x**2  &
        + 0.080595*x**3

      sigma = (yN*l_m/ms)**2*(mN/mh)**4/(4*pi)
      logf = y + dy - log10(sigma/cm2)
      xenon_limit = 10.d0**logf
      end function

!*************************************************************************
      function lm_max(ms,BFinvis)
        use types
        double precision lm_max
        double precision, intent(IN) :: ms, BFinvis
        lm_max = BFinvis / (1.d0 - BFinvis) * Gh0 * pi * 2.d0
        lm_max = lm_max / sqrt(mh*mh - 4.d0*ms*ms)
        lm_max = 4.d0 * mh/v0 * sqrt(lm_max)

      end function

!*************************************************************************
      function invisible_higgs_width(ms,lm)
        use types
        double precision invisible_higgs_width
        double precision, intent(IN) :: ms, lm
        invisible_higgs_width = lm*lm*v0*v0*sqrt(1.d0-4.d0*ms*ms/(mh*mh))/(32.d0*pi*mh)
        
        end function

!*************************************************************************
      function fractions(lm,ms)
      use types
      Type(BFset) :: fractions
      double precision, intent(IN) :: lm,ms
      double precision :: sqrts,Ghs,Ghmh,denom,sv_common,rh,rw,rz
      double precision :: sv,oldsv,newsv,sv_h,sv_w,sv_z,sv_u,sv_d,sv_s,sv_c
      double precision :: sv_b,sv_t,sv_e,sv_mu,sv_tau,sv_gamgam
      double precision :: total_calculated,sv_factorised
      
        call clear(fractions)

        rh = (mh/ms)**2               
        sqrts = 2*ms
        Ghmh = higgswid(mh)
        if (ms .le. mh*0.5) Ghmh = Ghmh + invisible_higgs_width(ms,lm)
               
        denom = (4-rh)**2 + rh*(Ghmh/ms)**2
        sv_common = (lm/ms)**2 / (8.d0*pi) / denom

        if (ms.gt.mh) then
          sv_h = sv_common * 2.d0 * ( rh*(lm/lh)*(1.d0-0.25d0*rh)/(rh-2.d0) + 1.d0 + rh*0.5d0 )**2 * sqrt(1.d0-rh)
        else
          sv_h = 0.d0
        endif
        if (ms.gt.mw) then
          rw = (mw/ms)**2
          sv_w = sv_common * rw*rw*(2.d0+(1.d0-2.d0/rw)**2)*sqrt(1.d0-rw)
        else
          sv_w = 0.d0
        endif
        if (ms.gt.mz) then 
          rz = (mz/ms)**2
          sv_z = sv_common * 0.5d0*rz*rz*(2.d0+(1.d0-2.d0/rz)**2)*sqrt(1.d0-rz)
        else 
          sv_z = 0.d0
        endif
        sv_u = sv_common*fermion_extra(mu,ms) !* ( 1.d0 + (1.5d0*log((mu/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )   !1-loop QCD correction; makes virtually no difference
        sv_d = sv_common*fermion_extra(md,ms) !* ( 1.d0 + (1.5d0*log((md/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
        sv_s = sv_common*fermion_extra(mst,ms)!* ( 1.d0 + (1.5d0*log((mst/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
        sv_c = sv_common*fermion_extra(mc,ms) !* ( 1.d0 + (1.5d0*log((mc/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
        sv_b = sv_common*fermion_extra(mb,ms) !* ( 1.d0 + (1.5d0*log((mb/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
        sv_t = sv_common*fermion_extra(mt,ms) !* ( 1.d0 + (1.5d0*log((mt/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
        sv_e = sv_common*fermion_extra(me,ms)/3.d0
        sv_mu = sv_common*fermion_extra(mmu,ms)/3.d0
        sv_tau = sv_common*fermion_extra(mtau,ms)/3.d0
        sv_gamgam = 0.d0

        total_calculated = sv_h + sv_w + sv_z + sv_u + sv_d + sv_s + sv_c + sv_b + sv_t + sv_e + sv_mu + sv_tau + sv_gamgam
        sv = total_calculated

        if (ms .lt. factorisation_maxmass .and. ms .gt. 0.5d0*mh) then
          Ghs = higgswid(sqrts)
          sv_factorised = 8.d0*(v0*lm/ms)**2*Ghs/sqrts**3.d0/denom + sv_h
          if (sv_factorised .gt. total_calculated) sv = sv_factorised
        endif

        fractions%h = sv_h/sv
        fractions%w = sv_w/sv
        fractions%z = sv_z/sv
        fractions%u = sv_u/sv
        fractions%d = sv_d/sv
        fractions%s = sv_s/sv
        fractions%c = sv_c/sv
        fractions%b = sv_b/sv
        fractions%t = sv_t/sv
        fractions%e = sv_e/sv
        fractions%mu = sv_mu/sv
        fractions%tau = sv_tau/sv
        fractions%gamgam = sv_gamgam/sv
        if (ms .le. 1.03d0*mz .and. ms .gt. 0.5d0*mh) then !Assume that the remainder goes into 3/4-body final states associated with the next threshold
          if (ms .le. 1.03d0*mw) then
            fractions%Wto4 = (sv - total_calculated)/sv
            fractions%Zto4 = 0.d0
          else
            fractions%Zto4 = (sv - total_calculated)/sv
            fractions%Wto4 = 0.d0
          endif
        else
          fractions%Zto4 = 0.d0
          fractions%Wto4 = 0.d0
        endif

        end function

!*************************************************************************
        function fermion_extra(mx,ms)
        double precision, intent(IN) :: mx,ms
        double precision :: fermion_extra, r

        if (ms .gt. mx) then
          r = (mx/ms)**2
          fermion_extra = 6.d0*r*(1.d0-r)**1.5d0
        else 
          fermion_extra = 0.d0
        endif

        end function

!*************************************************************************
      function sigmav(lm,ms,y)
      use types
      real*8 lm,ms,sv,sigmav,rw,rh,rz,Ghs,sqrts
      double precision, optional :: y
      double precision :: sv_u,sv_d,sv_s,sv_c,y_loc,denom,sv_common
      double precision :: sv_b,sv_t,sv_e,sv_mu,sv_tau,sv_gamgam,Ghmh

        if (present(y)) then
          y_loc = y
        else
          y_loc = 1.d0
        endif          

        rh = (mh/ms)**2
        sqrts = 2.d0*ms        
        Ghmh = higgswid(mh)
        if (ms .le. mh*0.5) Ghmh = Ghmh + invisible_higgs_width(ms,lm)
        
        if (ms .gt. mh .or. ms .ge. factorisation_maxmass) then
          denom = (4-rh)**2 + rh*(Ghmh/ms)**2
          sv_common = (lm/ms)**2 / (8.d0*pi) / denom
        endif


        if (ms .lt. factorisation_maxmass) then

          Ghs = higgswid(sqrts*sqrt(y_loc))
          sv = 8*(v0*lm/ms)**2*Ghs/sqrts**3.d0/sqrt(y_loc)
          sv = sv / ((4*y_loc-rh)**2 + rh*(Ghmh/ms)**2)
          
        else

          sv_u = sv_common*fermion_extra(mu,ms) !* ( 1.d0 + (1.5d0*log((mu/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )  !1-loop QCD correction; makes virtually no difference
          sv_d = sv_common*fermion_extra(md,ms) !* ( 1.d0 + (1.5d0*log((md/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
          sv_s = sv_common*fermion_extra(mst,ms)!* ( 1.d0 + (1.5d0*log((mst/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
          sv_c = sv_common*fermion_extra(mc,ms) !* ( 1.d0 + (1.5d0*log((mc/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
          sv_b = sv_common*fermion_extra(mb,ms) !* ( 1.d0 + (1.5d0*log((mb/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
          sv_t = sv_common*fermion_extra(mt,ms) !* ( 1.d0 + (1.5d0*log((mt/sqrts)**2) + 2.25d0)*4.d0*as/(3.d0*pi) )
          sv_e = sv_common*fermion_extra(me,ms)/3.d0
          sv_mu = sv_common*fermion_extra(mmu,ms)/3.d0
          sv_tau = sv_common*fermion_extra(mtau,ms)/3.d0
          sv_gamgam = 0.d0
          sv = sv_u + sv_d + sv_s + sv_c + sv_b + sv_t + sv_e + sv_mu + sv_tau + sv_gamgam

          if (ms.gt.mw) then
            rw = (mw/ms)**2
            sv = sv + sv_common * rw*rw*(2.d0+(1.d0-2.d0/rw)**2)*sqrt(1.d0-rw)
          endif
          if (ms.gt.mz) then 
            rz = (mz/ms)**2
            sv = sv + sv_common * 0.5d0*rz*rz*(2.d0+(1.d0-2.d0/rz)**2)*sqrt(1.d0-rz)
          endif

        endif

        if (ms.gt.mh) sv = sv + sv_common * 2.d0 * ( rh*(lm/lh)*(1.d0-0.25d0*rh)/(rh-2.d0) + &
         1.d0 + rh*0.5d0 )**2 * sqrt(1.d0-rh)

        sigmav = sv
        
      end function

!************************************************************************
      function relic_f_steigman(lm,ms)
      use types
      integer Nst
      parameter (Nst=37)
      real*8 lm,ms,sv,sv0_std,cm2,relic_f_steigman,logmdm(Nst),sig(Nst),sigspl(Nst),sv1,sv0
      parameter (cm2 = 1.d26/0.197**2, sv0_std=1.d-36*cm2)
      character (LEN=200) :: msstring, lmstring, command
      common/splsig/logmdm,sig,sigspl

        !if (ms .lt. mh*0.5d0) then
        
        !  write(msstring, '(e16.6)') ms
        !  write(lmstring, '(e16.6)') lm
        !  command = 'octave -qf kimmolib.m '//trim(msstring)//' '//trim(lmstring)//' > '//trim(tempfile)
        !  call system(command)
        !  open(1,file=tempfile,status='old')
        !  read(1,*) relic_f_steigman
        !  close(1)

        !else

          sv = sigmav(lm,ms)
          call splint(logmdm,sig,sigspl,Nst,log10(ms),sv1) 
          sv0 = sv1/3*sv0_std ! exact relic density cross section
          relic_f_steigman = sv0/sv

        !endif

      end function

!***********************************************************************
      function svrelo(lm,ms)
      integer Nst
      parameter (Nst=37)
      real*8 lm,mb,ms,mh,Gh0,sv0_std,cm2,pi,lh,mZ,mW,&
            rw,rh,rz,Ghs,sqrts,v0,svrelo,&
            logmdm(Nst),sig(Nst),sigspl(Nst),sv1,sv0,lm1,ms1,sv,&
            dv,rb,Ghtot,v1,v2,v_higgs,svrel1,svrel2,norm,beta_m,&
            norm1,norm2,vm
      parameter (mb=4.18d0, Gh0=4.07d-3, cm2 = 1.d26/0.197**2,lh=0.129d0,&
            pi=3.1415926536d0, mh=125.d0, sv0_std=1.d-36*cm2, &
            mW=80.385d0, mZ=91.188d0, v_higgs=246.22d0)
      common/splsig/logmdm,sig,sigspl
      common/lmms/lm1,ms1,rh,rb,Ghtot
      common/bmb/beta_m
      common/normb/norm

!      write(6,*) 'beta_m = ',beta_m
!      beta_m = 19.3d0 ! 32.9

      vm = sqrt(2.)/beta_m*sqrt(1+sqrt(1+beta_m**2))
      call qromb(norm_int,0.d0,vm,norm1)
      call qromo(norm_int,vm,1.d10,norm2,midinf)      
      norm = norm1 + norm2
      write(6,*) 'beta_m, norm = ',beta_m, norm

      lm1 = lm
      ms1 = ms
      rb = (mb/ms)**2
      rh = (mh/ms)**2
      Ghtot = Gh0
      if (ms.lt.mh/2) Ghtot = Ghtot &
            + (lm*v_higgs)**2/mh/(32*pi)*sqrt(1-4/rh)

      rh = (mh/ms)**2
      if (rh.gt.4) then
         v0 = sqrt(rh/4-1)
         dv = v0*(sqrt(1 +sqrt(rh)*Ghtot/ms/4/v0*2)-1)
         v1 = max(v0-500*dv,0.d0)
         v2 = v0+500*dv
!         call qromo(sv_int,v1,v0,svrel1,midinf)
!         call qromo(sv_int,v0,v2,svrel2,midinf)
         call qromb(sv_int,0.d0,v0,svrel1)
!         write(6,*) 'svrel1 = ', svrel1
!         call test_int(v1,v0)
!         stop

         call qromo(sv_int,v0,1.d10,svrel2,midinf)
!         write(6,*) 'svrel2 = ', svrel2
!         call test_int(v0,v2)
!         stop

      else
         v1 = 0.d0
         v2 = vm
         call qromo(sv_int,v2,1.d10,svrel1,midinf)
         call qromb(sv_int,v1,v2,svrel2)
!         write(6,*) 'svrel = ', svrel
!         call test_int(v1,v2)
!         stop

      endif      

         svrelo = svrel1 + svrel2      


      end function

!***********************************************************************
      function sv_int(v)
!      sv = sigma *vrel as a function of velocity of DM particles
!      Ghtot must be computed elsewhere before this is called
      integer Nst
      parameter (Nst=37)
      real*8 lm,mb,ms,mh,Gh0,sv,sv0_std,cm2,pi,lh,mZ,mW,&
            rw,rh,rz,Ghs,sqrts,v0,rb,&
            logmdm(Nst),sig(Nst),sigspl(Nst),sv1,sv0,v,Ghtot,gamma,&
            sv_int,beta_m,norm/0.01535/
      parameter (mb=4.18d0, Gh0=4.07d-3, cm2 = 1.d26/0.197**2,lh=0.129d0,&
            pi=3.1415926536d0, mh=125.d0, sv0_std=1.d-36*cm2, &
            mW=80.385d0, mZ=91.188d0,v0=246.22d0 )
      common/splsig/logmdm,sig,sigspl
      common/lmms/lm,ms,rh,rb,Ghtot
      common/bmb/beta_m
      common/normb/norm
      common/svb/sv

      gamma = 1+v**2

!      sv = lm**2/(8*pi*ms**2*gamma*((4*gamma-rh)**2 + rh*Ghtot/ms**2))
!     &      * 6*rb*(1-rb)**1.5d0

      sqrts = 2*ms
      Ghs = higgswid(sqrts*sqrt(gamma))
      sv = 8*(v0*lm/ms)**2*Ghs/sqrts**3.d0/gamma/&
        ((4*gamma-rh)**2 + rh*((Ghs+Ghtot-Gh0)/ms)**2)


      sv_int = v**2*sv*exp(-beta_m*(sqrt(gamma)-1))/norm

      end function

!***********************************************************************
      function relic_f_steigman2(lm,ms)
      integer Nst
      parameter (Nst=37)
      real*8 lm,mb,ms,mh,Gh0,sv0_std,cm2,pi,relic_f_steigman2,lh,mZ,mW,&
            rw,rh,rz,Ghs,sqrts,v0,svrel,&
            logmdm(Nst),sig(Nst),sigspl(Nst),sv1,sv0,lm1,ms1,sv,&
            dv,rb,Ghtot,v1,v2,v_higgs,svrel1,svrel2,norm,beta_m

      parameter (mb=4.18d0, Gh0=4.07d-3, cm2 = 1.d26/0.197**2,lh=0.129d0,&
            pi=3.1415926536d0, mh=125.d0, sv0_std=1.d-36*cm2, &
            mW=80.385d0, mZ=91.188d0, v_higgs=246.22d0)
      common/splsig/logmdm,sig,sigspl
      common/lmms/lm1,ms1,rh,rb,Ghtot
      common/svrb/svrel
      common/bmb/beta_m
      common/normb/norm

!      write(6,*) 'beta_m = ',beta_m
!      beta_m = 19.3d0 ! 32.9

      call qromb(norm_int,0.d0,20.d0/beta_m,norm)

      lm1 = lm
      ms1 = ms
      rb = (mb/ms)**2
      rh = (mh/ms)**2
      Ghtot = Gh0
      if (ms.lt.mh/2) Ghtot = Ghtot &
            + (lm*v_higgs)**2/mh/(32*pi)*sqrt(1-4/rh)

      rh = (mh/ms)**2
      if (rh.gt.4) then
         v0 = sqrt(rh/4-1)
         dv = v0*(sqrt(1 +sqrt(rh)*Ghtot/ms/4/v0*2)-1)
         v1 = max(v0-500*dv,0.d0)
         v2 = v0+500*dv
!         call qromo(sv_int,v1,v0,svrel1,midinf)
!         call qromo(sv_int,v0,v2,svrel2,midinf)
         call qromb(sv_int,v1,v0,svrel1)
!         write(6,*) 'svrel1 = ', svrel1
!         call test_int(v1,v0)
!         stop

         call qromb(sv_int,v0,v2,svrel2)
!         write(6,*) 'svrel2 = ', svrel2
!         call test_int(v0,v2)
!         stop


         svrel = svrel1 + svrel2      
      else
         v1 = 0.d0
         dv = 0.5*sqrt(sqrt(2*(4-rh)**2 + rh*(Ghtot/ms)**2)&
            -(4-rh))
         v2 = 10*dv
!         call qromo(sv_int,v1,v2,svrel,midinf)
         call qromb(sv_int,v1,v2,svrel)
!         write(6,*) 'svrel = ', svrel
!         call test_int(v1,v2)
!         stop

      endif      

      call splint(logmdm,sig,sigspl,Nst,log10(ms),sv1) 
      
      sv0 = sv1/3*sv0_std ! exact relic density cross section
      relic_f_steigman2 = sv0/svrel
      end function

!*************************************************************************
      function relic_f_old(lm,ms)
      real*8 lm,mb,ms,mh,Gh,sv,sv0,cm2,pi,relic_f_old,lh,mZ,mW, &
        rw,rh,rz
      parameter (mb=4.18d0, Gh=4.07d-3, cm2 = 1.d26/0.197**2,lh=0.129d0, &
        pi=3.1415926536d0, mh=125.d0, sv0=1.d-36*cm2,  &
        mW=80.385d0, mZ=91.188d0 )

      sv = 6*(lm*mb)**2/((4*ms**2-mh**2)**2 + (Gh*mh)**2)/(8*pi)
      rh = (mh/ms)**2
      if (ms.gt.mW) then
         rw = (mW/mS)**2
         sv = sv + 4*(lm*rw)**2/(4-rh)**2*(2+(1-2/rw)**2) * &
         sqrt(1-rw)/(32*pi*ms**2)
      endif

      if (ms.gt.mZ) then
         rz = (mZ/mS)**2
         sv = sv + 4*(lm*rz)**2/(4-rh)**2*(2+(1-2/rz)**2) * &
         sqrt(1-rz)/(64*pi*ms**2)
      endif
        
      if (ms.gt.mh) sv = sv + (lm/ms)**2/(64*pi)*( rh/(rh-2)*(lm/lh) + &
        (1+rh/2)/(1-rh/4) )**2*sqrt(1-rh)

      relic_f_old = sv0/sv
      end function


!*************************************************************************

      double precision function relic_f(lm,ms)
      double precision lm, ms, ms_shared, xf0, Yinf
      integer bessel_switch
      common/bsb/bessel_switch
      common/msb/ms_shared
      bessel_switch = 1
      ms_shared = ms
       
      xf0 = 20.d0
      call get_xf(lm,ms,xf0)
      call spline_svrel(lm,ms,xf0)
      call igr_boltz(xf0-20.d0,1.d6,Yinf)
      relic_f = 2.47d9*ms*Yinf

      end function


!*************************************************************************
      subroutine igr_boltz(x1,x2,Yinf)
      logical fail
      integer nvar/1/,nok,nbad,kmax,kount,NMAX,KMAXX,i
      real*8 x1,x2,Yinf,ystart(1),eps/1.d-6/,h1,hmin,dxsav,x12
      PARAMETER (NMAX=50,KMAXX=2000)
      REAL*8 xp(KMAXX), yp(NMAX,KMAXX)
      COMMON /path/ kmax,kount,dxsav,xp,yp


      x12 = x1


1     dxsav = 0.1d0
      kmax = KMAXX
      hmin = 0.d0
      h1 = 0.1d0
      ystart(1) = Yeq(x12)
!      write(6,*)  'x1,x2 = ',x12,x2

      call odeint(ystart,nvar,x12,x2,eps,h1,hmin,nok,nbad,boltz_deriv,rkqs,fail)

!      write(6,*) 'kount = ',kount
      
      if (kount.lt.15) then
         x12 = x12 + 1
         goto 1
      endif

      Yinf = yp(1,kount)

      return

      open(30,file='Y.dat',status='unknown')
      do i=1,kount
      write(30,*) log10(xp(i)),log10(yp(1,i))
      enddo
      close(30)

      write(6,*) 'Yinf = ',Yinf

      end subroutine

!*************************************************************************
      subroutine boltz_deriv(x,y,dydx)
      real*8 x(1), y(1), dydx(1), sigv, ms, mpl, gss, pi, c, msb
      parameter(mpl=1.22d19, gss=86.25, pi=3.1415926536d0, c = 2*pi**2/1.66/45*sqrt(gss)*mpl)
      common/msb/ms
      
      dydx(1) = -ms*c*svr_spl(x(1))*(y(1)**2 - Yeq(x(1))**2)/x(1)**2

!      write(6,*) dydx(1)

      end subroutine

!*************************************************************************
      function Yeq(x)
      real*8 Yeq, x, gss/86.25/,c/0.145/
      Yeq = c/gss*x**1.5d0*exp(-x)
      end function

!*************************************************************************
      subroutine get_xf(lm,ms,xf)
      integer i,Ni/1/,g/1/
      real*8 lm,ms,bm,bm0/17./,bm1/23./,dbm,svrel,fr,norm
      real*8 s0/0.267/,mpl/1.22d19/,n/1.6/,gs/86.25/,q,xf,lbm(0:1),lsv(0:1),ls0,xfold,pi/3.1415926536d0/
      common/bmb/bm
      common/svrb/svrel
      common/normb/norm

      xfold = -100.d0


      do while (abs(xf-xfold).gt.0.1)
      
         bm = xf
         xfold = xf
         fr = relic_f_steigman2(lm,ms)
         xf = log((0.6*svrel*ms*mpl*sqrt(bm))/(2*pi)**1.5d0*sqrt(gs))
!         write(6,*) 'xf = ',xf      
      
      enddo

      end subroutine

!*************************************************************************
      function svr_spl(xf)
      integer Ns
      parameter (Ns=100)
      real*8 xf, svr_spl,bm(Ns),svr(Ns),svrs(Ns)
      common/svrspl/bm,svr,svrs

      if (xf.lt.bm(1)) then
         svr_spl = svr(1)
      else if (xf.gt.bm(Ns)) then
         svr_spl = svr(Ns)
      else
         call splint(bm,svr,svrs,Ns,xf,svr_spl)
      endif
      end function

!*************************************************************************
      subroutine spline_svrel(lm,ms,xf0)
      integer i,Ns
      parameter (Ns=100)
      real*8  bm(Ns),svr(Ns),svrs(Ns),xf0,xf1/3.d0/,xf2,xf,dxf,&
       lm, ms,svmin, dummy,logxf1,logxf2,dlogxf,logxf
      common/bmb/xf
      common/svrspl/bm,svr,svrs
      common/svb/svmin

      logxf1 = log10(xf1)
      xf2 = xf0+10000.d0
      logxf2 = log10(xf2)
      dlogxf = (logxf2-logxf1)/(Ns-1)

!      dummy = sv_int(0.d0)

      
      do i=1,Ns
        logxf = logxf1 + (i-1)*dlogxf 
        xf = 10.d0**logxf
        bm(i) = xf
!        svr(i) = max(svrelo(lm,ms),svmin)
        svr(i) = svrel_new(lm,ms)
      enddo



!      open(30,file='spltest.dat',status='unknown')
!      do i=1,Ns
!        logxf = logxf1 + (i-1)*dlogxf 
!        write(30,*) (logxf), (svr(i))
!      enddo
!      close(30)
!      stop


      call spline(bm,svr,Ns,1.d31,1.d31,svrs)
      return

      open(30,file='spltest.dat',status='unknown')
      dxf = (xf2-xf1)/9999
      do i=1,10000
        xf = xf1 + (i-1)*dxf 
        write(30,*) log10(xf), log10(svr_spl(xf))
      enddo
      close(30)
      stop

      end subroutine

!*************************************************************************
      function Gamma_inv(lm,ms)
      real*8 Gamma_inv,lm,ms,mh/125./,v0/246./,pi/3.1415926536/,r

      r = 4*(ms/mh)**2
      if (r.lt.1) then
         Gamma_inv = 1e3* lm**2*v0**2/(32*pi*mh)*sqrt(1-r)
      else
         Gamma_inv = 0.d0
      endif
      end function

!***********************************************************************
      function norm_int(v)
      real*8 norm_int,v,beta_m
      common/bmb/beta_m

      norm_int = v**2*exp(-beta_m*(sqrt(1+v**2)-1))

      end function

!***********************************************************************
      function sv_int_new_old(y)
!      sv = sigma *vrel as a function of velocity of DM particles
!      Ghtot must be computed elsewhere before this is called
      integer Nst
      parameter (Nst=37)
      real*8 lm,mb,ms,mh,Gh0,sv,sv0_std,cm2,pi,lh,mZ,mW,&
            rw,rh,rz,Ghs,sqrts,v0,rb,&
            logmdm(Nst),sig(Nst),sigspl(Nst),sv1,sv0,v,Ghtot,gamma,&
            sv_int_new_old,beta_m,norm/0.01535/,y,sv2,num
      parameter (mb=4.18d0, Gh0=4.07d-3, cm2 = 1.d26/0.197**2,lh=0.129d0,&
            pi=3.1415926536d0, mh=125.d0, sv0_std=1.d-36*cm2, &
            mW=80.385d0, mZ=91.188d0,v0=246.22d0 )
      common/splsig/logmdm,sig,sigspl
      common/lmms/lm,ms,rh,rb,Ghtot
      common/bmb/beta_m
      common/normb/norm
      common/svb/sv
      common/svb2/sv2,num
      integer bs
      common/bsb/bs

      sqrts = 2*ms
      Ghs = higgswid(sqrts*sqrt(y))
      
      sv = 8*(v0*lm/ms)**2*Ghs/sqrts**3.d0/sqrt(y)/&
        ((4*y-rh)**2 + rh*((Ghs+Ghtot-Gh0)/ms)**2)
      rh = (mh/mS)**2
      if (rh.lt.1) sv = sv + (lm/ms)**2/(4.d0*pi)*sqrt(1.d0-rh)*&
       (lm/lh*rh*(1.d0-rh/4.d0)/(rh-2.d0)+ 1+rh/2.d0)**2 /&
       ((4.d0-rh)**2 + rh*((Ghs+Ghtot-Gh0)/ms)**2) 
!      sv2 = 1./ ((4*y-rh)**2 + rh*((Ghs+Ghtot-Gh0)/ms)**2)
!      num = Ghs/y

      sv_int_new_old = sv*sqrt(y*(y-1))*bessk1(2*beta_m*sqrt(y))/norm&
            * sqrt(y)*exp(bs*2*beta_m*(1-sqrt(y)))

      end function

!***********************************************************************
      function sv_int_new(y)
!      sv = sigma *vrel as a function of velocity of DM particles
      integer Nst
      parameter (Nst=37)
      real*8 lm,mb,ms,mh,Gh0,sv,sv0_std,cm2,pi,lh,mZ,mW,&
            rw,rh,rz,Ghs,sqrts,v0,rb,&
            logmdm(Nst),sig(Nst),sigspl(Nst),sv1,sv0,v,Ghtot,gamma,&
            sv_int_new,beta_m,norm/0.01535/,y,sv2,num
      parameter (mb=4.18d0, Gh0=4.07d-3, cm2 = 1.d26/0.197**2,lh=0.129d0,&
            pi=3.1415926536d0, mh=125.d0, sv0_std=1.d-36*cm2, &
            mW=80.385d0, mZ=91.188d0,v0=246.22d0 )
      common/splsig/logmdm,sig,sigspl
      common/lmms/lm,ms,rh,rb,Ghtot
      common/bmb/beta_m
      common/normb/norm
      common/svb/sv
      common/svb2/sv2,num
      integer bs
      common/bsb/bs

      sv = sigmav(lm,ms,y=y)
      sv_int_new = sv*sqrt(y*(y-1))*bessk1(2*beta_m*sqrt(y))/norm&
            * sqrt(y)*exp(bs*2*beta_m*(1-sqrt(y)))

      end function

!***********************************************************************
      function svrel_new(lm,ms)
      integer Nst
      parameter (Nst=37)
      real*8 lm,mb,ms,mh,Gh0,sv0_std,cm2,pi,lh,mZ,mW,&
            rw,rh,rz,Ghs,sqrts,v0,svrel_new,&
            logmdm(Nst),sig(Nst),sigspl(Nst),sv1,sv0,lm1,ms1,sv,&
            dv,rb,Ghtot,v1,v2,v_higgs,svrel1,svrel2,norm,beta_m,&
            y0,vm
      parameter (mb=4.18d0, Gh0=4.07d-3, cm2 = 1.d26/0.197**2,lh=0.129d0,&
            pi=3.1415926536d0, mh=125.d0, sv0_std=1.d-36*cm2, &
            mW=80.385d0, mZ=91.188d0, v_higgs=246.22d0)
      common/splsig/logmdm,sig,sigspl
      common/lmms/lm1,ms1,rh,rb,Ghtot
      common/bmb/beta_m
      common/normb/norm

1      format(i1,$)

      norm = bessk(2,beta_m)**2/(2*beta_m)
      lm1 = lm
      ms1 = ms
      rh = (mh/ms)**2
      Ghtot = Gh0
      if (ms.lt.mh/2) Ghtot = Ghtot + (lm*v_higgs)**2/mh/(32*pi)*sqrt(1-4/rh)

      rh = (mh/ms)**2
      if (rh.gt.4) then
         y0 = rh/4
         call qromb(sv_int_new,1.d0,y0,svrel1)
!      write(6,1) 1
         call qromo(sv_int_new,y0,1.d10,svrel2,midinf)
!      write(6,1) 2
      else
         vm = sqrt(2.)/beta_m*sqrt(1+sqrt(1+beta_m**2))
         y0 = 1 + vm**2
         call qromb(sv_int_new,1.d0,y0,svrel1)
!      write(6,1) 3
         call qromo(sv_int_new,y0,1.d10,svrel2,midinf)
!      write(6,1) 4
      endif      

      svrel_new = svrel1 + svrel2      

      end function

!*************************************************************************
      FUNCTION ran3(idum)
!     Returns a uniform random deviate between 0:0 and 1:0. Set idum to any negative value
!     to initialize or reinitialize the sequence.
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
!     According to Knuth, any large mbig, and any smaller (but still large) mseed can be substituted
!     for the above values.
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55) ! The value 55 is special and should not be modied; see
!     REAL mj,mk,ma(55) Knuth.
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then  ! Initialization.
         iff=1
         mj=abs(MSEED-abs(idum)) ! Initialize ma(55) using the seed idum and the large num
         mj=mod(mj,MBIG) ! ber mseed.
         ma(55)=mj
         mk=1
         do 11 i=1,54 ! Now initialize the rest of the table,
            ii=mod(21*i,55) ! in a slightly random order,
            ma(ii)=mk ! with numbers that are not especially random.
            mk=mj-mk
            if(mk.lt.MZ)mk=mk+MBIG
            mj=ma(ii)
11         enddo ! 11
         do 13 k=1,4 ! We randomize them by \warming up the generator."
         do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12         enddo ! 12
13         enddo ! 13
         inext=0 ! Prepare indices for our rst generated number.
         inextp=31 ! The constant 31 is special; see Knuth.
         idum=1
      endif
      inext=inext+1 ! Here is where we start, except on initialization. Increment
      if(inext.eq.56)inext=1 ! inext, wrapping around 56 to 1.
      inextp=inextp+1 ! Ditto for inextp.
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp) ! Now generate a new random number subtractively.
      if(mj.lt.MZ)mj=mj+MBIG ! Be sure that it is in range.
      ma(inext)=mj  ! Store it,
      ran3=mj*FAC  ! and output the derived uniform deviate.
      return
      end function

!*********************************************************************
!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=4000)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
      1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
      u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** &
      2)/6.
      return
      end subroutine


!*********************************************************************
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs, &
      rkqs,odeint_fail)
      logical fail, odeint_fail
      common/rkqsfail/fail
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      PARAMETER (MAXSTP=20000,NMAX=50,KMAXX=2000,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav/0/,dydx(NMAX),xp(KMAXX),y(NMAX), &
      yp(NMAX,KMAXX),yscal(NMAX)
      EXTERNAL rkqs,derivs
      COMMON /path/ kmax,kount,dxsav,xp,yp
      fail=.false.
      odeint_fail=.false.
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if (fail) return
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause &
      'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      write(6,*) 'too many steps in odeint'
      odeint_fail=.true.
      return
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
!U    USES derivs
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX), &
      ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, &
      B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., &
      B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, &
      B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512., &
      B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378., &
      C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648., &
      DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336., &
      DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+ &
      B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6* &
      ak6(i))
17    continue
      return
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      logical fail
      common/rkqsfail/fail
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
!U    USES derivs,rkck
      INTEGER i
      REAL*8 errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK, &
      ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1*h)then
          h=.1*h
        endif
        xnew=x+h
        if(xnew.eq.x)then
            fail=.true.
            return
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE qromo(func,a,b,ss,choose)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.e-7, JMAX=17, JMAXP=JMAX+1, K=5, KM=K-1)
!U    USES polint
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      pause 'too many steps in qromo'
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE midpnt(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)

          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end subroutine

!  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END SUBROUTINE
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
!     USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      return
      write(6,*) 'integral = ',ss
      write(6,*) 'a, b = ',a,b
      call test_int(a,b,func)
      pause 'too many steps in qromb'
      END SUBROUTINE
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      SUBROUTINE midinf(funk,aa,bb,s,n)
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 a,b,ddel,del,sum,tnm,func,x
      func(x)=funk(1./x)/x**2
      b=1./aa
      a=1./bb
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
!        write(6,*) 'a,b = ',a,b
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END SUBROUTINE
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      FUNCTION bessk0(x)
      REAL*8 bessk0,x
!     USES bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,&
      0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,&
      -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      integer bs
      common/bsb/bs
      if (x.le.2.0) then
        y=x*x/4.0
        bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*&
      (p6+y*p7))))))*exp(x*bs)
      else
        y=(2.0/x)
        bessk0=(exp(-x*(1-bs))/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
      q7))))))
      endif
      return
      END FUNCTION
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      FUNCTION bessk1(x)
      REAL*8 bessk1,x
!     USES bessi1
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,&
      -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,&
      0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      integer bs
      common/bsb/bs
      if (x.le.2.0) then
        y=x*x/4.0
        bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*&
      (p5+y*(p6+y*p7))))))*exp(x*bs)
      else
        y=2.0/x
        bessk1=(exp(-x*(1-bs))/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
      q7))))))
      endif
      return
      END FUNCTION
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      FUNCTION bessk(n,x)
      INTEGER n
      REAL*8 bessk,x
!     USES bessk0,bessk1
      INTEGER j
      REAL*8 bk,bkm,bkp,tox
      if (n.lt.2) pause 'bad argument n in bessk'
      tox=2.0/x
      bkm=bessk0(x)
      bk=bessk1(x)
      do 11 j=1,n-1
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
11    continue
      bessk=bk
      return
      END FUNCTION
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      FUNCTION bessi1(x)
      REAL*8 bessi1,x
      REAL*8 ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,&
      0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,&
      -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,&
      0.1787654d-1,-0.420059d-2/
      if (abs(x).lt.3.75) then
        y=(x/3.75)**2
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=3.75/ax
        bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
      (q7+y*(q8+y*q9))))))))
        if(x.lt.0.)bessi1=-bessi1
      endif
      return
      END FUNCTION
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      FUNCTION bessi0(x)
      REAL*8 bessi0,x
      REAL*8 ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,&
      1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,&
      0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,&
      -0.1647633d-1,0.392377d-2/
      if (abs(x).lt.3.75) then
        y=(x/3.75)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=abs(x)
        y=3.75/ax
        bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
      (q7+y*(q8+y*q9))))))))
      endif
      return
      END FUNCTION
!     (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.

end module jimlib
  
