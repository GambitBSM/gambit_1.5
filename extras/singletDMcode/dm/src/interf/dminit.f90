      subroutine dminit(InP, InN, DMIn)

      use parameters

      implicit none

      include'dssusy-f90.h'  
      include'dsidtag-f90.h' 
      include'dsprep-f90.h'
      include 'dshacom-f90.h'

      Type(Input_params), INTENT(IN) :: InP
      Type(Nuisance_params), INTENT(IN) :: InN
      Type(DM), INTENT(IN) :: DMIn

      integer i
 

      include '../ini/dmdir.h'      ! set DarkMatter root directory

!
!... internal fixed-for-ever values go here
!
!  knu=(nue,numu,nutau)   kl=(e,mu,tau)    kqu=(u,c,t)    kqd=(d,s,b)
!  ksqu=(~u1,~c1,~t1,~u2,~c2,~t2)    ksqd=(~d1,~s1,~b1,~d2,~s2,~b2)

!      data pname/'error','nu_e','e','nu_mu','mu','nu_tau','tau','u', &
!      'd','c','s','t','b','gamma','w','z','gluon','h1','h2','a','h+', &
!      's-nu_1','s-l_1','s-l_2','s-nu_2','s-l_3','s-l_4','s-nu_3',  &
!      's-l_5','s-l_6','s-qu_1','s-qu_2','s-qd_1','s-qd_2','s-qu_3',  &
!      's-qu_4','s-qd_3','s-qd_4','s-qu_5','s-qu_6','s-qd_5','s-qd_6', & 
!      'x0_1','x0_2','x0_3','x0_4','x+_1','x+_2','gluino','goldst0',  &
!      'goldst+'/


!... mathematical constants
      pi=3.141592653589793238d0
      gev2cm3s = 0.38937966d-27*3.d10     ! to convert GeV^-2 to cm^3 s^-1

!... standard model charges
      wiso3(knue)   =+0.5d0
      wiso3(ke)     =-0.5d0
      wiso3(knumu)  =+0.5d0
      wiso3(kmu)    =-0.5d0
      wiso3(knutau) =+0.5d0
      wiso3(ktau)   =-0.5d0
      wiso3(ku)     =+0.5d0
      wiso3(kd)     =-0.5d0
      wiso3(kc)     =+0.5d0
      wiso3(ks)     =-0.5d0
      wiso3(kt)     =+0.5d0
      wiso3(kb)     =-0.5d0
      echarg(knue)  =0.d0
      echarg(ke)    =-1.d0
      echarg(knumu) =0.d0
      echarg(kmu)   =-1.d0
      echarg(knutau)=0.d0
      echarg(ktau)  =-1.d0
      echarg(ku)    =+2.d0/3.d0
      echarg(kd)    =-1.d0/3.d0
      echarg(kc)    =+2.d0/3.d0
      echarg(ks)    =-1.d0/3.d0
      echarg(kt)    =+2.d0/3.d0
      echarg(kb)    =-1.d0/3.d0
      ncolor(knue)  =1.d0
      ncolor(ke)    =1.d0
      ncolor(knumu) =1.d0
      ncolor(kmu)   =1.d0
      ncolor(knutau)=1.d0
      ncolor(ktau)  =1.d0
      ncolor(ku)    =3.d0
      ncolor(kd)    =3.d0
      ncolor(kc)    =3.d0
      ncolor(ks)    =3.d0   
      ncolor(kt)    =3.d0
      ncolor(kb)    =3.d0
                     

!
!... default values go here
!
!... gauge coupling constants at the z scale
!....alphaem and alph3mz are nuisance parameters in mcmc

!... standard model masses
      mass(kgamma) =  MassData%gam
      mass(kgluon) =  MassData%g
      mass(kz)     =  MassData%Z
      mass(kw)     =  MassData%W
      mass(knue)   =  MassData%nu1
      mass(knumu)  =  MassData%mu
      mass(knutau) =  MassData%tau
! complete masses
      mass(ku) = MassData%u
      mass(kd) = MassData%d
      mass(ks) = MassData%s 
      mass(kc) = MassData%c
      mass(kb) = MassData%b     
      mass(kt) = InN%mtop 

!... program switches
!... the next 4 flags are related with original spectrum + bsg 
!... computation. Now they are useless....
      prtlevel = 0
!
!      incglue = .true.
!      incgaga = .true.
!      incgaz = .true.
      idtag = ' '
      luout = 6  ! unit where messages go

!... set-up defaults for modules
!... only parts on are init (ie ID stuff off)

      !pass veryting

      mx=InP%mx
       
      !we give sigmav in units of x 10^27
        
      sigmav = InP%sigmav*1d-27
      

!     sigmav = 0.38937966d-27*3.d10*wtot/(2.d0*mx**2) ! in cm^3/s
      wtot = sigmav * 2.d0*mx**2/(0.38937966d-27*3.d10) 

      prtial(12) = wtot * InP%br_zz
      prtial(13) = wtot * InP%br_ww       
      prtial(15) = wtot * InP%br_ee
      prtial(17) = wtot * InP%br_mumu
      prtial(19) = wtot * InP%br_tautau
      prtial(20) = wtot * InP%br_uubar
      prtial(21) = wtot * InP%br_ddbar
      prtial(22) = wtot * InP%br_ccbar       
      prtial(23) = wtot * InP%br_ssbar
      prtial(24) = wtot * InP%br_ttbar
      prtial(25) = wtot * InP%br_bbar
      prtial(26) = wtot * InP%br_gg       
      prtial(28) = wtot * InP%br_gamgam     
      prtial(29) = wtot * InP%br_zgam       
      prtial(30) = wtot * InP%br_hh

      do i=1,30
        sigv(i) = prtial(i)/wtot*sigmav
      enddo
  
      call dmhasetup  
      call dmibset('default')

!...default for new model checks
!      dsprepcalled=.false.

!...Warnings and error init.
 
      suwar = 0
      haerr = 0

      if (Debug) write(*,*) 'Initialization of DM complete.'

      return

      end subroutine dminit


