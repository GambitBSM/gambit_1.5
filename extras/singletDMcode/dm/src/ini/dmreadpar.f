      subroutine dsreadpar(unit)
      implicit none
      character line*80,var*80,val*80
      integer unit,i,j,dskillsp,dsival
      real*8 dsfval
      logical dslval
      character*80 dscval
      real*8 msq
      logical m1gut,m3gut
      include 'dssusy.h'
      include 'dsandwcom.h'
      include 'dsrncom.h'
      msq=0.d0
      m1gut=.false.
      m3gut=.false.
c
c... read parameters from file
c
      call dswrite(1,0,'reading darksusy parameters')
 10   continue
      read (unit,'(a)',end=99) line
      if (line(1:1).eq.'!'.or.line(1:1).eq.'#'.or.
     &  len(line).eq.0) goto 10
      j = dskillsp(line,line)
      i = index(line,'=')
      var = line(1:i-1)
      val = line(i+1:j)
      call dslowcase(var)
      call dswrite(2,0,'...setting '//var(1:i-1)//'='//val(1:j-i))
c... gauge coupling constants at the z scale
      if (var.eq.'s2thw') s2thw=dsfval(val)
      if (var.eq.'alphem') alphem=dsfval(val)
      if (var.eq.'alph3') alph3=dsfval(val)
c... standard model masses
      if (var.eq.'mass(kgamma)') mass(kgamma)=dsfval(val)
      if (var.eq.'mass(kgluon)') mass(kgluon)=dsfval(val)
      if (var.eq.'mass(kz)') mass(kz)=dsfval(val)
      if (var.eq.'mass(kw)') mass(kw)=dsfval(val)
      if (var.eq.'mass(knue)') mass(knue)=dsfval(val)
      if (var.eq.'mass(ke)') mass(ke)=dsfval(val)
      if (var.eq.'mass(knumu)') mass(knumu)=dsfval(val)
      if (var.eq.'mass(kmu)') mass(kmu)=dsfval(val)
      if (var.eq.'mass(knutau)') mass(knutau)=dsfval(val)
      if (var.eq.'mass(ktau)') mass(ktau)=dsfval(val)
      if (var.eq.'mass(ku)') mass(ku)=dsfval(val)
      if (var.eq.'mass(kd)') mass(kd)=dsfval(val)
      if (var.eq.'mass(kc)') mass(kc)=dsfval(val)
      if (var.eq.'mass(ks)') mass(ks)=dsfval(val)
      if (var.eq.'mass(kt)') mass(kt)=dsfval(val)
      if (var.eq.'mass(kb)') mass(kb)=dsfval(val)
c... standard model widths
      if (var.eq.'width(kgamma)') width(kgamma)=dsfval(val)
      if (var.eq.'width(kgluon)') width(kgluon)=dsfval(val)
      if (var.eq.'width(kz)') width(kz)=dsfval(val)
      if (var.eq.'width(kw)') width(kw)=dsfval(val)
      if (var.eq.'width(knue)') width(knue)=dsfval(val)
      if (var.eq.'width(ke)') width(ke)=dsfval(val)
      if (var.eq.'width(knumu)') width(knumu)=dsfval(val)
      if (var.eq.'width(kmu)') width(kmu)=dsfval(val)
      if (var.eq.'width(knutau)') width(knutau)=dsfval(val)
      if (var.eq.'width(ktau)') width(ktau)=dsfval(val)
      if (var.eq.'width(ku)') width(ku)=dsfval(val)
      if (var.eq.'width(kd)') width(kd)=dsfval(val)
      if (var.eq.'width(kc)') width(kc)=dsfval(val)
      if (var.eq.'width(ks)') width(ks)=dsfval(val)
      if (var.eq.'width(kt)') width(kt)=dsfval(val)
      if (var.eq.'width(kb)') width(kb)=dsfval(val)
c... cabibbo-kobayashi-maskawa mixing matrix
      if (var.eq.'ckms12') ckms12=dsfval(val)
      if (var.eq.'ckms23') ckms23=dsfval(val)
      if (var.eq.'ckms13') ckms13=dsfval(val)
      if (var.eq.'ckmdelta') ckmdelta=dsfval(val)
c... program switches
      if (var.eq.'higloop') higloop=dsival(val)
      if (var.eq.'prtlevel') prtlevel=dsival(val)
      if (var.eq.'neuloop') neuloop=dsival(val)
      if (var.eq.'bsgqcd') bsgqcd=dsival(val)
c... 1-loop annihilation 
      if (var.eq.'incglue') incglue=dslval(val)
      if (var.eq.'incgaga') incgaga=dslval(val)
      if (var.eq.'incgaz') incgaz=dslval(val)
c... supersymmetric parameters
      if (var.eq.'m1gut') m1gut=dslval(val)
      if (var.eq.'m3gut') m3gut=dslval(val)
      if (var.eq.'mu') mu=dsfval(val)
      if (var.eq.'m1') m1=dsfval(val)
      if (var.eq.'m2') m2=dsfval(val)
      if (var.eq.'m3') m3=dsfval(val)
      if (var.eq.'ma') ma=dsfval(val)
      if (var.eq.'tanbe') tanbe=dsfval(val)
      if (var.eq.'msquarks') msquarks=dsfval(val)
      if (var.eq.'msleptons') msleptons=dsfval(val)
      if (var.eq.'msq') msq=dsfval(val)
      if (var.eq.'mass2q(1)') mass2q(1)=dsfval(val)
      if (var.eq.'mass2u(1)') mass2u(1)=dsfval(val)
      if (var.eq.'mass2l(1)') mass2l(1)=dsfval(val)
      if (var.eq.'mass2e(1)') mass2e(1)=dsfval(val)
      if (var.eq.'mass2q(2)') mass2q(2)=dsfval(val)
      if (var.eq.'mass2u(2)') mass2u(2)=dsfval(val)
      if (var.eq.'mass2l(2)') mass2l(2)=dsfval(val)
      if (var.eq.'mass2e(2)') mass2e(2)=dsfval(val)
      if (var.eq.'mass2q(3)') mass2q(3)=dsfval(val)
      if (var.eq.'mass2u(3)') mass2u(3)=dsfval(val)
      if (var.eq.'mass2l(3)') mass2l(3)=dsfval(val)
      if (var.eq.'mass2e(3)') mass2e(3)=dsfval(val)
      if (var.eq.'asoftu(1)'.or.var.eq.'au') asoftu(1)=dsfval(val)
      if (var.eq.'asoftd(1)'.or.var.eq.'ad') asoftd(1)=dsfval(val)
      if (var.eq.'asoftu(2)'.or.var.eq.'ac') asoftu(2)=dsfval(val)
      if (var.eq.'asoftd(2)'.or.var.eq.'as') asoftd(2)=dsfval(val)
      if (var.eq.'asoftu(3)'.or.var.eq.'at') asoftu(3)=dsfval(val)
      if (var.eq.'asoftd(3)'.or.var.eq.'ab') asoftd(3)=dsfval(val)
c... relic density
      if (var.eq.'omtype') omtype=dsival(val)
      if (var.eq.'omfast') omfast=dsival(val)
c... set-up modules
      if (var.eq.'hmlabel') call dshmset(dscval(val))
      if (var.eq.'pblabel') call dspbset(dscval(val))
      if (var.eq.'ddsilabel') call dsddset('si',dscval(val))
      if (var.eq.'ddsdlabel') call dsddset('sd',dscval(val))
      if (var.eq.'eplabel') call dsepset(dscval(val))
      if (var.eq.'aclabel') call dsacset(dscval(val))
      goto 10
 99   continue
      if (m1gut) m1 = m2 * 5.0d0/3.0d0 * s2thw/(1.d0-s2thw)
      if (m3gut) m3 = m2 * (alph3*s2thw)/alphem
      if (msq.ne.0.d0) then
         mass2q(3) = msq**2
         mass2q(3) = msq**2
         mass2u(3) = msq**2
         mass2d(3) = 2.0d0*mass2q(3)-mass2u(3)
         mass2l(3) = mass2d(3)
         mass2e(3) = mass2d(3)
         do i=1,2
            mass2q(i) = mass2d(3)
            mass2u(i) = mass2d(3)
            mass2d(i) = mass2d(3)
            mass2l(i) = mass2d(3)
            mass2e(i) = mass2d(3)
         enddo
      endif
      return
      end
