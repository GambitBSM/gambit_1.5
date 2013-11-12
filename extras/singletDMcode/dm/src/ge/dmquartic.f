      subroutine dsquartic(a3,a2,a1,a0,z1,z2,z3,z4)
c_______________________________________________________________________
c  analytic solution of z^4 + a3 z^3 + a2 z^2 + a1 z + a0 = 0
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      real*8 a3,a2,a1,a0,z1,z2,z3,z4
      real*8 amax,b4,b3,b2,b1,b0,aa,bb,cc,rr,qq,theta,u,s1,s2,
     &  ap,p1,p2,aq,q1,q2,pi23,temp
c     2/3 pi
      parameter (pi23=2.0d0/3.0d0*3.1415926535897932384626d0)
      amax = 1.0d0
c      write(*,*) 'a3 = ',a3
c      write(*,*) 'a2 = ',a2
c      write(*,*) 'a1 = ',a1
c      write(*,*) 'a0 = ',a0
      if (abs(a3).gt.amax) amax=abs(a3)
      if (abs(a2).gt.amax) amax=abs(a2)
      if (abs(a1).gt.amax) amax=abs(a1)
      if (abs(a0).gt.amax) amax=abs(a0)
      b4 = 1.0d0/amax
      b3 = a3/amax
      b2 = a2/amax
      b1 = a1/amax
      b0 = a0/amax
      aa = -b2
      bb = b1*b3-4.0d0*b0*b4
      cc = 4.0d0*b0*b2*b4-b0*b3**2-b1**2*b4
      rr = (2.0d0*aa**3-9.0d0*aa*bb+27.0d0*cc)/54.0d0
      qq = (aa**2-3.0d0*bb)/9.0d0
      if (rr**2.le.qq**3) then
        temp = sqrt(qq**3)
        if (temp.eq.0.0d0) then
          theta = 0.0d0
        else
          theta = acos(rr/temp)/3.0d0
        endif
        u = -2.0d0*sqrt(qq)*cos(theta+pi23)-aa/3.0d0
        if (.not.( (u**2.ge.4.0d0*b0*b4) .and.
     &             (b3**2.ge.4.0d0*b4*(b2-u)) )) then
          u = -2.0d0*sqrt(qq)*cos(theta)-aa/3.0d0
        endif
        if (.not.( (u**2.ge.4.0d0*b0*b4) .and.
     &             (b3**2.ge.4.0d0*b4*(b2-u)) )) then
          u = -2.0d0*sqrt(qq)*cos(theta-pi23)-aa/3.0d0
        endif
        if (.not.( (u**2.ge.4.0d0*b0*b4) .and.
     &             (b3**2.ge.4.0d0*b4*(b2-u)) )) then
c ruiz            call dswrite(0,1,'suquartic: never get here')
          return
        endif
      else
        s1 = - sign( (abs(rr)+sqrt(rr**2-qq**3))**(1.0d0/3.0d0), rr)
        if (aa.eq.0.0d0) then
          s2 = 0.0d0
        else
          s2 = qq/aa
        endif
        u = (s1+s2)-aa/3.0d0
      endif
      ap = 0.5d0 * (b3 + sign( sqrt(b3**2-4.0d0*b4*(b2-u)), b3))
      p1 = ap/b4
      if (ap.eq.0.0d0) then
        p2 = 0.0d0
      else
        p2 = (b2-u)/ap
      endif
      aq = 0.5d0 * (u + sign( sqrt(u**2-4.0d0*b0*b4), u))
      if (sign(1.0d0,b3*u-2.0d0*b1*b4) .eq. sign(1.0d0,b3*u)) then
        q1 = aq/b4
        if (aq.eq.0.0d0) then
          q2 = 0.0d0
        else
          q2 = b0/aq
        endif
      else
        q2 = aq/b4
        if (aq.eq.0.0d0) then
          q1 = 0.0d0
        else
          q1 = b0/aq
        endif
      endif
      if (p1.eq.0.0d0) then
        z1 = -sqrt(-q1)
      else
        z1 = -0.5d0 * (p1 + sign( sqrt(p1**2-4.0d0*q1), p1))
      endif
      if (p2.eq.0.0d0) then
        z2 = -sqrt(-q2)
      else
        z2 = -0.5d0 * (p2 + sign( sqrt(p2**2-4.0d0*q2), p2))
      endif
      if (z1.eq.0.0d0) then
        z3 = 0.0d0
      else
        z3 = q1/z1
      endif
      if (z2.eq.0.0d0) then
        z4 = 0.0d0
      else
        z4 = q2/z2
      endif
      end
