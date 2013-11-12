      real*8 function dsbessei1(x)
c     exp(-|x|) i1(x)
      implicit none
      real*8 x,y,ax
      real*8 p1,p2,p3,p4,p5,p6,p7,
     &  q1,q2,q3,q4,q5,q6,q7,q8,q9
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      data p1,p2,p3,p4,p5,p6,p7/
     &  0.5d0,0.87890594d0,0.51498869d0,0.15084934d0,
     &  0.2658733d-1,0.301532d-2,0.32411d-3/
      data q1,q2,q3,q4,q5,q6,q7,q8,q9/
     &  0.39894228d0,-0.3988024d-1,-0.362018d-2,0.163801d-2,
     &  -0.1031555d-1,0.2282967d-1,-0.2895312d-1,
     &  0.1787654d-1,-0.420059d-2/
      
      ax=abs(x)
      if (ax.lt.3.75d0) then
        y=(x/3.75d0)**2
        dsbessei1=
     &    exp(-ax)*x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        y=3.75d0/ax
        dsbessei1=(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*
     &       (q8+y*q9))))))))/sqrt(ax)
        if (x.lt.0.0d0) dsbessei1=-dsbessei1
      endif
      
      return
      end
