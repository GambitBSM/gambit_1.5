*      -*- mode: fortran  -*-
      real*8 lpi,mw,mz,g2,g1,e,th,mmtau
      common/dsstauconst/lpi,mw,mz,g2,g1,e,th,mmtau

      REAL*8 mH1,mH2,mstau,mstn,msel,msnu,tb,al,mst,msn,mstau2,
     &  lma,atau,lmu,tht,mHp,msel2
      COMMON/dsstaupar/mH1,mH2,mstau,mstn,msel,msnu,tb,al,mst,msn,
     &     mstau2,lma,Atau,lmu,tht,mHp,msel2

      REAL*8 nevecs(4,4),nevals(4),mx(2),mp(2,2),mm(2,2)
      COMMON /dsinomix/nevecs,mm,mp
      COMMON /dsinomass/nevals,mx

      save /dsstauconst/,/dsstaupar/,/dsinomix/,/dsinomass/
