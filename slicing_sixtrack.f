
      subroutine stsld(star,cphi2,sphi2,sigzs,nsli,calpha,salpha)
      
Cf2py intent(inout) star      
Cf2py intent(in) cphi2 
Cf2py intent(in) sphi2  
Cf2py intent(in) sigzs
Cf2py intent(in) nsli
Cf2py intent(in) calpha
Cf2py intent(in) salpha

!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
!*******STSLD*********************************************************
!   makes longitudinal position of the strong slice for all slices
!*********************************************************************
!-----------------------------------------------------------------------
      implicit none
!+if crlibm
!+ca crlibco
!+ei
      integer i,nsli
      double precision bord,bord1,border,calpha,cphi,cphi2,gauinv,pi,   &
     &salpha,sigz,sigzs,sphi,sphi2,star,yy, half
!+ca parpro
!+ca parnum
      dimension star(3,nsli)
!-----------------------------------------------------------------------
      data border /8d0/
      save
!-----------------------------------------------------------------------

      half = 0.5

!+if crlibm
!      pi=4d0*atan_rn(1d0)
!+ei
!+if .not.crlibm
      pi=4d0*atan(1d0)
!+ei
      sigz=sigzs/cphi2
! DEFINE `STARRED' COORDINATES
!  BORD is longitudinal border star(3,mbea) is the barycenter of region
!  divided two borders.
      bord=+border
      do 101 i=nsli,1,-1
        yy=(1d0/dble(nsli))*dble(i-1)                                    !hr06
        if(i.ne.1) bord1=gauinv(yy)                                      !hr06
        if(i.eq.1) bord1=-1d0*border                                     !hr06
!+if crlibm
!        star(3,i)=(((exp_rn((-1d0*bord**2)*half)-                       &!hr06
!     &exp_rn((-1d0*bord1**2)*half))/sqrt(2d0*pi))*dble(nsli))*sigz       !hr06
!+ei
!+if .not.crlibm
       star(3,i)=(((exp((-1d0*bord**2)*half)-exp((-1d0*bord1**2)*half))/&!hr06
     &sqrt(2d0*pi))*dble(nsli))*sigz                                     !hr06
!+ei
        bord=bord1
        !JBG When doing slicing phi=0 for crab crossing
        ! star(1,i)=0.
        ! star(2,i)=0. 
        !JBG When doing slicing phi2 different tiltings of the strong beam
        star(1,i)=(star(3,i)*sphi2)*calpha
        star(2,i)=(star(3,i)*sphi2)*salpha  
        !star(1,i)=(star(3,i)*sphi)*calpha                                !hr06
        !star(2,i)=(star(3,i)*sphi)*salpha                                !hr06
 101  continue
      return
      end



      function gauinv(p0)
!GAUINV***********************************************
!  INVERSE OF (INTEGRATED) NORMAL DISTRIBUTION FUNCTION
!              1         X= Y
!     P(Y)=-----------* INTEGRAL EXP(-X**2/2) DX
!          SQRT(2*PI)    X= -INF
!     IF P(Y)=P0, THEN GAUINV(P0)=Y.
!        0 < P0 < 1 ,   -INF < Y < +INF
!  IF THIS ROUTINE IS USED TO CONVERT UNIFORM RANDOM NUMBERS TO
!  GAUSSIAN, MAXIMUM RELATIVE ERROR IN THE DISTRIBUTION FUNCTION
!  DP/DX=EXP(-X**2/2)/SQRT(2*PI) IS LESS THAN 0.640E-3 EVERYWHERE
!  IN THE RANGE  2**(-31) < P0 < 1-2**31.  (MINIMAX APPROXIMATION)
      implicit none
!+ca crcoall
!+if crlibm
!+ca crlibco
!+ei
      double precision a0,a1,a2,a3,b0,b1,b2,b3,b4,c0,c1,c2,c3,c4,d0,d1, &
     &d2,d3,d4,e0,e1,e2,e3,e4,f0,f1,f2,gauinv,p,p0,p1,p2,pp1,q,qq2,qq3, &
     &qq4,qq5,t
!-----------------------------------------------------------------------
      data pp1/0.334624883253d0/, qq2/0.090230446775d0/,                &
     &qq3/0.049905685242d0/, qq4/0.027852994157d0/,                     &
     &qq5/0.015645650215d0/
      data a3/ 4.5585614d+01/, a2/ 2.1635544d+00/, a1/ 2.7724523d+00/,  &
     &a0/ 2.5050240d+00/,                                               &
     &b4/ 4.0314354d+02/, b3/-2.7713713d+02/, b2/ 7.9731883d+01/,       &
     &b1/-1.4946512d+01/, b0/ 2.2157257d+00/,                           &
     &c4/ 4.1394487d+03/, c3/-1.5585873d+03/, c2/ 2.4648581d+02/,       &
     &c1/-2.4719139d+01/, c0/ 2.4335936d+00/,                           &
     &d4/ 4.0895693d+04/, d3/-8.5400893d+03/, d2/ 7.4942805d+02/,       &
     &d1/-4.1028898d+01/, d0/ 2.6346872d+00/,                           &
     &e4/ 3.9399134d+05/, e3/-4.6004775d+04/, e2/ 2.2566998d+03/,       &
     &e1/-6.8317697d+01/, e0/ 2.8224654d+00/
      data f0/-8.1807613d-02/, f1/-2.8358733d+00/, f2/ 1.4902469d+00/
      save
!-----------------------------------------------------------------------
      p=p0-0.5d0
      p1=abs(p)
      if(p1.ge.pp1) goto 120
      p2=p**2
      gauinv=(((a3*p2+a2)*p2+a1)*p2+a0)*p
      return
 120  q=0.5d0-p1
      if(q.le.qq2) goto 140
      gauinv=(((b4*q+b3)*q+b2)*q+b1)*q+b0
      goto 200
 140  if(q.le.qq3) goto 150
      gauinv=(((c4*q+c3)*q+c2)*q+c1)*q+c0
      goto 200
 150  if(q.le.qq4) goto 160
      gauinv=(((d4*q+d3)*q+d2)*q+d1)*q+d0
      goto 200
 160  if(q.le.qq5) goto 170
      gauinv=(((e4*q+e3)*q+e2)*q+e1)*q+e0
      goto 200
 170  if(q.le.0d0) goto 900
!+if crlibm
!      t=sqrt(-2d0*log_rn(q))
!+ei
!+if .not.crlibm
      t=sqrt(-2d0*log(q))
!+ei
      gauinv=(t+f0)+f1/(f2+t)                                            !hr06
 200  if(p.lt.0d0) gauinv=-1d0*gauinv                                    !hr06
      return
 900  write(*,910) p0
 910  format(' (FUNC.GAUINV) INVALID INPUT ARGUMENT ',1pd20.13)
!      call closeUnits
!+if cr
!      call abend('                                                  ')
!+ei
!+if .not.cr
!      stop
!+ei
      end
