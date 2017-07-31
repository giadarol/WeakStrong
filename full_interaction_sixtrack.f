c    I take the ERRF from pyecloud

      SUBROUTINE ERRF(XX, YY, WX, WY)
Cf2py intent(in)  XX   
Cf2py intent(in)  YY                                      
Cf2py intent(out) WX   
Cf2py intent(out) WY   
*----------------------------------------------------------------------*   
* Purpose:                                                             *   
*   Modification of WWERF, double precision complex error function,    *   
*   written at CERN by K. Koelbig.                                     *   
* Input:                                                               *   
*   XX, YY    (real)    Argument to CERF.                              *   
* Output:                                                              *   
*   WX, WY    (real)    Function result.                               *   
*----------------------------------------------------------------------*   
                                                                           
*---- Double precision version.                                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)                   
      PARAMETER         (MWFLT = 2, MREAL = 4)                             
      PARAMETER         (MCWRD = 4)                                        
      PARAMETER         (MCNAM = 16, MWNAM = MCNAM / MCWRD)                
      PARAMETER         (MCFIL = 80, MCRNG = 40, MCSTR = 80)               
                                                                           
      PARAMETER         (CC     = 1.12837 91670 9551D0)                    
      PARAMETER         (ONE    = 1.D0)                                    
      PARAMETER         (TWO    = 2.D0)                                    
      PARAMETER         (XLIM   = 5.33D0)                                  
      PARAMETER         (YLIM   = 4.29D0)                                  
      DIMENSION         RX(33), RY(33)                                     
                                                                           
      X = ABS(XX)                                                          
      Y = ABS(YY)                                                          
                                                                           
      IF (Y .LT. YLIM  .AND.  X .LT. XLIM) THEN                            
        Q  = (ONE - Y / YLIM) * SQRT(ONE - (X/XLIM)**2)                    
        H  = ONE / (3.2D0 * Q)                                             
        NC = 7 + INT(23.0*Q)                                               
        XL = H**(1 - NC)                                                   
        XH = Y + 0.5D0/H                                                   
        YH = X                                                             
        NU = 10 + INT(21.0*Q)                                              
        RX(NU+1) = 0.                                                      
        RY(NU+1) = 0.                                                      
                                                                           
        DO 10 N = NU, 1, -1                                                
          TX = XH + N * RX(N+1)                                            
          TY = YH - N * RY(N+1)                                            
          TN = TX*TX + TY*TY                                               
          RX(N) = 0.5D0 * TX / TN                                          
          RY(N) = 0.5D0 * TY / TN                                          
   10   CONTINUE                                                           
                                                                           
        SX = 0.                                                            
        SY = 0.                                                            
                                                                           
        DO 20 N = NC, 1, -1                                                
          SAUX = SX + XL                                                   
          SX = RX(N) * SAUX - RY(N) * SY                                   
          SY = RX(N) * SY + RY(N) * SAUX                                   
          XL = H * XL                                                      
   20   CONTINUE                                                           
                                                                           
        WX = CC * SX                                                       
        WY = CC * SY                                                       
      ELSE                                                                 
        XH = Y                                                             
        YH = X                                                             
        RX(1) = 0.                                                         
        RY(1) = 0.                                                         
                                                                           
        DO 30 N = 9, 1, -1                                                 
          TX = XH + N * RX(1)                                              
          TY = YH - N * RY(1)                                              
          TN = TX*TX + TY*TY                                               
          RX(1) = 0.5D0 * TX / TN                                          
          RY(1) = 0.5D0 * TY / TN                                          
   30   CONTINUE                                                           
                                                                           
        WX = CC * RX(1)                                                    
        WY = CC * RY(1)                                                    
      ENDIF                                                                
                                                                           
      IF(Y .EQ. 0.) WX = EXP(-X**2)                                        
      IF(YY .LT. 0.) THEN                                                  
        WX =   TWO * EXP(Y*Y-X*X) * COS(TWO*X*Y) - WX                      
        WY = - TWO * EXP(Y*Y-X*X) * SIN(TWO*X*Y) - WY                      
        IF(XX .GT. 0.) WY = -WY                                            
      ELSE                                                                 
        IF(XX .LT. 0.) WY = -WY                                            
      ENDIF                                                                
                                                                           
      END
      
      
      subroutine beamint(np,track,param,sigzs,bcu,ibb,ne,ibtyp,ibbc,    &
     &npart, nele, nbb, mbea, beam_expflag, pieni)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
!-----------------------------------------------------------------------
      implicit none
c  +ca crcoall
c  +if crlibm
c  +ca crlibco
c  +ei
      integer ibb,ibbc,ibtyp,ne,np,nsli,npart, nele, nbb, mbea,         &
     &beam_expflag
      double precision alpha,bcu,calpha,cphi,f,param,phi,salpha,sigzs,  
     &sphi,tphi,track,star,phi2,cphi2,sphi2,tphi2, pieni
c  +ca parpro
c  +ca parnum
      dimension track(6,npart)
      dimension param(nele,18),bcu(nbb,12)
      dimension star(3,1)
c  +ca parbeam_exp
      save
!-----------------------------------------------------------------------
      if (beam_expflag .eq. 0) then
         phi=param(ne,1)
         nsli=param(ne,2)
         alpha=param(ne,3)
         f=param(ne,4)/dble(nsli)
         phi2=param(ne,18)
      else if(beam_expflag .eq. 1) then
         alpha=param(ne,3)
         phi=param(ne,1)
         nsli=param(ne,2)
         !sepax=param(ne,4)     !Not actually used anywhere?
         !sepay=param(ne,5)     !Not actually used anywhere?
         f=param(ne,4)/dble(nsli)
         phi2=phi               !Note - phi2 is not a free parameter anymore
      else
         write(*,*) "ERROR in subroutine beamint"
         write(*,*)  "beam_expflag was", beam_expflag
         write(*,*)  " expected 0 or 1. This is a BUG!"
c        call prror(-1)
      endif

c  +if crlibm
c        sphi=sin_rn(phi)
c        sphi2=sin_rn(phi2)
c  +ei
c  +if .not.crlibm
      sphi=sin(phi)
      sphi2=sin(phi2)
c  +ei
c  +if crlibm
c        cphi=cos_rn(phi)
c        cphi2=cos_rn(phi2)
c  +ei
c  +if .not.crlibm
      cphi=cos(phi)
      cphi2=cos(phi2)
c  +ei
c  +if crlibm
c        tphi=tan_rn(phi)
c        tphi2=tan_rn(phi2)
c  +ei
c  +if .not.crlibm
      tphi=tan(phi)
      tphi2=tan(phi2)
c  +ei
c  +if crlibm
c        salpha=sin_rn(alpha)
c  +ei
c  +if .not.crlibm
      salpha=sin(alpha)
c  +ei
c  +if crlibm
c        calpha=cos_rn(alpha)
c  +ei
c  +if .not.crlibm
      calpha=cos(alpha)
c  +ei
!     define slices
c     call stsld(star,cphi2,sphi2,sigzs,nsli,calpha,salpha)
      call stsld(star,cphi2,sphi2,sigzs,nsli,calpha,salpha, mbea)
c     call boost(np,sphi,cphi,tphi,salpha,calpha,track)
      call boost(np,sphi,cphi,tphi,salpha,calpha,track, npart)
c     call sbc(np,star,cphi,cphi2,nsli,f,ibtyp,ibb,bcu,track,ibbc)
      call sbc(np,star,cphi,cphi2,nsli,f,ibtyp,ibb,bcu,track,ibbc,      &
     &mbea,npart,nbb, pieni)
c     call boosti(np,sphi,cphi,tphi,salpha,calpha,track)
      call boosti(np,sphi,cphi,tphi,salpha,calpha,track,npart)
      return
      end
      
      subroutine boost(np,sphi,cphi,tphi,salpha,calpha,track, npart)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
! BOOST Boost Operation ********************************************
!    P,Q,E are all normalized by P0
!-----------------------------------------------------------------------
      implicit none
c  +if crlibm
c  +ca crlibco
c  +ei
      integer i,np, npart
      double precision calpha,cphi,h,h1x,h1y,h1z,hd1,salpha,sphi,tphi,  &
     &track,x1,y1, one, half
c  +ca parpro
c  +ca parnum
      dimension track(6,npart)
      save
      
      one = 1.
      half = 0.5
!-----------------------------------------------------------------------
      do 1000 i=1,np
        h=(track(6,i)+one)-sqrt(((one+track(6,i))**2-                   &!hr06
     &track(2,i)**2)-track(4,i)**2)                                      !hr06
        track(6,i)=((track(6,i)-(calpha*tphi)*track(2,i))               &!hr06
     &-(track(4,i)*salpha)*tphi)+h*tphi**2                               !hr06
        track(2,i)=(track(2,i)-(tphi*h)*calpha)/cphi                     !hr06
        track(4,i)=(track(4,i)-(tphi*h)*salpha)/cphi                     !hr06
        hd1=sqrt(((one+track(6,i))**2-track(2,i)**2)-track(4,i)**2)      !hr06
        h1x=track(2,i)/hd1
        h1y=track(4,i)/hd1
        h1z=one-(one+track(6,i))/hd1
        x1=((calpha*tphi)*track(5,i)+(one+(calpha*sphi)*h1x)*track(1,i))&!hr06
     &+((track(3,i)*salpha)*sphi)*h1x                                    !hr06
        y1=((salpha*tphi)*track(5,i)+(one+(salpha*sphi)*h1y)*track(3,i))&!hr06
     &+((track(1,i)*calpha)*sphi)*h1y                                    !hr06
        track(5,i)=track(5,i)/cphi+h1z*((sphi*calpha)*track(1,i)        &!hr06
     &+(sphi*salpha)*track(3,i))                                         !hr06
        track(1,i)=x1
        track(3,i)=y1
 1000 continue
      return
      end
      subroutine sbc(np,star,cphi,cphi2,nsli,f,ibtyp,ibb,bcu,track,ibbc,&
     &mbea,npart,nbb, pieni)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
!**SBC ***Synchro-Beam for headon collision**********************
!  call BBF  (table) disabled
!****************************************************************
!-----------------------------------------------------------------------
      implicit none
c  +if crlibm
c  +ca crlibco
c  +ei
      integer i,ibb,ibbc,ibbc1,ibtyp,jsli,np,nsli,mbea,npart,nbb
      double precision bbf0,bbfx,bbfy,bbgx,bbgy,bcu,costh,costhp,cphi,  &
     &dum,f,s,sepx,sepx0,sepy,sepy0,sfac,sinth,sinthp,sp,star,sx,       &
     &sy,track,cphi2, four, one, half, zero, two, pieni
c  +ca parpro
c  +ca parnum
      dimension track(6,npart),bcu(nbb,12)
      dimension star(3,mbea),dum(13)
      save
      
      four = 4.
      one = 1.
      half = .5
      zero = 0.
      two = 2.
!-----------------------------------------------------------------------
      do 2000 jsli=1,nsli
        do 1000 i=1,np
          s=(track(5,i)-star(3,jsli))*half
          !write(*,*)'JBG - cphi2',cphi2
          sp=s/cphi2 
          dum(1)=(bcu(ibb,1)+(two*bcu(ibb,4))*sp)+bcu(ibb,6)*sp**2       !hr06
          dum(2)=(bcu(ibb,2)+(two*bcu(ibb,9))*sp)+bcu(ibb,10)*sp**2      !hr06
          dum(3)=(bcu(ibb,3)+(bcu(ibb,5)+bcu(ibb,7))*sp)+               &!hr06
     &bcu(ibb,8)*sp**2                                                   !hr06
          dum(4)=dum(1)-dum(2)
          dum(5)=dum(4)**2+four*dum(3)**2                                !hr06
          if(ibbc.eq.1.and.(abs(dum(4)).gt.pieni.and.                   &
     &abs(dum(5)).gt.pieni)) then
            ibbc1=1
            dum(5)=sqrt(dum(5))
         else
            ibbc1=0
          endif
        !JBG New set of canonical set of variables at the Col point (CP)
          sepx0=(track(1,i)+track(2,i)*s)-star(1,jsli)                   !hr06
          sepy0=(track(3,i)+track(4,i)*s)-star(2,jsli)                   !hr06
          if(ibbc1.eq.1) then
            sfac=one
            if(dum(4).lt.zero) sfac=-1d0*one                             !hr06
            dum(6)=(sfac*dum(4))/dum(5)                                  !hr06
            dum(7)=dum(1)+dum(2)
            costh=half*(one+dum(6))
            if(abs(costh).gt.pieni) then
              costh=sqrt(costh)
            else
              costh=zero
            endif
            sinth=half*(one-dum(6))
            if(abs(sinth).gt.pieni) then
              sinth=(-1d0*sfac)*sqrt(sinth)                              !hr06
            else
              sinth=zero
            endif
            if(dum(3).lt.zero) sinth=-1d0*sinth                          !hr06
            sy=sfac*dum(5)
            sx=(dum(7)+sy)*half
            sy=(dum(7)-sy)*half
            sepx=sepx0*costh+sepy0*sinth
            sepy=sepy0*costh-sepx0*sinth                                 !hr06
          else
            sx=dum(1)
            sy=dum(2)
            sepx=sepx0
            sepy=sepy0
          endif
          if(sx.gt.sy) then
            call bbf(sepx,sepy,sx,sy,bbfx,bbfy,bbgx,bbgy,ibtyp, pieni)
          else
            call bbf(sepy,sepx,sy,sx,bbfy,bbfx,bbgy,bbgx,ibtyp, pieni)
          endif
          bbfx=f*bbfx
          bbfy=f*bbfy
          bbgx=f*bbgx
          bbgy=f*bbgy
          if(ibbc1.eq.1) then
            dum(8)=two*((bcu(ibb,4)-bcu(ibb,9))+                        &!hr06
     &(bcu(ibb,6)-bcu(ibb,10))*sp)                                       !hr06
            dum(9)=(bcu(ibb,5)+bcu(ibb,7))+(two*bcu(ibb,8))*sp           !hr06
            dum(10)=(((dum(4)*dum(8)+(four*dum(3))*dum(9))/             &!hr06
     &dum(5))/dum(5))/dum(5)                                             !hr06
            dum(11)=sfac*(dum(8)/dum(5)-dum(4)*dum(10))
            dum(12)=(bcu(ibb,4)+bcu(ibb,9))+(bcu(ibb,6)+bcu(ibb,10))*sp  !hr06
      dum(13)=(sfac*((dum(4)*dum(8))*half+(two*dum(3))*dum(9)))/dum(5)   !hr06
            if(abs(costh).gt.pieni) then
              costhp=(dum(11)/four)/costh                                !hr06
            else
              costhp=zero
            endif
            if(abs(sinth).gt.pieni) then
              sinthp=((-1d0*dum(11))/four)/sinth                         !hr06
            else
              sinthp=zero
            endif
            track(6,i)=track(6,i)-                                      &!hr06
     &((((bbfx*(costhp*sepx0+sinthp*sepy0)+                             &!hr06
     &bbfy*(costhp*sepy0-sinthp*sepx0))+                                &!hr06
     &bbgx*(dum(12)+dum(13)))+bbgy*(dum(12)-dum(13)))/                  &!hr06
     &cphi)*half                                                         !hr06
            bbf0=bbfx
            bbfx=bbf0*costh-bbfy*sinth
            bbfy=bbf0*sinth+bbfy*costh
          else
            track(6,i)=track(6,i)-                                      &
     &(bbgx*(bcu(ibb,4)+bcu(ibb,6)*sp)+                                 &
     &bbgy*(bcu(ibb,9)+bcu(ibb,10)*sp))/cphi
          endif
          track(6,i)=track(6,i)-(bbfx*(track(2,i)-bbfx*half)+           &
     &bbfy*(track(4,i)-bbfy*half))*half
          track(1,i)=track(1,i)+s*bbfx
          track(2,i)=track(2,i)-bbfx
          track(3,i)=track(3,i)+s*bbfy
          track(4,i)=track(4,i)-bbfy
 1000   continue
 2000 continue
      return
      end
      subroutine boosti(np,sphi,cphi,tphi,salpha,calpha,track,npart)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
! BOOSTI **************inverse boost *****************
!-----------------------------------------------------------------------
      implicit none
c  +if crlibm
c  +ca crlibco
c  +ei
      integer i,np,npart
      double precision calpha,cphi,det,h1,h1d,h1x,h1y,h1z,salpha,sphi,  &
     &tphi,track,x1,y1,z1, one
c  +ca parpro
c  +ca parnum
      dimension track(6,npart)
      save
      
      one = 1.
!-----------------------------------------------------------------------
      do 1000 i=1,np
        h1d=sqrt(((one+track(6,i))**2-track(2,i)**2)-track(4,i)**2)      !hr06
        h1x=track(2,i)/h1d
        h1y=track(4,i)/h1d
        h1z=one-(one+track(6,i))/h1d
        h1=((track(6,i)+one)-sqrt(((one+track(6,i))**2-                 &!hr06
     &track(2,i)**2)-track(4,i)**2))*cphi**2                             !hr06
        det=one/cphi+tphi*((h1x*calpha+h1y*salpha)-h1z*sphi)             !hr06
        x1= (track(1,i)*(one/cphi+(salpha*(h1y-(h1z*salpha)*sphi))*tphi)&!hr06
     &+((track(3,i)*salpha)*tphi)*((h1z*calpha)*sphi-h1x))              &!hr06
     &-(track(5,i)*((calpha+((h1y*calpha)*salpha)*sphi)                 &!hr06
     &-(h1x*salpha**2)*sphi))*tphi                                       !hr06
        y1= (((track(1,i)*calpha)*tphi)*((h1z*salpha)*sphi-h1y)         &!hr06
     &+track(3,i)*(one/cphi+(calpha*(h1x-(h1z*calpha)*sphi))*tphi))     &!hr06
     &-(track(5,i)*(salpha-(h1y*calpha**2)*sphi                         &!hr06
     &+((h1x*calpha)*salpha)*sphi))*tphi                                 !hr06
        z1= (track(5,i)*((one+(h1x*calpha)*sphi)+(h1y*salpha)*sphi)     &!hr06
     &-((track(1,i)*h1z)*calpha)*sphi)-((track(3,i)*h1z)*salpha)*sphi    !hr06
        track(1,i)=x1/det
        track(3,i)=y1/det
        track(5,i)=z1/det
        track(6,i)=(track(6,i)+(calpha*sphi)*track(2,i))                &!hr06
     &+(salpha*sphi)*track(4,i)                                          !hr06
        track(2,i)=(track(2,i)+(calpha*sphi)*h1)*cphi                    !hr06
        track(4,i)=(track(4,i)+(salpha*sphi)*h1)*cphi                    !hr06
 1000 continue
      return
      end
      subroutine bbf(sepx,sepy,sigxx,sigyy,bbfx,bbfy,bbgx,bbgy,ibtyp,   &
     &pieni)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
!-----------------------------------------------------------------------
!**BBF   without using table ******************************************
! gives transverse (f_x and f_y) and longitudinal(g_x and g_y)
! beam-beam kicks except for the kinematical term (nr_e/\gamma)
! SIGXX is \Sigma
!**********************************************************************
      implicit none
c  +if crlibm
c  +ca crlibco
c  +ei
      integer ibtyp
      double precision arg1x,arg1y,arg2x,arg2y,bbfx,bbfy,bbgx,bbgy,     &
     &comfac,comfac2,const,expfac,fac,fac2,sepx,sepy,sigxx,sigxy,sigyy, &
     &sqrpi2,wx1,wx2,wy1,wy2,x,xxyy, pieni
c  +ca parpro
c  +ca parnum
      data sqrpi2/3.544907701811032d0/
      save
!-----------------------------------------------------------------------
      if(sigxx.eq.sigyy) then
        x=sepx**2+sepy**2
        xxyy=sigxx+sigyy
        const=0.0d0
        if(abs(xxyy).gt.pieni) const=x/xxyy
c  +if crlibm
c          expfac=exp_rn(-1d0*const)                                        !hr06
c  +ei
c  +if .not.crlibm
        expfac=exp(-1d0*const)                                           !hr06
c  +ei
        bbfx=0.0d0
        bbfy=0.0d0
        bbgx=0.0d0
        bbgy=0.0d0
        if(abs(x).gt.pieni) then
          bbfx=((2.0d0*sepx)*(1d0-expfac))/x                             !hr06
          bbfy=((2.0d0*sepy)*(1d0-expfac))/x                             !hr06
          comfac=sepy*bbfy-sepx*bbfx                                     !hr06
          comfac2=(abs(sigxx)+abs(sigyy))**2
          bbgx=(comfac+(((4d0*sepx**2)*const)/x)*expfac)/(2d0*x)         !hr06
          bbgy=((((4d0*sepy**2)*const)/x)*expfac-comfac)/(2d0*x)         !hr06
        endif
      else
        x=sepx**2/sigxx+sepy**2/sigyy
        fac2=2.d0*abs(sigxx-sigyy)
        fac=sqrt(fac2)
        const=sqrpi2/fac
        sigxy=sqrt(sigxx/sigyy)
        arg1x=abs(sepx/fac)
        arg1y=abs(sepy/fac)
        if(ibtyp.eq.0) call errf(arg1x,arg1y,wy1,wx1)
c~         if(ibtyp.eq.1) call wzsub(arg1x,arg1y,wy1,wx1)
        if(x.lt.100.d0) then
c  +if crlibm
c            expfac=exp_rn(-0.5d0*x)                                        !hr06
c  +ei
c  +if .not.crlibm
          expfac=exp(-0.5d0*x)                                           !hr06
c  +ei
          arg2x=arg1x/sigxy
          arg2y=arg1y*sigxy
          if(ibtyp.eq.0) call errf(arg2x,arg2y,wy2,wx2)
c~           if(ibtyp.eq.1) call wzsub(arg2x,arg2y,wy2,wx2)
          bbfx=const*(wx1-expfac*wx2)
          bbfy=const*(wy1-expfac*wy2)
          if(sepx.lt.0) bbfx=-1d0*bbfx                                   !hr06
          if(sepy.lt.0) bbfy=-1d0*bbfy                                   !hr06
          comfac=sepx*bbfx+sepy*bbfy
          bbgx=(-1d0*(comfac+2d0*(expfac/sigxy -1d0)))/fac2              !hr06
          bbgy= (comfac+2d0*(expfac*sigxy -1d0))/fac2                    !hr06
        else
          bbfx=const*wx1
          bbfy=const*wy1
          if(sepx.lt.0) bbfx=-1d0*bbfx                                   !hr06
          if(sepy.lt.0) bbfy=-1d0*bbfy                                   !hr06
          comfac=sepx*bbfx+sepy*bbfy
          bbgx=(-1d0*(comfac-2d0))/fac2                                  !hr06
          bbgy= -1d0*bbgx                                                !hr06
        endif
      endif
      return
      end
      subroutine stsld(star,cphi2,sphi2,sigzs,nsli,calpha,salpha, mbea)
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
c  +if crlibm
c  +ca crlibco
c  +ei
      integer i,nsli, mbea
      double precision bord,bord1,border,calpha,cphi,cphi2,gauinv,pi,   &
     &salpha,sigz,sigzs,sphi,sphi2,star,yy, half
c  +ca parpro
c  +ca parnum
      dimension star(3,mbea)
!-----------------------------------------------------------------------
      data border /8d0/
      save
      
      half = 0.5
!-----------------------------------------------------------------------
c  +if crlibm
c        pi=4d0*atan_rn(1d0)
c  +ei
c  +if .not.crlibm
      pi=4d0*atan(1d0)
c  +ei
      sigz=sigzs/cphi2
! DEFINE `STARRED' COORDINATES
!  BORD is longitudinal border star(3,mbea) is the barycenter of region
!  divided two borders.
      bord=+border
      do 101 i=nsli,1,-1
        yy=(1d0/dble(nsli))*dble(i-1)                                    !hr06
        if(i.ne.1) bord1=gauinv(yy)                                      !hr06
        if(i.eq.1) bord1=-1d0*border                                     !hr06
c  +if crlibm
c          star(3,i)=(((exp_rn((-1d0*bord**2)*half)-                       &!hr06
c       &exp_rn((-1d0*bord1**2)*half))/sqrt(2d0*pi))*dble(nsli))*sigz       !hr06
c  +ei
c  +if .not.crlibm

       star(3,i)=(((exp((-1d0*bord**2)*half)-exp((-1d0*bord1**2)*half))/&!hr06
     &sqrt(2d0*pi))*dble(nsli))*sigz                                     !hr06                                !hr06
c  +ei
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
c  +ca crcoall
c  +if crlibm
c  +ca crlibco
c  +ei
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
c  +if crlibm
c        t=sqrt(-2d0*log_rn(q))
c  +ei
c  +if .not.crlibm
      t=sqrt(-2d0*log(q))
c  +ei
      gauinv=(t+f0)+f1/(f2+t)                                            !hr06
 200  if(p.lt.0d0) gauinv=-1d0*gauinv                                    !hr06
      return
 900  write(*,910) p0
 910  format(' (FUNC.GAUINV) INVALID INPUT ARGUMENT ',1pd20.13)
      
      end
      
