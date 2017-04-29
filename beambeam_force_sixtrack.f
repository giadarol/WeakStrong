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
      
c     subroutine bbf(sepx,sepy,sigxx,sigyy,bbfx,bbfy,bbgx,bbgy,ibtyp)
      subroutine bbf(sepx,sepy,sigxx,sigyy,bbfx,bbfy,bbgx,bbgy)
Cf2py intent(in)  sepx  
Cf2py intent(in)  sepy
Cf2py intent(in)  sigxx
Cf2py intent(in)  sigyy
Cf2py intent(out) bbfx
Cf2py intent(out) bbfy
Cf2py intent(out) bbgx
Cf2py intent(out) bbgy
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

      integer ibtyp
      double precision arg1x,arg1y,arg2x,arg2y,bbfx,bbfy,bbgx,bbgy,     
     + comfac,comfac2,const,expfac,fac,fac2,sepx,sepy,sigxx,sigxy,sigyy, 
     + sqrpi2,wx1,wx2,wy1,wy2,x,xxyy, pieni

      data sqrpi2/3.544907701811032d0/
      save
!-----------------------------------------------------------------------
      
      pieni = 1e-10
      ibtyp = 0
      if(sigxx.eq.sigyy) then
        x=sepx**2+sepy**2
        xxyy=sigxx+sigyy
        const=0.0d0
        if(abs(xxyy).gt.pieni) const=x/xxyy

        expfac=exp(-1d0*const)                                           

        bbfx=0.0d0
        bbfy=0.0d0
        bbgx=0.0d0
        bbgy=0.0d0
        if(abs(x).gt.pieni) then
          bbfx=((2.0d0*sepx)*(1d0-expfac))/x                             
          bbfy=((2.0d0*sepy)*(1d0-expfac))/x                             
          comfac=sepy*bbfy-sepx*bbfx                                     
          comfac2=(abs(sigxx)+abs(sigyy))**2
          bbgx=(comfac+(((4d0*sepx**2)*const)/x)*expfac)/(2d0*x)         
          bbgy=((((4d0*sepy**2)*const)/x)*expfac-comfac)/(2d0*x)         
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
c        if(ibtyp.eq.1) call wzsub(arg1x,arg1y,wy1,wx1)
        if(x.lt.100.d0) then

          expfac=exp(-0.5d0*x)                                           

          arg2x=arg1x/sigxy
          arg2y=arg1y*sigxy
          if(ibtyp.eq.0) call errf(arg2x,arg2y,wy2,wx2)
c          if(ibtyp.eq.1) call wzsub(arg2x,arg2y,wy2,wx2)
          bbfx=const*(wx1-expfac*wx2)
          bbfy=const*(wy1-expfac*wy2)
          if(sepx.lt.0) bbfx=-1d0*bbfx                                   
          if(sepy.lt.0) bbfy=-1d0*bbfy                                   
          comfac=sepx*bbfx+sepy*bbfy
          bbgx=(-1d0*(comfac+2d0*(expfac/sigxy -1d0)))/fac2              
          bbgy= (comfac+2d0*(expfac*sigxy -1d0))/fac2                    
        else
          bbfx=const*wx1
          bbfy=const*wy1
          if(sepx.lt.0) bbfx=-1d0*bbfx                                   
          if(sepy.lt.0) bbfy=-1d0*bbfy                                   
          comfac=sepx*bbfx+sepy*bbfy
          bbgx=(-1d0*(comfac-2d0))/fac2                                  
          bbgy= -1d0*bbgx                                                
        endif
      endif
      return
      end
      
      
