*----------------------------------------------------------------------*
* COMPILE USING
*  f2py -m boost_sixtrack -c boost_sixtrack.f 


      subroutine boost(np,sphi,cphi,tphi,salpha,calpha,track)
      
Cf2py intent(in) np      
Cf2py intent(in) sphi
Cf2py intent(in) cphi
Cf2py intent(in) tphi
Cf2py intent(in) salpha
Cf2py intent(in) calpha
Cf2py intent(inout) track

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
      integer i,np
      real*8 calpha,cphi,h,h1x,h1y,h1z,hd1,salpha,sphi,tphi,x1,y1
      real*8 track(6,np)
!-----------------------------------------------------------------------
      do i=1,np
        h=(track(6,i)+1.)-sqrt(((1.+track(6,i))**2
     +    -track(2,i)**2)-track(4,i)**2)                                      
        track(6,i)=((track(6,i)-(calpha*tphi)*track(2,i))  
     +     -(track(4,i)*salpha)*tphi)+h*tphi**2                               
        track(2,i)=(track(2,i)-(tphi*h)*calpha)/cphi                     
        track(4,i)=(track(4,i)-(tphi*h)*salpha)/cphi                     
        hd1=sqrt(((1.+track(6,i))**2-track(2,i)**2)-track(4,i)**2)      
        h1x=track(2,i)/hd1
        h1y=track(4,i)/hd1
        h1z=1.-(1.+track(6,i))/hd1
        x1=((calpha*tphi)*track(5,i)+(1.+(calpha*sphi)*h1x)*track(1,i))
     +    +((track(3,i)*salpha)*sphi)*h1x                                    
        y1=((salpha*tphi)*track(5,i)+(1.+(salpha*sphi)*h1y)*track(3,i))
     +    +((track(1,i)*calpha)*sphi)*h1y                                    
        track(5,i)=track(5,i)/cphi+h1z*((sphi*calpha)*track(1,i)  
     +      +(sphi*salpha)*track(3,i))                                         
        track(1,i)=x1
        track(3,i)=y1
      end do
      end subroutine
      
      
      subroutine boosti(np,sphi,cphi,tphi,salpha,calpha,track)
Cf2py intent(in) np      
Cf2py intent(in) sphi
Cf2py intent(in) cphi
Cf2py intent(in) tphi
Cf2py intent(in) salpha
Cf2py intent(in) calpha
Cf2py intent(inout) track
      
      
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
! BOOSTI **************inverse boost *****************
!-----------------------------------------------------------------------
      implicit none
      integer i,np
      real*8 calpha,cphi,det,h1,h1d,h1x,h1y,h1z,salpha,sphi,tphi,
     +     x1,y1,z1
      real*8 track(6,np)

!-----------------------------------------------------------------------
      do 1000 i=1,np
        h1d=sqrt(((1.+track(6,i))**2-track(2,i)**2)-track(4,i)**2)      
        h1x=track(2,i)/h1d
        h1y=track(4,i)/h1d
        h1z=1.-(1.+track(6,i))/h1d
        h1=((track(6,i)+1.)-sqrt(((1.+track(6,i))**2-                 
     +      track(2,i)**2)-track(4,i)**2))*cphi**2                            
        det=1./cphi+tphi*((h1x*calpha+h1y*salpha)-h1z*sphi)            
        x1=(track(1,i)*(1./cphi+(salpha*(h1y-(h1z*salpha)*sphi))*tphi)
     +      +((track(3,i)*salpha)*tphi)*((h1z*calpha)*sphi-h1x))              
     +      -(track(5,i)*((calpha+((h1y*calpha)*salpha)*sphi)                 
     +      -(h1x*salpha**2)*sphi))*tphi                                      
        y1= (((track(1,i)*calpha)*tphi)*((h1z*salpha)*sphi-h1y)         
     +     +track(3,i)*(1./cphi+(calpha*(h1x-(h1z*calpha)*sphi))*tphi))    
     +      -(track(5,i)*(salpha-(h1y*calpha**2)*sphi                         
     +      +((h1x*calpha)*salpha)*sphi))*tphi                                
        z1= (track(5,i)*((1.+(h1x*calpha)*sphi)+(h1y*salpha)*sphi)     
     +   -((track(1,i)*h1z)*calpha)*sphi)-((track(3,i)*h1z)*salpha)*sphi   
        track(1,i)=x1/det
        track(3,i)=y1/det
        track(5,i)=z1/det
        track(6,i)=(track(6,i)+(calpha*sphi)*track(2,i))                
     +         +(salpha*sphi)*track(4,i)                                         
        track(2,i)=(track(2,i)+(calpha*sphi)*h1)*cphi                   
        track(4,i)=(track(4,i)+(salpha*sphi)*h1)*cphi                   
 1000 continue
      return
      end
