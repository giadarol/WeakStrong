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
