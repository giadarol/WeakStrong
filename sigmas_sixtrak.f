      subroutine propagate_sigma(sp, bcu, sinth, costh, sx, sy, 
     +costhp, sinthp, ds_sx, ds_sy)
     
Cf2py intent(in) sp
Cf2py intent(in) bcu
Cf2py intent(out) sinth
Cf2py intent(out) costh
Cf2py intent(out) sx
Cf2py intent(out) sy
Cf2py intent(out) costhp
Cf2py intent(out) sinthp
Cf2py intent(out) ds_sx
Cf2py intent(out) ds_sy

      implicit none
      
      integer ibbc1, ibbc
      
      
      real*8 sp, pieni, sinth, costh, sfac, sx, sy
      real*8 costhp, sinthp, ds_sx, ds_sy
      real*8 bcu(10), dum(13)
      
      pieni = 1e-10
      ibbc = 1

c~       integer i,ibb,ibbc,ibbc1,ibtyp,jsli,np,nsli
c~       double precision bbf0,bbfx,bbfy,bbgx,bbgy,bcu,costh,costhp,cphi,  
c~      +dum,f,s,sepx,sepx0,sepy,sepy0,sfac,sinth,sinthp,sp,star,sx,       
c~      +sy,track,cphi2

c~       dimension track(6,npart),bcu(nbb,12)
c~       dimension star(3,mbea),dum(13)
!-----------------------------------------------------------------------

          dum(1)=(bcu(1)+(2.*bcu(4))*sp)+bcu(6)*sp**2      
          dum(2)=(bcu(2)+(2.*bcu(9))*sp)+bcu(10)*sp**2     
          dum(3)=(bcu(3)+(bcu(5)+bcu(7))*sp)+               
     +bcu(8)*sp**2                                                  
          dum(4)=dum(1)-dum(2)
          dum(5)=dum(4)**2+4.*dum(3)**2                               
          if(ibbc.eq.1.and.(abs(dum(4)).gt.pieni.and.                   
     +abs(dum(5)).gt.pieni)) then
            ibbc1=1
            dum(5)=sqrt(dum(5))
          else
            ibbc1=1 !changed for testing
          endif
          
          if(ibbc1.eq.1) then
            sfac=1.
            if(dum(4).lt.0.) sfac=-1d0                            
            dum(6)=(sfac*dum(4))/dum(5)                                 
            dum(7)=dum(1)+dum(2)
            costh=0.5*(1.+dum(6))
            if(abs(costh).gt.pieni) then
              costh=sqrt(costh)
            else
              costh=0.
            endif
            sinth=0.5*(1.-dum(6))
            if(abs(sinth).gt.pieni) then
              sinth=(-1d0*sfac)*sqrt(sinth)                             
            else
              sinth=0.
            endif
            if(dum(3).lt.0.) sinth=-1d0*sinth                         
            sy=sfac*dum(5)
            sx=(dum(7)+sy)*0.5
            sy=(dum(7)-sy)*0.5                             
          else
            sx=dum(1)
            sy=dum(2)
          endif

          !if(ibbc1.eq.1) then
            dum(8)=2.*((bcu(4)-bcu(9))+                        
     +(bcu(6)-bcu(10))*sp)                                      
            dum(9)=(bcu(5)+bcu(7))+(2.*bcu(8))*sp          
            dum(10)=(((dum(4)*dum(8)+(4.*dum(3))*dum(9))/             
     +dum(5))/dum(5))/dum(5)                                            
            dum(11)=sfac*(dum(8)/dum(5)-dum(4)*dum(10))
            dum(12)=(bcu(4)+bcu(9))+(bcu(6)+bcu(10))*sp 
        dum(13)=(sfac*((dum(4)*dum(8))*0.5+(2.*dum(3))*dum(9)))/dum(5)   
            if(abs(costh).gt.pieni) then
              costhp=(dum(11)/4.)/costh                               
            else
              costhp=0.
            endif
            if(abs(sinth).gt.pieni) then
              sinthp=((-1d0*dum(11))/4.)/sinth                        
            else
              sinthp=0.
            endif
            
            ds_sx = (dum(12)+dum(13))
            ds_sy = (dum(12)-dum(13))

      end subroutine
