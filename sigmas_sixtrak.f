      subroutine propagate_sigma(sp, bcu, sinth, costh, sx, sy, 
     +costhp, sinthp, ds_sx, ds_sy, Sig_11, Sig_33, Sig_13)
     
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
Cf2py intent(out) Sig_11
Cf2py intent(out) Sig_33
Cf2py intent(out) Sig_13

      implicit none
      
      integer ibbc1, ibbc, jsli, i
      
      
      real*8 sp, pieni, sinth, costh, sfac, sx, sy
      real*8 costhp, sinthp, ds_sx, ds_sy, Sig_11, Sig_33, Sig_13
      real*8 bcu(10), dum(13)
      
      real*8 zero, one, two, half, four
      
      zero= 0.
      one = 1.
      half = 0.5
      four = 4.
      two = 2.
      
      pieni = 1e-10
      ibbc = 1

c~       integer i,ibb,ibbc,ibbc1,ibtyp,jsli,np,nsli
c~       double precision bbf0,bbfx,bbfy,bbgx,bbgy,bcu,costh,costhp,cphi,  
c~      +dum,f,s,sepx,sepx0,sepy,sepy0,sfac,sinth,sinthp,sp,star,sx,       
c~      +sy,track,cphi2

c~       dimension track(6,npart),bcu(nbb,12)
c~       dimension star(3,mbea),dum(13)
!-----------------------------------------------------------------------

      do 2000 jsli=1,1
        do 1000 i=1,1
         write(*,*) 'In!'
          !s=(track(5,i)-star(3,jsli))*half
          !write(*,*)'JBG - cphi2',cphi2
          !sp=s/cphi2 
          dum(1)=(bcu(1)+(two*bcu(4))*sp)+bcu(6)*sp**2      
          dum(2)=(bcu(2)+(two*bcu(9))*sp)+bcu(10)*sp**2     
          dum(3)=(bcu(3)+(bcu(5)+bcu(7))*sp)+               
     +bcu(8)*sp**2      
          Sig_11 = dum(1)
          Sig_33 = dum(2)
          Sig_13 = dum(3)
          dum(4)=dum(1)-dum(2)
          dum(5)=dum(4)**2+four*dum(3)**2 
          write(*,*) 'dum(4)'     
          write(*,*) dum(4)
          write(*,*) 'dum(5)'     
          write(*,*) dum(5)                         
          if(ibbc.eq.1.and.(abs(dum(4)).gt.pieni.and.                   
     +abs(dum(5)).gt.pieni**2)) then
            write(*,*) 'Coupling!'
            ibbc1=1
            dum(5)=sqrt(dum(5))
         else
            ibbc1=0
         endif
        !JBG New set of canonical set of variables at the Col point (CP)
          !sepx0=(track(1,i)+track(2,i)*s)-star(1,jsli)                   !hr06
          !sepy0=(track(3,i)+track(4,i)*s)-star(2,jsli)                   !hr06
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
            write(*,*) 'dum(5)'
            write(*,*) dum(5)
            sy=sfac*dum(5)
            sx=(dum(7)+sy)*half
            sy=(dum(7)-sy)*half
            !sepx=sepx0*costh+sepy0*sinth
            !sepy=sepy0*costh-sepx0*sinth                                 !hr06
          else
            sx=dum(1)
            sy=dum(2)
c~             sepx=sepx0
c~             sepy=sepy0
          endif
          !if(sx.gt.sy) then
          !  call bbf(sepx,sepy,sx,sy,bbfx,bbfy,bbgx,bbgy,ibtyp)
          !else
          !  call bbf(sepy,sepx,sy,sx,bbfy,bbfx,bbgy,bbgx,ibtyp)
          !endif
          !bbfx=f*bbfx
          !bbfy=f*bbfy
          !bbgx=f*bbgx
          !bbgy=f*bbgy
          if(ibbc1.eq.1) then
            dum(8)=two*((bcu(4)-bcu(9))+                        
     +(bcu(6)-bcu(10))*sp)                                      
            dum(9)=(bcu(5)+bcu(7))+(two*bcu(8))*sp          
            dum(10)=(((dum(4)*dum(8)+(four*dum(3))*dum(9))/             
     +dum(5))/dum(5))/dum(5)                                            
            dum(11)=sfac*(dum(8)/dum(5)-dum(4)*dum(10))
            dum(12)=(bcu(4)+bcu(9))+(bcu(6)+bcu(10))*sp 
      dum(13)=(sfac*((dum(4)*dum(8))*half+(two*dum(3))*dum(9)))/dum(5)  
            if(abs(costh).gt.pieni) then
              costhp=(dum(11)/four)/costh                                
            else
              costhp=zero
            endif
            if(abs(sinth).gt.pieni) then
              sinthp=((-1d0*dum(11))/four)/sinth                        
            else
              sinthp=zero
            endif
c~             track(6,i)=track(6,i)-                                      
c~      &((((bbfx*(costhp*sepx0+sinthp*sepy0)+                             
c~      &bbfy*(costhp*sepy0-sinthp*sepx0))+                                
c~      &bbgx*(dum(12)+dum(13)))+bbgy*(dum(12)-dum(13)))/                  
c~      &cphi)*half                                                        
c~             bbf0=bbfx
c~             bbfx=bbf0*costh-bbfy*sinth
c~             bbfy=bbf0*sinth+bbfy*costh
          else
c~             track(6,i)=track(6,i)-                                      &
c~      &(bbgx*(bcu(4)+bcu(6)*sp)+                                 &
c~      &bbgy*(bcu(9)+bcu(10)*sp))/cphi
          endif
c~           track(6,i)=track(6,i)-(bbfx*(track(2,i)-bbfx*half)+           &
c~      &bbfy*(track(4,i)-bbfy*half))*half
c~           track(1,i)=track(1,i)+s*bbfx
c~           track(2,i)=track(2,i)-bbfx
c~           track(3,i)=track(3,i)+s*bbfy
c~           track(4,i)=track(4,i)-bbfy
 1000   continue
 2000 continue
            ds_sx = (dum(12)+dum(13))
            ds_sy = (dum(12)-dum(13))

      end subroutine
