      program readoutput 
!********************************************************************   
!                                                                   *   
!  This program reads the unformatted file and prints readable      *   
!  numbers.                                                         *   
!                                                                   *   
!********************************************************************   
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension encm(0:idim1) 
      dimension x(0:idim1), v(0:idim1),a(0:idim1) 
      dimension q(idim), dq(idim), u(idim), ifleos(idim) 
      dimension deltam(idim), abar(idim), rho(idim) 
      dimension temp(idim), ye(idim), xp(idim), xn(idim) 
      dimension ynue(idim), ynueb(idim), ynux(idim) 
      dimension unue(idim), unueb(idim), unux(idim) 
      dimension ufreez(idim), pr(idim), u2(idim) 
      real ycc(idim,19) 
      dimension eta(idim) 
      dimension vsound(idim) 
      dimension f2ynue(idim),f2ynueb(idim),f2ynux(idim) 
      dimension f2unue(idim),f2unueb(idim),f2unux(idim) 
      dimension etanue(idim),etanueb(idim),etanux(idim) 
      dimension tempnue(idim), tempnueb(idim), tempnux(idim) 
      dimension dj(idim), steps(idim)
      logical te(idim), teb(idim), tx(idim) 
      character*1 sample,again
      character*1024 output,basename 
      character*32 dumpn 
      character(:), allocatable :: outname
      character*1024 infile 
      integer iskip, ndump, idump
      logical from_dump
!                                                                       
      double precision dm,press,enue,enueb,ks,ka,ksb,kab,               &
     &     ksx                                                          
!   
       if (command_argument_count() == 1) then
           print*, 'ERROR: Wrong number of arguments: [Input, Output, #Dumps]'
           call EXIT(0)
       else if (command_argument_count() == 2) then
           print*, 'ERROR: Wrong number of arguments: [Input, Output, #Dumps]'
           call EXIT(0)
       else if (command_argument_count() == 3) then
           call get_command_argument(number=1, value=infile)
           call get_command_argument(number=2, value=basename)
           call get_command_argument(number=3, value=dumpn)
           read(dumpn,*)  ndump
       else
          print*,' ==============================================='             
          print*, '  Reading from setup_readout'
          print*,' ==============================================='                       
          open(11,file='setup_readout') 
          read(11,*)
          read(11,10) infile
          read(11,*)
          read(11,*)
          read(11,10) basename   
          read(11,*)
          read(11,*)
          read(11,*) ndump 
   10 format(A) 
       endif
   
   
      !print*, "What's the data file name?"
      !read *, infile 
      open(42,file=infile,form='unformatted') 
!                                                                       
      pi43=3.14159265359*4.0/3.0 
   97 continue 
      !print *,'number of dumps?' 
      !read *, ndump 
!                                                                       
!--read data                                                            
!                                                                       
      !print *, 'output basename' 
      !read(*,fmt='(a)') basename 
      ibasenamelen = index(basename,' ')-1 
                                                                        
      do k=1,ndump         
         read(42) idump,nc,t,xmcore,rb,ftrape,ftrapb,ftrapx,             &
               shock_ind,shock_x,from_dump,rlumnue,rlumnueb,rlumnux,     &
               (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc), &
               (u(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),        &
               (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),          &
               (xp(i),i=1,nc),(xn(i),i=1,nc),(ifleos(i),i=1,nc),         &
               (ynue(i),i=1,nc),(ynueb(i),i=1,nc),(ynux(i),i=1,nc),      &
               (unue(i),i=1,nc),(unueb(i),i=1,nc),(unux(i),i=1,nc),      &
               (ufreez(i),i=1,nc),(pr(i),i=1,nc),(u2(i),i=1,nc),         &
               (dj(i),i=1,nc),                                           &
               (te(i),i=1,nc),(teb(i),i=1,nc),(tx(i),i=1,nc),            &
               (steps(i),i=1,nc),((ycc(i,j),j=1,nqn),i=1,nc)                  
!             (vsound(i),i=1,nc)                                        
!        
         !print*, 'rho(1)  ftrape  ftrapb  ftrapx'
         !print*, rho(1),ftrape,ftrapb,ftrapx 
         if (k.lt.40) then 
            imp=1 
!         else                                                          
!            imp=10                                                     
         end if 
         if (mod(k,imp).eq.0) then 
!            k1=36+int(k/imp)+1                                         
            k1=int(k/imp)+1 
!            print *, 'k1 ',k1 
            rhomax=0. 
!            print *, 'xmcore,t',xmcore,t 
            encm(0)=xmcore 
!xmcore-0.4107                                                          
            dke=0. 
            do i=1,nc 
               encm(i)=encm(i-1)+deltam(i) 
!pi43*(x(i)**3-x(i-1)**3)*rho(i)                                        
               vesc=dsqrt(2.*13.34*encm(i-1)/x(i)) 
               if ((v(i)-vesc).gt.0) then 
                  dke=dke+0.5*deltam(i)*(v(i)-vesc)**2 
               end if 
            end do 
!            print *,'dke ', dke 
            vmin=0. 
            do i=1,nc 
               if (v(i).lt.vmin) then 
                  vmin=v(i) 
                  xshock=x(i) 
               end if 
               rhomax=max(rho(i),rhomax) 
            end do 
!            write(72,103) t,vmin,xshock,rho(1),rhomax                              
            !write(dumpn,*) k1-1
            write(dumpn,*) idump
            allocate(character(LEN(TRIM(basename))+LEN(adjustl(dumpn))) :: outname)
            outname=basename(1:ibasenamelen)//'.'//adjustl(dumpn)
            open(69,file=outname) 
            
            tautot=0. 
            dk=0. 
            dene=0. 
            sumni=0. 
            sumti=0. 
            sumca=0. 
            sumsi=0. 
            sumc=0. 
            sumo=0. 
            sumne=0. 
            summg=0. 
            sums=0. 
            sumar=0. 
            sumcr=0. 
            sumzn=0. 
            sumfe=0. 
            iskip=0 
            write(69,*) 'Time [s]  R_shock [index]  R_shock [cm] nu_e_flux [foe/s]'
            write(69,108)10.d0*t, int(shock_ind), 1.d9*shock_x, 2.d-3*rlumnue
            write(69,*)'Cell  M_enclosed [M_sol]  Radius [cm]  Rho [g/cm^3]  Velocity [cm/s] &
                        & Ye  Pressure [g/cm/s^2] Temperature [K]'            
            do i=1,nc 
!               write(69,103)i,encm(i),x(i),rho(i),v(i),ye(i),          
!     $           vsound(i)                                             
!                                                                       
!eta(i),                                                                
!     $              temp(i),u2(i),u(i),ynue(i),ynueb(i),ynux(i)        
               if (encm(i).lt.199.) then 
                  if (i.eq.1) then 
                     dm=encm(1) 
                  else 
                     dm=encm(i)-encm(i-1) 
                  end if 
                  press=pr(i)*2.d22*dm 
                  if (ynue(i).gt.0) then 
                     enue=unue(i)/ynue(i)/96.44 
                  else 
                     enue=10.83 
                  end if 
                  if (ynueb(i).gt.0) then 
                     enueb=unueb(i)/ynueb(i)/96.44 
                  else 
                     enueb=16.8 
                  end if 
                  if (ynux(i).gt.0) then 
                     enux=unux(i)/ynux(i)/96.44 
                  else 
                     enux=25. 
                  end if 
                  ksx=0.2899/1.5799*1.512139d-20*enux**2*               &
     &                 2.d6*rho(i)                                      
                  ks=1.512139d-20*enue**2*2.d6*rho(i) 
                  ka=5.84973737d-20*enue**2*2.d6*rho(i) 
                  ksb=1.512139d-20*enueb**2*2.d6                        &
     &                 *rho(i)                                          
                  kab=5.84973737d-20*enueb**2*2.d6                      &
     &                 *rho(i)                                          
                  if (i.eq.1) then 
                     tautot=x(i)*(ks+ka) 
                  else 
                     tautot=tautot+(x(i)-x(i-1))*(ks+ka) 
                  end if 
                  vg=dsqrt(13.34*(encm(i))/x(i)) 
                  if (v(i).gt.vg) then 
                     dk=dk+0.5*(deltam(i))*                             &
     &                    (v(i)**2-vg**2)                               
                     if (ifleos(i).eq.1) then 
                        dene=dene+0.5*(deltam(i))*                      &
     &                       (v(i)**2-vg**2+2.*u(i))                    
                        if (u(i).lt.0.) then 
!                           print *, u(i),ifleos(i)                     
!                           stop                                        
                        end if 
                     elseif (ifleos(i).eq.2) then 
                        ui=u(i)+860. 
                        if (ui.lt.0.) then 
!                           print *, ui,ifleos(i)                       
!                           stop                                        
                        end if 
                        dene=dene+0.5*(deltam(i))*                      &
     &                       (v(i)**2-vg**2+2.*ui)                      
                     end if 
                  end if 
!                  if (i.gt.700) then                                   
!                     if (i.eq.701) write(69,*) t*10.                   
!                  if (encm(i).gt.1.512.and.encm(i).lt.11.0) then     
                  if (encm(i).gt.0..and.encm(i).lt.11.0) then
! Pressure: converting the units, you get a 2.d22 factor, so where did 1.d16 come from?

                     write(69,103)i,encm(i),1.d9*x(i),2.d6*rho(i),      &
!     &                    1.d8*v(i),ye(i),1.d16*pr(i)!,                  &
     &                    1.d8*v(i),ye(i),2.d22*pr(i),                  &
     &                    1.d9*temp(i)     
!     &                    1.d8*vsound(i)                                
!                  if (i.gt.1) then                                     
                     if (ufreez(i).lt.1.d-10) then 
                       dene2=dene2+                                     &
     &                       deltam(i)*(1.d16*(u(i))+v(i)**2)           
!                       write(69,105) encm(i),1.d9*x(i),                
!     $                       2.d6*rho(i),                              
!     $                       1.d16*u(i),                               
!$                       1.d16*(u(i)+860.),                             
!     $                    1.d16*pr(i),1.d8*v(i),1.d9*temp(i),ye(i),    
!     $                      5.165268d8*dsqrt(encm(i)/x(i))             
                     else 
                       dene2=dene2+                                     &
     &                       deltam(i)*(1.d16*u(i)-                     &
     &                       1.d14*ufreez(i)+v(i)**2-vg**2)             
!                        write(69,105) encm(i),1.d9*x(i),               
!     $                       2.d6*rho(i),1.d16*u(i),                   
!     $                       1.d16*pr(i),1.d8*v(i),1.d9*temp(i)        
!1.d14*ufreez(i),1.d16*pr(i),1.d8*v(i)                                  
                     end if 
                     write(70,106) t*10.,deltam(i),2.d6*rho(i)          &
     &                    ,1.d9*temp(i)                                 
                  else 
                     iskip=iskip+1 
                  end if 
                  if (v(i).gt.vg) then 
                     sumni=sumni+deltam(i)*56.d0*ycc(i,16) 
                     sumzn=sumzn+deltam(i)*60.d0*ycc(i,17) 
                     sumfe=sumfe+deltam(i)*52.d0*ycc(i,15) 
                     sumcr=sumcr+deltam(i)*48.d0*ycc(i,14) 
                     sumti=sumti+deltam(i)*44.d0*ycc(i,13) 
                     sumca=sumca+deltam(i)*40.d0*ycc(i,12) 
                     sumar=sumar+deltam(i)*36.d0*ycc(i,11) 
                     sums=sums+deltam(i)*32.d0*ycc(i,10) 
                     sumsi=sumsi+deltam(i)*28.d0*ycc(i,9) 
                     summg=summg+deltam(i)*24.d0*ycc(i,8) 
                     sumne=sumne+deltam(i)*20.d0*ycc(i,7) 
                     sumo=sumo+deltam(i)*16.d0*ycc(i,6) 
                     sumc=sumc+deltam(i)*12.d0*ycc(i,5) 
                  end if 
!                  if (encm(i).gt.1.4506.and.encm(i).lt.1.4509) then    
!                     write(71,103) t*10,encm(i),1.d9*x(i),1.d8*v(i),   
!     $                 1.d8*vg,2.d6*rho(i),                            
!     $                 temp(i),ye(i),abar(i),                          
!     $                 28.*ycc(i,9),56.*ycc(i,16)                      
!                  end if                                               
                  if (k.eq.nc) then 
!                     write(70,103)encm(i),1.d9*x(i),1.d8*v(i),         
!     $                    2.d6*rho(i),0.0862*eta(i)*temp(i),           
!     $              etanue(i)*tempnue(i),                              
!     $              etanueb(i)*tempnueb(i),                            
!     $              etanux(i)*tempnux(i),                              
!--for aimee                                                            
!     $              etanue(i),                                         
!     $              etanueb(i),                                        
!     $              etanux(i),                                         
!                                                                       
!     $                    xn(i),xp(i),pr(i),temp(i)                    
                  end if 
               end if 
            end do 
            !print *, 'Nickel',sumc,sumo,sumne,summg,sumsi,sums,         &
!     &           sumar,sumca,sumti,sumcr,sumfe,sumni,sumzn              
         end if 
 !        print *, 'energy',dk/50.,dene/50. 
         deallocate(outname)
         print*, 'Dump: ', idump
      end do 
  103 format(I5,1pe12.4,1pe14.6,22(1pe12.4)) 
  105 format(1pe12.4,1pe14.6,7(1pe12.4)) 
  106 format(4(1pe14.7)) 
                                                                        
  107 format(I3,24(1pe13.5)) 
  108 format(1pe12.4,I5,1pe12.4,1pe12.4) 

      print*,' ==============================================='       
      print*, '  Converted ',trim(adjustl(infile)),' to ',trim(adjustl(basename))
!                                                                       
      END                                           